#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <inttypes.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <zlib.h>
#include <time.h>

#include "seq_file.h"

static const char usage[] =
"usage: readsim [options] <out_base>\n"
" Simulate single base change sequencing errors. \n"
" Writes to <out_base>.1.fa.gz <out_base>.1.fa.gz\n"
"\n"
" Apply Error:\n"
"  -p <profile.fq> apply errors with distribution as seen in profile.fq\n"
" Simulate reads:\n"
"  -r <ref.fa>  sample reads from ref\n"
"  -t <t>       template size (insert size + 2*(read length)) [800]\n"
"  -v <v>       variance on template size as a proportion [0.1]\n"
"  -l <l>       read length [250]\n"
"  -d <d>       sequencing depth [1]\n"
" Load reads:\n"
"  -1 <in.1.fq> input reads\n"
"  -2 <in.2.fq> input reads\n";

#define MAX2(x,y) ((x) >= (y) ? (x) : (y))
#define MIN2(x,y) ((x) <= (y) ? (x) : (y))

// Sample a random number from normal distribution with [mean 0, stddev 1]
// From Heng Li
// http://en.wikipedia.org/wiki/Box_Muller_transform#Polar_form
// http://en.wikipedia.org/wiki/Marsaglia_polar_method
double ran_normal()
{ 
  static int iset = 0; 
  static double gset; 
  double fac, rsq, v1, v2; 
  if (iset == 0) {
    do { 
      v1 = 2.0 * drand48() - 1.0;
      v2 = 2.0 * drand48() - 1.0; 
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0 * log(rsq) / rsq); 
    gset = v1 * fac; 
    iset = 1;
    return v2 * fac;
  } else {
    iset = 0;
    return gset;
  }
}

void print_usage(const char *msg, const char *errfmt,  ...)
{
  if(errfmt != NULL) {
    fprintf(stderr, "Error: ");
    va_list argptr;
    va_start(argptr, errfmt);
    vfprintf(stderr, errfmt, argptr);
    va_end(argptr);

    if(errfmt[strlen(errfmt)-1] != '\n') fputc('\n', stderr);
  }

  fputs(msg, stderr);
  exit(EXIT_FAILURE);
}

#define die(fmt, ...) call_die(__FILE__, __LINE__, fmt, ##__VA_ARGS__)

void call_die(const char *file, int line, const char *fmt, ...)
__attribute__((format(printf, 3, 4)))
__attribute__((noreturn));

void call_die(const char *file, int line, const char *fmt, ...)
{
  va_list argptr;
  fflush(stdout);
  fprintf(stderr, "[%s:%i] Error: ", file, line);
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);
  if(*(fmt + strlen(fmt) - 1) != '\n') fputc('\n', stderr);
  exit(EXIT_FAILURE);
}

char test_file_readable(const char *file)
{
  FILE *fp = fopen(file, "r");
  if(fp == NULL) return 0;
  fclose(fp);
  return 1;
}

#define phred_to_prob(qual) powl(10, -(double)(qual)/10)

void sum_qual_scores(const read_t *read, int qoffset,
                     double *qprob_sum, size_t *counts, size_t len)
{
  size_t i, limit = len < read->qual.end ? len : read->qual.end;
  double qual;

  if(read->qual.end != read->seq.end) die("Invalid qual");

  for(i = 0; i < read->qual.end; i++) {
    if(read->qual.b[i] < qoffset) die("Invalid qual");
    if(read->qual.b[i] - qoffset > 50) die("Invalid qual");
    read->qual.b[i] -= qoffset;
  }

  for(i = 0; i < limit; i++) {
    if(toupper(read->seq.b[i]) != 'N') {
      qual = phred_to_prob(read->qual.b[i]);
      qprob_sum[i] += qual;
      counts[i]++;
    }
  }
}

size_t extend_arrays(size_t currlen, size_t newlen, double **qprob, size_t **counts)
{
  size_t i;
  *qprob = realloc(*qprob, newlen * sizeof(double));
  *counts = realloc(*counts, newlen * sizeof(size_t));
  if(*counts == NULL || *qprob == NULL) die("Out of memory");
  for(i = currlen; i < newlen; i++) { (*qprob)[i] = 0; (*counts)[i] = 0; }
  return newlen;
}

// Returns length of first read
void load_error_profile(seq_file_t **files, size_t num_files,
                        read_t **rptr, size_t *nreadsptr)
{
  int qoffset, minq, maxq, fmt;

  size_t i, f, longest_read = 0;
  size_t *counts = NULL;
  double *qprob = NULL;
  size_t sumlen = 0, meanlen = 0, maxlen = 0, minlen = SIZE_MAX;

  size_t rcap = 16, rlen = 0;
  read_t *r = malloc(rcap * sizeof(read_t));

  printf("Loading calibration reads...\n");

  for(f = 0; f < num_files; f++)
  {
    // Get FASTQ offset
    if((fmt = seq_guess_fastq_format(files[f], &minq, &maxq)) == -1)
      die("Cannot get quality score from: %s", files[f]->path);

    qoffset = FASTQ_OFFSET[fmt];
    printf(" reading: %s  [FASTQ offset: %i]\n", files[f]->path, qoffset);

    if(!seq_read_alloc(&r[rlen])) die("Out of memory");
    if(seq_read(files[f], &r[rlen]) <= 0)
      die("Empty profile file %s", files[f]->path);

    if(counts == NULL) {
      longest_read = r[rlen].seq.end;
      qprob = malloc(longest_read * sizeof(double));
      counts = calloc(longest_read, sizeof(size_t));
      if(counts == NULL || qprob == NULL) die("Out of memory");
      for(i = 0; i < longest_read; i++) qprob[i] = 0.0;
    }
    else if(r[rlen].seq.end > longest_read)
      longest_read = extend_arrays(longest_read, r[rlen].seq.end, &qprob, &counts);

    sum_qual_scores(&r[rlen], qoffset, qprob, counts, r[rlen].seq.end);
    sumlen += r[rlen].seq.end;
    maxlen = MAX2(maxlen, r[rlen].seq.end);
    minlen = MIN2(minlen, r[rlen].seq.end);
    rlen++;

    while(1) {
      if(rlen == rcap) { r = realloc(r, (rcap *= 2) * sizeof(read_t)); }
      if(!seq_read_alloc(&r[rlen])) die("Out of memory");
      if(seq_read(files[f], &r[rlen]) <= 0) { seq_read_dealloc(&r[rlen]); break; }
      if(r[rlen].seq.end > longest_read)
        longest_read = extend_arrays(longest_read, r[rlen].seq.end, &qprob, &counts);
      sum_qual_scores(&r[rlen], qoffset, qprob, counts, r[rlen].seq.end);
      sumlen += r[rlen].seq.end;
      maxlen = MAX2(maxlen, r[rlen].seq.end);
      minlen = MIN2(minlen, r[rlen].seq.end);
      rlen++;
    }
  }

  meanlen = sumlen / rlen;

  printf("Number of calibration reads: %zu\n", rlen);
  printf("MeanReadLen: %zu; MinReadLen: %zu; MaxReadLen: %zu\n",
         meanlen, minlen, maxlen);

  // Normalise qprob
  // double qmean = qtotal / len, err_rate = 0.01;
  // for(i = 0; i < len; i++)
  //   qprob[i] = err_rate * qprob[i] / qmean;

  if(rlen == 0) {
    printf("Counts: \nCountsMean: 0\n");
    printf("QProb: \nQProbMean: 0\n");
  }
  else
  {
    size_t sumcount;
    double sumprob;

    sumcount = counts[0];
    printf("Counts: %zu", counts[0]);
    for(i = 1; i < longest_read; i++) {
      printf(",%zu", counts[i]);
      sumcount += counts[i];
    }
    printf("\nCountsMean: %zu\n", sumcount / longest_read);

    sumprob = qprob[0];
    printf("QProb: %f", qprob[0] / counts[0]);
    for(i = 1; i < longest_read; i++) {
      printf(",%f", qprob[i] / counts[i]);
      sumprob += qprob[i];
    }
    printf("\nQProbMean: %f\n", sumprob / sumcount);
  }

  free(qprob);
  free(counts);

  *rptr = r;
  *nreadsptr = rlen;
}

static inline char dna_complement(char c)
{
  switch(c) {
    case 'a': return 't';
    case 'c': return 'g';
    case 'g': return 'c';
    case 't': return 'a';
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    case 'N': return 'N';
    case 'n': return 'n';
    default: die("Invalid base: '%c'", c);
  }
  die("Invalid base: '%c'", c);
}

static void dna_revcmp(char *dna, size_t len)
{
  if(len == 0) return;
  if(len == 1) { dna[0] = dna_complement(dna[0]); return; }
  size_t i, j;
  char tmp;
  for(i = 0, j = len-1; i <= j; i++, j--) {
    tmp = dna[i];
    dna[i] = dna_complement(dna[j]);
    dna[j] = dna_complement(tmp);
  }
}

// given base c and random i (0 <= i < 3)
static inline char mut(char c, int i)
{
  const char bases[] = "ACGTACGT";
  int b;
  switch(c) {
    case 'A': case 'a': b = 1; break;
    case 'C': case 'c': b = 2; break;
    case 'G': case 'g': b = 3; break;
    case 'T': case 't': b = 4; break;
    default: die("Invalid base: '%c'", c);
  }
  return bases[b+i];
}

void add_seq_error(char *seq, size_t seqlen, const read_t *r)
{
  size_t i, limit = seqlen < r->seq.end ? seqlen : r->seq.end;
  int rnd;
  for(i = 0; i < limit; i++) {
    if(r->seq.b[i] == 'N') seq[i] = 'N';
    else if((rnd = rand()) < phred_to_prob(r->qual.b[i]) * RAND_MAX)
      seq[i] = mut(seq[i], rnd % 3);
  }
}

// Load reads from a file, apply sequence error, dump
void mutate_reads(seq_file_t *sfile, gzFile gzout,
                  const read_t *rlist, size_t nreads)
{
  printf(" reading: %s\n", sfile->path);
  read_t r;
  size_t rndread;
  seq_read_alloc(&r);

  while(seq_read(sfile, &r) > 0) {
    rndread = drand48() * nreads;
    add_seq_error(r.seq.b, r.seq.end, &rlist[rndread]);
    gzprintf(gzout, "@%s\n%s\n+\n%s\n", r.name.b, r.seq.b, r.qual.b);
  }

  seq_read_dealloc(&r);
}

size_t rand_chrom(read_t *chroms, size_t nchroms, size_t totallen)
{
  uint64_t i, r = (((uint64_t)rand()) << 32) | rand();
  r %= totallen;
  for(i = 0; i < nchroms; i++)
    if(r < chroms[i].seq.end)
      return i;
  die("Shouldn't reach here");
}

// Returns ref genome size
size_t sim_reads(seq_file_t *reffile, gzFile out0, gzFile out1,
                 const read_t *rlist, size_t rlistlen,
                 size_t tlen, double tlen_stddev, size_t rlen, double depth)
{
  size_t i, chromcap = 16, nchroms, glen = 0, nreads, chr, pos0, pos1;
  read_t *chroms;

  chroms = malloc(chromcap * sizeof(read_t));
  nchroms = 0;

  // Load genome
  printf(" Loaded contigs:");
  while(1)
  {
    if(nchroms == chromcap) chroms = realloc(chroms, (chromcap *= 2) * sizeof(read_t));
    seq_read_alloc(&chroms[nchroms]);
    if(seq_read(reffile, &chroms[nchroms]) <= 0) { seq_read_dealloc(&chroms[nchroms]); break; }
    if(chroms[nchroms].seq.end < tlen) { seq_read_dealloc(&chroms[nchroms]); }
    else {
      seq_read_truncate_name(&chroms[nchroms]);
      printf(" %s[%zu]", chroms[nchroms].name.b, chroms[nchroms].seq.end);
      glen += chroms[nchroms].seq.end;
      nchroms++;
    }
  }
  printf("\n Genome size: %zu\n", glen);

  if(nchroms == 0) {
    die("No sequences long enough in ref genome file [min len: %zu]: %s",
        tlen, reffile->path);
  }

  // Sample
  nreads = (glen * depth) / (out1 == NULL ? rlen : (2 * rlen));
  char read0[rlen+1], read1[rlen+1];
  read0[rlen] = read1[rlen] = '\0';

  printf("Sampling %zu %sreads...\n", nreads,
         out1 == NULL ? "single " : "paired-end ");

  // Sample paired-end if out1 != NULL
  for(i = 0; i < nreads; i++)
  {
    chr = (nchroms == 1) ? 0 : rand_chrom(chroms, nchroms, glen);
    pos0 = drand48() * (chroms[chr].seq.end - (out1 == NULL ? rlen : tlen));
    pos1 = pos0;
    memcpy(read0, chroms[chr].seq.b+pos0, rlen);
    if(out1 != NULL) {
      pos1 = pos0 + tlen + ran_normal()*tlen*tlen_stddev - rlen;
      if(pos1 + rlen > chroms[chr].seq.end) pos1 = chroms[chr].seq.end-rlen;
      memcpy(read1, chroms[chr].seq.b+pos1, rlen);
    }
    if(rlistlen > 0) {
      add_seq_error(read0, rlen, &rlist[(int)(drand48() * rlistlen)]);
      if(out1 != NULL)
        add_seq_error(read1, rlen, &rlist[(int)(drand48() * rlistlen)]);
    }
    gzprintf(out0, ">r%zu:0:%s:%zu:%zu\n%.*s\n", i, chroms[chr].name.b,
                   pos0, pos1, (int)rlen, read0);
    if(out1 != NULL) {
      dna_revcmp(read1, rlen);
      gzprintf(out1, ">r%zu:1:%s:%zu:%zu\n%.*s\n", i, chroms[chr].name.b,
                     pos0, pos1, (int)rlen, read1);
    }
  }

  for(i = 0; i < nchroms; i++) seq_read_dealloc(&chroms[i]);
  free(chroms);

  return glen;
}

int main(int argc, char **argv)
{
  if(argc < 4) print_usage(usage, NULL);
  srand(time(NULL) + getpid());

  // Sample reads from ref
  char *refpath = NULL;
  int tlen = 800, rlen = 250, single_ended = 0;
  double depth = 1.0, tlen_stddev = 0.1; // stddev as proportion of tlen
  int optr = 0, optt = 0, optv = 0, optl = 0, optd = 0; // keeps track of values

  char *in0path = NULL, *in1path = NULL;

  char *profile_paths[argc];
  int i, num_profile_paths = 0;

  int c;
  while((c = getopt(argc, argv, "p:r:t:v:l:d:s1:2:")) >= 0) {
    switch (c) {
      case 'p': profile_paths[num_profile_paths++] = optarg; break;
      case 'r': refpath = optarg; optr++; break;
      case 't': tlen = atoi(optarg); optt++; break;
      case 'v': tlen_stddev = atof(optarg); optv++; break;
      case 'l': rlen = atoi(optarg); optl++; break;
      case 'd': depth = atof(optarg); optd++; break;
      case 's': single_ended = 1; break;
      case '1': in0path = optarg; break;
      case '2': in1path = optarg; break;
      default: die("Unknown option: %c", c);
    }
  }

  char *outbase = argv[optind];

  if(argc == optind) print_usage(usage, "Missing <out_base>");
  if(argc > optind + 1) print_usage(usage, "Too many args after %s", outbase);

  if(depth <= 0) print_usage(usage, "Depth [-d] cannot be <= 0");

  if(tlen_stddev < 0)
    print_usage(usage, "Template length standard deviation [-v] cannot be < 0");

  if((optt > 0 || optv > 0 || optl > 0 || optd > 0) && refpath == NULL)
    print_usage(usage, "Missing -r <in.fa>");

  if(optr > 1 || optt > 1 || optv > 1 || optl > 1 || optd > 1)
    print_usage(usage, "Duplicate args");

  if(in0path == NULL && in1path != NULL)
    print_usage(usage, "-2 <in> requires -1 <in>");

  if(in0path != NULL && in1path == NULL) {
    if(refpath == NULL) single_ended = 1;
    else if(!single_ended) print_usage(usage, "Missing -2 for paired-end output");
  }

  if(in0path != NULL && num_profile_paths == 0)
    print_usage(usage, "Need at least one -p <profile.fq.gz> to use -1 .. -2 ..");

  if(num_profile_paths == 0 && refpath == NULL)
    print_usage(usage, "Need one of -p or -r");

  seq_file_t *profile_sf[num_profile_paths];

  for(i = 0; i < num_profile_paths; i++)
    if((profile_sf[i] = seq_open(profile_paths[i])) == NULL)
      die("Cannot open file: %s", profile_paths[i]);

  size_t outlen = strlen(outbase), extlen = strlen(".1.fa.gz");
  char out0path[outlen+extlen+1], out1path[outlen+extlen+1];
  memcpy(out0path, outbase, outlen);
  memcpy(out1path, outbase, outlen);

  if(single_ended) strcpy(out0path+outlen, ".fa.gz");
  else {
    strcpy(out0path+outlen, ".1.fa.gz");
    strcpy(out1path+outlen, ".2.fa.gz");
  }

  gzFile gzout0 = NULL, gzout1 = NULL;
  seq_file_t *sf0 = NULL, *sf1 = NULL, *reffile = NULL;

  if(in0path != NULL && (sf0 = seq_open(in0path)) == NULL) die("Cannot read: %s", in0path);
  if(in1path != NULL && (sf1 = seq_open(in1path)) == NULL) die("Cannot read: %s", in1path);

  if(refpath != NULL)
  {
    if((reffile = seq_open(refpath)) == NULL) die("Cannot read: %s", refpath);
    if((gzout0 = gzopen(out0path, "w")) == NULL) die("Cannot open: %s", out0path);
    if(!single_ended && (gzout1 = gzopen(out1path, "w")) == NULL)
      die("Cannot open: %s", out1path);
  }

  // Generate profile from profile.fq
  size_t r, rlistlen = 0;
  read_t *rlist = NULL;

  // Load existing `profile' reads
  if(num_profile_paths > 0)
  {
    load_error_profile(profile_sf, num_profile_paths, &rlist, &rlistlen);

    for(i = 0; i < num_profile_paths; i++)
      seq_close(profile_sf[i]);
  }

  if(sf0 != NULL) printf("Adding error to input reads...\n");
  if(sf0 != NULL) { mutate_reads(sf0, gzout0, rlist, rlistlen); seq_close(sf0); }
  if(sf1 != NULL) {
    mutate_reads(sf1, single_ended ? gzout0 : gzout1, rlist, rlistlen);
    seq_close(sf1);
  }

  if(refpath != NULL)
  {
    printf("Sampling from %s\n", refpath);
    printf(" read length: %i\n", rlen);
    printf(" template length: %i\n", tlen);
    printf(" template stddev: %.2f * tlen = %.2f\n", tlen_stddev, tlen_stddev*tlen);
    printf(" sequencing depth: %.2f\n", depth);
    printf(" read pairs: %s\n", single_ended ? "no" : "yes");
    sim_reads(reffile, gzout0, gzout1, rlist, rlistlen,
              tlen, tlen_stddev, rlen, depth);
    seq_close(reffile);
  }

  if(gzout1 != NULL && gzout1 != NULL)
    printf("Wrote to: %s and %s\n", out0path, out1path);
  else if(gzout0 != NULL)
    printf("Wrote to: %s\n", out0path);

  if(gzout0 != NULL) gzclose(gzout0);
  if(gzout1 != NULL) gzclose(gzout1);

  if(rlist != NULL) {
    for(r = 0; r < rlistlen; r++) seq_read_dealloc(&rlist[r]);
    free(rlist);
  }

  return EXIT_SUCCESS;
}
