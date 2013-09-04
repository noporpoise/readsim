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
"usage: readsim [options] <profile.fq> <out.1.fq.gz> <out.2.fq.gz>\n"
" Simulate single base change sequencing errors\n"
"\n"
" Simulate reads:\n"
"  -r <ref.fa>  sample reads from ref\n"
"  -t <t>       template size (insert size + 2*(read length)) [800]\n"
"  -v <v>       variance on template size [100]\n"
"  -l <l>       read length [250]\n"
"  -d <d>       sequencing depth [1]\n"
" Load reads:\n"
"  -1 <in.1.fq> input reads\n"
"  -2 <in.2.fq> input reads\n";

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

// Returns length of first read
size_t load_error_profile(const char *path, read_t **rptr, size_t *nreadsptr)
{
  int qoffset, minq, maxq, fmt;

  size_t ncap = 16;
  read_t *r = malloc(ncap * sizeof(read_t));

  seq_file_t *file = seq_open(path);
  if(file == NULL) die("Cannot open file: %s", path);

  // Get FASTQ offset
  if((fmt = seq_guess_fastq_format(file, &minq, &maxq)) == -1)
    die("Cannot get quality score from: %s", path);

  qoffset = FASTQ_OFFSET[fmt];
  printf("FASTQ offset: %i\n", qoffset);

  if(!seq_read_alloc(&r[0])) die("Out of memory");
  if(seq_read(file, &r[0]) <= 0)
    die("Empty profile file %s", path);

  size_t i, len = r[0].seq.end, nreads;
  size_t sumcount = 0, *counts;
  double sumprob = 0, *qprob;

  qprob = malloc(len * sizeof(double));
  counts = calloc(len, sizeof(size_t));

  if(counts == NULL || qprob == NULL) die("Out of memory");

  for(i = 0; i < len; i++) qprob[i] = 0.0;

  sum_qual_scores(&r[0], qoffset, qprob, counts, len);

  for(nreads = 1; ; nreads++) {
    if(nreads == ncap) { r = realloc(r, (ncap *= 2) * sizeof(read_t)); }
    if(!seq_read_alloc(&r[nreads])) die("Out of memory");
    if(seq_read(file, &r[nreads]) <= 0) { seq_read_dealloc(&r[nreads]); break; }
    sum_qual_scores(&r[nreads], qoffset, qprob, counts, len);
  }

  printf("num profile reads: %zu\n", nreads);

  // Normalise qprob
  // double qmean = qtotal / len, err_rate = 0.01;
  // for(i = 0; i < len; i++)
  //   qprob[i] = err_rate * qprob[i] / qmean;

  // Print qprob
  printf("Counts: %zu", counts[0]);
  sumcount = counts[0];
  for(i = 1; i < len; i++) {
    printf(",%zu", counts[i]);
    sumcount += counts[i];
  }
  printf("\nCountsSum: %zu CountsMean: %zu\n", sumcount, sumcount / len);

  printf("QProb: %f", qprob[0] / counts[0]);
  for(i = 1; i < len; i++) {
    printf(",%f", qprob[i] / counts[i]);
    sumprob += qprob[i] / counts[i];
  }
  printf("\nQProbSum: %f QProbMean: %f\n", sumprob, sumprob / len);

  free(qprob);
  free(counts);

  *rptr = r;
  *nreadsptr = nreads;

  return len;
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
  read_t r;
  size_t rndread;
  seq_read_alloc(&r);

  while(seq_read(sfile, &r) > 0) {
    rndread = drand48() * (nreads-1);
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
                 size_t tlen, size_t tlen_stddev, size_t rlen, double depth)
{
  size_t i, chromcap = 16, nchroms, glen = 0, nreads, chr, pos0, pos1;
  read_t *chroms;

  chroms = malloc(chromcap * sizeof(read_t));
  nchroms = 0;

  while(1)
  {
    if(nchroms == chromcap) chroms = realloc(chroms, (chromcap *= 2) * sizeof(read_t));
    seq_read_alloc(&chroms[nchroms]);
    if(seq_read(reffile, &chroms[nchroms]) <= 0) { seq_read_dealloc(&chroms[nchroms]); break; }
    if(chroms[nchroms].seq.end < tlen) { seq_read_dealloc(&chroms[nchroms]); }
    else glen += chroms[nchroms++].seq.end;
  }

  if(nchroms == 0)
    die("No sequences long enough in ref genome file [min len: %zu]: %s",
        tlen, reffile->path);

  // Sample
  nreads = (glen * depth) / (2 * rlen);
  char read0[rlen+1], read1[rlen+1];
  read0[rlen] = read1[rlen] = '\0';

  printf("Sampling %zu read pairs\n", nreads);

  for(i = 0; i < nreads; i++)
  {
    chr = (nchroms == 1) ? 0 : rand_chrom(chroms, nchroms, glen);
    pos0 = drand48() * (chroms[chr].seq.end - tlen);
    pos1 = pos0 + tlen + ran_normal()*tlen_stddev - rlen;
    if(pos1 + rlen > chroms[chr].seq.end) pos1 = chroms[chr].seq.end-rlen;
    memcpy(read0, chroms[chr].seq.b+pos0, rlen);
    memcpy(read1, chroms[chr].seq.b+pos1, rlen);
    add_seq_error(read0, rlen, &rlist[(int)(drand48() * (rlistlen-1))]);
    add_seq_error(read1, rlen, &rlist[(int)(drand48() * (rlistlen-1))]);
    dna_revcmp(read1, rlen);
    gzprintf(out0, ">r%zu:0:%s:%zu:%zu\n%.*s\n", i, chroms[chr].name.b, pos0, pos1, (int)rlen, read0);
    gzprintf(out1, ">r%zu:1:%s:%zu:%zu\n%.*s\n", i, chroms[chr].name.b, pos0, pos1, (int)rlen, read1);
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
  int tlen = 800, tlen_stddev = 100, rlen = 250;
  double depth = 1;
  int optr = 0, optt = 0, optv = 0, optl = 0, optd = 0; // keeps track of values

  char *in0path = NULL, *in1path = NULL;

  int c;
  while((c = getopt(argc, argv, "r:t:v:l:d:1:2:")) >= 0) {
    switch (c) {
      case 'r': refpath = optarg; optr++; break;
      case 't': tlen = atoi(optarg); optt++; break;
      case 'v': tlen_stddev = atoi(optarg); optv++; break;
      case 'l': rlen = atoi(optarg); optl++; break;
      case 'd': depth = atol(optarg); optd++; break;
      case '1': in0path = optarg; break;
      case '2': in1path = optarg; break;
      default: die("Unknown option: %c", c);
    }
  }

  if(argc - optind != 3) print_usage(usage, "Missing args");

  if((optt > 0 || optv > 0 || optl > 0 || optd > 0) && refpath == NULL)
    print_usage(usage, "Missing -r <in.fa>");

  if(optr > 1 || optt > 1 || optv > 1 || optl > 1 || optd > 1)
    print_usage(usage, "Duplicate args");

  if((in0path == NULL) != (in1path == NULL))
    print_usage(usage, "Need both -1 <in> -2 <in>");

  char *out0path = argv[optind+1], *out1path = argv[optind+2];  
  gzFile gzout0 = NULL, gzout1 = NULL;
  seq_file_t *sf0 = NULL, *sf1 = NULL, *reffile = NULL;

  if(in0path != NULL)
  {
    printf("Reading from %s and %s\n", in0path, in1path);
    if((sf0 = seq_open(in0path)) == NULL) die("Cannot read: %s", in0path);
    if((sf1 = seq_open(in1path)) == NULL) die("Cannot read: %s", in1path);
  }

  if(refpath != NULL)
  {
    printf("Sampling from %s\n", refpath);
    printf(" read length: %i\n", rlen);
    printf(" template length: %i\n", tlen);
    printf(" template stddev: %i\n", tlen_stddev);
    printf(" sequencing depth: %f\n", depth);
    if((reffile = seq_open(refpath)) == NULL) die("Cannot read: %s", refpath);
    if((gzout0 = gzopen(out0path, "w")) == NULL) die("Cannot open: %s", out0path);
    if((gzout1 = gzopen(out1path, "w")) == NULL) die("Cannot open: %s", out1path);
  }

  // Generate profile from profile.fq
  size_t errlen, glen, rlistlen = 0;
  read_t *rlist;

  errlen = load_error_profile(argv[optind], &rlist, &rlistlen);

  if(sf0 != NULL) {
    mutate_reads(sf0, gzout0, rlist, rlistlen);
    mutate_reads(sf1, gzout1, rlist, rlistlen);
    seq_close(sf0);
    seq_close(sf1);
  }

  if(refpath != NULL) {
    glen = sim_reads(reffile, gzout0, gzout1, rlist, rlistlen,
                     tlen, tlen_stddev, rlen, depth);
    seq_close(reffile);
  }

  if(gzout0 != NULL) {
    gzclose(gzout0);
    gzclose(gzout1);
  }

  free(rlist);

  return EXIT_SUCCESS;
}
