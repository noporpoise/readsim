#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <inttypes.h>
#include <string.h>
#include <strings.h>
#include <limits.h>
#include <math.h>
#include <zlib.h>
#include <time.h>
#include <sys/time.h> // for seeding random

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
// "  -t <t>       template size (insert size + 2*(read length)) [800]\n"
// "  -v <v>       variance on template size as a proportion [0.1]\n"
"  -i <t>       insert size [250]\n"
"  -v <v>       variance on insert size as a proportion [0.2 => 50bp]\n"
"  -l <l>       read length [250]\n"
"  -d <d>       sequencing depth [1]\n"
" Load reads:\n"
"  -1 <in.1.fq> input reads\n"
"  -2 <in.2.fq> input reads\n";

#define MAX2(x,y) ((x) >= (y) ? (x) : (y))
#define MIN2(x,y) ((x) <= (y) ? (x) : (y))

#define ROUNDUP2POW(x) (0x1UL << (64 - __builtin_clzl(x)))

#define phred_to_prob(qual) powl(10, -(double)(qual)/10.0)
#define phred_to_lprob(qual) powl(10, -(long double)(qual)/10.0)

#define MAXQUAL 60
double qual_prob[MAXQUAL];

static void init_qual_prob()
{
  size_t i;
  for(i = 0; i < MAXQUAL; i++) qual_prob[i] = phred_to_prob(i);
}

// Sample a random number from normal distribution with [mean 0, stddev 1]
// N(mean,stddev) = mean + ran_normal()*stddev
// http://en.wikipedia.org/wiki/Box_Muller_transform#Polar_form
// http://en.wikipedia.org/wiki/Marsaglia_polar_method
// From Heng Li
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

typedef struct
{
  seq_file_t **files;
  int *fqoffsets;
  size_t num_files, curr, filesready;
  read_t read;
  size_t *errors, errors_len, errors_cap;
} FileList;

void filelist_alloc(FileList *flist, char **paths, size_t num)
{
  size_t i;
  flist->num_files = num;
  flist->curr = 0;
  flist->files = malloc(num * sizeof(seq_file_t*));
  flist->fqoffsets = malloc(num * sizeof(int));

  for(i = 0; i < num; i++) {
    if((flist->files[i] = seq_open(paths[i])) == NULL)
      die("Cannot open: %s", paths[i]);
    int min, max, fmt;
    fmt = seq_guess_fastq_format(flist->files[i], &min, &max);
    flist->fqoffsets[i] = FASTQ_OFFSET[fmt];
    printf(" profile: %s [offset: %i]\n", paths[i], FASTQ_OFFSET[fmt]);
  }

  seq_read_alloc(&flist->read);
  flist->filesready = 1;
  flist->errors_cap = 512;
  flist->errors_len = 0;
  flist->errors = calloc(flist->errors_cap, sizeof(size_t));
}

void filelist_dealloc(FileList *flist)
{
  size_t i;
  for(i = 0; i < flist->num_files; i++) seq_close(flist->files[i]);
  seq_read_dealloc(&flist->read);
  free(flist->files);
  free(flist->fqoffsets);
  free(flist->errors);
}

read_t* filelist_read(FileList *flist)
{
  read_t *r = &flist->read;
  size_t i; // i is number of file changes
  for(i = 0; seq_read(flist->files[flist->curr], r) <= 0 && i <= flist->num_files; i++)
  {
    flist->curr++;
    if(flist->curr == flist->num_files) { flist->curr = flist->filesready = 0; }
    if(!flist->filesready) {
      char path[PATH_MAX+1];
      strcpy(path, flist->files[flist->curr]->path);
      seq_close(flist->files[flist->curr]);
      seq_open(path);
    }
  }
  if(i > flist->num_files) die("All seq files empty");
  return r;
}

void filelist_extend_errarr(FileList *flist, size_t len)
{
  size_t i, oldcap = flist->errors_cap;
  flist->errors_len = MAX2(flist->errors_len, len);
  if(len > flist->errors_cap) {
    flist->errors_cap = ROUNDUP2POW(len);
    flist->errors = realloc(flist->errors, flist->errors_cap * sizeof(size_t));
    for(i = oldcap; i < flist->errors_cap; i++) flist->errors[i] = 0;
  }
}

void filelist_mean_err(FileList *flist)
{
  size_t i, f, maxread = 0, readcap = 512, newcap, carry;
  long double *sumprob = malloc(readcap * sizeof(long double));
  size_t *counts = malloc(readcap * sizeof(size_t));
  for(i = 0; i < readcap; i++) { sumprob[i] = 0; counts[i] = 0; }

  read_t *r = &flist->read;
  int fmt, fqoffset = 33, minq, maxq;
  for(f = flist->curr; f < flist->num_files; f++) {
    fmt = seq_guess_fastq_format(flist->files[f], &minq, &maxq);
    fqoffset = (fmt == -1 ? 33 : FASTQ_OFFSET[fmt]);
    while(seq_read(flist->files[f], r) > 0) {
      if(r->qual.end > readcap) {
        newcap = ROUNDUP2POW(r->qual.end);
        sumprob = realloc(sumprob, newcap * sizeof(long double));
        counts = realloc(counts, newcap * sizeof(size_t));
        for(i = readcap; i < newcap; i++) { sumprob[i] = 0; counts[i] = 0; }
        readcap = newcap;
      }
      counts[r->qual.end-1]++;
      for(i = 0; i < r->qual.end; i++)
        sumprob[i] += qual_prob[r->qual.b[i] - fqoffset];
      maxread = MAX2(maxread, r->qual.end);
    }
  }

  // Convert counts to cummulative (reverse)
  for(i = maxread-1, carry = 0; i != SIZE_MAX; i--) {
    carry += counts[i];
    counts[i] = carry;
  }

  for(i = 0; i < maxread; i++) {
    // printf(" %.8Lf/%zu", sumprob[i], counts[i]);
    printf(" %.2Lf", 100.0 * sumprob[i] / counts[i]);
  } printf("\n");

  free(counts);
  free(sumprob);
}

static inline char dna_complement(char c)
{
  switch(c) {
    case 'a': return 't'; case 'A': return 'T';
    case 'c': return 'g'; case 'C': return 'G';
    case 'g': return 'c'; case 'G': return 'C';
    case 't': return 'a'; case 'T': return 'A';
    case 'n': return 'n'; case 'N': return 'N';
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

void add_seq_error(char *seq, size_t seqlen, FileList *flist)
{
  const read_t *r = filelist_read(flist);
  int fqoffset = flist->fqoffsets[flist->curr];
  filelist_extend_errarr(flist, seqlen);
  size_t i, limit = MIN2(seqlen, r->seq.end);
  int rnd;
  for(i = 0; i < limit; i++) {
    if(toupper(seq[i]) == 'N' || toupper(r->seq.b[i]) == 'N') seq[i] = 'N';
    else if((rnd = rand()) < qual_prob[r->qual.b[i] - fqoffset] * RAND_MAX) {
      seq[i] = mut(seq[i], rnd % 3);
      flist->errors[i]++;
    }
  }
}

// Load reads from a file, apply sequence error, dump
// Return total number of bases
size_t mutate_reads(seq_file_t *sfile, gzFile gzout, FileList *flist)
{
  printf(" reading: %s\n", sfile->path);
  read_t r;
  seq_read_alloc(&r);
  size_t num_bases = 0;

  while(seq_read(sfile, &r) > 0) {
    add_seq_error(r.seq.b, r.seq.end, flist);
    gzprintf(gzout, "@%s\n%s\n+\n%s\n", r.name.b, r.seq.b, r.qual.b);
    num_bases += r.seq.end;
  }

  seq_read_dealloc(&r);
  return num_bases;
}

static size_t rand_chrom(read_t *chroms, size_t nchroms, size_t totallen)
{
  uint64_t i, r = (((uint64_t)rand()) << 32) | rand();
  r %= totallen;
  for(i = 0; i < nchroms; i++)
    if(r < chroms[i].seq.end)
      return i;
  die("Shouldn't reach here");
}

// Returns num of bases printed
size_t sim_reads(seq_file_t *reffile, gzFile out0, gzFile out1,
                 FileList *flist,
                 size_t insert, double insert_stddev, size_t rlen, double depth)
{
  size_t i, chromcap = 16, nchroms, glen = 0, nreads, chr, pos0, pos1, tlen;
  read_t *chroms;

  tlen = rlen + (out1 == NULL ? 0 : insert + rlen);

  chroms = malloc(chromcap * sizeof(read_t));
  nchroms = 0;

  // Load genome
  printf(" Loaded contigs:");
  while(1)
  {
    if(nchroms == chromcap) chroms = realloc(chroms, (chromcap*=2)*sizeof(read_t));
    seq_read_alloc(&chroms[nchroms]);
    if(seq_read(reffile, &chroms[nchroms]) <= 0)
    { seq_read_dealloc(&chroms[nchroms]); break; }
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
      pos1 = pos0 + rlen + insert + ran_normal()*insert_stddev;
      if(pos1 + rlen > chroms[chr].seq.end) pos1 = chroms[chr].seq.end-rlen;
      memcpy(read1, chroms[chr].seq.b+pos1, rlen);
    }
    if(flist != NULL) {
      add_seq_error(read0, rlen, flist);
      if(out1 != NULL)
        add_seq_error(read1, rlen, flist);
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

  size_t num_bases = nreads * rlen;
  if(out1 != NULL) num_bases *= 2;

  return num_bases;
}

static void seed_random()
{
  struct timeval time;
  gettimeofday(&time, NULL);
  srand((((time.tv_sec ^ getpid()) * 1000000) + time.tv_usec));
}

int main(int argc, char **argv)
{
  if(argc < 3) print_usage(usage, NULL);

  // Set up
  seed_random();
  init_qual_prob();

  // Sample reads from ref
  char *refpath = NULL;
  // int optt = 0, tlen = 800; double tlen_stddev = 0.1;
  int insert = 250, rlen = 250, single_ended = 0;
  double depth = 1.0, insert_stddev_prop = 0.2; // stddev as proportion of insert
  int optr = 0, opti = 0, optv = 0, optl = 0, optd = 0; // keeps track of values

  char *in0path = NULL, *in1path = NULL;

  char *profile_paths[argc];
  size_t num_profile_paths = 0, i, total_seq = 0;

  int c;
  while((c = getopt(argc, argv, "p:r:i:v:l:d:s1:2:")) >= 0) {
    switch (c) {
      case 'p': profile_paths[num_profile_paths++] = optarg; break;
      case 'r': refpath = optarg; optr++; break;
      // case 't': tlen = atoi(optarg); optt++; break;
      // case 'v': tlen_stddev = atof(optarg); optv++; break;
      case 'i': insert = atoi(optarg); opti++; break;
      case 'v': insert_stddev_prop = atof(optarg); optv++; break;
      case 'l': rlen = atoi(optarg); optl++; break;
      case 'd': depth = atof(optarg); optd++; break;
      case 's': single_ended = 1; break;
      case '1': in0path = optarg; break;
      case '2': in1path = optarg; break;
      default: die("Unknown option: %c", c);
    }
  }

  char *outbase = NULL;

  if(optind == argc) {}//print_usage(usage, "Missing <out_base>");
  else if(optind + 1 == argc) outbase = argv[optind];
  else if(optind + 1 < argc) print_usage(usage, "Too many args after %s", outbase);

  if(depth <= 0) print_usage(usage, "Depth [-d] cannot be <= 0");

  if(insert_stddev_prop < 0)
    print_usage(usage, "Insert length standard deviation [-v] cannot be < 0");

  if((opti > 0 || optv > 0 || optl > 0 || optd > 0) && refpath == NULL)
    print_usage(usage, "Missing -r <in.fa>");

  if(optr > 1 || opti > 1 || optv > 1 || optl > 1 || optd > 1)
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

  if(num_profile_paths == 0 && outbase == NULL)
    print_usage(usage, "More options required");

  // Profile reads
  FileList fliststore, *flist = NULL;
  if(num_profile_paths > 0) {
    flist = &fliststore;
    filelist_alloc(flist, profile_paths, num_profile_paths);
  }

  if(outbase == NULL)
  {
    // Summarise error profile in input
    filelist_mean_err(flist);
  }
  else
  {
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
    size_t total_seq = 0;

    if(in0path != NULL && (sf0 = seq_open(in0path)) == NULL) die("Cannot read: %s", in0path);
    if(in1path != NULL && (sf1 = seq_open(in1path)) == NULL) die("Cannot read: %s", in1path);

    if(refpath != NULL)
    {
      if((reffile = seq_open(refpath)) == NULL) die("Cannot read: %s", refpath);
      if((gzout0 = gzopen(out0path, "w")) == NULL) die("Cannot open: %s", out0path);
      if(!single_ended && (gzout1 = gzopen(out1path, "w")) == NULL)
        die("Cannot open: %s", out1path);
    }

    if(sf0 != NULL) {
      printf("Adding error to input reads...\n");
      total_seq += mutate_reads(sf0, gzout0, flist);
      seq_close(sf0);
    }
    if(sf1 != NULL) {
      total_seq += mutate_reads(sf1, single_ended ? gzout0 : gzout1, flist);
      seq_close(sf1);
    }

    if(refpath != NULL)
    {
      printf("Sampling from %s\n", refpath);
      printf(" read length: %i\n", rlen);
      printf(" insert length: %i\n", insert);
      printf(" insert stddev: %.2f * insert = %.2f\n",
             insert_stddev_prop, insert_stddev_prop*insert);
      printf(" sequencing depth: %.2f\n", depth);
      printf(" read pairs: %s\n", single_ended ? "no" : "yes");
      if(num_profile_paths == 0) printf(" sequencing errors: no\n");
      else {
        printf(" seq error files: %s", flist->files[0]->path);
        for(i = 1; i < num_profile_paths; i++)
          printf(",%s", flist->files[i]->path);
        printf("\n");
      }
      total_seq += sim_reads(reffile, gzout0, gzout1, flist,
                             insert, insert_stddev_prop*insert, rlen, depth);
      seq_close(reffile);
    }

    if(gzout1 != NULL && gzout1 != NULL)
      printf("Wrote %zu bases to: %s and %s\n", total_seq, out0path, out1path);
    else if(gzout0 != NULL)
      printf("Wrote %zu bases to: %s\n", total_seq, out0path);

    if(gzout0 != NULL) gzclose(gzout0);
    if(gzout1 != NULL) gzclose(gzout1);
  }

  if(flist != NULL)
  {
    // Print error distribution
    size_t err_total = 0;
    for(i = 0; i < flist->errors_len; i++) err_total += flist->errors[i];
    printf("Errors: %zu (%.2f%%)\n", err_total, (100.0*err_total) / total_seq);
    for(i = 0; i < flist->errors_len; i++) printf(" %zu", flist->errors[i]);
    printf("\n");

    filelist_dealloc(flist);
  }

  return EXIT_SUCCESS;
}
