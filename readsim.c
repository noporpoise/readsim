#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <inttypes.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <zlib.h>

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

// Sample a random number from normal distribution with [mean 0, variance 1]
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

void sum_qual_scores(long double *qprob_sum, size_t *counts, size_t len,
                     const read_t *read, int offset)
{
  size_t i, limit = len < read->qual.end ? len : read->qual.end;
  long double qual, qprob;
  if(read->qual.end != read->seq.end) die("Invalid qual");
  for(i = 0; i < limit; i++) {
    if(read->qual.b[i] < offset) die("Invalid qual");
    if(read->qual.b[i] - offset > 50) die("Invalid qual");
    if(read->seq.b[i] != 'N') {
      qual = read->qual.b[i] - offset;
      qprob = powl(10, -qual/10);
      qprob_sum[i] += qprob;
      counts[i]++;
    }
  }
}

size_t generate_qprob(const char *path,
                      long double **qprob_ptr, size_t **counts_ptr)
{
  int qoffset, minq, maxq, fmt;
  read_t r;

  seq_file_t *file = seq_open(path);
  seq_read_alloc(&r);

  // Get FASTQ offset
  if((fmt = seq_guess_fastq_format(file, &minq, &maxq)) == -1)
    die("Cannot get quality score from: %s", path);

  qoffset = FASTQ_OFFSET[fmt];

  if(seq_read(file, &r) <= 0)
    die("Empty profile file %s", path);

  size_t i, len = r.seq.end, num_reads;
  size_t *counts;
  long double *qprob;

  // qual_sum = calloc(len, sizeof(size_t));
  counts = calloc(len, sizeof(size_t));
  qprob = malloc(len * sizeof(long double));

  if(counts == NULL || qprob == NULL) die("Out of memory");

  for(i = 0; i < len; i++) qprob[i] = 0;

  sum_qual_scores(qprob, counts, len, &r, qoffset);

  for(num_reads = 1; seq_read(file, &r) > 0; num_reads++)
    sum_qual_scores(qprob, counts, len, &r, qoffset);

  seq_close(file);
  seq_read_dealloc(&r);

  printf("num profile reads: %zu\n", num_reads);

  long double qtotal = 0;
  for(i = 0; i < len; i++) {
    qprob[i] /= counts[i];
    qtotal += qprob[i];
  }

  // Normalise
  // long double qmean = qtotal / len, err_rate = 0.01;
  // for(i = 0; i < len; i++)
  //   qprob[i] = err_rate * qprob[i] / qmean;

  *qprob_ptr = qprob;
  *counts_ptr = counts;

  return len;
}

// given base c and random i (0 <= i <= 3)
char mut(char c, int i)
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

void add_seq_error(char *r, size_t rlen,
                   const long double *qprob, size_t *mcount, size_t errlen)
{
  size_t i, limit = rlen < errlen ? rlen : errlen;
  int rnd;
  for(i = 0; i < limit; i++) {
    if((rnd = rand()) < qprob[i] * RAND_MAX) {
      r[i] = mut(r[i], rnd % 3);
      mcount[i]++;
    }
  }
}

// Load reads from a file, apply sequence error, dump
void load_reads(const char *path, gzFile gzout,
                const long double *qprob, size_t *mcount, size_t len)
{
  seq_file_t *file;
  read_t r;
  int i;

  seq_read_alloc(&r);

  for(i = 0; i < 2; i++)
  {
    if((file = seq_open(path)) == NULL) die("Cannot read: %s", path);

    while(seq_read(file, &r) > 0) {
      add_seq_error(r.seq.b, r.seq.end, qprob, mcount, len);
      gzprintf(gzout, "@%s\n%s\n+\n%s\n", r.name.b, r.seq.b, r.qual.b);
    }

    seq_close(file);
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

void sim_reads(const char *refpath, gzFile out0, gzFile out1, long double *qprob,
               size_t *mcounts, size_t errlen, size_t tlen,
               size_t tlen_std_dev, size_t rlen, double depth)
{
  seq_file_t *file;
  size_t i, rcap, nchroms, glen = 0, nreads, chr, pos0, pos1;
  read_t *r;

  if((file = seq_open(refpath)) == NULL) die("Cannot read: %s", refpath);

  rcap = 16;
  r = malloc(rcap * sizeof(read_t));

  for(i = 0; ; i++) {
    if(i == rcap) {
      rcap *= 2;
      r = realloc(r, rcap * sizeof(read_t));
    }
    seq_read_alloc(&r[i]);
    if(seq_read(file, &r[i]) <= 0) break;
    glen += r[i].seq.end;
  }
  nchroms = i;

  if(nchroms == 0) die("No sequences in ref genome file: %s", refpath);

  // Sample
  nreads = (glen * depth) / (2 * rlen);
  char read0[rlen+1], read1[rlen+1];
  read0[rlen] = read1[rlen] = '\0';

  for(i = 0; i < nreads; i++)
  {
    chr = (nchroms == 1) ? 0 : rand_chrom(r, nchroms, glen);
    pos0 = rand() * ((double)(r[chr].seq.end-tlen) / RAND_MAX);
    pos1 = pos0 + tlen + ran_normal()*tlen_std_dev - rlen;
    if(pos1 + rlen > r[chr].seq.end) pos1 = r[chr].seq.end-rlen;
    memcpy(read0, r[chr].seq.b+pos0, rlen);
    memcpy(read1, r[chr].seq.b+pos1, rlen);
    add_seq_error(read0, rlen, qprob, mcounts, errlen);
    add_seq_error(read1, rlen, qprob, mcounts, errlen);
    gzprintf(out0, ">r%zu:%s:%zu:%zu\n.*s\n", i, (int)rlen, read0,
             r[chr].seq.b, pos0, pos1);
    gzprintf(out1, ">r%zu:%s:%zu:%zu\n.*s\n", i, (int)rlen, read1,
             r[chr].seq.b, pos0, pos1);
  }

  seq_close(file);

  for(i = 0; i < nchroms; i++) seq_read_dealloc(&r[i]);
  free(r);
}

int main(int argc, char **argv)
{
  if(argc < 4) print_usage(usage, NULL);

  // Sample reads from ref
  char *refpath = NULL;
  int tlen = 800, tlen_std_dev = 100, rlen = 250;
  double depth = 1;
  int optr = 0, optt = 0, optv = 0, optl = 0, optd = 0; // keeps track of values

  char *input0 = NULL, *input1 = NULL;

  int c;
  while((c = getopt(argc, argv, "r:t:l:d:1:2:")) >= 0) {
    switch (c) {
      case 'r': refpath = optarg; optr++; break;
      case 't': tlen = atoi(optarg); optt++; break;
      case 'v': tlen_std_dev = atoi(optarg); optv++; break;
      case 'l': rlen = atoi(optarg); optl++; break;
      case 'd': depth = atol(optarg); optd++; break;
      case '1': input0 = optarg; break;
      case '2': input1 = optarg; break;
    }
  }

  if(argc - optind < 3) print_usage(usage, NULL);

  if((optt > 0 || optv > 0 || optl > 0 || optd > 0) && refpath == NULL)
    print_usage(usage, "Missing -r <in.fa>");

  if(optr > 1 || optt > 1 || optv > 1 || optl > 1 || optd > 1)
    print_usage(usage, "Duplicate args");

  if((input0 == NULL) != (input1 == NULL))
    print_usage(usage, "Need both -1 <in> -2 <in>");

  char *out0path = argv[optind+1], *out1path = argv[optind+2];  
  gzFile out0, out1;

  if(input0 != NULL)
  {
    if(!test_file_readable(input0)) die("Cannot read: %s", input0);
    if(!test_file_readable(input1)) die("Cannot read: %s", input1);
  }

  if(refpath != NULL) {
    if((out0 = gzopen(out0path, "w")) == NULL) die("Cannot open: %s", out0path);
    if((out1 = gzopen(out1path, "w")) == NULL) die("Cannot open: %s", out1path);
  }

  // Generate profile from profile.fq

  // qprob[i] is the probability of making a seqn error at base i
  long double *qprob, sumprob;
  size_t i, errlen, *counts, sumcount;

  errlen = generate_qprob(argv[optind], &qprob, &counts);

  // Print qprob
  printf("QProb: %Lf", qprob[0]);
  sumprob = qprob[0];
  for(i = 1; i < errlen; i++) {
    printf(",%Lf", qprob[i]);
    sumprob += qprob[i];
  }
  printf("\nQProbSum: %Lf QProbMean: %Lf\n", sumprob, sumprob / errlen);

  printf("Counts: %zu", counts[0]);
  sumcount = counts[0];
  for(i = 1; i < errlen; i++) {
    printf(",%zu", counts[i]);
    sumcount += counts[i];
  }
  printf("\nCountsSum: %zu CountsMean: %zu\n", sumcount, sumcount / errlen);

  // Reuse mcounts for recording distribution of mutations
  memset(counts, 0, sizeof(*counts) * errlen);

  if(input0 != NULL)
  {
    load_reads(input0, out0, qprob, counts, errlen);
    load_reads(input1, out1, qprob, counts, errlen);
  }

  if(refpath != NULL) {
    sim_reads(refpath, out0, out1, qprob, counts, errlen,
              tlen, tlen_std_dev, rlen, depth);
  }

  if(out0 != NULL) {
    gzclose(out0);
    gzclose(out1);

    printf("SeqErr: %zu", counts[0]);
    sumcount = counts[0];
    for(i = 1; i < errlen; i++) {
      printf(",%zu", counts[i]);
      sumcount += counts[i];
    }
    printf("\nSeqErrSum: %zu SeqErrMean: %f\n", sumcount, (double)sumcount / errlen);
  }

  free(qprob);
  free(counts);

  return EXIT_SUCCESS;
}
