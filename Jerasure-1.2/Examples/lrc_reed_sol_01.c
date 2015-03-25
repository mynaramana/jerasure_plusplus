#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "jerasure.h"
#include "reed_sol.h"
#include "lrc_reed_sol.h"
#include "time.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

usage(char *s)
{
  fprintf(stderr, "usage: lrc_reed_sol_01 k m w l size- Vandermonde matrices in GF(2^w).\n");
  fprintf(stderr, "       \n");
  fprintf(stderr, "       k+m must be <= 2^w.  This simply prints out the \n");
  fprintf(stderr, "       Vandermonde matrix in GF(2^w), and then the distribution\n");
  fprintf(stderr, "       matrix that is constructed from it.  See [Plank-Ding-05] for\n");
  fprintf(stderr, "       information on how this construction proceeds\n");
  fprintf(stderr, "       \n");
  fprintf(stderr, "This demonstrates: reed_sol_extended_vandermonde_matrix()\n");
  fprintf(stderr, "                   reed_sol_big_vandermonde_coding_matrix()\n");
  fprintf(stderr, "                   reed_sol_vandermonde_coding_matrix()\n");
  fprintf(stderr, "                   jerasure_print_matrix()\n");
  if (s != NULL) fprintf(stderr, "%s\n", s);
  exit(1);
}

static void print_data_and_coding(int k, int m, int loc, int w, int size, 
				  char **data, char **coding) 
{
  int i, j, x;
  int n, sp;
  long l;

  if(k > m) n = k;
  else n = m;
  sp = size * 2 + size/(w/8) + 8;

  printf("%-*s%-*sLocalParity\n", sp, "Data", sp, "Codes");
  for(i = 0; i < n; i++) {
    if(i < k) {
      printf("D%-2d:", i);
      for(j=0;j< size; j+=(w/8)) { 
	printf(" ");
	for(x=0;x < w/8;x++){
	  printf("%02x", (unsigned char)data[i][j+x]);
	}
      }
      printf("    ");
    }
    else printf("%*s", sp, "");
    if(i < m) {
      printf("C%-2d:", i);
      for(j=0;j< size; j+=(w/8)) { 
	printf(" ");
	for(x=0;x < w/8;x++){
	  printf("%02x", (unsigned char)coding[i][j+x]);
	}
      }
      printf("    ");
    }
    if(i < loc) {
      printf("L%-2d:", i);
      for(j=0;j< size; j+=(w/8)) { 
	printf(" ");
	for(x=0;x < w/8;x++){
	  printf("%02x", (unsigned char)coding[i+m][j+x]);
	}
      }
    }
    printf("\n");
  }
  printf("\n");
}

int main(int argc, char **argv)
{
  long l;
  int k, w, i, j, m, loc;
  int *matrix;
  char **data, **coding;
  int *erasures, *erased;
  int temp_;
  char *temp, *var, *temp1;
  int size;
  if (argc != 6) usage(NULL);
  if (sscanf(argv[1], "%d", &k) == 0 || k <= 0) usage("Bad k");
  if (sscanf(argv[2], "%d", &m) == 0 || m <= 0) usage("Bad m");
  if (sscanf(argv[3], "%d", &w) == 0 || w <= 0 || w > 32) usage("Bad w");
  if (sscanf(argv[4], "%d", &loc) == 0 || loc > k || w < 0) usage("Bad loc");
  if (sscanf(argv[5], "%d", &size) == 0 || size < 0) usage("Bad size");
  if (w <= 30 && k + m > (1 << w)) usage("k + m is too big");

  matrix = lrc_rs_vandermonde_coding_matrix(k,m,w,loc);
  printf("Vandermonde Coding Matrix:\n\n");
  jerasure_print_matrix(matrix, m+loc, k, w);
  printf("\n");

  srand48(time(NULL));
  //srand48(0);
  data = talloc(char *, k);
  for (i = 0; i < k; i++) {
    data[i] = talloc(char, size);
    data[i+k] = talloc(char, size);
    for(j=0;j<size;j++){
      l = lrand48();
      memcpy(data[i]+j, &l, 1);
    }
  }

  coding = talloc(char *, m+loc);
  for (i = 0; i < m+loc; i++) {
    coding[i] = talloc(char, sizeof(long));
  }
  jerasure_matrix_encode(k, m+loc, w, matrix, data, coding, size);
  printf("Encoding Complete:\n\n");
  print_data_and_coding(k, m, loc, w, size, data, coding);
  
  int err = 1;
  
  erasures = talloc(int, 2);
  erased = talloc(int, k+m);
  for (i = 0; i < m+k; i++) erased[i] = 0;
  
  for (i = 0; i < err; ) {
    erasures[i] = lrand48()%k;
    if (erased[erasures[i]] == 0) {
      erased[erasures[i]] = 1;
      l = 0;
      for(j=0;j<size;j++){
	memcpy(data[erasures[i]]+j, &l, 1);
      }
      i++;
    }
  }
  erasures[i] = -1;
  printf("Erased Node %d\n",erasures[0]);
  print_data_and_coding( k, m, loc, w, size, data, coding);
  lrc_rs_vandermode_decode( k, m, w, loc, matrix, erasures, data, coding, size);
  print_data_and_coding( k, m, loc, w, size, data, coding);
}
