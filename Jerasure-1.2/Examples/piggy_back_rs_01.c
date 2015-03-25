#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "jerasure.h"
#include "reed_sol.h"
#include "piggyback_rs.h"
#include "time.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

usage(char *s)
{
  fprintf(stderr, "usage: piggy_back_01 k m w err size- Vandermonde matrices in GF(2^w).\n");
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

static void print_data_and_coding(int k, int m, int w, int size, 
				  char **data, char **coding) 
{
  int i, j, x, ind;
  int n, sp;
  long l;

  if(k > m) n = k;
  else n = m;
  sp = size * 2 + size/(w/8) + 8;

  printf("%-*sCoding\n", 2*sp-8, "Data");
  for(i = 0; i < n; i++) {
    if(i < k) {
      printf("D%-2d:", i);
      ind = i;
      while(ind <= i+k){
	for(j=0;j< size; j+=(w/8)) { 
	  printf(" ");
	  for(x=0;x < w/8;x++){
	    printf("%02x", (unsigned char)data[ind][j+x]);
	  }
	}
	ind = ind + k;
      }
      printf("    ");
    }
    else printf("%*s", sp, "");
    if(i < m) {
      printf("C%-2d:", i);
      ind = i;
      while(ind <= i+m){
	for(j=0;j< size; j+=(w/8)) { 
	  printf(" ");
	  for(x=0;x < w/8;x++){
	    printf("%02x", (unsigned char)coding[ind][j+x]);
	  }
	}
	ind = ind+m;
      }
    }
    printf("\n");
  }
  printf("\n");
}

int main(int argc, char **argv)
{
  long l;
  int k, w, i, j, m;
  int *matrix;
  char **data, **coding;
  int *erasures, *erased;
  int err;
  int size;
  
  if (argc != 6) usage(NULL);
  if (sscanf(argv[1], "%d", &k) == 0 || k <= 0) usage("Bad k");
  if (sscanf(argv[2], "%d", &m) == 0 || m <= 0) usage("Bad m");
  if (sscanf(argv[3], "%d", &w) == 0 || w <= 0 || w > 32) usage("Bad w");
  if (sscanf(argv[4], "%d", &err) == 0 ) usage("Bad err");
  if (sscanf(argv[5], "%d", &size) == 0 ) usage("Bad err");
  if (w <= 30 && k + m > (1 << w)) usage("k + m is too big");

  matrix = piggyback1_rs_vandermonde_coding_matrix(k,m,w);
  printf("Vandermonde Coding Matrix:\n\n");
  jerasure_print_matrix(matrix, 2*m, 2*k, w);
  printf("\n");

  srand48(time(NULL));
  //srand48(0);
  data = talloc(char *, 2*k);
  for (i = 0; i < k; i++) {
    data[i] = talloc(char, size);
    data[i+k] = talloc(char, size);
    for(j=0; j<size; j++){
      l = lrand48();
      memcpy(data[i]+j, &l, 1);
      l = lrand48();
      memcpy(data[i+k]+j, &l, 1);
    }
  }

  coding = talloc(char *, 2*m);
  for (i = 0; i < 2*m; i++) {
    coding[i] = talloc(char, size);
  }
  jerasure_matrix_encode(2*k, 2*m, w, matrix, data, coding, size);
  printf("Encoding Complete:\n\n");
  print_data_and_coding(k, m, w, size, data, coding);

  erasures = talloc(int, (2));
  erased = talloc(int, 2*(k+m));
  for (i = 0; i < 2*(m+k); i++) erased[i] = 0;
  l = 0;

  /*i = 0;
    while(i < 2*m){
    erasures[i] = lrand48()%k;
    if(erased[erasures[i]] == 0){
    erased[erasures[i]] = 1;
    erasures[i+1] = erasures[i]+k;
    erased[erasures[i+1]] = 1;
    i = i+2;
    }
    }

    l = 0;*/
  for (i = 0; i < 2*(err); ) {
    erasures[i] = lrand48()%k;
    if (erased[erasures[i]] == 0) {
      erased[erasures[i]] = 1;
      erasures[i+1] = erasures[i]+k;
      for(j = 0; j < size; j++){
	l = 0;
	memcpy(data[erasures[i]]+j, &l, 1);
	memcpy(data[erasures[i]+k]+j, &l, 1);
      }
      i=i+2;
    }
  }
  erasures[i] = -1;

  
  /*erased[erasures[0] + k + m] = 1;
    memcpy((erasures[0] < k) ? data[erasures[0]] : coding[erasures[0]-k], &l, sizeof(long));
    memcpy((erasures[0] < k) ? data[erasures[0]+k] : coding[erasures[0]-k+m], &l, sizeof(long));

    //erasures[2] = -1;*/

  printf("Erased D%d data device:\n\n",erasures[0]);
  print_data_and_coding(k, m, w, size, data, coding);
  
  piggyback1_rs_decode( k, m, w, matrix, erasures, data, coding, size);    

  // jerasure_matrix_decode(2*k, 2*m, w, matrix, 1, erasures, data, coding, sizeof(long));
  printf("After D%d data device recovered:\n\n",erasures[0]);
  print_data_and_coding(k, m, w, size, data, coding);
  return 0;
}
