#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "galois.h"
#include "jerasure.h"
#include "reed_sol.h"
#include "piggyback_rs.h"


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

int *lrc_rs_vandermonde_coding_matrix(int k, int m, int w, int l)
{
  int tmp;
  int *vdm, *dist;
  int i;
  vdm = reed_sol_vandermonde_coding_matrix(k, m, w);
  if (vdm == NULL) return NULL;
  dist = malloc((m+l)*k*sizeof(int));
  if (dist == NULL) {
    free(vdm);
    return NULL;
  }
  memset(dist+m*k, 0, l*k*sizeof(int));
  memcpy(dist, vdm, m*k*sizeof(int));
  int cur_set = 0;
  int cur_set_cnt = k/l;

  int prev_sets_cnt = 0;
  if((double)cur_set_cnt < ((double)k/(double)l)) cur_set_cnt++;
  while(cur_set_cnt){
    //printf("cur_Set_cnt = %d cur_set = %d\n",cur_set_cnt,cur_set);
    for(i=0;i<cur_set_cnt;i++){
      dist[m*k+cur_set*k+i+prev_sets_cnt] = 1;
    }

    prev_sets_cnt += cur_set_cnt;
    if( l - cur_set == 2)
      cur_set_cnt = k - prev_sets_cnt;
    else if( l- cur_set == 1)break;
    else{
      cur_set_cnt = (k-prev_sets_cnt)/(l- cur_set -1);
      if((double)cur_set_cnt < ((double)(k - prev_sets_cnt)/(double)(l-cur_set-1))) cur_set_cnt++;
    }
    cur_set++;
  }


  free(vdm);
  return dist;
}

int lrc_rs_vandermode_decode( int k, int m, int w, int l, int matrix, int* erasures, char** data, char** coding, int size){
  int* erased = jerasure_erasures_to_erased(k, m+l, erasures);

  int cur_set = 0;
  int cur_set_cnt = k/l;
  int prev_sets_cnt = 0;
  int i;
  int curset_eras_cnt;
  int erased_node;

  //print_data_and_coding(k, m, l,w,size,data,coding);
  if((double)cur_set_cnt < ((double)k/(double)l)) cur_set_cnt++;
  while(cur_set_cnt){
    curset_eras_cnt = 0;
    erased_node = -1;
    //printf("cur_set_cnt = %d cur_set = %d\n",cur_set_cnt,cur_set);
    /*check among the sets if any node is locally repairable*/
    for(i=0;i<cur_set_cnt;i++){
      if(erased[prev_sets_cnt+i] == 1){
	curset_eras_cnt++;
	erased_node = i;
	if(curset_eras_cnt > 1)break;
      }
    }
    //printf("cur_set_cnt = %d cur_set = %d erased_node = %d\n",cur_set_cnt,cur_set,erased_node);
    if(curset_eras_cnt == 1){/*perform local repair*/
      galois_region_xor( data[prev_sets_cnt+erased_node], coding[m+cur_set], data[prev_sets_cnt+erased_node],size);
      for(i=0;i<cur_set_cnt;i++){
	if(i!=erased_node)
	  galois_region_xor( data[prev_sets_cnt+erased_node], data[prev_sets_cnt+i], data[prev_sets_cnt+erased_node],size);
      }
      erased[erased_node+prev_sets_cnt] = 0;
    }
    prev_sets_cnt += cur_set_cnt;
    if( l - cur_set == 2)
      cur_set_cnt = k - prev_sets_cnt;
    else if( l- cur_set == 1)break;
    else{
      cur_set_cnt = (k-prev_sets_cnt)/(l- cur_set -1);
      if((double)cur_set_cnt < ((double)(k - prev_sets_cnt)/(double)(l-cur_set-1))) cur_set_cnt++;
    }
    cur_set++;
  }
  
}

int lrc_rs_vandermode_repair( int k, int m, int w, int l, int matrix, int* erasures, char** data, char** coding, int size, char* curdir, char* cs1, char* cs2, int n){
  int* erased = jerasure_erasures_to_erased(k, m+l, erasures);

  int cur_set = 0;
  int cur_set_cnt = k/l;
  int prev_sets_cnt = 0;
  int i;
  int curset_eras_cnt;
  int erased_node;

  //print_data_and_coding(k, m, l,w,size,data,coding);
  if((double)cur_set_cnt < ((double)k/(double)l)) cur_set_cnt++;
  while(cur_set_cnt){
    curset_eras_cnt = 0;
    erased_node = -1;
    //printf("cur_set_cnt = %d cur_set = %d\n",cur_set_cnt,cur_set);
    /*check among the sets if any node is locally repairable*/
    for(i=0;i<cur_set_cnt;i++){
      if(erased[prev_sets_cnt+i] == 1){
	curset_eras_cnt++;
	erased_node = i;
	if(curset_eras_cnt > 1)break;
      }
    }
    //printf("cur_set_cnt = %d cur_set = %d erased_node = %d\n",cur_set_cnt,cur_set,erased_node);
    if(curset_eras_cnt == 1){/*perform local repair*/
      galois_region_xor( data[prev_sets_cnt+erased_node], coding[m+cur_set], data[prev_sets_cnt+erased_node],size);
      for(i=0;i<cur_set_cnt;i++){
	if(i!=erased_node)
	  galois_region_xor( data[prev_sets_cnt+erased_node], data[prev_sets_cnt+i], data[prev_sets_cnt+erased_node],size);
      }
      erased[erased_node+prev_sets_cnt] = 0;
    }
    prev_sets_cnt += cur_set_cnt;
    if( l - cur_set == 2)
      cur_set_cnt = k - prev_sets_cnt;
    else if( l- cur_set == 1)break;
    else{
      cur_set_cnt = (k-prev_sets_cnt)/(l- cur_set -1);
      if((double)cur_set_cnt < ((double)(k - prev_sets_cnt)/(double)(l-cur_set-1))) cur_set_cnt++;
    }
    cur_set++;
  }
  
}

