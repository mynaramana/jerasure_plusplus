#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "galois.h"
#include "jerasure.h"
#include "reed_sol.h"
#include "piggyback_rs.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

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

int *piggyback1_rs_vandermonde_coding_matrix(int k, int m, int w)
{
  int tmp;
  int i, j, index;
  int *vdm, *dist;
  int q_s;
  int q_sl;
  vdm = reed_sol_vandermonde_coding_matrix(k, m, w);
  if (vdm == NULL) return NULL;
  dist = malloc(4*m*k*sizeof(int));
  if (dist == NULL) {
    free(vdm);
    return NULL;
  }
  memset(dist, 0, 4*m*k*sizeof(int));
  for( i = 0; i < k; i++){
    for (j = 0; j < m; j++) {
      dist[j*2*k+i] = vdm[j*k+i];
      dist[j*2*k+i+k+2*m*k] = vdm[j*k+i];
    }
  }
  if(m > 1){
    double q_sd = (double)k/(double)m + (double)(m-2)/(double)(2*m);
    q_s = (int)q_sd;
    if((double)q_s < q_sd) q_s++;//need to take upper bound here
    q_sl = k - (m-1)*q_s;
    //printf("t_1 = %d t_r = %d\n",q_s,q_sl);
  }
  for(j = 0; j < k; j++){
    dist[(m-1)*2*k +j+k ] = dist[(m-1)*2*k + j ];//((~dist[(m-1)*2*k + j ]) & 0xFF)+1;
    int row = j/q_s;
    int column = j%q_s;
    if( j < (m-1)*q_s)
      dist[2*(m+1)*k + row*2*k+ row*q_s+column] =  dist[(m-1)*2*k + j ];
    if(row == m-2)
      dist[(m-1)*2*k + j] = 0;
  }
  free(vdm);
  return dist;
}

int *piggyback2_rs_vandermonde_coding_matrix(int k, int m, int w)
{
  int tmp;
  int i, j, index;
  int *vdm, *dist;
  int q_s;
  int q_sl;
  vdm = reed_sol_vandermonde_coding_matrix(k, m, w);
  if (vdm == NULL) return NULL;
  dist = malloc(4*m*k*sizeof(int));
  if (dist == NULL) {
    free(vdm);
    return NULL;
  }
  memset(dist, 0, 4*m*k*sizeof(int));
  for( i = 0; i < k; i++){
    for (j = 0; j < m; j++) {
      dist[j*2*k+i] = vdm[j*k+i];
      dist[j*2*k+i+k+2*m*k] = vdm[j*k+i];
    }
  }
  if(m > 1){
    double q_sd = (double)k/(double)m + (double)(m-2)/(double)(2*m);
    q_s = (int)q_sd;
    if((double)q_s < q_sd) q_s++;//need to take upper bound here
    q_sl = k - (m-1)*q_s;
    //printf("t_1 = %d t_r = %d\n",q_s,q_sl);
  }
  for(j = 0; j < k; j++){
    dist[(m-1)*2*k +j+k ] = dist[(m-1)*2*k + j ];//((~dist[(m-1)*2*k + j ]) & 0xFF)+1;
    int row = j/q_s;
    int column = j%q_s;
    if( j < (m-1)*q_s)
      dist[2*(m+1)*k + row*2*k+ row*q_s+column] =  dist[(m-1)*2*k + j ];
    if(row == m-2)
      dist[(m-1)*2*k + j] = 0;
  }
  free(vdm);
  return dist;
}

void piggyback1_rs_decode( int k, int m, int w, int* matrix, int* erasures, char** data, char** coding, int size){
  char *temp, *var;
  int temp_;
  int i,j;
  double q_sd = (double)k/(double)m + (double)(m-2)/(double)(2*m);
  int q_s = (int)q_sd;
  if((double)q_s < q_sd) q_s++;//need to take upper bound here

  int set = erasures[0]/q_s;
  //printf("belongs to set %d\n",set);
  //print_data_and_coding(k, m, w, size, data, coding);
  /*first obtain b_e*/
  for(j = 0; j<k;j++){
    if(erasures[0] != j){
      galois_region_xor(data[erasures[0]+k],data[j+k],data[erasures[0]+k], size);
    }
  }
  galois_region_xor(data[erasures[0]+k],coding[m],data[erasures[0]+k], size);
  /*now obtain a_e*/
  if(set < m-1){
    jerasure_matrix_dotprod(k,w,&matrix[2*(m+1)*k+set*2*k+k],NULL, k+erasures[0], &data[k], data, size);
	
    for(j =0; j < k; j++){
      if(erasures[0] != j){		
        if(set == j/q_s){
		  
	  if(matrix[2*(m+1)*k+set*2*k+j] == 1){
	    galois_region_xor(data[erasures[0]],data[j],data[erasures[0]],size);
	  }else{
	    temp = talloc(char, size);
	    galois_w08_region_multiply( data[j],matrix[2*(m+1)*k+set*2*k+j],size,temp, 0);
	    galois_region_xor( temp, data[erasures[0]], data[erasures[0]], size);
	    free(temp);
	  }
	      
        }
      }
    }
    galois_region_xor(data[erasures[0]],coding[set+m+1],data[erasures[0]], size);
    for(i = 0; i < size; i++){
      temp_ = 0;
      memcpy(	&temp_, data[erasures[0]] + i,1);
      temp_ = galois_single_divide(temp_, matrix[2*(m+1)*k+set*2*k+erasures[0]], w);
      memcpy(data[erasures[0]]+i, &temp_, 1);
    }
    
  }else{
    var = talloc( char , size );
    memcpy(var, data[erasures[0]],size);
    for(j=1;j<m;j++){
      jerasure_matrix_dotprod(k,w,&matrix[2*m*k+j*2*k+k], NULL, k+erasures[0], &data[k], data, size);
      galois_region_xor( data[erasures[0]], var, var, size);
      if(j!=m-1)
	galois_region_xor( coding[j+m], var, var, size);
    }
    galois_region_xor( coding[m-1], var, var, size);
    for(j = set*q_s; j < k ; j++){
      if(j!=erasures[0]){
	temp = talloc(char, size);
	galois_w08_region_multiply( data[j],matrix[2*(m-1)*k+j],size,temp, 0);
	galois_region_xor( temp, var, var, size);
	free(temp);
      }
    }
    memcpy(data[erasures[0]], var, size);
    free(var);

    long val_ = 0;
    for(i = 0; i < size; i++){
      temp_ = 0;
      memcpy(	&temp_, data[erasures[0]] + i,1);
      temp_ = galois_single_divide(temp_, matrix[2*(m-1)*k+erasures[0]], w);
      val_ = val_ + (temp_ << (i*8));
    }
    memcpy(data[erasures[0]], &val_,size);
  }
	
}

void ReadDataFromFile(int ind, int DataBit, int FirstHalf, int size, char* data_ptr, char* curdir, char* cs1, char* cs2, int k, int n){
  char fname[100];
  char temp_str[10];
  int md;  
  FILE* fp;
	
  sprintf(temp_str, "%d", k);
  md = strlen(temp_str);

  if(DataBit){
    sprintf(fname, "%s/Coding/%s_k%0*d%s", curdir, cs1, md, ind, cs2);
  }else{
    sprintf(fname, "%s/Coding/%s_m%0*d%s", curdir, cs1, md, ind, cs2);
  }

  fp = fopen(fname, "rb");
  if(FirstHalf){
    fseek(fp, (n-1)*size*2, SEEK_SET);			
  }else{
    fseek(fp, (n-1)*size*2+size, SEEK_SET);			
  }
		
  fread(data_ptr, sizeof(char), size, fp);
  fclose(fp);

}
void piggyback1_rs_repair( int k, int m, int w, int* matrix, int* erasures, char** data, char** coding, int size, char* curdir, char* cs1, char* cs2, int n){
  char *temp, *var;
  int temp_;
  int i,j;
  double q_sd = (double)k/(double)m + (double)(m-2)/(double)(2*m);
  int q_s = (int)q_sd;
  FILE* fp;
  int md;
  char temp_str[10];
  char fname[100];
  
  sprintf(temp_str, "%d", k);
  md = strlen(temp_str);
  memset(data[erasures[0]], 0, size);
  memset(data[erasures[0]+k], 0, size);
	
  if((double)q_s < q_sd) q_s++;//need to take upper bound here

  int set = erasures[0]/q_s;
  //printf("belongs to set %d\n",set);
  //print_data_and_coding(k, m, w, size, data, coding);
  /*first obtain b_e*/
  for(j = 0; j<k;j++){
    if(erasures[0] != j){
      ReadDataFromFile(j+1,1, 0, size, data[j+k], curdir, cs1, cs2, k, n);
      galois_region_xor(data[erasures[0]+k],data[j+k],data[erasures[0]+k], size);
    }
  }
  ReadDataFromFile(1,0, 0, size, coding[m], curdir, cs1, cs2, k, n);
  galois_region_xor(data[erasures[0]+k],coding[m],data[erasures[0]+k], size);
  /*now obtain a_e*/
  if(set < m-1){
    jerasure_matrix_dotprod(k,w,&matrix[2*(m+1)*k+set*2*k+k],NULL, k+erasures[0], &data[k], data, size);
	
    for(j =0; j < k; j++){
      if(erasures[0] != j){		
        if(set == j/q_s){
	  ReadDataFromFile(j+1,1, 1, size, data[j], curdir, cs1, cs2, k, n);
	  if(matrix[2*(m+1)*k+set*2*k+j] == 1){
	    galois_region_xor(data[erasures[0]],data[j],data[erasures[0]],size);
	  }else{
	    temp = talloc(char, size);
	    galois_w08_region_multiply( data[j],matrix[2*(m+1)*k+set*2*k+j],size,temp, 0);
	    galois_region_xor( temp, data[erasures[0]], data[erasures[0]], size);
	    free(temp);
	  }
	      
        }
      }
    }
    ReadDataFromFile(set+2,0, 0, size, coding[set+m+1], curdir, cs1, cs2, k, n);
    sprintf(fname, "%s/Coding/%s_m%0*d%s", curdir, cs1, md, set+2, cs2);
    galois_region_xor(data[erasures[0]],coding[set+m+1],data[erasures[0]], size);
    for(i = 0; i < size; i++){
      temp_ = 0;
      memcpy(	&temp_, data[erasures[0]] + i,1);
      temp_ = galois_single_divide(temp_, matrix[2*(m+1)*k+set*2*k+erasures[0]], w);
      memcpy(data[erasures[0]]+i, &temp_, 1);
    }
    
  }else{
    var = talloc( char , size );
    memcpy(var, data[erasures[0]],size);
    for(j=1;j<m;j++){
      for(i = 0; i < k; i++){
	ReadDataFromFile(i+1, 1, 0, size, data[k+i], curdir, cs1, cs2, k, n);
      }
      jerasure_matrix_dotprod(k,w,&matrix[2*m*k+j*2*k+k], NULL, k+erasures[0], &data[k], data, size);
      galois_region_xor( data[erasures[0]], var, var, size);
      if(j!=m-1){
	ReadDataFromFile( j+1, 0, 0, size, coding[j+m], curdir, cs1, cs2, k, n);
	galois_region_xor( coding[j+m], var, var, size);
      }
    }
    ReadDataFromFile( j+1, 0, 1, size, coding[m-1], curdir, cs1, cs2, k, n);
    galois_region_xor( coding[m-1], var, var, size);
    for(j = set*q_s; j < k ; j++){
      if(j!=erasures[0]){
	temp = talloc(char, size);
	ReadDataFromFile( j+1, 1, 0, size, data[j], curdir, cs1, cs2, k, n);
	galois_w08_region_multiply( data[j],matrix[2*(m-1)*k+j],size,temp, 0);
	galois_region_xor( temp, var, var, size);
	free(temp);
      }
    }
    memcpy(data[erasures[0]], var, size);
    free(var);

    long val_ = 0;
    for(i = 0; i < size; i++){
      temp_ = 0;
      memcpy(	&temp_, data[erasures[0]] + i,1);
      temp_ = galois_single_divide(temp_, matrix[2*(m-1)*k+erasures[0]], w);
      val_ = val_ + (temp_ << (i*8));
    }
    memcpy(data[erasures[0]], &val_,size);
  }
	
}

