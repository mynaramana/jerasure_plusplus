extern int *lrc_rs_vandermonde_coding_matrix(int k,int m,int w,int l);
extern int lrc_rs_vandermode_decode(int k,int m,int w,int l,int matrix,int * erasures,char * * data,char * * coding,int size);
extern int lrc_rs_vandermode_repair(int k,int m,int w,int l,int matrix,int * erasures,char * * data,char * * coding,int size);
