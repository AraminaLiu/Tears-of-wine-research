/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/* swiped from numerical recipes (2nd) bandec.c, banbks.c */
#define swap(a,b) {double temp=(a);(a)=(b);(b)=temp;}
#define TINY 1.0e-20
#define max(a,b) ((a)>(b)? (a):(b))
#define min(a,b) ((b)>(a)? (a):(b))
void band2lu(double **a,int n,int m1,int m2,double **al,int *indx);
void bandsolve(double **a,int n,int m1,int m2,double **al,int *indx,double *b);
void bandmult(double **a,int n,int m1,int m2,double *x,double *b);

/*-------------------------------(Dynamic Memory allocation stuff)------*/
