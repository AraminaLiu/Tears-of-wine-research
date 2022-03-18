/*-------------------------------------MEMORY ALLOCATION----------*/
#include<stdio.h>
#include<stdlib.h>
#define myfree(p) {if(p){free(p);p=NULL;}}
double *myalloc(int);
int *myalloci(int);
double **matrix(int,int); 
/*-------------------------------(Dynamic Memory allocation stuff)------*/
double *myalloc(int n)
{
    double *t;
        if((t=(double *)calloc(n,sizeof(double)))==NULL)
        {fprintf(stderr,"Memory allocation error\n");exit(1);}
        return(t);
}
int *myalloci(int n)
{
        int *t;
        if((t=(int *)calloc(n,sizeof(int)))==NULL)
        {fprintf(stderr,"Memory allocation error\n");exit(1);}
        return(t);
}
 
/*--------------------------------------------------------------------------*/
/*-----------------------------------------------------------*/
/*-----------------------------------------------------------*/
double **matrix(int n,int m) /* double **a=matrix[0..n-1][0..m-1] */
{
        int i;
        double **matrx;
        matrx=(double **)calloc(n,sizeof(double*));
        if(matrx==NULL)
        {fprintf(stderr,"ERROR!\n");exit(1);}
        matrx[0]=myalloc(n*m);/* flat memory model */
        for(i=1;i<n;++i)
                matrx[i]=matrx[0]+i*m;
        return(matrx); /* return memory for matrix */
}
/*-------------------------------(Dynamic Memory allocation stuff)------*/
/*--------------------------------------------------------------------------*/
 


