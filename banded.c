/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
#include"banded.h"
#define fabs(x) ((x)>0.0 ?(x):-(x))


/* n=number of equations [0..N] -> N+1, etc... */

void band2lu(double **a,int n,int m1,int m2,double **al,int *indx)
{
        int i,j,k,l;
        int mm;
        double temp;

        mm=m1+m2+1;
        l=m1;
        for (i=1;i<=m1;i++)
        {
                for (j=m1+2-i;j<=mm;j++)
                        a[i-1][j-l-1]=a[i-1][j-1];
                l--;
                for (j=mm-l;j<=mm;j++)
                        a[i-1][j-1]=0.0;
        }
        l=m1;
        for (k=1;k<=n;k++)
        {
                temp=a[k-1][0];
                i=k;
                if (l < n) l++;
                for (j=k+1;j<=l;j++)
                {
                        if (fabs(a[j-1][0]) > fabs(temp))
                        {
                                temp=a[j-1][0];
                                i=j;
                        }
                }
                indx[k-1]=i;
                if (temp == 0.0)
                        a[k-1][0]=TINY;
                if (i != k)
                {
                        for (j=1;j<=mm;j++)
                                swap(a[k-1][j-1],a[i-1][j-1]);
                }
                for (i=k+1;i<=l;i++)
                {
                        temp=a[i-1][0]/a[k-1][0];
                        al[k-1][i-k-1]=temp;
                        for (j=2;j<=mm;j++)
                                a[i-1][j-2]=a[i-1][j-1]-temp*a[k-1][j-1];
                        a[i-1][mm-1]=0.0;
                }
        }
}

void bandsolve(double **a,int n,int m1,int m2,double **al,int *indx,double *b)
{
        int i,k,l;
        int mm;
        double temp;

        mm=m1+m2+1;
        l=m1;
        for (k=1;k<=n;k++)
        {
                i=indx[k-1];
                if (i != k)
                        swap(b[k-1],b[i-1]);
                if (l < n)
                        l++;
                for (i=k+1;i<=l;i++)
                        b[i-1] -= al[k-1][i-k-1]*b[k-1];
        }
        l=1;
        for (i=n;i>=1;i--)
        {
                temp=b[i-1];
                for (k=2;k<=l;k++)
                        temp -= a[i-1][k-1]*b[k+i-2];
                b[i-1]=temp/a[i-1][0];
                if (l < mm)
                        l++;
        }
}

void bandmult(double **a,int n,int m1,int m2,double *x,double *b)
{
        int i,j,k,tmploop;

        for (i=1;i<=n;i++)
        {
                k=i-m1-1;
                tmploop=min(m1+m2+1,n-k);
                b[i-1]=0.0;
                for (j=max(1,1-k);j<=tmploop;j++)
                {
                        b[i-1] += a[i-1][j-1]*x[j+k-1];
                }
        }
}
/*-----------------------------------------------------------*/



/*---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
