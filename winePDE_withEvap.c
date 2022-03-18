/*=====================================*/
/* PDE solver for tears of wine project with surface tension gradient, 
normal and tangential gravity, fourth order surface tension term, and 
curvature contributions on a conical substrate


h_t + (h^2-h^3-A*h^3/(L0+x)^2)_x = -(h^3 h_xxx)_x + D*(h^3h_x)_x + beta/(u + K)*(PP(u) - u_xx)

where PP(u) = eps^2/h^3 - eps^3/h^4 - p0


include evaporation effects in the moving film


This version uses meniscus boundary conditions

h -> -x/D for x -> -\infty (or dh/dx = -1/D)
h -> b for x -> infty

in numerics, let's set

h'(0) = -1/D; h''(0) = 0
h(L) = hR; h'(L) = 0;

Upwind scheme for the (h^2-h^3-A*h^3/(L0+x)^2+B*h^4/(L0+x)^3)_x term

Centered finite difference for the [h^3(D*h-h_xx)_x]_x term

compile with

cc winePDE_withEvap.c -lm -O3


2020.02.09

*/

/*=====================================*/
/*=====================================*/
#include<stdio.h>
#include<math.h>
#include"myalloc.h"
/*--------------------------------------------------------*/
#define sqr(x) ((x)*(x))
#define sgn(x) ((x)>=0.0 ? 1.0 : -1.0)
#define pi M_PI
/*------------------------------------------------------------------------*/
#define swap(a,b) {double temp=(a);(a)=(b);(b)=temp;}
#define TINY 1.0e-20
#include "banded.h"
//#include "banded.c"
/*------------------------------------------------------------------------*/

int N;
double L,dx,dt;
double err_tol;
double tend;
/*---------------------------*/



double *U;
double *U0;
double *U0_record;
double *U0_init;



double *rhs;
double **J;
double **J1;
int *indx;


/*-------------------------*/
double param_D;
double param_A;
double param_B;

double param_beta;
double param_K;
double param_p0;
double eps;

double L0;


// Left and Right Dirichlet boundary conditions
double uL;
double uR; 

/*-------------------------*/


void build_rhs_jac(double *u,double *u_old);
void memory_setup();
void run();
int timestep(double tt); 

FILE *out;
FILE *out1;




double tt;
double tnext,dtnext;

int BAD;

/*-------------------------*/

int main(int argc,char *argv[])
{
	int i,j;


	char text[100];
	double xx;
	double HL, HR; // left and right boundary conditions


/*-------------------------*/
 	out=fopen("out.dat","w");
	out1=fopen("out1.dat","w");

/*--------------------*/



	 N = 6000;
	 
	 L = 100.0;


     
//      L0 = 5.0;
//    	param_A = 3.443924815;  // with curvature effects	
  	param_A = 0.0; // without curvature effects
    
    /* Vuilleumier data: Experiment CI */
    // param_D = 0.43;
    // HR = 0.0353; // dimensional b = 2 micrometers (compressive - undercompressive)
	//param_D = 0.353;
	//HR = 0.0317;
	param_D = 0.639;
	HR = 0.031;
 
    //param_beta = 0.01; // seems to work
	//param_beta = 0.02;
    //param_beta = 0.001;
	//param_beta = 0.03;
	param_beta = 0.025;
	
    
    
    param_K = 0.1;
    param_p0 = -0.1;
    
    
    eps = HR;
 
/*--------------------*/

/*========================================================================*/
	dx=L/N;
	memory_setup();

/*====================== set ICs ================*/

    double shift;



    shift = 20;


    for(i=0;i<N+1;++i)
    {
        xx = i*dx-shift;
         
         /* IC from Evans & Munch Meniscus paper, sensitive to pairing with h = HR right-side boundary condition */
         
         if(xx < 0.0)
            U0_init[i] = pow(param_D,-3.0/2.0)*(exp(pow(param_D,0.5)*(xx))-pow(param_D,0.5)*(xx)-1)+HR;
        else    
            U0_init[i] = HR;
 
    }

    
    uR = HR;

	dt=1e-3;

 	tt = 0.0;
 	
    tnext= 0.0; 	
 	dtnext = 20.0;
//  	dtnext = 0.01;
	
	tend = 1000;
	
	printf("param_D = %g\t b=%g\t A=%g\n", param_D, HR, param_A);
	
	printf("|Hx(0)- (-1.0/D)| ~ %g\n", fabs((U0_init[2]-U0_init[0])/(2.0*dx)-(-1.0/param_D)));
 	
	run();

 

}

/*******************************************************************/


void run()
{

	int i,b_count=0;

	double umin,umin0,umax;
	int imin,imax; // umin and umax are global min/maximum; 
	
	double mass,ht,p;
	double func; //energy functional
	double flux;


	umin0 = 0.0;


	/*====================== Initial profile ================*/
	for(i=0;i<N+1;++i)
	{
		U0[i] = U0_init[i];	
		 fprintf(out,"%g\t %g\n",
			i*dx,U0[i]);
	}		
	 fprintf(out,"\n\n");
	 fflush(out);
	 


	/*====================== Evolution starts ================*/

	while(tt<tend)
	{
		BAD=0;
		if(timestep(tt))
		{
			tt+=dt;
			
			dt*=1.1;
			

 			dt = min(0.01,dt);

			b_count=0; /* reset badstep counter */

		}
		else
		{
			dt*=0.5;
			BAD=1;
			b_count++;
// 			printf("BAD %d\n", b_count);
			if(b_count>8) 
				exit(1);
			continue; /* restart loop from top */
		}
	
		if(tt>tnext)
		{	
			//umin0=umin;
			tnext+=dtnext;
			mass=0.0;
			for(i=0;i<N;++i)
				mass+=0.5*dx*(U[i]+U[i+1]);

		
			imax=0;
			umax=U[0];
			imin = 0;
			umin = U[0];
			for(i=0;i<N+1;++i)
			{
				if(U[i]>umax)
				{
					imax=i;
					umax=U[i];
				}
				
				if(U[i]<umin)
				{
					imin = i;
					umin = U[i];
				}
			}


	
				
			func = 0.0;
			for(i=0;i<N;++i)
			{
				ht = (U0_record[i]-U[i])/dt;
	

				fprintf(out,"%g\t %g\t %g\t %g\t %g\t %g %g\n",
				i*dx,U[i],(U[i+1]-U[i])/(dx),imax*dx,tt,umax,ht);
 
			}
		
			

			fprintf(out,"\n\n");
			fflush(out);
			
			fprintf(out1,"%g\t %g\t %g\t %g\t %g\t %g\t %g\n",
			tt,umin, imin*dx, umax,imax*dx,mass, param_D);
			fflush(out1);
					
			printf("%g \t %g\t %g \t %g \t %g\t %g\n",tt,dt,mass,umin, umax, param_D);

		}
	
	

	}

		fprintf(out1,"\n\n");
		fflush(out1);
		

}


int timestep(double tt) /* return 1=good, 0=bad */
{
    int i,j;
	double maxerr;
	double relax;

	for(i=0;i<N+1;++i)
		U[i]=U0[i];

	err_tol=1e-9;
	
	err_tol = 1e-7;
	
// 	err_tol = 1e-6;
	

	for(j=0;j<10;++j)
	{
		maxerr=0.0;

		build_rhs_jac(U,U0);
		

		band2lu(J,N+1,2,2,J1,indx);
		bandsolve(J,N+1,2,2,J1,indx,rhs);
		
		relax = 1.0;
		
// 		if (j < 2)
// 		    relax = 0.2;


		for(i=0;i<N+1;++i)
		{
			maxerr=max(maxerr,fabs(rhs[i]));
			U[i]+=relax*rhs[i];
//  			U[i]=fabs(U[i]); /* prevent numerical instability*/
		}
		

		if(maxerr<err_tol/10.0 && j >2)
			break;

	}
//   	printf("j=%d\t maxerr=%g dt=%g\n", j, maxerr,dt);
	if(maxerr<err_tol) 
	{
		for(i=0;i<N+1;++i)
		{
			U0_record[i] = U0[i];
			U0[i]=U[i];  /* store previous step */
		}
		return(1);
	}
	else
		return(0);

}




/*--------------------------------------------------------------------------*/

void build_rhs_jac(double *u,double *u_old)
{
	int i,j,ni;

	double u0,up,upp,um,umm,uold;
	double x0;

	#include"tempvars.h"

	for(i=0;i<N+1;++i) /* zero to start with */
	{	
		rhs[i]=0.0;
		for(j=0;j<5;++j)
			J[i][j]=0.0;
	}

	for(i=1;i<N;++i)
	{
		/*-------------------------*/
		/*-------------------------*/
		ni=N-i;

		u0=u[i];
		up=u[i+1];
		upp=u[i+2];
		um=u[i-1];
		umm=u[i-2];
		uold=u_old[i];
		
		x0 = i*dx;

		switch(i) /* ghost points for left symmetry BCs */
		{
			case 0:
					umm=upp + 1.0/param_D*4.0*dx;
					um=up + 1.0/param_D*2.0*dx;
		break;
			case 1:
					umm=u0 + 1.0/param_D*2.0*dx;
		break;
		}
		switch(ni) /* ghost points for right symmetry BCs */
		{
			case 0:
					upp=umm;
					up=um;
		break;
			case 1:
					upp=u0;
		break;
		}

		if(2.0*u0-3.0*u0*u0
		-3*param_A*u0*u0/pow((L0+x0),2.0)>0)
		{
			/*-------------------------*/
	 t5 = u0*u0;
      t8 = L0+x0;
      t9 = t8*t8;
      t14 = u0-um;
      t16 = 1/dx;
      t18 = u0+up;
      t21 = t18*t18*t18/8.0;
      t25 = u0+um;
      t28 = t25*t25*t25/8.0;
      t35 = 3.0*u0;
      t38 = dx*dx;
      t40 = 1/t38/dx;
      t48 = t5*u0;
      t54 = eps*eps;
      t58 = t5*t5;
      t70 = (u0-uold)/dt+(2.0*u0-3.0*t5-3.0*param_A*t5/t9)*t14*t16-param_D*(t21
*(up-u0)*t16-t28*t14*t16)*t16+(t21*(upp-3.0*up+t35-um)*t40-t28*(up-t35+3.0*um-
umm)*t40)*t16+2.0*t48*param_A/t9/t8+param_beta*(t54/t48-t54*eps/t58-param_p0-(
up-0.2E1*u0+um)/t38)/(u0+param_K);

			rhs[i]= t70;

			/*-------------------------*/


			/*--------------*/
		t1 = u0+um;
      t5 = dx*dx;
      t6 = t5*t5;
      t8 = t1*t1*t1/t6/8.0;


			J[i][0]=t8;

			/*--------------*/
	t2 = u0*u0;
      t6 = pow(L0+x0,2.0);
      t11 = 1/dx;
      t13 = u0+um;
      t14 = t13*t13/4.0;
      t20 = t13*t14/2.0;
      t25 = u0+up;
      t29 = dx*dx;
      t31 = 1/t29/dx;
      t48 = -(2.0*u0-3.0*t2-3.0*param_A*t2/t6)*t11-param_D*(-3.0/2.0*t14*(u0-um
)*t11+t20*t11)*t11+(-t25*t25*t25*t31/8.0-3.0/2.0*t14*(up-3.0*u0+3.0*um-umm)*t31
-3.0*t20*t31)*t11-param_beta/t29/(u0+param_K);

			J[i][1]=t48;

			/*--------------*/
	t4 = L0+x0;
      t5 = t4*t4;
      t6 = 1/t5;
      t10 = u0-um;
      t12 = 1/dx;
      t15 = u0*u0;
      t17 = param_A*t15;
      t22 = u0+up;
      t23 = t22*t22/4.0;
      t29 = t23*t22/2.0;
      t31 = u0+um;
      t32 = t31*t31/4.0;
      t37 = t32*t31/2.0;
      t43 = 3.0*u0;
      t46 = dx*dx;
      t48 = 1/t46/dx;
      t66 = eps*eps;
      t67 = t15*t15;
      t68 = 1/t67;
      t71 = t66*eps;
      t76 = 1/t46;
      t80 = u0+param_K;
      t92 = t80*t80;
      t95 = 1/dt+(-6.0*param_A*u0*t6-6.0*u0+2.0)*t10*t12+(-3.0*t17*t6-3.0*t15+
2.0*u0)*t12-param_D*(3.0/2.0*t23*(up-u0)*t12-t29*t12-3.0/2.0*t32*t10*t12-t37*
t12)*t12+(3.0/2.0*t23*(upp-3.0*up+t43-um)*t48+3.0*t29*t48-3.0/2.0*t32*(up-t43+
3.0*um-umm)*t48+3.0*t37*t48)*t12+6.0*t17/t5/t4+param_beta*(-3.0*t66*t68+4.0*t71
/t67/u0+0.2E1*t76)/t80-param_beta*(t66/t15/u0-t71*t68-param_p0-(up-0.2E1*u0+um)
*t76)/t92;


			J[i][2]=t95;

			/*--------------*/
	 t1 = u0+up;
      t2 = t1*t1/4.0;
      t5 = 1/dx;
      t9 = t1*t2/2.0;
      t18 = dx*dx;
      t20 = 1/t18/dx;
      t25 = u0+um;
      t37 = -param_D*(3.0/2.0*t5*t2*(up-u0)+t9*t5)*t5+(3.0/2.0*t2*(upp-3.0*up+
3.0*u0-um)*t20-3.0*t9*t20-t25*t25*t25*t20/8.0)*t5-param_beta/t18/(u0+param_K);


			J[i][3]=t37;

			/*--------------*/
	t1 = u0+up;
      t5 = dx*dx;
      t6 = t5*t5;
      t8 = t1*t1*t1/t6/8.0;

			J[i][4]=t8;

			/*-------------------------*/
			/*-------------------------*/
		}
		else
		{
			/*-------------------------*/
	t5 = u0*u0;
      t8 = L0+x0;
      t9 = t8*t8;
      t14 = up-u0;
      t16 = 1/dx;
      t18 = u0+up;
      t21 = t18*t18*t18/8.0;
      t24 = u0+um;
      t27 = t24*t24*t24/8.0;
      t35 = 3.0*u0;
      t38 = dx*dx;
      t40 = 1/t38/dx;
      t48 = t5*u0;
      t54 = eps*eps;
      t58 = t5*t5;
      t70 = (u0-uold)/dt+(2.0*u0-3.0*t5-3.0*param_A*t5/t9)*t14*t16-param_D*(t21
*t14*t16-t27*(u0-um)*t16)*t16+(t21*(upp-3.0*up+t35-um)*t40-t27*(up-t35+3.0*um-
umm)*t40)*t16+2.0*t48*param_A/t9/t8+param_beta*(t54/t48-t54*eps/t58-param_p0-(
up-0.2E1*u0+um)/t38)/(u0+param_K);


			rhs[i]= t70;

			/*-------------------------*/


			/*--------------*/
	t1 = u0+um;
      t5 = dx*dx;
      t6 = t5*t5;
      t8 = t1*t1*t1/t6/8.0;

			J[i][0]=t8;

			/*--------------*/
	t1 = u0+um;
      t2 = t1*t1/4.0;
      t5 = 1/dx;
      t9 = t1*t2/2.0;
      t14 = u0+up;
      t18 = dx*dx;
      t20 = 1/t18/dx;
      t37 = -param_D*(-3.0/2.0*t5*t2*(u0-um)+t9*t5)*t5+(-t14*t14*t14*t20/8.0
-3.0/2.0*t2*(up-3.0*u0+3.0*um-umm)*t20-3.0*t9*t20)*t5-param_beta/t18/(u0+
param_K);

			J[i][1]=t37;

			/*--------------*/
	t4 = L0+x0;
      t5 = t4*t4;
      t6 = 1/t5;
      t10 = up-u0;
      t12 = 1/dx;
      t15 = u0*u0;
      t17 = param_A*t15;
      t22 = u0+up;
      t23 = t22*t22/4.0;
      t28 = t23*t22/2.0;
      t30 = u0+um;
      t31 = t30*t30/4.0;
      t37 = t31*t30/2.0;
      t43 = 3.0*u0;
      t46 = dx*dx;
      t48 = 1/t46/dx;
      t66 = eps*eps;
      t67 = t15*t15;
      t68 = 1/t67;
      t71 = t66*eps;
      t76 = 1/t46;
      t80 = u0+param_K;
      t92 = t80*t80;
      t95 = 1/dt+(-6.0*param_A*u0*t6-6.0*u0+2.0)*t10*t12-(-3.0*t17*t6-3.0*t15+
2.0*u0)*t12-param_D*(3.0/2.0*t23*t10*t12-t28*t12-3.0/2.0*(u0-um)*t31*t12-t37*
t12)*t12+(3.0/2.0*t23*(upp-3.0*up+t43-um)*t48+3.0*t28*t48-3.0/2.0*t31*(up-t43+
3.0*um-umm)*t48+3.0*t37*t48)*t12+6.0*t17/t5/t4+param_beta*(-3.0*t66*t68+4.0*t71
/t67/u0+0.2E1*t76)/t80-param_beta*(t66/t15/u0-t71*t68-param_p0-(up-0.2E1*u0+um)
*t76)/t92;


			J[i][2]=t95;

			/*--------------*/
	t2 = u0*u0;
      t6 = pow(L0+x0,2.0);
      t11 = 1/dx;
      t13 = u0+up;
      t14 = t13*t13/4.0;
      t20 = t13*t14/2.0;
      t29 = dx*dx;
      t31 = 1/t29/dx;
      t36 = u0+um;
      t48 = (2.0*u0-3.0*t2-3.0*param_A*t2/t6)*t11-param_D*(3.0/2.0*t14*(up-u0)*
t11+t20*t11)*t11+(3.0/2.0*t14*(upp-3.0*up+3.0*u0-um)*t31-3.0*t20*t31-t36*t36*
t36*t31/8.0)*t11-param_beta/t29/(u0+param_K);


			J[i][3]=t48;

			/*--------------*/
	t1 = u0+up;
      t5 = dx*dx;
      t6 = t5*t5;
      t8 = t1*t1*t1/t6/8.0;

			J[i][4]=t8;

			/*-------------------------*/
			/*-------------------------*/

		}
	}
	
	
	rhs[0] = (u[1]-u[0])/dx + 1.0/param_D;

    J[0][2] = - 1.0/dx;
    J[0][3] = 1.0/dx; 

	
	rhs[1] = (u[2]-u[0])/(2.0*dx) + 1.0/param_D;
	J[1][3] = 1.0/(2.0*dx);
	J[1][1] = -1.0/(2.0*dx);
	

	rhs[N] = u[N]-uR;
	J[N][2] = 1.0;
	J[N-1][2]+=J[N-1][4];

	for(i=0;i<N+1;++i)
		rhs[i]*= -1.0;

}

void memory_setup()
{
        U=myalloc(N+1);
        U0=myalloc(N+1);
        U0_record=myalloc(N+1);
        U0_init=myalloc(N+1);

        rhs=myalloc(N+1);
        J=matrix(N+1,5);
        J1=matrix(N+1,5);
        indx=myalloci(N+1);
}




