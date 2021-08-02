#include <iostream>
#include <fstream>
#include <math.h>
#include "Lambert.h"
//#include <cilk/cilk.h>

double PI  = 3.14159265358979323846264338327950288419716939937510;  
double PI2 = 2.*PI;

/*// ------------------------------------------------------------------------ //
// This source file is part of the 'ESA Advanced Concepts Team's			    //
// Space Mechanics Toolbox' software.                                           //
//                                                                              //
// The source files are for research use only,                                  //
// and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.         //
//                                                                              //
// Copyright (c) 2004-2007 European Space Agency                                //
// ------------------------------------------------------------------------   //


/*
 This routine implements a new algorithm that solves Lambert's problem. The
 algorithm has two major characteristics that makes it favorable to other
 existing ones.

   1) It describes the generic orbit solution of the boundary condition
   problem through the variable X=log(1+cos(alpha/2)). By doing so the
   graphs of the time of flight become defined in the entire real axis and
   resembles a straight line. Convergence is granted within few iterations
   for all the possible geometries (except, of course, when the transfer
   angle is zero). When multiple revolutions are considered the variable is
   X=tan(cos(alpha/2)*pi/2).

   2) Once the orbit has been determined in the plane, this routine
   evaluates the velocity vectors at the two points in a way that is not
   singular for the transfer angle approaching to pi (Lagrange coefficient
   based methods are numerically not well suited for this purpose).

   As a result Lambert's problem is solved (with multiple revolutions
   being accounted for) with the same computational effort for all
   possible geometries. The case of near 180 transfers is also solved
   efficiently.

   We note here that even when the transfer angle is exactly equal to pi
   the algorithm does solve the problem in the plane (it finds X), but it
   is not able to evaluate the plane in which the orbit lies. A solution
   to this would be to provide the direction of the plane containing the
   transfer orbit from outside. This has not been implemented in this
   routine since such a direction would depend on which application the
   transfer is going to be used in.

Inputs:
           r1=Position vector at departure (column) 
           r2=Position vector at arrival (column, same units as r1)
           t=Transfer time (scalar)
           mu=gravitational parameter (scalar, units have to be
           consistent with r1,t units) 
           lw=1 if long way is chosen    
           branch = 1 if the left branch is selected in a problem where N
           is not 0 (multirevolution)
           N=number of revolutions

Outputs:
           v1=Velocity at departure        (consistent units)
           v2=Velocity at arrival
           a=semi major axis of the solution
           p=semi latus rectum of the solution 
           theta=transfer angle in rad
           iter=number of iteration made by the newton solver (usually 6)

*/

void LambertI (const double *r1_in, const double *r2_in, double t, const double mu, const int lw, int N, const int branch, // INPUT
	           double *v1, double *v2, double &a, double &p, double &theta, int &iter) // OUTPUT

{
	using std::cout;
	using std::endl;
	const double tolerance=1e-11;
	int i, i_count;
	double r1[3], r2[3],r2_vers[3], alfa, beta, psi, eta, eta2, sigma1;
	double ih_dum[3], ih[3], dum[3];
	double	V, T, R=0.0, r2_mod = 0.0, dot_prod = 0.0,err,
	c,		        // non-dimensional chord
    s,		        // non dimesnional semi-perimeter
    am,		        // minimum energy ellipse semi major axis
    lambda,	        // lambda parameter defined in Battin's Book
	x1, x2, y1, y2,
	x, x_new=0,y_new, vr1, vt1, vt2, vr2;

	//printf("%s", "Start Lambert solver \n");
	
	for (i = 0; i < 3; i++)
    {
      r1[i] = r1_in[i];
      r2[i] = r2_in[i];
	  R += r1[i]*r1[i];
    }

	R = sqrt(R);
	V = sqrt(mu/R);
	T = R/V;

	// working with non-dimensional radii and time-of-flight
	t /= T;
	for (i = 0;i <3;i++)  // r1 dimension is 3
    {
		r1[i] /= R;
		r2[i] /= R;
		r2_mod += r2[i]*r2[i];
    }
	// Evaluation of the relevant geometry parameters in non dimensional units
	r2_mod = sqrt(r2_mod);
	
	for (i = 0;i < 3;i++)
      dot_prod += (r1[i] * r2[i]);

	theta = acos(dot_prod/r2_mod);
	if (lw)	theta=2*acos(-1.0)-theta;

	c = sqrt(1 + r2_mod*(r2_mod - 2.0 * cos(theta)));
	s = (1 + r2_mod + c)/2.0;
	am = s/2.0;
	lambda = sqrt (r2_mod) * cos (theta/2.0)/s;

double inn1, inn2;
if (N==0)
	{
	x1=log(0.4767);
	x2=log(1.5233);
	y1=log(x2tof(-.5233,s,c,lw, N))-log(t);
	y2=log(x2tof(.5233,s,c,lw, N))-log(t);

	// Newton iterations
	err=1;
	i_count=0;
	while ((err>tolerance) && (y1 != y2))
    {
		i_count++;
		x_new=(x1*y2-y1*x2)/(y2-y1);
		y_new=log(x2tof(exp(x_new)-1,s,c,lw, N))-log(t); 
		x1=x2;
		y1=y2;
		x2=x_new;
		y2=y_new;
		err = fabs(x1-x_new);
    }
	iter = i_count;
	x = exp(x_new)-1; 
	}
else
{
	if (branch == 1)   // left branch	  
	{
		inn1=-0.5234;
		inn2=-0.2234;
	}
	else			   // right branch
	{
		inn1=0.7234;
		inn2=0.5234;
	}
	x1=tan(inn1*acos(-1.0)/2);
    x2=tan(inn2*acos(-1.0)/2);
    y1=x2tof(inn1,s,c,lw,N)-t;
    y2=x2tof(inn2,s,c,lw,N)-t;
    err=1;
    i=0; int imax = 30;
  // Newton Iteration
	while ( (err>tolerance) && (y1 != y2) && i<imax )
	{
		i++;
        x_new=(x1*y2-y1*x2)/(y2-y1);
        y_new=x2tof(atan(x_new)*2/acos(-1.0),s,c,lw,N)-t;
        x1=x2;
        y1=y2;
        x2=x_new;
        y2=y_new;
        err=abs(x1-x_new);
	}
	x=atan(x_new)*2/acos(-1.0);

	if (i==imax) { 	iter = -1; }
	else { iter = i; }
}
	// The solution has been evaluated in terms of log(x+1) or tan(x*pi/2), we
    // now need the conic. As for transfer angles near to pi the lagrange
    // coefficient technique goes singular (dg approaches a zero/zero that is
    // numerically bad) we here use a different technique for those cases. When
    // the transfer angle is exactly equal to pi, then the ih unit vector is not
    // determined. The remaining equations, though, are still valid.
	a = am/(1 - x*x);			  // solution semimajor axis
	// psi evaluation
	if (x < 1)  // ellipse
    {
		beta = 2 * asin (sqrt( (s-c)/(2*a) ));
		if (lw) beta = -beta;
		alfa=2*acos(x);
		psi=(alfa-beta)/2;
		eta2=2*a*pow(sin(psi),2)/s;
		eta=sqrt(eta2);
    }
	else       // hyperbola
    {
		beta = 2*asinh(sqrt((c-s)/(2*a)));
		if (lw) beta = -beta;
		alfa = 2*acosh(x);
		psi = (alfa-beta)/2;
		eta2 = -2 * a * pow(sinh(psi),2)/s;
		eta = sqrt(eta2);
    }
	// parameter of the solution
	p = ( r2_mod / (am * eta2) ) * pow (sin (theta/2),2);
	sigma1 = (1/(eta * sqrt(am)) )* (2 * lambda * am - (lambda + x * eta));
	getCrossProductVectors(r1,r2,ih_dum);
	vers(ih_dum,ih) ;

	if (lw)
    {
		for (i = 0; i < 3;i++)
			ih[i]= -ih[i];
    }

	vr1 = sigma1;
	vt1 = sqrt(p);
	getCrossProductVectors(ih,r1,dum);
   
	for (i = 0;i < 3 ;i++)
		v1[i] = vr1 * r1[i] + vt1 * dum[i];
   
	vt2 = vt1 / r2_mod;
	vr2 = -vr1 + (vt1 - vt2)/tan(theta/2);
	
	vers(r2,r2_vers);
	getCrossProductVectors(ih,r2_vers,dum);
	for (i = 0;i < 3 ;i++)
		v2[i] = vr2 * r2[i] / r2_mod + vt2 * dum[i];
 
	for (i = 0;i < 3;i++)
    {
		v1[i] *= V;
		v2[i] *= V;
    }
	a *= R;
	p *= R;
	/*
	cout<<"a = "<<a<<endl;
	cout<<"p = "<<p<<endl;
	cout<<"v1[0] = "<<v1[0]<<endl;
	cout<<"v1[1] = "<<v1[1]<<endl;
	cout<<"v1[2] = "<<v1[2]<<endl;
	cout<<"v1    = "<<sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])<<endl<<endl;
	cout<<"v2[0] = "<<v2[0]<<endl;;
	cout<<"v2[1] = "<<v2[1]<<endl;
	cout<<"v2[2] = "<<v2[2]<<endl<<endl;
	cout<<"v2    = "<<sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2])<<endl<<endl;
	cout<<"a = "<<a<<endl;
	cout<<"p = "<<p<<endl;
	cout<<"theta = "<<theta*180/acos(-1.0)<<endl;
	cout<<"iter = "<<iter<<endl;
	*/
}	   


void GetVidMin(double rp, double vp, double dt,             // INPUT
			   double VidMin,								// OUTPUT
			   double *workd){								// working place 

//       workd[23]

		 workd[0] = 0.5*vp*dt;
		 workd[5] = sqrt(rp*rp+workd[0]*workd[0]);
		 workd[6] = vp;
		 workd[7] = workd[0]/workd[5];
		 workd[8] = rp/workd[5];

		 workd[9] = workd[6]*workd[7];
         workd[10] = workd[6]*workd[8];
         workd[11] = dt*dt;
         workd[12] = (-2.0*workd[9]+6.0*workd[5]/dt)/dt;
         workd[13] = -2.0*workd[10]/dt;
         workd[15] = (6.0*workd[9]-12.0*workd[5]/dt)/workd[11];
         workd[16] = 6.0*workd[10]/workd[11];
         workd[14] = workd[12]*workd[12]+workd[13]*workd[13];
         workd[17] = workd[15]*workd[15]+workd[16]*workd[16];
         workd[18] = sqrt(workd[17]);
         workd[19] = 2.0*(workd[12]*workd[15]+workd[13]*workd[16])/workd[17];
         workd[20] = workd[14]/workd[17];
         workd[1] = dt+0.5*workd[19];
         workd[2] = sqrt(workd[20]+dt*(workd[19]+dt));
         workd[3] = 0.5*workd[20]-0.125*workd[19]*workd[19];
         workd[4] = sqrt(workd[20]);
         workd[21] = 0.5*workd[1]*workd[2]+workd[3]*log(workd[1]+workd[2]);
         workd[22] = 0.25*workd[19]*workd[4]+workd[3]*log(0.5*workd[19]+workd[4]);
         VidMin = workd[18]*(workd[21]-workd[22]);

}


double asinh (double x) { return log(x + sqrt(x*x + 1)); };
double acosh (double x) { return log(x + sqrt(x*x - 1)); };
double x2tof(const double &x,const double &s,const double &c,const int lw, int N)
{
double am, a, alfa, beta;
//%Subfunction that evaluates the time of flight as a function of x
	am = s/2;
	a = am/(1-x*x);
	if (x < 1)//ellipse
    {
		beta = 2 * asin (sqrt((s - c)/(2*a)));
		if (lw) beta = -beta;
		alfa = 2 * acos(x);
    }
	
	else  //hyperbola
    
	{
		alfa = 2 * acosh(x);
		beta = 2 * asinh(sqrt ((s - c)/(-2 * a)));
		if (lw) beta = -beta;
    }

	if (a > 0)
    {
      return (a * sqrt (a)* ( (alfa - sin(alfa)) - (beta - sin(beta)) +N*2*acos(-1.0) ));
    }
	else
    {
      return (-a * sqrt(-a)*( (sinh(alfa) - alfa) - ( sinh(beta) - beta)) );
    }

}

void getCrossProductVectors(const double * vectorFirst, const double * vectorSecond, double *crossProductVector)
{
	crossProductVector[0] = (vectorFirst[1] * vectorSecond[2] - vectorFirst[2] * vectorSecond[1]);
	crossProductVector[1] = (vectorFirst[2] * vectorSecond[0] - vectorFirst[0] * vectorSecond[2]);
	crossProductVector[2] = (vectorFirst[0] * vectorSecond[1] - vectorFirst[1] * vectorSecond[0]);
}

void vers(const double *V_in, double *Ver_out)
{
	double v_mod = 0;
	int i;

	for (i = 0;i < 3;i++)
    {
		v_mod += V_in[i]*V_in[i];
	}
	
	double sqrtv_mod = sqrt(v_mod);

	for (i = 0;i < 3;i++)
	{
		Ver_out[i] = V_in[i]/sqrtv_mod;
	}
}


static double sign  (double a, double b)
{
  return (b < 0.0)? -fabs(a) : fabs(a);

} /* sign */
void vector_length (double *vec_1, double *vec_2, double &vec_length) {

	vec_length = 0;
	for (int i=0;i<3;i++){	
		vec_length +=vec_1[i]*vec_2[i];
	}
	vec_length = sqrt(vec_length);
}
void vector_dot    (double *vec_1, double *vec_2, double &vec_dot){

	vec_dot = 0;
	for (int i=0;i<3;i++){	
		vec_dot +=vec_1[i]*vec_2[i];
	}
}
void sigma (double u, double &out){

	out=4./3.+u*(2./5.+u*(3./14.+u*(5./36.+u*(35./352.+u*945./12480.))));

}
void dsigma(double u, double &out){

	out=2./5.+u*(3./7.+u*(5./12.+u*(35./88.+u*(945./2496.+u*10395./28800.))));

}
void Tx(double &x,double &q, int N, double &T, double &DTDx){
	double K,out1,out2,E,y,z,f,g,d;
// Lancaster's Algorithm
// initialize variables
	T=0;
	DTDx=0;
	y=0;
	z=0;
	f=0;
	g=0;
	d=0;
	K=q*q;
	E=x*x-1;
// if abs(x-1)<1e-2
// формула для производной не действительна в случае, если х=0 и k=1 и при х=1
// нужно переходить к формуле
// DtDx=2*x*(q*K**2*sigma(-K*E)-sigma(-E)), где
// sigma=sum(n*an*u**(n-1),n=1,ifinity)
// an=1...3...5...(2n-1)/2**(n-2)*(2*n+3)*n!

	if (abs(x-1) < 1.e-2) {
		sigma(-E,out1);
		sigma(-K*E,out2);
// T=sigma(-E)-q*K*sigma(-K*E);
		T=out1-q*K*out2;
		dsigma(-K*E,out1);
		dsigma(-E,out2);
		DTDx=2*x*(q*K*K*out1-out2);
	}
	else {
		y=sqrt(abs(E));
		z=sqrt(1.0+K*E);
		if (z==0) {
			z=2.2251-308;
		}
		f=y*(z-q*x);
		g=x*z-q*E;
		if (E<0) {
			d=N*PI+acos(g/sqrt(g*g+f*f));
		}
		else {
			d=log(f+g);
		}
		T=2*(x-q*z-d/y)/E;
		DTDx=(4-4*q*K*x/z-3*x*T)/E;
	}
}
int NorbT_1(double &q, double &T, int N, double &yy){
	double x,TxN,DTDxN1,DTDxN2,DTDxN3;
	double xx,zz,xN1,xN2,xN3,par1,par2;
//  yy is index of q that satisfies possibility of N revs
	yy=1.0;
	zz=0;
	xN1=0;
//  min energy time (good initial guess)
	Tx(xN1,q,N,TxN,DTDxN1);
//  N+1 rev possible for not aa
	if (T<TxN) {
		yy=1;
		zz=0;
		xN2=0.1;
		while (zz==0) {
			Tx(xN2,q,N,TxN,DTDxN2);
			if(T<TxN) {
				if (abs(DTDxN2)>1e-7) {
					zz=0;
					DTDxN3=DTDxN2-DTDxN1;
					xN3=DTDxN2*xN1-DTDxN1*xN2;
					par1=xN3*1.e-50;
					if ((par1-DTDxN3)>0) {
						yy=0;
						zz=1;
						return 0;
					}
// next x
					xN3=xN3/DTDxN3;
// Check found in MASL code
					par1=abs(xN2-xN3);
					par2=abs(xN1-xN3);
					if (par1<par2) {
						xN1=xN2;
						DTDxN1=DTDxN2;
						xN2=xN3;
					}
					else {
						yy=0;
						zz=1;
					}
				}
				else {
					yy=0;
					zz=1;
				}
			}
			else {
				yy=0;
				zz=1;
			}
		}
	}
	return 0;
}


