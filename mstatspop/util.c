/*
 * util.c
 *
 *  Created on: Apr 18, 2012
 *      Author: gvera
 *
 *  Functions collected though previous versions of divpop or mstatspop... still growing.
 */


/* Optional available! but be careful!
 *
 * #define MIN(a,b)  ((a) < (b) ? (a) : (b))
 * it breaks with: min(3,++i)
 */

#include "util.h"

#define PI 3.14159265359

double min(int a, int b)
{
	return((a < b)? (a):(b));
}

double max(int a, int b)
{
	return((a > b)? (a):(b));
}

double gammln(double zz)
{
	/*Based on Numerical Recipes in C, Press et al. 1992. p. 213. and on 
	 Lanczos 1964 J. SIAM Numer. Anal. Ser. B, Vol. 1 pp 86-96.*/
	
	/*gamma distribution for a z integer*/
	double loggammaz;
	double z,logg,h,sumc;
	static double gamma = 5.0;
	static double c0 =  1.000000000178;
	static double c1 = 76.180091729406;
	static double c2 = 86.505320327112;
	static double c3 = 24.014098222230;
	static double c4 =  1.231739516140;
	static double c5 =  0.001208580030;
	static double c6 =  0.000005363820;
	
	if(zz <= 0.) {
		puts("Error gamma");
		return (double)-10000.;
	}
	
	z = (double)zz;
	h = (double)sqrt(2. * PI);
	sumc = c0 + c1/(z+1.) - c2/(z+2.) + c3/(z+3.)  - c4/(z+4.) + c5/(z+5.) - c6/(z+6.);
	logg = (z + 0.5)*(double)log((double)(z + gamma + 0.5)) - (z + gamma + 0.5);
	loggammaz = log((double)h);
	loggammaz += logg + log((double)sumc);
	loggammaz -= log((double)z);
	
	return (double)loggammaz;
}


double largebinomialdist(double pp, double n) 
{
	/*Based on Numerical Recipes in C, Press et al. 1992 and on
	 Fishman 1979 J. American Statistical Association Vol. 74, No. 366, pp 418-423*/
	
	double ran1(void);	
	double p,np;
	int N;
	double A,B,C,D,V,s;
	double g,plog,pclog,sq,angle,y,em,tt;
	double gammln(double);	
	
	if(pp > 0.5) p = (double)1.-pp;
	else p = pp;
	
	np = n * p;
	
	if(n==0) {
		puts("Error bindist");
		return (double)-10000.;
	}
	if(p==(double)0) {
		if(pp > 0.5) return (double)n;
		return (double)0;
	}
	
	if(np < (double)10) {
		/*Rejection Method: BI Algorithm*/
		s = (double)1- p;
		A = (double)1;
		B = p/s;
		C = ((double)n+(double)1)*B;
		D = A;
		N = 0;
		V = ran1()/(double)pow(s,(double)n);
		while(V > A) {
			N++;
			D *= (C/(double)N - B);
			A += D;
			if(N > n) break;
		}
	}
	else { /*Rejection method with a Lorentzian comparison distribution*/
		g = gammln(n+1.);
		plog  = (double)log(p);
		pclog = (double)log((1.0 - p));
		sq = (double)sqrt(2.0*np*(1.0 - p));
		do {
			do {
				angle = PI*ran1();
				y = tan(angle);
				em=sq*y+np;
			} while(em < 0.0 || em >= (n + 1.0));
			em = floor(em);
			tt = 1.2*sq*(1.0+y*y)*exp(g-gammln(em+1.0) - gammln(n-em+1.0)+em*plog+(n-em)*pclog);
		} while(ran1() > tt);
		N = (int)em;
	}
	
	if(pp > 0.5) N = (int)n - N;
	return (double)N;
}

