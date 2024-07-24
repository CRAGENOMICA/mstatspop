/*
 * File:   freq_stats.h
 * Author: gvera
 *
 * Created on April 17, 2012, 2:46 PM
 */
 
#ifndef FREQ_STATS_H
#define	FREQ_STATS_H

#include "common.h"
#include "ran1.h"

#ifdef	__cplusplus
extern "C" {
#endif

struct covar_ij {
	int i1;
	int i2;
	int n1;
	int n2;
	int n12;
	double Covij;
};
	
int calc_freqstats(int npops, int *nsam, char *matrix_pol,long int length, struct stats *statistics,
				   int outgroup_presence,int force_outgroup,int include_unknown, int n_ccov, int H1frq);

void init_coef(double *p,int sample_size);

double tajima_d(double k_, long int S_, double *coef_taj);

double fl_f( long int fr1, long int S, double pi, double *coef);
double fl_d( long int fr1,long int S, double *coef);
double fay_wu_normalized2(int n,double thetaL,double thetaw,double S,double *coef,double pi);
double E_zeng(int n,double thetaL,double thetaw,double S,double *coef);
double fl_d2(int sample_size,long int fr1w,long int S, double *coef);
double fl_f2(int sample_size,long int fr1w, long int S, double pi, double *coef);

double freqtesto_achaz(int sample_size,long int *fr,int singleton,double *w1,double *w2);
double freqtestn_achaz(int sample_size,long int *fr,int singleton,double *w1,double *w2);

double an(int n);
double a2n(int n);
double bn(int n,int i);
double sigmaii(int n,int i);
double sigmaij(int n,int i,int j);
/*double omegai(int n,int i,double *w1,double *w2);*/
double omegai(int n,int i,double *w1,double *w2,double sumw1,double sumw2);
double psii(int n,int i);
double rhoii(int n,int i);
double rhoij(int n,int i,int j);
/*double omegain(int n,int i,double *w1,double *w2);*/
double omegain(int n,int i,double *w1,double *w2,double sumw1,double sumw2);

float cov_missing(int, int, int, int, int, float);
int ominx_tajD(double **ominx, int **nx, int *no, double any, int nsam, long int length, long int Sc);
int ominx_Dfl(double **ominx, int **nx, int *no, double any, int nsam, long int length, long int Sco);
int ominx_Ffl(double **ominx, int **nx, int *no, int nsam, long int length, long int Sco);
int ominx_Hnfw(double **ominx, int **nx, int *no, int nsam, long int length, long int Sco);
int ominx_Ezeng(double **ominx, int **nx, int *no, double any, int nsam, long int length, long int Sco);
int ominx_Yachaz(double **ominx, int **nx, int *no, double any, int nsam, long int length, long int Sco);
int ominx_FH(double **ominx, int **nx, int *no, double any, int nsam, long int length, long int Sco);
int ominx_tajDnooutg(double **ominx, int **nx, double any, int nsam, long int length, long int Sc);
int ominx_Dfl2(double **ominx, int **nx, double any, int nsam, long int length, long int Sc);
int ominx_Ffl2(double **ominx, int **nx, int nsam, long int length, long int Sc);
int ominx_Yachaz2(double **ominx, int **nx, double any, int nsam, long int length, long int Sc);
int ominx_Optimal_unfolded(double **ominx, int **nx, int *no, double any, int nsam, long int length, long int Sco, double *mean_freqsptr);
int ominx_Optimal_folded(double **ominx, int **nx, double any, int nsam, long int length, long int S, double *mean_freqsptr);
double freqtest_outg_missing(double **ominx,long int **eix,int **nx, int *no, double thetaw, int nsam, 
							 long int length, int n_ccov, struct covar_ij **Covi1i2, long int *count_cov,int singleton, 
							 double *bnxp, double *theta_square);
double freqtest_noutg_missing(double **ominx,long int **eix,int **nx, double thetaw, int nsam, 
							  long int length, int n_ccov, struct covar_ij **Covi1i2, long int *count_cov,int singleton, 
							  double *bnxp, double *theta_square);
float watterson_variance(float theta, int l, int* nx, int* ny,int* nxy, float* var_d, float* var_h);
long int maxd(double, double);
double ran1(void);
    
/*	
	double compute_Omega_Y( int n, long int *Xi_array, char norm );
	void compute_theta2estimator_omega(  int n,  long int *Xi_array, double *w, double *theta, double *theta2  );
	double compute_T( int n, long int *Xi_array, double *w1, double *w2, char norm, double theta, double theta2 );
	double * compute_HarmonicSums( int n );
	double beta_FU1995( int i, double *HarmonicSums, int n  );
	double sigma_ii_FU1995( int i, double *HarmonicSums, int n );
	double sigma_ij_FU1995( int i, int j, double *HarmonicSums, int n );
*/
#ifdef	__cplusplus
}
#endif

#endif	/* FREQ_STATS_H */

