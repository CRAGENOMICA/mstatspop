/*
 * common.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gvera
 */
 
#ifndef COMMON_H_
#define COMMON_H_

#ifdef	__cplusplus
extern "C" {
#endif

#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define MSTATSPOP "\nmstatspop v.0.1beta (20170224)\n" \
		   "Sebastian E. Ramos-Onsins, Luca Ferretti, Emanuele Raineri, Giacomo Marmorini, William Burgos, Joan Jene and Gonzalo Vera\n" \
		   "Variability Analyses of multiple populations: " \
		   "Calculation and estimation of statistics and neutrality tests.\n"


#define MSP_MAX_FILENAME			(unsigned long) 1024 /**< @brief Maximum Filename Length allowed */
#define MSP_MAX_GFF_WORDLEN         (unsigned long) 20
#define MSP_GENETIC_CODETYPE_LEN	(unsigned long) 50	/* e.g. "Nuclear universal" */
#define MSP_GENCODE_COMBINATIONS    (unsigned long) 64 /* 4^3 */
#define MSP_GFF_CRITERIA_MSG_LEN    (unsigned long) 20 /* e.g. "MIN" */
#define MSP_MAX_FILELINE_LEN		(unsigned long) 102400

#define SAMPLE_LARGE 4000

/* For compatibility for some old compilers */
#ifndef NULL
#define NULL	0
#endif

/* # TODO: Hacer un ENUM!
-o [output: 0:extended, 10: complete extended, 1:single line, 2:single line freq spectrum, 3:single line joint freq distrib.\n");
	printf("\t				4:single line pairwise distribution 5:single line frequency variant per line\n");
	*/

struct stats {
	long int *Sanc;
	double *piw;
	double *pia;
	double *piT;
	double *piant;
	double *piTnt;
	double *fst;
	double *piwHKY;
	double *piaHKY;
	double *piTHKY;
	double *fstHKY;
	double *fst1all;
	double *Gst;
	double *hapw;
	double *hapa;
	double *hapT;
	double *fsth;
	double *fsth1all;
		
	double fstALL;
	double fsthALL;
	double GstALL;

	double *S;
	double *So;
	double *thetaS;
	double *thetaSo;
	double *thetaT;
	double *thetaTo;
	double *thetaTHKY;
	double *thetaFL;
	double *thetaFW;
	double *thetaL;
	double *thetaSA;
	double *thetaTA;
	double *K;
	double *KHKY;
	double *Dtaj;
	double *Dfl;
	double *Ffl;
	double *Hnfw;
	double *Ez;
	double *Yach;
	double *R2;
    double *Fs;
    double *FH;

	long int *Rm;
	double *ZnA;

	long int **freq;
	int nh;
	int *nhpop;
	long int **freqh;

	double *length;
	double *length2;
    double **lengthamng;
    double **lengthamng_outg;
	double total_length;
	double total_real_length;

	double total_svratio;
	double *total_tcga;
	double **tcga;
	double ***sv;
	double ***svT;
	long int nmhits;

	double **H1freq;
	double *thetaH1;
	double **H0freq;
	double *thetaH0;

	double *ToH0_ii;
	double *ToH0_00;
	double *To_ii;
	double *To_00;
	double *To_i0;
	double *To_Qc_ii;
	double *To_Qw_ii;
	double *To_Lc_ii;

	double *mdsd; /*mismatch distribution: standard deviation*/
	double *mdg1;/*mismatch distribution: skewness*/
	double *mdg2;/*mismatch distribution: kurtosis*/
	double **mdw; /*the complete mismatch distribution*/

	double **linefreq; /*frequency of variants (calculated within population) that contain each line*/
	
	double *anx;
	double *bnx;
	double *anxo;
	double *bnxo;
	
	double **R2p;
};

struct probs {
	long int *i;
	long int *ih;
	long int *igh;
	long int *i1;
	long int *ih1;
	long int *niteri;
	long int *niterih;
	long int *niterigh;
	long int *niteri1;
	long int *niterih1;
	long int iall;
	long int ihall;
	long int ighall;
	long int niteriall;
	long int niterihall;
	long int niterighall;
};

struct test_rinf {
	int n;
	double *w;
};

struct test_r0 {
	int n;
	double *w;
	double *dw;
	double *d2w;
	double theta;
	long l;
	double var;
	double dvar;
	double d2var;
};

struct test_quad_wc {
	int n;
	double theta;
	long l;
	double gamma;
	double *w;
	double **w2;

};

struct test_quad {
	int n;
	double *w;
	double theta;
	long l;
	double **w2;
};


struct test_lc {
	int n;
	double *w;
	double theta;
	long l;
	double gamma;
};

#ifdef	__cplusplus
}
#endif

#endif /* COMMON_H_ */
