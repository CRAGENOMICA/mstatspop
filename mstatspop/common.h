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
#include "zutil.h"
#include "zindex.h"


// #define STRINGIFY(x) #x
// #define TOSTRING(x) STRINGIFY(x)
#define FULL_VERSION "v." VERSION_NUMBER " (" BUILD_NUMBER ")"

//uncomment FULL VERSION for compiling using build.sh!∫
#define MSTATSPOP "#mstatspop " //FULL_VERSION "\n" \
"Sebastian E. Ramos-Onsins, Luca Ferretti, Emanuele Raineri, Giacomo Marmorini, William Burgos, Joan Jene,  Gonzalo Vera and Ahmed Hafez\n" \
"Variability Analyses of multiple populations: " \
"Calculation and estimation of statistics and neutrality tests.\n"
#define MSTATSPOPVERSION "version 0.1-20251102\n"
#define MSTATSPOPTITLE "Variability Analyses of multiple populations: " \
"Calculation and estimation of statistics and neutrality tests.\n"
#define MSTATSPOPAUTHORS "Sebastian E. Ramos-Onsins, Luca Ferretti, Emanuele Raineri, Giacomo Marmorini, William Burgos, Joan Jene,  Gonzalo Vera and Ahmed Hafez\n"


#define MSP_MAX_FILENAME			(unsigned long) 1024 /**< @brief Maximum Filename Length allowed */
#define MSP_MAX_GFF_WORDLEN         (unsigned long) 20
#define MSP_GENETIC_CODETYPE_LEN	(unsigned long) 50	/* e.g. "Nuclear universal" */
#define MSP_GENCODE_COMBINATIONS    (unsigned long) 64 /* 4^3 */
#define MSP_GFF_CRITERIA_MSG_LEN    (unsigned long) 20 /* e.g. "MIN" */
#define MSP_MAX_FILELINE_LEN		(unsigned long) 102400
#define MSP_MAX_NAME                (unsigned long) 1024
#define MSP_MAX_COL_LINE            (unsigned long) 32767*2
#define SAMPLE_LARGE                (unsigned long) 4000

/* For compatibility for some old compilers */
#ifndef NULL
#define NULL	0
#endif

/* # TODO: Hacer un ENUM!
 -o [output: 0:extended, 10: complete extended, 1:single line, 2:single line freq spectrum, 3:single line joint freq distrib.\n");
 printf("\t				4:single line pairwise distribution 5:single line frequency variant per line\n");
 */

// possible formats for the input file
#define FASTA_FORMAT 0
#define NBRF_FORMAT 0
#define MS_FORMAT 1
#define MS_X_FORMAT 2
#define TFA_FORMAT 3

/* o Output type, from 0 to 11 - TODO: define*/

//  0 (extended),
#define OUTPUT_EXTENDED 0

//  1 (single line/window)
#define OUTPUT_SINGLE_LINE_WINDOW 1

//  2 (single line SFS/window)
#define OUTPUT_SINGLE_LINE_SFS_WINDOW 2

//  3 (dadi-like format)
#define OUTPUT_DADI_LIKE 3

//  4 (single line pairwise distribution)
#define OUTPUT_SINGLE_LINE_PAIRWISE_DISTRIBUTION 4

//  5 (single line freq. variant per line/window)
#define OUTPUT_SINGLE_LINE_FREQ_VARIANT_PER_LINE_WINDOW 5

//  6 (SNP genotype matrix)
#define OUTPUT_SNP_GENOTYPE_MATRIX 6

//  7 (SweepFiinder format -only first pop-)
#define OUTPUT_SWEEP_FINDER 7

//  8 (single line/window: Frequency of each haplotype in the populations)
#define OUTPUT_SINGLE_LINE_FREQ_HAPLOTYPE_POP 8

//  9 (single line/window: Frequency of variants per line and population)
#define OUTPUT_SINGLE_LINE_FREQ_VARIANT_POP 9

//  92 (single line/window: -rSFS- Frequency of variants per population relative to all)
#define OUTPUT_SINGLE_LINE_FREQ_VARIANT_POP_RELATIVE 92

//  10 (full extended)
#define OUTPUT_FULL_EXTENDED 10

typedef struct
{
    const char *format_arg; // format of the input file [fasta,ms,ms_x,tfa]
    int formatfile;         // format type of the input file [0,1,2,3]
    int output;             /* TODO: tipo de output, de 0 a 11 */
    char ploidy[2];         // ploidy, 1: haploid, 2: diploid
    long int niterdata;     // number of iterations (ms) or permutations (others)
    char file_in[MSP_MAX_FILENAME] /* input file */;
    char file_out[MSP_MAX_FILENAME]; /* output file */
    char file_mas[MSP_MAX_FILENAME]; /* mask file */
    char file_GFF[MSP_MAX_FILENAME]; /* GFF file */
    
    /* a Name of file of the Alternative frequency spectrum */
    /* Only with optimal tests (in GSL libs) */
    char file_H1f[MSP_MAX_FILENAME];
    
    /* a Name of file of the NULL frequency spectrum */
    /* Only with optimal tests (in GSL libs) */
    char file_H0f[MSP_MAX_FILENAME];
    /* 1: missing values allowed , 0: excluded missing positions */
    int include_unknown;
    int include_rsfs;
    //int include_fstall;
    
    /* variables defining more data*/
    int npops; // TODO args.npops
    /* Number of samples for each population, each element a population */
    int *vint_perpop_nsam; /* old nsam */
    
    /* Sum of all nsam of all populations */
    int int_total_nsam; /* old nsamtot*/
    
    int coordfile;
    
    /* number of iterations (ms) or permutations (others) */
    long int niter;
    
    /* GFF variables */
    int gfffiles;
    /*int 	observed_data	= 0;*/
    char subset_positions[MSP_MAX_GFF_WORDLEN];
    char code_name[MSP_GENETIC_CODETYPE_LEN];
    char genetic_code[MSP_GENCODE_COMBINATIONS + 1];
    char criteria_transcript[MSP_GFF_CRITERIA_MSG_LEN];
    
    long int length; /* Sequence length */
    int kind_length; // TODO args.kind_length
    
    float ms_svratio;        // TODO args.ms_svratio
    int outgroup_presence;   // TODO args.outgroup_presence
    int force_outgroup;      // TODO args.force_outgroup
    float freq_revert;       // TODO args.freq_revert
    double freq_missing_ms;  // TODO args.freq_missing_ms
    int location_missing_ms; // TODO args.location_missing_ms
    
    /*R2_p*/
    int *r2i_ploidies; // TODO args.r2i_ploidies
    /*Covariances in missing values*/
    int n_ccov;                                              // TODO args.n_ccov
    int int_total_nsam_order;                                /*number of samples. value for ordering the samples in populations in case they are not consecutive*/
    int *sort_nsam; /*vector with the position in the file*/ // TODO : args.sort_nsam
    long int slide;                                          // TODO args.slide
    long int first_slide;                                    // TODO args.first_slide
    long int window;                                         // TODO args.window
    char file_wps[MSP_MAX_FILENAME];                         // TODO args.file_wps
    char file_Wcoord[MSP_MAX_FILENAME];                      // TODO args.file_Wcoord
    int Physical_length;                                     // TODO args.Physical_length
    int mask_print;                                          // TODO args.mask_print
    char file_chr_name_all[MSP_MAX_NAME];                    // TODO :: args.file_chr_name_all
    long int nseed;
    int argc;
} mstatspop_args_t;


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
    
    double **linefreq; /*frequency of variants (in relation to total population, within population) that contain each line*/
    double **popfreq; /*frequency of variants (in relation to all populations) that contain each pop*/
    
    double *anx;
    double *bnx;
    double *anxo;
    double *bnxo;
    
    double **R2p;
    
    /* relative theta stats */
    double *rS;
    double *rSo;
    double *rthetaS;
    double *rthetaT;
    double *rthetaSo;
    double *rthetaTo;
    double *rthetaFL;
    double *rthetaFW;
    double *rthetaL;
    double *rthetaSA;
    double *rthetaTA;
    double *rDtaj;//rthetaT - rthetaS;
    double *rDfl;//rthetaS - rthetaFL;
    double *rFfl;//rthetaT - rthetaFL;
    double *rHnfw;//rthetaT - rthetaFW;
    double *rEz;//rthetaL - rthetaS;
    double *rYach;
    double *rFH;
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

struct rweights {
    /*weightss*/
    double **ww;
    double **wt;
    double **wfw;
    double **wfl;
    double **wl;
    double **wwA;
    double **wtA;
    double **wfwA;
    double **wflA;
    double **wlA;
    double **wphi_i;
    double ***wpsi_ij;
};

#ifdef	__cplusplus
}
#endif

#endif /* COMMON_H_ */
