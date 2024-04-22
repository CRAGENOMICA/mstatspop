/*
 *  mstatspop.h
 *  mstatspop
 *
 *  Created by Sebastian Ramos-Onsins on 11/08/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef MSTATSPOP_H
#define	MSTATSPOP_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "common.h"

#include "calcR2.h"
#include "calcFs.h"
#include "calc_Toptimal_tests.h"
#include "freq_stats.h"
#include "missing_freqs.h"
#include "jointfreqdist.h"
#include "mismatch.h"
#include "permute.h"
#include "print_output.h"
#include "ran1.h"
#include "util.h"
#include "zutil.h"
#include "zindex.h"
#include "get_msdata.h"    
#include "sancestral.h"
#include "fstcalc.h"
#include "get_obsdata.h"
#include "get_tfadata.h"
/*introduce the FASTA file and the number of populations (with sample sizes)*/
// int get_obsdata( FILE *,SGZip *,FILE *,SGZip *,FILE *, SGZip *,FILE *,char *,int ,char *,char *,char **,
// 						long int **,long int **, long int *, long int *, double *,
// 						long int *,int,char *,int *,int,double *,double *,int,
// 						double *,double **,long int **,long int *,int,int,char *,
// 						double *,double *,double *,double *,double *,double *, double *,double *, 
// 						double *, double *, double **,double **,int *, char **,char *,int);
						
/*introduce the ms file and the number of populations (with sample sizes) and length*/			
    // int get_msdata( FILE *file_input,SGZip *,
	//  // FILE *file_logerr,SGZip *file_logerr_gz,
    //                char **matrix_pol, long int **matrix_freq, long int **matrix_pos,
    //                long int *length_seg, int *nsamuser, int npops, int nsamtot, long int length, long int *nmhits,
    //                int *matrix_mask, float *vector_mask, float svratio, double **vector_priors, int *npriors,
    //                long int **matrix_sv, int outgroup_presence, int force_outgroup, float freq_revert ,double *sum_sam,
    //                double *nsites1_pop, double *nsites1_pop_outg, int formatfile,double *nsites2_pop,double *nsites2_pop_outg,
    //                double *nsites3_pop,double *nsites3_pop_outg,
    //                double *anx, double *bnx,double *anxo, double *bnxo,double **lengthamng,double **lengthamng_outg,int include_unknown,
    //                char *file_mas, double freq_missing_ms, int kind_length,
    //                double **sum_sam_mask, double *length_mask, long int *length_mask_real, double *missratio, int location_missing_ms, int *sort_nsam);


/*calculate exclusive, fixed, shared among populations*/
// int calc_sxsfss(	int ,int *, char *,long int *, long int, struct stats *, 
// 						long int *, int,int); 



/*calculate fst, piw, pia*/
// int calc_piwpiafst(int, int, int, int *, char *, long int, struct stats *,
// 						long int *,int, int);

/*calculate fstH, Hapw, Hapa*/
// int calc_hwhafsth (int, int *, char *, long int, struct stats *);	

// int read_coordinates(FILE *file_wcoor, SGZip *file_wcoor_gz, FILE *file_output, SGZip *file_output_gz, FILE *file_logerr, SGZip *file_logerr_gz, long int **wgenes, long int *nwindows,char *chr_name);

int get_tfadata(FILE *file_output,
    SGZip *file_output_gz,
	FILE *file_input,
	SGZip *input_gz,
    struct SGZIndex *index_input,
    char *file_wps,
    FILE *file_ws,
    SGZip *file_ws_gz,
    struct SGZIndex *index_w,
    //FILE *file_logerr,
    //SGZip *file_logerr_gz,
    char **matrix_pol,
	long int **matrix_freq,
	long int **matrix_pos,
	long int *length,
	long int *length_seg,
	double *length_al,
	long int *length_al_real,
	int mainargc,
	int *nsamuser,
	int npops,
	double *svratio,
	double *missratio,
	int include_unknown,
	double *sum_sam,
	double **tcga,
	long int **matrix_sv,
	long int *nmhits,
	int output,
	int outgroup_presence,
	double *nsites1_pop,
	double *nsites1_pop_outg,
	double *nsites2_pop,
	double *nsites2_pop_outg,
	double *nsites3_pop,
	double *nsites3_pop_outg,
	double *anx,
	double *bnx,
	double *anxo,
	double *bnxo,
    double **lengthamng,
    double **lengthamng_outg,
	int *sort_nsam,
	long int *wgenes,
	long int nwindows,
    long int first_slide,
	long int slide,/**/
	long int window,/**/
	int Physical_length,
	long int *li/**/,
    int *npriors,
    double **vector_priors,
    char **matrix_pol_tcga,
    char *chr_name,
    int first
);
void usage(void);
int read_index_file(char *chr_name_all, unsigned long *nscaffolds,char ***chr_name_array,char ***chr_length_array);
    
#ifdef	__cplusplus
}
#endif

#endif	/* MSTATSPOP_H */
