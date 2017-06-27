/*
 * print_output.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gvera
 */

#ifndef PRINT_OUTPUT_H_
#define PRINT_OUTPUT_H_

#include "common.h"
#include "util.h"
#include "missing_freqs.h"

#ifdef	__cplusplus
extern "C" {
#endif

#define TO_NEW	0
#define DENCOV_CORRECTION	1
#define MATRIXCOV	0
#define LUCA_CRs	1

/* prints out ALL the results... */
int print_output( int mainargc,int npops,int *nsam,
					FILE *file_out, SGZip *file_out_gz,
                    char *file_input, char *file_output,
					int gfffiles, char *file_GFF, char *subset_positions,
					char *code_name, char *genetic_code,
					long int length, long int length_seg,
					double length_al,long int length_al_real,
					struct stats *statistics, struct probs *piter,
					long int niter, long int *sites_matrix, char *ploidy,
					double svratio, double missratio, int include_unknown,
					long int *matrix_pos,
					double **jfd, int **nfd, int output, int H1frq,int H0frq,
					long int nseed, char *file_H1f, char *file_H0f,
					double *vector_priors, int npriors, int formatfile,
					int outgroup_presence,int force_outgroup,double freq_missing_ms,
					double *nsites1_pop, double *nsites1_pop_outg,
				    double *nsites2_pop, double *nsites2_pop_outg,
					double *nsites3_pop, double *nsites3_pop_outg,
					long int niterdata, char *matrix_pol, int *r2i_ploidies,
                    char *matrix_pol_tcga,char *chr_name);

#ifdef	__cplusplus
}
#endif



#endif /* PRINT_OUTPUT_H_ */
