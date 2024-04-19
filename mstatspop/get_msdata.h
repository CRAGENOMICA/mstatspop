/*
 * get_msdata.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gvera
 */

#ifndef GET_MSDATA_H_
#define GET_MSDATA_H_

#include "common.h"
#include "util.h"
#include "ran1.h"
#include "zutil.h"

#ifdef	__cplusplus
extern "C" {
#endif

int get_msdata( FILE *file_input,SGZip *input_gz,FILE *file_logerr,SGZip *file_logerr_gz,
		char **matrix_pol, long int **matrix_freq, long int **matrix_pos,
		long int *length_seg, int *nsamuser, int npops, int nsamtot, long int length, long int *nmhits,
		int *matrix_mask, float *vector_mask, float svratio, double **vector_priors, int *npriors,
		long int **matrix_sv, int outgroup_presence, int force_outgroup, float freq_revert ,double *sum_sam,
		double *nsites1_pop, double *nsites1_pop_outg, int formatfile,double *nsites2_pop,double *nsites2_pop_outg,
		double *nsites3_pop,double *nsites3_pop_outg,
		double *anx, double *bnx,double *anxo, double *bnxo,double **lengthamng,double **lengthamng_outg,int include_unknown,
		char *file_mas, double freq_missing_ms, int kind_length,
		double **sum_sam_mask, double *length_mask, long int *length_mask_real, double *missratio, int location_missing_ms, int *sort_nsam);

#ifdef	__cplusplus
}
#endif


#endif /* GET_MSDATA_H_ */
