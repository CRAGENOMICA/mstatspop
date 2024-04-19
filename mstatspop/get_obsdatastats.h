/*
 * get_obsdatastats.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gvera
 */

#ifndef GET_OBSDATASTATS_H_
#define GET_OBSDATASTATS_H_

#include "common.h"

#ifdef	__cplusplus
extern "C" {
#endif

#include "zutil.h"

int get_obsstats(FILE *file_output,SGZip *file_output_gz,
                 FILE *file_mask,
                 // FILE *file_logerr,SGZip *file_logerr_gz,
                 int n_samp, long int n_site,
				 long int *n_realsite,char **names,char *DNA_matr,double *matrix_sizepos,
				 double *matrix_segrpos,char **matrix_pol,long int **matrix_freq,
				 long int **matrix_pos,double *length_al,long int *length_seg,
				 int *nsamuser,int npops,double *svratio,double *missratio, 
				 int include_unknown,double *sum_sam,double **tcga,long int **matrix_sv,
				 long int *nmhits,int output,char *ploidy,int outgroup_presence,
				 double *nsites1_pop, double *nsites1_pop_outg,
				 double *nsites2_pop,double *nsites2_pop_outg,double *nsites3_pop,double *nsites3_pop_outg, 
				 double *anx, double *bnx,double *anxo, double *bnxo,
				 double **lengthamng, double **lengthamng_outg, long int *mhitbp,char **matrix_pol_tcga,long int formatfile);

#ifdef	__cplusplus
}
#endif



#endif /* GET_OBSDATASTATS_H_ */
