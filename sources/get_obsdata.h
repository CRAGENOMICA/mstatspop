/*
 * get_obsdata.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gvera
 */

#ifndef GET_OBSDATA_H_
#define GET_OBSDATA_H_

#include "common.h"
#include "usegff.h"
#include "get_obsdatastats.h"
#include "zutil.h"

#ifdef	__cplusplus
extern "C" {
#endif

int get_obsdata(FILE *file_output,
		        FILE *file_input, SGZip *input_gz,
				FILE *file_mask,
				char *name_fileinputgff,int gfffiles,char *subset_positions,
				char *genetic_code,char **matrix_pol,long int **matrix_freq,
				long int **matrix_pos, long int *length, long int *length_seg,
				double *length_al, long int *length_al_real, int mainargc,
				char *ploidy,int *nsamuser,int npops,double *svratio,
				double *missratio,int include_unknown, double *sum_sam,
				double **tcga, long int **matrix_sv,long int *nmhits,
				int output,int outgroup_presence,char *criteria_transcript,
				double *nsites1_pop, double *nsites1_pop_outg,
				double *nsites2_pop,double *nsites2_pop_outg,double *nsites3_pop,double *nsites3_pop_outg,
				double *anx, double *bnx,double *anxo, double *bnxo,
				double **lengthamng, double **lengthamng_outg,int *sort_nsam,char **matrix_pol_tcga
);

int var_char(FILE *file_input,SGZip *input_gz,long int *count,int *c,int *n_sam,long int *n_sit,int *nseq,int *maxsam,char ***names,char **DNA_matr,
	long int *n_site,int excludelines,char *name_excluded,int *n_excl,int includelines,char *name_ingroups,char *name_outgroup,int outgroup,int nsamuser_eff,char *ploidy);

int assigna(FILE *file_input,SGZip *input_gz,int *c,int *nseq,int *maxsam,char ***names);

#ifdef	__cplusplus
}
#endif



#endif /* GET_OBSDATA_H_ */
