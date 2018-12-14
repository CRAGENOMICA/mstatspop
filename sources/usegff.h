/*
 * usegff.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gvera
 */

#ifndef USEGFF_H_
#define USEGFF_H_

#include "common.h"
#define SIZE_ROW 10000

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdio.h> /* para BUFSIZ??*/

/**
 * todo: documentarrrrr, ya casi esta
 */
struct valuesgff
{
	char filename[256];
	char source[256];
	char feature[256];
	char strand[1];
	long int start;
	long int end;
	char score[256];
	char frame[1];
	char seqname[256];
	char gene_id[256];
	char transcript_id[256];
};
/*
 name_fileinputgff: name of the file gff/gtf <- defined before
 subset_positions: string with the 'feature' value to filter: in case coding region use synonymous, nonsynonymous, silent. Also noncoding. <- defined before
 genetic_code: a 64 length char vector with the genetic code <- defined before
 matrix_sizepos: vector with 1/0 n_site values indicating whether use or not that position <- defined before but HERE MODIFIED
 n_samp: number of samples in data <- defined before. In case NGS include 2 samples per (diploid) individual or two per population in case pools.
 n_site: total number of positions <- defined before
 DNA_matr: all data in a matrix <- defined before. In case NGS include a matrix with n_site x 2 seqs per individual (or 2 per population in case pools)
 matrix_segrpos: matrix with only the variants <- defined before but HERE MODIFIED
 file_output: name of the file_output <- defined before
 mainargc: in fact NOT USED
 include_unknown: flag to say if include missing values or not <- defined before
 criteria_transcripts: char vector wit "max" or "min" values <- defined before
*/

int use_gff(char *name_fileinputgff,char *subset_positions,char *genetic_code,
			double *matrix_sizepos,int n_samp,long int n_site,char *DNA_matr,
			double *matrix_segrpos,FILE *file_output,SGZip *file_output_gz,int mainargc,
            FILE *file_logerr, SGZip *file_logerr_gz,int include_unknown,
			char *criteria_transcripts, int type_output, long int *nmhits, long int *mhitbp,
			int outgroup_presence, int nsamoutg,char *chr_name, int first);

int tripletnsamp(char *cod3n,char *DNA_matr,char strand,double *cmat,
					int n_samp,long int n_site,long int end,long int ii,
					FILE *file_output,SGZip *file_output_gz/*,int mainargc*/,int include_unknown,int type_output,
					long int *nmhits, long int *mhitbp, int outgroup_presence, int nsamoutg,
                    FILE *file_logerr, SGZip *file_logerr_gz);
	
int comp_trcpt_id(const void *a,const void *b);
int comp_start_id(const void *a,const void *b);
int comp_end_id(const void *a,const void *b);
int comp_gene_id(const void *a,const void *b);
/*do a function that read the genetic code and add the degeneration to each position at each codon*/
int function_do_nfold_triplets(int n_fold[64][3], char *genetic_code, char tripletsN[64][3]);

#ifdef	__cplusplus
}
#endif



#endif /* USEGFF_H_ */
