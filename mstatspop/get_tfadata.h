//
//  get_tfadata.h
//  xcode_project
//
//  Created by Sebastian Ramos-Onsins on 09/08/15.
//
//

#ifndef __xcode_project__get_tfadata__
#define __xcode_project__get_tfadata__

#include <stdio.h>

#include "common.h"
#include "usegff.h"
#include "get_obsdatastats.h"

#ifdef	__cplusplus
extern "C" {
#endif
    
#define MAXLEN 1000

int read_coordinates(FILE *file_wcoor, SGZip *file_wcoor_gz, FILE *file_output, SGZip *file_output_gz, FILE *file_logerr, SGZip *file_logerr_gz, long int **wgenes, long int *nwindows,char *chr_name);
int read_weights_positions_file(FILE *file_ws, SGZip *file_ws_gz, struct SGZIndex *index_w,FILE *file_logerr,SGZip *file_logerr_gz, double **wP, double **wPV, double **wV, long int *wlimit_end,long int init_site, double *window_size, long int *n_sitesw, int weight_window, char *chr_name,int first, long int length);
int function_read_tfasta(FILE *file_input,SGZip *input_gz,struct SGZIndex *index_input,FILE *file_logerr,SGZip *file_logerr_gz,long int init_site,long int end_site,int *n_sam, long int *n_site, char ***names, char **DNA_matr,char **matrix_pol_tcga,char *chr_name,int first,long int length);
int transform_beg_chr(char *ID, char *chr_name, long int beg,int nchstr, int count0s);
int check_comment(char *c, FILE *file_input, SGZip *file_input_gz);
    
#ifdef	__cplusplus
}
#endif

#endif /* defined(__xcode_project__get_tfadata__) */
