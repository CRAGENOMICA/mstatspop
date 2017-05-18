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
#include "zutil.h"

#ifdef	__cplusplus
extern "C" {
#endif

int read_coordinates(FILE *file_wcoor, FILE *file_output,long int **wgenes, long int *nwindows);
int read_weights_positions_file(FILE *file_ws, SGZip *file_ws_gz,FILE *file_output, float **wP, float **wPV, float **wV,long int *wlimit_end);
int read_weights_file(FILE *file_es, FILE *file_output, float **wV, long int **Pp, long int *nV, long int *welimit_end);
int function_read_tfasta(FILE *file_input,SGZip *input_gz,long int init_site,long int end_site,int *n_sam, long int *n_site, char ***names, char **DNA_matr,char **matrix_pol_tcga);

#ifdef	__cplusplus
}
#endif

#endif /* defined(__xcode_project__get_tfadata__) */
