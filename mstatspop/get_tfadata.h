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
#include "tfasta.h"

#ifdef __cplusplus
extern "C"
{
#endif

#define MAXLEN 1000

    int read_coordinates(
        FILE *file_wcoor,
        SGZip *file_wcoor_gz,
        FILE *file_output,
        SGZip *file_output_gz,
        // FILE *file_logerr, SGZip *file_logerr_gz,
        long int **wgenes,
        long int *nwindows,
        char *chr_name);
    int read_weights_positions_file(
        wtfasta_file *wtfasta,
        //FILE *file_ws,
        //SGZip *file_ws_gz,
        //struct SGZIndex *index_w,
        // FILE *file_logerr, SGZip *file_logerr_gz,
        double **wP,
        double **wPV,
        double **wV,
        long int *wlimit_end,
        long int init_site,
        double *window_size,
        long int *n_sitesw,
        int weight_window,
        char *chr_name,
        int first,
        long int length);

    // Prepare to be Removed
    // int function_read_tfasta_2(
    //     FILE *file_input,
    //     SGZip *input_gz,
    //     struct SGZIndex *index_input,
    //     // FILE *file_logerr, SGZip *file_logerr_gz,
    //     long int init_site,
    //     long int end_site,
    //     int *n_sam,
    //     long int *n_site,
    //     char ***names,
    //     char **DNA_matr,
    //     char **matrix_pol_tcga,
    //     char *chr_name,
    //     int first,
    //     long int length);

    // int function_read_tfasta(
    //     FILE *file_input,
    //     SGZip *file_input_gz,
    //     struct SGZIndex *index_input,
    //     // FILE *file_logerr,SGZip *file_logerr_gz,
    //     long int init_site,
    //     long int end_site,
    //     int *n_sam,
    //     long long *n_site,
    //     char ***names,
    //     char **DNA_matr,
    //     char **matrix_pol_tcga,
    //     char *chr_name,
    //     int first,
    //     long int length);

    int transform_beg_chr(char *ID, char *chr_name, long int beg, int nchstr, int count0s);
    int check_comment(char *c, FILE *file_input, SGZip *file_input_gz);

    int get_tfadata(
        FILE *file_output,
        SGZip *file_output_gz,
        tfasta_file *tfasta,
        // FILE *file_input,
        // SGZip *input_gz,
        // struct SGZIndex *index_input,
        wtfasta_file *wtfasta,
        //char *file_wps,
        //FILE *file_ws,
        //SGZip *file_ws_gz,
        //struct SGZIndex *index_w,
        // FILE *file_logerr,
        // SGZip *file_logerr_gz,
        char **matrix_pol,
        long int **matrix_freq,
        long int **matrix_pos,
        long int *length,
        long int *length_seg,
        double *length_al,
        long int *length_al_real,
        // int mainargc,
        // int *vint_perpop_nsam,
        // int npops,
        double *svratio,
        double *missratio,
        // int include_unknown,
        double *sum_sam,
        double **tcga,
        long int **matrix_sv,
        long int *nmhits,
        // int output,
        // int outgroup_presence,
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
        // long int first_slide,
        // long int slide,/**/
        // long int window,/**/
        // int Physical_length,
        long int *li /**/,
        int *npriors,
        double **vector_priors,
        char **matrix_pol_tcga,
        char *chr_name,
        int first,
        mstatspop_args_t *args);

#ifdef __cplusplus
}
#endif

#endif /* defined(__xcode_project__get_tfadata__) */
