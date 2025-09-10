/*
 *  mstatspop.c
 *  mstatspop
 *
 *  First version created by Sebastian Ramos-Onsins on 28/07/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "mstatspop.h"
#include "log.h"
#include "tfasta.h"

// // possible formats for the input file
// #define FASTA_FORMAT 0
// #define NBRF_FORMAT 0
// #define MS_FORMAT 1
// #define MS_X_FORMAT 2
// #define TFA_FORMAT 3

// /* o Output type, from 0 to 11 - TODO: define*/

// //  0 (extended),
// #define OUTPUT_EXTENDED 0

// //  1 (single line/window)
// #define OUTPUT_SINGLE_LINE_WINDOW 1

// //  2 (single line SFS/window)
// #define OUTPUT_SINGLE_LINE_SFS_WINDOW 2

// //  3 (dadi-like format)
// #define OUTPUT_DADI_LIKE 3

// //  4 (single line pairwise distribution)
// #define OUTPUT_SINGLE_LINE_PAIRWISE_DISTRIBUTION 4

// //  5 (single line freq. variant per line/window)
// #define OUTPUT_SINGLE_LINE_FREQ_VARIANT_PER_LINE_WINDOW 5

// //  6 (SNP genotype matrix)
// #define OUTPUT_SNP_GENOTYPE_MATRIX 6

// //  7 (SweepFiinder format -only first pop-)
// #define OUTPUT_SWEEP_FINDER 7

// //  8 (single line/window: Frequency of each haplotype in the populations)
// #define OUTPUT_SINGLE_LINE_FREQ_HAPLOTYPE_POP 8

// //  9 (single line/window: Frequency of variants per line and population)
// #define OUTPUT_SINGLE_LINE_FREQ_VARIANT_POP 9

// //  92 (single line/window: -rSFS- Frequency of variants per population relative to all)
// #define OUTPUT_SINGLE_LINE_FREQ_VARIANT_POP_RELATIVE 92

// //  10 (full extended)
// #define OUTPUT_FULL_EXTENDED 10

// typedef struct
// {
//   const char *format_arg; // format of the input file [fasta,ms,ms_x,tfa]
//   int formatfile;         // format type of the input file [0,1,2,3]
//   int output;             /* TODO: tipo de output, de 0 a 11 */
//   char ploidy[2];         // ploidy, 1: haploid, 2: diploid
//   long int niterdata;     // number of iterations (ms) or permutations (others)
//   char file_in[MSP_MAX_FILENAME] /* input file */;
//   char file_out[MSP_MAX_FILENAME]; /* output file */
//   char file_mas[MSP_MAX_FILENAME]; /* mask file */
//   char file_GFF[MSP_MAX_FILENAME]; /* GFF file */

//   /* a Name of file of the Alternative frequency spectrum */
//   /* Only with optimal tests (in GSL libs) */
//   char file_H1f[MSP_MAX_FILENAME];

//   /* a Name of file of the NULL frequency spectrum */
//   /* Only with optimal tests (in GSL libs) */
//   char file_H0f[MSP_MAX_FILENAME];
//   /* 1: missing values allowed , 0: excluded missing positions */
//   int include_unknown;

//   /* variables defining more data*/
//   int npops; // TODO args.npops
//   /* Number of samples for each population, each element a population */
//   int *vint_perpop_nsam; /* old nsam */

//   /* Sum of all nsam of all populations */
//   int int_total_nsam; /* old nsamtot*/

//   int coordfile;

//   /* number of iterations (ms) or permutations (others) */
//   long int niter;

//   /* GFF variables */
//   int gfffiles;
//   /*int 	observed_data	= 0;*/
//   char subset_positions[MSP_MAX_GFF_WORDLEN];
//   char code_name[MSP_GENETIC_CODETYPE_LEN];
//   char genetic_code[MSP_GENCODE_COMBINATIONS + 1];
//   char criteria_transcript[MSP_GFF_CRITERIA_MSG_LEN];

//   long int length; /* Sequence length */
//   int kind_length; // TODO args.kind_length

//   float ms_svratio;        // TODO args.ms_svratio
//   int outgroup_presence;   // TODO args.outgroup_presence
//   int force_outgroup;      // TODO args.force_outgroup
//   float freq_revert;       // TODO args.freq_revert
//   double freq_missing_ms;  // TODO args.freq_missing_ms
//   int location_missing_ms; // TODO args.location_missing_ms

//   /*R2_p*/
//   int *r2i_ploidies; // TODO args.r2i_ploidies
//   /*Covariances in missing values*/
//   int n_ccov;                                              // TODO args.n_ccov
//   int int_total_nsam_order;                                /*number of samples. value for ordering the samples in populations in case they are not consecutive*/
//   int *sort_nsam; /*vector with the position in the file*/ // TODO : args.sort_nsam
//   long int slide;                                          // TODO args.slide
//   long int first_slide;                                    // TODO args.first_slide
//   long int window;                                         // TODO args.window
//   char file_wps[MSP_MAX_FILENAME];                         // TODO args.file_wps
//   char file_Wcoord[MSP_MAX_FILENAME];                      // TODO args.file_Wcoord
//   int Physical_length;                                     // TODO args.Physical_length
//   int mask_print;                                          // TODO args.mask_print
//   char file_chr_name_all[MSP_MAX_NAME];                    // TODO :: args.file_chr_name_all
//   long int nseed;

// } mstatspop_args_t;

// function to init  mstatspop_args_t with default values
static inline void init_args(mstatspop_args_t *args)
{
    // set all values to 0
    memset(args, 0, sizeof(mstatspop_args_t));
    /* variables defining more data*/
    args->npops = 0;

    /* Number of samples for each population, each element a population */
    args->vint_perpop_nsam = NULL; /* old nsam */
    /* S->m of all nsam of all populations */
    args->int_total_nsam = 0; /* old nsamtot*/
    /* n->mber of iterations (ms) or permutations (others) */
    args->niter = 0;

    // i->t include_unknown;
    args->output = OUTPUT_EXTENDED; /* TODO: tipo de output, de 0 a 11 */
    args->formatfile = 0;

    memset(args->file_in, 0, MSP_MAX_FILENAME);
    memset(args->file_out, 0, MSP_MAX_FILENAME);
    memset(args->file_mas, 0, MSP_MAX_FILENAME);
    memset(args->file_GFF, 0, MSP_MAX_FILENAME);
    memset(args->file_H1f, 0, MSP_MAX_FILENAME);
    memset(args->file_H0f, 0, MSP_MAX_FILENAME);

    /* GFF variables */
    args->gfffiles = 0; // TODO args.gfffiles

    args->length = 0; /* Sequence length */
    args->kind_length = 0;

    args->nseed = 1234;

    // int outgroup_presence;  // TODO args.outgroup_presence
    // int force_outgroup; // TODO args.force_outgroup
    args->freq_revert = 0.;        // TODO args.freq_revert
    args->freq_missing_ms = 0.;    // TODO args.freq_missing_ms
    args->location_missing_ms = 0; // TODO args.location_missing_ms

    /*R2_p*/
    args->r2i_ploidies = 0; // TODO args.r2i_ploidies

    /*Covariances in missing values*/
    args->n_ccov = 0; // TODO args.n_ccov

    // TODO args.int_total_nsam_order
    args->int_total_nsam_order = 0;                               /*number of samples. value for ordering the samples in populations in case they are not consecutive*/
    args->sort_nsam = 0; /*vector with the position in the file*/ // TODO : args.sort_nsam

    /* ms */
    args->ms_svratio = 0.5; // TODO args.ms_svratio
    args->slide = 0;        // TODO args.slide
    args->first_slide = 0;  // TODO args.first_slide
    args->window = 0;       // TODO args.window
    // char file_wps[MSP_MAX_FILENAME]; // TODO args.file_wps
    // char file_Wcoord[MSP_MAX_FILENAME]; // TODO args.file_Wcoord
    args->Physical_length = 1; // TODO args.Physical_length

    memset(args->subset_positions, 0, MSP_MAX_GFF_WORDLEN);
    memset(args->code_name, 0, MSP_GENETIC_CODETYPE_LEN);
    memset(args->genetic_code, 0, MSP_GENCODE_COMBINATIONS + 1);
    memset(args->criteria_transcript, 0, MSP_GFF_CRITERIA_MSG_LEN);

    memset(args->file_wps, 0, MSP_MAX_FILENAME);
    memset(args->file_Wcoord, 0, MSP_MAX_FILENAME);

    /***** DEFAULTS *****/
    // refactoring: move all args to a struct

    args->formatfile = TFA_FORMAT;            /* 3 tfasta format*/
    args->output = OUTPUT_SINGLE_LINE_WINDOW; /* 1 single line statistics per window*/
    args->length = 0;                         /*undefined*/
    /*------------------*/
    args->npops = 0;                /*undefined*/
    args->outgroup_presence = 0;    /*no outgroup*/
    args->include_unknown = 0;      /*excluded missing positions*/
    args->include_rsfs = 0;         /*analysis using the relative site frequency spectrum (rSFS)*/
    //args->include_fstall = 0;       /*if 0, analysis of Fst if less or equal than 6 pops*/
    args->int_total_nsam_order = 0; /*undefined*/
    /*r2i_ploidies=0;*/             /*already defined*/
    args->window = 0;
    args->slide = 0;
    args->Physical_length = 1;
    args->niter = 0;
    args->nseed = 123456; // TODO args.nseed
    args->niterdata = 1;
    strcpy(args->file_mas, "-1\0"); /* -1 = No file provided*/
    args->ms_svratio = 0.5;
    args->force_outgroup = 0;
    args->freq_revert = 0;
    args->ploidy[0] = '1'; /* 1:haploid, 2:diploid, and so on...NO, next Ver */
    args->ploidy[1] = '\0';
    strcpy(args->criteria_transcript, "long\0");
    args->coordfile = 0;
    args->mask_print = 0;
    args->first_slide = 0;
}

/**
 * \brief main function
 * \details More details
 * @param argc
 * @param argv
 * @return
 */
/* [input_file] [#_pops] [#samples_pop1] ... [#samples_popN] [GFF_file] */
int main(int argc, const char *argv[])
{

    // struct to hold all command line arguments
    mstatspop_args_t args;
    // set default values
    init_args(&args);
    args.argc = argc;

    /*Statistics*/
    struct stats *statistics = 0;
    struct stats *stats_iter = 0;
    struct probs *piter = 0;
    struct rweights *rw = 0;
    
    int w, x, y, z, yy /*,ii*/;
    long int i;
    char *f;
    // char ploidy[2];

    /* Population index */
    int pop_index = 0;

    /* name and data files */
    // char file_in[ MSP_MAX_FILENAME];
    // char file_out[MSP_MAX_FILENAME];
    // char file_mas[MSP_MAX_FILENAME];
    // char file_GFF[MSP_MAX_FILENAME];
    // char file_H1f[MSP_MAX_FILENAME];
    // char file_H0f[MSP_MAX_FILENAME];
    char file_log[MSP_MAX_FILENAME];

    // currently used with ms and fasta formats
    FILE *file_input = 0;
    SGZip file_input_gz;

    FILE *error_log_file = 0; // file_logerr

    FILE *file_output = stdout;
    SGZip file_output_gz;

    FILE *file_mask = 0;
    FILE *file_H1freq = 0;
    FILE *file_H0freq = 0;

    // SGZip file_logerr_gz;

    /*int 	observed_data	= 0;*/
    // char 	subset_positions[ MSP_MAX_GFF_WORDLEN ];	 // TODO args.subset_positions
    // char 	code_name[ MSP_GENETIC_CODETYPE_LEN ];  // TODO args.code_name
    // char 	genetic_code[ MSP_GENCODE_COMBINATIONS +1]; // TODO args.genetic_code
    // char 	criteria_transcript[ MSP_GFF_CRITERIA_MSG_LEN ]; // TODO args.criteria_transcript

    /* Alignment data variables */
    char *matrix_pol = 0;         /*the matrix of SNPs alignment data*/
    char *matrix_pol_tcga = 0;    /*the matrix of SNPs alignment data with four nucleotides*/
    long int *matrix_freq = 0;    /*frequency of the segregating site*/
    long int *matrix_pos = 0;     /*position of the segregating site*/
    /*long int *matrix_GC = 0; */ /*vector with 1 if GC exists in these position and the contiguous  */
    long int *matrix_sv = 0;      /*vector with x  for each position, where x */
    /* equals: 1 transition, 2 transversion, 0 nothing
     * (biallelic only, no multiple hits considered) */
    long int *sites_matrix = 0;
    long int length_seg = 0;
    double length_al = 0;
    long int length_al_real;

    double svratio;
    double missratio = 0.;

    double *sum_sam; /*length of each sequence excluding non-tcga*/
    double **tcga;
    long int nmhits = 0;
    int *matrix_mask = 0;   /*mask matrix for ms inputs*/
    float *vector_mask = 0; /*value of each position: eg. in case synonymous, the position can count less than 1*/

    /* joint frequency distribution */
    double **jfd; /*frequency of each segregating site per population*/
    int **nfd;    /*number of samples for each segregating site and population*/

    /* permutation test */
    char *matrix_perm = 0; /*the matrix of SNPs alignment data permutated*/
    int nsam2[3];          /*the sample size for each pop in permutation test. Include outg*/
    int psam2[3];          /*the position of the first sample in matrix_pol of each pop for perm. test. Include outg*/

    /* Optimal tests */
    int H1frq, H0frq;
    int npf = 0;
    double **freqspH1 = 0;
    double **freqspH0 = 0;
    double *thetaH1 = 0;
    double *thetaH0 = 0;
    char *cad, *cad1;
    int len;

    // long int niterdata;
    long int li, li2;
    int flaghky;

    /* read mask_file */
    char c[2], vvm[11];
    int vli;
    long int n;

    double *sum_sam_mask = 0;
    double length_mask = 0;
    long int length_mask_real = 0;
    int maxnsam;

    int arg;

    char *el;
    double *vector_priors = 0;
    int npriors = 0;

    //pop analysis
    double *nsites1_pop;
    double *nsites2_pop;
    double *nsites3_pop;
    double *nsites1_pop_outg;
    double *nsites2_pop_outg;
    double *nsites3_pop_outg;
    double *anx, *bnx;
    double *anxo, *bnxo;
    double **lengthamng;
    double **lengthamng_outg;
    //rSFS analysis
    double *rnsites1_pop=0;
    double *rnsites2_pop=0;
    double *rnsites3_pop=0;
    double *rnsites1_pop_outg=0;
    double *rnsites2_pop_outg=0;
    double *rnsites3_pop_outg=0;

    int sort_index; /*used in the loop*/

    /*tfasta windows and weights*/
    long int *wgenes;  /*init and end coordinates*/
    long int nwindows; /*number of fragments*/

    // TODO :: weights files need to be converted and indexed same was as tfasta file. DONE
    // FILE *file_ws = 0;
    //SGZip file_ws_gz;

    FILE *file_wcoor = 0;
    SGZip file_wcoor_gz;

    // struct SGZIndex index_input;
    struct SGZIndex index_w;

    // int mask_print; // TODO args.mask_print

    unsigned long nscaffolds;
    // char file_chr_name_all[ MSP_MAX_NAME]; // TODO :: args.file_chr_name_all
    char *chr_name = 0;
    char chr_name_all[MSP_MAX_NAME];
    char **chr_name_array;
    char **chr_length_array;
    int first = 0;

    memset(chr_name_all, 0, MSP_MAX_NAME);

    memset(c, 0, 2);

    /*Default values*/
    /*file_out            = stdout;*/ /*already defined as stdout*/
    H1frq = 0;
    H0frq = 0;

    /********************/

    /* STEP 1: PARSING COMMAND LINE ------------------------------------------ */

    if (argc <= 1)
    {
        // No arguments given, print help and exit with error
        usage();
        exit(1);
    }

    // set program name as mstatspop
    const char *program_name = "mstatspop";
#ifdef DEBUG
    log_set_level(LOG_TRACE);
#else
    log_set_level(LOG_INFO);
#endif
    log_trace("Parsing command line arguments");

    log_start(program_name, argc, argv);

    arg = 1;
    while (arg < argc)
    {
        if (argv[arg][0] != '-')
        {
            if (argv[arg][0] == '>')
                break;
            // printf(" argument should be -%s ?\n", argv[arg]);
            log_error("argument should be -%s ?", argv[arg]);
            usage();
            exit(1);
        }

        switch (argv[arg][1])
        {
            case 'f': /* f FORMAT */
                arg++;

                args.format_arg = argv[arg];
                /* if format is not recognized, error and out */
                if (strcmp(argv[arg], "fasta") != 0 &&
                    strcmp(argv[arg], "nbrf") != 0 &&
                    strcmp(argv[arg], "ms") != 0 &&
                    strcmp(argv[arg], "tfa") != 0)
                {
                    // printf("\n Error: the argument -f has only the choices 'fasta', 'tfa' or 'ms'.");
                    log_error("Error: the argument -f has only the choices 'fasta', 'tfa' or 'ms'.");
                    usage();
                    exit(1);
                }
                else if (strcmp(argv[arg],
                                "fasta") == 0 || strcmp(argv[arg],
                                                        "nbrf") == 0)
                {
                    args.formatfile = FASTA_FORMAT; // 0;
                    args.niterdata = 1;
                }
                else if (strcmp(argv[arg], "ms") == 0)
                {
                    args.formatfile = MS_FORMAT; // 1;
                }
                else if (strcmp(argv[arg], "ms_x") == 0)
                {
                    args.formatfile = MS_X_FORMAT; //  2;
                }
                else if (strcmp(argv[arg], "tfa") == 0)
                {
                    args.formatfile = TFA_FORMAT; //  3;
                    args.niterdata = 1;
                }
                break;

            case 'i': /* i Input File, el path */
                arg++;
                strcpy(args.file_in, argv[arg]);
                break;

            case 'T': /* i output File, el path */
                arg++;
                strcpy(args.file_out, argv[arg]);
                if ((file_output = fopen(args.file_out, "w")) == 0)
                { /*zipped not available yet*/
                    // fprintf(stdout,"\n It is not possible to write in the output file %s\n", file_out);
                    log_error("It is not possible to write in the output file %s\n",
                              args.file_out);
                    exit(1);
                }
                strcpy(file_log, args.file_out);
                strcat(file_log, ".log");
                if ((error_log_file = fopen(file_log, "w")) == 0)
                {
                    // fprintf(stdout,"\n It is not possible to write the log file %s.", file_log);
                    log_error("It is not possible to write the log file %s.",
                              file_log);
                    exit(1);
                }
                // TODO :: set log level from command line
                log_add_fp(error_log_file, LOG_DEBUG);
                break;
            case 'o': /* o Output type, from 0 to 11 - TODO: define*/
                arg++;
                args.output = (int)atoi(argv[arg]);
                // TODO :: validate output type of any defined value
                if (args.output < 0 || args.output > 100)
                {
                    // printf("\n Error in -o argument: only values between 0 and 10 are allowed.");
                    log_error("Error in -o argument: only values between 0 and 10 are allowed.");
                    usage();
                    exit(1);
                }
                break;

            case 'p': /* p Ploidy, 1: haploid, 2: diploid */
                arg++;
                args.ploidy[0] = argv[arg][0];
                if (args.ploidy[0] != '1' && args.ploidy[0] != '2')
                {
                    // printf("\n Error in -p argument: only the values 1 or 2 are allowed.");
                    log_error("Error in -p argument: only the values 1 or 2 are allowed.");
                    usage();
                    exit(1);
                }
                break;

            case 'u': /* u Missing Values Allowed or not, 1: allowed */
                arg++;
                args.include_unknown = (int)atoi(argv[arg]);

                if (args.include_unknown != 0 && args.include_unknown != 1)
                {
                    // printf("\n Error in -u argument: only the values 0 or 1 are allowed.");
                    log_error("Error in -u argument: only the values 0 or 1 are allowed.");
                    usage();
                    exit(1);
                }
                break;

            case 'R': /* R performs rSFS analysis. (default) 0: do not perform. 1: performs */
                arg++;
                args.include_rsfs = (int)atoi(argv[arg]);
            
                if (args.include_rsfs != 0 && args.include_rsfs != 1)
                {
                    // printf("\n Error in -R argument: only the values 0 or 1 are allowed.");
                    log_error("Error in -R argument: only the values 0 or 1 are allowed.");
                    usage();
                    exit(1);
                }
                break;
            /*
            case 'D': // S Calculate Fst if more than 6 pop. (default) 0: do calculate,  1: calculate
                arg++;
                args.include_fstall = (int)atoi(argv[arg]);
            
                if (args.include_fstall != 0 && args.include_fstall != 1)
                {
                    // printf("\n Error in -F argument: only the values 0 or 1 are allowed.");
                    log_error("Error in -F argument: only the values 0 or 1 are allowed.");
                    usage();
                    exit(1);
                }
                break;
            */
            case 'm': /* m Mask file, if -1, no file included */
                arg++;
                strcpy(args.file_mas, argv[arg]);
                break;

            case 's': /* s Seed, positive value */
                arg++;
                args.nseed = (long int)atol(argv[arg]);
                break;

            case 't': /* t number of permutations, in case of fasta format */
                arg++;
                args.niter = (long int)atol(argv[arg]);
                break;

            case 'r': /* number of iterations, in case of ms */ /*NO or tfasta */
                arg++;
                args.niterdata = (long int)atol(argv[arg]);
                break;

            case 'v': /* v ratio transitions/transversions, only in ms */
                arg++;
                args.ms_svratio = (double)atof(argv[arg]);
                break;

            case 'l': /* l length of the sequence, only for ms format */
                arg++;
                args.length = (long int)atol(argv[arg]); // actual usage might be different, as we pass it around
                break;

            case 'k': /* kind of l length of the sequence (0,1,2,3), only for ms format */
                arg++;
                args.kind_length = (int)atoi(argv[arg]);
                break;

            case 'G': /* G outgroup: included=1, non included=0 */
                arg++;
                args.outgroup_presence = (int)atoi(argv[arg]);
                if (args.outgroup_presence != 0 && args.outgroup_presence != 1)
                {
                    // printf("\n Error in -G argument: only the values 0 or 1 are allowed.");
                    log_error("Error in -G argument: only the values 0 or 1 are allowed.");
                    usage();
                    exit(1);
                }
                break;
            case 'F': /* F forced_outgroup=1 */
                arg++;
                args.force_outgroup = (int)atoi(argv[arg]);
                break;

            case 'q': /* frequency of revert mutation when forcing outgroup, only in ms format */
                arg++;
                args.freq_revert = (double)atof(argv[arg]);
                break;

            case 'N': /* N number of populations, Warning!, followed by
                       N numbers indicating the sample size of each population
                       */
                arg++;
                args.npops = atoi(argv[arg]);
                if ((args.vint_perpop_nsam = (int *)calloc((unsigned long)args.npops,
                                                           sizeof(int))) == 0)
                {
                    // printf("Error allocating memory");
                    log_fatal("Error allocating memory, vint_perpop_nsam");
                    exit(1);
                }
                args.int_total_nsam = 0;
                for (pop_index = 0; pop_index < args.npops; pop_index++)
                {
                    arg++;
                    args.vint_perpop_nsam[pop_index] = atoi(argv[arg]);
                    args.int_total_nsam += args.vint_perpop_nsam[pop_index];
                }
                break;

            case 'O': /* O the order of each individual in the original data, Warning!, followed by
                       O numbers indicating the order (0 is the first) in case samples are not consecutive. Only for fasta data!
                       */
                arg++;
                args.int_total_nsam_order = atoi(argv[arg]);
                if ((args.sort_nsam = (int *)calloc((unsigned long)args.int_total_nsam_order,
                                                    sizeof(int))) == 0)
                {
                    // printf("Error allocating memory");
                    log_fatal("Error allocating memory, sort_nsam");
                    exit(1);
                }
                for (sort_index = 0; sort_index < args.int_total_nsam_order; sort_index++)
                {
                    arg++;
                    args.sort_nsam[sort_index] = atoi(argv[arg]);
                }
                break;

            case 'g': /* g GFF file name, AND more words
                       2nd : synonymous, nonsynonymous, silent or whatever
                       3rd : selected genetic code name or "Other"
                       next 64th's : in case of 'Other', provide 64 IUPAC letters of each codon.
                       * Check order.
                       */
                arg++;
                strcpy(args.file_GFF, argv[arg]);
                arg++;
                strcpy(args.subset_positions, argv[arg]);

                args.gfffiles = 1;

                /* Go if known coding option - */
                if ((strcmp(args.subset_positions, "synonymous") == 0 ||
                     strcmp(args.subset_positions, "nonsynonymous") == 0 ||
                     strcmp(args.subset_positions, "0-fold") == 0 ||
                     strcmp(args.subset_positions, "2-fold") == 0 ||
                     strcmp(args.subset_positions, "3-fold") == 0 ||
                     strcmp(args.subset_positions, "4-fold") == 0 ||
                     strcmp(args.subset_positions, "silent") == 0))
                {
                    arg++;
                    strcpy(args.code_name, argv[arg]);

                    if (strcmp(args.code_name, "Other") == 0)
                    {
                        for (x = 0; x < 64; x++)
                        {
                            arg++;
                            if (argv[arg][0] == '-')
                            {
                                // printf("\n Error in -g argument: In case use \"Other\", include the genetic code of the 64 aa values.");
                                log_error("Error in -g argument: In case use \"Other\", include the genetic code of the 64 aa values.");
                                usage();
                                exit(1);
                            }
                            args.genetic_code[x] = atoi(argv[arg]);
                        }
                    }
                }
                break;

            case 'A': /* a Name of file of the Alternative frequency spectrum */
                /* Only with optimal tests (in GSL libs) */
                arg++;
                strcpy(args.file_H1f, argv[arg]);
                H1frq = 1;
                break;

            case 'S': /* a Name of file of the NULL frequency spectrum */
                /* Only with optimal tests (in GSL libs) */
                arg++;
                strcpy(args.file_H0f, argv[arg]);
                H0frq = 1;
                break;

            case 'c': /* c Criteria used for analyzing the transcripts */
                /* Basically, max or min */
                arg++;
                strcpy(args.criteria_transcript, argv[arg]);
                if (strcmp(args.criteria_transcript, "max") != 0 &&
                    strcmp(args.criteria_transcript, "min") != 0 &&
                    strcmp(args.criteria_transcript, "first") != 0 &&
                    strcmp(args.criteria_transcript, "long") != 0)
                {
                    // printf("\n Error: the argument -c has only the choices 'max', 'min', 'first' or 'long'.");
                    log_error("Error: the argument -c has only the choices 'max', 'min', 'first' or 'long'.");
                    usage();
                    exit(1);
                }
                break;

            case 'x': /* proportion of missing values , only in ms format with the option -u 1 */
                arg++;
                args.freq_missing_ms = (double)atof(argv[arg]);
                break;

            case 'y': /* number of comparisons for the calculation of covariances in neutrality tests
                       using missing values (0 if not desired) */
                arg++;
                args.n_ccov = (int)atol(argv[arg]);
                break;

            case 'M': /* in case ms_e, column location of missing values at prior row (default 3), only in ms_e format with the option -u 1 */
                arg++;
                args.location_missing_ms = (double)atoi(argv[arg]);
                break;

            case 'P': /*Calculation of R2_p: first value is the number of values to include, next are the ploidies to consider. ex: -P 6 1 2 4 8 16 64*/
                arg++;
                if ((args.r2i_ploidies = (int *)malloc((unsigned long)(atoi(argv[arg]) + 1) * sizeof(int))) == NULL)
                {
                    // printf("\nError: memory not reallocated. mstatspop.c.00 \n");
                    log_fatal("Error: memory not reallocated, r2i_ploidies mstatspop.c.00");
                    exit(1);
                }
                x = 0;
                args.r2i_ploidies[x] = atoi(argv[arg]);
                while (x < args.r2i_ploidies[0])
                {
                    x++;
                    arg++;
                    args.r2i_ploidies[x] = atoi(argv[arg]);
                }
                break;
            case 'Y': /* physical length or effective length (only valid positions) */
                arg++;
                args.Physical_length = (int)atoi(argv[arg]);
                if (args.Physical_length != 0 && args.Physical_length != 1)
                {
                    // printf("\n Error in -l argument: only the values 0 or 1 are allowed.");
                    log_error("Error in -l argument: only the values 0 or 1 are allowed.");
                    usage();
                    exit(1);
                }
                break;
            case 'w': /* window size */
                arg++;
                args.window = (long int)atol(argv[arg]);
                break;
            case 'z': /* slide size */
                arg++;
                args.slide = (long int)atol(argv[arg]);
                break;
            case 'W': /* file with the coordinates of each window [init end](overwrite options -w and -s)*/
                arg++;
                strcpy(args.file_Wcoord, argv[arg]);
                args.coordfile = 1;
                break;
            case 'E': /*file with the weight for each position */
                arg++;
                strcpy(args.file_wps, argv[arg]);
                break;
            case 'K': /*output a mask file or not */
                arg++;
                args.mask_print = (int)atoi(argv[arg]);
                break;
            case 'Z': /* first slide size */
                arg++;
                args.first_slide = (long int)atol(argv[arg]);
                break;
            case 'n': /* name of the file with scaffold(s) and total lengths to analyze*/
                arg++;
                strcpy(args.file_chr_name_all, argv[arg]);
                // strcpy( chr_name_all, argv[arg] );
                break;
            case 'h': /* h HEEEEEEEEEEELPPPPPPPPPPPPPP!!! */
                usage();
                exit(0);
                break;
        }
        arg++;
    } // end while (arg < argc)



    // validate arguments
    // validate input tfasta file and index file
    tfasta_file tfasta;

    // data structure to hold tfasta weights file and index
    wtfasta_file *wtfasta;
    wtfasta = NULL;
    // allocat mem for wtfasta structure
 
    if (args.formatfile == TFA_FORMAT)
    {
        log_info("Input tfasta format : running tfasta mode.");
        log_info("Validating tfasta file %s", args.file_in);
        // validate tfasta file exists and is readable and in TFAv2.0 format
        if (access(args.file_in, R_OK) != 0)
        {
            // fprintf(file_logerr, "\nError: the file %s does not exist or is not readable.\n", file_in);
            log_error("Error: the file %s does not exist or is not readable.",
                      args.file_in);
            exit(1);
        }
        memset(&tfasta, 0, sizeof(tfasta));
        int ret_status = init_tfasta_file(&tfasta, args.file_in);
        if (ret_status == TFA_ERROR)
        {
            log_error("Can not initialize tfasta file %s", args.file_in);
            exit(1);
        }
        if(ret_status == TFA_INVALID_FORMAT)
        {
            log_error("The file %s is not in TFAv2.0 format", args.file_in);
            exit(1);
        }
        if(ret_status == TFA_INVALID_INDEX)
        {
            log_error("The index file for %s is invalid", args.file_in);
            exit(1);
        }
        // if we pass here, the file is valid and we can use it
        // validate other parameters for tfasta format
        log_info("tfasta file %s is valid", args.file_in);
        if (args.file_wps[0] != '\0')
        {
            log_info("Using tfasta weights file %s", args.file_wps);
            log_info("Validating tfasta weights file %s", args.file_wps);
            // validate weights file exists and is readable
            if (access(args.file_wps, R_OK) != 0)
            {
                log_error("Error: the tfasta weights file %s does not exist or is not readable.",
                          args.file_wps);
                exit(1);
            }

            if ((wtfasta = (wtfasta_file *)malloc(sizeof(wtfasta_file))) == NULL)
            {
                // fprintf(file_logerr, "\nError: memory not reallocated. mstatspop.c.00 \n");
                log_fatal("Error: can not allocate memory for wtfasta_file structure.");
                exit(1);
            }

            // validate weights file is in TFAv2.0 format
            // memset the structure to 0 
            memset(wtfasta, 0, sizeof(wtfasta_file));
            int ret_status = init_wtfasta_file(wtfasta, args.file_wps);
            if (ret_status == TFA_ERROR)
            {
                log_error("Can not initialize tfasta weights file %s",
                          args.file_wps);
                exit(1);
            }
            if(ret_status == TFA_INVALID_FORMAT)
            {
                log_error("The tfasta weights file %s is not in wTFAv2.0 format",
                          args.file_wps);
                exit(1);
            }
            if(ret_status == TFA_INVALID_INDEX)
            {
                log_error("The index file for %s is invalid", args.file_wps);
                exit(1);
            }
            log_info("tfasta weights file %s is valid", args.file_wps);
        }
    }

    /*few filters*/
    // if (!(args.formatfile == 1 || args.formatfile == 2) && args.force_outgroup == 1)
    if (!(args.formatfile == MS_FORMAT || args.formatfile == MS_X_FORMAT) && args.force_outgroup == 1)
    {
        // fprintf(file_logerr,"\nError. The option -F 1 is only compatible with -f 'ms'.\n");
        log_error("Error. The option -F 1 is only compatible with -f 'ms'.");
        exit(1);
    }
    if (args.npops == 0)
    {
        // fprintf(file_logerr,"\nError. The option -N must be included.\n");
        log_error("Error. The option -N must be included.");
        exit(1);
    }
    // if (args.window <= 0 && args.formatfile == 3 && args.coordfile == 0)
    if (args.window <= 0 && args.formatfile == TFA_FORMAT && args.coordfile == 0)
    {
        // fprintf(file_logerr,"\nError. The option -w or -W must be included with option -f tfa.\n");
        log_error("Error. The option -w or -W must be included with option -f tfa.");
        exit(1);
    }
    if (args.slide == 0 && args.window > 0 && args.coordfile == 0)
        args.slide = args.window;
    // if (args.slide <= 0 && args.formatfile == 3 && args.coordfile == 0)
    if (args.slide <= 0 && args.formatfile == TFA_FORMAT && args.coordfile == 0)
    {
        // fprintf(file_logerr,"\nError. The value at option -z (slide) must be larger than 0\n");
        log_error("Error. The value at option -z (slide) must be larger than 0");
        exit(1);
    }
    // if (args.first_slide < 0 && args.formatfile == 3)
    if (args.first_slide < 0 && args.formatfile == TFA_FORMAT)
    {
        // fprintf(file_logerr,"\nError. The value at option -Z (first slide) must be larger than or equal to 0.\n");
        log_error("Error. The value at option -Z (first slide) must be larger than or equal to 0.");
        exit(1);
    }
    /*
     if(length == 0 && formatfile == 1) {
     // fprintf(file_logerr,"\nError. length (-l option) must be defined with ms input file.\n");
     log_error("Error. length (-l option) must be defined with ms input file.");
     exit(1);
     }*/
    // if (args.formatfile == 0 && args.niterdata > 1)
    if (args.formatfile == FASTA_FORMAT && args.niterdata > 1)
    {
        // fprintf(file_logerr,"\nError. The option -f fasta does not accept the option -r\n");
        log_error("Error. The option -f fasta does not accept the option -r");
        exit(1);
    }
    if (args.formatfile == TFA_FORMAT && args.niterdata > 1)
    {
        // fprintf(file_logerr,"\nError. The option -f tfa does not accept the option -r\n");
        log_error("Error. The option -f tfa does not accept the option -r");
        exit(1);
    }
    if (args.include_unknown == 1)
        args.niter = 0;

    if (strcmp(args.file_chr_name_all, "") == 0)
    {
        // fprintf(file_logerr,"\nError: the file name containing the scaffold(s) (option -n) must be defined\n");
        log_error("Error: the file name containing the scaffold(s) (option -n) must be defined");
        exit(1);
    }
    /*
     if(file_Wcoord[0]!=0 && (slide > 0 && window>0)) {
     // fprintf(file_logerr,"\n the option -W (coordinates file) is incompatible with definitions of -w and -z ");
     log_error("the option -W (coordinates file) is incompatible with definitions of -w and -z ");
     exit(1);
     }
     */
    /* STEP 2: READING FILE INFO --------------------------------------------- */

    /* Opening files */
    // legacy code,
    if (args.formatfile != TFA_FORMAT)
    {
        if (args.file_in[0] == '\0')
        {
            // if no input file is provided, read from stdin
            file_input = stdin;
            file_input_gz.file_compressed = 0;
            // if (args.formatfile == 3)
            // if (args.formatfile == TFA_FORMAT)
            // {
            //   /*tfa must be a file because it needs the index file*/
            //   // fprintf(file_logerr,"\nError: tfa format file needs indexation. \nstdout is not available using this option. \n");
            //   log_error("Error: tfa format file needs indexation. \nstdout is not available using this option.");
            //   exit(1);
            // }
        }
        else
        {
            // if input file is provided, open it
            if ((file_input = fzopen(args.file_in, "r", &file_input_gz)) == 0)
            {
                // fprintf(file_logerr,"\n It is not possible to open the input file %s.", file_in);
                log_error("It is not possible to open the input file %s.",
                          args.file_in);
                exit(1);
            }
            // if tfasta format, load index file as well
            // if (args.formatfile == 3)
            // if (args.formatfile == TFA_FORMAT)
            // { /*tfasta file*/
            //   load_index_from_file(file_input_gz.index_file_name, &index_input);
            // }
        }

        if ((f = (char *)malloc((unsigned long)BUFSIZ * 10)) == NULL)
        {
            // fprintf(file_logerr,"\nError: memory not reallocated. main.4 \n");
            log_fatal("Error: memory not reallocated, main.4");
            exit(1);
        }
        /* Definition of a File Stream Buffer, for buffered IO */
        setbuf(file_input, f);
    }

    /*separate all values of the list chr_name_all in chr_name_array: */
    if (read_index_file(args.file_chr_name_all,
                        &nscaffolds,
                        &chr_name_array,
                        &chr_length_array))
    {
        // printf("Error reading the scaffold names file %s\n",file_chr_name_all);
        log_error("Error reading the scaffold names file %s\n",
                  args.file_chr_name_all);
        exit(1);
    }

    // if (args.formatfile != 3 && nscaffolds > 1)
    if (args.formatfile != TFA_FORMAT && nscaffolds > 1)
    {
        // printf("Error: it only possible to read one scaffold with option -f fasta or -f ms.\n");
        log_error("Error: it only possible to read one scaffold with option -f fasta or -f ms.");
        exit(1);
    }
    if (args.length == 0)
    {
        args.length = atol(chr_length_array[0]); /*for ms files*/
    }
    /*separate all values of the list chr_name_all in chr_name_array: */
    /* Only do the list if input and output is tfa*/
    /*
     nscaffolds = 1;
     if(formatfile == 3) {
     chr_name_array = (char **)calloc(nscaffolds,sizeof(char *));
     chr_name_array[0] = (char *)calloc(1,sizeof(MSP_MAX_NAME));
     j=0;
     while(chr_name_all[j] != '\0') {
     k=0;
     while(chr_name_all[j] != ',' && chr_name_all[j] != '\0' && j < MSP_MAX_NAME) {
     chr_name_array[nscaffolds-1][k] = chr_name_all[j];
     j++; k++;
     }
     if(chr_name_all[j] == ',') {
     nscaffolds += 1;
     chr_name_array = (char **)realloc(chr_name_array,nscaffolds*sizeof(char *));
     chr_name_array[nscaffolds-1] = (char *)calloc(1,sizeof(MSP_MAX_NAME));
     j++;
     }
     }
     }
     else {
     chr_name_array = (char **)calloc(1,sizeof(char *));
     chr_name_array[0] = (char *)calloc(1,sizeof(MSP_MAX_NAME));
     strcpy(chr_name_array[0],chr_name_all);
     }
     */
    /*if( (formatfile == 1 || formatfile == 2) && output == 0 && niter > 1) output = 1;*/
    /*if( (formatfile == 1 || formatfile == 2)) length_al = length;*/
    if (args.formatfile == FASTA_FORMAT)
        strcpy(args.file_mas, args.file_in);
    if (args.formatfile == TFA_FORMAT && ((args.slide == 0 && args.window == 0) && args.file_Wcoord[0] == '\0'))
        strcpy(args.file_mas, args.file_in);

    /* TODO: por eso la suma nunca es 2, o 0 o 1... */
    if (args.outgroup_presence == 1)
    {
        args.force_outgroup = 0;
    }
    else if (args.outgroup_presence == 0)
    {
        /* Al cargar los datos, -N, hemos rellenado el vector */
        /* TODO: Check for errors in realloc */
        if ((args.vint_perpop_nsam = (int *)realloc(args.vint_perpop_nsam,
                                                    (args.npops + 1) * sizeof(int))) == 0)
        {
            // fprintf(file_logerr,"Error allocating memory");
            log_fatal("Error allocating memory, vint_perpop_nsam");
            exit(1);
        }

        args.vint_perpop_nsam[args.npops] = 1;
        args.int_total_nsam += 1;
        args.npops++;
    }

    if (args.force_outgroup == 1)
        args.outgroup_presence = 0;
    if (args.include_unknown == 1 && args.ploidy[0] == '2' && args.n_ccov > 0)
        args.n_ccov = 0; /*warning?*/

    /* -------------------------------------------------------------------- */
    /* -a alternative frequency spectrum: sum=S */

    if (H1frq == 1 && args.include_rsfs == 0)
    {
        npf = args.npops - 1; /* TODO: Mirar de cambiar por npops - 1...*/

        if ((file_H1freq = fopen(args.file_H1f, "r")) == 0)
        {
            // fprintf(file_logerr,"\n Error opening Alternative frequency spectrum file %s",file_H1f);
            log_error("Error opening Alternative frequency spectrum file %s",
                      args.file_H1f);
            exit(1);
        }

        /* Solo lee la primera linea, y un máximo de caracteres */
        if (!feof(file_H1freq))
        {
            /* TODO: check for errors */
            /* TODO: Cambiar para leer longitudes indeterminadas, no MAX_LEN
             */
            /* antes 1026 y 1024 */
            if ((cad = (char *)calloc(MSP_MAX_FILELINE_LEN, sizeof(char))) == 0)
            {
                // fprintf(file_logerr,"Error allocating memory");
                log_fatal("Error allocating memory, cad");
                exit(1);
            }

            fgets(cad, MSP_MAX_FILELINE_LEN, file_H1freq);
        }
        else
        {
            // fprintf(file_logerr,"\n Error reading Alternative frequency spectrum file %s", file_H1f);
            log_error("Error reading Alternative frequency spectrum file %s",
                      args.file_H1f);
            exit(1);
        }

        if ((freqspH1 = (double **)calloc((unsigned long)npf,
                                          sizeof(double *))) == 0)
        {
            // fprintf(file_logerr,"Error allocating memory");
            log_fatal("Error allocating memory, freqspH1");
            exit(1);
        }

        if ((freqspH0 = (double **)calloc((unsigned long)npf,
                                          sizeof(double *))) == 0)
        {
            // fprintf(file_logerr,"Error allocating memory");
            log_fatal("Error allocating memory, freqspH0");
            exit(1);
        }

        if ((thetaH1 = (double *)calloc((unsigned long)npf,
                                        sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"Error allocating memory");
            log_fatal("Error allocating memory, thetaH1");
            exit(1);
        }

        if ((thetaH0 = (double *)calloc((unsigned long)npf,
                                        sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"Error allocating memory");
            log_fatal("Error allocating memory, thetaH0");
            exit(1);
        }

        for (x = 0; x < npf; x++)
        {
            if ((freqspH1[x] = (double *)calloc((unsigned long)args.vint_perpop_nsam[x],sizeof(double))) == 0)
            {
                // fprintf(file_logerr,"Error allocating memory");
                log_fatal("Error allocating memory, freqspH1[x]");
                exit(1);
            }

            if ((freqspH0[x] = (double *)calloc((unsigned long)args.vint_perpop_nsam[x],sizeof(double))) == 0)
            {
                // fprintf(file_logerr,"Error allocating memory");
                log_fatal("Error allocating memory, freqspH0[x]");
                exit(1);
            }

            if (!feof(file_H1freq))
            {
                fgets(cad, MSP_MAX_FILELINE_LEN, file_H1freq);
                cad1 = cad; /* TODO: REVISAR, copiando punteros! */
            }
            else
            {
                // fprintf(file_logerr,"\n  Error reading Alternative frequency spectrum file %s, line %d",file_H1f,x+2);
                log_error("Error reading Alternative frequency spectrum file %s, line %d",
                          args.file_H1f,
                          x + 2);
                exit(1);
            }

            for (y = 1; y < args.vint_perpop_nsam[x]; y++)
            {
                if ((len = (int)strlen(cad1)) > 0)
                {
                    freqspH1[x][y] = (double)atof(cad1);
                    cad1 = strchr(cad1, '\t');
                    while (cad1 && (*cad1 == ' ' || *cad1 == '\t'))
                    {
                        cad1 += 1;
                    }
                }
            }
            thetaH1[x] = (double)atof(cad1);
        }
        free(cad);
        fclose(file_H1freq);
    }
    else
    {
        H1frq = 0;
    }

    /* -------------------------------------------------------------------- */
    /* -n null frequency spectrum: sum=S */
    /* TODO: Es igual que la anterior !!!!!!!!
     * REFACTOR IT!!!!
     * */
    if ((H0frq == 1 && H1frq == 1) && args.include_rsfs == 0)
    {
        if ((file_H0freq = fopen(args.file_H0f, "r")) == 0)
        {
            // fprintf(file_logerr,"\n Error opening NULL frequency spectrum file %s",file_H0f);
            log_error("Error opening NULL frequency spectrum file %s",
                      args.file_H0f);
            exit(1);
        }
        if (!feof(file_H0freq))
        {
            if ((cad = (char *)calloc(1026, sizeof(char))) == 0)
            {
                // fprintf(file_logerr,"Error allocating memory");
                log_fatal("Error allocating memory, cad");
                exit(1);
            }

            fgets(cad, 1024, file_H1freq);
        }
        else
        {
            // fprintf(file_logerr,"\n Error reading NULL frequency spectrum file %s",file_H0f);
            log_error("Error reading NULL frequency spectrum file %s",
                      args.file_H0f);
            exit(1);
        }

        npf = args.npops - 1;

        for (x = 0; x < npf; x++)
        {
            if (!feof(file_H0freq))
            {
                fgets(cad, 1024, file_H0freq);
                cad1 = cad;
            }
            else
            {
                // fprintf(file_logerr,"\n  Error reading NULL frequency spectrum file %s, line %d",file_H0f,x+2);
                log_error("Error reading NULL frequency spectrum file %s, line %d",
                          args.file_H0f,
                          x + 2);
                exit(1);
            }
            for (y = 1; y < args.vint_perpop_nsam[x]; y++)
            {
                if ((len = (int)strlen(cad1)) > 0)
                {
                    freqspH0[x][y] = (double)atof(cad1);
                    cad1 = strchr(cad1, '\t');
                    while (cad1 && (*cad1 == ' ' || *cad1 == '\t'))
                        cad1 += 1;
                }
            }
            thetaH0[x] = (double)atof(cad1);
        }
        free(cad);
        fclose(file_H0freq);
    }
    else
    {
        H0frq = 0;
        if (H1frq)
        {
            npf = args.npops - 1;
            for (x = 0; x < npf; x++)
            {
                for (y = 1; y < args.vint_perpop_nsam[x]; y++)
                {
                    freqspH0[x][y] = 1. / (double)y;
                }
            }
        }
    }

    /* -------------------------------------------------------------------- */
    /*GFF files*/

    /*observed_data = 0;*/

    if (args.gfffiles == 1 &&
        (strcmp(args.subset_positions, "synonymous") == 0 ||
         strcmp(args.subset_positions, "nonsynonymous") == 0 ||
         strcmp(args.subset_positions, "0-fold") == 0 ||
         strcmp(args.subset_positions, "2-fold") == 0 ||
         strcmp(args.subset_positions, "3-fold") == 0 ||
         strcmp(args.subset_positions, "4-fold") == 0 ||
         strcmp(args.subset_positions, "silent") == 0))
    {
        /* TODO: Reordenar, esto ya esta definido en el gff_data.c !
         */
        if (strcmp(args.code_name, "Nuclear_Universal") == 0)
        {
            args.genetic_code[0] = 'F';
            args.genetic_code[1] = 'F';
            args.genetic_code[2] = 'L';
            args.genetic_code[3] = 'L';
            args.genetic_code[4] = 'S';
            args.genetic_code[5] = 'S';
            args.genetic_code[6] = 'S';
            args.genetic_code[7] = 'S';
            args.genetic_code[8] = 'Y';
            args.genetic_code[9] = 'Y';
            args.genetic_code[10] = '*';
            args.genetic_code[11] = '*';
            args.genetic_code[12] = 'C';
            args.genetic_code[13] = 'C';
            args.genetic_code[14] = '*';
            args.genetic_code[15] = 'W';
            args.genetic_code[16] = 'L';
            args.genetic_code[17] = 'L';
            args.genetic_code[18] = 'L';
            args.genetic_code[19] = 'L';
            args.genetic_code[20] = 'P';
            args.genetic_code[21] = 'P';
            args.genetic_code[22] = 'P';
            args.genetic_code[23] = 'P';
            args.genetic_code[24] = 'H';
            args.genetic_code[25] = 'H';
            args.genetic_code[26] = 'Q';
            args.genetic_code[27] = 'Q';
            args.genetic_code[28] = 'R';
            args.genetic_code[29] = 'R';
            args.genetic_code[30] = 'R';
            args.genetic_code[31] = 'R';
            args.genetic_code[32] = 'I';
            args.genetic_code[33] = 'I';
            args.genetic_code[34] = 'I';
            args.genetic_code[35] = 'M';
            args.genetic_code[36] = 'T';
            args.genetic_code[37] = 'T';
            args.genetic_code[38] = 'T';
            args.genetic_code[39] = 'T';
            args.genetic_code[40] = 'N';
            args.genetic_code[41] = 'N';
            args.genetic_code[42] = 'K';
            args.genetic_code[43] = 'K';
            args.genetic_code[44] = 'S';
            args.genetic_code[45] = 'S';
            args.genetic_code[46] = 'R';
            args.genetic_code[47] = 'R';
            args.genetic_code[48] = 'V';
            args.genetic_code[49] = 'V';
            args.genetic_code[50] = 'V';
            args.genetic_code[51] = 'V';
            args.genetic_code[52] = 'A';
            args.genetic_code[53] = 'A';
            args.genetic_code[54] = 'A';
            args.genetic_code[55] = 'A';
            args.genetic_code[56] = 'D';
            args.genetic_code[57] = 'D';
            args.genetic_code[58] = 'E';
            args.genetic_code[59] = 'E';
            args.genetic_code[60] = 'G';
            args.genetic_code[61] = 'G';
            args.genetic_code[62] = 'G';
            args.genetic_code[63] = 'G';
        }
        else if (strcmp(args.code_name, "mtDNA_Drosophila") == 0)
        {
            args.genetic_code[0] = 'F';
            args.genetic_code[1] = 'F';
            args.genetic_code[2] = 'L';
            args.genetic_code[3] = 'L';
            args.genetic_code[4] = 'S';
            args.genetic_code[5] = 'S';
            args.genetic_code[6] = 'S';
            args.genetic_code[7] = 'S';
            args.genetic_code[8] = 'Y';
            args.genetic_code[9] = 'Y';
            args.genetic_code[10] = '*';
            args.genetic_code[11] = '*';
            args.genetic_code[12] = 'C';
            args.genetic_code[13] = 'C';
            args.genetic_code[14] = 'W';
            args.genetic_code[15] = 'W';
            args.genetic_code[16] = 'L';
            args.genetic_code[17] = 'L';
            args.genetic_code[18] = 'L';
            args.genetic_code[19] = 'L';
            args.genetic_code[20] = 'P';
            args.genetic_code[21] = 'P';
            args.genetic_code[22] = 'P';
            args.genetic_code[23] = 'P';
            args.genetic_code[24] = 'H';
            args.genetic_code[25] = 'H';
            args.genetic_code[26] = 'Q';
            args.genetic_code[27] = 'Q';
            args.genetic_code[28] = 'R';
            args.genetic_code[29] = 'R';
            args.genetic_code[30] = 'R';
            args.genetic_code[31] = 'R';
            args.genetic_code[32] = 'I';
            args.genetic_code[33] = 'I';
            args.genetic_code[34] = 'M';
            args.genetic_code[35] = 'M';
            args.genetic_code[36] = 'T';
            args.genetic_code[37] = 'T';
            args.genetic_code[38] = 'T';
            args.genetic_code[39] = 'T';
            args.genetic_code[40] = 'N';
            args.genetic_code[41] = 'N';
            args.genetic_code[42] = 'K';
            args.genetic_code[43] = 'K';
            args.genetic_code[44] = 'S';
            args.genetic_code[45] = 'S';
            args.genetic_code[46] = 'S';
            args.genetic_code[47] = 'S';
            args.genetic_code[48] = 'V';
            args.genetic_code[49] = 'V';
            args.genetic_code[50] = 'V';
            args.genetic_code[51] = 'V';
            args.genetic_code[52] = 'A';
            args.genetic_code[53] = 'A';
            args.genetic_code[54] = 'A';
            args.genetic_code[55] = 'A';
            args.genetic_code[56] = 'D';
            args.genetic_code[57] = 'D';
            args.genetic_code[58] = 'E';
            args.genetic_code[59] = 'E';
            args.genetic_code[60] = 'G';
            args.genetic_code[61] = 'G';
            args.genetic_code[62] = 'G';
            args.genetic_code[63] = 'G';
        }
        else if (strcmp(args.code_name, "mtDNA_Mammals") == 0)
        {
            args.genetic_code[0] = 'F';
            args.genetic_code[1] = 'F';
            args.genetic_code[2] = 'L';
            args.genetic_code[3] = 'L';
            args.genetic_code[4] = 'S';
            args.genetic_code[5] = 'S';
            args.genetic_code[6] = 'S';
            args.genetic_code[7] = 'S';
            args.genetic_code[8] = 'Y';
            args.genetic_code[9] = 'Y';
            args.genetic_code[10] = '*';
            args.genetic_code[11] = '*';
            args.genetic_code[12] = 'C';
            args.genetic_code[13] = 'C';
            args.genetic_code[14] = 'W';
            args.genetic_code[15] = 'W';
            args.genetic_code[16] = 'L';
            args.genetic_code[17] = 'L';
            args.genetic_code[18] = 'L';
            args.genetic_code[19] = 'L';
            args.genetic_code[20] = 'P';
            args.genetic_code[21] = 'P';
            args.genetic_code[22] = 'P';
            args.genetic_code[23] = 'P';
            args.genetic_code[24] = 'H';
            args.genetic_code[25] = 'H';
            args.genetic_code[26] = 'Q';
            args.genetic_code[27] = 'Q';
            args.genetic_code[28] = 'R';
            args.genetic_code[29] = 'R';
            args.genetic_code[30] = 'R';
            args.genetic_code[31] = 'R';
            args.genetic_code[32] = 'I';
            args.genetic_code[33] = 'I';
            args.genetic_code[34] = 'M';
            args.genetic_code[35] = 'M';
            args.genetic_code[36] = 'T';
            args.genetic_code[37] = 'T';
            args.genetic_code[38] = 'T';
            args.genetic_code[39] = 'T';
            args.genetic_code[40] = 'N';
            args.genetic_code[41] = 'N';
            args.genetic_code[42] = 'K';
            args.genetic_code[43] = 'K';
            args.genetic_code[44] = 'S';
            args.genetic_code[45] = 'S';
            args.genetic_code[46] = '*';
            args.genetic_code[47] = '*';
            args.genetic_code[48] = 'V';
            args.genetic_code[49] = 'V';
            args.genetic_code[50] = 'V';
            args.genetic_code[51] = 'V';
            args.genetic_code[52] = 'A';
            args.genetic_code[53] = 'A';
            args.genetic_code[54] = 'A';
            args.genetic_code[55] = 'A';
            args.genetic_code[56] = 'D';
            args.genetic_code[57] = 'D';
            args.genetic_code[58] = 'E';
            args.genetic_code[59] = 'E';
            args.genetic_code[60] = 'G';
            args.genetic_code[61] = 'G';
            args.genetic_code[62] = 'G';
            args.genetic_code[63] = 'G';
        }
        else if (strcmp(args.code_name, "Other") == 0)
        {
            ; /* TODO: Con verbose, dejar rastro para debuggar */
            /* Con control de errores, comprobar que es correcto (IUPAC,etc) */
        }
        else
        {
            // fprintf(file_logerr," %s: Unknown code, sorry", code_name);
            log_error(" %s: Unknown code, sorry", args.code_name);
            exit(1);
        }
    }

    /*ordering data: in case O is not a flag included*/
    if (args.int_total_nsam_order > 0 && args.int_total_nsam_order + !args.outgroup_presence != args.int_total_nsam)
    {
        // fprintf(file_logerr,"Error: the number of samples defined in -N and -O are different");
        log_error("Error: the number of samples defined in -N and -O are different");
        exit(1);
    }
    if (args.int_total_nsam_order == 0)
    {
        args.int_total_nsam_order = args.int_total_nsam - !args.outgroup_presence;
        if ((args.sort_nsam = (int *)calloc((unsigned long)args.int_total_nsam,
                                            sizeof(int))) == 0)
        {
            // fprintf(file_logerr,"Error allocating memory");
            log_fatal("Error allocating memory, sort_nsam");
            exit(1);
        }
        for (sort_index = 0; sort_index < args.int_total_nsam; sort_index++)
        {
            args.sort_nsam[sort_index] = sort_index;
        }
    }

    /* -------------------------------------------------------------------- */
    /* -------------------------------------------------------------------- */
    /* -------------------------------------------------------------------- */

    // if(args.output == 0 || args.output == 10) {
    if (args.output == OUTPUT_EXTENDED || args.output == OUTPUT_FULL_EXTENDED)
    {

        fprintf(file_output, MSTATSPOP );
        fprintf(file_output, "\nmstatspop ");
        for (x = 1; x < arg; x++)
        {
            fprintf(file_output, "%s ", argv[x]);
        }
        fprintf(file_output, "\n\n");
        fprintf(file_output,
                "\n****************************************************************************\n");
        fprintf(file_output,
                "*  NUCLEOTIDE VARIABILITY, NEUTRALITY TEST AND POPULATION DIFFERENTIATION  *\n");
        fprintf(file_output,
                "****************************************************************************\n");
    }

    /* STEP 3: And finally, some calculations ---------------------------------*/
    init_seed1(args.nseed);

    /*alloc memory for lengths of populations*/
    if ((nsites1_pop = (double *)calloc((unsigned long)args.npops,
                                        sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, nsites1_pop");
        exit(1);
    }
    if ((nsites2_pop = (double *)calloc((unsigned long)args.npops,
                                        sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, nsites2_pop");
        exit(1);
    }
    if ((nsites3_pop = (double *)calloc((unsigned long)args.npops,
                                        sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, nsites3_pop");
        exit(1);
    }
    if ((nsites1_pop_outg = (double *)calloc((unsigned long)args.npops,
                                             sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, nsites1_pop_outg");
        exit(1);
    }
    if ((nsites2_pop_outg = (double *)calloc((unsigned long)args.npops,
                                             sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, nsites2_pop_outg");
        exit(1);
    }
    if ((nsites3_pop_outg = (double *)calloc((unsigned long)args.npops,
                                             sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, nsites3_pop_outg");
        exit(1);
    }
    if ((rnsites1_pop = (double *)calloc((unsigned long)args.npops,
                                        sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, nsites1_pop");
        exit(1);
    }
    if ((rnsites2_pop = (double *)calloc((unsigned long)args.npops,
                                        sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, nsites2_pop");
        exit(1);
    }
    if ((rnsites3_pop = (double *)calloc((unsigned long)args.npops,
                                        sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, nsites3_pop");
        exit(1);
    }
    if ((rnsites1_pop_outg = (double *)calloc((unsigned long)args.npops,
                                             sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, nsites1_pop_outg");
        exit(1);
    }
    if ((rnsites2_pop_outg = (double *)calloc((unsigned long)args.npops,
                                             sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, nsites2_pop_outg");
        exit(1);
    }
    if ((rnsites3_pop_outg = (double *)calloc((unsigned long)args.npops,
                                             sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, nsites3_pop_outg");
        exit(1);
    }
    if ((anx = (double *)calloc((unsigned long)args.npops,
                                sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, anx");
        exit(1);
    }
    if ((bnx = (double *)calloc((unsigned long)args.npops,
                                sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, bnx");
        exit(1);
    }
    if ((anxo = (double *)calloc((unsigned long)args.npops,
                                 sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, anxo");
        exit(1);
    }
    if ((bnxo = (double *)calloc((unsigned long)args.npops,
                                 sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, bnxo");
        exit(1);
    }
    if ((lengthamng = (double **)calloc((unsigned long)args.npops,
                                        sizeof(double *))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, lengthamng");
        exit(1);
    }
    if ((lengthamng_outg = (double **)calloc((unsigned long)args.npops,
                                             sizeof(double *))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, lengthamng_outg");
        exit(1);
    }
    for (x = 0; x < args.npops; x++)
    {
        if ((lengthamng[x] = (double *)calloc((unsigned long)args.npops,
                                              sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"Error allocating memory");
            log_fatal("Error allocating memory, lengthamng[x]");
            exit(1);
        }
        if ((lengthamng_outg[x] = (double *)calloc((unsigned long)args.npops,
                                                   sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"Error allocating memory");
            log_fatal("Error allocating memory, lengthamng_outg[x]");
            exit(1);
        }
    }

    /* introduce the data and mask file (last if necessary) */

    if (args.mask_print == 1 && (args.formatfile == FASTA_FORMAT || (args.formatfile == TFA_FORMAT && ((args.slide == 0 && args.window == 0) && args.file_Wcoord[0] == '\0'))))
    {
        // if no input file is provided, create a default mask file
        if (args.file_in[0] == '\0')
        {
            sprintf(args.file_mas, "file_npop%d_nsam%d",
                    args.npops - !args.outgroup_presence,
                    args.int_total_nsam - !args.outgroup_presence);
            if (args.gfffiles == 1)
            {
                strcat(args.file_mas, "_");
                strcat(args.file_mas, args.subset_positions);
                strcat(args.file_mas, "_");
                strcat(args.file_mas, args.criteria_transcript);
            }
            if (!args.include_unknown)
                strcat(args.file_mas, "_ExcludeMissingVariantsmhits");
            else
                strcat(args.file_mas, "_IncludeMissingVariantsmhits");
            if (args.outgroup_presence == 0)
                strcat(args.file_mas, "_NOoutg");
            if (args.outgroup_presence == 1)
                strcat(args.file_mas, "_outg");
            if (args.ploidy[0] == '1')
                strcat(args.file_mas, "_ploidy1");
            if (args.ploidy[0] == '2')
                strcat(args.file_mas, "_ploidy2");
            strcat(args.file_mas, "_MASK.txt");

            if ((file_mask = fopen(args.file_mas, "w")) == 0)
            {
                // fprintf(file_logerr,"Error in mask file %s.",file_out);
                log_error("Error in mask file %s.", args.file_out);
                exit(1);
            }
            /*file_mask = stderr;*/ /* TODO: My GOD!!!! */
            /* TODO: formatfile=0 from pipeline => solo fasta o nbrf=> generara un
             * fichero de mascara, que no sirve pa na => lo tiramos a la basura =>
             * el stderr no es la basura!!!! Es para ERRORES!!! FIX IT.*/
        }
        else
        {
            el = strrchr(args.file_mas, '.');
            *el = '\0';
            // TODO :: FIXME
            sprintf(args.file_mas, "%s_npop%d_nsam%d", args.file_mas,
                    args.npops - !args.outgroup_presence,
                    args.int_total_nsam - !args.outgroup_presence);
            if (args.gfffiles == 1)
            {
                strcat(args.file_mas, "_");
                strcat(args.file_mas, args.subset_positions);
                strcat(args.file_mas, "_");
                strcat(args.file_mas, args.criteria_transcript);
            }
            if (!args.include_unknown)
                strcat(args.file_mas, "_ExcludeMissingVariantsmhits");
            else
                strcat(args.file_mas, "_IncludeMissingVariantsmhits");
            if (args.outgroup_presence == 0)
                strcat(args.file_mas, "_NOoutg");
            if (args.outgroup_presence == 1)
                strcat(args.file_mas, "_outg");
            if (args.ploidy[0] == '1')
                strcat(args.file_mas, "_ploidy1");
            if (args.ploidy[0] == '2')
                strcat(args.file_mas, "_ploidy2");
            strcat(args.file_mas, "_MASK.txt");

            if ((file_mask = fopen(args.file_mas, "w")) == 0)
            {
                // fprintf(file_logerr,"Error in mask file %s.",file_out);
                log_error("Error in mask file %s.", args.file_out);
                exit(1);
            }
        }
    }

    /* diploid or haploid */
    if (args.ploidy[0] == '2' && args.formatfile == FASTA_FORMAT)
    {
        args.int_total_nsam *= 2;

        for (x = 0; x < args.npops; x++)
        {
            args.vint_perpop_nsam[x] *= 2;
        }
        if (args.outgroup_presence == 0)
        {
            args.int_total_nsam -= 1;
            args.vint_perpop_nsam[args.npops - 1] -= 1;
        }
    }

    /*calloc struct stats*/
    /* TODO: Fix error checking */
    statistics = 0;
    if ((statistics = (struct stats *)calloc(1, sizeof(struct stats))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, statistics");
        exit(1);
    }
    /* length of each sequence, excluding non-tcga */
    if ((sum_sam = (double *)calloc(args.int_total_nsam + (!args.outgroup_presence),
                                    sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, sum_sam");
        exit(1);
    }
    /*tcga content for each sample*/
    if ((tcga = (double **)calloc(args.int_total_nsam + (!args.outgroup_presence),
                                  sizeof(double *))) == 0)
    {
        // fprintf(file_logerr,"Error allocating memory");
        log_fatal("Error allocating memory, tcga");
        exit(1);
    }
    for (x = 0; x < args.int_total_nsam + (!args.outgroup_presence); x++)
    {
        if ((tcga[x] = (double *)calloc(4, sizeof(double))) == 0)
        { /*tcga content*/
            // fprintf(file_logerr,"Error allocating memory");
            log_fatal("Error allocating memory, tcga[x]");
            exit(1);
        }
    }
    
    //DEFINE WEIGHTS IN CASE CALCULATE rSFS
    if(args.include_rsfs) {
        
        maxnsam = 0;
        for(x=0;x<args.npops-1;x++) {
            if(maxnsam < args.vint_perpop_nsam[x]) maxnsam = args.vint_perpop_nsam[x];
        }

        if ((rw = (struct rweights *)calloc(1, sizeof(struct rweights))) == 0)
        {
            // fprintf(file_logerr,"Error allocating memory");
            log_fatal("Error allocating memory, rweights");
            exit(1);
        }
        rw->ww   = (double **) calloc(args.int_total_nsam,sizeof(double *));
        rw->wt   = (double **) calloc(args.int_total_nsam,sizeof(double *));
        rw->wfl  = (double **) calloc(args.int_total_nsam,sizeof(double *));
        rw->wfw  = (double **) calloc(args.int_total_nsam,sizeof(double *));
        rw->wl   = (double **) calloc(args.int_total_nsam,sizeof(double *));
        rw->wwA  = (double **) calloc(args.int_total_nsam,sizeof(double *));
        rw->wtA  = (double **) calloc(args.int_total_nsam,sizeof(double *));
        rw->wflA = (double **) calloc(args.int_total_nsam,sizeof(double *));
        rw->wfwA = (double **) calloc(args.int_total_nsam,sizeof(double *));
        rw->wlA  = (double **) calloc(args.int_total_nsam,sizeof(double *));
        
        rw->wphi_i  = (double **) calloc(args.int_total_nsam,sizeof(double *));
        rw->wpsi_ij = (double ***) calloc(args.int_total_nsam,sizeof(double **));
        
        for(x=0;x<args.int_total_nsam;x++) {
            rw->ww[x]   = (double *) calloc(x+1,sizeof(double));
            rw->wt[x]   = (double *) calloc(x+1,sizeof(double));
            rw->wfl[x]  = (double *) calloc(x+1,sizeof(double));
            rw->wfw[x]  = (double *) calloc(x+1,sizeof(double));
            rw->wl[x]   = (double *) calloc(x+1,sizeof(double));
            rw->wwA[x]  = (double *) calloc(x+1,sizeof(double));
            rw->wtA[x]  = (double *) calloc(x+1,sizeof(double));
            rw->wflA[x] = (double *) calloc(x+1,sizeof(double));
            rw->wfwA[x] = (double *) calloc(x+1,sizeof(double));
            rw->wlA[x]  = (double *) calloc(x+1,sizeof(double));
            rw->wphi_i[x]  = (double *) calloc(x+1,sizeof(double));
            
            if(x>1) {
                if(args.outgroup_presence+args.force_outgroup==1) {
                    /*calculate the weights for unfolded*/
                    weights_unfolded(rw->ww[x], rw->wt[x],  rw->wfw[x], rw->wfl[x], rw->wl[x], x, 0.0);
                    /*calculate the weights for unfolded excluding singletons*/
                    weights_unfolded(rw->wwA[x], rw->wtA[x], rw->wfwA[x], rw->wflA[x], rw->wlA[x], x, 1.0/(double)x);
                }
                else {
                    /*calculate the weights for folded*/
                    weights_folded(rw->ww[x], rw->wt[x],  rw->wfl[x], rw->wl[x], rw->wphi_i[x], x, 0.0);
                    /*calculate the weights for folded excluding singletons*/
                    weights_folded(rw->wwA[x], rw->wtA[x], rw->wflA[x], rw->wlA[x], rw->wphi_i[x], x, 1.0/(double)x);
                }
            }
            
            rw->wpsi_ij[x] = (double **) calloc(x+1,sizeof(double *));
            for(z=0;z<x+1;z++) {
                rw->wpsi_ij[x][z] = (double *) calloc(maxnsam+1,sizeof(double));
                if(x>1) {
                    if(args.outgroup_presence+args.force_outgroup==1)
                        weights_unfolded_wpsi(rw->wpsi_ij[x][z], x, z, maxnsam);
                    else
                        weights_folded_wpsi(rw->wpsi_ij[x][z], x, z, maxnsam);
                }
            }
        }

        if ((statistics[0].rS = (double *)calloc(1 * args.npops,
                                                     sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rS");
            exit(1);
        }

        if ((statistics[0].rSo = (double *)calloc(1 * args.npops,
                                                     sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rSo");
            exit(1);
        }

        if ((statistics[0].rthetaS = (double *)calloc(1 * args.npops,
                                                     sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rthetaS");
            exit(1);
        }

        if ((statistics[0].rthetaSo = (double *)calloc(1 * args.npops,
                                                      sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rthetaSo");
            exit(1);
        }

        if ((statistics[0].rthetaT = (double *)calloc(1 * args.npops,
                                                     sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rthetaT");
            exit(1);
        }

        if ((statistics[0].rthetaTo = (double *)calloc(1 * args.npops,
                                                      sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rthetaTo");
            exit(1);
        }

        if ((statistics[0].rthetaFL = (double *)calloc(1 * args.npops,
                                                      sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rthetaFL");
            exit(1);
        }

        if ((statistics[0].rthetaFW = (double *)calloc(1 * args.npops,
                                                      sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rthetaFW");
            exit(1);
        }

        if ((statistics[0].rthetaL = (double *)calloc(1 * args.npops,
                                                     sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rthetaL");
            exit(1);
        }

        if ((statistics[0].rthetaSA = (double *)calloc(1 * args.npops,
                                                      sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rthetaSA");
            exit(1);
        }

        if ((statistics[0].rthetaTA = (double *)calloc(1 * args.npops,
                                                      sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rthetaTA");
            exit(1);
        }


        if ((statistics[0].rDtaj = (double *)calloc(1 * args.npops,
                                                   sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rDtaj");
            exit(1);
        }

        if ((statistics[0].rDfl = (double *)calloc(1 * args.npops,
                                                  sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rDfl");
            exit(1);
        }

        if ((statistics[0].rFfl = (double *)calloc(1 * args.npops,
                                                  sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rFfl");
            exit(1);
        }

        if ((statistics[0].rHnfw = (double *)calloc(1 * args.npops,
                                                   sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rHnfw");
            exit(1);
        }

        if ((statistics[0].rEz = (double *)calloc(1 * args.npops,
                                                 sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rEz");
            exit(1);
        }

        if ((statistics[0].rYach = (double *)calloc(1 * args.npops,
                                                   sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rYach");
            exit(1);
        }
        if ((statistics[0].rFH = (double *)calloc(1 * args.npops,
                                                 sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].rFH");
            exit(1);
        }
    }
    
    /* ************* READ FASTA DATA AND WRITE FILE_MASK ************ */
    if (args.formatfile == FASTA_FORMAT) // == 0
    {
        /*data from fasta*/
        chr_name = chr_name_array[0];
        if (args.length == 0)
            args.length = atol(chr_length_array[0]);
        first = 0;

        if (get_obsdata(
                        file_output,
                        &file_output_gz,
                        file_input,
                        &file_input_gz,
                        // file_logerr,&file_logerr_gz,
                        file_mask,
                        // args.file_GFF,
                        // args.gfffiles,
                        // args.subset_positions,
                        // args.genetic_code,
                        &matrix_pol,
                        &matrix_freq,
                        &matrix_pos,
                        &args.length,
                        &length_seg,
                        &length_al,
                        &length_al_real,
                        // argc,
                        // args.ploidy,
                        // args.vint_perpop_nsam,
                        // args.npops,
                        &svratio,
                        &missratio,
                        // args.include_unknown,
                        sum_sam,
                        tcga,
                        &matrix_sv,
                        &nmhits,
                        // args.output,
                        // args.outgroup_presence,
                        // args.criteria_transcript,
                        nsites1_pop,
                        nsites1_pop_outg,
                        nsites2_pop,
                        nsites2_pop_outg,
                        nsites3_pop,
                        nsites3_pop_outg,
                        rnsites1_pop,
                        rnsites1_pop_outg,
                        rnsites2_pop,
                        rnsites2_pop_outg,
                        rnsites3_pop,
                        rnsites3_pop_outg,
                        anx,
                        bnx,
                        anxo,
                        bnxo,
                        lengthamng,
                        lengthamng_outg,
                        // args.sort_nsam,
                        &matrix_pol_tcga,
                        chr_name,
                        first,
                        &args) == 0)
        {
            // fprintf(file_logerr, "Error processing input data.\n");
            log_error("Error processing input data.");
            exit(1);
        }

        fzclose(file_input, &file_input_gz);
        if (file_mask)
            fclose(file_mask);

        args.niterdata = 1;
    }

    /*DEFINE PARAMETERS FOR MS MASK FILE AND READ DATA IF DEFINED*/
    // if (args.formatfile == 1 || args.formatfile == 2)
    if (args.formatfile == MS_FORMAT || args.formatfile == MS_X_FORMAT)
    {                            /* MASK FILE MS FORMAT => Ponerlo en una funcion */
        args.niter = 0;            /*permutation tests invalidated*/
        svratio = args.ms_svratio; /*-10000*/
        length_al_real = args.length;
        length_al = length_al_real;

        if ((vector_mask = (float *)calloc((unsigned int)args.length,
                                           sizeof(float))) == 0)
        {
            // fprintf(file_logerr,"Error allocating memory");
            log_fatal("Error allocating memory, vector_mask");
            exit(1);
        }

        if ((matrix_mask = (int *)calloc((unsigned int)(args.int_total_nsam + 1) * args.length,
                                         sizeof(int))) == 0)
        {
            // fprintf(file_logerr,"Error allocating memory");
            log_fatal("Error allocating memory, matrix_mask");
            exit(1);
        }

        /*if(include_unknown) {*/
        /*if(file_mas[0] == '-'  && file_mas[1] == '1') {*/ /*all positions are accepted*/ /*
                                                                                            for(li=0;li<length;li++) vector_mask[li] = 1.0;
                                                                                            if(freq_missing_ms) {
                                                                                            for(li=0;li<(int_total_nsam)*length;li++) {
                                                                                            if(rand() > freq_missing_ms) matrix_mask[li] = 0;
                                                                                            else matrix_mask[li] = -1;
                                                                                            }
                                                                                            }
                                                                                            else {
                                                                                            for(li=0;li<(int_total_nsam)*length;li++) matrix_mask[li] = 0;
                                                                                            }
                                                                                            }
                                                                                            else {*/
        if (args.file_mas[0] != '-')                                                       /*file_mask defined*/
        {
            if ((file_mask = fopen(args.file_mas, "r")) == 0)
            {
                // fprintf(file_logerr,"\n  It is not possible to open the input mask file %s.",file_mas);
                log_error("It is not possible to open the input mask file %s.",
                          args.file_mas);
                exit(1);
            }
            li = 0;
            n = 0;
            vli = 0;
            *c = fgetc(file_mask);
            while (*c != 0 && *c != -1 && n < args.int_total_nsam - !args.outgroup_presence + 1 && li < args.length)
            {
                while (*c != 10 && *c != 13 && *c != 0 && *c != -1)
                {
                    if (n == 0)
                    {
                        vvm[vli] = *c;
                        while ((*c = fgetc(file_mask)) != 32 && *c != 10 && *c != 13 && vli < 10)
                        {
                            vli++;
                            vvm[vli] = *c;
                        }
                        vvm[vli + 1] = (char)0;
                        vector_mask[(unsigned int)li] = (float)atof(vvm); /*value of the position, usually 1 except for noncounting values or Syn/Nsyn (between 0 and 1)*/
                        for (vli = 0; vli < 10; vli++)
                            vvm[vli] = (char)0;
                        vli = 0;
                    }
                    else
                    {
                        matrix_mask[(unsigned int)(n - 1) * (unsigned int)args.length + (unsigned int)li] = (int)atoi(c) - 1; /*in case not defined file, all values are zero, in case file exist, normal is zero, missing is -1*/
                        if ((int)atoi(c) == 0 &&
                            args.include_unknown == 0)
                            vector_mask[(unsigned int)li] = 0.; /*not counting missing values when not included*/
                    }
                    if (*c == 10 || *c == 13 || *c == 0 || *c == -1)
                        break;
                    while ((*c = fgetc(file_mask)) == 32)
                        ;
                    li++;
                }
                if (*c == 10 || *c == 13 || li > args.length)
                {
                    if (li > args.length)
                    {
                        // fprintf(file_logerr,"\n  Error: Length of rows in mask file %s are longer than defined (row %ld is %ld > %ld). ",file_mas, n, li, length);
                        log_error("Error: Length of rows in mask file %s are longer than defined (row %ld is %ld > %ld). ",
                                  args.file_mas,
                                  n,
                                  li,
                                  args.length);
                        exit(1);
                    }
                    n++;
                    li = 0;
                    *c = fgetc(file_mask);
                }
            }
            fclose(file_mask);

            if (args.outgroup_presence == 0)
            { /*the "cryptic" outgroup is all 1*/
                for (li = 0; li < args.length; li++)
                {
                    matrix_mask[(unsigned int)(args.int_total_nsam) * (unsigned int)args.length + (unsigned int)li] = 0;
                }
            }
            /*}*/
            /*}*/
            /*if(include_unknown) {*/
            if ((sum_sam_mask = (double *)calloc(args.int_total_nsam,
                                                 sizeof(double))) == 0)
            {
                // fprintf(file_logerr,"Error allocating memory");
                log_fatal("Error allocating memory, sum_sam_mask");
                exit(1);
            }
            /*}*/
            length_mask = 0.;
            length_mask_real = 0;
            li2 = 0;
            for (li = 0; li < args.length; li++)
            {
                if (vector_mask[li] > 0.)
                {
                    if (args.npops > 1)
                    {
                        x = 0;
                        for (n = 0; n < args.int_total_nsam - args.vint_perpop_nsam[args.npops - 1]; n++)
                        {
                            li2 += 1;
                            if (matrix_mask[(unsigned int)n * (unsigned int)args.length + (unsigned int)li] == 0)
                            {
                                sum_sam_mask[n] += 1.; /*1 sum, 0 no sum*/
                                x += 1;
                            }
                            else
                                missratio += 1.;
                        }
                    }
                    else
                        x = 1;
                    y = 0;
                    for (n = args.int_total_nsam - args.vint_perpop_nsam[args.npops - 1]; n < args.int_total_nsam; n++)
                    {
                        if (args.npops == 1)
                            li2 += 1;
                        if (matrix_mask[(unsigned int)n * (unsigned int)args.length + (unsigned int)li] == 0)
                        {
                            sum_sam_mask[n] += 1.; /*1 sum, 0 no sum*/
                            y += 1;
                        }
                        else
                        {
                            if (args.npops == 1)
                                missratio += 1.;
                        }
                    }
                    if (x > 0 && y > 0)
                    {
                        length_mask_real += (long int)1;              /*length where the outgroup and one of the rest lines exists*/
                        length_mask += vector_mask[(unsigned int)li]; /*length where the outgroup and one of the rest lines exists*/
                    }
                    else
                    {
                        if (x == 0 || (y == 0 && x > 0))
                        {
                            missratio -= ((args.int_total_nsam - args.vint_perpop_nsam[args.npops - 1]) - x);
                            li2 -= (args.int_total_nsam - args.vint_perpop_nsam[args.npops - 1]);
                        }
                        if (y == 0 && args.npops == 1)
                        {
                            missratio -= (args.int_total_nsam);
                            li2 -= (args.int_total_nsam);
                        }
                    }
                }
            }
            if (li2)
                missratio = (double)missratio / (double)li2;
            else
                missratio = (double)1;
            length_al = length_mask;
            length_al_real = length_mask_real;
            /*
             li=0;
             for(n=0;n<int_total_nsam;n++) {
             li += sum_sam_mask[n];
             }
             missratio = 1.0 - (double)li2/(double)(li2+missratio);
             */
            /*}*/
        }
    }
    /* MASK FILE MS FORMAT => Ponerlo en una funcion HASTA AQUI ^^^^ */

    /*calloc pointers in struct*/
    if ((statistics[0].Sanc = (long int *)calloc(4 * args.npops,
                                                 sizeof(long int))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].Sanc");
        exit(1);
    }
    if ((statistics[0].piw = (double *)calloc(1 * args.npops,
                                              sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].piw");
        exit(1);
    }
    if ((statistics[0].piwHKY = (double *)calloc(1 * args.npops,
                                                 sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].piwHKY");
        exit(1);
    }
    if ((statistics[0].hapw = (double *)calloc(1 * args.npops,
                                               sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].hapw");
        exit(1);
    }


    if(args.include_rsfs == 0) {
        
        if ((statistics[0].pia = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                  sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].pia");
            exit(1);
        }
        if ((statistics[0].piT = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                  sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].piT");
            exit(1);
        }
        
        if ((statistics[0].piant = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                    sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].piant");
            exit(1);
        }
        if ((statistics[0].piTnt = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                    sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].piTnt");
            exit(1);
        }
        
        if ((statistics[0].fst = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                  sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].fst");
            exit(1);
        }
        
        if ((statistics[0].piaHKY = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                     sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].piaHKY");
            exit(1);
        }
        if ((statistics[0].piTHKY = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                     sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].piTHKY");
            exit(1);
        }
        
        if ((statistics[0].fstHKY = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                     sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].fstHKY");
            exit(1);
        }
        
        if ((statistics[0].fst1all = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                      sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].fst1all");
            exit(1);
        }
        
        if ((statistics[0].hapa = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                   sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].hapa");
            exit(1);
        }
        
        if ((statistics[0].hapT = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                   sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].hapT");
            exit(1);
        }
        
        if ((statistics[0].fsth = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                   sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].fsth");
            exit(1);
        }
        
        if ((statistics[0].fsth1all = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                       sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].fsth1all");
            exit(1);
        }
        
        if ((statistics[0].Gst = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                  sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].Gst");
            exit(1);
        }
    }
    if ((statistics[0].S = (double *)calloc(1 * args.npops,
                                            sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].S");
        exit(1);
    }

    if ((statistics[0].So = (double *)calloc(1 * args.npops,
                                             sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].So");
        exit(1);
    }

    if ((statistics[0].thetaS = (double *)calloc(1 * args.npops,
                                                 sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].thetaS");
        exit(1);
    }

    if ((statistics[0].thetaSo = (double *)calloc(1 * args.npops,
                                                  sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].thetaSo");
        exit(1);
    }

    if ((statistics[0].thetaT = (double *)calloc(1 * args.npops,
                                                 sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].thetaT");
        exit(1);
    }

    if ((statistics[0].thetaTo = (double *)calloc(1 * args.npops,
                                                  sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].thetaTo");
        exit(1);
    }

    if ((statistics[0].thetaTHKY = (double *)calloc(1 * args.npops,
                                                    sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].thetaTHKY");
        exit(1);
    }

    if ((statistics[0].thetaFL = (double *)calloc(1 * args.npops,
                                                  sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].thetaFL");
        exit(1);
    }

    if ((statistics[0].thetaFW = (double *)calloc(1 * args.npops,
                                                  sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].thetaFW");
        exit(1);
    }

    if ((statistics[0].thetaL = (double *)calloc(1 * args.npops,
                                                 sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].thetaL");
        exit(1);
    }

    if ((statistics[0].thetaSA = (double *)calloc(1 * args.npops,
                                                  sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].thetaSA");
        exit(1);
    }

    if ((statistics[0].thetaTA = (double *)calloc(1 * args.npops,
                                                  sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].thetaTA");
        exit(1);
    }

    if ((statistics[0].K = (double *)calloc(1 * args.npops,
                                            sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].K");
        exit(1);
    }

    if ((statistics[0].KHKY = (double *)calloc(1 * args.npops,
                                               sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].KHKY");
        exit(1);
    }

    if ((statistics[0].Dtaj = (double *)calloc(1 * args.npops,
                                               sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].Dtaj");
        exit(1);
    }

    if ((statistics[0].Dfl = (double *)calloc(1 * args.npops,
                                              sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].Dfl");
        exit(1);
    }

    if ((statistics[0].Ffl = (double *)calloc(1 * args.npops,
                                              sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].Ffl");
        exit(1);
    }

    if ((statistics[0].Hnfw = (double *)calloc(1 * args.npops,
                                               sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].Hnfw");
        exit(1);
    }

    if ((statistics[0].Ez = (double *)calloc(1 * args.npops,
                                             sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].Ez");
        exit(1);
    }

    if ((statistics[0].Yach = (double *)calloc(1 * args.npops,
                                               sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].Yach");
        exit(1);
    }
    if ((statistics[0].FH = (double *)calloc(1 * args.npops,
                                             sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].FH");
        exit(1);
    }

    if ((statistics[0].R2 = (double *)calloc(1 * args.npops,
                                             sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].R2");
        exit(1);
    }

    if ((statistics[0].Fs = (double *)calloc(1 * args.npops,
                                             sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].Fs");
        exit(1);
    }

    if ((statistics[0].nhpop = (int *)calloc(1 * args.npops, sizeof(int))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].nhpop");
        exit(1);
    }

    if ((statistics[0].length = (double *)calloc(1 * args.npops,
                                                 sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].length");
        exit(1);
    }
    if ((statistics[0].length2 = (double *)calloc(1 * args.npops,
                                                  sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].length2");
        exit(1);
    }
    if ((statistics[0].lengthamng = (double **)calloc(args.npops,
                                                      sizeof(double *))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].lengthamng");
        exit(1);
    }
    if ((statistics[0].lengthamng_outg = (double **)calloc(args.npops,
                                                           sizeof(double *))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].lengthamng_outg");
        exit(1);
    }
    for (x = 0; x < args.npops; x++)
    {
        if ((statistics[0].lengthamng[x] = (double *)calloc(args.npops,
                                                            sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].lengthamng[x]");
            exit(1);
        }
        if ((statistics[0].lengthamng_outg[x] = (double *)calloc(args.npops,
                                                                 sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].lengthamng_outg[x]");
            exit(1);
        }
    }
    if ((statistics[0].total_tcga = (double *)calloc(4, sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].total_tcga");
        exit(1);
    }

    if ((statistics[0].tcga = (double **)calloc(1 * args.npops,
                                                sizeof(double *))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].tcga");
        exit(1);
    }

    if ((statistics[0].sv = (double ***)calloc(1 * args.npops,
                                               sizeof(double **))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].sv");
        exit(1);
    }
    if ((statistics[0].svT = (double ***)calloc(1 * args.npops,
                                                sizeof(double **))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].svT");
        exit(1);
    }

    if ((statistics[0].freq = (long int **)calloc(1 * args.npops,
                                                  sizeof(long int *))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].freq");
        exit(1);
    }

    if ((statistics[0].freqh = (long int **)calloc(1 * args.npops,
                                                   sizeof(long int *))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].freqh");
        exit(1);
    }

    if ((statistics[0].ToH0_ii = (double *)calloc(1 * args.npops,
                                                  sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].ToH0_ii");
        exit(1);
    }

    if ((statistics[0].To_ii = (double *)calloc(1 * args.npops,
                                                sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].To_ii");
        exit(1);
    }

    if ((statistics[0].To_00 = (double *)calloc(1 * args.npops,
                                                sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].To_00");
        exit(1);
    }

    if ((statistics[0].To_i0 = (double *)calloc(1 * args.npops,
                                                sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].To_i0");
        exit(1);
    }

    if ((statistics[0].ToH0_00 = (double *)calloc(1 * args.npops,
                                                  sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].ToH0_00");
        exit(1);
    }

    if ((statistics[0].To_Qc_ii = (double *)calloc(1 * args.npops,
                                                   sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].To_Qc_ii");
        exit(1);
    }

    if ((statistics[0].To_Qw_ii = (double *)calloc(1 * args.npops,
                                                   sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].To_Qw_ii");
        exit(1);
    }

    if ((statistics[0].To_Lc_ii = (double *)calloc(1 * args.npops,
                                                   sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].To_Lc_ii");
        exit(1);
    }

    if ((statistics[0].Rm = (long int *)calloc(1 * args.npops,
                                               sizeof(long int))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].Rm");
        exit(1);
    }

    if ((statistics[0].ZnA = (double *)calloc(1 * args.npops,
                                              sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].ZnA");
        exit(1);
    }

    if ((statistics[0].mdsd = (double *)calloc(1 * args.npops,
                                               sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].mdsd");
        exit(1);
    }

    if ((statistics[0].mdg1 = (double *)calloc(1 * args.npops,
                                               sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].mdg1");
        exit(1);
    }

    if ((statistics[0].mdg2 = (double *)calloc(1 * args.npops,
                                               sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].mdg2");
        exit(1);
    }
    if ((statistics[0].anx = (double *)calloc(1 * args.npops,
                                              sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].anx");
        exit(1);
    }
    if ((statistics[0].bnx = (double *)calloc(1 * args.npops,
                                              sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].bnx");
        exit(1);
    }
    if ((statistics[0].anxo = (double *)calloc(1 * args.npops,
                                               sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].anxo");
        exit(1);
    }
    if ((statistics[0].bnxo = (double *)calloc(1 * args.npops,
                                               sizeof(double))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].bnxo");
        exit(1);
    }
    if ((statistics[0].mdw = (double **)calloc(1 * args.npops,
                                               sizeof(double *))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].mdw");
        exit(1);
    }

    if ((statistics[0].linefreq = (double **)calloc(args.int_total_nsam,
                                                    sizeof(double *))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].linefreq");
        exit(1);
    }
    if ((statistics[0].popfreq = (double **)calloc(1 * args.npops,
                                                   sizeof(double *))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].popfreq");
        exit(1);
    }

    if (H1frq)
    {
        statistics[0].H1freq = freqspH1; /*pointer to Alternative Frequency Spectrum*/
        statistics[0].H0freq = freqspH0; /*pointer to Null Frequency Spectrum*/
        statistics[0].thetaH1 = thetaH1; /*pointer to Alternative theta*/
        statistics[0].thetaH0 = thetaH0; /*pointer to Null theta*/
    }

    for (x = 0; x < args.npops; x++)
    {
        if ((statistics[0].tcga[x] = (double *)calloc(4, sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].tcga");
            exit(1);
        }

        if ((statistics[0].freq[x] = (long int *)calloc(args.vint_perpop_nsam[x],
                                                        sizeof(long int))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].freq");
            exit(1);
        }

        if ((statistics[0].freqh[x] = (long int *)calloc(args.int_total_nsam,
                                                         sizeof(long int))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].freqh");
            exit(1);
        }

        if ((statistics[0].sv[x] = (double **)calloc(1 * args.npops,
                                                     sizeof(double *))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].sv");
            exit(1);
        }
        if ((statistics[0].svT[x] = (double **)calloc(1 * args.npops,
                                                      sizeof(double *))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].svT");
            exit(1);
        }

        if ((statistics[0].mdw[x] = (double *)calloc((args.vint_perpop_nsam[x] * (args.vint_perpop_nsam[x] - 1)) / 2,
                                                     sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].mdw");
            exit(1);
        }

        for (y = 0; y < args.npops; y++)
        {
            if ((statistics[0].sv[x][y] = (double *)calloc(2,
                                                           sizeof(double))) == 0)
            {
                // fprintf(file_logerr,"\n  Error allocating memory.");
                log_fatal("Error allocating memory, statistics[0].sv");
                exit(1);
            }
            if ((statistics[0].svT[x][y] = (double *)calloc(2,
                                                            sizeof(double))) == 0)
            {
                // fprintf(file_logerr,"\n  Error allocating memory.");
                log_fatal("Error allocating memory, statistics[0].svT");
                exit(1);
            }
      
        }
        if ((statistics[0].popfreq[x] = (double *)calloc(args.int_total_nsam + 1,
                                                         sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].popfreq");
            exit(1);
        }
    }
    for (x = 0; x < args.int_total_nsam; x++)
    {
        if ((statistics[0].linefreq[x] = (double *)calloc(args.int_total_nsam + 1,
                                                          sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].linefreq");
            exit(1);
        }
    }
    /*in case r2i_ploidies is undefined*/
    if (args.r2i_ploidies == 0)
    {
        if ((args.r2i_ploidies = (int *)calloc(2, sizeof(int))) == 0)
        {
            // fprintf(file_logerr,"\nError: memory not reallocated. mstatspop.c.00 \n");
            log_fatal("Error: memory not reallocated, r2i_ploidies mstatspop.c.00 \n");
            exit(1);
        }
        args.r2i_ploidies[0] = 1;
        if (args.ploidy[0] == '1')
            args.r2i_ploidies[1] = 1;
        if (args.ploidy[0] == '2')
            args.r2i_ploidies[1] = 2;
    }
    /*allocate R2p*/
    if ((statistics[0].R2p = (double **)calloc(args.r2i_ploidies[0],
                                               sizeof(double *))) == 0)
    {
        // fprintf(file_logerr,"\n  Error allocating memory.");
        log_fatal("Error allocating memory, statistics[0].R2p");
        exit(1);
    }
    for (x = 0; x < args.r2i_ploidies[0]; x++)
    {
        if ((statistics[0].R2p[x] = (double *)calloc(1 * args.npops,
                                                     sizeof(double))) == 0)
        {
            // fprintf(file_logerr,"\n  Error allocating memory.");
            log_fatal("Error allocating memory, statistics[0].R2p");
            exit(1);
        }
    }

    /* Aqui el PROCESAMIENTO DE DATOS (y leer los fragmentos de datos (MS o TFA)*/
    int chr_length_array_size = 0;
    long int processed_sites = 0;
    for (first = 0; first < nscaffolds; first++)
    {
    
        // if (args.formatfile == 3)
        if (args.formatfile == TFA_FORMAT)
        {
            chr_name = chr_name_array[first];
            log_info("Processing scaffold %s ....", chr_name);
            args.length = atol(chr_length_array[first]);
            chr_length_array_size = atol(chr_length_array[first]);
            processed_sites = 0;
            wgenes = 0;
            nwindows = 0;
            /*read the file for weigth for positions, if included*/
            if (args.file_Wcoord[0] != '\0')
            {
                if ((file_wcoor = fzopen(args.file_Wcoord,
                                         "r",
                                         &file_wcoor_gz)) == 0)
                {
                    // fprintf(file_logerr,"\n It is not possible to open the coordinates file %s\n", file_Wcoord);
                    log_fatal("It is not possible to open the coordinates file %s\n",
                              args.file_Wcoord);
                    exit(1);
                }
                if (
                    read_coordinates(
                                     file_wcoor,
                                     &file_wcoor_gz,
                                     file_output,
                                     &file_output_gz,
                                     // file_logerr,
                                     // &file_logerr_gz,
                                     &wgenes,
                                     &nwindows,
                                     chr_name) == 0)
                {
                    exit(1);
                }
                args.window = -1;
                args.slide = -1;
                fzclose(file_wcoor, &file_wcoor_gz);
            }
            // Replaced with the usage of wtfasta_file structure
            // if (args.file_wps[0] != '\0' && first == 0)
            // {
            //   if ((file_ws = fzopen(args.file_wps, "r", &file_ws_gz)) == 0)
            //   {
            //     // fprintf(file_logerr,"\n It is not possible to open the weights file %s\n", file_wps);
            //     log_fatal("It is not possible to open the weights file %s\n", args.file_wps);
            //     exit(1);
            //   }
            //   load_index_from_file(file_ws_gz.index_file_name, &index_w);
            // }
        }
        li = 0;
        while (li < args.niterdata)
        {
            flaghky = 1; /* TODO: corrector de multiple hits... en ms format
                          es tratado de manera diferente */

            /*read the ms format for each iteration*/
            /* FUNCTION TO READ MS FILES*/
            // if (args.formatfile == 1 || args.formatfile == 2)
            if (args.formatfile == MS_FORMAT || args.formatfile == MS_X_FORMAT)
            {
                /*read ms file*/
                if (get_msdata(
                               file_input, &file_input_gz,
                               // file_logerr,&file_logerr_gz,
                               &matrix_pol, &matrix_freq, &matrix_pos,
                               &length_seg, args.vint_perpop_nsam, args.npops, args.int_total_nsam,
                               args.length, &nmhits, matrix_mask, vector_mask, args.ms_svratio,
                               &vector_priors, &npriors, &matrix_sv,
                               args.outgroup_presence, args.force_outgroup, args.freq_revert, sum_sam,
                               nsites1_pop, nsites1_pop_outg, args.formatfile,
                               nsites2_pop, nsites2_pop_outg, nsites3_pop, nsites3_pop_outg,
                               rnsites1_pop,
                               rnsites1_pop_outg,
                               rnsites2_pop,
                               rnsites2_pop_outg,
                               rnsites3_pop,
                               rnsites3_pop_outg,
                               anx, bnx, anxo, bnxo, lengthamng, lengthamng_outg,
                               args.include_unknown, args.file_mas, args.freq_missing_ms, args.kind_length,
                               &sum_sam_mask, &length_mask, &length_mask_real, &missratio, args.location_missing_ms, args.sort_nsam))
                {
                    // fprintf(file_logerr,"\nError processing ms data.\n");
                    log_fatal("Error processing ms data.\n");
                    exit(1);
                }
                /*
                 length_al = 0.;
                 for(li2=0;li2<length;li2++) {
                 length_al += vector_mask[li2];
                 }
                 */
                flaghky = 0;
                for (x = 0; x < args.int_total_nsam; x++)
                {
                    for (y = 0; y < args.int_total_nsam + 1; y++)
                    {
                        statistics[0].linefreq[x][y] = 0.;
                    }
                }
                for (x = 0; x < args.npops; x++)
                {
                    for (y = 0; y < args.int_total_nsam + 1; y++)
                    {
                        statistics[0].popfreq[x][y] = 0.;
                    }
                }

                if (args.include_unknown)
                {
                    if (!(args.file_mas[0] == '-' && args.file_mas[1] == '1'))
                    {
                        for (x = 0; x < args.int_total_nsam; x++)
                            sum_sam[x] = sum_sam_mask[x] /*- nmhits*/;
                        length_al = length_mask - nmhits;
                    }
                    else
                    { /*here is counting the positions with 0 and 1 (not gaps[8 or 9] and mhits)*/
                        /*in case including gaps in ms format, it must be included the total length region!!! (not only variants!!)*/
                        /*if(freq_missing_ms > 0.) {*/
                        for (x = 0; x < args.int_total_nsam; x++)
                            sum_sam[x] = (sum_sam_mask[x] /*- nmhits*/);
                        length_al = length_mask - nmhits;
                        /*}
                         else {
                         length_al = length_mask - nmhits;
                         }
                         */
                    }
                }
                else
                {
                    for (x = 0; x < args.int_total_nsam; x++)
                    {
                        if (args.length - nmhits > 0)
                            sum_sam[x] = args.length - nmhits;
                        else
                            sum_sam[x] = 0;
                    }
                    if (args.length - nmhits > 0)
                        length_al = args.length - nmhits;
                    else
                        length_al = 0;
                }
                z = 0;
                for (x = 0; x < args.npops; x++)
                {
                    for (y = z; y < z + args.vint_perpop_nsam[x]; y++)
                    {
                        statistics[0].length[x] = 0;
                        for (w = 0; w < 4; w++)
                        {
                            statistics[0].total_tcga[w] = 0;
                            statistics[0].tcga[x][w] = 0;
                        }
                    }
                    z += args.vint_perpop_nsam[x];
                }
            }

            /* FUNCTION TO READ TRANSPOSED FASTA (TFA) FILES*/
            // if (args.formatfile == 3)
            if (args.formatfile == TFA_FORMAT)
            {
                if (get_tfadata(
                                file_output, 
                                &file_output_gz,
                                &tfasta,
                                // file_input, 
                                // &file_input_gz, 
                                // &index_input,
                                wtfasta,
                                // args.file_wps, 
                                // file_ws, 
                                //&file_ws_gz, 
                                //&index_w,
                                // file_logerr,&file_logerr_gz,
                                &matrix_pol, 
                                &matrix_freq,
                                &matrix_pos, 
                                &args.length, 
                                &length_seg, 
                                &length_al,
                                &length_al_real, 
                                //args.argc, 
                                //args.vint_perpop_nsam, 
                                //args.npops,
                                &svratio, 
                                &missratio, 
                                //args.include_unknown, 
                                sum_sam, 
                                tcga,
                                &matrix_sv, 
                                &nmhits, 
                                //args.output, 
                                //args.outgroup_presence,
                                nsites1_pop, 
                                nsites1_pop_outg,
                                nsites2_pop, 
                                nsites2_pop_outg,
                                nsites3_pop, 
                                nsites3_pop_outg,
                                rnsites1_pop,
                                rnsites1_pop_outg,
                                rnsites2_pop,
                                rnsites2_pop_outg,
                                rnsites3_pop,
                                rnsites3_pop_outg,
                                anx,
                                bnx,
                                anxo, 
                                bnxo,
                                lengthamng, 
                                lengthamng_outg, 
                                args.sort_nsam,
                                wgenes, 
                                nwindows,
                                //args.first_slide, 
                                //args.slide, 
                                //args.window,
                                //args.Physical_length,
                                &li, 
                                &npriors, 
                                &vector_priors,
                                &matrix_pol_tcga,
                                chr_name, 
                                first, 
                                &args) == 0)
                {
                    /*printf("End processing input tfa data.\n");*/
                    break;
                    /*exit(1);*/
                }

                // report progress how many bases have been processed
                processed_sites += args.length;
        
                // only report progress 
                if (processed_sites % (chr_length_array_size / 100) == 0)
                {
                    // calculate the percentage of processed sites
                    double processed_sites_percentage = (double)processed_sites / chr_length_array_size * 100;
                    // report process
                    log_info("Processing scaffold %s %ld/%ld sites (%.2f%%)",
                             chr_name,
                             processed_sites,
                             chr_length_array_size,
                             processed_sites_percentage);
                }
                // calculate the percentage of processed sites
                // double processed_sites_percentage = (double)processed_sites / chr_length_array_size * 100;
                // // report process 
                // log_info("Processing scaffold %s %ld/%ld sites (%.2f%%)", processed_sites, chr_length_array_size, processed_sites_percentage);


                /*When the tfa file is finished then li=0, otherwise li=-1 and the loop continue*/
                z = 0;
                for (x = 0; x < args.npops; x++)
                {
                    for (y = z; y < z + args.vint_perpop_nsam[x]; y++)
                    {
                        statistics[0].length[x] = 0;
                        for (w = 0; w < 4; w++)
                        {
                            statistics[0].total_tcga[w] = 0;
                            statistics[0].tcga[x][w] = 0;
                        }
                    }
                    z += args.vint_perpop_nsam[x];
                }
            }

            /*calculate number of effective nucleotides per population and tcga frequencies*/
            z = 0;
            for (x = 0; x < args.npops; x++)
            {
                /*
                 nsites1_pop[x] -= nmhits;
                 nsites1_pop_outg[x] -= nmhits;
                 */
                if (args.outgroup_presence)
                    statistics[0].length2[x] = nsites2_pop_outg[x];
                else
                    statistics[0].length2[x] = nsites2_pop[x];
        
                statistics[0].anx[x] = anx[x];
                statistics[0].bnx[x] = bnx[x];
                statistics[0].anxo[x] = anxo[x];
                statistics[0].bnxo[x] = bnxo[x];

                for (y = z; y < z + args.vint_perpop_nsam[x]; y++)
                {
                    statistics[0].length[x] += sum_sam[y];
                    for (w = 0; w < 4; w++)
                    {
                        statistics[0].total_tcga[w] += tcga[y][w];
                        statistics[0].tcga[x][w] += tcga[y][w];
                    }
                }
                z += args.vint_perpop_nsam[x];
            }

            /*include lengthamng and calculations*/
            for (x = 0; x < args.npops - !args.outgroup_presence; x++)
            {
                for (y = 0; y < args.npops - !args.outgroup_presence; y++)
                {
                    statistics[0].lengthamng[x][y] = lengthamng[x][y];
                    statistics[0].lengthamng_outg[x][y] = lengthamng_outg[x][y];
                }
            }
            for (x = 0; x < args.int_total_nsam; x++)
            {
                for (y = 0; y < args.int_total_nsam + 1; y++)
                {
                    statistics[0].linefreq[x][y] = 0.;
                }
            }
            for (x = 0; x < args.npops; x++)
            {
                for (y = 0; y < args.int_total_nsam + 1; y++)
                {
                    statistics[0].popfreq[x][y] = 0.;
                }
            }

            sites_matrix = (long int *)calloc(4 * (length_seg + 1) * args.npops,
                                              sizeof(long int));
            jfd = (double **)calloc(args.npops, sizeof(double *));
            nfd = (int **)calloc(args.npops, sizeof(int *));

            for (x = 0; x < args.npops; x++)
            {
                jfd[x] = (double *)calloc(length_seg, sizeof(double));
                nfd[x] = (int *)calloc(length_seg, sizeof(int));
                for (y = 1; y < args.vint_perpop_nsam[x]; y++)
                {
                    statistics[0].freq[x][y] = 0;
                }
                statistics[0].mdsd[x] = -10000;
            }
            statistics[0].total_length = length_al;
            statistics[0].total_real_length = length_al_real;
            statistics[0].total_svratio = svratio;
            statistics[0].nmhits = nmhits;

            if(args.include_rsfs==0) {
                /*calculate statistics ------------------------------------------------ */
                if (calc_sxsfss(
                                args.npops,
                                args.vint_perpop_nsam,
                                matrix_pol,
                                matrix_pos,
                                length_seg,
                                statistics,
                                sites_matrix,
                                args.outgroup_presence,
                                args.force_outgroup) == 0)
                {
                    // fprintf(file_logerr,"\nError in calc_sxsfss function.1.\n");
                    log_error("Error in calc_sxsfss function.1.");
                    exit(1);
                }
                if (jointfreqdist(
                                  args.npops,
                                  args.vint_perpop_nsam,
                                  matrix_pol,
                                  matrix_pos,
                                  length_seg,
                                  statistics,
                                  sites_matrix,
                                  jfd,
                                  nfd,
                                  args.outgroup_presence,
                                  args.force_outgroup,
                                  args.output) == 0)
                {
                    // fprintf(file_logerr,"\nError in jointfreqdist function.1.\n");
                    log_error("Error in jointfreqdist function.1.");
                    exit(1);
                }
                if (calc_piwpiafst(
                                   flaghky,
                                   args.formatfile,
                                   args.npops,
                                   args.vint_perpop_nsam,
                                   matrix_pol,
                                   length_seg,
                                   statistics,
                                   matrix_sv,
                                   args.outgroup_presence,
                                   args.force_outgroup) == 0)
                {
                    // fprintf(file_logerr,"\nError in calc_piwpiafst function.1.\n");
                    log_error("Error in calc_piwpiafst function.1.");
                    exit(1);
                }
                if (calc_freqstats(args.npops, args.vint_perpop_nsam, matrix_pol, length_seg,
                                   statistics, args.outgroup_presence, args.force_outgroup, args.include_unknown, args.n_ccov, H1frq) == 0)
                {
                    // fprintf(file_logerr,"\nError in freqstats function.1.\n");
                    log_error("Error in freqstats function.1.");
                    exit(1);
                }
                if (H1frq && args.include_unknown == 0)
                {
                    if (calc_Toptimal_tests(args.npops,
                                            args.vint_perpop_nsam,
                                            statistics) == 0)
                    {
                        // fprintf(file_logerr,"\nError in calc_Toptimal_tests function.1.\n");
                        log_error("Error in calc_Toptimal_tests function.1.");
                        exit(1);
                    }
                }
                if (calcR2(args.npops,
                           args.vint_perpop_nsam,
                           matrix_pol,
                           length_seg,
                           statistics,
                           args.ploidy) == 0)
                {
                    // fprintf(file_logerr,"\nError in calc_R2 function.1.\n");
                    log_error("Error in calc_R2 function.1.");
                    exit(1);
                }
                if (calcR2p(args.npops,
                            args.vint_perpop_nsam,
                            matrix_pol,
                            length_seg,
                            statistics,
                            sum_sam,
                            args.r2i_ploidies,
                            args.outgroup_presence + args.force_outgroup) == 0)
                {
                    // fprintf(file_logerr,"\nError in calc_R2 function.1.\n");
                    log_error("Error in calc_R2 function.1.");
                    exit(1);
                }
            
                /* Calculate statistics for haplotypes */
                if (args.int_total_nsam < SAMPLE_LARGE)
                {
                    /*if(include_unknown==0) {*/
                    if (calc_mismatch(args.npops,
                                      args.vint_perpop_nsam,
                                      matrix_pol,
                                      length_seg,
                                      statistics,
                                      args.ploidy,
                                      args.outgroup_presence + args.force_outgroup) == 0)
                    {
                        // fprintf(file_logerr,"\nError in calc_mismatch function.1.\n");
                        log_error("Error in calc_mismatch function.1.");
                        exit(1);
                    }
                    /*}*/
                    if (args.include_unknown == 0 && args.ploidy[0] == '1')
                    {
                        if (calc_hwhafsth(args.npops,
                                          args.vint_perpop_nsam,
                                          matrix_pol,
                                          length_seg,
                                          statistics) == 0)
                        {
                            // fprintf(file_logerr,"\nError in calc_hwhafsth function.1.\n");
                            log_error("Error in calc_hwhafsth function.1.");
                            exit(1);
                        }
                        if (calcFs(args.npops,
                                   args.vint_perpop_nsam,
                                   statistics) == 0)
                        { /*need calc_freqstats and calc_hwhafsth to be calculated*/
                            // fprintf(file_logerr,"\nError in calc_Fs function.1.\n");
                            log_error("Error in calc_Fs function.1.");
                            exit(1);
                        }
                    }
                }
            
                if (file_output)
                {
                    /*
                     fprintf(file_output,"Done.\n");
                     fflush(file_output);
                     fprintf(file_logerr,"Done.\n");
                     */
                    fflush(stdout);
                }
            
                /* Un monton de cosas casi-repes que volvemos a utilizar! */
                /* PERMUTATION */
                /* Only if not MS format */
                if (args.niter && args.npops > 2 /*one is outgroup, forced or not*/ && args.include_unknown == 0)
                {
                    /*calloc pointers in structs for permutation test*/
                    if ((stats_iter = (struct stats *)calloc(1,
                                                             sizeof(struct stats))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter");
                        exit(1);
                    }
                
                    if ((stats_iter[0].piw = (double *)calloc(1 * args.npops,
                                                              sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].piw");
                        exit(1);
                    }
                
                    if ((stats_iter[0].pia = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                              sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].pia");
                        exit(1);
                    }
                    if ((stats_iter[0].piT = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                              sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].piT");
                        exit(1);
                    }
                
                    if ((stats_iter[0].piant = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                                sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].piant");
                        exit(1);
                    }
                    if ((stats_iter[0].piTnt = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                                sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].piTnt");
                        exit(1);
                    }
                
                    if ((stats_iter[0].fst = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                              sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].fst");
                        exit(1);
                    }
                
                    if ((stats_iter[0].piwHKY = (double *)calloc(1 * args.npops,
                                                                 sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].piwHKY");
                        exit(1);
                    }
                    if ((stats_iter[0].thetaTHKY = (double *)calloc(1 * args.npops,
                                                                    sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].thetaTHKY");
                        exit(1);
                    }
                
                    if ((stats_iter[0].piaHKY = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                                 sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].piaHKY");
                        exit(1);
                    }
                    if ((stats_iter[0].piTHKY = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                                 sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].piTHKY");
                        exit(1);
                    }
                
                    if ((stats_iter[0].fstHKY = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                                 sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].fstHKY");
                        exit(1);
                    }
                
                    if ((stats_iter[0].fst1all = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                                  sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].fst1all");
                        exit(1);
                    }
                
                    if ((stats_iter[0].hapw = (double *)calloc(1 * args.npops,
                                                               sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].hapw");
                        exit(1);
                    }
                
                    if ((stats_iter[0].hapa = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                               sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].hapa");
                        exit(1);
                    }
                
                    if ((stats_iter[0].hapT = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                               sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].hapT");
                        exit(1);
                    }
                
                    if ((stats_iter[0].fsth = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                               sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].fsth");
                        exit(1);
                    }
                
                    if ((stats_iter[0].fsth1all = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                                   sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].fsth1all");
                        exit(1);
                    }
                
                    if ((stats_iter[0].Gst = (double *)calloc((args.npops * (args.npops - 0)) / 2,
                                                              sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].Gst");
                        exit(1);
                    }
                
                    if ((stats_iter[0].K = (double *)calloc(1 * args.npops,
                                                            sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].K");
                        exit(1);
                    }
                
                    if ((stats_iter[0].KHKY = (double *)calloc(1 * args.npops,
                                                               sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].KHKY");
                        exit(1);
                    }
                
                    stats_iter[0].freq = 0;
                    if ((stats_iter[0].sv = (double ***)calloc(1 * args.npops,
                                                               sizeof(double **))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].sv");
                        exit(1);
                    }
                    if ((stats_iter[0].svT = (double ***)calloc(1 * args.npops,
                                                                sizeof(double **))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].svT");
                        exit(1);
                    }
                
                    if ((stats_iter[0].nhpop = (int *)calloc(1 * args.npops,
                                                             sizeof(int))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].nhpop");
                        exit(1);
                    }
                
                    if ((stats_iter[0].freqh = (long int **)calloc(1 * args.npops,
                                                                   sizeof(long int *))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].freqh");
                        exit(1);
                    }
                
                    if ((stats_iter[0].tcga = (double **)calloc(1 * args.npops,
                                                                sizeof(double *))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].tcga");
                        exit(1);
                    }
                
                    if ((stats_iter[0].length = (double *)calloc(1 * args.npops,
                                                                 sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].length");
                        exit(1);
                    }
                    if ((stats_iter[0].length2 = (double *)calloc(1 * args.npops,
                                                                  sizeof(double))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].length2");
                        exit(1);
                    }
                    if ((stats_iter[0].lengthamng = (double **)calloc(args.npops,
                                                                      sizeof(double *))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].lengthamng");
                        exit(1);
                    }
                    if ((stats_iter[0].lengthamng_outg = (double **)calloc(args.npops,
                                                                           sizeof(double *))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, stats_iter[0].lengthamng_outg");
                        exit(1);
                    }
                    for (x = 0; x < args.npops; x++)
                    {
                        if ((stats_iter[0].lengthamng[x] = (double *)calloc(args.npops,
                                                                            sizeof(double))) == 0)
                        {
                            // fprintf(file_logerr,"\n  Error allocating memory.");
                            log_fatal("Error allocating memory, stats_iter[0].lengthamng");
                            exit(1);
                        }
                        if ((stats_iter[0].lengthamng_outg[x] = (double *)calloc(args.npops,
                                                                                 sizeof(double))) == 0)
                        {
                            // fprintf(file_logerr,"\n  Error allocating memory.");
                            log_fatal("Error allocating memory, stats_iter[0].lengthamng_outg");
                            exit(1);
                        }
                        if ((stats_iter[0].tcga[x] = (double *)calloc(4,
                                                                      sizeof(double))) == 0)
                        {
                            // fprintf(file_logerr,"\n  Error allocating memory.");
                            log_fatal("Error allocating memory, stats_iter[0].tcga");
                            exit(1);
                        }
                    
                        if ((stats_iter[0].freqh[x] = (long int *)calloc(args.int_total_nsam,
                                                                         sizeof(long int))) == 0)
                        {
                            // fprintf(file_logerr,"\n  Error allocating memory.");
                            log_fatal("Error allocating memory, stats_iter[0].freqh");
                            exit(1);
                        }
                    
                        if ((stats_iter[0].sv[x] = (double **)calloc(1 * args.npops,
                                                                     sizeof(double *))) == 0)
                        {
                            // fprintf(file_logerr,"\n  Error allocating memory.");
                            log_fatal("Error allocating memory, stats_iter[0].sv");
                            exit(1);
                        }
                        if ((stats_iter[0].svT[x] = (double **)calloc(1 * args.npops,
                                                                      sizeof(double *))) == 0)
                        {
                            // fprintf(file_logerr,"\n  Error allocating memory.");
                            log_fatal("Error allocating memory, stats_iter[0].svT");
                            exit(1);
                        }
                    
                        for (y = 0; y < args.npops; y++)
                        {
                            if ((stats_iter[0].sv[x][y] = (double *)calloc(2,
                                                                           sizeof(double))) == 0)
                            {
                                // fprintf(file_logerr,"\n  Error allocating memory.");
                                log_fatal("Error allocating memory, stats_iter[0].sv");
                                exit(1);
                            }
                            if ((stats_iter[0].svT[x][y] = (double *)calloc(2,
                                                                            sizeof(double))) == 0)
                            {
                                // fprintf(file_logerr,"\n  Error allocating memory.");
                                log_fatal("Error allocating memory, stats_iter[0].svT");
                                exit(1);
                            }
                        }
                    }
                    for (x = 0; x < args.npops; x++)
                    {
                        for (w = 0; w < 4; w++)
                        {
                            stats_iter[0].tcga[x][w] = statistics[0].tcga[x][w];
                        }
                    }
                    stats_iter[0].total_length = length_al;
                    stats_iter[0].total_svratio = svratio;
                    stats_iter[0].nmhits = nmhits;
                
                    if ((piter = (struct probs *)calloc(1,
                                                        sizeof(struct probs))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, piter");
                        exit(1);
                    }
                
                    if ((piter[0].i1 = (long int *)calloc(1 * args.npops,
                                                          sizeof(long int))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, piter[0].i1");
                        exit(1);
                    }
                
                    if ((piter[0].ih1 = (long int *)calloc(1 * args.npops,
                                                           sizeof(long int))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, piter[0].ih1");
                        exit(1);
                    }
                
                    if ((piter[0].i = (long int *)calloc((args.npops * (args.npops - 0)) / 2,
                                                         sizeof(long int))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, piter[0].i");
                        exit(1);
                    }
                
                    if ((piter[0].ih = (long int *)calloc((args.npops * (args.npops - 0)) / 2,
                                                          sizeof(long int))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, piter[0].ih");
                        exit(1);
                    }
                
                    if ((piter[0].igh = (long int *)calloc((args.npops * (args.npops - 0)) / 2,
                                                           sizeof(long int))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, piter[0].igh");
                        exit(1);
                    }
                
                    if ((piter[0].niteri1 = (long int *)calloc(1 * args.npops,
                                                               sizeof(long int))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, piter[0].niteri1");
                        exit(1);
                    }
                
                    if ((piter[0].niterih1 = (long int *)calloc(1 * args.npops,
                                                                sizeof(long int))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, piter[0].niterih1");
                        exit(1);
                    }
                
                    if ((piter[0].niteri = (long int *)calloc((args.npops * (args.npops - 0)) / 2,
                                                              sizeof(long int))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, piter[0].niteri");
                        exit(1);
                    }
                
                    if ((piter[0].niterih = (long int *)calloc((args.npops * (args.npops - 0)) / 2,
                                                               sizeof(long int))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, piter[0].niterih");
                        exit(1);
                    }
                
                    if ((piter[0].niterigh = (long int *)calloc((args.npops * (args.npops - 0)) / 2,
                                                                sizeof(long int))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, piter[0].niterigh");
                        exit(1);
                    }
                
                    piter[0].iall = 0;
                    piter[0].ihall = 0;
                    piter[0].ighall = 0;
                    piter[0].niteriall = 0;
                    piter[0].niterihall = 0;
                    piter[0].niterighall = 0;
                
                    if ((matrix_perm = (char *)calloc(length_seg * args.int_total_nsam,
                                                      sizeof(char))) == 0)
                    {
                        // fprintf(file_logerr,"\n  Error allocating memory.");
                        log_fatal("Error allocating memory, matrix_perm");
                        exit(1);
                    }
                
                    // if(file_output && (args.output == 0 || args.output == 10)) {
                    if (file_output && (args.output == OUTPUT_EXTENDED || args.output == OUTPUT_FULL_EXTENDED))
                    {
                        /*fprintf(file_output,"Calculating permutation tests...\n");*/
                        fflush(file_output);
                        /*
                         fprintf(file_logerr,"Calculating permutation test...\n");
                         fflush(stdout);
                         */
                    }
                
                    /*permute samples in pops:*/
                    if (args.npops > 2)
                    {
                        /* TODO: OJO que necesita inicializar las poblaciones de otra manera, segun ...*/
                        for (i = 0; i < args.niter; i++)
                        {
                            /*assign the (2) groups that will be included in permutation test*/
                            nsam2[0] = args.vint_perpop_nsam[0];
                            nsam2[1] = args.int_total_nsam - args.vint_perpop_nsam[args.npops - 1] - args.vint_perpop_nsam[0];
                            /*nsam2[2] = vint_perpop_nsam[npops-1];*/
                        
                            psam2[0] = 0;
                            psam2[1] = args.vint_perpop_nsam[0];
                            /*psam2[2] = psam2[1] + nsam2[1];*/
                        
                            if (permute(matrix_pol,
                                        length_seg,
                                        args.int_total_nsam,
                                        matrix_perm,
                                        nsam2,
                                        psam2,
                                        args.npops,
                                        args.vint_perpop_nsam[args.npops - 1],
                                        args.int_total_nsam - args
                                        .vint_perpop_nsam[args.npops - 1]) == 0)
                            {
                                // fprintf(file_logerr,"\nError in permute function.\n");
                                log_error("Error in permute function.");
                                exit(1);
                            }
                            /*HERE WE ASSUME THAT THE LENGTH SIZE ARE THE SAME THAN FOR ORIGINAL SAMPLES!! (OK FOR NO MISSING)*/
                            for (x = 0; x < args.npops; x++)
                            {
                                stats_iter[0]
                                    .length[x] = statistics[0].length[x];
                                stats_iter[0]
                                    .length2[x] = statistics[0].length2[x];
                                for (yy = 0; yy < args.npops; yy++)
                                {
                                    stats_iter[0]
                                        .lengthamng[x][yy] = statistics[0].lengthamng[x][yy];
                                    stats_iter[0]
                                        .lengthamng_outg[x][yy] = statistics[0].lengthamng_outg[x][yy];
                                }
                            }
                        
                            if (calc_piwpiafst(0, 0, args.npops, args.vint_perpop_nsam, matrix_perm,
                                               length_seg, stats_iter, matrix_sv,
                                               args.outgroup_presence, args.force_outgroup) == 0)
                            {
                                // fprintf(file_logerr,"\nError in calc_piwpiafst function.2.\n");
                                log_error("Error in calc_piwpiafst function.2.");
                                exit(1);
                            }
                            if (args.ploidy[0] == '1')
                            {
                                if (calc_hwhafsth(args.npops,
                                                  args.vint_perpop_nsam,
                                                  matrix_perm,
                                                  length_seg,
                                                  stats_iter) == 0)
                                {
                                    // fprintf(file_logerr,"\nError in calc_hwhafsth function.2.\n");
                                    log_error("Error in calc_hwhafsth function.2.");
                                    exit(1);
                                }
                            }
                        
                            z = 0;
                            for (x = 0; x < args.npops - 1; x++)
                            {
                                if (stats_iter[0].fst1all[x] != -10000 && statistics[0].fst1all[x] != -10000)
                                {
                                    if (statistics[0].fst1all[x] <= stats_iter[0].fst1all[x])
                                        piter[0].i1[x]++;
                                    piter[0].niteri1[x]++;
                                    /*fprintf(file_output,"%.3f\t",stats_iter[0].fst1all[x]);*/
                                }
                                if (args.ploidy[0] == '1')
                                {
                                    if (stats_iter[0].fsth1all[x] != -10000 && statistics[0].fsth1all[x] != -10000)
                                    {
                                        if (statistics[0].fsth1all[x] <= stats_iter[0].fsth1all[x])
                                            piter[0].ih1[x]++;
                                        piter[0].niterih1[x]++;
                                    }
                                }
                            }
                        
                            /*fprintf(file_output,"\n");*/
                            if (stats_iter[0].fstALL != -10000 && statistics[0].fstALL != -10000)
                            {
                                if (statistics[0].fstALL <= stats_iter[0].fstALL)
                                    piter[0].iall++;
                                piter[0].niteriall++;
                            }
                        
                            if (args.ploidy[0] == '1')
                            {
                                if (stats_iter[0].fsthALL != -10000 && statistics[0].fsthALL != -10000)
                                {
                                    if (statistics[0].fsthALL <= stats_iter[0].fsthALL)
                                        piter[0].ihall++;
                                    piter[0].niterihall++;
                                }
                                if (stats_iter[0].GstALL != -10000 && statistics[0].GstALL != -10000)
                                {
                                    if (statistics[0].GstALL <= stats_iter[0].GstALL)
                                        piter[0].ighall++;
                                    piter[0].niterighall++;
                                }
                            }
                        }
                        // if(file_output && (args.output == 0 || args.output == 10)) {
                        if (file_output && (args.output == OUTPUT_EXTENDED || args.output == OUTPUT_FULL_EXTENDED))
                        {
                            /*fprintf(file_output,"Permutation test one vs all Done.\n");*/
                            fflush(file_output);
                            /*
                             fprintf(file_logerr,"Permutation test one vs all Done.\n");
                             fflush(stdout);
                             */
                        }
                        /* permute pairs of pops */
                        for (i = 0; i < args.niter; i++)
                        {
                            z = 0;
                            psam2[0] = 0;
                            psam2[1] = 0;
                        
                            for (x = 0; x < args.npops - 1; x++)
                            {
                                nsam2[0] = args.vint_perpop_nsam[x];
                                psam2[1] = psam2[0] + args.vint_perpop_nsam[x];
                                for (y = x + 1; y < args.npops - 0; y++)
                                {
                                    nsam2[1] = args.vint_perpop_nsam[y];
                                    psam2[2] = psam2[1] + args.vint_perpop_nsam[y];
                                    nsam2[2] = args.vint_perpop_nsam[args.npops - 1];
                                    if (permute(matrix_pol,
                                                length_seg,
                                                args.int_total_nsam,
                                                matrix_perm,
                                                nsam2,
                                                psam2,
                                                args.npops,
                                                args
                                                .vint_perpop_nsam[args.npops - 1],
                                                args.int_total_nsam - args
                                                .vint_perpop_nsam[args.npops - 1]) == 0)
                                    {
                                        // fprintf(file_logerr,"\nError in permute function.\n");
                                        log_error("Error in permute function.");
                                        exit(1);
                                    }
                                    /*HERE WE ASSUME THAT THE LENGTH SIZE ARE THE SAME THAN FOR ORIGINAL SAMPLES!! (OK FOR NO MISSING)*/
                                    stats_iter[0]
                                        .length[0] = statistics[0].length[x];
                                    stats_iter[0]
                                        .length2[0] = statistics[0].length2[x];
                                    stats_iter[0]
                                        .length[1] = statistics[0].length[y];
                                    stats_iter[0]
                                        .length2[1] = statistics[0].length2[y];
                                    stats_iter[0]
                                        .lengthamng[0][1] = statistics[0].lengthamng[x][y];
                                    stats_iter[0]
                                        .lengthamng_outg[0][1] = statistics[0].lengthamng_outg[x][y];
                                
                                    if (calc_piwpiafst(0,
                                                       0,
                                                       2 + 1 /*outg*/,
                                                       nsam2,
                                                       matrix_perm,
                                                       length_seg,
                                                       stats_iter,
                                                       matrix_sv,
                                                       args.outgroup_presence,
                                                       args.force_outgroup) == 0)
                                    {
                                        // fprintf(file_logerr,"\nError in calc_piwpiafst function.2.\n");
                                        log_error("Error in calc_piwpiafst function.2.");
                                        exit(1);
                                    }
                                    if (args.int_total_nsam < SAMPLE_LARGE)
                                    {
                                        if (args.ploidy[0] == '1')
                                        {
                                            if (calc_hwhafsth(2 + 1,
                                                              nsam2,
                                                              matrix_perm,
                                                              length_seg,
                                                              stats_iter) == 0)
                                            {
                                                // fprintf(file_logerr,"\nError in calc_hwhafsth function.2.\n");
                                                log_error("Error in calc_hwhafsth function.2.");
                                                exit(1);
                                            }
                                        }
                                    }
                                    if (stats_iter[0].fst[0] != -10000 && statistics[0].fst[z] != -10000)
                                    {
                                        if (statistics[0].fst[z] <= stats_iter[0].fst[0])
                                            piter[0].i[z]++;
                                        piter[0].niteri[z]++;
                                    }
                                    if (args.ploidy[0] == '1')
                                    {
                                        if (stats_iter[0].fsth[0] != -10000 && statistics[0].fsth[z] != -10000)
                                        {
                                            if (statistics[0].fsth[z] <= stats_iter[0].fsth[0])
                                                piter[0].ih[z]++;
                                            piter[0].niterih[z]++;
                                        }
                                        if (stats_iter[0].Gst[0] != -10000 && statistics[0].Gst[z] != -10000)
                                        {
                                            if (statistics[0].Gst[z] <= stats_iter[0].Gst[0])
                                                piter[0].igh[z]++;
                                            piter[0].niterigh[z]++;
                                        }
                                    }
                                    psam2[1] += args.vint_perpop_nsam[y]; /*!outgroup_presence at loops*/
                                    z++;
                                }
                                psam2[0] += args.vint_perpop_nsam[x];
                            }
                        }
                    }
                    if (file_output)
                    {
                        /*
                         fprintf(file_output,"Done.\n");
                         fflush(file_output);
                         printf("Done.\n");
                         */
                        fflush(stdout);
                    }
                }
            }
            else {
                //functions to calculate rSFS statistics
                if(args.include_rsfs) {
                    if (calc_rfreqstats(args.npops, args.vint_perpop_nsam, matrix_pol, length_seg,
                                       statistics, rw, args.outgroup_presence, args.force_outgroup, args.include_unknown) == 0)
                    {
                        // fprintf(file_logerr,"\nError in rfreqstats function.1.\n");
                        log_error("Error in rfreqstats function.1.");
                        exit(1);
                    }
                }
            }

            /*print results*/
            /* TODO: mas elegante el tema de force_outgroup */
            if (
                print_output(
                             //argc,
                             //args.npops,
                             //args.vint_perpop_nsam,
                             file_output,
                             &file_output_gz,
                             //args.file_in,
                             //args.file_out,
                             //args.gfffiles,
                             //args.file_GFF, 
                             //args.subset_positions, 
                             //args.code_name,
                             //args.genetic_code, 
                             //args.length, 
                             length_seg, 
                             length_al,
                             length_al_real, 
                             statistics, 
                             piter, 
                             //args.niter,
                             sites_matrix, 
                             args.ploidy, 
                             svratio, 
                             missratio,
                             //args.include_unknown, 
                             matrix_pos, 
                             jfd, 
                             nfd, 
                             //args.output,
                             H1frq, 
                             H0frq, 
                             //args.nseed, 
                             //args.file_H1f, 
                             //args.file_H0f, 
                             vector_priors,
                             npriors, 
                             //args.formatfile,
                             //args.outgroup_presence + args.force_outgroup, 
                             //args.force_outgroup, 
                             //args.freq_missing_ms,
                             nsites1_pop, 
                             nsites1_pop_outg, 
                             nsites2_pop, 
                             nsites2_pop_outg, 
                             nsites3_pop, 
                             nsites3_pop_outg,
                             rnsites1_pop,
                             rnsites1_pop_outg,
                             rnsites2_pop,
                             rnsites2_pop_outg,
                             rnsites3_pop,
                             rnsites3_pop_outg,
                             li + 1,
                             matrix_pol, 
                             //args.r2i_ploidies, 
                             matrix_pol_tcga, 
                             chr_name,
                             &args) == 0)
            {

                // fprintf(file_logerr,"\nSorry. Error in printing function.\n");
                log_error("Sorry. Error in printing function.");
                exit(1);
            }
            /* TODO: Check cleaning the house */
            log_trace("Cleaning matrix_perm");
            free(matrix_pol);
            // if (!(args.formatfile == 1 || args.formatfile == 2))
            if (!(args.formatfile == MS_FORMAT || args.formatfile == MS_X_FORMAT)) {
                log_trace("Cleaning matrix_pol_tcga");
                free(matrix_pol_tcga);
            }
            log_trace("Cleaning matrix_freq");
            free(matrix_freq);
            log_trace("Cleaning matrix_pos");
            free(matrix_pos);
            /*free(matrix_GC);*/
            log_trace("Cleaning matrix_sv");
            free(matrix_sv);
            log_trace("Cleaning sites_matrix");
            free(sites_matrix);

            for (x = 0; x < args.npops; x++)
            {
                log_trace("Cleaning jfd[%d]",x);
                free(jfd[x]);
                log_trace("Cleaning nfd[%d]",x);
                free(nfd[x]);
            }
            log_trace("Cleaning jfd");
            free(jfd);
            log_trace("Cleaning nfd");
            free(nfd);
            // if (/*include_unknown && */ args.file_mas[0] == '-' && (args.formatfile == 1 || args.formatfile == 2))
            if (/*include_unknown && */ args.file_mas[0] == '-' && (args.formatfile == MS_FORMAT || args.formatfile == MS_X_FORMAT))
            {
                log_trace("Cleaning sum_sam_mask");
                free(sum_sam_mask);
            }
            li++;
        } // end of while( li < args.niterdata)
        // if (args.formatfile == 3)
        if (args.formatfile == TFA_FORMAT)
        {
            log_info("Processing scaffold %s Done", chr_name);
            if (file_wcoor) {
                log_trace("cleaning wgenes");
                free(wgenes);
            }
        }
    }
    // if (file_ws)
    //   fzclose(file_ws, &file_ws_gz);
  

    if(wtfasta) {
        log_trace("Cleaning tfasta weight file");
        close_wtfasta_file(wtfasta);
    }
    
  
    // if (/*include_unknown && */ args.file_mas[0] != '-' && (args.formatfile == 1 || args.formatfile == 2))
    if (/*include_unknown && */ args.file_mas[0] != '-' && (args.formatfile == MS_FORMAT || args.formatfile == MS_X_FORMAT))
    {
        log_trace("Clearning sum_sam_mask");
        free(sum_sam_mask);
    }
    if (args.formatfile > 0)
    {
        log_trace("Cleaning file_input and file_input_gz");
        fzclose(file_input, &file_input_gz);
    }
    log_trace("cleaning file_output");
    if (file_output)
        fclose(file_output);

    for (x = 0; x < nscaffolds; x++)
    {
        log_trace("cleaning chr_name_array[%d]",x);
        free(chr_name_array[x]);
        log_trace("cleaning chr_length_array[%d]",x);
        free(chr_length_array[x]);
    }
  
    log_trace("cleaning chr_name_array");
    free(chr_name_array);
    log_trace("cleaning chr_length_array");
    free(chr_length_array);

    log_trace("cleaning args.chr_name");
    free(args.sort_nsam);

    log_trace("cleaning nsites1_pop");
    free(nsites1_pop);
  
    log_trace("cleaning nsites2_pop");
    free(nsites2_pop);

    log_trace("cleaning nsites3_pop");
    free(nsites3_pop);

    log_trace("cleaning nsites1_pop_outg");
    free(nsites1_pop_outg);

    log_trace("cleaning nsites2_pop_outg");
    free(nsites2_pop_outg);

    log_trace("cleaning nsites3_pop_outg");
    free(nsites3_pop_outg);

    log_trace("cleaning rnsites1_pop");
    free(rnsites1_pop);
  
    log_trace("cleaning rnsites2_pop");
    free(rnsites2_pop);

    log_trace("cleaning rnsites3_pop");
    free(rnsites3_pop);

    log_trace("cleaning rnsites1_pop_outg");
    free(rnsites1_pop_outg);

    log_trace("cleaning rnsites2_pop_outg");
    free(rnsites2_pop_outg);

    log_trace("cleaning rnsites3_pop_outg");
    free(rnsites3_pop_outg);

    log_trace("cleaning anx");
    free(anx);

    log_trace("cleaning bnx");
    free(bnx);

    log_trace("cleaning anxo");
    free(anxo);

    log_trace("cleaning bnxo");
    free(bnxo);
    for (x = 0; x < args.npops; x++)
    {
        log_trace("cleaning lengthamng[%d]",x);
        free(lengthamng[x]);
        log_trace("cleaning lengthamng_outg[%d]",x);
        free(lengthamng_outg[x]);
    }
    log_trace("cleaning lengthamng");
    free(lengthamng);
    log_trace("cleaning lengthamng_outg");
    free(lengthamng_outg);

    log_trace("cleaning args.vint_perpop_nsam");
    free(args.vint_perpop_nsam);

    if (H1frq && args.include_rsfs == 0)
    {
        for (x = 0; x < npf; x++)
        {
            log_trace("cleaning freqspH1[%d]",x);
            free(freqspH1[x]);
        }
        log_trace("cleaning freqspH1");
        free(freqspH1);
        log_trace("cleaning thetaH1");
        free(thetaH1);
        for (x = 0; x < npf; x++){
            log_trace("cleaning freqspH0[%d]",x);
            free(freqspH0[x]);
        }
        log_trace("cleaning freqspH0");
        free(freqspH0);
        log_trace("cleaning thetaH0");
        free(thetaH0);
    }

    log_trace("cleaning statistics[0].Sanc");
    free(statistics[0].Sanc);
    log_trace("cleaning statistics[0].piw");
    free(statistics[0].piw);
    log_trace("cleaning statistics[0].piwHKY");
    free(statistics[0].piwHKY);
    log_trace("cleaning statistics[0].hapw");
    free(statistics[0].hapw);

    if(args.include_rsfs == 0) {
        log_trace("cleaning statistics[0].pia");
        free(statistics[0].pia);
        log_trace("cleaning statistics[0].piT");
        free(statistics[0].piT);
        log_trace("cleaning statistics[0].piant");
        free(statistics[0].piant);
        log_trace("cleaning statistics[0].piTnt");
        free(statistics[0].piTnt);
        log_trace("cleaning statistics[0].fst");
        free(statistics[0].fst);
        log_trace("cleaning statistics[0].piaHKY");
        free(statistics[0].piaHKY);
        log_trace("cleaning statistics[0].piTHKY");
        free(statistics[0].piTHKY);
        log_trace("cleaning statistics[0].fstHKY");
        free(statistics[0].fstHKY);
        log_trace("cleaning statistics[0].fst1all");
        free(statistics[0].fst1all);
        log_trace("cleaning statistics[0].hapa");
        free(statistics[0].hapa);
        log_trace("cleaning statistics[0].hapT");
        free(statistics[0].hapT);
        log_trace("cleaning statistics[0].fsth");
        free(statistics[0].fsth);
        log_trace("cleaning statistics[0].fsth1all");
        free(statistics[0].fsth1all);
        log_trace("cleaning statistics[0].Gst");
        free(statistics[0].Gst);
    }
    
    log_trace("cleaning statistics[0].S");
    free(statistics[0].S);
    log_trace("cleaning statistics[0].thetaS");
    free(statistics[0].thetaS);
    log_trace("cleaning statistics[0].thetaT");
    free(statistics[0].thetaT);
    log_trace("cleaning statistics[0].So");
    free(statistics[0].So);
    log_trace("cleaning statistics[0].thetaSo");
    free(statistics[0].thetaSo);
    log_trace("cleaning statistics[0].thetaTo");
    free(statistics[0].thetaTo);
    log_trace("cleaning statistics[0].thetaTHKY");
    free(statistics[0].thetaTHKY);
    log_trace("cleaning statistics[0].thetaFL");
    free(statistics[0].thetaFL);
    log_trace("cleaning statistics[0].thetaFW");
    free(statistics[0].thetaFW);
    log_trace("cleaning statistics[0].thetaL");
    free(statistics[0].thetaL);
    log_trace("cleaning statistics[0].thetaSA");
    free(statistics[0].thetaSA);
    log_trace("cleaning statistics[0].thetaTA");
    free(statistics[0].thetaTA);
    log_trace("cleaning statistics[0].K");
    free(statistics[0].K);
    log_trace("cleaning statistics[0].KHKY");
    free(statistics[0].KHKY);
    log_trace("cleaning statistics[0].Dtaj");
    free(statistics[0].Dtaj);
    log_trace("cleaning statistics[0].Dfl");
    free(statistics[0].Dfl);
    log_trace("cleaning statistics[0].Ffl");
    free(statistics[0].Ffl);
    log_trace("cleaning statistics[0].Hnfw");
    free(statistics[0].Hnfw);
    log_trace("cleaning statistics[0].Ez");
    free(statistics[0].Ez);
    log_trace("cleaning statistics[0].Yach");
    free(statistics[0].Yach);
    log_trace("cleaning statistics[0].FH");
    free(statistics[0].FH);
    log_trace("cleaning statistics[0].R2");
    free(statistics[0].R2);
    log_trace("cleaning statistics[0].Fs");  
    free(statistics[0].Fs);
    log_trace("cleaning statistics[0].nhpop");
    free(statistics[0].nhpop);
    log_trace("cleaning statistics[0].length");
    free(statistics[0].length);
    log_trace("cleaning statistics[0].total_length");
    free(statistics[0].total_tcga);
    log_trace("cleaning statistics[0].ToH0_ii");
    free(statistics[0].ToH0_ii);
    log_trace("cleaning statistics[0].To_ii");
    free(statistics[0].To_ii);
    log_trace("cleaning statistics[0].To_00");
    free(statistics[0].To_00);
    log_trace("cleaning statistics[0].To_i0");
    free(statistics[0].To_i0);
    log_trace("cleaning statistics[0].ToH0_00");
    free(statistics[0].ToH0_00);
    log_trace("cleaning statistics[0].To_Qc_ii");
    free(statistics[0].To_Qc_ii);
    log_trace("cleaning statistics[0].To_Qw_ii");
    free(statistics[0].To_Qw_ii);
    log_trace("cleaning statistics[0].To_Lc_ii");
    free(statistics[0].To_Lc_ii);


    log_trace("cleaning statistics[0].Rm");
    free(statistics[0].Rm);
    log_trace("cleaning statistics[0].ZnA");
    free(statistics[0].ZnA);

    log_trace("cleaning statistics[0].mdsd");
    free(statistics[0].mdsd);

    log_trace("cleaning statistics[0].mdg1");
    free(statistics[0].mdg1);

    log_trace("cleaning statistics[0].mdg2");
    free(statistics[0].mdg2);

    log_trace("cleaning statistics[0].anx");
    free(statistics[0].anx);

    log_trace("cleaning statistics[0].bnx");
    free(statistics[0].bnx);

    log_trace("cleaning statistics[0].anxo");
    free(statistics[0].anxo);

    log_trace("cleaning statistics[0].bnxo");
    free(statistics[0].bnxo);

    for (x = 0; x < args.npops; x++)
    {
        log_trace("cleaning statistics[0].freq[%d]",x);
        free(statistics[0].freq[x]);
        log_trace("cleaning statistics[0].freqh[%d]",x);
        free(statistics[0].freqh[x]);
        log_trace("cleaning statistics[0].tcga[%d]",x);
        free(statistics[0].tcga[x]);
        log_trace("cleaning statistics[0].mdw[%d]",x);
        free(statistics[0].mdw[x]);
        log_trace("cleaning statistics[0].lengthamng[%d]",x);
        free(statistics[0].lengthamng[x]);
        log_trace("cleaning statistics[0].lengthamng_outg[%d]",x);
        free(statistics[0].lengthamng_outg[x]);
        for (y = 0; y < args.npops; y++)
        {
            log_trace("cleaning statistics[0].sv[%d][%d]",x,y);
            free(statistics[0].sv[x][y]);
            log_trace("cleaning statistics[0].svT[%d][%d]",x,y);
            free(statistics[0].svT[x][y]);
        }
        log_trace("cleaning statistics[0].sv[%d]",x);
        free(statistics[0].sv[x]);
        log_trace("cleaning statistics[0].svT[%d]",x);
        free(statistics[0].svT[x]);
    }

    log_trace("cleaning statistics[0].length2");
    free(statistics[0].length2);
    log_trace("cleaning statistics[0].lengthamng");
    free(statistics[0].lengthamng);
    log_trace("cleaning statistics[0].lengthamng_outg");
    free(statistics[0].lengthamng_outg);

    for (x = 0; x < args.int_total_nsam; x++)
    {
        log_trace("cleaning statistics[0].linefreq[%d]",x);
        free(statistics[0].linefreq[x]);
    }
    for (x = 0; x < args.npops; x++)
    {
        log_trace("cleaning statistics[0].popfreq[%d]",x);
        free(statistics[0].popfreq[x]);
    }
    for (x = 0; x < args.r2i_ploidies[0]; x++)
    {
        log_trace("cleaning statistics[0].R2p[%d]",x);
        free(statistics[0].R2p[x]);
    }

    log_trace("cleaning statistics[0].R2p");
    free(statistics[0].R2p);

    log_trace("cleaning statistics[0].freq");
    free(statistics[0].freq);

    log_trace("cleaning statistics[0].freqh");
    free(statistics[0].freqh);

    log_trace("cleaning statistics[0].tcga");
    free(statistics[0].tcga);

    log_trace("cleaning statistics[0].sv");
    free(statistics[0].sv);

    log_trace("cleaning statistics[0].svT");
    free(statistics[0].svT);

    log_trace("cleaning statistics[0].mdw");
    free(statistics[0].mdw);

    log_trace("cleaning statistics[0].linefreq");
    free(statistics[0].linefreq);

    log_trace("cleaning statistics[0].popfreq");
    free(statistics[0].popfreq);

    log_trace("cleaning statistics");
    free(statistics);

    log_trace("cleaning sum_sam");
    free(sum_sam);
    for (x = 0; x < args.int_total_nsam + (!args.outgroup_presence); x++) {
        log_trace("cleaning tcga[%d]",x);  
        free(tcga[x]);
    }
    log_trace("cleaning tcga");
    free(tcga);
    // if (!(args.formatfile == 0 || (args.formatfile == 3 && ((args.slide == 0 && args.window == 0) && args.file_Wcoord[0] == '\0'))))
    if (!(args.formatfile == FASTA_FORMAT || (args.formatfile == TFA_FORMAT && ((args.slide == 0 && args.window == 0) && args.file_Wcoord[0] == '\0'))))

    {
        log_trace("cleaning matrix_mask");
        free(matrix_mask);
        log_trace("cleaning vector_mask");
        free(vector_mask);
    }

    if (args.niter && args.npops > 2)
    {
        log_trace("cleaning matrix_perm");
        free(matrix_perm);

        log_trace("cleaning stats_iter[0].piw");
        free(stats_iter[0].piw);

        log_trace("cleaning stats_iter[0].pia");
        free(stats_iter[0].pia);

        log_trace("cleaning stats_iter[0].piT");
        free(stats_iter[0].piT);

        log_trace("cleaning stats_iter[0].piant");
        free(stats_iter[0].piant);

        log_trace("cleaning stats_iter[0].piTnt");
        free(stats_iter[0].piTnt);

        log_trace("cleaning stats_iter[0].fst");
        free(stats_iter[0].fst);

        log_trace("cleaning stats_iter[0].piwHKY");
        free(stats_iter[0].piwHKY);

        log_trace("cleaning stats_iter[0].thetaTHKY");
        free(stats_iter[0].thetaTHKY);

        log_trace("cleaning stats_iter[0].piaHKY");
        free(stats_iter[0].piaHKY);

        log_trace("cleaning stats_iter[0].piTHKY");
        free(stats_iter[0].piTHKY);

        log_trace("cleaning stats_iter[0].fstHKY");
        free(stats_iter[0].fstHKY);

        log_trace("cleaning stats_iter[0].fst1all");
        free(stats_iter[0].fst1all);

        log_trace("cleaning stats_iter[0].hapw");
        free(stats_iter[0].hapw);

        log_trace("cleaning stats_iter[0].hapa");
        free(stats_iter[0].hapa);

        log_trace("cleaning stats_iter[0].hapT");
        free(stats_iter[0].hapT);

        log_trace("cleaning stats_iter[0].fsth");
        free(stats_iter[0].fsth);

        log_trace("cleaning stats_iter[0].fsth1all");
        free(stats_iter[0].fsth1all);

        log_trace("cleaning stats_iter[0].Gst");
        free(stats_iter[0].Gst);

        log_trace("cleaning stats_iter[0].K");
        free(stats_iter[0].K);

        log_trace("cleaning stats_iter[0].KHKY");
        free(stats_iter[0].KHKY);

        log_trace("cleaning stats_iter[0].length");
        free(stats_iter[0].length);

        log_trace("cleaning stats_iter[0].length2");
        free(stats_iter[0].length2);

        for (x = 0; x < args.npops; x++)
        {
            log_trace("cleaning stats_iter[0].freqh[%d]",x);
            free(stats_iter[0].freqh[x]);
            log_trace("cleaning stats_iter[0].tcga[%d]",x);
            free(stats_iter[0].tcga[x]);
            for (y = 0; y < args.npops; y++)
            {
                log_trace("cleaning stats_iter[0].sv[%d][%d]",x,y);
                free(stats_iter[0].sv[x][y]);
                log_trace("cleaning stats_iter[0].svT[%d][%d]",x,y);
                free(stats_iter[0].svT[x][y]);
            }
            log_trace("cleaning stats_iter[0].sv[%d]",x);
            free(stats_iter[0].sv[x]);
            log_trace("cleaning stats_iter[0].svT[%d]",x);
            free(stats_iter[0].svT[x]);
            log_trace("cleaning stats_iter[0].lengthamng[%d]",x);
            free(stats_iter[0].lengthamng[x]);
            log_trace("cleaning stats_iter[0].lengthamng_outg[%d]",x);
            free(stats_iter[0].lengthamng_outg[x]);
        }
        log_trace("cleaning stats_iter[0].freqh");
        free(stats_iter[0].freqh);
        log_trace("cleaning stats_iter[0].nhpop");
        free(stats_iter[0].nhpop);
        log_trace("cleaning stats_iter[0].tcga");
        free(stats_iter[0].tcga);
        log_trace("cleaning stats_iter[0].sv");
        free(stats_iter[0].sv);
        log_trace("cleaning stats_iter[0].svT");
        free(stats_iter[0].svT);
        log_trace("cleaning stats_iter[0].lengthamng");
        free(stats_iter[0].lengthamng);
        log_trace("cleaning stats_iter[0].lengthamng_outg");
        free(stats_iter[0].lengthamng_outg);

        log_trace("cleaning stats_iter");
        free(stats_iter);

        log_trace("cleaning piter[0].i1");
        free(piter[0].i1);

        log_trace("cleaning piter[0].ih1");
        free(piter[0].ih1);

        log_trace("cleaning piter[0].i");
        free(piter[0].i);

        log_trace("cleaning piter[0].ih");
        free(piter[0].ih);

        log_trace("cleaning piter[0].igh");
        free(piter[0].igh);

        log_trace("cleaning piter[0].niteri1");
        free(piter[0].niteri1);
        log_trace("cleaning piter[0].niterih1");
        free(piter[0].niterih1);

        log_trace("cleaning piter[0].niteri");
        free(piter[0].niteri);

        log_trace("cleaning piter[0].niterih");
        free(piter[0].niterih);

        log_trace("cleaning piter[0].niterigh");
        free(piter[0].niterigh);

        if (npriors) {
            log_trace("cleaning vector_priors");
            free(vector_priors);
        }
        log_trace("cleaning piter");
        free(piter);
    }
    //free rSFS stats
    if(args.include_rsfs) {
        log_trace("cleaning rSFS weights");
        for(x=0;x<args.int_total_nsam;x++) {
            free(rw->ww[x]);
            free(rw->wt[x]);
            free(rw->wfw[x]);
            free(rw->wfl[x]);
            free(rw->wl[x]);
            free(rw->wwA[x]);
            free(rw->wtA[x]);
            free(rw->wfwA[x]);
            free(rw->wflA[x]);
            free(rw->wlA[x]);
            free(rw->wphi_i[x]);
            for(z=0;z<x+1;z++)
                free(rw->wpsi_ij[x][z]);
            free(rw->wpsi_ij[x]);
        }
        free(rw->ww);
        free(rw->wt);
        free(rw->wfw);
        free(rw->wfl);
        free(rw->wl);
        free(rw->wwA);
        free(rw->wtA);
        free(rw->wfwA);
        free(rw->wflA);
        free(rw->wlA);
        free(rw->wphi_i);
        free(rw->wpsi_ij);
        
        log_trace("cleaning rSFS rweights pointer");
        free(rw);
        
        log_trace("cleaning statistics[0].rS");
        free(statistics[0].rS);
        log_trace("cleaning statistics[0].rSo");
        free(statistics[0].rSo);
        log_trace("cleaning statistics[0].rthetaS");
        free(statistics[0].rthetaS);
        log_trace("cleaning statistics[0].rthetaT");
        free(statistics[0].rthetaT);
        log_trace("cleaning statistics[0].rthetaSo");
        free(statistics[0].rthetaSo);
        log_trace("cleaning statistics[0].rthetaTo");
        free(statistics[0].rthetaTo);
        log_trace("cleaning statistics[0].rthetaFL");
        free(statistics[0].rthetaFL);
        log_trace("cleaning statistics[0].rthetaFW");
        free(statistics[0].rthetaFW);
        log_trace("cleaning statistics[0].rthetaL");
        free(statistics[0].rthetaL);
        log_trace("cleaning statistics[0].rthetaSA");
        free(statistics[0].rthetaSA);
        log_trace("cleaning statistics[0].rthetaTA");
        free(statistics[0].rthetaTA);
        log_trace("cleaning statistics[0].rDtaj");
        free(statistics[0].rDtaj);
        log_trace("cleaning statistics[0].rDfl");
        free(statistics[0].rDfl);
        log_trace("cleaning statistics[0].rFfl");
        free(statistics[0].rFfl);
        log_trace("cleaning statistics[0].rHnfw");
        free(statistics[0].rHnfw);
        log_trace("cleaning statistics[0].rEz");
        free(statistics[0].rEz);
        log_trace("cleaning statistics[0].rYach");
        free(statistics[0].rYach);
        log_trace("cleaning statistics[0].rFH");
        free(statistics[0].rFH);

    }
    
    // fix issue : free f if formatfile is not tfa
    if (args.formatfile != TFA_FORMAT)
    {
        log_trace("cleaning f");
        free(f);
    }

    log_trace("cleaning args.r2i_ploidies");
    free(args.r2i_ploidies);

    // if(file_logerr) fprintf(file_logerr,"\nProgram Ended\n");
    log_info("Program Ended");
    if (error_log_file)
        fclose(error_log_file);

    exit(0);
}

void usage(void)
{
    printf( MSTATSPOP );
    printf("Flags:\n");
    printf("      -f [input format file: ms, fasta OR tfa (gz file indexed)]\n");
    printf("      -i [path and name of the input file]\n");
    printf("      -o [output format file: 0 (extended),\n");
    printf("                              1 (single line/window),\n");
    printf("                              2 (single line SFS/window),\n");
    printf("                              3 (dadi-like format),\n");
    printf("                              4 (single line pairwise distribution)\n");
    printf("                              5 (single line freq. variant per line/window)\n");
    printf("                              6 (SNP genotype matrix)\n");
    printf("                              7 (SweepFinder-like format -only first pop-)\n");
    printf("                              8 (single line/window: Frequency of each haplotype in the populations)\n");
    printf("                              9 (single line/window: Frequency of variants per line and population)\n");
    //printf("                              92 (single line/window: -rSFS- Frequency of variants per population relative to all)\n");
    printf("                             10 (full extended)]\n");
    printf("      -N [#_pops] [#samples_pop1] ... [#samples_popN]\n");
    // printf("      -n [name of a single scaffold to analyze. For tfa can be a list separated by commas(ex. -n chr1,chr2,chr3]\n");
    printf("      -n [name of the file containing the name(s) of scaffold(s) and their length (separated by a tab), one per line (ex. fai file)]\n");
    printf("      -T [path and name of the output file]. DEFAULT stdout.\n");
    printf("   OPTIONAL GENERAL PARAMETERS:\n");
    printf("      -G [outgroup (0/1)] (last population). DEFAULT 0.\n");
    printf("      -u [include unknown positions (0/1)].  DEFAULT 0.\n");
    printf("      -R [performs analysis using only rSFS (0/1)].  DEFAULT 0.\n");
    // printf("      -D [calculate Fst if more than 6 pops (0/1)].  DEFAULT 0 (no more than 6).\n");
    printf("      -A [Alternative Spectrum File (Only for Optimal Test): alternative_spectrum for each population (except outg)\n");
    printf("          File format: (average absolute values) header plus fr(0,1) fr(0,2) ... fr(0,n-1) theta(0)/nt,\n");
    printf("          fr(1,1) fr(1,2) ... fr(1,n-1) theta(1)/nt...; -u 1 not allowed yet]\n");
    printf("      -S [Null Spectrum File (only if -A is defined): null_spectrum for each population (except outg).\n");
    printf("          (average absolute values) header plus fr(0,1) fr(0,2) ... fr(0,n-1) theta(0)/nt,\n");
    printf("          fr(1,1) fr(1,2) ... fr(1,n-1) theta(1)/nt...]. DEFAULT SNM.\n");
    // printf("      -P [Only for Calculation of R2_p: first value is the number of values to include, \n");
    // printf("                       next are the number of lines to consider. ex: -P 6 1 2 4 8 16 64]\n");
    printf("    Optional Parameters for fasta and tfa input files:\n");
    printf("      -O [#_nsam] [number order of first sample, number 0 is the first sample] [second sample] ...etc. up to nsamples.\n");
    printf("         DEFAULT current order.\n");
    printf("      -t [# permutations per window (H0: Fst=0). Only available with option -u 0]. DEFAULT 0.\n");
    printf("      -s [seed]. DEFAULT 123456.\n");
    printf("   PARAMETERS FOR TFASTA INPUT (-f tfa): 'SLIDING WINDOW ANALYSIS OF EMPIRICAL DATA'\n");
    printf("      -w [window size].\n");
    printf("        OR\n");
    printf("      -W [file with the coordinates of each window [scaffold init end] (instead options -w and -z).\n");
    printf("         DEFAULT one whole window.\n");
    printf("    Optional:\n");
    printf("      -z [slide size (must be a positive value)]. DEFAULT window size.\n");
    printf("      -Z [first window size displacement [for comparing overlapped windows])]. DEFAULT 0.\n");
    printf("      -Y [define window lengths in 'physical' positions (1) or in 'effective' positions (0)]. DEFAULT 1.\n");
    printf("      -E [input file with weights for positions:\n");
    printf("         include three columns with a header,\n");
    printf("         first the physical positions (1...end),\n");
    printf("         second the weight for positions and\n");
    printf("         third a boolean weight for the variant (eg. syn variant in nsyn counts is 0.000)].\n");
    printf("         DEFAULT all 1.000\n");
    printf("   PARAMETERS FOR MS INPUT (-f ms):'SIMULATION ANALYSIS OF A SINGLE REGION'\n");
    // printf("      -l [length]\n");//to eliminate!!
    printf("    Optional:\n");
    printf("      -r [# ms iterations]. DEFAULT 1.\n");
    printf("      -m [include mask_filename] DEFAULT -1 (all positions included).\n");
    printf("         [mask_file format: 1st row with 'length' weights, next sample rows x lengths: missing 0, sequenced 1)].\n");
    printf("         DEFAULT no mask.\n");
    printf("      -v [ratio transitions/transversions]. DEFAULT 0.5.\n");
    printf("      -F [force analysis to include outgroup (0/1) (0 in ms means ancestral)]. DEFAULT 0.\n");
    printf("      -q [frequency of reverted mutation] (only with -F 1). DEFAULT 0.\n");
    printf("   PARAMETERS FOR FASTA INPUT (-f fasta): 'WHOLE REGION ANALYSIS'\n");
    printf("    Optional:\n");
    printf("      -p [Number of lineages per sequence (1/2)]. DEFAULT 1.\n");
    printf("      -g [GFF_file]\n");
    printf("         [add also: coding,noncoding,synonymous,nonsynonymous,silent, others (or whatever annotated)]\n");
    printf("         [if 'synonymous', 'nonsynonymous', 'silent' add: Genetic_Code: Nuclear_Universal,mtDNA_Drosophila,mtDNA_Mammals,Other]\n");
    printf("         [if 'Other', introduce the code for the 64 triplets in the order UUU UUC UUA UUG ... etc.].\n");
    printf("         DEFAULT no annotation.\n");
    printf("      -c [in case use coding regions, criteria to consider transcripts (max/min/first/long)]. DEFAULT long.\n");
    printf("      -K [make a MASK file with the valid positions for this fasta. Useful for running ms simulations (1/0)]. DEFAULT 0.\n");
    printf("   HELP:\n");
    printf("      -h [help and exit]\n");
    /*printf("\t-y [in case missing: number of comparisons for approximate calculation of covariance. 0 if rough approach]\n");*/ /*to eliminate*/
    /*
     printf("\t-T [count only transitions (not for ms format)\n"); NOT DONE
     printf("\t-V [count only transversions (not for ms format)\n"); NOT DONE
     printf("\t-G [count only G/C mutations (not for ms format)\n"); NOT DONE
     printf("\t-A [count only A/T mutations (not for ms format)\n"); NOT DONE
     */
}

int read_index_file(char *chr_name_all,
                    unsigned long *nscaffolds,
                    char ***chr_name_array,
                    char ***chr_length_array)
{

    FILE *file_scaffolds;
    char *buf;
    int c;
    int k;

    *nscaffolds = 1;
    chr_name_array[0] = (char **)calloc(*nscaffolds, sizeof(char *));
    chr_name_array[0][0] = (char *)calloc(MSP_MAX_NAME, sizeof(char));
    chr_length_array[0] = (char **)calloc(*nscaffolds, sizeof(char *));
    chr_length_array[0][0] = (char *)calloc(MSP_MAX_NAME, sizeof(char));

    if (!(file_scaffolds = fopen(chr_name_all, "r")))
    {
        // printf("Error opening the scaffold names file %s\n",chr_name_all);
        log_error("Error opening the scaffold names file %s", chr_name_all);
        return (1);
    }
    if (!(buf = (char *)malloc(BUFSIZ)))
    {
        // puts("\nError: Not enough memory to read the scaffold names file.\n");
        log_error("Error: Not enough memory to read the scaffold names file.");
        return (1);
    }
    setbuf(file_scaffolds, buf);
    c = fgetc(file_scaffolds);
    while (c != EOF)
    {
        k = 0;
        chr_name_array[0][*nscaffolds - 1][k] = c;
        k++;
        while ((c = fgetc(file_scaffolds)) != 9 && c != 10 && c != 13 && c != -1 && c != 0 && k < MSP_MAX_NAME - 1)
        {
            chr_name_array[0][*nscaffolds - 1][k] = c;
            k++;
        }
        chr_name_array[0][*nscaffolds - 1][k] = '\0';
        if (c != 9 && c != 32)
        {
            // printf("Error reading the scaffold names file %s:\n scaffold (%s) without length information.\n",chr_name_all, chr_name_array[0][*nscaffolds-1]);
            log_error("Error reading the scaffold names file %s:\n scaffold (%s) without length information.",
                      chr_name_all,
                      chr_name_array[0][*nscaffolds - 1]);
            free(buf);
            return (1);
        }
        do
        {
            c = fgetc(file_scaffolds);
        } while (!(c != 9 && c != 32 && c != 10 && c != 13 && c != -1 && c != EOF));
        if (c == EOF)
        {
            // printf("Error reading the scaffold names file %s:\n scaffold (%s) without length information.\n",chr_name_all, chr_name_array[0][*nscaffolds-1]);
            log_error("Error reading the scaffold names file %s:\n scaffold (%s) without length information.",
                      chr_name_all,
                      chr_name_array[0][*nscaffolds - 1]);
            free(buf);
            return (1);
        }
        k = 0;
        chr_length_array[0][*nscaffolds - 1][k] = c;
        k++;
        while ((c = fgetc(file_scaffolds)) != 9 && c != 32 && c != 10 && c != 13 && c != -1 && c != 0 && k < MSP_MAX_NAME - 1)
        {
            chr_length_array[0][*nscaffolds - 1][k] = c;
            k++;
        }
        chr_length_array[0][*nscaffolds - 1][k] = '\0';
        /*check next line, if exist*/
        if (c == 32 || c == 9)
        {
            do
            {
                c = fgetc(file_scaffolds);
            } while (c != 10 && c != 13 && c != EOF);
        }
        while (c == 10 || c == 13)
        {
            c = fgetc(file_scaffolds);
        }
        if (c == EOF)
            break;
        /*if exist, prepare new row in arrays*/
        *nscaffolds += 1;
        chr_name_array[0] = (char **)realloc(chr_name_array[0],
                                             *nscaffolds * sizeof(char *));
        chr_name_array[0][*nscaffolds - 1] = (char *)calloc(MSP_MAX_NAME,
                                                            sizeof(char));
        chr_length_array[0] = (char **)realloc(chr_length_array[0],
                                               *nscaffolds * sizeof(char *));
        chr_length_array[0][*nscaffolds - 1] = (char *)calloc(MSP_MAX_NAME,
                                                              sizeof(char));
    }
    free(buf);
    fclose(file_scaffolds);

    return (0);
}
