/*
 *  mstatspop.c
 *  mstatspop
 *
 *  First version created by Sebastian Ramos-Onsins on 28/07/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "mstatspop.h"

/**
 * \brief main function
 * \details More details
 * @param argc
 * @param argv
 * @return
 */
/* [input_file] [#_pops] [#samples_pop1] ... [#samples_popN] [GFF_file] */
int main(int argc, const char * argv[])
{
	/*Statistics*/
	struct stats *statistics=0;
	struct stats *stats_iter=0;
	struct probs *piter=0;
	int w,x,y,z,yy	/*,ii*/;
	long int i;
	char *f;
	char ploidy[2];
	
	/* variables defining more data*/
	int npops = 0;
	
	/* Number of samples for each population, each element a population */
	int *	vint_perpop_nsam = NULL;	/* old nsam */
	/* Population index */
	int pop_index = 0;
	/* Sum of all nsam of all populations */
	int int_total_nsam = 0; /* old nsamtot*/	
	/* number of iterations (ms) or permutations (others) */
	long int niter = 0;
	
	int include_unknown;
	int output = 0; 				/* TODO: tipo de output, de 0 a 11 */										 
	int formatfile = 0;
    int coordfile;
	
	/* name and data files */
	char file_in[ MSP_MAX_FILENAME];
	char file_out[MSP_MAX_FILENAME];
	char file_mas[MSP_MAX_FILENAME];
	char file_GFF[MSP_MAX_FILENAME];
	char file_H1f[MSP_MAX_FILENAME];
	char file_H0f[MSP_MAX_FILENAME];
    char file_log[MSP_MAX_FILENAME];
	
	memset( file_in,  0, MSP_MAX_FILENAME);
	memset( file_out, 0, MSP_MAX_FILENAME);
	memset( file_mas, 0, MSP_MAX_FILENAME);
	memset( file_GFF, 0, MSP_MAX_FILENAME);
	memset( file_H1f, 0, MSP_MAX_FILENAME);
	memset( file_H0f, 0, MSP_MAX_FILENAME);
			
	FILE *file_input	= 	0;
    FILE *file_logerr   =   stdout;
	FILE *file_output	=	stdout;
	FILE *file_mask	    =	0;
	FILE *file_H1freq	=	0;
	FILE *file_H0freq	=	0;

    SGZip file_input_gz;
    SGZip file_output_gz;
    SGZip file_logerr_gz;

    /* GFF variables */
	int 	gfffiles			= 0;
	/*int 	observed_data	= 0;*/
	char 	subset_positions[ MSP_MAX_GFF_WORDLEN ];		
	char 	code_name[ MSP_GENETIC_CODETYPE_LEN ];
	char 	genetic_code[ MSP_GENCODE_COMBINATIONS +1];
	char 	criteria_transcript[ MSP_GFF_CRITERIA_MSG_LEN ];
	
	/* Alignment data variables */
    char *matrix_pol = 0; 		/*the matrix of SNPs alignment data*/
    char *matrix_pol_tcga = 0; 		/*the matrix of SNPs alignment data with four nucleotides*/
	long int *matrix_freq = 0;	/*frequency of the segregating site*/
	long int *matrix_pos = 0; 	/*position of the segregating site*/
	/*long int *matrix_GC = 0; */	/*vector with 1 if GC exists in these position and the contiguous  */
	long int *matrix_sv = 0; 	/*vector with x  for each position, where x */
										/* equals: 1 transition, 2 transversion, 0 nothing
										 * (biallelic only, no multiple hits considered) */
	long int *sites_matrix=0;
	long int length_seg=0;
	long int length = 0;		/* Sequence length */
	double length_al=0;
	long int length_al_real;
	int kind_length = 0;
	
	double svratio;
	double missratio=0.;
	
	double *	sum_sam; /*length of each sequence excluding non-tcga*/
	double **	tcga;
	long int 	nmhits=0;
	int *		matrix_mask=0; /*mask matrix for ms inputs*/
	float *		vector_mask=0; /*value of each position: eg. in case synonymous, the position can count less than 1*/
	
	/* joint frequency distribution */
	double **jfd; /*frequency of each segregating site per population*/
	int    **nfd; /*number of samples for each segregating site and population*/
	
	/* permutation test */
	char *matrix_perm = 0; /*the matrix of SNPs alignment data permutated*/ 
	int nsam2[3];		/*the sample size for each pop in permutation test. Include outg*/
	int psam2[3];		/*the position of the first sample in matrix_pol of each pop for perm. test. Include outg*/
	
	/* Optimal tests */
	int H1frq,H0frq;
	int npf=0;
	double **	freqspH1	=	0;
	double **	freqspH0	=	0;
	double *		thetaH1	=	0;
	double *		thetaH0	=	0;
	char *cad, *cad1;
	int len;
	
	long int nseed = 1234;
	long int niterdata;
	long int li,li2;
	int flaghky;

	/* read mask_file */
	char c[2],vvm[11];
	int vli;
	long int n;

	double *sum_sam_mask=0;
    double length_mask=0;
    long int length_mask_real=0;
	
	int arg;
	
	/* ms */
	float ms_svratio = 0.5;
	char *el;
	double *vector_priors=0;
	int npriors=0;
	int outgroup_presence;
	int force_outgroup;
	float freq_revert = 0.;
	double freq_missing_ms = 0.;
	int location_missing_ms =  0;
	
	/**/
	double *nsites1_pop;
	double *nsites2_pop;
	double *nsites3_pop;
	double *nsites1_pop_outg;
	double *nsites2_pop_outg;
	double *nsites3_pop_outg;
	double *anx,*bnx;
	double *anxo,*bnxo;
    double **lengthamng;
    double **lengthamng_outg;
	
	
	/*R2_p*/
	int *r2i_ploidies=0;
	
	/*Covariances in missing values*/
	int n_ccov = 0;
	
	int int_total_nsam_order=0; /*number of samples. value for ordering the samples in populations in case they are not consecutive*/
	int *sort_nsam=0; /*vector with the position in the file*/
	int sort_index;/*used in the loop*/
	
	/*tfasta windows and weights*/
	long int *wgenes;/*init and end coordinates*/
	long int nwindows;/*number of fragments*/

    long int slide=0;
    long int first_slide=0;
	long int window=0;
    char file_wps[MSP_MAX_FILENAME];
	char file_Wcoord[MSP_MAX_FILENAME];
    FILE *file_ws   	= 	0;
    SGZip file_ws_gz;

	FILE *file_wcoor    =   0;
    SGZip file_wcoor_gz;
	
    struct SGZIndex index_input;
    struct SGZIndex index_w;
    
    int Physical_length=1;
    int mask_print;
    
    int nscaffolds;
    char *chr_name = 0;
    char chr_name_all[ MSP_MAX_NAME];
    char **chr_name_array;
    int j,k;
    int first = 0;
    
    memset( chr_name_all, 0, MSP_MAX_NAME);
    memset( subset_positions, 0, 	MSP_MAX_GFF_WORDLEN );
    memset( code_name, 0, 			MSP_GENETIC_CODETYPE_LEN );
    memset( genetic_code, 0, 		MSP_GENCODE_COMBINATIONS +1);
    memset( criteria_transcript, 0, MSP_GFF_CRITERIA_MSG_LEN );
    
    memset( file_wps,    0, MSP_MAX_FILENAME);
    memset( file_Wcoord, 0, MSP_MAX_FILENAME);
    memset( c,0,2);
    
	/***** DEFAULTS *****/
    formatfile          = 3;/*tfasta format*/
    output              = 1;/*single line statistics per window*/
    length              = 0; /*undefined*/
    /*------------------*/
    npops 				= 0;/*undefined*/
    outgroup_presence   = 0;/*no outgroup*/
    include_unknown 	= 0;/*excluded missing positions*/
    int_total_nsam_order= 0;/*undefined*/
    /*file_out            = stdout;*//*already defined as stdout*/
    H1frq 				= 0;
    H0frq 				= 0;
    /*r2i_ploidies=0;*//*already defined*/
    window              = 0;
    slide               = 0;
    Physical_length     = 1;
    niter               = 0;
    nseed 				= 123456;
    niterdata 			= 1;
    strcpy( file_mas,"-1\0"); /* -1 = No file provided*/
    ms_svratio 			= 0.5;
    force_outgroup 	    = 0;
    freq_revert         = 0;
    ploidy[0]           = '1';/* 1:haploid, 2:diploid, and so on...NO, next Ver */
    ploidy[1]           = '\0';
    strcpy( criteria_transcript,"long\0");
    coordfile           = 0;
    mask_print          = 0;
    first_slide         = 0;
    /********************/

	/* STEP 1: PARSING COMMAND LINE ------------------------------------------ */
	if( argc > 1 ) 
	{
		arg = 1;
		while(arg < argc)
		{
			if( argv[arg][0] != '-' ) 
			{
                if(argv[arg][0] == '>')
                    break;
                printf(" argument should be -%s ?\n", argv[arg]);
				usage();
				exit(1);
			}
			
			switch (argv[arg][1])
			{
				case 'f' : /* f FORMAT */
					arg++;
					/* if format is not recognized, error and out */
					if( 	strcmp( argv[arg], "fasta" ) != 0 && 
							strcmp( argv[arg], "nbrf") !=0 && 
							strcmp( argv[arg], "ms")!=0 &&
							strcmp( argv[arg], "tfa")!=0 )
					{
						printf("\n Error: the argument -f has only the choices 'fasta', 'tfa' or 'ms'.");
						usage();
						exit(1);
					}
					else if( strcmp(argv[arg],"fasta")==0 || strcmp(argv[arg],"nbrf")==0) { 
						formatfile = 0;
                        niterdata = 1;
					}
					else if(strcmp(argv[arg],"ms")==0) { 
						formatfile = 1;
					}
					else if(strcmp(argv[arg],"ms_x")==0) {
						formatfile = 2;
					}
					else if(strcmp(argv[arg],"tfa")==0) {
						formatfile = 3;
					}
					
					break;
					
                case 'i': /* i Input File, el path */
                    arg++;
                    strcpy( file_in, argv[arg] );
                    break;
                
                case 'T': /* i output File, el path */
                    arg++;
                    strcpy( file_out, argv[arg] );
                    if( (file_output = fopen( file_out, "w")) == 0) { /*zipped not available yet*/
                        fprintf(stdout,"\n It is not possible to write in the output file %s\n", file_out);
                        exit(1);
                    }
                    strcpy(file_log, file_out);
                    strcat(file_log,".log");
                    if( (file_logerr = fzopen( file_log, "w", &file_logerr_gz)) == 0) {
                        fprintf(stdout,"\n It is not possible to write the log file %s.", file_log);
                        exit(1);
                    }
                    break;
				case 'o' : /* o Output type, from 0 to 11 - TODO: define*/
					arg++;
					output = (int)atoi( argv[arg] );
					if(output < 0 || output > 10) {
						printf("\n Error in -o argument: only values between 0 and 10 are allowed.");
						usage();
						exit(1);
					}
					break;
					
				case 'p' : /* p Ploidy, 1: haploid, 2: diploid */
					arg++;
					ploidy[0] = argv[arg][0];
					if(ploidy[0] != '1' && ploidy[0] != '2') {
						printf("\n Error in -p argument: only the values 1 or 2 are allowed.");
						usage();
						exit(1);
					}
					break;
					
				case 'u' : /* u Missing Values Allowed or not, 1: allowed */
					arg++;
					include_unknown = (int)atoi( argv[arg] );
					
					if(include_unknown != 0 && include_unknown != 1) {
						printf("\n Error in -u argument: only the values 0 or 1 are allowed.");
						usage();
						exit(1);
					}
					break;
					
				case 'm' : /* m Mask file, if -1, no file included */
					arg++;
					strcpy(file_mas,argv[arg]);
					break;
					
				case 's' : /* s Seed, positive value */
					arg++;
					nseed = (long int)atol(argv[arg]);
					break;
					
				case 't' : /* t number of permutations, in case of fasta format */
					arg++;
					niter = (long int)atol( argv[arg] );
					break;
 
                case 'r' : /* number of iterations, in case of ms or tfasta */
                    arg++;
                    niterdata = (long int)atol( argv[arg] );
                    break;
					
				case 'v' : /* v ratio transitions/transversions, only in ms */
					arg++;
					ms_svratio = (double)atof(argv[arg]);
					break;
					
				case 'l' : /* l length of the sequence, only for ms format */
					arg++;
					length = (long int)atol(argv[arg]);
					break;
					
				case 'k' : /* kind of l length of the sequence (0,1,2,3), only for ms format */
					arg++;
					kind_length = (int)atoi(argv[arg]);
					break;
					
				case 'G' : /* G outgroup: included=1, non included=0 */
					arg++;
					outgroup_presence = (int)atoi(argv[arg]);					
					if(outgroup_presence != 0 && outgroup_presence != 1) {
						printf("\n Error in -G argument: only the values 0 or 1 are allowed.");
						usage();
						exit(1);
					}
					break;
				case 'F' : /* F forced_outgroup=1 */
					arg++;
					force_outgroup = (int)atoi(argv[arg]);
					break;
					
				case 'q' : /* frequency of revert mutation when forcing outgroup, only in ms format */
					arg++;
					freq_revert = (double)atof(argv[arg]);
					break;
					
				case 'N' : /* N number of populations, Warning!, followed by
									N numbers indicating the sample size of each population
								*/
					arg++;
					npops = atoi(argv[arg]);
					if((vint_perpop_nsam = (int *) calloc( (unsigned long)npops, sizeof(int) )) == 0) {
						printf("Error allocating memory");
						exit(1);
					}
					int_total_nsam = 0;
					for( pop_index = 0; pop_index < npops; pop_index++ ) 
					{
						arg++;
						vint_perpop_nsam[ pop_index ] = atoi( argv[arg] );
						int_total_nsam += vint_perpop_nsam[ pop_index ];
					}					
					break;
					
				case 'O' : /* O the order of each individual in the original data, Warning!, followed by
							O numbers indicating the order (0 is the first) in case samples are not consecutive. Only for fasta data!
							*/
					arg++;
					int_total_nsam_order = atoi(argv[arg]);
					if((sort_nsam = (int *) calloc( (unsigned long)int_total_nsam_order, sizeof(int) )) == 0) {
						printf("Error allocating memory");
						exit(1);
					}
					for( sort_index = 0; sort_index < int_total_nsam_order; sort_index++ ) 
					{
						arg++;
						sort_nsam[ sort_index ] = atoi( argv[arg] );
					}					
					break;
					
				case 'g': /* g GFF file name, AND more words 
									2nd : synonymous, nonsynonymous, silent or whatever
									3rd : selected genetic code name or "Other"
									next 64th's : in case of 'Other', provide 64 IUPAC letters of each codon. 
									* Check order. 
								*/
					arg++;
					strcpy( file_GFF, argv[arg]);					
					arg++;
					strcpy( subset_positions, argv[arg] );
					
					gfffiles = 1;	
					
					/* Go if known coding option - */
					if( (strcmp(subset_positions,"synonymous")==0 ||
                         strcmp(subset_positions,"nonsynonymous")==0 ||
                         strcmp(subset_positions,"0-fold")==0 ||
                         strcmp(subset_positions,"2-fold")==0 ||
                         strcmp(subset_positions,"3-fold")==0 ||
                         strcmp(subset_positions,"4-fold")==0 ||
                         strcmp(subset_positions,"silent")==0))
					{
						arg++;
						strcpy( code_name,argv[arg] ); 

						if( strcmp(code_name, "Other")==0 ) 
						{
							for(x=0;x<64;x++) {
								arg++;
								if(argv[arg][0] == '-') {
									printf("\n Error in -g argument: In case use \"Other\", include the genetic code of the 64 aa values.");
									usage();
									exit(1);
								}
								genetic_code[x] = atoi(argv[arg]);
							}
						}
					}
					break;
					
				case 'A': /* a Name of file of the Alternative frequency spectrum */
							 /* Only with optimal tests (in GSL libs) */
					arg++;
					strcpy( file_H1f, argv[arg] );					
					H1frq = 1;
					break;
					
				case 'S': /* a Name of file of the NULL frequency spectrum */
							 /* Only with optimal tests (in GSL libs) */
					arg++;
					strcpy( file_H0f, argv[arg] );					
					H0frq = 1;
					break;
					
				case 'c' : /* c Criteria used for analyzing the transcripts */
								/* Basically, max or min */
					arg++;
					strcpy(criteria_transcript,argv[arg]);
					if( 	strcmp( criteria_transcript, "max")!=0  && 
							strcmp( criteria_transcript, "min")!=0  && 
							strcmp( criteria_transcript, "first")!=0   && 
							strcmp( criteria_transcript, "long")!=0  ) 
					{
						printf("\n Error: the argument -c has only the choices 'max', 'min', 'first' or 'long'.");
						usage();
						exit(1);
					}
					break;
					
				case 'x' : /* proportion of missing values , only in ms format with the option -u 1 */
					arg++;
					freq_missing_ms = (double)atof(argv[arg]);
					break;
					
				case 'y' : /* number of comparisons for the calculation of covariances in neutrality tests
							using missing values (0 if not desired) */
					arg++;
					n_ccov = (int)atol( argv[arg] );
					break;
					
				case 'M' : /* in case ms_e, column location of missing values at prior row (default 3), only in ms_e format with the option -u 1 */
					arg++;
					location_missing_ms = (double)atoi(argv[arg]);
					break;
					
				case 'P' : /*Calculation of R2_p: first value is the number of values to include, next are the ploidies to consider. ex: -P 6 1 2 4 8 16 64*/
					arg++;
					if( (r2i_ploidies = (int *)malloc((unsigned long)(atoi(argv[arg])+1)*sizeof(int))) == NULL ) {
						printf("\nError: memory not reallocated. mstatspop.c.00 \n");
						exit(1);
					}
					x=0;
					r2i_ploidies[x] = atoi(argv[arg]);
					while(x < r2i_ploidies[0]) {
						x++;
						arg++;
						r2i_ploidies[x] = atoi(argv[arg]);
					}
					break;
				case 'Y' : /* physical length or effective length (only valid positions) */
					arg++;
					Physical_length = (int)atoi(argv[arg]);
					if(Physical_length != 0 && Physical_length != 1) {
						printf("\n Error in -l argument: only the values 0 or 1 are allowed.");
						usage();
						exit(1);
					}
					break;
				case 'w' : /* window size */
					arg++;
					window = (long int)atol(argv[arg]);
					break;
				case 'z' : /* slide size */
					arg++;
					slide = (long int)atol(argv[arg]);
					break;					
				case 'W' : /* file with the coordinates of each window [init end](overwrite options -w and -s)*/
					arg++;
					strcpy(file_Wcoord, argv[arg] );
                    coordfile = 1;
					break;
                case 'E' : /*file with the weight for each position */
                    arg++;
                    strcpy(file_wps, argv[arg]);
                    break;
                case 'K' : /*output a mask file or not */
                    arg++;
                    mask_print = (int)atoi( argv[arg] );
                    break;
                case 'Z' : /* first slide size */
                    arg++;
                    first_slide = (long int)atol(argv[arg]);
                    break;					
                case 'n' : /* name of the scaffold to analyze*/
                    arg++;
                    strcpy( chr_name_all, argv[arg] );
                    break;
				case 'h' : /* h HEEEEEEEEEEELPPPPPPPPPPPPPP!!! */
					usage();
					exit(0);
					break;
			}
			arg++;
		}
		
		/*few filters*/
		if(!(formatfile == 1 || formatfile == 2) && force_outgroup == 1) {
			fzprintf(file_logerr,&file_logerr_gz,"\nError. The option -F 1 is only compatible with -f 'ms'.\n");
			exit(1);
		}
        if(npops == 0) {
            fzprintf(file_logerr,&file_logerr_gz,"\nError. The option -N must be included.\n");
            exit(1);
        }
        if(window <=0 && formatfile == 3 && coordfile == 0) {
            fzprintf(file_logerr,&file_logerr_gz,"\nError. The option -w or -W must be included with option -f tfa.\n");
            exit(1);
        }
        if(slide == 0 && window > 0 && coordfile == 0)
            slide = window;
		if(slide <= 0 && formatfile == 3 && coordfile == 0) {
			fzprintf(file_logerr,&file_logerr_gz,"\nError. The value at option -z (slide) must be larger than 0\n");
			exit(1);
		}
        if(first_slide < 0 && formatfile == 3) {
            fzprintf(file_logerr,&file_logerr_gz,"\nError. The value at option -Z (first slide) must be larger than or equal to 0.\n");
            exit(1);
        }
        if(length == 0 && formatfile == 1) {
            fzprintf(file_logerr,&file_logerr_gz,"\nError. length (-l option) must be defined with ms input file.\n");
            exit(1);
        }
        if(formatfile == 0 && niterdata > 1) {
            fzprintf(file_logerr,&file_logerr_gz,"\nError. The option ");
            exit(1);
        }
        if(include_unknown == 1)
            niter = 0;
		
        if(strcmp(chr_name_all,"") == 0) {
               fzprintf(file_logerr,&file_logerr_gz,"\nError: the name of the scaffold (option -n) must be defined\n");
               exit(1);
        }
        /*
        if(file_Wcoord[0]!=0 && (slide > 0 && window>0)) {
            fzprintf(file_logerr,&file_logerr_gz,"\n the option -W (coordinates file) is incompatible with definitions of -w and -z ");
            exit(1);
        }
        */
        /* STEP 2: READING FILE INFO --------------------------------------------- */
		
		/* Opening files */
		if( file_in[0] == '\0' ) {
			file_input = stdin;
			file_input_gz.file_compressed = 0;
            if(formatfile == 3) {
                /*tfa must be a file because it needs the index file*/
                fzprintf(file_logerr,&file_logerr_gz,"\nError: tfa format file needs indexation. \nstdout is not available using this option. \n");
                exit(1);
            }
		}
		else {
			if( (file_input = fzopen( file_in, "r", &file_input_gz)) == 0) {
				fzprintf(file_logerr,&file_logerr_gz,"\n It is not possible to open the input file %s.", file_in);
				exit(1);
			}
            if(formatfile == 3) {/*tfasta file*/
                load_index_from_file(file_input_gz.index_file_name, &index_input);
            }
		}
		
		if( (f = (char *)malloc((unsigned long)BUFSIZ*10)) == NULL ) {
			fzprintf(file_logerr,&file_logerr_gz,"\nError: memory not reallocated. main.4 \n");
			exit(1);
		}
		/* Definition of a File Stream Buffer, for buffered IO */
		setbuf(file_input,f);
				
        /*separate all values of the list chr_name_all in chr_name_array: */
        /* Only do the list if input and output is tfa*/
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
                
        /*if( (formatfile == 1 || formatfile == 2) && output == 0 && niter > 1) output = 1;*/
        /*if( (formatfile == 1 || formatfile == 2)) length_al = length;*/
		if( formatfile == 0) strcpy(file_mas,file_in);		
		if( formatfile == 3 && ((slide == 0 && window == 0) && file_Wcoord[0]=='\0')) 
			strcpy(file_mas,file_in);		
		
		/* TODO: por eso la suma nunca es 2, o 0 o 1... */
		if( outgroup_presence == 1) {
			 force_outgroup = 0;
		}		
		else if(outgroup_presence == 0) 
		{
			/* Al cargar los datos, -N, hemos rellenado el vector */
			/* TODO: Check for errors in realloc */
			if((vint_perpop_nsam = (int *)realloc( vint_perpop_nsam, (npops+1)*sizeof(int) )) == 0) {
				fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
				exit(1);
			}
			
			vint_perpop_nsam[npops] = 1;
			int_total_nsam += 1;
			npops++;
		}
		
		if( force_outgroup == 1) outgroup_presence = 0;
		if( include_unknown == 1 && ploidy[0] == '2' && n_ccov > 0) 
			n_ccov = 0; /*warning?*/
		
		/* -------------------------------------------------------------------- */
		/* -a alternative frequency spectrum: sum=S */
		
		if(H1frq == 1) 
		{
			npf = npops - 1; /* TODO: Mirar de cambiar por npops - 1...*/
						
			if((file_H1freq = fopen(file_H1f, "r")) == 0) {
				fzprintf(file_logerr,&file_logerr_gz,"\n Error opening Alternative frequency spectrum file %s",file_H1f);
				exit(1);
			}
			
			/* Solo lee la primera linea, y un máximo de caracteres */
			if( ! feof( file_H1freq ) ) 
			{
				/* TODO: check for errors */
				/* TODO: Cambiar para leer longitudes indeterminadas, no MAX_LEN
				 */
				/* antes 1026 y 1024 */
				if((cad = (char *) calloc( MSP_MAX_FILELINE_LEN, sizeof(char))) == 0) {
					fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
					exit(1);
				}
				
				fgets( cad, MSP_MAX_FILELINE_LEN, file_H1freq); 
			}
			else {
				fzprintf(file_logerr,&file_logerr_gz,"\n Error reading Alternative frequency spectrum file %s", file_H1f);
				exit(1);
			}

			if((freqspH1 	= (double **)	calloc( (unsigned long)npf, sizeof(double *))) == 0) {
				fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
				exit(1);
			}
			
			if((freqspH0 	= (double **)	calloc( (unsigned long)npf, sizeof(double *))) == 0) {
				fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
				exit(1);
			}
			
			if((thetaH1 		= (double *)	calloc( (unsigned long)npf, sizeof(double))) == 0) {
				fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
				exit(1);
			}
			
			if((thetaH0 		= (double *)	calloc( (unsigned long)npf, sizeof(double))) == 0) {
				fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
				exit(1);
			}
			
			
			for(x=0;x<npf;x++) 
			{
				if((freqspH1[x] = (double *) calloc((unsigned long)vint_perpop_nsam[x], sizeof(double))) == 0) {
					fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
					exit(1);
				}
				
				if((freqspH0[x] = (double *) calloc((unsigned long)vint_perpop_nsam[x], sizeof(double))) == 0) {
					fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
					exit(1);
				}
					
				
				if( !feof(file_H1freq)) {
					fgets( cad, MSP_MAX_FILELINE_LEN, file_H1freq);
					cad1 = cad;	/* TODO: REVISAR, copiando punteros! */
				}
				else {
					fzprintf(file_logerr,&file_logerr_gz,"\n  Error reading Alternative frequency spectrum file %s, line %d",file_H1f,x+2);
					exit(1);
				}
				
				for( y=1; y < vint_perpop_nsam[x]; y++) 
				{
					if( (len = (int)strlen(cad1)) > 0) {
						freqspH1[x][y] = (double) atof( cad1 );
						cad1 = strchr(cad1,'\t');
						while( cad1 && ( *cad1 == ' '|| *cad1 == '\t') ) {
							cad1 += 1;
						}
					}
				}
				thetaH1[x] = (double)atof(cad1);
			}
			free(cad);
			fclose(file_H1freq);
		}
		else {
			H1frq = 0;
		}
		
		/* -------------------------------------------------------------------- */
		/* -n null frequency spectrum: sum=S */ 
		/* TODO: Es igual que la anterior !!!!!!!!
		 * REFACTOR IT!!!!
		 * */
		if(H0frq == 1 && H1frq == 1) 
		{
			if((file_H0freq = fopen(file_H0f,"r")) == 0) {
				fzprintf(file_logerr,&file_logerr_gz,"\n Error opening NULL frequency spectrum file %s",file_H0f);
				exit(1);
			}
			if(!feof(file_H0freq)) {
				if((cad = (char *)calloc(1026,sizeof(char))) == 0) {
					fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
					exit(1);
				}
				
				fgets(cad,1024,file_H1freq);
			}
			else {
				fzprintf(file_logerr,&file_logerr_gz,"\n Error reading NULL frequency spectrum file %s",file_H0f);
				exit(1);
			}
			
			npf = npops-1;
			
			for(x=0;x<npf;x++) {
				if(!feof(file_H0freq)) {
					fgets(cad,1024,file_H0freq);
					cad1 = cad;
				}
				else {
					fzprintf(file_logerr,&file_logerr_gz,"\n  Error reading NULL frequency spectrum file %s, line %d",file_H0f,x+2);
					exit(1);
				}
				for(y=1;y<vint_perpop_nsam[x];y++) {
					if((len = (int)strlen(cad1)) > 0) {
						freqspH0[x][y] = (double)atof(cad1);
						cad1 = strchr(cad1,'\t');
						while(cad1 && (*cad1 == ' '|| *cad1 == '\t')) 
							cad1 += 1;
					}
				}
				thetaH0[x] = (double)atof(cad1);
			}
			free(cad);
			fclose(file_H0freq);
		}
		else {
			H0frq = 0;
			if(H1frq) {
				npf = npops-1;
				for(x=0;x<npf;x++) {
					for(y=1;y<vint_perpop_nsam[x];y++) {
						freqspH0[x][y] = 1./(double)y;
					}
				}
			}
		}
	
		/* -------------------------------------------------------------------- */
		/*GFF files*/
	
		/*observed_data = 0;*/
	
		if( gfffiles == 1 &&
			(	strcmp(subset_positions,"synonymous")==0 || 
                strcmp(subset_positions,"nonsynonymous")==0 ||
                strcmp(subset_positions,"0-fold")==0 ||
                strcmp(subset_positions,"2-fold")==0 ||
                strcmp(subset_positions,"3-fold")==0 ||
                strcmp(subset_positions,"4-fold")==0 ||
				strcmp(subset_positions,"silent")==0)
			) 
		{	
			/* TODO: Reordenar, esto ya esta definido en el gff_data.c !
			 */
			if( strcmp(code_name, "Nuclear_Universal") == 0 ) 
			{
				genetic_code[0] = 'F';
				genetic_code[1] = 'F';
				genetic_code[2] = 'L';
				genetic_code[3] = 'L';
				genetic_code[4] = 'S';
				genetic_code[5] = 'S';
				genetic_code[6] = 'S';
				genetic_code[7] = 'S';
				genetic_code[8] = 'Y';
				genetic_code[9] = 'Y';
				genetic_code[10] = '*';
				genetic_code[11] = '*';
				genetic_code[12] = 'C';
				genetic_code[13] = 'C';
				genetic_code[14] = '*';
				genetic_code[15] = 'W';
				genetic_code[16] = 'L';
				genetic_code[17] = 'L';
				genetic_code[18] = 'L';
				genetic_code[19] = 'L';
				genetic_code[20] = 'P';
				genetic_code[21] = 'P';
				genetic_code[22] = 'P';
				genetic_code[23] = 'P';
				genetic_code[24] = 'H';
				genetic_code[25] = 'H';
				genetic_code[26] = 'Q';
				genetic_code[27] = 'Q';
				genetic_code[28] = 'R';
				genetic_code[29] = 'R';
				genetic_code[30] = 'R';
				genetic_code[31] = 'R';
				genetic_code[32] = 'I';
				genetic_code[33] = 'I';
				genetic_code[34] = 'I';
				genetic_code[35] = 'M';
				genetic_code[36] = 'T';
				genetic_code[37] = 'T';
				genetic_code[38] = 'T';
				genetic_code[39] = 'T';
				genetic_code[40] = 'N';
				genetic_code[41] = 'N';
				genetic_code[42] = 'K';
				genetic_code[43] = 'K';
				genetic_code[44] = 'S';
				genetic_code[45] = 'S';
				genetic_code[46] = 'R';
				genetic_code[47] = 'R';
				genetic_code[48] = 'V';
				genetic_code[49] = 'V';
				genetic_code[50] = 'V';
				genetic_code[51] = 'V';
				genetic_code[52] = 'A';
				genetic_code[53] = 'A';
				genetic_code[54] = 'A';
				genetic_code[55] = 'A';
				genetic_code[56] = 'D';
				genetic_code[57] = 'D';
				genetic_code[58] = 'E';
				genetic_code[59] = 'E';
				genetic_code[60] = 'G';
				genetic_code[61] = 'G';
				genetic_code[62] = 'G';
				genetic_code[63] = 'G';
			}
			else if(strcmp(code_name,"mtDNA_Drosophila")==0) 
			{
				genetic_code[0] = 'F';
				genetic_code[1] = 'F';
				genetic_code[2] = 'L';
				genetic_code[3] = 'L';
				genetic_code[4] = 'S';
				genetic_code[5] = 'S';
				genetic_code[6] = 'S';
				genetic_code[7] = 'S';
				genetic_code[8] = 'Y';
				genetic_code[9] = 'Y';
				genetic_code[10] = '*';
				genetic_code[11] = '*';
				genetic_code[12] = 'C';
				genetic_code[13] = 'C';
				genetic_code[14] = 'W';
				genetic_code[15] = 'W';
				genetic_code[16] = 'L';
				genetic_code[17] = 'L';
				genetic_code[18] = 'L';
				genetic_code[19] = 'L';
				genetic_code[20] = 'P';
				genetic_code[21] = 'P';
				genetic_code[22] = 'P';
				genetic_code[23] = 'P';
				genetic_code[24] = 'H';
				genetic_code[25] = 'H';
				genetic_code[26] = 'Q';
				genetic_code[27] = 'Q';
				genetic_code[28] = 'R';
				genetic_code[29] = 'R';
				genetic_code[30] = 'R';
				genetic_code[31] = 'R';
				genetic_code[32] = 'I';
				genetic_code[33] = 'I';
				genetic_code[34] = 'M';
				genetic_code[35] = 'M';
				genetic_code[36] = 'T';
				genetic_code[37] = 'T';
				genetic_code[38] = 'T';
				genetic_code[39] = 'T';
				genetic_code[40] = 'N';
				genetic_code[41] = 'N';
				genetic_code[42] = 'K';
				genetic_code[43] = 'K';
				genetic_code[44] = 'S';
				genetic_code[45] = 'S';
				genetic_code[46] = 'S';
				genetic_code[47] = 'S';
				genetic_code[48] = 'V';
				genetic_code[49] = 'V';
				genetic_code[50] = 'V';
				genetic_code[51] = 'V';
				genetic_code[52] = 'A';
				genetic_code[53] = 'A';
				genetic_code[54] = 'A';
				genetic_code[55] = 'A';
				genetic_code[56] = 'D';
				genetic_code[57] = 'D';
				genetic_code[58] = 'E';
				genetic_code[59] = 'E';
				genetic_code[60] = 'G';
				genetic_code[61] = 'G';
				genetic_code[62] = 'G';
				genetic_code[63] = 'G';
			}
			else if( strcmp(code_name,"mtDNA_Mammals") == 0 ) 
			{
				genetic_code[0] = 'F';
				genetic_code[1] = 'F';
				genetic_code[2] = 'L';
				genetic_code[3] = 'L';
				genetic_code[4] = 'S';
				genetic_code[5] = 'S';
				genetic_code[6] = 'S';
				genetic_code[7] = 'S';
				genetic_code[8] = 'Y';
				genetic_code[9] = 'Y';
				genetic_code[10] = '*';
				genetic_code[11] = '*';
				genetic_code[12] = 'C';
				genetic_code[13] = 'C';
				genetic_code[14] = 'W';
				genetic_code[15] = 'W';
				genetic_code[16] = 'L';
				genetic_code[17] = 'L';
				genetic_code[18] = 'L';
				genetic_code[19] = 'L';
				genetic_code[20] = 'P';
				genetic_code[21] = 'P';
				genetic_code[22] = 'P';
				genetic_code[23] = 'P';
				genetic_code[24] = 'H';
				genetic_code[25] = 'H';
				genetic_code[26] = 'Q';
				genetic_code[27] = 'Q';
				genetic_code[28] = 'R';
				genetic_code[29] = 'R';
				genetic_code[30] = 'R';
				genetic_code[31] = 'R';
				genetic_code[32] = 'I';
				genetic_code[33] = 'I';
				genetic_code[34] = 'M';
				genetic_code[35] = 'M';
				genetic_code[36] = 'T';
				genetic_code[37] = 'T';
				genetic_code[38] = 'T';
				genetic_code[39] = 'T';
				genetic_code[40] = 'N';
				genetic_code[41] = 'N';
				genetic_code[42] = 'K';
				genetic_code[43] = 'K';
				genetic_code[44] = 'S';
				genetic_code[45] = 'S';
				genetic_code[46] = '*';
				genetic_code[47] = '*';
				genetic_code[48] = 'V';
				genetic_code[49] = 'V';
				genetic_code[50] = 'V';
				genetic_code[51] = 'V';
				genetic_code[52] = 'A';
				genetic_code[53] = 'A';
				genetic_code[54] = 'A';
				genetic_code[55] = 'A';
				genetic_code[56] = 'D';
				genetic_code[57] = 'D';
				genetic_code[58] = 'E';
				genetic_code[59] = 'E';
				genetic_code[60] = 'G';
				genetic_code[61] = 'G';
				genetic_code[62] = 'G';
				genetic_code[63] = 'G';
			}
			else if( strcmp(code_name,"Other") == 0 ) {
				; /* TODO: Con verbose, dejar rastro para debuggar */
				/* Con control de errores, comprobar que es correcto (IUPAC,etc) */
			}
			else {
					fzprintf(file_logerr,&file_logerr_gz," %s: Unknown code, sorry", code_name);
					exit(1);
			}	
		}
		
		/*ordering data: in case O is not a flag included*/
		if(int_total_nsam_order > 0 && int_total_nsam_order+!outgroup_presence != int_total_nsam) {
			fzprintf(file_logerr,&file_logerr_gz,"Error: the number of samples defined in -N and -O are different");
			exit(1);
		}
		if(int_total_nsam_order == 0) {
			int_total_nsam_order = int_total_nsam-!outgroup_presence;
			if((sort_nsam = (int *) calloc( (unsigned long)int_total_nsam, sizeof(int) )) == 0) {
				fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
				exit(1);
			}
			for( sort_index = 0; sort_index < int_total_nsam; sort_index++ ) {
				sort_nsam[ sort_index ] = sort_index;
			}
		}
		
		/* -------------------------------------------------------------------- */
        /* -------------------------------------------------------------------- */
		/* -------------------------------------------------------------------- */
		
		if(output == 0 || output == 10) {
			fprintf(file_output,MSTATSPOP);
			fprintf(file_output,"\nmstatspop ");
			for(x=1;x<arg;x++) {
				fprintf(file_output,"%s ",argv[x]);
			}
			fprintf(file_output,"\n\n");
			fprintf(file_output,"\n****************************************************************************\n");	
			fprintf(file_output,  "*  NUCLEOTIDE VARIABILITY, NEUTRALITY TEST AND POPULATION DIFFERENTIATION  *\n");
			fprintf(file_output,  "****************************************************************************\n");
		}
	}
	else {
		usage();
		exit(0);
	}
		
	/* STEP 3: And finally, some calculations ---------------------------------*/
	init_seed1(nseed);

	/*alloc memory for lengths of populations*/
	if((nsites1_pop 		= (double *)	calloc( (unsigned long)npops, sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
		exit(1);
	}
	if((nsites2_pop 	= (double *)	calloc( (unsigned long)npops, sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
		exit(1);
	}
	if((nsites3_pop 	= (double *)	calloc( (unsigned long)npops, sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
		exit(1);
	}
	if((nsites1_pop_outg = (double *)	calloc( (unsigned long)npops, sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
		exit(1);
	}
	if((nsites2_pop_outg = (double *)	calloc( (unsigned long)npops, sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
		exit(1);
	}
	if((nsites3_pop_outg = (double *)	calloc( (unsigned long)npops, sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
		exit(1);
	}
	if((anx = (double *)	calloc( (unsigned long)npops, sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
		exit(1);
	}
	if((bnx = (double *)	calloc( (unsigned long)npops, sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
		exit(1);
	}
	if((anxo = (double *)	calloc( (unsigned long)npops, sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
		exit(1);
	}
	if((bnxo = (double *)	calloc( (unsigned long)npops, sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
		exit(1);
	}
    if((lengthamng = (double **)	calloc( (unsigned long)npops, sizeof(double *))) == 0) {
        fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
        exit(1);
    }
    if((lengthamng_outg = (double **)	calloc( (unsigned long)npops, sizeof(double *))) == 0) {
        fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
        exit(1);
    }
    for(x=0;x<npops;x++) {
        if((lengthamng[x] = (double *) calloc( (unsigned long)npops, sizeof(double))) == 0) {
            fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
            exit(1);
        }
        if((lengthamng_outg[x] = (double *) calloc( (unsigned long)npops, sizeof(double))) == 0) {
            fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
            exit(1);
        }
    }
	
	/* introduce the data and mask file (last if necessary) */
	
	if(mask_print==1 && (formatfile == 0 || ( formatfile == 3 && ((slide == 0 && window == 0) && file_Wcoord[0]=='\0'))))
	{
		if( file_in[0] == '\0' ) {
            sprintf( file_mas, "file_npop%d_nsam%d",
                    npops-!outgroup_presence,
                    int_total_nsam-!outgroup_presence);
            if(gfffiles == 1) {
                strcat(file_mas,"_");
                strcat(file_mas,subset_positions);
                strcat(file_mas,"_");
                strcat(file_mas,criteria_transcript);
            }
            if(!include_unknown) strcat(file_mas,"_ExcludeMissingVariantsmhits");
            else strcat(file_mas,"_IncludeMissingVariantsmhits");
            if(outgroup_presence==0) strcat(file_mas,"_NOoutg");
            if(outgroup_presence==1) strcat(file_mas,"_outg");
            if(ploidy[0]=='1') strcat(file_mas,"_ploidy1");
            if(ploidy[0]=='2') strcat(file_mas,"_ploidy2");
            strcat(file_mas,"_MASK.txt");
			
            if((file_mask = fopen(file_mas,"w")) == 0) {
				fzprintf(file_logerr,&file_logerr_gz,"Error in mask file %s.",file_out);
				exit(1);
			}			
			/*file_mask = stderr;*/ /* TODO: My GOD!!!! */
			/* TODO: formatfile=0 from pipeline => solo fasta o nbrf=> generara un 
			 * fichero de mascara, que no sirve pa na => lo tiramos a la basura =>
			 * el stderr no es la basura!!!! Es para ERRORES!!! FIX IT.*/
		}
		else 
		{
            el = strrchr(file_mas,'.');
            *el = '\0';
            sprintf( file_mas, "%s_npop%d_nsam%d", file_mas,
                    npops-!outgroup_presence,
                    int_total_nsam-!outgroup_presence);
            if(gfffiles == 1) {
                strcat(file_mas,"_");
                strcat(file_mas,subset_positions);
                strcat(file_mas,"_");
                strcat(file_mas,criteria_transcript);
            }
            if(!include_unknown) strcat(file_mas,"_ExcludeMissingVariantsmhits");
            else strcat(file_mas,"_IncludeMissingVariantsmhits");
            if(outgroup_presence==0) strcat(file_mas,"_NOoutg");
            if(outgroup_presence==1) strcat(file_mas,"_outg");
            if(ploidy[0]=='1') strcat(file_mas,"_ploidy1");
            if(ploidy[0]=='2') strcat(file_mas,"_ploidy2");
            strcat(file_mas,"_MASK.txt");

            if((file_mask = fopen(file_mas,"w")) == 0) {
				fzprintf(file_logerr,&file_logerr_gz,"Error in mask file %s.",file_out);
				exit(1);
			}
		}
	}

	/* diploid or haploid */
	if(ploidy[0]=='2' && formatfile == 0) 
	{
		int_total_nsam *= 2;
		
		for(x=0;x<npops;x++) {
			vint_perpop_nsam[x] *= 2;
		}
		if(outgroup_presence == 0) {
			int_total_nsam -= 1;
			vint_perpop_nsam[npops-1] -= 1;
		}
	}

	/*calloc struct stats*/
	/* TODO: Fix error checking */
	statistics = 0;
	if((statistics = (struct stats *)calloc(1,sizeof(struct stats))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
		exit(1);
	}
	/* length of each sequence, excluding non-tcga */
	if((sum_sam = 	(double *) calloc(int_total_nsam+(!outgroup_presence), sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
		exit(1);
	}
	/*tcga content for each sample*/
	if((tcga = (double **) calloc(int_total_nsam+(!outgroup_presence), sizeof(double *))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
		exit(1);
	}
	for( x=0; x<int_total_nsam+(!outgroup_presence); x++) 
	{
		if((tcga[x] = (double *) calloc(4,sizeof(double ))) == 0) { /*tcga content*/
			fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
			exit(1);
		}
	}

	/******** READ FASTA DATA AND WRITE FILE_MASK *************/
	if( formatfile == 0 ) 
	{
		/*data from fasta*/
        chr_name = chr_name_array[0];
        first = 0;
        
		if (get_obsdata(file_output,&file_output_gz,
				        file_input, &file_input_gz,
                        file_logerr,&file_logerr_gz,
						file_mask,
						file_GFF, gfffiles,
						subset_positions, genetic_code, &matrix_pol, &matrix_freq,
						&matrix_pos, &length, &length_seg, &length_al,
						&length_al_real, argc, ploidy, vint_perpop_nsam,
						npops, &svratio, &missratio, include_unknown,
						sum_sam,tcga, &matrix_sv, &nmhits, output,
						outgroup_presence, criteria_transcript,nsites1_pop,nsites1_pop_outg,
						nsites2_pop,nsites2_pop_outg,nsites3_pop,nsites3_pop_outg,
						anx,bnx,anxo,bnxo,lengthamng,lengthamng_outg,
						sort_nsam,&matrix_pol_tcga,chr_name,first) == 0)
		{
			fzprintf(file_logerr,&file_logerr_gz,"Error processing input data.\n");
			exit(1);
		}
		
		fzclose(file_input, &file_input_gz);
		if(file_mask) fclose(file_mask);

		niterdata = 1;
	}
	
    /*DEFINE PARAMETERS FOR MS MASK FILE AND READ DATA IF DEFINED*/
	if( formatfile == 1 || formatfile == 2 )  
	{	/* MASK FILE MS FORMAT => Ponerlo en una funcion */
		niter = 0; /*permutation tests invalidated*/
		svratio = ms_svratio;/*-10000*/
        length_al_real = length;
        length_al = length_al_real;
	
		if((vector_mask  = (float *)calloc((unsigned int)length,sizeof(float))) == 0) {
			fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
			exit(1);
		}
		
		if((matrix_mask = (int *)calloc((unsigned int)(int_total_nsam+1)*length,sizeof(int))) == 0) {
			fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
			exit(1);
		}
		
		/*if(include_unknown) {*/
		/*if(file_mas[0] == '-'  && file_mas[1] == '1') {*/ /*all positions are accepted*//*
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
		if(file_mas[0] != '-') /*file_mask defined*/
        {
			if((file_mask = fopen (file_mas,"r")) == 0) {
				fzprintf(file_logerr,&file_logerr_gz,"\n  It is not possible to open the input mask file %s.",file_mas);
				exit(1);
			}
			li=0;n=0;vli=0;
			*c = fgetc(file_mask);
			while(*c != 0 && *c != -1 && n < int_total_nsam-!outgroup_presence+1 && li < length) {
				while(*c != 10 && *c != 13 && *c != 0 && *c != -1) {
					if(n==0) {
						vvm[vli] = *c;
						while((*c = fgetc(file_mask)) != 32 && *c != 10 && *c != 13 && vli < 10) {
							vli++;
							vvm[vli] = *c;
						}
						vvm[vli+1] = (char)0; vector_mask[(unsigned int)li] = (float)atof(vvm); /*value of the position, usually 1 except for noncounting values or Syn/Nsyn (between 0 and 1)*/
 						for(vli=0;vli<10;vli++) vvm[vli] = (char)0;
						vli = 0;
					}
					else {
						matrix_mask[(unsigned int)(n-1)*(unsigned int)length+(unsigned int)li] = (int)atoi(c) - 1;/*in case not defined file, all values are zero, in case file exist, normal is zero, missing is -1*/
						if((int)atoi(c) == 0 &&
                           include_unknown == 0)
                            vector_mask[(unsigned int)li] = 0.; /*not counting missing values when not included*/
					}
					while((*c = fgetc(file_mask)) == 32);
					li++;
				}
				if(*c == 10 || *c == 13 || li > length) {
					if(li > length) {
						fzprintf(file_logerr,&file_logerr_gz,"\n  Error: Length of rows in mask file %s are longer than defined (row %ld is %ld > %ld). ",file_mas, n, li, length);
						exit(1);
					}
					n++;
					li = 0;
					*c = fgetc(file_mask);
				}
			}
			fclose(file_mask);
			
			if(outgroup_presence == 0) {/*the "cryptic" outgroup is all 1*/
				for(li=0;li<length;li++) {
					matrix_mask[(unsigned int)(int_total_nsam)*(unsigned int)length+(unsigned int)li] = 0;
				}
			}
		/*}*/
		/*}*/
			/*if(include_unknown) {*/
				if((sum_sam_mask = (double *)calloc(int_total_nsam,sizeof(double))) == 0) {
					fzprintf(file_logerr,&file_logerr_gz,"Error allocating memory");
					exit(1);
				}
            /*}*/
            length_mask = 0.;
            length_mask_real = 0;
            li2 = 0;
            for(li=0;li<length;li++) {
                if(vector_mask[li] > 0.) {
                    if(npops > 1) {
                        x = 0;
                        for(n=0;n<int_total_nsam-vint_perpop_nsam[npops-1];n++) {
                            li2 += 1;
                            if(matrix_mask[(unsigned int)n*(unsigned int)length+(unsigned int)li] == 0) {
                                sum_sam_mask[n] += 1.; /*1 sum, 0 no sum*/
                                x += 1;
                            }
                            else missratio += 1.;
                        }
                    }
                    else x = 1;
                    y = 0;
                    for(n=int_total_nsam-vint_perpop_nsam[npops-1];n<int_total_nsam;n++) {
                        if(npops == 1) li2 += 1;
                        if(matrix_mask[(unsigned int)n*(unsigned int)length+(unsigned int)li] == 0) {
                            sum_sam_mask[n] += 1.; /*1 sum, 0 no sum*/
                            y += 1;
                        }
                        else {
                            if(npops == 1) missratio += 1.;
                        }
                    }
                    if(x > 0 && y > 0) {
                        length_mask_real += (long int)1; /*length where the outgroup and one of the rest lines exists*/
                        length_mask += vector_mask[(unsigned int)li]; /*length where the outgroup and one of the rest lines exists*/
                    }
                    else {
                         if(x==0 || (y==0 && x>0)) {
                            missratio -= ((int_total_nsam - vint_perpop_nsam[npops-1]) - x);
                            li2 -= (int_total_nsam - vint_perpop_nsam[npops-1]);
                        }
                        if(y==0 && npops == 1) {
                            missratio -= (int_total_nsam);
                            li2 -= (int_total_nsam);
                        }
                    }
                }
            }
            if(li2) missratio = (double)missratio/(double)li2;
            else missratio = (double)1;
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

	if((statistics[0].Sanc = (long int *)calloc(4*npops,sizeof(long int))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	if((statistics[0].piw  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].pia  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	if((statistics[0].piT  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].piant  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	if((statistics[0].piTnt  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].fst  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
    if((statistics[0].piwHKY  = (double *)calloc(1*npops,sizeof(double))) == 0) {
        fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
        exit(1);
    }
	if((statistics[0].piaHKY  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	if((statistics[0].piTHKY  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].fstHKY  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].fst1all  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].hapw = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].hapa = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].hapT = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].fsth = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].fsth1all = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].Gst = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	

	if((statistics[0].S  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].So  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].thetaS  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].thetaSo  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].thetaT  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].thetaTo  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].thetaTHKY  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].thetaFL  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].thetaFW  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].thetaL  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].thetaSA  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].thetaTA  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].K     = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].KHKY    = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].Dtaj  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].Dfl  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].Ffl  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].Hnfw  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].Ez  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
    if((statistics[0].Yach  = (double *)calloc(1*npops,sizeof(double))) == 0) {
        fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
        exit(1);
    }
    if((statistics[0].FH  = (double *)calloc(1*npops,sizeof(double))) == 0) {
        fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
        exit(1);
    }
	
	if((statistics[0].R2  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].Fs  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].nhpop = (int *)calloc(1*npops,sizeof(int))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].length = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	if((statistics[0].length2 = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
    if((statistics[0].lengthamng  = (double **)calloc(npops,sizeof(double *))) == 0) {
        fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
        exit(1);
    }
    if((statistics[0].lengthamng_outg  = (double **)calloc(npops,sizeof(double *))) == 0) {
        fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
        exit(1);
    }
	for(x=0;x<npops;x++) {
        if((statistics[0].lengthamng[x]  = (double *)calloc(npops,sizeof(double))) == 0) {
            fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
            exit(1);
        }
        if((statistics[0].lengthamng_outg[x]  = (double *)calloc(npops,sizeof(double))) == 0) {
            fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
            exit(1);
        }
	}
	if((statistics[0].total_tcga = (double *)calloc(4,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].tcga = (double **)calloc(1*npops,sizeof(double *))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].sv = (double ***)calloc(1*npops,sizeof(double **))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	if((statistics[0].svT = (double ***)calloc(1*npops,sizeof(double **))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].freq   = (long int **)calloc(1*npops,sizeof(long int *))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].freqh  = (long int **)calloc(1*npops,sizeof(long int *))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].ToH0_ii  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].To_ii  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].To_00  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].To_i0  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].ToH0_00  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].To_Qc_ii  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].To_Qw_ii  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].To_Lc_ii  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	
	if((statistics[0].Rm  = (long int *)calloc(1*npops,sizeof(long int))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].ZnA  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	

	if((statistics[0].mdsd  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].mdg1  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	if((statistics[0].mdg2  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	if((statistics[0].anx  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	if((statistics[0].bnx  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	if((statistics[0].anxo  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	if((statistics[0].bnxo  = (double *)calloc(1*npops,sizeof(double))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	if((statistics[0].mdw   = (double **)calloc(1*npops,sizeof(double *))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	

	if((statistics[0].linefreq   = (double **)calloc(int_total_nsam,sizeof(double *))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	
	
	if(H1frq) {
		statistics[0].H1freq = freqspH1; /*pointer to Alternative Frequency Spectrum*/
		statistics[0].H0freq = freqspH0; /*pointer to Null Frequency Spectrum*/
		statistics[0].thetaH1 = thetaH1; /*pointer to Alternative theta*/
		statistics[0].thetaH0 = thetaH0; /*pointer to Null theta*/
	}
	
	for(x=0;x<npops;x++) {
		if((statistics[0].tcga [x]   = (double *)calloc(4,sizeof(double))) == 0) {
			fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
			exit(1);
		}
		
		if((statistics[0].freq [x]   = (long int *)calloc(vint_perpop_nsam[x],sizeof(long int))) == 0) {
			fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
			exit(1);
		}
		
		if((statistics[0].freqh[x]   = (long int *)calloc(int_total_nsam,sizeof(long int))) == 0) {
			fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
			exit(1);
		}
		
		if((statistics[0].sv[x]      = (double **)calloc(1*npops,sizeof(double *))) == 0) {
			fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
			exit(1);
		}
		if((statistics[0].svT[x]      = (double **)calloc(1*npops,sizeof(double *))) == 0) {
			fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
			exit(1);
		}
		
		if((statistics[0].mdw[x]     = (double *)calloc((vint_perpop_nsam[x]*(vint_perpop_nsam[x]-1))/2,sizeof(double))) == 0) {
			fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
			exit(1);
		}
		
		for(y=0;y<npops;y++) {
			if((statistics[0].sv[x][y] = (double *)calloc(2,sizeof(double))) == 0) {
				fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
				exit(1);
			}
			if((statistics[0].svT[x][y] = (double *)calloc(2,sizeof(double))) == 0) {
				fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
				exit(1);
			}
		}
	}
	for(x=0;x<int_total_nsam;x++) {
		if((statistics[0].linefreq[x]= (double *)calloc(int_total_nsam+1,sizeof(double))) == 0) {
			fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
			exit(1);
		}
		
	}
	/*in case r2i_ploidies is undefined*/
	if(r2i_ploidies == 0) {
		if( (r2i_ploidies = (int *)calloc(2,sizeof(int)))== 0) {
			fzprintf(file_logerr,&file_logerr_gz,"\nError: memory not reallocated. mstatspop.c.00 \n");
			exit(1);
		}
		r2i_ploidies[0] = 1;
		if(ploidy[0]=='1') r2i_ploidies[1] = 1;
		if(ploidy[0]=='2') r2i_ploidies[1] = 2;
	}
	/*allocate R2p*/
	if((statistics[0].R2p= (double **)calloc(r2i_ploidies[0],sizeof(double *))) == 0) {
		fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
		exit(1);
	}
	for(x=0;x<r2i_ploidies[0];x++) {
		if((statistics[0].R2p[x]= (double *)calloc(1*npops,sizeof(double))) == 0) {
			fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
			exit(1);
		}
	}
	
	/* Aqui el PROCESAMIENTO DE DATOS (y leer los fragmentos de datos (MS o TFA)*/
    for(first=0;first<nscaffolds;first++) {
        if(formatfile == 3) {
            chr_name = chr_name_array[first];
            wgenes = 0;
            nwindows = 0;
            /*read the file for weigth for positions, if included*/
            if( file_Wcoord[0] != '\0') {
                if( (file_wcoor = fzopen( file_Wcoord, "r", &file_wcoor_gz)) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n It is not possible to open the coordinates file %s\n", file_Wcoord);
                    exit(1);
                }
                if(read_coordinates(file_wcoor,&file_wcoor_gz,file_output,&file_output_gz,file_logerr,&file_logerr_gz,&wgenes, &nwindows,chr_name) == 0) {
                    exit(1);
                }
                window = -1;
                slide = -1;
                fzclose(file_wcoor, &file_wcoor_gz);
            }
            if( file_wps[0] != '\0' && first == 0) {
                if( (file_ws = fzopen( file_wps, "r", &file_ws_gz)) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n It is not possible to open the weights file %s\n", file_wps);
                    exit(1);
                }
                load_index_from_file(file_ws_gz.index_file_name, &index_w);
            }
        }
        li = 0;
        while(li < niterdata)
        {
            flaghky = 1; /* TODO: corrector de multiple hits... en ms format
                                         es tratado de manera diferente */
            
            /*read the ms format for each iteration*/
            /* FUNCTION TO READ MS FILES*/
            if(formatfile == 1 || formatfile == 2) 
            {
                /*read ms file*/
                if(get_msdata( file_input,&file_input_gz,file_logerr,&file_logerr_gz,
                                    &matrix_pol,&matrix_freq,&matrix_pos,
                                    &length_seg,vint_perpop_nsam,npops,int_total_nsam,
                                    length,&nmhits,matrix_mask,vector_mask,ms_svratio,
                                    &vector_priors,&npriors,&matrix_sv,
                                    outgroup_presence,force_outgroup,freq_revert,sum_sam,
                                    nsites1_pop, nsites1_pop_outg,formatfile,
                                    nsites2_pop,nsites2_pop_outg,nsites3_pop,nsites3_pop_outg,anx,bnx,anxo,bnxo,lengthamng,lengthamng_outg,
                                    include_unknown,file_mas,freq_missing_ms,kind_length,
                                    &sum_sam_mask,&length_mask,&length_mask_real,&missratio,location_missing_ms,sort_nsam)) {
                    fzprintf(file_logerr,&file_logerr_gz,"\nError processing ms data.\n");
                    exit(1);
                }
                /*
                length_al = 0.;
                for(li2=0;li2<length;li2++) {
                    length_al += vector_mask[li2];
                }
                */
                flaghky = 0;
                for(x=0;x<int_total_nsam;x++) {
                    for(y=0;y<int_total_nsam+1;y++) {
                        statistics[0].linefreq[x][y] = 0.;
                    }
                }
                
                if(include_unknown) {
                    if(!(file_mas[0] == '-'  && file_mas[1] == '1')) {
                        for(x=0;x<int_total_nsam;x++) 
                            sum_sam[x] = sum_sam_mask[x] /*- nmhits*/ ;
                        length_al = length_mask- nmhits;
                    }
                    else {/*here is counting the positions with 0 and 1 (not gaps[8 or 9] and mhits)*/
                        /*in case including gaps in ms format, it must be included the total length region!!! (not only variants!!)*/
                        /*if(freq_missing_ms > 0.) {*/
                            for(x=0;x<int_total_nsam;x++) 
                                sum_sam[x] = (sum_sam_mask[x] /*- nmhits*/);
                        length_al = length_mask- nmhits;
                        /*}
                        else {
                            length_al = length_mask - nmhits;
                        }
                        */
                    }
                }
                else {
                    for(x=0;x<int_total_nsam;x++) {
                        if (length - nmhits > 0)
                            sum_sam[x] = length - nmhits;
                        else
                            sum_sam[x] = 0;
                    }
                    if(length - nmhits > 0)
                        length_al = length - nmhits;
                    else
                        length_al = 0;
                }
                z=0;
                for(x=0;x<npops;x++) {
                    for(y=z;y<z+vint_perpop_nsam[x];y++) {
                        statistics[0].length[x] = 0;
                        for(w=0;w<4;w++) {
                            statistics[0].total_tcga[w] = 0;
                            statistics[0].tcga[x][w] = 0;
                        }
                    }
                    z += vint_perpop_nsam[x];
                }
            }

            /* FUNCTION TO READ TRANSPOSED FASTA (TFA) FILES*/
            if(formatfile == 3) { 
                if( get_tfadata(file_output,&file_output_gz,
                                file_input,&file_input_gz,&index_input,
                                file_wps,file_ws,&file_ws_gz,&index_w,
                                file_logerr,&file_logerr_gz,
                                &matrix_pol,&matrix_freq,
                                &matrix_pos,&length,&length_seg,&length_al, 
                                &length_al_real,argc,vint_perpop_nsam,npops, 
                                &svratio,&missratio,include_unknown,sum_sam,tcga, 
                                &matrix_sv,&nmhits,output,outgroup_presence, 
                                nsites1_pop,nsites1_pop_outg,
                                nsites2_pop,nsites2_pop_outg,
                                nsites3_pop,nsites3_pop_outg,
                                anx,bnx,anxo,bnxo,
                                lengthamng,lengthamng_outg,sort_nsam,
                                wgenes,nwindows,
                                first_slide,slide,window,
                                Physical_length,
                                &li,&npriors,&vector_priors,
                                &matrix_pol_tcga,
                                chr_name,first) == 0)
                {
                    /*printf("End processing input tfa data.\n");*/
                    break;
                    /*exit(1);*/
                }
                /*When the tfa file is finished then li=0, otherwise li=-1 and the loop continue*/
                z=0;
                for(x=0;x<npops;x++) {
                    for(y=z;y<z+vint_perpop_nsam[x];y++) {
                        statistics[0].length[x] = 0;
                        for(w=0;w<4;w++) {
                            statistics[0].total_tcga[w] = 0;
                            statistics[0].tcga[x][w] = 0;
                        }
                    }
                    z += vint_perpop_nsam[x];
                }
            }
            
            /*calculate number of effective nucleotides per population and tcga frequencies*/
            z=0;
            for(x=0;x<npops;x++) {
                /*
                nsites1_pop[x] -= nmhits;
                nsites1_pop_outg[x] -= nmhits;
                */
                if(outgroup_presence)
                    statistics[0].length2[x] = nsites2_pop_outg[x];
                else
                    statistics[0].length2[x] = nsites2_pop[x];
                statistics[0].anx[x] = anx[x];
                statistics[0].bnx[x] = bnx[x];
                statistics[0].anxo[x] = anxo[x];
                statistics[0].bnxo[x] = bnxo[x];
                
                for(y=z;y<z+vint_perpop_nsam[x];y++) {
                    statistics[0].length[x] += sum_sam[y];
                        for(w=0;w<4;w++) {
                        statistics[0].total_tcga[w] += tcga[y][w];
                        statistics[0].tcga[x][w] += tcga[y][w];
                    }
                }
                z += vint_perpop_nsam[x];
            }
            
            /*include lengthamng and calculations*/
            for(x=0;x<npops-!outgroup_presence;x++) {
                for(y=0;y<npops-!outgroup_presence;y++) {
                    statistics[0].lengthamng[x][y] = lengthamng[x][y];
                    statistics[0].lengthamng_outg[x][y] = lengthamng_outg[x][y];
                }
            }
            for(x=0;x<int_total_nsam;x++) {
                for(y=0;y<int_total_nsam+1;y++) {
                    statistics[0].linefreq[x][y] = 0.;
                }
            }
            
            sites_matrix = (long int *)calloc(4*(length_seg+1)*npops,sizeof(long int));
            jfd = 			(double **)calloc(npops,sizeof(double *));
            nfd = 			(int **)calloc(npops,sizeof(int *));
            
            for(x=0;x<npops;x++) {
                jfd[x] = (double *)calloc(length_seg,sizeof(double));
                nfd[x] = (int *)calloc(length_seg,sizeof(int));
                for(y=1;y<vint_perpop_nsam[x];y++) {
                    statistics[0].freq[x][y] = 0;
                }
                statistics[0].mdsd[x] = -10000;
            }
            statistics[0].total_length = length_al;
            statistics[0].total_real_length = length_al_real;
            statistics[0].total_svratio = svratio;
            statistics[0].nmhits = nmhits;

            /*calculate statistics ------------------------------------------------ */
            if( calc_sxsfss( npops,vint_perpop_nsam,matrix_pol,matrix_pos,
                                    length_seg,statistics,sites_matrix,outgroup_presence,force_outgroup) == 0) {
                fzprintf(file_logerr,&file_logerr_gz,"\nError in calc_sxsfss function.1.\n");
                exit(1);
            }
            if( jointfreqdist(npops,vint_perpop_nsam,matrix_pol,matrix_pos,
                                    length_seg,statistics,sites_matrix,jfd,nfd,outgroup_presence,force_outgroup) == 0) {
                fzprintf(file_logerr,&file_logerr_gz,"\nError in jointfreqdist function.1.\n");
                exit(1);
            }
            if( calc_piwpiafst(flaghky,formatfile,npops,vint_perpop_nsam,
                                        matrix_pol,length_seg,statistics,matrix_sv,outgroup_presence,force_outgroup) == 0) {
                fzprintf(file_logerr,&file_logerr_gz,"\nError in calc_piwpiafst function.1.\n");
                exit(1);
            }
            if( calc_freqstats(npops,vint_perpop_nsam,matrix_pol,length_seg,
                                        statistics,outgroup_presence,force_outgroup,include_unknown,n_ccov,H1frq) == 0) {
                fzprintf(file_logerr,&file_logerr_gz,"\nError in freqstats function.1.\n");
                exit(1);
            }
            if(H1frq && include_unknown == 0) 
            {
                if(calc_Toptimal_tests(npops,vint_perpop_nsam,statistics) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\nError in calc_Toptimal_tests function.1.\n");
                    exit(1);
                }
            }
            
            if(calcR2(npops,vint_perpop_nsam,matrix_pol,length_seg,statistics,ploidy)==0) {
                fzprintf(file_logerr,&file_logerr_gz,"\nError in calc_R2 function.1.\n");
                exit(1);
            }	
            
            if(calcR2p(npops,vint_perpop_nsam,matrix_pol,length_seg,statistics,sum_sam,r2i_ploidies)==0) {
                fzprintf(file_logerr,&file_logerr_gz,"\nError in calc_R2 function.1.\n");
                exit(1);
            }	
            
            /* Calculate statistics for haplotypes */
            if(int_total_nsam < SAMPLE_LARGE) {
                /*if(include_unknown==0) {*/
                if(calc_mismatch(npops,vint_perpop_nsam,matrix_pol,length_seg,statistics,ploidy) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\nError in calc_mismatch function.1.\n");
                    exit(1);
                }
                /*}*/
                if(include_unknown == 0 &&  ploidy[0] == '1')
                    {
                    if(calc_hwhafsth(npops,vint_perpop_nsam,matrix_pol,length_seg,statistics) == 0) {
                        fzprintf(file_logerr,&file_logerr_gz,"\nError in calc_hwhafsth function.1.\n");
                        exit(1);
                    }
                    if(calcFs(npops,vint_perpop_nsam,statistics)==0) { /*need calc_freqstats and calc_hwhafsth to be calculated*/
                        fzprintf(file_logerr,&file_logerr_gz,"\nError in calc_Fs function.1.\n");
                        exit(1);
                    }
                }
            }
            
            if(file_output) {
                /*
                fprintf(file_output,"Done.\n");
                fflush(file_output);
                fzprintf(file_logerr,&file_logerr_gz,"Done.\n");
                 */
                fflush(stdout);
            }
            
            /* Un monton de cosas casi-repes que volvemos a utilizar! */
            /* PERMUTATION */
            /* Only if not MS format */
            if(niter && npops > 2 /*one is outgroup, forced or not*/ && include_unknown == 0)
            {			
                /*calloc pointers in structs for permutation test*/
                if((stats_iter = (struct stats *)calloc(1,sizeof(struct stats))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }

                if((stats_iter[0].piw  = (double *)calloc(1*npops,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].pia  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                if((stats_iter[0].piT  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].piant  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                if((stats_iter[0].piTnt  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].fst  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].piwHKY  = (double *)calloc(1*npops,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                if((stats_iter[0].thetaTHKY  = (double *)calloc(1*npops,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].piaHKY  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                if((stats_iter[0].piTHKY  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].fstHKY  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].fst1all  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].hapw = (double *)calloc(1*npops,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].hapa = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].hapT = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].fsth = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].fsth1all = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].Gst  = (double *)calloc((npops*(npops-0))/2,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].K  = (double *)calloc(1*npops,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].KHKY  = (double *)calloc(1*npops,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                stats_iter[0].freq  = 0;
                if((stats_iter[0].sv = (double ***)calloc(1*npops,sizeof(double **))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                if((stats_iter[0].svT = (double ***)calloc(1*npops,sizeof(double **))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].nhpop  = (int *)calloc(1*npops,sizeof(int))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].freqh  = (long int **)calloc(1*npops,sizeof(long int *))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].tcga  = (double **)calloc(1*npops,sizeof(double *))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((stats_iter[0].length = (double *)calloc(1*npops,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                if((stats_iter[0].length2 = (double *)calloc(1*npops,sizeof(double))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                if((stats_iter[0].lengthamng  = (double **)calloc(npops,sizeof(double *))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                if((stats_iter[0].lengthamng_outg  = (double **)calloc(npops,sizeof(double *))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                for(x=0;x<npops;x++) {
                    if((stats_iter[0].lengthamng[x]  = (double *)calloc(npops,sizeof(double))) == 0) {
                        fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                        exit(1);
                    }
                    if((stats_iter[0].lengthamng_outg[x]  = (double *)calloc(npops,sizeof(double))) == 0) {
                        fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                        exit(1);
                    }
                    if((stats_iter[0].tcga [x]  = (double *)calloc(4,sizeof(double))) == 0) {
                        fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                        exit(1);
                    }
                    
                    if((stats_iter[0].freqh[x]  = (long int *)calloc(int_total_nsam,sizeof(long int))) == 0) {
                        fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                        exit(1);
                    }
                    
                    if((stats_iter[0].sv[x]     = (double **)calloc(1*npops,sizeof(double *))) == 0) {
                        fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                        exit(1);
                    }
                    if((stats_iter[0].svT[x]     = (double **)calloc(1*npops,sizeof(double *))) == 0) {
                        fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                        exit(1);
                    }
                    
                    for(y=0;y<npops;y++) {
                        if((stats_iter[0].sv[x][y] = (double *)calloc(2,sizeof(double))) == 0) {
                            fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                            exit(1);
                        }
                        if((stats_iter[0].svT[x][y] = (double *)calloc(2,sizeof(double))) == 0) {
                            fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                            exit(1);
                        }
                    }
                    
                }
                for(x=0;x<npops;x++) {
                    for(w=0;w<4;w++) {
                        stats_iter[0].tcga[x][w] = statistics[0].tcga[x][w];
                    }
                }
                stats_iter[0].total_length = length_al;
                stats_iter[0].total_svratio = svratio;
                stats_iter[0].nmhits = nmhits;
                
                if((piter = (struct probs *)calloc(1,sizeof(struct probs))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }			
                
                if((piter[0].i1  = (long int *)calloc(1*npops,sizeof(long int))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((piter[0].ih1 = (long int *)calloc(1*npops,sizeof(long int))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((piter[0].i	  = (long int *)calloc((npops*(npops-0))/2,sizeof(long int))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((piter[0].ih  = (long int *)calloc((npops*(npops-0))/2,sizeof(long int))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((piter[0].igh  = (long int *)calloc((npops*(npops-0))/2,sizeof(long int))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                
                if((piter[0].niteri1  = (long int *)calloc(1*npops,sizeof(long int))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((piter[0].niterih1 = (long int *)calloc(1*npops,sizeof(long int))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((piter[0].niteri   = (long int *)calloc((npops*(npops-0))/2,sizeof(long int))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((piter[0].niterih  = (long int *)calloc((npops*(npops-0))/2,sizeof(long int))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if((piter[0].niterigh  = (long int *)calloc((npops*(npops-0))/2,sizeof(long int))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                
                piter[0].iall	    = 0;
                piter[0].ihall      = 0;
                piter[0].ighall     = 0;
                piter[0].niteriall  = 0;
                piter[0].niterihall = 0;	
                piter[0].niterighall= 0;
                
                if((matrix_perm = (char *)calloc(length_seg*int_total_nsam,sizeof(char))) == 0) {
                    fzprintf(file_logerr,&file_logerr_gz,"\n  Error allocating memory.");
                    exit(1);
                }
                
                if(file_output && (output == 0 || output == 10)) {
                    /*fprintf(file_output,"Calculating permutation tests...\n");*/
                    fflush(file_output);
                    /*
                    fzprintf(file_logerr,&file_logerr_gz,"Calculating permutation test...\n");
                    fflush(stdout);
                    */
                }

                /*permute samples in pops:*/
                if(npops > 2)
                {
                    /* TODO: OJO que necesita inicializar las poblaciones de otra manera, segun ...*/
                    for(i=0;i<niter;i++) 
                    {				
                        /*assign the (2) groups that will be included in permutation test*/
                        nsam2[0] = vint_perpop_nsam[0];
                        nsam2[1] = int_total_nsam - vint_perpop_nsam[npops-1] - vint_perpop_nsam[0];
                        /*nsam2[2] = vint_perpop_nsam[npops-1];*/
                        
                        psam2[0] = 0;
                        psam2[1] = vint_perpop_nsam[0];
                        /*psam2[2] = psam2[1] + nsam2[1];*/
                        
                        if( permute(matrix_pol,length_seg,int_total_nsam,matrix_perm,nsam2,psam2,npops,vint_perpop_nsam[npops-1],int_total_nsam-vint_perpop_nsam[npops-1]) == 0) {
                            fzprintf(file_logerr,&file_logerr_gz,"\nError in permute function.\n");
                            exit(1);
                        }
                        /*HERE WE ASSUME THAT THE LENGTH SIZE ARE THE SAME THAN FOR ORIGINAL SAMPLES!! (OK FOR NO MISSING)*/
                        for(x=0;x<npops;x++) {
                            stats_iter[0].length[x]  = statistics[0].length[x];
                            stats_iter[0].length2[x] = statistics[0].length2[x];
                            for(yy=0;yy<npops;yy++) {
                                stats_iter[0].lengthamng[x][yy] = statistics[0].lengthamng[x][yy];
                                stats_iter[0].lengthamng_outg[x][yy] = statistics[0].lengthamng_outg[x][yy];
                            }
                        }

                        if( calc_piwpiafst(0,0,npops,vint_perpop_nsam,matrix_perm,
                                           length_seg,stats_iter,matrix_sv,
                                           outgroup_presence,force_outgroup) == 0) {
                            fzprintf(file_logerr,&file_logerr_gz,"\nError in calc_piwpiafst function.2.\n");
                            exit(1);
                        }
                        if( ploidy[0] == '1') {
                            if(calc_hwhafsth(npops,vint_perpop_nsam,matrix_perm,length_seg,stats_iter) == 0) {
                                fzprintf(file_logerr,&file_logerr_gz,"\nError in calc_hwhafsth function.2.\n");
                                exit(1);
                            }
                        }
                        
                        z=0;
                        for(x=0;x<npops-1;x++) 
                        {
                            if(stats_iter[0].fst1all[x] != -10000 && statistics[0].fst1all[x] != -10000) {
                                if(statistics[0].fst1all[x] <= stats_iter[0].fst1all[x]) 
                                    piter[0].i1[x]++;
                                piter[0].niteri1[x]++;
                                /*fprintf(file_output,"%.3f\t",stats_iter[0].fst1all[x]);*/
                            }
                            if(ploidy[0] == '1'){
                                if(stats_iter[0].fsth1all[x] != -10000 && statistics[0].fsth1all[x] != -10000) {
                                    if(statistics[0].fsth1all[x] <= stats_iter[0].fsth1all[x]) 
                                        piter[0].ih1[x]++;
                                    piter[0].niterih1[x]++;
                                }
                            }
                        }
                        
                        /*fprintf(file_output,"\n");*/
                        if(stats_iter[0].fstALL != -10000 && statistics[0].fstALL != -10000) {
                            if(statistics[0].fstALL <= stats_iter[0].fstALL) 
                                piter[0].iall++;
                            piter[0].niteriall++;
                        }
                        
                        if(ploidy[0] == '1')
                        {
                            if(stats_iter[0].fsthALL != -10000 && statistics[0].fsthALL != -10000) {
                                if(statistics[0].fsthALL <= stats_iter[0].fsthALL) 
                                    piter[0].ihall++;
                                piter[0].niterihall++;				
                            }
                            if(stats_iter[0].GstALL != -10000 && statistics[0].GstALL != -10000) {
                                if(statistics[0].GstALL <= stats_iter[0].GstALL) 
                                    piter[0].ighall++;
                                piter[0].niterighall++;				
                            }
                        }
                    }
                    if(file_output && (output == 0 || output == 10)) {
                         /*fprintf(file_output,"Permutation test one vs all Done.\n");*/
                         fflush(file_output);
                         /*
                         fzprintf(file_logerr,&file_logerr_gz,"Permutation test one vs all Done.\n");
                         fflush(stdout);
                         */
                    }
                    /* permute pairs of pops */
                    for(i=0;i<niter;i++)
                    {
                        z = 0;
                        psam2[0] = 0;
                        psam2[1] = 0;
                        
                        for(x=0;x<npops-1;x++)
                        {
                            nsam2[0] = vint_perpop_nsam[x];
                            psam2[1] = psam2[0] + vint_perpop_nsam[x];
                            for(y=x+1;y<npops-0;y++) {
                                nsam2[1] = vint_perpop_nsam[y];
                                psam2[2] = psam2[1] + vint_perpop_nsam[y];
                                nsam2[2] = vint_perpop_nsam[npops-1];
                                if( permute(matrix_pol,length_seg,int_total_nsam,matrix_perm,nsam2,psam2,npops,vint_perpop_nsam[npops-1],int_total_nsam-vint_perpop_nsam[npops-1])== 0) {
                                    fzprintf(file_logerr,&file_logerr_gz,"\nError in permute function.\n");
                                    exit(1);
                                }
                                /*HERE WE ASSUME THAT THE LENGTH SIZE ARE THE SAME THAN FOR ORIGINAL SAMPLES!! (OK FOR NO MISSING)*/
                                stats_iter[0].length[0]  = statistics[0].length[x];
                                stats_iter[0].length2[0] = statistics[0].length2[x];
                                stats_iter[0].length[1]  = statistics[0].length[y];
                                stats_iter[0].length2[1] = statistics[0].length2[y];
                                stats_iter[0].lengthamng[0][1] = statistics[0].lengthamng[x][y];
                                stats_iter[0].lengthamng_outg[0][1] = statistics[0].lengthamng_outg[x][y];
                                
                                if( calc_piwpiafst(0,0,2+1/*outg*/,nsam2,matrix_perm,length_seg,stats_iter,matrix_sv,outgroup_presence,force_outgroup) == 0) {
                                    fzprintf(file_logerr,&file_logerr_gz,"\nError in calc_piwpiafst function.2.\n");
                                    exit(1);
                                }
                                if(int_total_nsam < SAMPLE_LARGE) {
                                    if( ploidy[0] == '1') {
                                        if(calc_hwhafsth(2+1,nsam2,matrix_perm,length_seg,stats_iter) == 0) {
                                            fzprintf(file_logerr,&file_logerr_gz,"\nError in calc_hwhafsth function.2.\n");
                                            exit(1);
                                        }
                                    }
                                }
                                if( stats_iter[0].fst[0] != -10000 && statistics[0].fst[z] != -10000) {
                                    if(statistics[0].fst[z] <= stats_iter[0].fst[0])
                                        piter[0].i[z]++;
                                    piter[0].niteri[z]++;
                                }
                                if( ploidy[0] == '1' )
                                {
                                    if(stats_iter[0].fsth[0] != -10000 && statistics[0].fsth[z] != -10000) {
                                        if(statistics[0].fsth[z] <= stats_iter[0].fsth[0])
                                            piter[0].ih[z]++;
                                        piter[0].niterih[z]++;
                                    }
                                    if(stats_iter[0].Gst[0] != -10000 && statistics[0].Gst[z] != -10000) {
                                        if(statistics[0].Gst[z] <= stats_iter[0].Gst[0])
                                            piter[0].igh[z]++;
                                        piter[0].niterigh[z]++;
                                    }
                                }
                                psam2[1] += vint_perpop_nsam[y]; /*!outgroup_presence at loops*/
                                z++;
                            }
                            psam2[0] += vint_perpop_nsam[x];
                        }
                    }
                }
                if(file_output) {
                    /*
                    fprintf(file_output,"Done.\n");
                    fflush(file_output);
                    printf("Done.\n");
                     */
                    fflush(stdout);
                }
            }
            /*print results*/
            /* TODO: mas elegante el tema de force_outgroup */
            if(print_output( argc,npops,vint_perpop_nsam,file_output,&file_output_gz,file_in,file_out,
                                    gfffiles,file_GFF,subset_positions,code_name,
                                    genetic_code,length,length_seg,length_al,
                                    length_al_real,statistics,piter,niter,
                                    sites_matrix,ploidy,svratio,missratio,
                                    include_unknown,matrix_pos,jfd,nfd,output,
                                    H1frq,H0frq,nseed,file_H1f,file_H0f,vector_priors,
                                    npriors,formatfile, 
                                    outgroup_presence+force_outgroup, force_outgroup,freq_missing_ms,
                                    nsites1_pop,nsites1_pop_outg,nsites2_pop,nsites2_pop_outg,nsites3_pop,nsites3_pop_outg,
                                    li+1, matrix_pol,r2i_ploidies,matrix_pol_tcga,chr_name) == 0)
            {

                fzprintf(file_logerr,&file_logerr_gz,"\nSorry. Error in printing function.\n");
                exit(1);
            }
            /* TODO: Check cleaning the house */
            free(matrix_pol);
            if(!(formatfile == 1 || formatfile == 2)) free(matrix_pol_tcga);
            free(matrix_freq);
            free(matrix_pos);
            /*free(matrix_GC);*/
            free(matrix_sv);
            free(sites_matrix);

            for(x=0;x<npops;x++) {
                free(jfd[x]);
                free(nfd[x]);
            }
            free(jfd);
            free(nfd);
            if(/*include_unknown && */file_mas[0] == '-' && (formatfile==1 || formatfile==2))
                free(sum_sam_mask);
            
            li++;
        }
        if(formatfile == 3) {
            if(file_wcoor) free(wgenes);
        }
    }
    if(file_ws) fzclose(file_ws, &file_ws_gz);
    
    if(/*include_unknown && */file_mas[0] != '-' && (formatfile==1 || formatfile==2))
        free(sum_sam_mask);
    
    if(formatfile > 0) {
    	fzclose(file_input, &file_input_gz);
    }
	if(file_output) fclose(file_output);
	
    for(x=0;x<nscaffolds;x++)
        free(chr_name_array[x]);
    free(chr_name_array);
    
	free(sort_nsam);
	free(nsites1_pop);
	free(nsites2_pop);
	free(nsites3_pop);
	free(nsites1_pop_outg);
	free(nsites2_pop_outg);
	free(nsites3_pop_outg);
	free(anx);
	free(bnx);
	free(anxo);
	free(bnxo);
    for(x=0;x<npops;x++) {
        free(lengthamng[x]);
        free(lengthamng_outg[x]);
    }
    free(lengthamng);
    free(lengthamng_outg);
	free(vint_perpop_nsam);
	
	if(H1frq) {
		for(x=0;x<npf;x++)
			free(freqspH1[x]);
		free(freqspH1);
		free(thetaH1);
		for(x=0;x<npf;x++)
			free(freqspH0[x]);
		free(freqspH0);
		free(thetaH0);
	}
	
	free(statistics[0].Sanc);
	free(statistics[0].piw);
	free(statistics[0].pia);
	free(statistics[0].piT);
	free(statistics[0].piant);
	free(statistics[0].piTnt);
	free(statistics[0].fst);
    free(statistics[0].piwHKY);
	free(statistics[0].piaHKY);
	free(statistics[0].piTHKY);
	free(statistics[0].fstHKY);
	free(statistics[0].fst1all);
	free(statistics[0].hapw);
	free(statistics[0].hapa);
	free(statistics[0].hapT);
	free(statistics[0].fsth);
	free(statistics[0].fsth1all);
	free(statistics[0].Gst);
	
	free(statistics[0].S);
	free(statistics[0].thetaS);
	free(statistics[0].thetaT);
	free(statistics[0].So);
	free(statistics[0].thetaSo);
	free(statistics[0].thetaTo);
	free(statistics[0].thetaTHKY);
	free(statistics[0].thetaFL);
	free(statistics[0].thetaFW);
	free(statistics[0].thetaL);
	free(statistics[0].thetaSA);
	free(statistics[0].thetaTA);
	free(statistics[0].K);
	free(statistics[0].KHKY);
	free(statistics[0].Dtaj);
	free(statistics[0].Dfl);
	free(statistics[0].Ffl);
	free(statistics[0].Hnfw);
	free(statistics[0].Ez);
    free(statistics[0].Yach);
    free(statistics[0].FH);
	free(statistics[0].R2);
	free(statistics[0].Fs);
	free(statistics[0].nhpop);
	free(statistics[0].length);
	free(statistics[0].total_tcga);
	free(statistics[0].ToH0_ii);
	free(statistics[0].To_ii);
	free(statistics[0].To_00);
	free(statistics[0].To_i0);
	free(statistics[0].ToH0_00);
	free(statistics[0].To_Qc_ii);
	free(statistics[0].To_Qw_ii);
	free(statistics[0].To_Lc_ii);
	
	free(statistics[0].Rm);
	free(statistics[0].ZnA);

	free(statistics[0].mdsd);
	free(statistics[0].mdg1);
	free(statistics[0].mdg2);

	free(statistics[0].anx);
	free(statistics[0].bnx);
	free(statistics[0].anxo);
	free(statistics[0].bnxo);

	for(x=0;x<npops;x++) {
		free(statistics[0].freq[x]);
		free(statistics[0].freqh[x]);
		free(statistics[0].tcga[x]);
		free(statistics[0].mdw[x]);
        free(statistics[0].lengthamng[x]);
        free(statistics[0].lengthamng_outg[x]);
		for(y=0;y<npops;y++) {
			free(statistics[0].sv[x][y]);
			free(statistics[0].svT[x][y]);
		}
		free(statistics[0].sv[x]);
		free(statistics[0].svT[x]);
	}
    free(statistics[0].length2);
    free(statistics[0].lengthamng);
    free(statistics[0].lengthamng_outg);

	for(x=0;x<int_total_nsam;x++) {free(statistics[0].linefreq[x]);}
	for(x=0;x<r2i_ploidies[0];x++) {free(statistics[0].R2p[x]);}
	free(statistics[0].R2p);
	
	free(statistics[0].freq);
	free(statistics[0].freqh);
	free(statistics[0].tcga);
	free(statistics[0].sv);
	free(statistics[0].svT);
	free(statistics[0].mdw);
	free(statistics[0].linefreq);
	free(statistics);
	free(sum_sam);
    for( x=0; x<int_total_nsam+(!outgroup_presence); x++)
        free(tcga[x]);
	free(tcga);
	if(!(formatfile == 0 || ( formatfile == 3 && ((slide == 0 && window == 0) && file_Wcoord[0]=='\0')))) {
		free(matrix_mask);
		free(vector_mask);
	}
	
	if(niter && npops > 2) {
		free(matrix_perm);
        
		free(stats_iter[0].piw);
        free(stats_iter[0].pia);
        free(stats_iter[0].piT);
        free(stats_iter[0].piant);
        free(stats_iter[0].piTnt);
		free(stats_iter[0].fst);
        free(stats_iter[0].piwHKY);
        free(stats_iter[0].thetaTHKY);
		free(stats_iter[0].piaHKY);
		free(stats_iter[0].piTHKY);
		free(stats_iter[0].fstHKY);
		free(stats_iter[0].fst1all);
		free(stats_iter[0].hapw);
		free(stats_iter[0].hapa);
		free(stats_iter[0].hapT);
		free(stats_iter[0].fsth);
		free(stats_iter[0].fsth1all);
		free(stats_iter[0].Gst);
		free(stats_iter[0].K);
		free(stats_iter[0].KHKY);
        free(stats_iter[0].length);
        free(stats_iter[0].length2);
		
		for(x=0;x<npops;x++) {
			free(stats_iter[0].freqh[x]);
			free(stats_iter[0].tcga[x]);
			for(y=0;y<npops;y++) {
				free(stats_iter[0].sv[x][y]);
				free(stats_iter[0].svT[x][y]);
			}
			free(stats_iter[0].sv[x]);
            free(stats_iter[0].svT[x]);
            free(stats_iter[0].lengthamng[x]);
            free(stats_iter[0].lengthamng_outg[x]);
		}
		free(stats_iter[0].freqh);
		free(stats_iter[0].nhpop);
		free(stats_iter[0].tcga);
		free(stats_iter[0].sv);
		free(stats_iter[0].svT);
        free(stats_iter[0].lengthamng);
        free(stats_iter[0].lengthamng_outg);

		free(stats_iter);
		
		free(piter[0].i1);
		free(piter[0].ih1);
		free(piter[0].i);
		free(piter[0].ih);
		free(piter[0].igh);
		
		free(piter[0].niteri1);
		free(piter[0].niterih1);	
		free(piter[0].niteri);
		free(piter[0].niterih);
		free(piter[0].niterigh);
		
		if(npriors) free(vector_priors);
		
		free(piter);
	}
	free(f);
    free(r2i_ploidies);
    
    fzprintf(file_logerr,&file_logerr_gz,"\nProgram Ended\n");
    fzclose(file_logerr, &file_logerr_gz);

	exit(0);
}

void usage(void)
{
    printf(MSTATSPOP);
    printf("Flags:\n");
    printf("      -f [input format file: fasta, tfa, ms]\n");
    printf("      -i [path and name of the input file]\n");
    printf("      -o [output format file: 0 (extended),\n");
    printf("                              1 (single line/window),\n");
    printf("                              2 (single line SFS/window),\n");
    printf("                              3 (dadi-like format),\n");
    printf("                              4 (single line pairwise distribution)\n");
    printf("                              5 (single line freq. variant per line/window)\n");
    printf("                              6 (SNP genotype matrix)\n");
    printf("                             10 (full extended)]\n");
    printf("      -N [#_pops] [#samples_pop1] ... [#samples_popN]\n");
    printf("      -n [name of a single scaffold to analyze. For tfa can be a list separated by commas(ex. -n chr1,chr2,chr3]\n");
    printf("      -T [path and name of the output file]. DEFAULT stdout.\n");
    printf("   OPTIONAL GENERAL PARAMETERS:\n");
    printf("      -G [outgroup (0/1)] (last population). DEFAULT 0.\n");
    printf("      -u [include unknown positions (0/1)].  DEFAULT 0.\n");
    printf("      -A [Alternative Spectrum File (Only for Optimal Test): alternative_spectrum for each population (except outg)\n");
    printf("          File format: (average absolute values) header plus fr(0,1) fr(0,2) ... fr(0,n-1) theta(0)/nt,\n");
    printf("          fr(1,1) fr(1,2) ... fr(1,n-1) theta(1)/nt...]\n");
    printf("      -S [Null Spectrum File (only if -a is defined): null_spectrum for each population (except outg).\n");
    printf("          (average absolute values) header plus fr(0,1) fr(0,2) ... fr(0,n-1) theta(0)/nt,\n");
    printf("          fr(1,1) fr(1,2) ... fr(1,n-1) theta(1)/nt...]. DEFAULT SNM.\n");
    printf("      -P [Only for Calculation of R2_p: first value is the number of values to include, \n");
    printf("                       next are the number of lines to consider. ex: -P 6 1 2 4 8 16 64]\n");
    printf("    Optional Parameters for fasta and tfa input files:\n");
    printf("      -O [#_nsam] [number order of first sample, number 0 is the first sample] [second sample] ...etc. up to nsamples.\n");
    printf("         DEFAULT current order.\n");
    printf("      -t [# permutations per window (H0: Fst=0). Only available with option -u 0]. DEFAULT 0.\n");
    printf("      -s [seed]. DEFAULT 123456.\n");
    printf("   PARAMETERS FOR TFASTA INPUT (-f tfa): 'SLIDING WINDOW ANALYSIS OF EMPIRICAL DATA'\n");
    printf("      -w [window size].\n");
    printf("    Optional:\n");
    printf("      -z [slide size (must be a positive value)]. DEFAULT window size.\n");
    printf("      -Z [first window size displacement [for comparing overlapped windows])]. DEFAULT 0.\n");
    printf("      -Y [define window lengths in 'physical' positions (1) or in 'effective' positions (0)]. DEFAULT 1.\n");
    printf("      -W [file with the coordinates of each window [init end] (instead options -w and -z).\n");
    printf("         DEFAULT one whole window.\n");
    printf("      -E [input file with weights for positions:\n");
    printf("         include three columns with a header,\n");
    printf("         first the physical positions (1...end),\n");
    printf("         second the weight for positions and\n");
    printf("         third a boolean weight for the variant (eg. syn variant in nsyn counts is 0.000)].\n");
    printf("         DEFAULT all 1.000\n");
    printf("   PARAMETERS FOR MS INPUT (-f ms):'SIMULATION ANALYSIS OF A SINGLE REGION'\n");
    printf("      -l [length]\n");
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
    /*printf("\t-y [in case missing: number of comparisons for approximate calculation of covariance. 0 if rough approach]\n");*//*to eliminate*/
    /*
     printf("\t-T [count only transitions (not for ms format)\n"); NOT DONE
     printf("\t-V [count only transversions (not for ms format)\n"); NOT DONE
     printf("\t-G [count only G/C mutations (not for ms format)\n"); NOT DONE
     printf("\t-A [count only A/T mutations (not for ms format)\n"); NOT DONE
     */
}
