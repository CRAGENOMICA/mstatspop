//
//  get_tfadata.c
//  xcode_project
//
//  Created by Sebastian Ramos-Onsins on 09/08/15.
//
//

#include "get_tfadata.h"
#include "log.h"
#include "tfasta.h"
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
    // long int slide,  /**/
    // long int window, /**/
    // int Physical_length,
    long int *li /**/,
    int *npriors,
    double **vector_priors,
    char **matrix_pol_tcga,
    char *chr_name,
    int first,
    mstatspop_args_t *args)
{
  static bool first_time = true;
  // static tfasta_file tfasta;
  // if(first_time){
  //   memset(&tfasta, 0, sizeof(tfasta));
  //   int ret_status = init_tfasta_file(&tfasta, args->file_in);
  //   if(ret_status != 1){
  //     log_fatal("Error initializing tfasta file");
  //     return 0;
  //   }
  //   first_time = false;
  // }

  char *DNA_matr2 = 0;
  char **names2 = 0; /* limit of the name of the lines to 50 characters. be careful ... */

  /*	double *matrix_sizepos = 0; */   /*size of each position, 0 not include, >0 include, Syn/Nsyn positions are allowed*/
  /*    double *matrix_segrpos = 0; */ /*always 1 except in those biallelic positions that are not desired (e.g., choose syn but we do not want nsyn)*/

  long int *mhitbp; /*vector of positions where mhits are*/
  /* long int count; */
  long int xx;
  /* int c; */
  int n_sam /*,ns*/;
  /*long int n_sit;*/
  /* int nseq; */
  /* int maxsam; */
  /* int n_excl; */
  int n_samp;
  long long n_site;
  char ploidy[2];

  int x;
  /*    static long int maxsites = 0;*/

  int nsamtot;
  int nsamuser_eff;
  int flag_change_sort; /*in case the order of samples is not consecutive*/

  // Usage of DNA_matr:
  // DNA_matr is a static character array used to store the DNA sequences read from a tfasta file.
  // It is initialized to 0 and its size is determined dynamically based on the number of samples and sites.
  // The DNA sequences are stored in a contiguous manner, with each sequence occupying n_site characters.
  // The total size of DNA_matr is n_site * n_samp characters.

  // Infered structure of DNA_matr:
  // DNA_matr is a 2-dimensional character array, where each row represents a DNA sequence and each column represents a site.
  // The structure can be visualized as follows:
  // DNA_matr[0][0]  DNA_matr[0][1]  DNA_matr[0][2]  ...  DNA_matr[0][n_site-1]
  // DNA_matr[1][0]  DNA_matr[1][1]  DNA_matr[1][2]  ...  DNA_matr[1][n_site-1]
  // ...
  // DNA_matr[n_samp-1][0]  DNA_matr[n_samp-1][1]  DNA_matr[n_samp-1][2]  ...  DNA_matr[n_samp-1][n_site-1]

  // Usage of names:
  // names is a static character array used to store the names of the samples.
  // It is initialized to 0 and its size is determined dynamically based on the number of samples.
  // Each name is limited to 50 characters.
  // The names are stored as an array of character pointers, where each pointer points to a character array representing a name.
  // The total size of names is 128 * 50 characters.

  // Infered structure of names:
  // names is a 1-dimensional character array, where each element represents a name.
  // The structure can be visualized as follows:
  // names[0]  names[1]  names[2]  ...  names[127]

  
  static char *DNA_matr = 0;
  
  /* limit of the name of the lines to 50 characters. be careful ... */
  // static char **names = 0; 
  
  static char chr_name2[MSP_MAX_NAME];

  static long int init_site = 1;
  static long wc = 0;
  long int beg, end;
  double wl;

  /*tfasta windows and weights*/
  /*wV: weight at variant (effect sizes)*/ /*not yet functional although we can recover*/
  /*wP: weight for each position*/
  /*wPV: weight for the variant at each position*/
  long int wlimit_end = 0; /*initial value*/
  long int n_sitesw = 0;

  double *wV;
  double *wP;
  double *wPV;
  double window_size;
  int weight_window;

  if ((wP = (double *)calloc(1000, sizeof(double))) == 0)
  {
    // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.14 \n");
    log_fatal("Error: memory not reallocated, wP get_obsdata.14");
    return (0);
  }
  if ((wPV = (double *)calloc(1000, sizeof(double))) == 0)
  {
    // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.25 \n");
    log_fatal("Error: memory not reallocated, wPV get_obsdata.25");
    return (0);
  }
  if ((wV = (double *)calloc(1000, sizeof(double))) == 0)
  {
    // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.23 \n");
    log_fatal("Error: memory not reallocated, wV get_obsdata.23");
    return (0);
  }

  ploidy[0] = '1';
  ploidy[1] = '\0';
  /* count = 0; */
  /* c = 0; */
  n_sam = 0;
  /*n_sit = 0;*/
  /* nseq  = 0; */
  /* maxsam= 128; */
  n_samp = 0;
  n_site = 0;
  /* n_excl= 0; */

  nsamtot = 0;
  for (x = 0; x < args->npops; x++)
  {
    nsamtot += args->vint_perpop_nsam[x];
  }
  nsamuser_eff = (nsamtot - !args->outgroup_presence);

  // allocate memory for names moved to tfasta struct & initialized in init_tfasta_file
  // if (names == 0)
  // { /* only initialize once. Check */
  //   if ((names = (char **)calloc(128, sizeof(char *))) == 0)
  //   {
  //     // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.1 \n");
  //     log_fatal("Error: memory not reallocated, names get_obsdata.1");
  //     free(wP);
  //     free(wPV);
  //     free(wV);
  //     return (0);
  //   }
  //   for (x = 0; x < 128; x++)
  //   {
  //     if ((names[x] = (char *)calloc(50, sizeof(char))) == 0)
  //     {
  //       // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.2 \n");
  //       log_fatal("Error: memory not reallocated, names get_obsdata.2");
  //       free(wP);
  //       free(wPV);
  //       free(wV);
  //       return (0);
  //     }
  //   }
  // }
  if (DNA_matr == 0) 
  {
    if ((DNA_matr = (char *)calloc(10000, sizeof(char))) == 0)
    {
      // for (x = 0; x < 128; x++)
      //   free(names[x]);
      // free(names);

      // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.3 \n");
      log_fatal("Error: memory not reallocated, DNA_matr get_obsdata.3");
      free(wP);
      free(wPV);
      free(wV);
      return (0);
    }
  }
  if (strcmp(chr_name2, chr_name) != 0) // A hard reset for changing chromosomes
  {
    init_site = 1;
    wc = 0;
    strcpy(chr_name2, chr_name);
  }

  if (init_site == 1)
    init_site += args->first_slide;

  /*FIND THE WINDOW OF POSITIONS TO ANALYZE*/
  if (nwindows == 0)
  {
    if (args->Physical_length == 1)
    {
      weight_window = 0; /*physical positions*/
      window_size = (double)args->window;
      beg = init_site;                    /*the initial position of the window*/
      end = init_site + args->window - 1; /*the final position of the window*/
      init_site += args->slide;           /*init for next window (static variable)*/
      if (wtfasta != NULL)
      {
        if (read_weights_positions_file(
					wtfasta,
          //file_ws, 
          //file_ws_gz, 
          //index_w,
          // file_logerr,file_logerr_gz,
          &wP, 
          &wPV, 
          &wV, 
          &wlimit_end,
          beg, 
          &window_size, 
          &n_sitesw, 
          weight_window,
          chr_name, 
          first, 
          *length) == 0)
        {
          // fprintf(file_logerr,"Error processing weighting file %s\n", file_wps);
          log_error("Error processing weighting file %s", wtfasta->wtfasta_fname);
          free(wP);
          free(wPV);
          free(wV);
          exit(1);
        }
      }
      else
      {
        n_sitesw = end - beg + 1 + 1; /*the position 0 does not count: we translate 1 right*/
        if ((wP = (double *)realloc((double *)wP, n_sitesw * sizeof(double))) == 0)
        {
          // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.14 \n");
          log_fatal("Error: memory not reallocated, wP get_obsdata.14");
          free(wP);
          free(wPV);
          free(wV);
          return (0);
        }
        if ((wPV = (double *)realloc((double *)wPV, n_sitesw * sizeof(double))) == 0)
        {
          // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.25 \n");
          log_fatal("Error: memory not reallocated, wPV get_obsdata.25");
          free(wP);
          free(wPV);
          free(wV);
          return (0);
        }
        if ((wV = (double *)realloc((double *)wV, n_sitesw * sizeof(double))) == 0)
        {
          // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.23 \n");
          log_fatal("Error: memory not reallocated, wV get_obsdata.23");
          free(wP);
          free(wPV);
          free(wV);
          return (0);
        }
        for (xx = 0; xx < n_sitesw; xx++)
        {
          wP[xx] = wPV[xx] = wV[xx] = (float)1;
        }
      }
    }
    else
    {                    /*if(Physical_length == 0), that is weighting windows*/
      weight_window = 1; /*weight positions*/
      window_size = (double)args->window;
      beg = init_site; /*the initial position of the window*/
      if (wtfasta == 0)
      {
        // fprintf(file_logerr,"Error: no weights file\n");
        log_error("Error: no weights file");
        free(wP);
        free(wPV);
        free(wV);
        exit(1);
      }
      if (read_weights_positions_file(
				wtfasta,
				//file_ws, file_ws_gz, index_w,
                                      // file_logerr,file_logerr_gz,
                                      &wP, &wPV, &wV, &wlimit_end,
                                      beg, &window_size, &n_sitesw, weight_window,
                                      chr_name, first, *length) == 0)
      {
        // fprintf(file_logerr,"Error processing weights file\n");
        log_error("Error processing weights file");
        free(wP);
        free(wPV);
        free(wV);
        exit(1);
      }
      end = wlimit_end; /*the final position of the window (static variable)*/

      if (args->slide < args->window)
      {
        wl = 0.0;
        xx = 0;
        while (wl < args->slide && xx < n_sitesw)
        {
          wl += wP[xx] * wPV[xx];
          xx++;
        }
        init_site = beg + xx;
      }
      else
      {
        if (args->slide == args->window)
        {
          init_site = end + 1;
        }
        else
        {
          init_site = end + 1;
          window_size = args->slide - args->window;
          if (read_weights_positions_file(
						wtfasta,
						//file_ws, file_ws_gz, index_w,
                                          // file_logerr,file_logerr_gz,
                                          &wP, &wPV, &wV, &wlimit_end,
                                          init_site, &window_size, &n_sitesw, weight_window,
                                          chr_name, first, *length) == 0)
          {
            // fprintf(file_logerr,"Error processing weights file\n");
            log_error("Error processing weights file");
            free(wP);
            free(wPV);
            free(wV);
            exit(1);
          }
          init_site = wlimit_end; /*init for next window*/
        }
      }
    }
  }
  else
  { /*(nwindows>0)*/    /*that is, using coordinate positions*/
    beg = wgenes[wc++]; /*the initial position of the window*/
    end = wgenes[wc++]; /*the final position of the window*/
    weight_window = 0;  /*physical positions*/
    window_size = (double)end - beg + 1;
    if (wtfasta != NULL)
    {
      if (read_weights_positions_file(
				wtfasta,
				//file_ws, file_ws_gz, index_w,
                                      // file_logerr,file_logerr_gz,
                                      &wP, &wPV, &wV, &wlimit_end,
                                      beg, &window_size, &n_sitesw, weight_window,
                                      chr_name, first, *length) == 0)
      {
        // fprintf(file_logerr,"Error processing weights file\n");
        log_error("Error processing weights file");
        free(wP);
        free(wPV);
        free(wV);
        exit(1);
      }
    }
    else
    {
      n_sitesw = end - beg + 1 + 1; /*the position 0 does not count: we translate 1 right*/
      if ((wP = (double *)realloc((double *)wP, n_sitesw * sizeof(double))) == 0)
      {
        // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.14 \n");
        log_fatal("Error: memory not reallocated, wP get_obsdata.14");
        free(wP);
        free(wPV);
        free(wV);
        return (0);
      }
      if ((wPV = (double *)realloc((double *)wPV, n_sitesw * sizeof(double))) == 0)
      {
        // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.25 \n");
        log_fatal("Error: memory not reallocated, wPV get_obsdata.25");
        free(wP);
        free(wPV);
        free(wV);
        return (0);
      }
      if ((wV = (double *)realloc((double *)wV, n_sitesw * sizeof(double))) == 0)
      {
        // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.23 \n");
        log_fatal("Error: memory not reallocated, wV get_obsdata.23");
        free(wP);
        free(wPV);
        free(wV);
        return (0);
      }
      for (xx = 0; xx < n_sitesw; xx++)
      {
        wP[xx] = wPV[xx] = wV[xx] = (float)1;
      }
    }
  }

  if (*vector_priors == 0)
  {
    *npriors = 2;
    if ((*vector_priors = (double *)calloc((long int)*npriors, sizeof(double))) == 0)
    {
      // fprintf(file_logerr,"Error: memory not allocated. get_tfadata.01");
      log_fatal("Error: memory not allocated, vector_priors get_tfadata.01");
      free(wP);
      free(wPV);
      free(wV);
      return (0);
    }
  }
  vector_priors[0][0] = (double)beg;
  vector_priors[0][1] = (double)end;

  /*READ TFASTA FILE*/
  /*define the init and the end site first! use slide and window if necessary: use also the weights if necessary*/
  /*detect the end of the file! (*li=0,otherwise *li=-1)*/

  // if (
  //     (
  //         x = function_read_tfasta(
  //             file_input,
  //             input_gz,
  //             index_input,
  //             // file_logerr,file_logerr_gz,
  //             beg,
  //             end,
  //             &n_sam,
  //             &n_site,
  //             &names,
  //             &DNA_matr,
  //             matrix_pol_tcga,
  //             chr_name,
  //             first,
  //             *length)) == 0)
  if(
    read_tfasta_DNA(tfasta, chr_name, beg, end, &n_sam, &n_site, &DNA_matr) == 0
  )
  {
    // fprintf(file_logerr,"Unable reading tfasta file\n");
    log_error("Unable reading tfasta file");
    // for (x = 0; x < n_sam; x++)
    //   free(tfasta.names[x]);
    // free(tfasta.names);
    close_tfasta_file(tfasta);
    free(wP);
    free(wPV);
    free(wV);
    free(DNA_matr);
    exit(1);
  }

  /*check that wP, wPV and wV (n_sitesw) have at least n_site positions!!! (if not, add more positions to the vector with zero values)*/
  if (n_site > n_sitesw)
  {
    if ((wP = (double *)realloc((float *)wP, n_site * sizeof(double))) == 0)
    {
      // fprintf(file_output,"\nError: memory not reallocated. get_obsdata.14 \n");
      log_fatal("Error: memory not reallocated, wP get_obsdata.14");
      exit(1);
    }
    if ((wPV = (double *)realloc((float *)wPV, n_site * sizeof(double))) == 0)
    {
      // fprintf(file_output,"\nError: memory not reallocated. get_obsdata.25 \n");
      log_fatal("Error: memory not reallocated, wPV get_obsdata.25");
      exit(1);
    }
    if ((wV = (double *)realloc((float *)wV, n_site * sizeof(double))) == 0)
    {
      // fprintf(file_output,"\nError: memory not reallocated. get_obsdata.23 \n");
      log_fatal("Error: memory not reallocated, wV get_obsdata.23");
      exit(1);
    }
    for (xx = n_sitesw; xx < n_site; xx++)
    {
      wP[xx] = wPV[xx] = wV[xx] = (double)0;
    }
    n_sitesw = n_site;
  }
  // li == 0 means that the end of the file has been reached
  // li == -1 means still more data to read !
  // x == -1 means that the end of the file has been reached
  if (x == -1 || (nwindows > 0 && wc == 2 * nwindows))
    *li = 0;
  else
    *li = -1;

  n_samp = n_sam;
  *length = n_site;

  if (n_samp < nsamuser_eff)
  {
    free(wP);
    free(wPV);
    free(wV);
    return (0);
  }
  if (n_samp == 0 || n_site == 0)
  {
    free(wP);
    free(wPV);
    free(wV);
    return (0);
  }
  else
  {
    /*modify the order of samples using option flag O*/
    flag_change_sort = 0;
    for (x = 0; x < nsamuser_eff; x++)
    {
      if (sort_nsam[x] != x)
      {
        flag_change_sort = 1;
        break;
      }
    }
    if (flag_change_sort == 1)
    {
      /*define duplicated matr*/
      if ((DNA_matr2 = (char *)calloc(n_site * (long long)n_samp, sizeof(char))) == 0)
      {
        // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.23d \n");
        log_fatal("Error: memory not reallocated, DNA_matr2 get_obsdata.23d");
        // for (x = 0; x < n_samp; x++)
        //   free(names[x]);
        // free(names);
        close_tfasta_file(tfasta);
        free(wP);
        free(wPV);
        free(wV);
        free(DNA_matr);
        return (0);
      }
      if ((names2 = (char **)calloc(n_samp, sizeof(char *))) == 0)
      {
        // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.1s2 \n");
        log_fatal("Error: memory not reallocated, names2 get_obsdata.1s2");
        // for (x = 0; x < n_samp; x++)
        //   free(names[x]);
        // free(names);
        close_tfasta_file(tfasta);
        free(wP);
        free(wPV);
        free(wV);
        free(DNA_matr);
        return (0);
      }
      for (x = 0; x < n_samp; x++)
      {
        if ((names2[x] = (char *)calloc(50, sizeof(char))) == 0)
        {
          // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.22 \n");
          log_fatal("Error: memory not reallocated, names2 get_obsdata.22");
          // for (x = 0; x < n_samp; x++)
          //   free(names[x]);
          // free(names);
          close_tfasta_file(tfasta);
          free(wP);
          free(wPV);
          free(wV);
          free(DNA_matr);
          return (0);
        }
        /*copy duplicated data*/
        // strncpy(names2[x], names[x], 50);
        strncpy(names2[x], tfasta->names[x], 50);
      }
      strncpy(DNA_matr2 + (long long)n_site * (long long)0, DNA_matr + (long long)n_site * (long long)0, (long long)n_site * n_samp);
      // for(x=0;x<n_samp;x++) {
      //     strncpy(DNA_matr2+(long long)n_site*(long long)x,DNA_matr+(long long)n_site*(long long)x,n_site);
      // }
      /*end define and duplicate*/

      /*include data in *DNA_matr and in *names[] in the correct order*/
      for (x = 0; x < nsamuser_eff; x++)
      {
        strncpy(DNA_matr + (long long)n_site * (long long)x, DNA_matr2 + (long long)n_site * (long long)sort_nsam[x], n_site);
        // strncpy(names[x], names2[sort_nsam[x]], 50);
        strncpy(tfasta->names[x], tfasta->names[sort_nsam[x]], 50);
      }
      /*delete duplicated matr*/
      for (x = 0; x < n_samp; x++)
        free(names2[x]);
      free(names2);
      names2 = 0;
      free(DNA_matr2);

      /*erase lines no used*/
      if (nsamuser_eff > n_samp)
      {
        // fprintf(file_logerr,"Error: too low samples in the file according to defined in -N flag.\n");
        log_error("Error: too low samples in the file according to defined in -N flag.");
        // for (x = 0; x < n_sam; x++)
        //   free(names[x]);
        // free(names);
        close_tfasta_file(tfasta);
        free(wP);
        free(wPV);
        free(wV);
        free(DNA_matr);
        return (0);
      }
    }
    /*end option flag O*/
    if (nsamuser_eff > 32167)
    {
      // fprintf(file_logerr,"Error: too much samples. Only 32167 samples per loci are allowed.\n");
      log_error("Error: too much samples. Only 32167 samples per loci are allowed.");
      // for (x = 0; x < n_sam; x++)
      //   free(names[x]);
      // free(names);
      close_tfasta_file(tfasta);
      free(wP);
      free(wPV);
      free(wV);
      free(DNA_matr);
      return (0);
    }
    /*init matrix_sizepos*/
    /*
            if(matrix_sizepos == 0) {
                if((matrix_sizepos = (double *)malloc(n_site*sizeof(double))) == 0) {
                    fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.2");
                    for(x=0;x<n_samp;x++) free(names[x]); free(names);
                    free(DNA_matr);
                    return(0);
                }
                if((matrix_segrpos = (double *)malloc(n_site*sizeof(double))) == 0) {
                    fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.2");
                    for(x=0;x<n_samp;x++) free(names[x]); free(names);
                    free(DNA_matr);
                    free(matrix_sizepos);
                    return(0);
                }
                maxsites = n_site;
            }
            else{
                if(n_site > maxsites) {
                    if((matrix_sizepos = (double *)realloc(matrix_sizepos,n_site*sizeof(double))) == 0) {
                        fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.2b");
                        for(x=0;x<n_samp;x++) free(names[x]); free(names);
                        free(DNA_matr);
                        return(0);
                    }
                    if((matrix_segrpos = (double *)realloc(matrix_segrpos,n_site*sizeof(double))) == 0) {
                        fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.2b");
                        for(x=0;x<n_samp;x++) free(names[x]); free(names);
                        free(DNA_matr);
                        free(matrix_sizepos);
                        return(0);
                    }
                    maxsites = n_site;
                }
            }
            for(xx=0;xx<n_site;xx++) {
                matrix_sizepos[xx] = (double)1;
                matrix_segrpos[xx] = (double)1;
            }
    */
    if (args->outgroup_presence == 0)
    {
      if ((DNA_matr = (char *)realloc(DNA_matr, (long long)n_site * (nsamuser_eff + !args->outgroup_presence) * sizeof(char))) == 0)
      {
        // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.23a \n");
        log_fatal("Error: memory not reallocated, DNA_matr get_obsdata.23a");
        // for (x = 0; x < n_samp; x++)
        //   free(names[x]);
        // free(names);
        close_tfasta_file(tfasta);
        free(wP);
        free(wPV);
        free(wV);
        free(DNA_matr);
        //                free(matrix_sizepos);
        //                free(matrix_segrpos);
        return (0);
      }
      int ns;
      int f1, f2, f3, f4, f0;
      int nf1, nf2, nf3, nf4, nfm, nf0, countf;
      for (xx = 0; xx < n_site; xx++)
      {
        /*define f1 for a given nt*/
        ns = nfm = 0;
        while (ns < nsamuser_eff - 1 && DNA_matr[(long long)n_site * (unsigned long)ns + xx] > '4')
        {
          ns++;
          nfm++;
        }
        if (ns < nsamuser_eff - 1)
        {
          f1 = (int)DNA_matr[(long long)n_site * (unsigned long)ns + xx];
          nf1 = 1;
          /*count the freqs of all nt*/
          nf2 = nf3 = nf4 = 0;
          f2 = f3 = f4 = 0;
          while (ns < nsamuser_eff - 1)
          {
            ns++;
            if ((int)DNA_matr[(long long)n_site * (unsigned long)ns + xx] == f1)
              nf1 += 1;
            if (DNA_matr[(long long)n_site * (unsigned long)ns + xx] > '4')
              nfm += 1;
            if ((int)DNA_matr[(long long)n_site * (unsigned long)ns + xx] != f1 &&
                DNA_matr[(long long)n_site * (unsigned long)ns + xx] < '5')
            {
              if (f2 == 0 || (int)DNA_matr[(long long)n_site * (unsigned long)ns + xx] == f2)
              {
                f2 = (int)DNA_matr[(long long)n_site * (unsigned long)ns + xx];
                nf2 += 1;
              }
              else
              {
                if (f3 == 0 || (int)DNA_matr[(long long)n_site * (unsigned long)ns + xx] == f3)
                {
                  f3 = (int)DNA_matr[(long long)n_site * (unsigned long)ns + xx];
                  nf3 += 1;
                }
                else
                {
                  if (f4 == 0 || (int)DNA_matr[(long long)n_site * (unsigned long)ns + xx] == f4)
                  {
                    f4 = (int)DNA_matr[(long long)n_site * (unsigned long)ns + xx];
                    nf4 += 1;
                  }
                }
              }
            }
          }
          /*choose one of the two highest frequencies, not the lowest if more than 2 (give '5')*/
          countf = 1;
          f0 = f1;
          nf0 = nf1;
          if (nf0 == nf2)
          {
            countf += 1;
          }
          if (nf0 < nf2)
          {
            f0 = f2;
            nf0 = nf2;
            countf = 1;
          }
          if (nf0 == nf3)
          {
            countf += 1;
          }
          if (nf0 < nf3)
          {
            f0 = f3;
            nf0 = nf3;
            countf = 1;
          }
          if (nf0 == nf4)
          {
            countf += 1;
          }
          if (nf0 < nf4)
          {
            f0 = f4;
            nf0 = nf4;
            countf = 1;
          }

          if (countf < 3)
            DNA_matr[(unsigned long)n_site * (long long)(nsamuser_eff) + xx] = (char)f0;
          else
            DNA_matr[(unsigned long)n_site * (long long)(nsamuser_eff) + xx] = '5';
        }
        else
        {
          DNA_matr[(unsigned long)n_site * (long long)(nsamuser_eff) + xx] = '5';
        }
        // ns = 0;
        // while(ns < nsamuser_eff-1 && DNA_matr[(long long)n_site*(unsigned long)ns+xx] > '4') ns++;
        // DNA_matr[(unsigned long)n_site*(long long)(nsamuser_eff)+xx] = DNA_matr[(unsigned long)n_site*(unsigned long)ns+xx];
      }
      nsamuser_eff += 1;
      /*we forced the invented outgroup without gaps or uncertainties, if possible*/
      // for(xx=0;xx<n_site;xx++) {
      //     int ns = 0;
      //     while(ns < nsamuser_eff-1 && DNA_matr[(long long)n_site*(unsigned long)ns+xx] > '4') ns++;
      //     DNA_matr[(unsigned long)n_site*(long long)(nsamuser_eff)+xx] = DNA_matr[(unsigned long)n_site*(unsigned long)ns+xx];
      // }
      /*strncpy(DNA_matr+(unsigned long)n_site*(unsigned long)(nsamuser_eff),DNA_matr+(unsigned long)n_site*(unsigned long)(nsamuser_eff-1),n_site);*/
      // nsamuser_eff += 1;
    }

    //        if(wP!=0) {
    /*define the weights*/
    //            for(n_sit=0;n_sit<n_site;n_sit++) {
    //                matrix_sizepos[n_sit] =  wP[0][/*beg+*/n_sit/*-1*/];
    //                matrix_segrpos[n_sit] = wPV[0][/*beg+*/n_sit/*-1*/] /* * wV[nsit-1]*/;
    //            }
    //        }
    /*define variables for mhits*/
    *nmhits = 0;
    if ((mhitbp = (long int *)calloc(n_site, sizeof(long int))) == 0)
    {
      // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.6");
      log_fatal("Error: memory not reallocated, mhitbp get_obsstat.6");
      // for (x = 0; x < n_sam; x++)
      //   free(names[x]);
      // free(names);
      close_tfasta_file(tfasta);
      free(wP);
      free(wPV);
      free(wV);
      free(DNA_matr);
      return (0);
    }
    /*function to analyze all data*/
    if (get_obsstats(
            file_output,
            file_output_gz,
            (FILE *)0,
            // file_logerr,file_logerr_gz,
            nsamuser_eff,
            n_site,
            length_al_real,
            // names /*2*/,
            tfasta->names,
            DNA_matr /*2*/,
            wP /*matrix_sizepos*/,
            wPV /*matrix_segrpos*/,
            matrix_pol,
            matrix_freq,
            matrix_pos,
            length_al,
            length_seg,
            // args->vint_perpop_nsam,
            // args->npops,
            svratio,
            missratio,
            // args->include_unknown,
            sum_sam,
            tcga,
            matrix_sv,
            nmhits,
            // args->output,
            // args->ploidy,
            // args->outgroup_presence,
            nsites1_pop,
            nsites1_pop_outg,
            nsites2_pop,
            nsites2_pop_outg,
            nsites3_pop,
            nsites3_pop_outg,
            anx,
            bnx,
            anxo,
            bnxo,
            lengthamng,
            lengthamng_outg,
            mhitbp,
            matrix_pol_tcga,
            (long int)beg,
            args) == 0)
    {
      // for (x = 0; x < n_sam; x++)
      //   free(names[x]);
      // free(names);
      close_tfasta_file(tfasta);
      // names = 0;
      free(wP);
      free(wPV);
      free(wV);
      free(DNA_matr);
      /*free(DNA_matr2);*/
      /*          free(matrix_sizepos);
                  free(matrix_segrpos);
      */
      free(mhitbp);
      return (0);
    }
    /*free(names2);
free(DNA_matr2);*/
    /*for(x=0;x<n_sam;x++) free(names[x]); free(names);*/
    /*		free(matrix_sizepos);
        free(matrix_segrpos);
    */
    free(mhitbp);
    free(wP);
    free(wPV);
    free(wV);
  }
  return (1);
}

/*very careful!!!!!. The coordinates file is very strict in shape: name \t start \t end \n !!!!!!!!*/
int read_coordinates(
  FILE *file_wcoor, 
  SGZip *file_wcoor_gz, 
  FILE *file_output, 
  SGZip *file_output_gz,
  // FILE *file_logerr, SGZip *file_logerr_gz,
  long int **wgenes, 
  long int *nwindows, 
  char *chr_name)
{

  char *valn = 0;
  char c;
  int x;
  long int xx;
  long int dd;
  double ee;
  int inside_chr;

  /*printf("\nReading coordinates file...");*/
  fflush(stdout);
  // fprintf(file_logerr,"\nReading coordinates file...");
  log_info("Reading coordinates file...");

  if ((valn = (char *)calloc(100, sizeof(char))) == 0)
  {
    // fprintf(file_logerr,"\nError: memory not reallocated. read_coordinates.00 \n");
    log_fatal("Error: memory not reallocated, valn read_coordinates.00");
    return (0);
  }
  if ((*wgenes = (long int *)calloc(10000, sizeof(long int))) == 0)
  {
    // fprintf(file_logerr,"\nError: memory not reallocated. read_coordinates.0 \n");
    log_fatal("Error: memory not reallocated, wgenes read_coordinates.0");
    free(valn);
    return (0);
  }

  *nwindows = 0;
  c = fzgetc(file_wcoor, file_wcoor_gz);

  if (check_comment(&c, file_wcoor, file_wcoor_gz) == 0)
  {
    // fprintf(file_logerr,"\nWarning: no coordinates assigned. \n");
    log_warn("Warning: no coordinates assigned.");
    free(valn);
    return (0);
  }
  inside_chr = 0;
  /*
   if(c==0 || c==-1 || c < 1 || c=='\xff' || c=='\xfe') {
   fprintf(file_logerr,"\nWarning: no coordinates assigned. \n");
   free(*wgenes);
   file_wcoor=0;
   return(0);
   }
   else {*/
  xx = 0;
  while (c != 0 && c != -1 && c != '\xff' && c != '\xfe')
  {
    /*now keep all values: two columns, only numbers*/
    /*first column is the name of the scaffold*/
    while (c == 32 || c == 9 || c == 13 || c == 10)
      c = fzgetc(file_wcoor, file_wcoor_gz);
    if (check_comment(&c, file_wcoor, file_wcoor_gz) == 0)
    {
      c = 0;
      break;
    }
    while (c == 32 || c == 9 || c == 13 || c == 10)
      c = fzgetc(file_wcoor, file_wcoor_gz);
    x = 0;
    while (c != 32 && c != 9 && c != 13 && c != 10 && c != 0 && c != -1 && x < 100 && c != '\xff' && c != '\xfe' && c > 0)
    {
      valn[x] = c;
      c = fgetc(file_wcoor);
      x++;
    }
    valn[x] = '\0'; /*scaffold name*/
    if (c == 0 || c == -1 || c < 1 || c == '\xff' || c == '\xfe')
    {
      c = 0;
      break;
    }
    while (strcmp(chr_name, valn) != 0)
    { /*discard those rows having different scaffold than chr_name*/
      if (inside_chr == 1)
      {
        inside_chr = -1; /*we passed target region*/
        break;
      }
      while (c != 13 && c != 10 && c != 0 && c != -1 && c != '\xff' && c != '\xfe')
        c = fzgetc(file_wcoor, file_wcoor_gz);
      if (check_comment(&c, file_wcoor, file_wcoor_gz) == 0)
      {
        free(valn);
        *nwindows = (xx) / 2;
        return (1);
      }
      while (c == 32 || c == 9 || c == 13 || c == 10)
        c = fzgetc(file_wcoor, file_wcoor_gz);
      x = 0;
      while (c != 32 && c != 9 && c != 13 && c != 10 && c != 0 && c != -1 && x < 100 && c != '\xff' && c != '\xfe' && c > 0)
      {
        valn[x] = c;
        c = fgetc(file_wcoor);
        x++;
      }
      valn[x] = '\0'; /*scaffold name*/
      if (c == 0 || c == -1 || c < 1 || c == '\xff' || c == '\xfe')
      {
        c = 0;
        break;
      }
    }
    if (c == -1 || c == 0)
      break;
    if (inside_chr == -1)
      break; /*we go out*/
    inside_chr = 1;
    /*KEEP POSITIONS (first initial position, then end, next region and so on)*/
    while (c == 32 || c == 9 || c == 13 || c == 10)
      c = fzgetc(file_wcoor, file_wcoor_gz);
    if (check_comment(&c, file_wcoor, file_wcoor_gz) == 0)
    {
      c = 0;
      free(valn);
      return (0);
    }
    while (c == 32 || c == 9 || c == 13 || c == 10)
      c = fzgetc(file_wcoor, file_wcoor_gz);
    x = 0;
    while (c != 32 && c != 9 && c != 13 && c != 10 && c != 0 && c != -1 && x < 100 && c != '\xff' && c != '\xfe' && c > 0)
    {
      valn[x] = c;
      c = fzgetc(file_wcoor, file_wcoor_gz);
      x++;
    }
    valn[x] = '\0';
    wgenes[0][xx] = (long int)atof(valn);

    xx++;
    while (c == 32 || c == 9 || c == 13 || c == 10)
      c = fgetc(file_wcoor);
    if (check_comment(&c, file_wcoor, file_wcoor_gz) == 0)
    {
      c = 0;
      free(valn);
      return (0);
    }
    while (c == 32 || c == 9 || c == 13 || c == 10)
      c = fzgetc(file_wcoor, file_wcoor_gz);
    x = 0;
    while (c != 32 && c != 9 && c != 13 && c != 10 && c != 0 && c != -1 && x < 100 && c != '\xff' && c != '\xfe' && c > 0)
    {
      valn[x] = c;
      c = fgetc(file_wcoor);
      x++;
    }
    valn[x] = '\0';
    wgenes[0][xx] = (long int)round((double)atof(valn));

    xx++;
    dd = (long int)floor((double)xx / (double)10000);
    ee = (double)xx / (double)10000;
    if (dd == ee)
    {
      if ((*wgenes = realloc(*wgenes, ((long int)(dd + 1) * (long int)10000 * (long int)sizeof(long int)))) == 0)
      {
        puts("Error: realloc error read_coordinates.1\n");
        free(valn); /*free(*wgenes);*/
        return (0);
      }
    }
  }
  if (xx == 0)
  {
    // fprintf(file_logerr,"\nError: no coordinates assigned, read_coordinates.2 \n");
    log_error("Error: no coordinates assigned, read_coordinates.2");
    /*free(*wgenes);
     file_wcoor=0;*/
    return (0);
  }
  *nwindows = (xx) / 2;
  /*}*/
  free(valn);
  return 1;
}




/**
 * Reads weights and positions from a file and populates the provided variables.
 *
 * @param wtfasta The wtfasta_file structure to store the weights and positions.
 * @param file_ws The file pointer to the weights and positions file.
 * @param file_ws_gz The SGZip structure for compressed weights and positions file.
 * @param index_w The SGZIndex structure for indexing the weights and positions file.
 * @param wP Pointer to the array to store the weights - weight for each position.
 * @param wPV Pointer to the array to store  weights for the variant at each position.
 * @param wV Pointer to the array to store weight at variant (effect sizes) not yet functional although we can recover.
 * @param wlimit_end Pointer to the variable to store the end limit of weights.
 * @param init_site The initial site value.
 * @param window_size Pointer to the variable to store the window size.
 * @param n_sitesw Pointer to the variable to store the number of sites.
 * @param weight_window The weight window value.
 * @param chr_name Pointer to the character array to store the chromosome name.
 * @param first The first value.
 * @param length The length value.
 * @return Returns 1 on success, 0 on failure.
 */
int read_weights_positions_file(
		wtfasta_file *wtfasta,
		//FILE *file_ws,
		//SGZip *file_ws_gz,
		//struct SGZIndex *index_w,
		// FILE *file_logerr,SGZip *file_logerr_gz,
		double **wP,
		double **wPV,
		double **wV,
		long int *wlimit_end,
		long int init_site,
		double *window_size,
		long int *n_sitesw,
		int weight_window,
		char *chr_name,
		int first, /* Not used */
		long int length /* Not used */)
{

	double curr_window = 0;
	char *valn = 0;
	long int dd;
	long int xx;
	double ee;
	int x;
	int y;

	static char c[1];
	static char line[MSP_MAX_COL_LINE];
	static char line2[MSP_MAX_COL_LINE];
	static int col = 0;
	static int col2 = 0;
	static long int position = 0;
	static char *ID = 0;
	static int count0s, nchstr;
	static long int row_num = -1; /* fzseek: Set to -1 if you want to search by ID. */
																/*Or set NULL to the ID if you want to seach by position..*/

	/*read weights file*/
	/*keep wP, wPV, wV (not used yet) for all positions, do not need to keep positions: all are correlative*/
	// THIS should be done in the main function, not here
	if (wtfasta == NULL)
	{
		return (0);
	}
	// construct and char string with chr_name and init_site and end_site
  // construct a region string in the format chr_name:init_site-end_site
  // allocate memory for the region string dynamically
  // 25 is the length of the string ":1234567890-1234567890"
  // leading to a total of 12 chars for the number which is enough for the length of a chromosome for up to 10^12 == 1 trillion bases
  // 999999999999-999999999999
  char *region = (char *)malloc(strlen(chr_name) + 25);

	long int end_site = init_site + *window_size - 1;

  sprintf(region, "%s:%ld-%ld", chr_name, init_site, end_site);
	 hts_itr_t *iter = tbx_itr_querys(wtfasta->tbx, region);
  if (iter == NULL)
  {
    // it is possible that the region is not found in the index
		log_error("Failed to parse region: %s", chr_name);
    return 0;
  }
	int expected_sites = end_site - init_site + 1;
	
	// for wP, wPV and wV free the memory if it is already allocated
	if (*wP != NULL)
	{
		free(*wP);
		*wP = NULL;
	}
	if (*wPV != NULL)
	{
		free(*wPV);
		*wPV = NULL;
	}
	if (*wV != NULL)
	{
		free(*wV);
		*wV = NULL;
	}
	// allocate memory for wP, wPV and wV
	*wP = (double *)malloc(expected_sites * sizeof(double));
	*wPV = (double *)malloc(expected_sites * sizeof(double));
	*wV = (double *)malloc(expected_sites * sizeof(double));
 


	kstring_t str = {0, 0, NULL};
  const char *delim = ":\t\n";
	int n_sites = 0;
	while (tbx_itr_next(wtfasta->fp, wtfasta->tbx, iter, &str) >= 0) {

		if(n_sites >= expected_sites) {
      // need to reallocate memory for the matrix, if our calculations are correct, this should not happen
      log_debug("Reallocating memory for wP, wPV and wV");
    }

		// if line start with # then it is a comment
    if (str.s[0] == '#')
    {
      // skip comments :: not going to happen Here
      // log_debug("Comment : %s", str.s);
      continue;
    }
		char *cc = strtok(str.s, delim);
    int col = 0;

		// by default wV[0][n_sites] = 1.0
		wV[0][n_sites] = 1.0;

    while (cc != NULL)
    {
			// col 0 is the name of the sequence
      if(col == 0) {

      }
			// col 1 is the position
      if(col == 1) {
        position = atol(cc); // need to check on position make sure we are going in a sequential order
      }
			if(col == 2) { // wP
				wP[0][n_sites] = (double)atof(cc);
			}
			if(col == 3) { // wPV
				wPV[0][n_sites] = (double)atof(cc);
			}
			// 
			if(col == 4) { // optional wV
				wV[0][n_sites] = (double)atof(cc); 			
			}
			// increment the column
      col++;
      // get the next token
      cc = strtok(NULL, delim);
			
		}

		if (weight_window)
			curr_window += wP[0][xx] * wV[0][xx];
		else
			curr_window += 1.0;
		// increment the site
    n_sites += 1;
    
	}
  // if *n_site == 0, then the region is not found in the file
  if(n_sites == 0) {
   	*wlimit_end = init_site;
		*window_size = 0;
    return 1;
  }

	*wlimit_end = init_site + n_sites;
	*window_size = curr_window;
	*n_sitesw = n_sites;

	return (1);

	// if ((valn = (char *)calloc(100, sizeof(char))) == 0)
	// {
	// 	// fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.34 \n");
	// 	log_fatal("Error: memory not reallocated, valn get_obsdata.34");
	// 	return (0);
	// }
	// if (ID == 0)
	// {
	// 	if ((ID = (char *)calloc(100, sizeof(char))) == 0)
	// 	{
	// 		// fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.34 \n");
	// 		log_fatal("Error: memory not reallocated, ID get_obsdata.34");
	// 		free(valn);
	// 		return (0);
	// 	}
	// }
	// if (position == 0)
	// { /*assign the value of count0s and nchrstr*/
	// 	/*pass scaffold name:*/
	// 	*c = fzgetc(file_ws, file_ws_gz);
	// 	if (check_comment(c, file_ws, file_ws_gz) == 0)
	// 	{
	// 		free(valn);
	// 		free(ID);
	// 		return (0);
	// 	}
	// 	if (!(*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe'))
	// 		*c = fzgetc(file_ws, file_ws_gz);

	// 	while (*c == 10 || *c == 13 || *c == 9 || *c == 32)
	// 		*c = fzgetc(file_ws, file_ws_gz);
	// 	col = 0;
	// 	while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != 58 && *c != '\xff' && *c != '\xfe')
	// 	{
	// 		line[col] = *c;
	// 		col++;
	// 		*c = fzgetc(file_ws, file_ws_gz);
	// 	}
	// 	if (*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe')
	// 	{
	// 		// fprintf(file_logerr,"\nError: no weights assigned for %s scaffold \n",chr_name);
	// 		log_error("Error: no weights assigned for %s scaffold", chr_name);
	// 		free(valn);
	// 		return (0);
	// 	}
	// 	line[col] = '\0';
	// 	/*pass position:*/
	// 	if (*c == 58)
	// 	{
	// 		*c = fzgetc(file_ws, file_ws_gz);
	// 		y = 0;
	// 		count0s = 0;
	// 		col2 = 0;
	// 		while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != '\xff' && *c != '\xfe')
	// 		{
	// 			/*count 0s at left*/
	// 			if (*c == '0' && y == 0)
	// 				count0s += 1;
	// 			else
	// 				y += 1;
	// 			/** in case count0s is > 0, count0s+y is the total number of characters
	// 			 in the string to search for positions (add zeros at left) */

	// 			line2[col2] = *c;
	// 			col2++;
	// 			*c = fzgetc(file_ws, file_ws_gz);
	// 		}
	// 		nchstr = y + count0s; /*the number of characters in the string*/

	// 		if (*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe')
	// 		{
	// 			// fprintf(file_logerr,"\nError: no weights assigned for %s scaffold \n",chr_name);
	// 			log_error("Error: no weights assigned for %s scaffold", chr_name);
	// 			free(valn);
	// 			return (0);
	// 		}
	// 		line2[col2] = '\0';
	// 		position = atol(line2);
	// 	}
	// 	else
	// 	{
	// 		// fprintf(file_logerr,"\nError: no position assigned for %s scaffold at position %ld \nHave you included additional comments?\n",chr_name,init_site);
	// 		log_error("Error: no position assigned for %s scaffold at position %ld \nHave you included additional comments?", chr_name, init_site);
	// 		free(valn);
	// 		return (0);
	// 	}
	// }

	// /*Search for the value beg into ID: it is transformed to a string with the scaffold name */
	// /* and perhaps zeros at left if count0s > 0*/
	// if (transform_beg_chr(ID, chr_name, init_site, nchstr, count0s) != 1)
	// {
	// 	// fprintf(file_logerr,"Error transforming beg into string.\n");
	// 	log_error("Error transforming beg into string.");
	// 	free(valn);
	// 	return (0);
	// }

	// row_num = -1;
	// if (fzseekNearest(file_ws, file_ws_gz, index_w, ID, MAXLEN, &row_num) != GZ_OK)
	// { //==GZ_ERROR_DATA_FILE?
	// 	if (init_site <= length)
	// 		// fprintf(file_logerr,"ID not found in the weights file (or not indexed position): %s\n",ID);
	// 		log_error("ID not found in the weights file (or not indexed position): %s", ID);
	// 	/*no position found. Assume the file window is finished*/
	// 	free(valn);
	// 	*wlimit_end = init_site;
	// 	*window_size = 0;
	// 	return (1);
	// }

	// /*get scaffold name: perhaps is not chr_name*/
	// *c = fzgetc(file_ws, file_ws_gz);
	// while (*c == 10 || *c == 13 || *c == 9 || *c == 32)
	// 	*c = fzgetc(file_ws, file_ws_gz);
	// col = 0;
	// while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != 58 && *c != '\xff' && *c != '\xfe')
	// {
	// 	line[col] = *c;
	// 	col++;
	// 	*c = fzgetc(file_ws, file_ws_gz);
	// }
	// if (*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe')
	// {
	// 	// fprintf(file_logerr,"\nError: no weights assigned for %s scaffold \n",chr_name);
	// 	log_error("Error: no weights assigned for %s scaffold", chr_name);
	// 	free(valn);
	// 	return (0);
	// }
	// line[col] = '\0';
	// /*get position: perhaps is not init_site*/
	// if (*c == 58)
	// {
	// 	*c = fzgetc(file_ws, file_ws_gz);
	// 	y = 0;
	// 	count0s = 0;
	// 	col2 = 0;
	// 	while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != '\xff' && *c != '\xfe')
	// 	{
	// 		/*count 0s at left*/
	// 		if (*c == '0' && y == 0)
	// 			count0s += 1;
	// 		else
	// 			y += 1;
	// 		/** in case count0s is > 0, count0s+y is the total number of characters
	// 		 in the string to search for positions (add zeros at left) */

	// 		line2[col2] = *c;
	// 		col2++;
	// 		*c = fzgetc(file_ws, file_ws_gz);
	// 	}
	// 	nchstr = y + count0s; /*the number of characters in the string*/

	// 	if (*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe')
	// 	{
	// 		// fprintf(file_logerr,"\nError: no weights assigned for %s scaffold \n",chr_name);
	// 		log_error("Error: no weights assigned for %s scaffold", chr_name);
	// 		free(valn);
	// 		return (0);
	// 	}
	// 	line2[col2] = '\0';
	// 	position = atol(line2);
	// }
	// else
	// {
	// 	// fprintf(file_logerr,"\nError: no position assigned for %s scaffold at position %ld \nHave you included additional comments?\n",chr_name,init_site);
	// 	log_error("Error: no position assigned for %s scaffold at position %ld \nHave you included additional comments?", chr_name, init_site);
	// 	free(valn);
	// 	return (0);
	// }

	// if (check_comment(c, file_ws, file_ws_gz) == 0)
	// {
	// 	free(valn);
	// 	return (0);
	// }
	// if (!(*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe'))
	// 	*c = fzgetc(file_ws, file_ws_gz);

	// xx = 0; /*xx=0*/ /* VERY IMPORTANT: SITES FOR WEIGHTS START FROM 0 AND NOT FROM 1 !*/
	// curr_window = 0;
	// col = 0;
	// while (strcmp(line, chr_name) == 0 && curr_window < window_size[0])
	// {
	// 	/*Weight Position*/
	// 	while (*c == 32 || *c == 9) /* || *c == 13 || *c == 10)*/
	// 		*c = fzgetc(file_ws, file_ws_gz);
	// 	if (*c == 0 || *c == -1)
	// 		break;
	// 	x = 0;
	// 	while (*c != 32 && *c != 9 && *c != 13 && *c != 10 && *c != 0 && *c != -1 && x < 100)
	// 	{
	// 		valn[x] = *c;
	// 		*c = fzgetc(file_ws, file_ws_gz);
	// 		x++;
	// 	}
	// 	valn[x] = '\0';
	// 	wP[0][xx] = (double)atof(valn);

	// 	/*Weight Variant*/
	// 	while (*c == 32 || *c == 9) /* || *c == 13 || *c == 10)*/
	// 		*c = fzgetc(file_ws, file_ws_gz);
	// 	if (*c == 0 || *c == -1)
	// 		break;
	// 	x = 0;
	// 	while (*c != 32 && *c != 9 && *c != 13 && *c != 10 && *c != 0 && *c != -1 && x < 100)
	// 	{
	// 		valn[x] = *c;
	// 		*c = fzgetc(file_ws, file_ws_gz);
	// 		x++;
	// 	}
	// 	valn[x] = '\0';
	// 	wPV[0][xx] = (double)atof(valn);

	// 	while (*c == 32 || *c == 9)
	// 		*c = fzgetc(file_ws, file_ws_gz);
	// 	if (!(*c == 13 || *c == 10 || *c == 0 || *c == -1 || *c == '#'))
	// 	{
	// 		/*Effect size*/
	// 		while (*c == 32 || *c == 9) /* || *c == 13 || *c == 10)*/
	// 			*c = fzgetc(file_ws, file_ws_gz);
	// 		if (*c == 0 || *c == -1)
	// 			break;
	// 		x = 0;
	// 		while (*c != 32 && *c != 9 && *c != 13 && *c != 10 && *c != 0 && *c != -1 && x < 100)
	// 		{
	// 			valn[x] = *c;
	// 			*c = fzgetc(file_ws, file_ws_gz);
	// 			x++;
	// 		}
	// 		valn[x] = '\0';
	// 		wV[0][xx] = (double)atof(valn);
	// 	}
	// 	else
	// 	{
	// 		wV[0][xx] = 1.0; /*if undefined the value is 1.0 for all*/
	// 	}

	// 	if (check_comment(c, file_ws, file_ws_gz) == 0)
	// 	{
	// 		free(valn);
	// 		break;
	// 	}

	// 	/*sum window sizes and positions*/
	// 	if (weight_window)
	// 		curr_window += wP[0][xx] * wV[0][xx];
	// 	else
	// 		curr_window += 1.0;
	// 	xx++;

	// 	/*collect a new position (chr:pos)*/
	// 	if (*c < 0 || *c == 10 || *c == 13 || *c == -1 || *c == 0)
	// 	{
	// 		if (check_comment(c, file_ws, file_ws_gz) == 0)
	// 		{
	// 			free(valn);
	// 			*wlimit_end = init_site + xx /* - 1*/;
	// 			*window_size = curr_window;
	// 			*n_sitesw = xx - 1;
	// 			return (1);
	// 		}
	// 		col = 0;
	// 		while (*c == 10 || *c == 13 || *c == 9 || *c == 32)
	// 		{
	// 			*c = fzgetc(file_ws, file_ws_gz);
	// 			col++;
	// 			if (col >= MSP_MAX_COL_LINE)
	// 			{
	// 				// fprintf(file_logerr,"Register too large: position after %s:%ld, size %d\n ",chr_name,position,col);
	// 				log_error("Register too large: position after %s:%ld, size %d", chr_name, position, col);
	// 				free(valn);
	// 				*wlimit_end = init_site + xx - 1;
	// 				*window_size = curr_window;
	// 				*n_sitesw = xx - 1;
	// 				return (1);
	// 			}
	// 		}
	// 		while (check_comment(c, file_ws, file_ws_gz) == 0)
	// 		{
	// 			free(valn);
	// 			*wlimit_end = init_site + xx - 1;
	// 			*window_size = curr_window;
	// 			*n_sitesw = xx - 1;
	// 			return (1);
	// 		}
	// 		col = 0;
	// 		while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != '\xff' && *c != '\xfe' && *c != 58 && *c > 0 && col < MSP_MAX_COL_LINE - 1)
	// 		{
	// 			line[col] = *c;
	// 			col++;
	// 			*c = fzgetc(file_ws, file_ws_gz);
	// 		}
	// 		if (col >= MSP_MAX_COL_LINE)
	// 		{
	// 			// fprintf(file_logerr,"Register too large: position after %s:%ld, size %d \n",chr_name,position,col);
	// 			log_error("Register too large: position after %s:%ld, size %d", chr_name, position, col);
	// 			free(valn);
	// 			*wlimit_end = init_site + xx - 1;
	// 			*window_size = curr_window;
	// 			*n_sitesw = xx - 1;
	// 			return (1);
	// 		}
	// 		line[col] = '\0';
	// 		if (check_comment(c, file_ws, file_ws_gz) == 0)
	// 		{
	// 			free(valn);
	// 			return (0);
	// 		}
	// 		col = 0;
	// 		if (*c == 58)
	// 		{
	// 			while (*c == 10 || *c == 13 || *c == 58)
	// 				*c = fzgetc(file_ws, file_ws_gz);
	// 			col2 = 0;
	// 			while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != '\xff' && *c != '\xfe' && *c > 0 && col < MSP_MAX_COL_LINE - 1)
	// 			{
	// 				line2[col2] = *c;
	// 				col2++;
	// 				*c = fzgetc(file_ws, file_ws_gz);
	// 			}
	// 			if (col2 >= MSP_MAX_COL_LINE)
	// 			{
	// 				// fprintf(file_logerr,"Register too large: position after %s:%ld, size %d\n ",chr_name,position,col2);
	// 				log_error("Register too large: position after %s:%ld, size %d", chr_name, position, col2);
	// 				free(valn);
	// 				*wlimit_end = init_site + xx - 1;
	// 				*window_size = curr_window;
	// 				*n_sitesw = xx - 1;
	// 				return (1);
	// 			}
	// 			col2 = 0;
	// 			position = atol(line2);
	// 			if (check_comment(c, file_ws, file_ws_gz) == 0)
	// 			{
	// 				free(valn);
	// 				*wlimit_end = init_site + xx - 1 - 1;
	// 				*window_size = curr_window;
	// 				*n_sitesw = xx - 1;
	// 				return (1);
	// 			}
	// 		}
	// 	}
	// 	dd = (long int)floor((double)xx / (double)1000);
	// 	ee = (double)xx / (double)1000;
	// 	if (dd == ee)
	// 	{
	// 		if ((*wP = realloc(wP[0], ((long int)(dd + 1) * (long int)1000 * sizeof(double)))) == 0)
	// 		{
	// 			file_ws = 0;
	// 			*wP = 0;
	// 			*wPV = 0;
	// 			*wV = 0;
	// 			free(*wP);
	// 			free(*wPV);
	// 			free(*wV);
	// 			free(valn);
	// 			// fprintf(file_logerr,"Error: realloc error get_obsdata.11\n");
	// 			log_fatal("Error: realloc error , wP get_obsdata.11");
	// 			return (0);
	// 		}
	// 		if ((*wPV = realloc(wPV[0], ((long int)(dd + 1) * (long int)1000 * sizeof(double)))) == 0)
	// 		{
	// 			file_ws = 0;
	// 			*wP = 0;
	// 			*wPV = 0;
	// 			*wV = 0;
	// 			free(*wP);
	// 			free(*wPV);
	// 			free(*wV);
	// 			free(valn);
	// 			// fprintf(file_logerr,"Error: realloc error get_obsdata.11\n");
	// 			log_fatal("Error: realloc error , wPV get_obsdata.11");
	// 			return (0);
	// 		}
	// 		if ((*wV = realloc(wV[0], ((long int)(dd + 1) * (long int)1000 * sizeof(double)))) == 0)
	// 		{
	// 			file_ws = 0;
	// 			*wP = 0;
	// 			*wPV = 0;
	// 			*wV = 0;
	// 			free(*wP);
	// 			free(*wPV);
	// 			free(*wV);
	// 			free(valn);
	// 			// fprintf(file_logerr,"Error: realloc error get_obsdata.11\n");
	// 			log_fatal("Error: realloc error , wV get_obsdata.11");
	// 			return (0);
	// 		}
	// 	}
	// } // end while (strcmp(line, chr_name) == 0 && curr_window < window_size[0]

	// free(valn);
	// *wlimit_end = init_site + xx;
	// *window_size = curr_window;
	// *n_sitesw = xx;

	// return (1);
}



// int read_weights_positions_file_2(
// 		wtfasta_file *wtfasta,
// 		FILE *file_ws,
// 		SGZip *file_ws_gz,
// 		struct SGZIndex *index_w,
// 		// FILE *file_logerr,SGZip *file_logerr_gz,
// 		double **wP,
// 		double **wPV,
// 		double **wV,
// 		long int *wlimit_end,
// 		long int init_site,
// 		double *window_size,
// 		long int *n_sitesw,
// 		int weight_window,
// 		char *chr_name,
// 		int first,
// 		long int length)
// {

// 	double curr_window = 0;
// 	char *valn = 0;
// 	long int dd;
// 	long int xx;
// 	double ee;
// 	int x;
// 	int y;

// 	static char c[1];
// 	static char line[MSP_MAX_COL_LINE];
// 	static char line2[MSP_MAX_COL_LINE];
// 	static int col = 0;
// 	static int col2 = 0;
// 	static long int position = 0;
// 	static char *ID = 0;
// 	static int count0s, nchstr;
// 	static long int row_num = -1; /* fzseek: Set to -1 if you want to search by ID. */
// 																/*Or set NULL to the ID if you want to seach by position..*/

// 	/*read weights file*/
// 	/*keep wP, wPV, wV (not used yet) for all positions, do not need to keep positions: all are correlative*/
// 	if (file_ws == 0)
// 	{
// 		return (0);
// 	}
// 	if ((valn = (char *)calloc(100, sizeof(char))) == 0)
// 	{
// 		// fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.34 \n");
// 		log_fatal("Error: memory not reallocated, valn get_obsdata.34");
// 		return (0);
// 	}
// 	if (ID == 0)
// 	{
// 		if ((ID = (char *)calloc(100, sizeof(char))) == 0)
// 		{
// 			// fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.34 \n");
// 			log_fatal("Error: memory not reallocated, ID get_obsdata.34");
// 			free(valn);
// 			return (0);
// 		}
// 	}
// 	if (position == 0)
// 	{ /*assign the value of count0s and nchrstr*/
// 		/*pass scaffold name:*/
// 		*c = fzgetc(file_ws, file_ws_gz);
// 		if (check_comment(c, file_ws, file_ws_gz) == 0)
// 		{
// 			free(valn);
// 			free(ID);
// 			return (0);
// 		}
// 		if (!(*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe'))
// 			*c = fzgetc(file_ws, file_ws_gz);

// 		while (*c == 10 || *c == 13 || *c == 9 || *c == 32)
// 			*c = fzgetc(file_ws, file_ws_gz);
// 		col = 0;
// 		while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != 58 && *c != '\xff' && *c != '\xfe')
// 		{
// 			line[col] = *c;
// 			col++;
// 			*c = fzgetc(file_ws, file_ws_gz);
// 		}
// 		if (*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe')
// 		{
// 			// fprintf(file_logerr,"\nError: no weights assigned for %s scaffold \n",chr_name);
// 			log_error("Error: no weights assigned for %s scaffold", chr_name);
// 			free(valn);
// 			return (0);
// 		}
// 		line[col] = '\0';
// 		/*pass position:*/
// 		if (*c == 58)
// 		{
// 			*c = fzgetc(file_ws, file_ws_gz);
// 			y = 0;
// 			count0s = 0;
// 			col2 = 0;
// 			while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != '\xff' && *c != '\xfe')
// 			{
// 				/*count 0s at left*/
// 				if (*c == '0' && y == 0)
// 					count0s += 1;
// 				else
// 					y += 1;
// 				/** in case count0s is > 0, count0s+y is the total number of characters
// 				 in the string to search for positions (add zeros at left) */

// 				line2[col2] = *c;
// 				col2++;
// 				*c = fzgetc(file_ws, file_ws_gz);
// 			}
// 			nchstr = y + count0s; /*the number of characters in the string*/

// 			if (*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe')
// 			{
// 				// fprintf(file_logerr,"\nError: no weights assigned for %s scaffold \n",chr_name);
// 				log_error("Error: no weights assigned for %s scaffold", chr_name);
// 				free(valn);
// 				return (0);
// 			}
// 			line2[col2] = '\0';
// 			position = atol(line2);
// 		}
// 		else
// 		{
// 			// fprintf(file_logerr,"\nError: no position assigned for %s scaffold at position %ld \nHave you included additional comments?\n",chr_name,init_site);
// 			log_error("Error: no position assigned for %s scaffold at position %ld \nHave you included additional comments?", chr_name, init_site);
// 			free(valn);
// 			return (0);
// 		}
// 	}

// 	/*Search for the value beg into ID: it is transformed to a string with the scaffold name */
// 	/* and perhaps zeros at left if count0s > 0*/
// 	if (transform_beg_chr(ID, chr_name, init_site, nchstr, count0s) != 1)
// 	{
// 		// fprintf(file_logerr,"Error transforming beg into string.\n");
// 		log_error("Error transforming beg into string.");
// 		free(valn);
// 		return (0);
// 	}

// 	row_num = -1;
// 	if (fzseekNearest(file_ws, file_ws_gz, index_w, ID, MAXLEN, &row_num) != GZ_OK)
// 	{ //==GZ_ERROR_DATA_FILE?
// 		if (init_site <= length)
// 			// fprintf(file_logerr,"ID not found in the weights file (or not indexed position): %s\n",ID);
// 			log_error("ID not found in the weights file (or not indexed position): %s", ID);
// 		/*no position found. Assume the file window is finished*/
// 		free(valn);
// 		*wlimit_end = init_site;
// 		*window_size = 0;
// 		return (1);
// 	}

// 	/*get scaffold name: perhaps is not chr_name*/
// 	*c = fzgetc(file_ws, file_ws_gz);
// 	while (*c == 10 || *c == 13 || *c == 9 || *c == 32)
// 		*c = fzgetc(file_ws, file_ws_gz);
// 	col = 0;
// 	while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != 58 && *c != '\xff' && *c != '\xfe')
// 	{
// 		line[col] = *c;
// 		col++;
// 		*c = fzgetc(file_ws, file_ws_gz);
// 	}
// 	if (*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe')
// 	{
// 		// fprintf(file_logerr,"\nError: no weights assigned for %s scaffold \n",chr_name);
// 		log_error("Error: no weights assigned for %s scaffold", chr_name);
// 		free(valn);
// 		return (0);
// 	}
// 	line[col] = '\0';
// 	/*get position: perhaps is not init_site*/
// 	if (*c == 58)
// 	{
// 		*c = fzgetc(file_ws, file_ws_gz);
// 		y = 0;
// 		count0s = 0;
// 		col2 = 0;
// 		while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != '\xff' && *c != '\xfe')
// 		{
// 			/*count 0s at left*/
// 			if (*c == '0' && y == 0)
// 				count0s += 1;
// 			else
// 				y += 1;
// 			/** in case count0s is > 0, count0s+y is the total number of characters
// 			 in the string to search for positions (add zeros at left) */

// 			line2[col2] = *c;
// 			col2++;
// 			*c = fzgetc(file_ws, file_ws_gz);
// 		}
// 		nchstr = y + count0s; /*the number of characters in the string*/

// 		if (*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe')
// 		{
// 			// fprintf(file_logerr,"\nError: no weights assigned for %s scaffold \n",chr_name);
// 			log_error("Error: no weights assigned for %s scaffold", chr_name);
// 			free(valn);
// 			return (0);
// 		}
// 		line2[col2] = '\0';
// 		position = atol(line2);
// 	}
// 	else
// 	{
// 		// fprintf(file_logerr,"\nError: no position assigned for %s scaffold at position %ld \nHave you included additional comments?\n",chr_name,init_site);
// 		log_error("Error: no position assigned for %s scaffold at position %ld \nHave you included additional comments?", chr_name, init_site);
// 		free(valn);
// 		return (0);
// 	}

// 	if (check_comment(c, file_ws, file_ws_gz) == 0)
// 	{
// 		free(valn);
// 		return (0);
// 	}
// 	if (!(*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe'))
// 		*c = fzgetc(file_ws, file_ws_gz);

// 	xx = 0; /*xx=0*/ /* VERY IMPORTANT: SITES FOR WEIGHTS START FROM 0 AND NOT FROM 1 !*/
// 	curr_window = 0;
// 	col = 0;
// 	while (strcmp(line, chr_name) == 0 && curr_window < window_size[0])
// 	{
// 		/*Weight Position*/
// 		while (*c == 32 || *c == 9) /* || *c == 13 || *c == 10)*/
// 			*c = fzgetc(file_ws, file_ws_gz);
// 		if (*c == 0 || *c == -1)
// 			break;
// 		x = 0;
// 		while (*c != 32 && *c != 9 && *c != 13 && *c != 10 && *c != 0 && *c != -1 && x < 100)
// 		{
// 			valn[x] = *c;
// 			*c = fzgetc(file_ws, file_ws_gz);
// 			x++;
// 		}
// 		valn[x] = '\0';
// 		wP[0][xx] = (double)atof(valn);

// 		/*Weight Variant*/
// 		while (*c == 32 || *c == 9) /* || *c == 13 || *c == 10)*/
// 			*c = fzgetc(file_ws, file_ws_gz);
// 		if (*c == 0 || *c == -1)
// 			break;
// 		x = 0;
// 		while (*c != 32 && *c != 9 && *c != 13 && *c != 10 && *c != 0 && *c != -1 && x < 100)
// 		{
// 			valn[x] = *c;
// 			*c = fzgetc(file_ws, file_ws_gz);
// 			x++;
// 		}
// 		valn[x] = '\0';
// 		wPV[0][xx] = (double)atof(valn);

// 		while (*c == 32 || *c == 9)
// 			*c = fzgetc(file_ws, file_ws_gz);
// 		if (!(*c == 13 || *c == 10 || *c == 0 || *c == -1 || *c == '#'))
// 		{
// 			/*Effect size*/
// 			while (*c == 32 || *c == 9) /* || *c == 13 || *c == 10)*/
// 				*c = fzgetc(file_ws, file_ws_gz);
// 			if (*c == 0 || *c == -1)
// 				break;
// 			x = 0;
// 			while (*c != 32 && *c != 9 && *c != 13 && *c != 10 && *c != 0 && *c != -1 && x < 100)
// 			{
// 				valn[x] = *c;
// 				*c = fzgetc(file_ws, file_ws_gz);
// 				x++;
// 			}
// 			valn[x] = '\0';
// 			wV[0][xx] = (double)atof(valn);
// 		}
// 		else
// 		{
// 			wV[0][xx] = 1.0; /*if undefined the value is 1.0 for all*/
// 		}

// 		if (check_comment(c, file_ws, file_ws_gz) == 0)
// 		{
// 			free(valn);
// 			break;
// 		}

// 		/*sum window sizes and positions*/
// 		if (weight_window)
// 			curr_window += wP[0][xx] * wV[0][xx];
// 		else
// 			curr_window += 1.0;
// 		xx++;

// 		/*collect a new position (chr:pos)*/
// 		if (*c < 0 || *c == 10 || *c == 13 || *c == -1 || *c == 0)
// 		{
// 			if (check_comment(c, file_ws, file_ws_gz) == 0)
// 			{
// 				free(valn);
// 				*wlimit_end = init_site + xx /* - 1*/;
// 				*window_size = curr_window;
// 				*n_sitesw = xx - 1;
// 				return (1);
// 			}
// 			col = 0;
// 			while (*c == 10 || *c == 13 || *c == 9 || *c == 32)
// 			{
// 				*c = fzgetc(file_ws, file_ws_gz);
// 				col++;
// 				if (col >= MSP_MAX_COL_LINE)
// 				{
// 					// fprintf(file_logerr,"Register too large: position after %s:%ld, size %d\n ",chr_name,position,col);
// 					log_error("Register too large: position after %s:%ld, size %d", chr_name, position, col);
// 					free(valn);
// 					*wlimit_end = init_site + xx - 1;
// 					*window_size = curr_window;
// 					*n_sitesw = xx - 1;
// 					return (1);
// 				}
// 			}
// 			while (check_comment(c, file_ws, file_ws_gz) == 0)
// 			{
// 				free(valn);
// 				*wlimit_end = init_site + xx - 1;
// 				*window_size = curr_window;
// 				*n_sitesw = xx - 1;
// 				return (1);
// 			}
// 			col = 0;
// 			while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != '\xff' && *c != '\xfe' && *c != 58 && *c > 0 && col < MSP_MAX_COL_LINE - 1)
// 			{
// 				line[col] = *c;
// 				col++;
// 				*c = fzgetc(file_ws, file_ws_gz);
// 			}
// 			if (col >= MSP_MAX_COL_LINE)
// 			{
// 				// fprintf(file_logerr,"Register too large: position after %s:%ld, size %d \n",chr_name,position,col);
// 				log_error("Register too large: position after %s:%ld, size %d", chr_name, position, col);
// 				free(valn);
// 				*wlimit_end = init_site + xx - 1;
// 				*window_size = curr_window;
// 				*n_sitesw = xx - 1;
// 				return (1);
// 			}
// 			line[col] = '\0';
// 			if (check_comment(c, file_ws, file_ws_gz) == 0)
// 			{
// 				free(valn);
// 				return (0);
// 			}
// 			col = 0;
// 			if (*c == 58)
// 			{
// 				while (*c == 10 || *c == 13 || *c == 58)
// 					*c = fzgetc(file_ws, file_ws_gz);
// 				col2 = 0;
// 				while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != '\xff' && *c != '\xfe' && *c > 0 && col < MSP_MAX_COL_LINE - 1)
// 				{
// 					line2[col2] = *c;
// 					col2++;
// 					*c = fzgetc(file_ws, file_ws_gz);
// 				}
// 				if (col2 >= MSP_MAX_COL_LINE)
// 				{
// 					// fprintf(file_logerr,"Register too large: position after %s:%ld, size %d\n ",chr_name,position,col2);
// 					log_error("Register too large: position after %s:%ld, size %d", chr_name, position, col2);
// 					free(valn);
// 					*wlimit_end = init_site + xx - 1;
// 					*window_size = curr_window;
// 					*n_sitesw = xx - 1;
// 					return (1);
// 				}
// 				col2 = 0;
// 				position = atol(line2);
// 				if (check_comment(c, file_ws, file_ws_gz) == 0)
// 				{
// 					free(valn);
// 					*wlimit_end = init_site + xx - 1 - 1;
// 					*window_size = curr_window;
// 					*n_sitesw = xx - 1;
// 					return (1);
// 				}
// 			}
// 		}
// 		dd = (long int)floor((double)xx / (double)1000);
// 		ee = (double)xx / (double)1000;
// 		if (dd == ee)
// 		{
// 			if ((*wP = realloc(wP[0], ((long int)(dd + 1) * (long int)1000 * sizeof(double)))) == 0)
// 			{
// 				file_ws = 0;
// 				*wP = 0;
// 				*wPV = 0;
// 				*wV = 0;
// 				free(*wP);
// 				free(*wPV);
// 				free(*wV);
// 				free(valn);
// 				// fprintf(file_logerr,"Error: realloc error get_obsdata.11\n");
// 				log_fatal("Error: realloc error , wP get_obsdata.11");
// 				return (0);
// 			}
// 			if ((*wPV = realloc(wPV[0], ((long int)(dd + 1) * (long int)1000 * sizeof(double)))) == 0)
// 			{
// 				file_ws = 0;
// 				*wP = 0;
// 				*wPV = 0;
// 				*wV = 0;
// 				free(*wP);
// 				free(*wPV);
// 				free(*wV);
// 				free(valn);
// 				// fprintf(file_logerr,"Error: realloc error get_obsdata.11\n");
// 				log_fatal("Error: realloc error , wPV get_obsdata.11");
// 				return (0);
// 			}
// 			if ((*wV = realloc(wV[0], ((long int)(dd + 1) * (long int)1000 * sizeof(double)))) == 0)
// 			{
// 				file_ws = 0;
// 				*wP = 0;
// 				*wPV = 0;
// 				*wV = 0;
// 				free(*wP);
// 				free(*wPV);
// 				free(*wV);
// 				free(valn);
// 				// fprintf(file_logerr,"Error: realloc error get_obsdata.11\n");
// 				log_fatal("Error: realloc error , wV get_obsdata.11");
// 				return (0);
// 			}
// 		}
// 	} // end while (strcmp(line, chr_name) == 0 && curr_window < window_size[0]

// 	free(valn);
// 	*wlimit_end = init_site + xx;
// 	*window_size = curr_window;
// 	*n_sitesw = xx;

// 	return (1);
// }

// old function to read tfasta files
// n_sam initialized once when first is 0
// names initialized once when first is 0
// the actual output of this function is the DNA_matr, n_site
// return 0 if error
// return 1 if success 
// return -1 if end of file
// int function_read_tfasta_2(
//     FILE *file_input,
//     SGZip *file_input_gz,
//     struct SGZIndex *index_input,
//     // FILE *file_logerr,SGZip *file_logerr_gz,
//     long int init_site,
//     long int end_site,
//     int *n_sam,
//     long int *n_site,
//     char ***names,
//     char **DNA_matr,
//     char **matrix_pol_tcga,
//     char *chr_name,
//     int first,
//     long int length)
// {
//   static char c[1];
//   char *cc;
//   static char line[MSP_MAX_COL_LINE];
//   static char line2[MSP_MAX_COL_LINE];
//   static int col = 0;
//   static int col2 = 0;
//   static int nseq = 0;
//   int x;
//   static int maxsam = 128;
//   static long int position = 0;
//   long int end_position = end_site;
//   long int dd, count, xx;
//   double ee;
//   char *DNA_matr2;
//   static long int row_num = -1; /* fzseek: Set to -1 if you want to search by ID. */
//                                 /*Or set NULL to the ID if you want to seach by position..*/
//   static char *ID = 0;
//   /*long int f_num = 0;*/
//   int y;
//   static int count0s, nchstr;
//   /*static int chr_defined=0;*/
//   /*static int position_defined=0;*/

//   if ((DNA_matr2 = (char *)calloc(1e6, sizeof(char))) == 0)
//   {
//     // fprintf(file_logerr,"\nError: memory not reallocated. get_tfadata.23d \n");
//     log_fatal("Error: memory not reallocated, DNA_matr2 get_tfadata.23d");
//     free(DNA_matr2);
//     return (0);
//   }
//   if (ID == 0)
//   {
//     if ((ID = (char *)calloc(100, sizeof(char))) == 0)
//     {
//       // fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.34 \n");
//       log_fatal("Error: memory not reallocated, ID get_obsdata.34");
//       free(DNA_matr2);
//       return (0);
//     }
//   }

//   if (position == 0 && first == 0)
//   { /*to allow sliding windows non-overlapped*/
//     /*if names and number of samples are defined then skip the definition of names*/
//     *c = fzgetc(file_input, file_input_gz);
//     while (*c == 9 || *c == 10 || *c == 13)
//       *c = fzgetc(file_input, file_input_gz);
//     if ((*c == -1 || *c == 0 || *c == '\xff' || *c == '\xfe'))
//     {
//       // fprintf(file_logerr,"\nError: no sequence assigned for %s scaffold \n",chr_name);
//       log_error("Error: no sequence assigned for %s scaffold", chr_name);
//       free(DNA_matr2);
//       return 0;
//     }
//     line[col] = *c;
//     col++;

//     while (*c == '#')
//     {
//       /*parse comments*/
//       while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != '\xff' && *c != '\xfe')
//       {
//         *c = fzgetc(file_input, file_input_gz);
//         line[col] = *c;
//         col++;
//       }
//       line[col] = '\0';
//       if ((cc = strstr(line, "#NAMES:")) != 0)
//       {
//         /*collect names*/
//         nseq = 0;
//         cc = strtok(line, ">\n\r ");
//         while (cc != NULL)
//         {
//           if (strstr(cc, "#NAMES:") == 0)
//           {
//             strcpy(names[0][nseq], cc);
//             nseq++;
//             if (nseq == maxsam)
//             {
//               maxsam += 128;
//               if (maxsam > 32767)
//               {
//                 // fprintf(file_logerr,"\n Sorry, no more than 32767 samples are allowed.");
//                 log_error("Sorry, no more than 32767 samples are allowed.");
//                 free(DNA_matr2);
//                 return 0;
//               }
//               if ((*names = (char **)realloc(*names, maxsam * sizeof(char *))) == 0)
//               {
//                 // fprintf(file_logerr,"\nError: memory not reallocated. assigna.1 \n");
//                 log_fatal("Error: memory not reallocated, names assigna.1");
//                 free(DNA_matr2);
//                 return (0);
//               }
//               for (x = nseq; x < maxsam; x++)
//               {
//                 if ((names[0][x] = (char *)calloc(50, sizeof(char))) == 0)
//                 {
//                   // fprintf(file_logerr,"\nError: memory not reallocated. assigna.2 \n");
//                   log_fatal("Error: memory not reallocated, names assigna.2");
//                   free(DNA_matr2);
//                   return (0);
//                 }
//               }
//             }
//           }
//           cc = strtok(NULL, ">\n\r ");
//         }
//       }
//       col = 0;
//       *c = fzgetc(file_input, file_input_gz);
//       line[col] = *c;
//       col++;
//     }

//     /*assign count0s and nchrstr*/
//     /*first pass scaffold name:*/
//     col = 0;
//     while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != 58 && *c != '\xff' && *c != '\xfe')
//     {
//       line[col] = *c;
//       col++;
//       *c = fzgetc(file_input, file_input_gz);
//     }
//     if (*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe')
//     {
//       free(DNA_matr2);
//       return (-1);
//     }
//     line[col] = '\0';
//     /*second: pass position value. Count 0s*/
//     if (*c == 58)
//     {
//       *c = fzgetc(file_input, file_input_gz);
//       y = 0;
//       count0s = 0;
//       col2 = 0;
//       while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != '\xff' && *c != '\xfe')
//       {
//         // count 0s at left
//         if (*c == '0' && y == 0)
//           count0s += 1;
//         else
//           y += 1;
//         // in case count0s is > 0, count0s+y is the total number of characters
//         // in the string to search for positions (add zeros at left)

//         line2[col2] = *c;
//         col2++;
//         *c = fzgetc(file_input, file_input_gz);
//       }
//       nchstr = y + count0s; // the number of characters in the string

//       if (*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe')
//       {
//         free(DNA_matr2);
//         return (-1);
//       }
//       line2[col2] = '\0';
//       position = atol(line2);
//     }
//   }
//   n_sam[0] = nseq;
//   /*
//   if(check_comment(c,file_input,file_input_gz) == 0) {
//       free(DNA_matr2);
//       return(-1);
//   }
//   */
//   /*Search for the value beg into ID: it is transformed to a string with the scaffold name */
//   /* and perhaps zeros at left if count0s > 0*/
//   if (transform_beg_chr(ID, chr_name, init_site, nchstr, count0s) != 1)
//   {
//     // fprintf(file_logerr,"Error transforming beg into string.\n");
//     log_error("Error transforming beg into string.");
//     free(DNA_matr2);
//     return (0);
//   }
//   /*f_num = position;*/
//   row_num = -1;

//   if (fzseekNearest(file_input, file_input_gz, index_input, ID, MAXLEN, &row_num) != GZ_OK)
//   { //==GZ_ERROR_DATA_FILE?
//     if (init_site <= length)
//       // fprintf(file_logerr,"ID not found in the tfa file (or not indexed position): %s\n",ID);
//       log_error("ID not found in the tfa file (or not indexed position): %s", ID);
//     /*no position found. Assume the file window is finished*/
//     free(DNA_matr2);
//     return (1);
//   }
//   /*f_num = row_num;*/

//   /*get scaffold name*/
//   *c = fzgetc(file_input, file_input_gz);
//   while (*c == '#')
//   {
//     while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != '\xff' && *c != '\xfe')
//       *c = fzgetc(file_input, file_input_gz);
//     *c = fzgetc(file_input, file_input_gz);
//   }
//   col = 0;
//   while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != 58 && *c != '\xff' && *c != '\xfe')
//   {
//     line[col] = *c;
//     col++;
//     *c = fzgetc(file_input, file_input_gz);
//   }
//   if (*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe')
//   {
//     /*fprintf(file_logerr,"\nError: no sequence assigned for %s scaffold \n",chr_name);*/
//     free(DNA_matr2);
//     return (-1);
//   }
//   line[col] = '\0';
//   /*get position*/
//   if (*c == 58)
//   {
//     *c = fzgetc(file_input, file_input_gz);
//     y = 0;
//     count0s = 0;
//     col2 = 0;
//     while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != '\xff' && *c != '\xfe')
//     {
//       /*count 0s at left*/
//       if (*c == '0' && y == 0)
//         count0s += 1;
//       else
//         y += 1;
//       /** in case count0s is > 0, count0s+y is the total number of characters
//        in the string to search for positions (add zeros at left) */

//       line2[col2] = *c;
//       col2++;
//       *c = fzgetc(file_input, file_input_gz);
//     }
//     nchstr = y + count0s; /*the number of characters in the string*/

//     if (*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe')
//     {
//       /*fprintf(file_logerr,"\nError: no sequence assigned for %s scaffold \n",chr_name);*/
//       free(DNA_matr2);
//       return (-1);
//     }
//     line2[col2] = '\0';
//     position = atol(line2);
//   }
//   else
//   {
//     /*fprintf(file_logerr,"\nError: no position assigned for %s scaffold at row %ld \n",chr_name,row_num);*/
//     free(DNA_matr2);
//     return (-1);
//   }

//   if (check_comment(c, file_input, file_input_gz) == 0)
//   {
//     free(DNA_matr2);
//     return (-1);
//   }
//   if (!(*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe'))
//     *c = fzgetc(file_input, file_input_gz);

//   *n_site = 0;
//   col = 0;
//   count = 0;
//   while (strcmp(line, chr_name) == 0 && position <= end_position)
//   {
//     /*chr_defined = 0;*/
//     /*position_defined = 0;*/
//     switch (*c)
//     {
//     case 'T':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '1';
//       col += 1;
//       count += 1;
//       break;
//     case 't':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '1';
//       col += 1;
//       count += 1;
//       break;
//     case 'U':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '1';
//       col += 1;
//       count += 1;
//       break;
//     case 'u':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '1';
//       col += 1;
//       count += 1;
//       break;
//     case 'C':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '2';
//       col += 1;
//       count += 1;
//       break;
//     case 'c':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '2';
//       col += 1;
//       count += 1;
//       break;
//     case 'G':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '3';
//       col += 1;
//       count += 1;
//       break;
//     case 'g':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '3';
//       col += 1;
//       count += 1;
//       break;
//     case 'A':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '4';
//       col += 1;
//       count += 1;
//       break;
//     case 'a':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '4';
//       col += 1;
//       count += 1;
//       break;
//     case 0:
//       break;
//     case -1:
//       break;
//     case 10:
//       break;
//     case 13:
//       break;
//     case 32:
//       break;
//     case 9:
//       break;
//     case 'N':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '5';
//       col += 1;
//       count += 1;
//       break;
//     case 'n':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '5';
//       col += 1;
//       count += 1;
//       break;
//     case '?':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '5';
//       col += 1;
//       count += 1;
//       break;
//     case '-':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '6';
//       col += 1;
//       count += 1;
//       break;
//       /*gaps are converted to N*/
//     case 'W':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = 'w';
//       col += 1;
//       count += 1;
//       break;
//     case 'w':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = 'w';
//       col += 1;
//       count += 1;
//       break;
//     case 'M':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = 'm';
//       col += 1;
//       count += 1;
//       break;
//     case 'm':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = 'm';
//       col += 1;
//       count += 1;
//       break;
//     case 'R':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = 'r';
//       col += 1;
//       count += 1;
//       break;
//     case 'r':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = 'r';
//       col += 1;
//       count += 1;
//       break;
//     case 'Y':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = 'y';
//       col += 1;
//       count += 1;
//       break;
//     case 'y':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = 'y';
//       col += 1;
//       count += 1;
//       break;
//     case 'K':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = 'k';
//       col += 1;
//       count += 1;
//       break;
//     case 'k':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = 'k';
//       col += 1;
//       count += 1;
//       break;
//     case 'S':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = 's';
//       col += 1;
//       count += 1;
//       break;
//     case 's':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = 's';
//       col += 1;
//       count += 1;
//       break;
//     case 'b':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '5';
//       col += 1;
//       count += 1;
//       break;
//     case 'B':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '5';
//       col += 1;
//       count += 1;
//       break;
//     case 'd':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '5';
//       col += 1;
//       count += 1;
//       break;
//     case 'D':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '5';
//       col += 1;
//       count += 1;
//       break;
//     case 'h':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '5';
//       col += 1;
//       count += 1;
//       break;
//     case 'H':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '5';
//       col += 1;
//       count += 1;
//       break;
//     case 'v':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '5';
//       col += 1;
//       count += 1;
//       break;
//     case 'V':
//       DNA_matr2[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '5';
//       col += 1;
//       count += 1;
//       break;
//     default:
//       // fprintf(file_logerr,"Unexpected value in tfa file: position %ld, sample %d \n ",position,col);
//       //  fprintf(file_logerr,"%c\n",*c);
//       log_error("Unexpected value in tfa file: position %ld, sample %d \n%c", position, col, *c);
//       free(DNA_matr2);
//       return (-1);
//       break;
//     }
//     if (check_comment(c, file_input, file_input_gz) == 0)
//     {
//       break;
//     }
//     *c = fzgetc(file_input, file_input_gz);
//     if (*c < 0 || *c == 10 || *c == 13 || *c == -1 || *c == 0 || *c == '\xff' || *c == '\xfe')
//     {
//       if (check_comment(c, file_input, file_input_gz) == 0)
//       {
//         if (col == *n_sam || col == 0)
//         {
//           *n_site += 1;
//           break;
//         }
//         else
//         {
//           free(DNA_matr2);
//           return (-1);
//         }
//       }
//       while (*c == 10 || *c == 13)
//         *c = fzgetc(file_input, file_input_gz);
//       if (check_comment(c, file_input, file_input_gz) == 0)
//       {
//         break;
//       }
//       *n_site += 1;
//       /*if(position != *n_site) {
//           printf("check here");
//       }*/
//       /*read new line*/
//       col = 0;
//       while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c != '\xff' && *c != '\xfe' && *c != 58 && *c > 0 && col < MSP_MAX_COL_LINE - 1)
//       {
//         line[col] = *c;
//         col++;
//         *c = fzgetc(file_input, file_input_gz);
//       }
//       line[col] = '\0';
//       /*chr_defined = 1;*/
//       if (check_comment(c, file_input, file_input_gz) == 0)
//       {
//         break;
//       }
//       col = 0;
//       if (*c == 58)
//       {
//         while (*c == 10 || *c == 13 || *c == 58)
//           *c = fzgetc(file_input, file_input_gz);
//         col2 = 0;
//         while (*c != 0 && *c != -1 && *c != 10 && *c != 13 && *c != 9 && *c != 32 && *c > 0)
//         {
//           line2[col2] = *c;
//           col2++;
//           *c = fzgetc(file_input, file_input_gz);
//         }
//         line2[col2] = '\0';
//         col2 = 0;
//         position = atol(line2);
//         if (end_site == -1)
//           end_position = position + 1;
//         /*position_defined = 1;*/
//         if (check_comment(c, file_input, file_input_gz) == 0)
//         {
//           free(DNA_matr2);
//           return (-1);
//         }
//       }
//       /*else {
//           printf("\nError: Problem reading line %s:%ld\n",line,position);
//       }*/
//     }
//     /*realloc DNAmatr if nnecessary*/
//     dd = (long int)floor((double)count / (double)1e6);
//     ee = (double)count / (double)1e6;
//     if ((double)dd == ee)
//     {
//       if ((DNA_matr2 = (char *)realloc((char *)DNA_matr2, ((long int)dd + (long int)1) * (long int)1e6 * sizeof(char))) == 0)
//       {
//         // fprintf(file_logerr,"Error: realloc error varchar.1\n");
//         log_fatal("Error: realloc error, DNA_matr2 varchar.1");
//         return (0);
//       }
//     }
//   }
//   /**n_site += 1;*/
//   /**/
//   if (*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe')
//   {
//     if (col == *n_sam || col == 0)
//       *n_site += 1;
//     else
//       return (-1);
//   }
//   /**/
//   /*return n_site,n_sam,DNA matrix and names in pointers*/
//   /*
//   if ((*DNA_matr = (char *)realloc((char *) *DNA_matr,n_site[0]*(long long)n_sam[0]*sizeof(char))) == 0) {
//       fprintf(file_logerr,"\nError: memory not reallocated. get_tfadata.23d \n");
//       free(DNA_matr2);
//       return(0);
//   }
//   */
//   free(*DNA_matr);
//   if ((*DNA_matr = (char *)calloc(n_site[0] * (long long)n_sam[0], sizeof(char))) == 0)
//   {
//     // fprintf(file_logerr,"\nError: memory not reallocated. get_tfadata.23d2 \n");
//     log_fatal("Error: memory not reallocated, DNA_matr get_tfadata.23d2");
//     free(DNA_matr2);
//     return (0);
//   }
//   for (x = 0; x < n_sam[0]; x++)
//   {
//     for (xx = 0; xx < n_site[0]; xx++)
//     { /*transpose */
//       DNA_matr[0][(((long long)n_site[0] * (unsigned long)x) + (unsigned long)xx)] =
//           DNA_matr2[(((long long)n_sam[0] * (unsigned long)xx) + (unsigned long)x)];
//     }
//   }
//   free(DNA_matr2);

//   if (check_comment(c, file_input, file_input_gz) == 0)
//   {
//     return (-1);
//   }
//   return (1);
// }

// Deprecated function
int transform_beg_chr(char *ID, char *chr_name, long int beg, int nchstr, int count0s)
{
  /*Search for the value beg: it may be transformed to a string with zeros at left if count0s > 0*/
  /*first count the number of digits in beg and transform beg in a string*/
  long int zi = beg;
  double zr;
  int p10 /*,p10i*/;
  int x /*,nnchstr*/;
  char ID2[100];

  memset(ID2, 0, 100);

  p10 = 0;
  do
  {
    p10 += 1;
    zr = (double)zi / (double)pow((double)10, (double)p10);
  } while (zr > 1.0);

  if (count0s > 0)
  { /*transformed to a string with zeros at left if count0s > 0*/
    /*p10 is the number of numbers after the zeros*/
    for (x = nchstr - 1; x >= p10; x--)
    {
      ID2[nchstr - 1 - x] = '0';
    }
    /*nnchstr = nchstr;*/
  } /*
   else {
       nnchstr = p10;
   }
   */
  sprintf(ID2, "%s%ld", ID2, beg);

  /*
  p10i = p10-1;
  for(x=0;x<p10;x++) {
      zr = (double)zi/(double)pow((double)10,(double)p10i);
      zr = floor(zr);
      zi -=  zr*(double)pow((double)10,(double)p10i);
      ID2[nnchstr-1-p10i] = zr+48;
      p10i -= 1;
  }
  ID2[nnchstr] = '\0';
  */

  strcpy(ID, chr_name);
  strcat(ID, ":");
  strcat(ID, ID2);

  return (1);
}
int check_comment(char *c, FILE *file_input, SGZip *file_input_gz)
{
  if (*c == '#')
  {
    /*parse comments*/
    while (*c != 10 && *c != 13 && *c != 0 && *c != -1 && *c != '\xff' && *c != '\xfe')
    {
      *c = fzgetc(file_input, file_input_gz);
    }
    if (*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe')
      return (0);
  }
  if (*c == 0 || *c == -1 || *c == '\xff' || *c == '\xfe')
    return (0);

  return (1);
}

// #include <htslib/tbx.h>
// #include <htslib/kseq.h>













// // n_sam initialized once when first is 0
// // names initialized once when first is 0
// // the actual output of this function is the DNA_matr, n_site
// // return 0 if error
// // return 1 if success 
// // return -1 if end of file : since we are using an index and we are reading a specific region of the file, we should not reach the end of the file
// // we can not reach the end of the file, we should always read the region we are interested in
// // or retun n_site = 0 if the region we are interested in is not in the file
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
//     long int length)
// {
//   // n_sam : number of samples
//   // n_site : number of sites, sites correspond to positions in the DNA sequences
//   // names : names of the samples
//   // DNA_matr : matrix of DNA sequences
//   // matrix_pol_tcga :  ?? not used here ??
//   // chr_name : name of the chromosome / scaffold to get data for
//   // first : First time we call this function
//   // length : length of the chromosome / scaffold

//   // this get the data for chr_name sequence , start from init_site to end_site
//   // could end earlier if the length of the sequence is smaller than end_site

//   // return 0 for error, 1 for success but did not reach end of sequences  and -1 for end of sequence

//   // DNA_matr structure :
//   // DNA_matr is a matrix of DNA sequences,
//   // each row is a site and each column is a variant value for a sample
//   // example :
//   // Nucleotides : A, C, G, T, N, -
//   // Sample 1 : A, C, G, T, N, -
//   // ACTGAAAAAAAAAACCA
//   // ACCCCCCCCCCCCCCCC
//   // means that the first site is A for sample 1 and C for sample 2 , and T for sample 3 and G for sample 4 , etc
//   // then the second site is A for sample 1 and C for sample 2 , and C for sample 3 and C for sample 4 , etc
//   // the matrix is transposed later, so each row is a sample and each column is a site

//   // some static variables to keep track of the position in the file
//   // TODO : check if we can remove them
//   // sample increment to reallocate memory for samples names

//   static tfasta_file tfasta;
//   memset(&tfasta, 0, sizeof(tfasta_file));
//   char *tfasta_fname = file_input_gz->file_name;
  
  
//   static long int position = 1;

//   // names :: is pre allocated to 128 samples
//   // each sample name is allocated to 50 characters

//   char *c;
  
//   // // return 0 for error, 1 for success
//   // // Open the file with tabix index
//   // tbx_t *tbx = tbx_index_load(filename);
//   // if (tbx == NULL)
//   // {
//   //   // fprintf(stderr, "Failed to load index for: %s\n", filename);
//   //   log_error("Failed to load index for: %s", filename);
//   //   return 0;
//   // }
//   // htsFile *fp = hts_open(filename, "r");
//   // if (fp == NULL)
//   // {
//   //   // fprintf(stderr, "Failed to open file: %s\n", filename);
//   //   log_error("Failed to open file: %s", filename);
//   //   tbx_destroy(tbx);
//   //   return 0;
//   // }

//   // if first time we call this function
//   // TODO :: this should be moved to the main function and not here

//   init_tfasta_file(&tfasta, tfasta_fname);
  
//   // copy the number of samples and names to the output variables
//   *n_sam = tfasta.n_sam;
//   *names = tfasta.names;
//   *n_site = 0;
//   if (first == 0)
//   {
    
//     // // read and fill samples names
//     // kstring_t str = {0, 0, NULL};
//     // int len;
//     // while ((len = hts_getline(fp, KS_SEP_LINE, &str)) >= 0)
//     // {
//     //   if (str.s[0] != '#')
//     //     break; // Stop at the first non-header line
//     //   printf("%s\n", str.s);
//     //   // if line start #NAMES: then parse the names
//     //   if (strstr(str.s, "#NAMES:") != 0)
//     //   {
//     //     // collect names
//     //     int nseq = 0;
//     //     char *cc = strtok(str.s, ">\n\r ");
//     //     while (cc != NULL)
//     //     {
//     //       if (strstr(cc, "#NAMES:") == 0)
//     //       {
//     //         strcpy(names[0][nseq], cc);
//     //         nseq++;
//     //         if (nseq == maxsam)
//     //         {
//     //           maxsam += NSAM_INC;
//     //           if (maxsam > 32767)
//     //           {
//     //             // fprintf(file_logerr,"\n Sorry, no more than 32767 samples are allowed.");
//     //             log_error("Sorry, no more than 32767 samples are allowed.");
//     //             // free(DNA_matr2);
//     //             return 0;
//     //           }
//     //           if ((*names = (char **)realloc(*names, maxsam * sizeof(char *))) == 0)
//     //           {
//     //             // fprintf(file_logerr,"\nError: memory not reallocated. assigna.1 \n");
//     //             log_fatal("Error: memory not reallocated, names assigna.1");
//     //             // free(DNA_matr2);
//     //             return (0);
//     //           }
//     //           for (int x = nseq; x < maxsam; x++)
//     //           {
//     //             if ((names[0][x] = (char *)calloc(50, sizeof(char))) == 0)
//     //             {
//     //               // fprintf(file_logerr,"\nError: memory not reallocated. assigna.2 \n");
//     //               log_fatal("Error: memory not reallocated, names assigna.2");
//     //               // free(DNA_matr2);
//     //               return (0);
//     //             }
//     //           }
//     //         }
//     //       }
//     //       cc = strtok(NULL, ">\n\r ");
//     //       if (cc != NULL) printf("%s\n", cc); // for debug
//     //       log_debug("%s", cc);
//     //     }
//     //     *n_sam = nseq;
//     //   }
//     // }
//     // if (len < 0)
//     // {
//     //   // fprintf(stderr, "Failed to read the first line of the file: %s\n", filename);
//     //   log_error("Failed to read the first line of the file: %s", filename);
//     //   hts_close(fp);
//     //   tbx_destroy(tbx);
//     //   return 0;
//     // }
//   }

//   // construct and char string with chr_name and init_site and end_site
//   // construct a region string in the format chr_name:init_site-end_site
//   // allocate memory for the region string dynamically
//   // 25 is the length of the string ":1234567890-1234567890"
//   // leading to a total of 12 chars for the number which is enough for the length of a chromosome for up to 10^12 == 1 trillion bases
//   // 999999999999-999999999999
//   char *region = (char *)malloc(strlen(chr_name) + 25);

//   sprintf(region, "%s:%ld-%ld", chr_name, init_site, end_site);

//   hts_itr_t *iter = tbx_itr_querys(tfasta.tbx, region);
//   if (iter == NULL)
//   {
//     // it is possible that the region is not found in the index
//     fprintf(stderr, "Failed to parse region: %s\n", chr_name);
//     // hts_close(tfasta.fp);
//     // tbx_destroy(tfasta.tbx);
//     // reset the position to 0
    
//     position = 1;
//     return 0;
//   }

//   // allocate
//   // DNA_matr2 is a temporary matrix to store the data
//   char *DNA_matr2;
//   // allocate memory for the matrix by muliple of tfasta.n_sam by init n positions
//   // we expect that the number of sites will be in the range of init_site and end_site
//   int expected_sites = end_site - init_site + 1;
//   long long DNA_matr2_size = tfasta.n_sam * expected_sites;
//   if ((DNA_matr2 = (char *)calloc(DNA_matr2_size, sizeof(char))) == 0)
//   {
//     // fprintf(file_logerr,"\nError: memory not reallocated. get_tfadata.23d \n");
//     log_fatal("Error: memory not reallocated, DNA_matr2 get_tfadata.23d");
//     free(DNA_matr2);
//     // close the file and free memory
//     close_tfasta_file(&tfasta);
//     return (0);
//   }
  


//   // tbx_itr_querys will return an iterator to the region in the file
//   // it only return what we asked for, so we need to iterate over the iterator to get the data

//   kstring_t str = {0, 0, NULL};
//   const char *delim = ":\t\n";
//   // keep track DNA_matr2 size
//   int count = 0;
//   while (tbx_itr_next(tfasta.fp, tfasta.tbx, iter, &str) >= 0)
//   {
//     // if line start with # then it is a comment
//     if (str.s[0] == '#')
//     {
//       // skip comments :: not going to happen Here
//       // log_debug("Comment : %s", str.s);
//       continue;
//     }

//     // log_debug("TFA at site %d  : %s", *n_site ,  str.s);
//     // tokenize the line to get the data
//     char *cc = strtok(str.s, delim);
//     int col = 0;
//     while (cc != NULL)
//     {
//       // log_debug("col %d  : %s", col,  cc);

//       // col 0 is the name of the sequence
//       if(col == 0) {

//       }

//       // col 1 is the position
//       if(col == 1) {
//         position = atol(cc);
//       }
      
//       // col 2 is the nucleotides per sample 
//       if(col == 2) {
//         // TODO : fill the DNA_matr matrix

//         // DOES we require that strlen(cc) == n_sam ??

//         // iterate over the nucleotides and fill the matrix
//         for (int i = 0; i < strlen(cc); i++)
//         {
//           // log_debug("Nucleotide %d  : %c", i,  cc[i]);
//           // fill the matrix
//           // DNA_matr[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '1';
//           int DNA_matr2_index = (((long long)tfasta.n_sam * (unsigned long)*n_site) + (unsigned long)i);
//           char dna_char =  get_DNA_char(&cc[i]);
//           if(dna_char == -1) {
//             log_error("Unexpected value in tfa file: position %ld, sample %d \n%c", position, i, cc[i]);
//             free(DNA_matr2);
//             return (-1);
//           }
//           if(dna_char > 0) {
//             DNA_matr2[DNA_matr2_index] = dna_char;
//             count++;
//           }
//         }

//       } 

//       // in case we have more than 3 columns
//       // we can ignore them for now
//       if(col > 2) {
//         log_debug("Ignoring unsupported data column %d  : %s", col,  cc);
//       }

//       // increment the column
//       col++;
//       // get the next token
//       cc = strtok(NULL, delim);
//     }
//     // increment the site
//     *n_site += 1;
//     if(*n_site > expected_sites) {
//       // need to reallocate memory for the matrix, if our calculations are correct, this should not happen
//       log_debug("Reallocating memory for DNA_matr2");
//     }
//   }
//   // if *n_site == 0, then the region is not found in the file
//   if(*n_site == 0) {
   
//     return 1;
//   }

//   free(region);
//   free(str.s);

//   // possible side effect, memory was allocated for DNA_matr2 before
//   free(*DNA_matr);
//   // actual size of the matrix as n_site can be less than expected_sites
//   // but it should not be more than expected_sites
//   DNA_matr2_size = tfasta.n_sam * *n_site;
//   if ((*DNA_matr = (char *)calloc(DNA_matr2_size, sizeof(char))) == 0)
//   {
//     // fprintf(file_logerr,"\nError: memory not reallocated. get_tfadata.23d2 \n");
//     log_fatal("Error: memory not reallocated, DNA_matr get_tfadata.23d2");
//     free(DNA_matr2);
//     return (0);
//   }

//   for (unsigned long x = 0; x < tfasta.n_sam; x++) {
//     for (unsigned long xx = 0; xx < *n_site; xx++) { /*transpose */
//       (*DNA_matr)[((*n_site * x) + xx)] =
//           DNA_matr2[((tfasta.n_sam * xx) + x)];
//     }
//   }
//   free(DNA_matr2);


//   return (1);
// }