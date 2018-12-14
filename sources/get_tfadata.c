//
//  get_tfadata.c
//  xcode_project
//
//  Created by Sebastian Ramos-Onsins on 09/08/15.
//
//

#include "get_tfadata.h"

int get_tfadata(FILE *file_output,
                SGZip *file_output_gz,
                FILE *file_input,
                SGZip *input_gz,
                struct SGZIndex *index_input,
                char *file_wps,
                FILE *file_ws,
                SGZip *file_ws_gz,
                struct SGZIndex *index_w,
                FILE *file_logerr,
                SGZip *file_logerr_gz,
                char **matrix_pol,
                long int **matrix_freq,
                long int **matrix_pos,
                long int *length,
                long int *length_seg,
                double *length_al,
                long int *length_al_real,
                int mainargc,
                int *nsamuser,
                int npops,
                double *svratio,
                double *missratio,
                int include_unknown,
                double *sum_sam,
                double **tcga,
                long int **matrix_sv,
                long int *nmhits,
                int output,
                int outgroup_presence,
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
                long int first_slide,
                long int slide,/**/
                long int window,/**/
                int Physical_length,
                long int *li/**/,
                int *npriors,
                double **vector_priors,
                char **matrix_pol_tcga,
                char *chr_name,
                int first
    )
{
    char *DNA_matr2 = 0;
    char **names2   = 0; /* limit of the name of the lines to 50 characters. be careful ... */
    
/*	double *matrix_sizepos = 0; *//*size of each position, 0 not include, >0 include, Syn/Nsyn positions are allowed*/
/*    double *matrix_segrpos = 0; *//*always 1 except in those biallelic positions that are not desired (e.g., choose syn but we do not want nsyn)*/
    
    long int *mhitbp; /*vector of positions where mhits are*/
    /* long int count; */
    long int xx;
    /* int c; */
    int n_sam/*,ns*/;
    /*long int n_sit;*/
    /* int nseq; */
    /* int maxsam; */
    /* int n_excl; */
    int n_samp;
    long int n_site;
	char ploidy[2];
    
    int x;
/*    static long int maxsites = 0;*/
    
    int nsamtot;
    int nsamuser_eff;
    int flag_change_sort; /*in case the order of samples is not consecutive*/
    
    static char *DNA_matr = 0;
    static char **names   = 0; /* limit of the name of the lines to 50 characters. be careful ... */
    static char chr_name2[MSP_MAX_NAME];
	
    static long int init_site = 1;
	static long wc = 0;
	long int beg,end;
	double wl;

    /*tfasta windows and weights*/
    /*wV: weight at variant (effect sizes)*//*not yet functional although we can recover*/
    /*wP: weight for each position*/
    /*wPV: weight for the variant at each position*/
    long int wlimit_end = 0; /*initial value*/
    long int n_sitesw=0;
    
    double *wV;
    double *wP;
    double *wPV;
    double window_size;
    int weight_window;
    
    if((wP = (double *)calloc(1000,sizeof(double))) == 0) {
        fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.14 \n");
        return(0);
    }
    if((wPV = (double *)calloc(1000,sizeof(double))) == 0) {
        fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.25 \n");
        return(0);
    }
    if((wV = (double *)calloc(1000,sizeof(double))) == 0) {
        fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.23 \n");
        return(0);
    }
    
    ploidy[0]='1';
    ploidy[1]='\0';
    /* count = 0; */
    /* c = 0; */
    n_sam = 0;
    /*n_sit = 0;*/
    /* nseq  = 0; */
    /* maxsam= 128; */
    n_samp= 0;
    n_site= 0;
    /* n_excl= 0; */
    
    nsamtot = 0;
    for(x=0;x<npops;x++) {
        nsamtot += nsamuser[x];
    }    
    nsamuser_eff = (nsamtot- !outgroup_presence) ;
    
    if(names == 0) { /* only initialize once. Check */
        if((names = (char **)calloc(128,sizeof(char *))) == 0) {
            fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.1 \n");
            free(wP);free(wPV);free(wV);
            return(0);
        }
        for(x=0;x<128;x++) {
            if((names[x] = (char *)calloc(50,sizeof(char))) == 0) {
                fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.2 \n");
                free(wP);free(wPV);free(wV);
                return(0);
            }
        }
        if((DNA_matr = (char *)calloc(10000,sizeof(char))) == 0) {
            for(x=0;x<128;x++) free(names[x]); free(names);
            fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.3 \n");
            free(wP);free(wPV);free(wV);
            return(0);
        }
    }
    if(strcmp(chr_name2, chr_name) != 0) {
        init_site = 1;
        wc = 0;
        strcpy(chr_name2,chr_name);
    }
    
    if(init_site == 1)
        init_site += first_slide;
    
	/*FIND THE WINDOW OF POSITIONS TO ANALYZE*/
    if(nwindows==0) {
		if(Physical_length == 1) {
            weight_window = 0; /*physical positions*/
            window_size = (double)window;
			beg = init_site;/*the initial position of the window*/
			end = init_site + window - 1;/*the final position of the window*/
			init_site += slide;/*init for next window (static variable)*/
            if(file_ws != 0) {
                if(read_weights_positions_file(file_ws,file_ws_gz,index_w,
                                               file_logerr,file_logerr_gz,
                                               &wP,&wPV,&wV,&wlimit_end,
                                               beg,&window_size,&n_sitesw,weight_window,
                                               chr_name,first,*length) == 0) {
                    fprintf(file_logerr,"Error processing weighting file %s\n", file_wps);
                    free(wP);free(wPV);free(wV);
                    exit(1);
                }
            }
            else {
                n_sitesw = end-beg+1+1; /*the position 0 does not count: we translate 1 right*/
                if((wP = (double *)realloc((double *)wP,n_sitesw*sizeof(double))) == 0) {
                    fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.14 \n");
                    free(wP);free(wPV);free(wV);
                    return(0);
                }
                if((wPV = (double *)realloc((double *)wPV,n_sitesw*sizeof(double))) == 0) {
                    fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.25 \n");
                    free(wP);free(wPV);free(wV);
                    return(0);
                }
                if((wV = (double *)realloc((double *)wV,n_sitesw*sizeof(double))) == 0) {
                    fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.23 \n");
                    free(wP);free(wPV);free(wV);
                    return(0);
                }
                for(xx=0;xx<n_sitesw;xx++) {
                    wP[xx] = wPV[xx] = wV[xx] = (float)1;
                }
            }
		}
		else {/*if(Physical_length == 0), that is weighting windows*/
            weight_window = 1; /*weight positions*/
            window_size = (double)window;
			beg = init_site;/*the initial position of the window*/
            if(file_ws == 0) {
                fprintf(file_logerr,"Error: no weights file\n");
                free(wP);free(wPV);free(wV);
                exit(1);
            }
            if(read_weights_positions_file(file_ws,file_ws_gz,index_w,
                                           file_logerr,file_logerr_gz,
                                           &wP,&wPV,&wV,&wlimit_end,
                                           beg,&window_size,&n_sitesw,weight_window,
                                           chr_name,first,*length) == 0) {
                fprintf(file_logerr,"Error processing weights file\n");
                free(wP);free(wPV);free(wV);
                exit(1);
            }
            end = wlimit_end;/*the final position of the window (static variable)*/
            
            if(slide < window) {
                wl = 0.0;
                xx = 0;
                while (wl < slide && xx < n_sitesw) {
                    wl += wP[xx] * wPV[xx];
                    xx++;
                }
                init_site = beg + xx;
            }
            else {
                if(slide == window) {
                    init_site = end + 1;
                }
                else {
                    init_site = end + 1;
                    window_size = slide - window;
                    if(read_weights_positions_file(file_ws,file_ws_gz,index_w,
                                                   file_logerr,file_logerr_gz,
                                                   &wP,&wPV,&wV,&wlimit_end,
                                                   init_site,&window_size,&n_sitesw,weight_window,
                                                   chr_name,first,*length) == 0) {
                        fprintf(file_logerr,"Error processing weights file\n");
                        free(wP);free(wPV);free(wV);
                        exit(1);
                    }
                    init_site = wlimit_end;/*init for next window*/
                }
            }
		}
	}
	else {/*(nwindows>0)*//*that is, using coordinate positions*/
		beg = wgenes[wc++];/*the initial position of the window*/
		end = wgenes[wc++];/*the final position of the window*/
        weight_window = 0; /*physical positions*/
        window_size = (double)end-beg+1;
        if(file_ws != 0) {
            if(read_weights_positions_file(file_ws,file_ws_gz,index_w,
                                           file_logerr,file_logerr_gz,
                                           &wP,&wPV,&wV,&wlimit_end,
                                           beg,&window_size,&n_sitesw,weight_window,
                                           chr_name,first,*length) == 0) {
                fprintf(file_logerr,"Error processing weights file\n");
                free(wP);free(wPV);free(wV);
                exit(1);
            }
        }
        else {
            n_sitesw = end-beg+1+1; /*the position 0 does not count: we translate 1 right*/
            if((wP = (double *)realloc((double *)wP,n_sitesw*sizeof(double))) == 0) {
                fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.14 \n");
                free(wP);free(wPV);free(wV);
                return(0);
            }
            if((wPV = (double *)realloc((double *)wPV,n_sitesw*sizeof(double))) == 0) {
                fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.25 \n");
                free(wP);free(wPV);free(wV);
                return(0);
            }
            if((wV = (double *)realloc((double *)wV,n_sitesw*sizeof(double))) == 0) {
                fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.23 \n");
                free(wP);free(wPV);free(wV);
                return(0);
            }
            for(xx=0;xx<n_sitesw;xx++) {
                wP[xx] = wPV[xx] = wV[xx] = (float)1;
            }
        }
    }

    if(*vector_priors==0) {
        *npriors = 2;
        if((*vector_priors = (double *)calloc((long int)*npriors,sizeof(double)))==0) {
            fprintf(file_logerr,"Error: memory not allocated. get_tfadata.01");
            free(wP);free(wPV);free(wV);
            return(0);
        }
    }
    vector_priors[0][0] = (double)beg;
    vector_priors[0][1] = (double)end;
    
    /*READ TFASTA FILE*/
	/*define the init and the end site first! use slide and window if necessary: use also the weights if necessary*/	
	/*detect the end of the file! (*li=0,otherwise *li=-1)*/
	
    if((x=function_read_tfasta(file_input,input_gz,index_input,file_logerr,file_logerr_gz,beg, end,&n_sam,&n_site,&names,&DNA_matr,matrix_pol_tcga,chr_name,first,*length))==0) {
        fprintf(file_logerr,"Unable reading tfasta file\n");
        for(x=0;x<n_sam;x++) free(names[x]); free(names);
        free(wP);free(wPV);free(wV);
        free(DNA_matr);
        exit(1);
    }
    
    /*check that wP, wPV and wV (n_sitesw) have at least n_site positions!!! (if not, add more positions to the vector with zero values)*/
    if(n_site > n_sitesw) {
        if((wP = (double *)realloc((float *)wP,n_site*sizeof(double))) == 0) {
            fprintf(file_output,"\nError: memory not reallocated. get_obsdata.14 \n");
            exit(1);
        }
        if((wPV = (double *)realloc((float *)wPV,n_site*sizeof(double))) == 0) {
            fprintf(file_output,"\nError: memory not reallocated. get_obsdata.25 \n");
            exit(1);
        }
        if((wV = (double *)realloc((float *)wV,n_site*sizeof(double))) == 0) {
            fprintf(file_output,"\nError: memory not reallocated. get_obsdata.23 \n");
            exit(1);
        }
        for(xx=n_sitesw;xx<n_site;xx++) {
            wP[xx] = wPV[xx] = wV[xx] = (double)0;
        }
        n_sitesw = n_site;
    }
    
 
    if(x==-1 || (nwindows > 0 && wc == 2*nwindows))
        *li=0;
	else
        *li=-1;
	   
    n_samp = n_sam;
    *length = n_site;
    
    if(n_samp < nsamuser_eff) {
        free(wP);free(wPV);free(wV);
        return(0);
    }
    if(n_samp == 0 || n_site == 0) {
        free(wP);free(wPV);free(wV);
        return(0);
    }
    else {
        /*modify the order of samples using option flag O*/
        flag_change_sort = 0;
        for(x=0;x<nsamuser_eff;x++) {
            if(sort_nsam[x] != x) {
                flag_change_sort = 1;
                break;
            }
        }
        if(flag_change_sort == 1) {
            /*define duplicated matr*/
            if ((DNA_matr2 = (char *)calloc(n_site*(long long)n_samp,sizeof(char))) == 0) {
                fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.23d \n");
                for(x=0;x<n_samp;x++) free(names[x]); free(names);
                free(wP);free(wPV);free(wV);
                free(DNA_matr);
                return(0);
            }
            if((names2 = (char **)calloc(n_samp,sizeof(char *))) == 0) {
                fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.1s2 \n");
                for(x=0;x<n_samp;x++) free(names[x]); free(names);
                free(wP);free(wPV);free(wV);
                free(DNA_matr);
                return(0);
            }
            for(x=0;x<n_samp;x++) {
                if((names2[x] = (char *)calloc(50,sizeof(char))) == 0) {
                    fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.22 \n");
                    for(x=0;x<n_samp;x++) free(names[x]); free(names);
                    free(wP);free(wPV);free(wV);
                    free(DNA_matr);
                    return(0);
                }
                /*copy duplicated data*/
                strncpy(names2[x],names[x],50);
            }
            strncpy(DNA_matr2+(long long)n_site*(long long)0,DNA_matr+(long long)n_site*(long long)0,(long long)n_site*n_samp);
            //for(x=0;x<n_samp;x++) {
            //    strncpy(DNA_matr2+(long long)n_site*(long long)x,DNA_matr+(long long)n_site*(long long)x,n_site);
            //}
            /*end define and duplicate*/
            
            /*include data in *DNA_matr and in *names[] in the correct order*/
            for(x=0;x<nsamuser_eff;x++) {
                strncpy(DNA_matr+(long long)n_site*(long long)x,DNA_matr2+(long long)n_site*(long long)sort_nsam[x],n_site);
                strncpy(names[x],names2[sort_nsam[x]],50);
            }
            /*delete duplicated matr*/
            for(x=0;x<n_samp;x++) free(names2[x]);
            free(names2); names2 = 0;
            free(DNA_matr2);
            
            /*erase lines no used*/
            if(nsamuser_eff > n_samp) {
                fprintf(file_logerr,"Error: too low samples in the file according to defined in -N flag.\n");
                for(x=0;x<n_sam;x++) free(names[x]); free(names);
                free(wP);free(wPV);free(wV);
                free(DNA_matr);
                return(0);
            }
        }
        /*end option flag O*/
        if(nsamuser_eff > 32167) {
            fprintf(file_logerr,"Error: too much samples. Only 32167 samples per loci are allowed.\n");
            for(x=0;x<n_sam;x++) free(names[x]); free(names);
            free(wP);free(wPV);free(wV);
            free(DNA_matr);
            return(0);
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
        if(outgroup_presence == 0) {
            if ((DNA_matr = (char *)realloc(DNA_matr,(long long)n_site*(nsamuser_eff+!outgroup_presence)*sizeof(char))) == 0) {
                fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.23a \n");
                for(x=0;x<n_samp;x++) free(names[x]); free(names);
                free(wP);free(wPV);free(wV);
                free(DNA_matr);
//                free(matrix_sizepos);
//                free(matrix_segrpos);
                return(0);
            }
            int ns;
            int f1,f2,f3,f4,f0;
            int nf1,nf2,nf3,nf4,nfm,nf0,countf;
            for(xx=0;xx<n_site;xx++) {
                /*define f1 for a given nt*/
                ns=nfm=0;
                while(ns < nsamuser_eff-1 && DNA_matr[(long long)n_site*(unsigned long)ns+xx] > '4') {ns++;nfm++;}
                if(ns<nsamuser_eff-1) {
                    f1 = (int)DNA_matr[(long long)n_site*(unsigned long)ns+xx]; nf1 = 1;
                    /*count the freqs of all nt*/
                    nf2=nf3=nf4=0;
                    f2=f3=f4=0;
                    while(ns < nsamuser_eff-1)
                    {ns++;
                        if((int)DNA_matr[(long long)n_site*(unsigned long)ns+xx]==f1) nf1 += 1;
                        if(DNA_matr[(long long)n_site*(unsigned long)ns+xx]>'4') nfm += 1;
                        if((int)DNA_matr[(long long)n_site*(unsigned long)ns+xx]!=f1 &&
                           DNA_matr[(long long)n_site*(unsigned long)ns+xx]<'5') {
                            if(f2==0 || (int)DNA_matr[(long long)n_site*(unsigned long)ns+xx]==f2) {
                                f2 = (int)DNA_matr[(long long)n_site*(unsigned long)ns+xx];
                                nf2 +=1;
                            }
                            else {
                                if(f3==0 || (int)DNA_matr[(long long)n_site*(unsigned long)ns+xx]==f3) {
                                    f3 = (int)DNA_matr[(long long)n_site*(unsigned long)ns+xx];
                                    nf3 +=1;
                                }
                                else {
                                    if(f4==0 || (int)DNA_matr[(long long)n_site*(unsigned long)ns+xx]==f4) {
                                        f4 = (int)DNA_matr[(long long)n_site*(unsigned long)ns+xx];
                                        nf4 +=1;
                                    }
                                }
                            }
                        }
                    }
                    /*choose one of the two highest frequencies, not the lowest if more than 2 (give '5')*/
                    countf = 1;
                    f0 = f1; nf0 = nf1;
                    if(nf0 == nf2) {countf += 1;}
                    if(nf0  < nf2) {f0 = f2; nf0 = nf2; countf = 1;}
                    if(nf0 == nf3) {countf += 1;}
                    if(nf0  < nf3) {f0 = f3; nf0 = nf3; countf = 1;}
                    if(nf0 == nf4) {countf += 1;}
                    if(nf0  < nf4) {f0 = f4; nf0 = nf4; countf = 1;}
                    
                    if(countf < 3) DNA_matr[(unsigned long)n_site*(long long)(nsamuser_eff)+xx] = (char)f0;
                    else
                        DNA_matr[(unsigned long)n_site*(long long)(nsamuser_eff)+xx] = '5';
                }
                else {
                    DNA_matr[(unsigned long)n_site*(long long)(nsamuser_eff)+xx] = '5';
                }
                //ns = 0;
                //while(ns < nsamuser_eff-1 && DNA_matr[(long long)n_site*(unsigned long)ns+xx] > '4') ns++;
                //DNA_matr[(unsigned long)n_site*(long long)(nsamuser_eff)+xx] = DNA_matr[(unsigned long)n_site*(unsigned long)ns+xx];
            }
            nsamuser_eff += 1;
            /*we forced the invented outgroup without gaps or uncertainties, if possible*/
            //for(xx=0;xx<n_site;xx++) {
            //    int ns = 0;
            //    while(ns < nsamuser_eff-1 && DNA_matr[(long long)n_site*(unsigned long)ns+xx] > '4') ns++;
            //    DNA_matr[(unsigned long)n_site*(long long)(nsamuser_eff)+xx] = DNA_matr[(unsigned long)n_site*(unsigned long)ns+xx];
            //}
            /*strncpy(DNA_matr+(unsigned long)n_site*(unsigned long)(nsamuser_eff),DNA_matr+(unsigned long)n_site*(unsigned long)(nsamuser_eff-1),n_site);*/
            //nsamuser_eff += 1;
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
        if((mhitbp = (long int *) calloc (n_site, sizeof(long int))) == 0) {
            fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.6");
            for(x=0;x<n_sam;x++) free(names[x]); free(names);
            free(wP);free(wPV);free(wV);
            free(DNA_matr);
            return(0);
        }
        /*function to analyze all data*/
        if(get_obsstats(file_output,file_output_gz,(FILE *)0,file_logerr,file_logerr_gz,
                        nsamuser_eff,n_site,length_al_real,names/*2*/,
                        DNA_matr/*2*/,wP/*matrix_sizepos*/,wPV/*matrix_segrpos*/,
                        matrix_pol,matrix_freq,matrix_pos,length_al,length_seg,
                        nsamuser,npops,svratio,missratio,include_unknown,
                        sum_sam,tcga,matrix_sv,nmhits,output,ploidy,outgroup_presence,
                        nsites1_pop,nsites1_pop_outg,
                        nsites2_pop,nsites2_pop_outg,nsites3_pop,nsites3_pop_outg,
                        anx,bnx,anxo,bnxo,lengthamng,lengthamng_outg,mhitbp,
                        matrix_pol_tcga,(long int)beg) == 0) {
            for(x=0;x<n_sam;x++) free(names[x]); free(names);names = 0;
            free(wP);free(wPV);free(wV);
            free(DNA_matr);
            /*free(DNA_matr2);*/
/*          free(matrix_sizepos);
            free(matrix_segrpos);
*/          free(mhitbp);
            return(0);
        }
		/*free(names2);
		free(DNA_matr2);*/
        /*for(x=0;x<n_sam;x++) free(names[x]); free(names);*/
/*		free(matrix_sizepos);
		free(matrix_segrpos);
*/
        free(mhitbp);
        free(wP);free(wPV);free(wV);
    }
    return(1);
}

/*very careful!!!!!. The coordinates file is very strict in shape: name \t start \t end \n !!!!!!!!*/
int read_coordinates(FILE *file_wcoor, SGZip *file_wcoor_gz, FILE *file_output, SGZip *file_output_gz, FILE *file_logerr, SGZip *file_logerr_gz, long int **wgenes, long int *nwindows,char *chr_name) {
    
    char *valn=0;
    char c;
    int x;
    long int xx;
    long int dd;
    double ee;
    int inside_chr;
    
    /*printf("\nReading coordinates file...");*/
    fflush(stdout);
    fprintf(file_logerr,"\nReading coordinates file...");
    
    if((valn = (char *)calloc(100,sizeof(char))) == 0) {
        fprintf(file_logerr,"\nError: memory not reallocated. read_coordinates.00 \n");
        return(0);
    }
    if((*wgenes = (long int *)calloc(10000,sizeof(long int))) == 0) {
        fprintf(file_logerr,"\nError: memory not reallocated. read_coordinates.0 \n");
        free(valn);
        return(0);
    }
    
    *nwindows = 0;
    c = fzgetc(file_wcoor, file_wcoor_gz);
    
    if(check_comment(&c,file_wcoor,file_wcoor_gz) == 0) {
        fprintf(file_logerr,"\nWarning: no coordinates assigned. \n");
        free(valn);
        return(0);
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
    xx=0;
    while (c != 0 && c != -1 && c!='\xff' && c!='\xfe') {
        /*now keep all values: two columns, only numbers*/
        /*first column is the name of the scaffold*/
        while(c == 32 || c == 9 || c == 13 || c == 10)
            c = fzgetc(file_wcoor, file_wcoor_gz);
        if(check_comment(&c,file_wcoor, file_wcoor_gz) == 0) {
            c=0;break;
        }
        while(c == 32 || c == 9 || c == 13 || c == 10)
            c = fzgetc(file_wcoor, file_wcoor_gz);
        x=0;
        while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100  && c!='\xff' && c!='\xfe' && c>0) {
            valn[x] = c;
            c = fgetc(file_wcoor);
            x++;
        }
        valn[x] = '\0';/*scaffold name*/
        if(c==0 || c==-1 || c < 1  || c=='\xff' || c=='\xfe') {
            c=0;break;
        }
        while(strcmp(chr_name,valn) != 0) { /*discard those rows having different scaffold than chr_name*/
            if(inside_chr == 1) {
                inside_chr = -1; /*we passed target region*/
                break;
            }
            while(c != 13 && c != 10 && c!=0 && c!=-1  && c!='\xff' && c!='\xfe')
                c = fzgetc(file_wcoor, file_wcoor_gz);
            if(check_comment(&c,file_wcoor, file_wcoor_gz) == 0){
                free(valn);*nwindows = (xx)/2;return(1);
            }
            while(c == 32 || c == 9 || c == 13 || c == 10)
                c = fzgetc(file_wcoor, file_wcoor_gz);
            x=0;
            while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100  && c!='\xff' && c!='\xfe' && c>0) {
                valn[x] = c;
                c = fgetc(file_wcoor);
                x++;
            }
            valn[x] = '\0';/*scaffold name*/
            if(c==0 || c==-1 || c < 1  || c=='\xff' || c=='\xfe') {
                c=0;break;
            }
        }
        if(c==-1 || c == 0)
            break;
        if(inside_chr == -1)
            break; /*we go out*/
        inside_chr = 1;
        /*KEEP POSITIONS (first initial position, then end, next region and so on)*/
        while(c == 32 || c == 9 || c == 13 || c == 10)
            c = fzgetc(file_wcoor, file_wcoor_gz);
        if(check_comment(&c,file_wcoor, file_wcoor_gz) == 0) {
            c=0;free(valn);return(0);
        }
        while(c == 32 || c == 9 || c == 13 || c == 10)
            c = fzgetc(file_wcoor, file_wcoor_gz);
        x=0;
        while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100 && c!='\xff' && c!='\xfe' && c>0) {
            valn[x] = c;
            c = fzgetc(file_wcoor, file_wcoor_gz);
            x++;
        }
        valn[x] = '\0';
        wgenes[0][xx] = (long int)atof(valn);
        
        xx++;
        while(c == 32 || c == 9 || c == 13 || c == 10)
            c = fgetc(file_wcoor);
        if(check_comment(&c,file_wcoor, file_wcoor_gz) == 0) {
            c=0;free(valn);return(0);
        }
        while(c == 32 || c == 9 || c == 13 || c == 10)
            c = fzgetc(file_wcoor, file_wcoor_gz);
        x=0;
        while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100 && c!='\xff' && c!='\xfe' && c>0) {
            valn[x] = c;
            c = fgetc(file_wcoor);
            x++;
        }
        valn[x] = '\0';
        wgenes[0][xx] = (long int)round((double)atof(valn));
        
        xx++;
        dd = (long int)floor((double)xx/(double)10000);
        ee = (double)xx/(double)10000;
        if(dd==ee) {
            if((*wgenes = realloc(*wgenes,((long int)(dd+1)*(long int)10000*(long int)sizeof(long int)))) == 0) {
                puts("Error: realloc error read_coordinates.1\n");
                free(valn);/*free(*wgenes);*/
                return(0);
            }
        }
    }
    if(xx == 0) {
        fprintf(file_logerr,"\nError: no coordinates assigned, read_coordinates.2 \n");
        /*free(*wgenes);
         file_wcoor=0;*/
        return(0);
    }
    *nwindows = (xx)/2;
    /*}*/
    free(valn);
    return 1;
}

int read_weights_positions_file(FILE *file_ws, SGZip *file_ws_gz, struct SGZIndex *index_w,FILE *file_logerr,SGZip *file_logerr_gz, double **wP, double **wPV, double **wV, long int *wlimit_end,long int init_site, double *window_size, long int *n_sitesw, int weight_window, char *chr_name,int first, long int length) {
    
    double curr_window=0;
    char *valn=0;
    long int dd;
    long int xx;
    double ee;
    int x;
    int y;

    static char c[1];
    static char line[MSP_MAX_COL_LINE];
    static char line2[MSP_MAX_COL_LINE];
    static int col=0;
    static int col2=0;
    static long int position=0;
    static char *ID = 0;
    static int count0s,nchstr;
    static long int row_num = -1;  /* fzseek: Set to -1 if you want to search by ID. */
                                    /*Or set NULL to the ID if you want to seach by position..*/

    /*read weights file*/
    /*keep wP, wPV, wV (not used yet) for all positions, do not need to keep positions: all are correlative*/
    if(file_ws != 0) {
        if((valn = (char *)calloc(100,sizeof(char))) == 0) {
            fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.34 \n");
            return(0);
        }
        if(ID==0) {
            if((ID = (char *)calloc(100,sizeof(char))) == 0) {
                fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.34 \n");
                free(valn);
                return(0);
            }
        }
        if(position == 0) { /*assign the value of count0s and nchrstr*/
            /*pass scaffold name:*/
            *c = fzgetc(file_ws, file_ws_gz);
            while(*c==10 || *c==13 || *c == 9 || *c == 32)
                *c = fzgetc(file_ws, file_ws_gz);
            col = 0;
            while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c != 9 && *c != 32 && *c != 58 && *c!='\xff' && *c!='\xfe') {
                line[col] = *c;
                col++;
                *c = fzgetc(file_ws, file_ws_gz);
            }
            if(*c == 0 || *c==-1  || *c=='\xff' || *c=='\xfe') {
                fprintf(file_logerr,"\nError: no weights assigned for %s scaffold \n",chr_name);
                free(valn);
                return(0);
            }
            line[col] = '\0';
            /*pass position:*/
            if(*c==58) {
                *c = fzgetc(file_ws, file_ws_gz);
                y=0; count0s=0;col2 = 0;
                while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c != 9 && *c != 32 && *c!='\xff' && *c!='\xfe') {
                    /*count 0s at left*/
                    if(*c == '0' && y == 0)
                        count0s += 1;
                    else y += 1;
                    /** in case count0s is > 0, count0s+y is the total number of characters
                     in the string to search for positions (add zeros at left) */
                    
                    line2[col2] = *c;
                    col2++;
                    *c = fzgetc(file_ws, file_ws_gz);
                }
                nchstr = y+count0s; /*the number of characters in the string*/
                
                if(*c == 0 || *c==-1  || *c=='\xff' || *c=='\xfe') {
                    fprintf(file_logerr,"\nError: no weights assigned for %s scaffold \n",chr_name);
                    free(valn);
                    return(0);
                }
                line2[col2] = '\0';
                position = atol(line2);
            }
            else {
                fprintf(file_logerr,"\nError: no position assigned for %s scaffold at position %ld \nHave you included additional comments?\n",chr_name,init_site);
                free(valn);
                return(0);
            }
        }
 
        /*Search for the value beg into ID: it is transformed to a string with the scaffold name */
        /* and perhaps zeros at left if count0s > 0*/
        if(transform_beg_chr(ID,chr_name,init_site,nchstr,count0s) != 1) {
            fprintf(file_logerr,"Error transforming beg into string.\n");
            free(valn);
            return(0);
        }
        
        row_num = -1;
        if(fzseekNearest(file_ws, file_ws_gz,index_w, ID, MAXLEN, &row_num) != GZ_OK) { //==GZ_ERROR_DATA_FILE?
            if(init_site <= length)
                fprintf(file_logerr,"ID not found in the weights file (or not indexed position): %s\n",ID);
            /*no position found. Assume the file window is finished*/
            free(valn);
            *wlimit_end = init_site;
            *window_size = 0;
            return(1);
        }

        /*get scaffold name: perhaps is not chr_name*/
        *c = fzgetc(file_ws, file_ws_gz);
        while(*c==10 || *c==13 || *c == 9 || *c == 32)
            *c = fzgetc(file_ws, file_ws_gz);
        col = 0;
        while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c != 9 && *c != 32 && *c != 58 && *c!='\xff' && *c!='\xfe') {
            line[col] = *c;
            col++;
            *c = fzgetc(file_ws, file_ws_gz);
        }
        if(*c == 0 || *c==-1  || *c=='\xff' || *c=='\xfe') {
            fprintf(file_logerr,"\nError: no weights assigned for %s scaffold \n",chr_name);
            free(valn);
            return(0);
        }
        line[col] = '\0';
        /*get position: perhaps is not init_site*/
        if(*c==58) {
            *c = fzgetc(file_ws, file_ws_gz);
            y=0; count0s=0;col2 = 0;
            while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c != 9 && *c != 32 && *c!='\xff' && *c!='\xfe') {
                /*count 0s at left*/
                if(*c == '0' && y == 0)
                    count0s += 1;
                else y += 1;
                /** in case count0s is > 0, count0s+y is the total number of characters
                 in the string to search for positions (add zeros at left) */
                
                line2[col2] = *c;
                col2++;
                *c = fzgetc(file_ws, file_ws_gz);
            }
            nchstr = y+count0s; /*the number of characters in the string*/
            
            if(*c == 0 || *c==-1  || *c=='\xff' || *c=='\xfe') {
                fprintf(file_logerr,"\nError: no weights assigned for %s scaffold \n",chr_name);
                free(valn);
                return(0);
            }
            line2[col2] = '\0';
            position = atol(line2);
        }
        else {
            fprintf(file_logerr,"\nError: no position assigned for %s scaffold at position %ld \nHave you included additional comments?\n",chr_name,init_site);
            free(valn);
            return(0);
        }
        
        if(check_comment(c,file_ws,file_ws_gz) == 0) {
            free(valn);
            return(0);
        }
        if(!(*c == 0 || *c==-1 || *c=='\xff' || *c=='\xfe'))
            *c = fzgetc(file_ws, file_ws_gz);
        

        xx = 0;/*xx=0*/ /* VERY IMPORTANT: SITES FOR WEIGHTS START FROM 0 AND NOT FROM 1 !*/
        curr_window = 0;
        col = 0;
        while (strcmp(line, chr_name) == 0 && curr_window < window_size[0]) {
            /*Weight Position*/
            while(*c == 32 || *c == 9)/* || *c == 13 || *c == 10)*/
                *c = fzgetc(file_ws,file_ws_gz);
            if(*c==0 || *c==-1)
                break;
            x=0;
            while(*c != 32 && *c != 9 && *c != 13 && *c != 10 && *c!=0 && *c!=-1 && x < 100) {
                valn[x] = *c;
                *c = fzgetc(file_ws,file_ws_gz);
                x++;
            }
            valn[x] = '\0';
            wP[0][xx] = (double)atof(valn);
            
            /*Weight Variant*/
            while(*c == 32 || *c == 9)/* || *c == 13 || *c == 10)*/
                *c = fzgetc(file_ws,file_ws_gz);
            if(*c==0 || *c==-1)
                break;
            x=0;
            while(*c != 32 && *c != 9 && *c != 13 && *c != 10 && *c!=0 && *c!=-1 && x < 100) {
                valn[x] = *c;
                *c = fzgetc(file_ws,file_ws_gz);
                x++;
            }
            valn[x] = '\0';
            wPV[0][xx] = (double)atof(valn);
            
            while(*c == 32 || *c == 9) *c = fzgetc(file_ws,file_ws_gz);
            if(!(*c == 13 || *c == 10 || *c == 0 || *c == -1 || *c == '#')) {
                /*Effect size*/
                while(*c == 32 || *c == 9)/* || *c == 13 || *c == 10)*/
                    *c = fzgetc(file_ws,file_ws_gz);
                if(*c==0 || *c==-1)
                    break;
                x=0;
                while(*c != 32 && *c != 9 && *c != 13 && *c != 10 && *c!=0 && *c!=-1 && x < 100) {
                    valn[x] = *c;
                    *c = fzgetc(file_ws,file_ws_gz);
                    x++;
                }
                valn[x] = '\0';
                wV[0][xx] = (double)atof(valn);
            }
            else {
                wV[0][xx] = 1.0; /*if undefined the value is 1.0 for all*/
            }

            if(check_comment(c,file_ws,file_ws_gz) == 0) {
                free(valn);
                break;
            }
            
            /*sum window sizes and positions*/
            if(weight_window)
                curr_window += wP[0][xx] * wV[0][xx];
            else
                curr_window += 1.0;
            xx++;
            
            /*collect a new position (chr:pos)*/
            if(*c < 0 || *c==10 || *c==13 || *c == -1 || *c == 0) {
                if(check_comment(c,file_ws,file_ws_gz) == 0) {
                    free(valn);
                    *wlimit_end = init_site + xx/* - 1*/;
                    *window_size = curr_window;
                    *n_sitesw = xx-1;
                    return(1);
                }
                col = 0;
                while(*c==10 || *c==13 || *c == 9 || *c == 32) {
                    *c = fzgetc(file_ws, file_ws_gz);
                    col++;
                    if(col >= MSP_MAX_COL_LINE) {
                        fprintf(file_logerr,"Register too large: position after %s:%ld, size %d\n ",chr_name,position,col);
                        free(valn);
                        *wlimit_end = init_site + xx - 1;
                        *window_size = curr_window;
                        *n_sitesw = xx-1;
                        return(1);
                    }
                }
                while(check_comment(c,file_ws,file_ws_gz) == 0) {
                    free(valn);
                    *wlimit_end = init_site + xx - 1;
                    *window_size = curr_window;
                    *n_sitesw = xx-1;
                    return(1);
                }
                col = 0;
                while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c != 9 && *c != 32 \
                      && *c!='\xff' && *c!='\xfe' && *c != 58 \
                      && *c>0  && col < MSP_MAX_COL_LINE-1) {
                    line[col] = *c;
                    col++;
                    *c = fzgetc(file_ws,file_ws_gz);
                }
                if(col >= MSP_MAX_COL_LINE) {
                    fprintf(file_logerr,"Register too large: position after %s:%ld, size %d \n",chr_name,position,col);
                    free(valn);
                    *wlimit_end = init_site + xx - 1;
                    *window_size = curr_window;
                    *n_sitesw = xx-1;
                    return(1);
                }
                line[col] = '\0';
                if(check_comment(c,file_ws,file_ws_gz) == 0) {
                    free(valn);
                    return(0);
                }
                col = 0;
                if(*c == 58) {
                    while(*c==10 || *c==13 || *c==58) *c = fzgetc(file_ws,file_ws_gz);
                    col2 = 0;
                    while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c != 9 && *c != 32 \
                          && *c!='\xff' && *c!='\xfe' \
                          && *c>0 && col < MSP_MAX_COL_LINE-1) {
                        line2[col2] = *c;
                        col2++;
                        *c = fzgetc(file_ws,file_ws_gz);
                    }
                    if(col2 >= MSP_MAX_COL_LINE) {
                        fprintf(file_logerr,"Register too large: position after %s:%ld, size %d\n ",chr_name,position,col2);
                        free(valn);
                        *wlimit_end = init_site + xx - 1;
                        *window_size = curr_window;
                        *n_sitesw = xx-1;
                        return(1);
                    }
                    col2 = 0;
                    position = atol(line2);
                    if(check_comment(c,file_ws,file_ws_gz) == 0) {
                        free(valn);
                        *wlimit_end = init_site + xx - 1 - 1;
                        *window_size = curr_window;
                        *n_sitesw = xx-1;
                        return(1);
                    }
                }
            }
            dd = (long int)floor((double)xx/(double)1000);
            ee = (double)xx/(double)1000;
            if(dd==ee) {
                if((*wP = realloc(wP[0],((long int)(dd+1)*(long int)1000*sizeof(double)))) == 0) {
                    file_ws=0;*wP=0;*wPV=0;*wV=0;
                    free(*wP);free(*wPV);free(*wV);
                    free(valn);
                    fprintf(file_logerr,"Error: realloc error get_obsdata.11\n");
                    return(0);
                }
                if((*wPV = realloc(wPV[0],((long int)(dd+1)*(long int)1000*sizeof(double)))) == 0) {
                    file_ws=0;*wP=0;*wPV=0;*wV=0;
                    free(*wP);free(*wPV);free(*wV);
                    free(valn);
                    fprintf(file_logerr,"Error: realloc error get_obsdata.11\n");
                    return(0);
                }
                if((*wV = realloc(wV[0],((long int)(dd+1)*(long int)1000*sizeof(double)))) == 0) {
                    file_ws=0;*wP=0;*wPV=0;*wV=0;
                    free(*wP);free(*wPV);free(*wV);
                    free(valn);
                    fprintf(file_logerr,"Error: realloc error get_obsdata.11\n");
                    return(0);
                }
            }
        }
        
        free(valn);
        *wlimit_end = init_site + xx;
        *window_size = curr_window;
        *n_sitesw = xx;
    }
    else return(0);
    
    return(1);
}

int function_read_tfasta(FILE *file_input,SGZip *file_input_gz,struct SGZIndex *index_input,FILE *file_logerr,SGZip *file_logerr_gz,long int init_site,long int end_site,int *n_sam, long int *n_site, char ***names, char **DNA_matr,char **matrix_pol_tcga,char *chr_name,int first, long int length)
{
    static char c[1];
    char *cc;
    static char line[MSP_MAX_COL_LINE];
    static char line2[MSP_MAX_COL_LINE];
    static int col=0;
    static int col2=0;
    static int nseq=0;
    int x;
    static int maxsam=128;
    static long int position=0;
    long int end_position=end_site;
    long int dd,count,xx;
    double ee;
    char *DNA_matr2;
    static long int row_num = -1;  /* fzseek: Set to -1 if you want to search by ID. */
                                    /*Or set NULL to the ID if you want to seach by position..*/
    static char *ID = 0;
    /*long int f_num = 0;*/
    int y;
    static int count0s,nchstr;
    /*static int chr_defined=0;*/
    /*static int position_defined=0;*/
    
    if ((DNA_matr2 = (char *)calloc(1e6,sizeof(char))) == 0) {
        fprintf(file_logerr,"\nError: memory not reallocated. get_tfadata.23d \n");
        free(DNA_matr2);
        return(0);
    }
    if(ID==0) {
        if((ID = (char *)calloc(100,sizeof(char))) == 0) {
            fprintf(file_logerr,"\nError: memory not reallocated. get_obsdata.34 \n");
            free(DNA_matr2);
            return(0);
        }
    }

    if(position == 0 && first == 0) {/*to allow sliding windows non-overlapped*/
        /*if names and number of samples are defined then skip the definition of names*/
        *c = fzgetc(file_input, file_input_gz);
        while(*c == 9 || *c == 10 || *c == 13)
            *c = fzgetc(file_input, file_input_gz);
        if((*c == -1 || *c == 0 || *c=='\xff' || *c=='\xfe')) {
            fprintf(file_logerr,"\nError: no sequence assigned for %s scaffold \n",chr_name);
            free(DNA_matr2);
            return 0;
        }
        line[col] = *c;
        col++;
        
        while(*c=='#') {
            /*parse comments*/
            while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c!='\xff' && *c!='\xfe') {
                *c = fzgetc(file_input, file_input_gz);
                line[col] = *c;
                col++;
            }
            line[col] = '\0';
            if((cc=strstr(line, "#NAMES:")) != 0) {
                /*collect names*/
                nseq = 0;
                cc = strtok (line,">\n\r ");
                while (cc != NULL)
                {
                    if(strstr(cc,"#NAMES:") == 0) {
                        strcpy(names[0][nseq],cc);
                        nseq++;
                        if(nseq == maxsam) {
                            maxsam += 128;
                            if(maxsam > 32767) {
                                fprintf(file_logerr,"\n Sorry, no more than 32767 samples are allowed.");
                                free(DNA_matr2);
                                return 0;
                            }
                            if ((*names = (char **)realloc(*names,maxsam*sizeof(char *))) == 0) {
                                fprintf(file_logerr,"\nError: memory not reallocated. assigna.1 \n");
                                free(DNA_matr2);
                                return(0);
                            }
                            for(x=nseq;x<maxsam;x++) {
                                if ((names[0][x] = (char *)calloc(50,sizeof(char))) == 0) {
                                    fprintf(file_logerr,"\nError: memory not reallocated. assigna.2 \n");
                                    free(DNA_matr2);
                                    return(0);
                                }
                            }
                        }
                    }
                    cc = strtok(NULL, ">\n\r ");
                }
            }
            col=0;
            *c = fzgetc(file_input, file_input_gz);
            line[col] = *c;
            col++;
        }
        
        /*assign count0s and nchrstr*/
        /*first pass scaffold name:*/
        col = 0;
        while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c != 9 && *c != 32 && *c != 58 && *c!='\xff' && *c!='\xfe') {
            line[col] = *c;
            col++;
            *c = fzgetc(file_input, file_input_gz);
        }
        if(*c == 0 || *c==-1  || *c=='\xff' || *c=='\xfe') {
            free(DNA_matr2);
            return(-1);
        }
        line[col] = '\0';
        /*second: pass position value. Count 0s*/
        if(*c==58) {
            *c = fzgetc(file_input, file_input_gz);
            y=0; count0s=0;col2 = 0;
            while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c != 9 && *c != 32 && *c!='\xff' && *c!='\xfe') {
                //count 0s at left
                if(*c == '0' && y == 0)
                    count0s += 1;
                else y += 1;
                // in case count0s is > 0, count0s+y is the total number of characters
                // in the string to search for positions (add zeros at left)
                
                line2[col2] = *c;
                col2++;
                *c = fzgetc(file_input, file_input_gz);
            }
            nchstr = y+count0s; //the number of characters in the string
            
            if(*c == 0 || *c==-1  || *c=='\xff' || *c=='\xfe') {
                free(DNA_matr2);
                return(-1);
            }
            line2[col2] = '\0';
            position = atol(line2);
        }
    }
    n_sam[0] = nseq;
    /*
    if(check_comment(c,file_input,file_input_gz) == 0) {
        free(DNA_matr2);
        return(-1);
    }
    */
    /*Search for the value beg into ID: it is transformed to a string with the scaffold name */
    /* and perhaps zeros at left if count0s > 0*/
    if(transform_beg_chr(ID,chr_name,init_site,nchstr,count0s) != 1) {
        fprintf(file_logerr,"Error transforming beg into string.\n");
        free(DNA_matr2);
        return(0);
    }
    /*f_num = position;*/
    row_num = -1;
    
    if(fzseekNearest(file_input, file_input_gz,index_input, ID, MAXLEN, &row_num) != GZ_OK) { //==GZ_ERROR_DATA_FILE?
        if(init_site <= length)
            fprintf(file_logerr,"ID not found in the tfa file (or not indexed position): %s\n",ID);
        /*no position found. Assume the file window is finished*/
        free(DNA_matr2);
        return(1);
    }
    /*f_num = row_num;*/

    /*get scaffold name*/
    *c = fzgetc(file_input, file_input_gz);
    while(*c=='#') {
        while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c!='\xff' && *c!='\xfe')
            *c = fzgetc(file_input, file_input_gz);
        *c = fzgetc(file_input, file_input_gz);
    }
    col = 0;
    while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c != 9 && *c != 32 && *c != 58 && *c!='\xff' && *c!='\xfe') {
        line[col] = *c;
        col++;
        *c = fzgetc(file_input, file_input_gz);
    }
    if(*c == 0 || *c==-1  || *c=='\xff' || *c=='\xfe') {
        /*fprintf(file_logerr,"\nError: no sequence assigned for %s scaffold \n",chr_name);*/
        free(DNA_matr2);
        return(-1);
    }
    line[col] = '\0';
    /*get position*/
    if(*c==58) {
        *c = fzgetc(file_input, file_input_gz);
        y=0; count0s=0;col2 = 0;
        while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c != 9 && *c != 32 && *c!='\xff' && *c!='\xfe') {
            /*count 0s at left*/
            if(*c == '0' && y == 0)
                count0s += 1;
            else y += 1;
            /** in case count0s is > 0, count0s+y is the total number of characters
             in the string to search for positions (add zeros at left) */
            
            line2[col2] = *c;
            col2++;
            *c = fzgetc(file_input, file_input_gz);
        }
        nchstr = y+count0s; /*the number of characters in the string*/
        
        if(*c == 0 || *c==-1  || *c=='\xff' || *c=='\xfe') {
            /*fprintf(file_logerr,"\nError: no sequence assigned for %s scaffold \n",chr_name);*/
            free(DNA_matr2);
            return(-1);
        }
        line2[col2] = '\0';
        position = atol(line2);
    }
    else {
        /*fprintf(file_logerr,"\nError: no position assigned for %s scaffold at row %ld \n",chr_name,row_num);*/
        free(DNA_matr2);
        return(-1);
    }
    
    if(check_comment(c,file_input,file_input_gz) == 0) {
        free(DNA_matr2);
        return(-1);
    }
    if(!(*c == 0 || *c==-1 || *c=='\xff' || *c=='\xfe'))
        *c = fzgetc(file_input, file_input_gz);

    *n_site = 0;
    col = 0;
    count=0;
    while (strcmp(line, chr_name) == 0 && position <= end_position) {
        /*chr_defined = 0;*/
        /*position_defined = 0;*/
        switch(*c) {
            case 'T':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '1';
                col += 1;
                count += 1;
                break;
            case 't':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '1';
                col += 1;
                count += 1;
                break;
            case 'U':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '1';
                col += 1;
                count += 1;
                break;
            case 'u':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '1';
                col += 1;
                count += 1;
                break;
            case 'C':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '2';
                col += 1;
                count += 1;
                break;
            case 'c':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '2';
                col += 1;
                count += 1;
                break;
            case 'G':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '3';
                col += 1;
                count += 1;
                break;
            case 'g':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '3';
                col += 1;
                count += 1;
                break;
            case 'A':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '4';
                col += 1;
                count += 1;
                break;
            case 'a':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '4';
                col += 1;
                count += 1;
                break;
            case 0:
                break;
            case -1:
                break;
            case 10:
                break;
            case 13:
                break;
            case 32:
                break;
            case 9:
                break;
            case 'N':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '5';
                col += 1;
                count += 1;
                break;
            case 'n':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '5';
                col += 1;
                count += 1;
                break;
            case '?':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '5';
                col += 1;
                count += 1;
                break;
            case '-':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '6';
                col += 1;
                count += 1;
                break;
                /*gaps are converted to N*/
            case 'W':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = 'w';
                col += 1;
                count += 1;
                break;
            case 'w':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = 'w';
                col += 1;
                count += 1;
                break;
            case 'M':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = 'm';
                col += 1;
                count += 1;
                break;
            case 'm':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = 'm';
                col += 1;
                count += 1;
                break;
            case 'R':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = 'r';
                col += 1;
                count += 1;
                break;
            case 'r':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = 'r';
                col += 1;
                count += 1;
                break;
            case 'Y':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = 'y';
                col += 1;
                count += 1;
                break;
            case 'y':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = 'y';
                col += 1;
                count += 1;
                break;
            case 'K':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = 'k';
                col += 1;
                count += 1;
                break;
            case 'k':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = 'k';
                col += 1;
                count += 1;
                break;
            case 'S':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = 's';
                col += 1;
                count += 1;
                break;
            case 's':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = 's';
                col += 1;
                count += 1;
                break;
            case 'b':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '5';
                col += 1;
                count += 1;
                break;
            case 'B':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '5';
                col += 1;
                count += 1;
                break;
            case 'd':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '5';
                col += 1;
                count += 1;
                break;
            case 'D':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '5';
                col += 1;
                count += 1;
                break;
            case 'h':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '5';
                col += 1;
                count += 1;
                break;
            case 'H':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '5';
                col += 1;
                count += 1;
                break;
            case 'v':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '5';
                col += 1;
                count += 1;
                break;
            case 'V':
                DNA_matr2[(((long long)nseq*(unsigned long)*n_site)+(unsigned long)col)] = '5';
                col += 1;
                count += 1;
                break;
            default:
                fprintf(file_logerr,"Unexpected value in tfa file: position %ld, sample %d \n ",position,col);
                fprintf(file_logerr,"%c\n",*c);
                free(DNA_matr2);
                return(-1);
                break;
        }
        if(check_comment(c,file_input,file_input_gz) == 0) {
            break;
        }
        *c = fzgetc(file_input, file_input_gz);
        if(*c < 0 || *c==10 || *c==13 || *c == -1 || *c == 0 || *c=='\xff' || *c=='\xfe')  {
            if(check_comment(c,file_input,file_input_gz) == 0) {
                if(col == *n_sam || col == 0) {
                    *n_site += 1;
                    break;
                }
                else {
                    free(DNA_matr2);
                    return(-1);}
            }
            while(*c==10 || *c==13)
                *c = fzgetc(file_input, file_input_gz);
            if(check_comment(c,file_input,file_input_gz) == 0) {
                break;
            }
            *n_site += 1;
            /*if(position != *n_site) {
                printf("check here");
            }*/
            /*read new line*/
            col = 0;
            while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c != 9 && *c != 32 \
                  && *c!='\xff' && *c!='\xfe' && *c != 58 \
                  && *c>0 && col < MSP_MAX_COL_LINE-1) {
                line[col] = *c;
                col++;
                *c = fzgetc(file_input, file_input_gz);
            }
            line[col] = '\0';
            /*chr_defined = 1;*/
            if(check_comment(c,file_input,file_input_gz) == 0) {
                break;
            }
            col = 0;
            if(*c == 58) {
                while(*c==10 || *c==13 || *c==58) *c = fzgetc(file_input, file_input_gz);
                col2 = 0;
                while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c != 9 && *c != 32 && *c>0) {
                    line2[col2] = *c;
                    col2++;
                    *c = fzgetc(file_input, file_input_gz);
                }
                line2[col2] = '\0';
                col2 = 0;
                position = atol(line2);
                if(end_site == -1) end_position = position + 1;
                /*position_defined = 1;*/
                if(check_comment(c,file_input,file_input_gz) == 0) {
                    free(DNA_matr2);
                    return(-1);
                }
            }
            /*else {
                printf("\nError: Problem reading line %s:%ld\n",line,position);
            }*/
        }
        /*realloc DNAmatr if nnecessary*/
        dd = (long int)floor((double)count/(double)1e6);
        ee = (double)count/(double)1e6;
        if((double)dd == ee) {
            if((DNA_matr2 = (char *)realloc((char *) DNA_matr2,((long int)dd+(long int)1)*(long int)1e6*sizeof(char))) == 0) {
                fprintf(file_logerr,"Error: realloc error varchar.1\n");
                return(0);
            }
        }
    }
    /**n_site += 1;*/
    /**/
    if(*c==0 || *c==-1 || *c =='\xff' || *c=='\xfe') {
        if(col == *n_sam || col == 0) *n_site += 1;
        else return(-1);
    }
    /**/
    /*return n_site,n_sam,DNA matrix and names in pointers*/
    /*
    if ((*DNA_matr = (char *)realloc((char *) *DNA_matr,n_site[0]*(long long)n_sam[0]*sizeof(char))) == 0) {
        fprintf(file_logerr,"\nError: memory not reallocated. get_tfadata.23d \n");
        free(DNA_matr2);
        return(0);
    }
    */
    free(*DNA_matr);
    if((*DNA_matr = (char *)calloc(n_site[0]*(long long)n_sam[0],sizeof(char))) == 0) {
        fprintf(file_logerr,"\nError: memory not reallocated. get_tfadata.23d2 \n");
        free(DNA_matr2);
        return(0);
    }
    for(x=0;x<n_sam[0];x++) {
        for(xx=0;xx<n_site[0];xx++) { /*transpose */
            DNA_matr[0][(((long long)n_site[0]*(unsigned long)x)+(unsigned long)xx)] =
            DNA_matr2 [(((long long)n_sam[0]*(unsigned long)xx)+(unsigned long)x)];
        }
    }
    free(DNA_matr2);
    
    if(check_comment(c,file_input,file_input_gz) == 0) {
        return(-1);
    }
    return(1);
}


int transform_beg_chr(char *ID, char *chr_name,long int beg,int nchstr, int count0s)
{
    /*Search for the value beg: it may be transformed to a string with zeros at left if count0s > 0*/
    /*first count the number of digits in beg and transform beg in a string*/
    long int zi=beg;
    double zr;
    int p10/*,p10i*/;
    int x/*,nnchstr*/;
    char ID2[100];
    
    memset(ID2,0,100);
    
    p10=0;
    do{
        p10 += 1;
        zr = (double)zi/(double)pow((double)10,(double)p10);
    }while(zr > 1.0);
    
    if(count0s > 0) {/*transformed to a string with zeros at left if count0s > 0*/
        /*p10 is the number of numbers after the zeros*/
        for(x=nchstr-1;x>=p10;x--) {
            ID2[nchstr-1-x] = '0';
        }
        /*nnchstr = nchstr;*/
    }/*
    else {
        nnchstr = p10;
    }
    */
    sprintf(ID2,"%s%ld",ID2,beg);
    
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
    
    strcpy(ID,chr_name);
    strcat(ID,":");
    strcat(ID,ID2);
    
    return(1);
}
int  check_comment(char *c, FILE *file_input, SGZip *file_input_gz) {
    if(*c=='#') {
        /*parse comments*/
        while(*c!=10 && *c!=13 && *c != 0 && *c!=-1 && *c!='\xff' && *c!='\xfe') {
            *c = fzgetc(file_input, file_input_gz);
        }
        if(*c==0 || *c==-1 || *c =='\xff' || *c=='\xfe')
            return(0);
    }
    if(*c==0 || *c==-1 || *c =='\xff' || *c=='\xfe')
        return(0);
    
    return(1);
}


