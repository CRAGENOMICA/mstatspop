//
//  get_tfadata.c
//  xcode_project
//
//  Created by Sebastian Ramos-Onsins on 09/08/15.
//
//

#include "get_tfadata.h"

int get_tfadata(FILE *file_output,
                FILE *file_input,
                SGZip *input_gz,
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
                double **length_amng_outg,
                int *sort_nsam,
                
                float *wV,
                long int *Pp,
                long int nV,
				float *wP,
				float *wPV,
                long int *wgenes,
                long int nwindows,
                long int first_slide,
                long int slide,/**/
                long int window,/**/
				int Physical_length,
				long int *li/**/,
                int *npriors,
                double **vector_priors,
                long int wlimit_end,
                long int welimit_end,
                char **matrix_pol_tcga
				)
{
    char *DNA_matr2 = 0;
    char **names2   = 0; /* limit of the name of the lines to 50 characters. be careful ... */
    
	double *matrix_sizepos = 0; /*size of each position, 0 not include, >0 include, Syn/Nsyn positions are allowed*/
    double *matrix_segrpos = 0; /*always 1 except in those biallelic positions that are not desired (e.g., choose syn but we do not want nsyn)*/
    
    long int *mhitbp; /*vector of positions where mhits are*/
    /* long int count; */
    long int xx;
    /* int c; */
    int n_sam/*,ns*/;
    long int n_sit;
    /* int nseq; */
    /* int maxsam; */
    /* int n_excl; */
    int n_samp;
    long int n_site;
	char ploidy[2];
    
    int x;
    static long int maxsites = 0;
    
    int nsamtot;
    int nsamuser_eff;
    int flag_change_sort; /*in case the order of samples is not consecutive*/
    
    static char *DNA_matr = 0;
    static char **names   = 0; /* limit of the name of the lines to 50 characters. be careful ... */
	
    static long int init_site = 1;
	static long wc = 0;
	long int lx,beg,end;
	double wl;

    ploidy[0]='1';
    ploidy[1]='\0';
    /* count = 0; */
    /* c = 0; */
    n_sam = 0;
    n_sit = 0;
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
            fprintf(file_output,"\nError: memory not reallocated. get_obsdata.1 \n");
            return(0);
        }
        for(x=0;x<128;x++) {
            if((names[x] = (char *)calloc(50,sizeof(char))) == 0) {
                fprintf(file_output,"\nError: memory not reallocated. get_obsdata.2 \n");
                return(0);
            }
        }
        if((DNA_matr = (char *)calloc(10000,sizeof(char))) == 0) {
            for(x=0;x<128;x++) free(names[x]); free(names);
            fprintf(file_output,"\nError: memory not reallocated. get_obsdata.3 \n");
            return(0);
        }
    }
    
    init_site += first_slide;
    
	/*FIND THE WINDOW OF POSITIONS TO ANALYZE*/
	if(nwindows==0) {
		if(Physical_length == 1) {
			beg = init_site;/*the initial position of the window*/
			end = init_site + window - 1;/*the final position of the window*/
			init_site += slide;/*init for next window (static variable)*/
		}
		else {/*if(Physical_length == 0)*/
			beg = init_site;/*the initial position of the window*/
			end = beg;
			wl=0.0;
			while(wl<window && end <= wlimit_end) {
				wl += wP[end-1] * wPV[end-1];
				end += 1;
			}
			end -= 1;/*the final position of the window (static variable)*/
			
			lx = end;
			wl=0.0;
			while(wl<slide && lx <= wlimit_end) {
				wl += wP[lx-1] * wPV[lx-1];
				lx += 1;
			}
			init_site = lx;/*init for next window*/
		}
	}
	else {/*(nwindows>0)*//*that is, using coordinate positions*/
		beg = wgenes[wc++];/*the initial position of the window*/
		end = wgenes[wc++];/*the final position of the window*/
	}
   if(*vector_priors==0) {
        *npriors = 2;
        if((*vector_priors = (double *)calloc((long int)*npriors,sizeof(double)))==0) {
            puts("Error: memory not allocated. get_tfadata.01");
            return(1);
        }
    }
    vector_priors[0][0] = (double)beg;
    vector_priors[0][1] = (double)end;
    
    /*READ TFASTA FILE*/
	/*define the init and the end site first! use slide and window if necessary: use also the weights if necessary*/	
	/*detect the end of the file! (*li=0,otherwise *li=-1)*/
	
    if((x=function_read_tfasta(file_input,input_gz,beg,end,&n_sam,&n_site,&names,&DNA_matr,matrix_pol_tcga))==0) {
        printf("Unable reading tfasta file\n");
        exit(1);
    }
    
    if(x==-1) *li=0;
	else *li=-1;
	   
    n_samp = n_sam;
    *length = n_site;
    
    if(n_samp < nsamuser_eff)
        return(0);
    if(n_samp == 0 || n_site == 0)
        return(0);
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
                fprintf(file_output,"\nError: memory not reallocated. get_obsdata.23d \n");
                for(x=0;x<n_samp;x++) free(names[x]); free(names);
                free(DNA_matr);
                return(0);
            }
            if((names2 = (char **)calloc(n_samp,sizeof(char *))) == 0) {
                fprintf(file_output,"\nError: memory not reallocated. get_obsdata.1s2 \n");
                for(x=0;x<n_samp;x++) free(names[x]); free(names);
                free(DNA_matr);
                return(0);
            }
            for(x=0;x<n_samp;x++) {
                if((names2[x] = (char *)calloc(50,sizeof(char))) == 0) {
                    fprintf(file_output,"\nError: memory not reallocated. get_obsdata.22 \n");
                    for(x=0;x<n_samp;x++) free(names[x]); free(names);
                    free(DNA_matr);
                    return(0);
                }
                /*copy duplicated data*/
                strncpy(names2[x],names[x],50);
            }
            strncpy(DNA_matr2+(long long)n_site*(long long)0,DNA_matr+(long long)n_site*(long long)0,(long long)n_site*n_samp);
            /*end define and duplicate*/
            
            /*include data in *DNA_matr and in *names[] in the correct order*/
            for(x=0;x<nsamuser_eff;x++) {
                strncpy(DNA_matr+(long long)n_site*(long long)x,DNA_matr2+(long long)n_site*(long long)sort_nsam[x],n_site);
                strncpy(names[x],names2[sort_nsam[x]],50);
            }
            /*delete duplicated matr*/
            for(x=0;x<n_samp;x++) free(names2[x]);
            free(names2);
            free(DNA_matr2);
            
            /*erase lines no used*/
            if(nsamuser_eff > n_samp) {
                fprintf(file_output,"Error: too low samples in the file according to defined in -N flag.\n");
                for(x=0;x<n_sam;x++) free(names[x]); free(names);
                free(DNA_matr);
                return(0);
            }
        }
        /*end option flag O*/
        if(nsamuser_eff > 32167) {
            fprintf(file_output,"Error: too much samples. Only 32167 samples per loci are allowed.\n");
            return(0);
        }
        /*init matrix_sizepos*/
        if(matrix_sizepos == 0) {
            if((matrix_sizepos = (double *)malloc(n_site*sizeof(double))) == 0) {
                fprintf(file_output,"Error: memory not reallocated. get_obsstat.2");
                for(x=0;x<n_samp;x++) free(names[x]); free(names);
                free(DNA_matr);
                return(0);
            }
            if((matrix_segrpos = (double *)malloc(n_site*sizeof(double))) == 0) {
                fprintf(file_output,"Error: memory not reallocated. get_obsstat.2");
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
                    fprintf(file_output,"Error: memory not reallocated. get_obsstat.2b");
                    for(x=0;x<n_samp;x++) free(names[x]); free(names);
                    free(DNA_matr);
                    return(0);
                }
                if((matrix_segrpos = (double *)realloc(matrix_segrpos,n_site*sizeof(double))) == 0) {
                    fprintf(file_output,"Error: memory not reallocated. get_obsstat.2b");
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
        if(outgroup_presence == 0) {
            if ((DNA_matr = (char *)realloc(DNA_matr,(long long)n_site*(nsamuser_eff+!outgroup_presence)*sizeof(char))) == 0) {
                fprintf(file_output,"\nError: memory not reallocated. get_obsdata.23a \n");
                for(x=0;x<n_samp;x++) free(names[x]); free(names);
                free(DNA_matr);
                free(matrix_sizepos);
                free(matrix_segrpos);
                return(0);
            }
            /*we forced the invented outgroup without gaps or uncertainties, if possible*/
            strncpy(DNA_matr+(unsigned long)n_site*(unsigned long)(nsamuser_eff),DNA_matr+(unsigned long)n_site*(unsigned long)(nsamuser_eff-1),n_site);
            nsamuser_eff += 1;
        }
        
        if(wP!=0) {
            /*define the weights*/
            for(n_sit=0;n_sit<n_site;n_sit++) {
                matrix_sizepos[n_sit] =  wP[beg+n_sit-1];
                matrix_segrpos[n_sit] = wPV[beg+n_sit-1]/* * wV[nsit-1]*/;
            }
        }
        /*define variables for mhits*/
        *nmhits = 0;
        if((mhitbp = (long int *) calloc (n_site, sizeof(long int))) == 0) {
            fprintf(file_output,"Error: memory not reallocated. get_obsstat.6");
            return(0);
        }
        /*function to analyze all data*/
        if(get_obsstats(file_output,'\0',nsamuser_eff,n_site,length_al_real,names/*2*/,DNA_matr/*2*/,matrix_sizepos,matrix_segrpos,
                        matrix_pol,matrix_freq,matrix_pos,length_al,length_seg,nsamuser,npops,svratio,missratio,include_unknown,
                        sum_sam,tcga,matrix_sv,nmhits,output,ploidy,outgroup_presence,nsites1_pop,nsites1_pop_outg,
                        nsites2_pop,nsites2_pop_outg,nsites3_pop,nsites3_pop_outg,anx,bnx,anxo,bnxo,lengthamng,length_amng_outg,mhitbp,matrix_pol_tcga,(long int)beg) == 0) {
            for(x=0;x<n_sam;x++) free(names[x]); free(names);
            free(DNA_matr);
            /*free(DNA_matr2);*/
            free(matrix_sizepos);
            free(matrix_segrpos);
            free(mhitbp);
            return(0);
        }
		/*free(names2);
		free(DNA_matr2);*/
		free(matrix_sizepos);
		free(matrix_segrpos);
        free(mhitbp);
    }
    return(1);
}

int read_coordinates(FILE *file_wcoor, FILE *file_output,long int **wgenes, long int *nwindows) {
    
    char *valn=0;
    int c,x;
    long int xx;
    long int dd;
    double ee;
    long int prevwin = 0;
    
    /*
    printf("\nReading coordinates file...");
    fflush(stdout);
    *//*
    fprintf(file_output,"\nReading coordinates file...");
    fflush(file_output);
    */
    if((valn = (char *)calloc(100,sizeof(char))) == 0) {
        fprintf(file_output,"\nError: memory not reallocated. read_coordinates.00 \n");
        return(0);
    }
    if((*wgenes = (long int *)calloc(1000,sizeof(long int))) == 0) {
        fprintf(file_output,"\nError: memory not reallocated. read_coordinates.0 \n");
        return(0);
    }
    
    *nwindows = 0;
    c = fgetc(file_wcoor);
    if(c=='#') {
        while(c != 13 && c != 10 && c != 0 && c != -1)
            c = fgetc(file_wcoor); /*exclude header*/
    }
    while(c == 13 || c == 10)
        c = fgetc(file_wcoor);
    if(c==0 || c==-1) {
        fprintf(file_output,"\nError: no coordinates assigned, read_coordinates.2 \n");
        free(*wgenes);
        file_wcoor=0;
        return(0);
    }
    else {
        xx=0;
        while (c != 0 && c != -1) {
            /*now keep all values: three columns, name and two columns of numbers*/
            /*first column is the name of the scaffold*/
            while(c == 32 || c == 9 || c == 13 || c == 10)
                c = fgetc(file_wcoor);
            if(c==-1)
                break;
            x=0;
            while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100) {
                valn[x] = c;
                c = fgetc(file_wcoor);
                x++;
            }
            valn[x] = '\0';/*scaffold name*/
            /*KEEP POSITIONS (first initial position, then end, next region and so on)*/
            while(c == 32 || c == 9 || c == 13 || c == 10)
                c = fgetc(file_wcoor);
            if(c==-1)
                break;
            x=0;
            while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100) {
                valn[x] = c;
                c = fgetc(file_wcoor);
                x++;
            }
            valn[x] = '\0';
            wgenes[0][xx] = (long int)round((double)atof(valn));
            if(wgenes[0][xx] <= prevwin) {
                printf("\nError: file with coordinates has overlapped or unsorted window positions: %ld.\n",wgenes[0][xx]);
                exit(1);
            }
            prevwin = wgenes[0][xx];
            
            xx++;
            while(c == 32 || c == 9 || c == 13 || c == 10)
                c = fgetc(file_wcoor);
            if(c==-1)
                break;
            x=0;
            while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100) {
                valn[x] = c;
                c = fgetc(file_wcoor);
                x++;
            }
            valn[x] = '\0';
            wgenes[0][xx] = (long int)round((double)atof(valn));
            if(wgenes[0][xx] <= prevwin) {
                printf("\nError: file with coordinates has overlapped or unsorted window positions: %ld.\n",wgenes[0][xx]);
                exit(1);
            }
            prevwin = wgenes[0][xx];

            xx++;
            dd = (long int)floor((double)xx/(double)1000);
            ee = (double)xx/(double)1000;
            if(dd==ee) {
                if((*wgenes = realloc(*wgenes,((long int)(dd+1)*(long int)1000*sizeof(long int)))) == 0) {
                    puts("Error: realloc error read_coordinates.1\n");
                    free(*wgenes);
                    return(0);
                }    
            }
        }
        if(xx == 0) {
            fprintf(file_output,"\nError: no coordinates assigned, read_coordinates.2 \n");
            free(*wgenes);
            file_wcoor=0;
            return(0);
        }
        *nwindows = (xx)/2;
    }		
    free(valn);
    return 1;
}

int read_weights_positions_file(FILE *file_ws, SGZip *file_ws_gz, FILE *file_output, float **wP, float **wPV, float **wV, long int *wlimit_end) {
    
    long int position;
    char *valn=0;
    long int dd;
    double ee;
    int c;
    int x;
    long int xx;
    
    /*read weights file*/
    /*fprintf(file_output,"\nReading weight file...");*/
    fflush(file_output);
    /*
    printf("\nReading weight file...");
    fflush(stdout);
    */
    
    /*keep wP, wPV, wV (not used yet) for all positions, do not need to keep positions: all are correlative*/
    if(file_ws != 0) {
        if((*wP = (float *)calloc(1000,sizeof(float))) == 0) {
            fprintf(file_output,"\nError: memory not reallocated. get_obsdata.14 \n");
            return(0);
        }
        if((*wPV = (float *)calloc(1000,sizeof(float))) == 0) {
            fprintf(file_output,"\nError: memory not reallocated. get_obsdata.25 \n");
            return(0);
        }
        if((*wV = (float *)calloc(1000,sizeof(float))) == 0) {
            fprintf(file_output,"\nError: memory not reallocated. get_obsdata.23 \n");
            return(0);
        }
        if((valn = (char *)calloc(100,sizeof(char))) == 0) {
            fprintf(file_output,"\nError: memory not reallocated. get_obsdata.34 \n");
            return(0);
        }
        c = fzgetc(file_ws,file_ws_gz);
        if(c=='#') {
            while(c != 13 && c != 10 && c != 0 && c != -1)
                c = fzgetc(file_ws,file_ws_gz); /*exclude header*/
        }
        while(c == 13 || c == 10)
            c = fzgetc(file_ws,file_ws_gz);
        if(c==0 || c==-1) {
            file_ws=0;
            *wP=0;
            *wPV=0;
            *wV=0;
            free(*wP);
            free(*wPV);
            free(*wV);
            fprintf(file_output,"\nError: no weights assigned \n");
            return(0);
        }
        else {
            /*now keep all values: three or four columns, only numbers or decimals*/
            xx=0;
            while (c != 0 && c != -1) {
                /*POSITION*/
                while(c == 32 || c == 9 || c == 13 || c == 10)
                    c = fzgetc(file_ws,file_ws_gz);
                if(c==0 || c==-1)
                    break;
                x=0;
                while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100) {
                    valn[x] = c;
                    c = fzgetc(file_ws,file_ws_gz);
                    x++;
                }
                valn[x] = '\0';
                position = atol(valn);
                if(xx < position-1) {
                    xx = position-1;
                    if((*wP = realloc(*wP,((long int)(xx+1000)*sizeof(float)))) == 0) {
                        file_ws=0;*wP=0;*wPV=0;*wV=0;
                        free(*wP);free(*wPV);free(*wV);
                        puts("Error: realloc error get_obsdata.11\n");
                        return(0);
                    }
                    if((*wPV = realloc(*wPV,((long int)(xx+1000)*sizeof(float)))) == 0) {
                        file_ws=0;*wP=0;*wPV=0;*wV=0;
                        free(*wP);free(*wPV);free(*wV);
                        puts("Error: realloc error get_obsdata.11\n");
                        return(0);
                    }
                    if((*wV = realloc(*wV,((long int)(xx+1000)*sizeof(float)))) == 0) {
                        file_ws=0;*wP=0;*wPV=0;*wV=0;
                        free(*wP);free(*wPV);free(*wV);
                        puts("Error: realloc error get_obsdata.11\n");
                        return(0);
                    }
                }
                
                /*Weight Position*/
                while(c == 32 || c == 9)/* || c == 13 || c == 10)*/
                    c = fzgetc(file_ws,file_ws_gz);
                if(c==0 || c==-1)
                    break;
                x=0;
                while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100) {
                    valn[x] = c;
                    c = fzgetc(file_ws,file_ws_gz);
                    x++;
                }
                valn[x] = '\0';
                wP[0][xx] = (double)atof(valn);
                
                /*Weight Variant*/
                while(c == 32 || c == 9)/* || c == 13 || c == 10)*/
                    c = fzgetc(file_ws,file_ws_gz);
                if(c==0 || c==-1)
                    break;
                x=0;
                while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100) {
                    valn[x] = c;
                    c = fzgetc(file_ws,file_ws_gz);
                    x++;
                }
                valn[x] = '\0';
                wPV[0][xx] = (double)atof(valn);
                
                while(c == 32 || c == 9) c = fzgetc(file_ws,file_ws_gz);
                if(!(c == 13 || c == 10 || c == 0 || c == -1)) {
                    /*Effect size*/
                    while(c == 32 || c == 9)/* || c == 13 || c == 10)*/
                        c = fzgetc(file_ws,file_ws_gz);
                    if(c==0 || c==-1)
                        break;
                    x=0;
                    while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100) {
                        valn[x] = c;
                        c = fzgetc(file_ws,file_ws_gz);
                        x++;
                    }
                    valn[x] = '\0';
                    wV[0][xx] = (double)atof(valn);
                }
                else {
                    wV[0][xx] = 1.0; /*if undefined the value is 1.0 for all*/
                }
                
                xx++;
                dd = (long int)floor((double)xx/(double)1000);
                ee = (double)xx/(double)1000;
                if(dd==ee) {
                    if((*wP = realloc(*wP,((long int)(dd+1)*(long int)1000*sizeof(float)))) == 0) {
                        file_ws=0;*wP=0;*wPV=0;*wV=0;
                        free(*wP);free(*wPV);free(*wV);
                        puts("Error: realloc error get_obsdata.11\n");
                        return(0);
                    }
                    if((*wPV = realloc(*wPV,((long int)(dd+1)*(long int)1000*sizeof(float)))) == 0) {
                        file_ws=0;*wP=0;*wPV=0;*wV=0;
                        free(*wP);free(*wPV);free(*wV);
                        puts("Error: realloc error get_obsdata.11\n");
                        return(0);
                    }
                    if((*wV = realloc(*wV,((long int)(dd+1)*(long int)1000*sizeof(float)))) == 0) {
                        file_ws=0;*wP=0;*wPV=0;*wV=0;
                        free(*wP);free(*wPV);free(*wV);
                        puts("Error: realloc error get_obsdata.11\n");
                        return(0);
                    }
                }
            }
        }
        free(valn);
        *wlimit_end = xx;
    }
    else return(0);
    
    return(1);
}

int read_weights_file(FILE *file_es, FILE *file_output, float **wV, long int **Pp, long int *nV, long int *welimit_end) {
    
    long int *effsz_site=0; /*positions*/
    float    *effsz_wght=0; /*effect sizes*/
    long int tot_effszP = 0; /*total values*/
    char *valn=0;
    long int dd;
    double ee;
    int c;
    int x;
    long int xx;
    
    /*printf("\nReading Effect sizes file...");*/
    fflush(stdout);
    
    if(file_es != 0) {
        if((effsz_site = (long int *)calloc(1000,sizeof(long int))) == 0) {
            fprintf(file_output,"\nError: memory not reallocated. get_obsdata.1 \n");
            return(0);
        }
        if((effsz_wght = (float *)calloc(1000,sizeof(float))) == 0) {
            fprintf(file_output,"\nError: memory not reallocated. get_obsdata.2 \n");
            return(0);
        }
        if((valn = (char *)calloc(100,sizeof(char))) == 0) {
            fprintf(file_output,"\nError: memory not reallocated. get_obsdata.3 \n");
            return(0);
        }
        
        c = fgetc(file_es);
        while(c != 13 && c != 10 && c != 0 && c != -1)
            c = fgetc(file_es); /*exclude header*/
        while(c == 13 || c == 10)
            c = fgetc(file_es);
        if(c==0 || c==-1) {
            file_es=0;
            *wV=0;
            free(effsz_site);
            free(effsz_wght);
            fprintf(file_output,"\nError: no effect sizes assigned \n");
            return(0);
        }
        else {
            /*now keep all values: two columns, only numbers or decimals*/
            xx=0;
            while (c != 0 && c != -1) {
                /*POSITION*/
                while(c == 32 || c == 9 || c == 13 || c == 10)
                    c = fgetc(file_es);
                if(c==-1)
                    break;
                x=0;
                while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100) {
                    valn[x] = c;
                    c = fgetc(file_es);
                    x++;
                }
                valn[x] = '\0';
                effsz_site[xx] = (double)atof(valn);
                
                /*Effect size*/
                while(c == 32 || c == 9 || c == 13 || c == 10)
                    c = fgetc(file_es);
                if(c==-1)
                    break;
                x=0;
                while(c != 32 && c != 9 && c != 13 && c != 10 && c!=0 && c!=-1 && x < 100) {
                    valn[x] = c;
                    c = fgetc(file_es);
                    x++;
                }
                valn[x] = '\0';
                effsz_wght[xx] = (double)atof(valn);
                
                xx++;
                dd = (long int)floor((double)xx/(double)1000);
                ee = (double)xx/(double)1000;
                if(dd==ee) {
                    if((effsz_site = (long int *)realloc((long int *)effsz_site,((long int)(dd+1)*(long int)1000*sizeof(long int)))) == 0) {
                        file_es=0;*wV=0;
                        free(effsz_site);free(effsz_wght);
                        puts("Error: realloc error get_obsdata.11\n");
                        return(0);
                    }
                    if((effsz_wght = (float *)realloc((float *)effsz_wght,((long int)(dd+1)*(long int)1000*sizeof(float)))) == 0) {
                        file_es=0;*wV=0;
                        free(effsz_site);free(effsz_wght);
                        puts("Error: realloc errorget_obsdatavarchar.12\n");
                        return(0);
                    }
                }
            }
            if(effsz_site[xx]== 0) xx--;
            tot_effszP = xx;/*total number of values added*/
            
            /*now we need to assign this values to the variant positions here observed and include them in wV vector*/
            if(xx>0) {
                if((*Pp = (long int *)calloc(xx,sizeof(long int))) == 0) {
                    file_es=0;*wV=0;
                    free(effsz_site);free(effsz_wght);
                    fprintf(file_output,"\nError: memory not reallocated. get_obsdata.3 \n");
                    return(0);
                }
                if((*wV = (float *)calloc(xx,sizeof(float))) == 0) {
                    file_es=0;*wV=0;
                    free(effsz_site);free(effsz_wght);
                    fprintf(file_output,"\nError: memory not reallocated. get_obsdata.2 \n");
                    return(0);
                }
                for( x = 0; x < (int)xx; x++ ) {
                    Pp[0][x] = effsz_site[x];
                    wV[0][x] = effsz_wght[x];
                }
                *nV = tot_effszP;
            }
            free(effsz_site);
            free(effsz_wght);
        }
        free(valn);
        *welimit_end = xx;
    }
    /*end collecting data for effect sizes*/
    return 1;
}

int function_read_tfasta(FILE *file_input,SGZip *input_gz,long int init_site,long int end_site,int *n_sam, long int *n_site, char ***names, char **DNA_matr,char **matrix_pol_tcga)
{
    char c[1];
    char *cc;
    char line[32767*2];
    int col=0;
    static int nseq=0;
    int x;
    static int maxsam=128;
    static long int position=0;
    long int end_position=end_site;
    long int dd,ee,count,xx;
    char *DNA_matr2;
    
    c[0] = '\t'; //! added
    if ((DNA_matr2 = (char *)calloc(10000,sizeof(char))) == 0) {
        printf("\nError: memory not reallocated. get_tfadata.23d \n");
        free(DNA_matr2);
        return(0);
    }

    if(position == 0) {/*to allow sliding windows non-overlapped*/
        /*if names and number of samples are defined then skip the definition of names*/
        *c = fzgetc(file_input, input_gz);
        if((*c==10 || *c==13 || *c == -1 || *c == 0 || *c=='\xff' || *c=='\xfe'))
            return -1;
        line[col] = *c;
        col++;
        
        while(*c=='#') {
            /*parse comments*/
            while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c!='\xff' && *c!='\xfe') {
                *c = fzgetc(file_input, input_gz);
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
                                puts("\n Sorry, no more than 32767 samples are allowed.");
                                return 0;
                            }
                            if ((*names = (char **)realloc(*names,maxsam*sizeof(char *))) == 0) {
                                puts("\nError: memory not reallocated. assigna.1 \n");
                                return(0);
                            }
                            for(x=nseq;x<maxsam;x++) {
                                if ((names[0][x] = (char *)calloc(50,sizeof(char))) == 0) {
                                    puts("\nError: memory not reallocated. assigna.2 \n");
                                    return(0);
                                }
                            }
                        }
                    }
                    cc = strtok(NULL, ">\n\r ");
                }
                n_sam[0] = nseq;
            }
            col=0;
            *c = fzgetc(file_input, input_gz);
            line[col] = *c;
            col++;
        }
        
        /*include in DNAmatrix the positions from the init_site to the end_site (if defined)*/
        col = 0;
        while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c != 9 && *c != 32 && *c!='\xff' && *c!='\xfe') {
            line[col] = *c;
            col++;
            *c = fzgetc(file_input, input_gz);
        }
        if(*c == 0 || *c==-1 || *c=='\xff' || *c=='\xfe') return(-1);
        
        line[col] = '\0';
        position = atol(line);
    }
    if(end_site == -1) end_position = init_site + 1;
    while(position < init_site) {
        while(!(*c==10 || *c==13 || *c == -1 || *c == 0 || *c=='\xff' || *c=='\xfe')) *c = fzgetc(file_input, input_gz);
        if(*c == 0 || *c==-1 || *c=='\xff' || *c=='\xfe') return(-1);
        if(*c==10 || *c==13) {
            while(*c==10 || *c==13) *c = fzgetc(file_input, input_gz);
            col = 0;
            while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c != 9 && *c != 32 && *c!='\xff' && *c!='\xfe') {
                line[col] = *c;
                col++;
                *c = fzgetc(file_input, input_gz);
            }
            line[col] = '\0';
            position = atol(line);
            if(end_site == -1) end_position = position + 1;
        }
    }
    if(!(*c == 0 || *c==-1 || *c=='\xff' || *c=='\xfe'))
        *c = fzgetc(file_input, input_gz);

    *n_site = 0;
    col = 0;
    count=0;
    while (position <= end_position) {
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
                printf("Unexpected value in tfa file: position %ld, sample %d: ",position,col);
                printf("%c\n",*c);
                return(0);
                break;
        }
        if(*c == -1 || *c == 0) break;
        *c = fzgetc(file_input, input_gz);
        if(*c==10 || *c==13 || *c == -1 || *c == 0) {
            while(*c==10 || *c==13) *c = fzgetc(file_input, input_gz);
            if(col != nseq) {
                printf("\nError: The number of columns are not coincident in all positions: %ld.\n",position);
            }
            n_sam[0] = col;
            col = 0;
            while(*c != 0 && *c!=-1 && *c!=10 && *c!=13 && *c != 9 && *c != 32) {
                line[col] = *c;
                col++;
                *c = fzgetc(file_input, input_gz);
            }
            line[col] = '\0';
            col = 0;
            position = atol(line);
            *n_site += 1;
            if(end_site == -1) end_position = position + 1;
        }
        /*realloc DNAmatr if nnecessary*/
        dd = (long int)floor((double)count/(double)10000);
        ee = (double)count/(double)10000;
        if(dd == ee) {
            if((DNA_matr2 = (char *)realloc((char *) DNA_matr2,((long int)dd+(long int)1)*(long int)10000*sizeof(char))) == 0) {
                puts("Error: realloc error varchar.1\n");
                return(0);
            }
        }
    }
    /*return n_site,n_sam,DNA matrix and names in pointers*/
    /*
    if ((*DNA_matr = (char *)realloc((char *) *DNA_matr,n_site[0]*(long long)n_sam[0]*sizeof(char))) == 0) {
        printf("\nError: memory not reallocated. get_tfadata.23d \n");
        free(DNA_matr2);
        return(0);
    }
    */
    free(*DNA_matr);
    if((*DNA_matr = (char *)calloc(n_site[0]*(long long)n_sam[0],sizeof(char))) == 0) {
        printf("\nError: memory not reallocated. get_tfadata.23d2 \n"); free(DNA_matr2);
        return(0);
    }
    for(x=0;x<n_sam[0];x++) {
        for(xx=0;xx<n_site[0];xx++) { /*transpose */
            DNA_matr[0][(((long long)n_site[0]*(unsigned long)x)+(unsigned long)xx)] =
            DNA_matr2 [(((long long)n_sam[0]*(unsigned long)xx)+(unsigned long)x)];
        }
    }
    free(DNA_matr2);
    
    return(1);
}
