/*
 *  get_obsdata.c
 *  MuLoNeTests
 *
 *  Created by sebas on Sun Feb 23 2003.
 *
 */

#include "get_obsdata.h"

int get_obsdata(FILE *file_output,SGZip *file_output_gz,
				FILE *file_input, SGZip *input_gz,
                FILE *file_logerr,SGZip *file_logerr_gz,
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
				double **lengthamng, double **lenghtamng_outg,int *sort_nsam,char **matrix_pol_tcga,
                char *chr_name,int first)
{    
    static char *DNA_matr = 0;
    static char *DNA_matr2 = 0;
    static char **names = 0; /* limit of the name of the lines to 50 characters. be careful ... */
    static char **names2 = 0; /* limit of the name of the lines to 50 characters. be careful ... */
	static double *matrix_sizepos = 0; /*size of each position, 0 not include, >0 include, Syn/Nsyn positions are allowed*/
	static double *matrix_segrpos = 0; /*always 1 except in those biallelic positions that are not desired (e.g., choose syn but we do not want nsyn)*/

	long int *mhitbp; /*vector of positions where mhits are*/
    long int count,xx;
    int c;
    int n_sam,ns;
    long int n_sit;
    int nseq;
    int maxsam;
	int n_excl;
    int n_samp;
    long int n_site;
        
    int x;
	static long int maxsites = 0;
	
	/*all these variables are zero because are not used here*/
	int excludelines = 0;
	char *name_excluded = 0;
	int includelines = 0;
	char *name_ingroups = 0;
	char *name_outgroup = 0;
	int outgroup = 0;
	int nsamtot;
	int nsamuser_eff;
	int flag_change_sort; /*in case the order of samples is not consecutive*/
		
    count = 0;
    c = 0;
    n_sam = 0;
    n_sit = 0;
    nseq  = 0;
    maxsam= 128;
    n_samp= 0;
    n_site= 0;
	n_excl= 0;
	
	nsamtot = 0;
	for(x=0;x<npops;x++) {
		nsamtot += nsamuser[x];
	}
	
	nsamuser_eff = (nsamtot- !outgroup_presence)/atoi(ploidy) ;
	
    if(names == 0) { /* only initialize once. Check */
        if((names = (char **)calloc(128,sizeof(char *))) == 0) {
            fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. get_obsdata.1 \n");
            return(0);
        }
        for(x=0;x<128;x++) {
            if((names[x] = (char *)calloc(50,sizeof(char))) == 0) {
                fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. get_obsdata.2 \n");
                return(0);
            }
        }    
        if((DNA_matr = (char *)calloc(10000,sizeof(char))) == 0) {
			for(x=0;x<128;x++) free(names[x]); free(names);
            fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. get_obsdata.3 \n");
            return(0);
        }
    }
    
    c = fzgetc(file_input, input_gz);
    while (c != 0 && c != -1 /*&& n_sam < nsamuser_eff*/) {
        while(c == 32 || c == 9 || c == 13 || c == 10 || c == '*') c = fzgetc(file_input, input_gz);
        n_sit = 0;
        if(!(var_char(file_input,input_gz,file_logerr,file_logerr_gz,&count,&c,&n_sam,&n_sit,&nseq,&maxsam,&names,&DNA_matr,&n_site,excludelines,name_excluded,&n_excl,includelines,name_ingroups,name_outgroup,outgroup,nsamuser_eff,ploidy)))
            return 0;
        if(n_sam == 32767) {
            fzprintf(file_logerr,file_logerr_gz,"Only 32167 samples per loci are allowed.\n");
            break;
        }
    }
    n_samp = n_sam;
	*length = n_site;
	
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
			fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. get_obsdata.23d \n");
			for(x=0;x<n_samp;x++) free(names[x]); free(names);
			free(DNA_matr);
			return(0);
		}
		if((names2 = (char **)calloc(n_samp,sizeof(char *))) == 0) {
			fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. get_obsdata.1s2 \n");
			for(x=0;x<n_samp;x++) free(names[x]); free(names);
			free(DNA_matr);
			return(0);
		}
		for(x=0;x<n_samp;x++) {
			if((names2[x] = (char *)calloc(50,sizeof(char))) == 0) {
				fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. get_obsdata.22 \n");
				for(x=0;x<n_samp;x++) free(names[x]); free(names);
				free(DNA_matr);
				return(0);
			}
			/*copy duplicated data*/
			strncpy(DNA_matr2+(long long)n_site*(long long)x,DNA_matr+(long long)n_site*(long long)x,n_site);			
			strncpy(names2[x],names[x],50);
		}
		/*end define and duplicate*/
		
		/*include data in *DNA_matr and in *names[] in the correct order*/
		for(x=0;x<nsamuser_eff;x++) {
			strncpy(DNA_matr+(long long)n_site*(long long)x,DNA_matr2+(long long)n_site*(long long)sort_nsam[x],n_site);			
			strncpy(names[x],names2[sort_nsam[x]],50);
		}
		/*delete duplicated matr*/
		for(x=0;x<n_samp;x++) free(names2[x]); free(names2);
		free(DNA_matr2);
		
		/*erase lines no used*/
		if(nsamuser_eff > n_samp) {
			fzprintf(file_logerr,file_logerr_gz,"Error: too low samples in the file according to defined in -N flag.\n");
			for(x=0;x<n_sam;x++) free(names[x]); free(names);
			free(DNA_matr);
			return(0);
		}
	}
	/*end option flag O*/
    
	if(n_samp * atoi(ploidy) < nsamtot-!outgroup_presence) return(0);
    if(n_samp == 0 || n_site == 0) return(0);
    else {
        if(nsamuser_eff > 32767) {
            fzprintf(file_logerr,file_logerr_gz,"Error: too much samples. Only 32767 samples per loci are allowed.\n");
            return(0);
        }
		/*init matrix_sizepos*/
		if(matrix_sizepos == 0) {
			if((matrix_sizepos = (double *)malloc(n_site*sizeof(double))) == 0) {
				fzprintf(file_logerr,file_logerr_gz,"Error: memory not reallocated. get_obsstat.2"); 
				for(x=0;x<n_sam;x++) free(names[x]); free(names);
				free(DNA_matr);
				return(0);
			}
			if((matrix_segrpos = (double *)malloc(n_site*sizeof(double))) == 0) {
				fzprintf(file_logerr,file_logerr_gz,"Error: memory not reallocated. get_obsstat.2"); 
				for(x=0;x<n_sam;x++) free(names[x]); free(names);
				free(DNA_matr);
				free(matrix_sizepos);
				return(0);
			}
			maxsites = n_site;
		}
		else{
			if(n_site > maxsites) {
				if((matrix_sizepos = (double *)realloc(matrix_sizepos,n_site*sizeof(double))) == 0) {
					fzprintf(file_logerr,file_logerr_gz,"Error: memory not reallocated. get_obsstat.2b"); 
					for(x=0;x<n_sam;x++) free(names[x]); free(names);
					free(DNA_matr);
					return(0);
				}
				if((matrix_segrpos = (double *)realloc(matrix_segrpos,n_site*sizeof(double))) == 0) {
					fzprintf(file_logerr,file_logerr_gz,"Error: memory not reallocated. get_obsstat.2b"); 
					for(x=0;x<n_sam;x++) free(names[x]); free(names);
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
		if(ploidy[0]=='2'){
			if ((DNA_matr2 = (char *)calloc(n_site*(long long)(nsamuser_eff)*2,sizeof(char))) == 0) {
				fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. get_obsdata.23 \n");
				for(x=0;x<n_sam;x++) free(names[x]); free(names);
				free(DNA_matr);
				free(matrix_sizepos);
				free(matrix_segrpos);
				return(0);
			}
			if((names2 = (char **)calloc((nsamuser_eff)*2,sizeof(char *))) == 0) {
				fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. get_obsdata.12 \n");
				return(0);
			}
			for(x=0;x<(nsamuser_eff)*2;x++) {
				if((names2[x] = (char *)calloc(50,sizeof(char))) == 0) {
					fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. get_obsdata.22 \n");
					return(0);
				}
			}
			for(x=0;x<nsamuser_eff;x++) {
				strncpy(DNA_matr2+(long long)n_site*(long long)x*2,DNA_matr+(long long)n_site*(long long)x,n_site);
				strncpy(DNA_matr2+(long long)n_site*(long long)(x*2+1),DNA_matr+(long long)n_site*(long long)x,n_site);	
				
				strncpy(names2[x*2],names[x],50);
				strncpy(names2[x*2+1],names[x],50);
				strncat(names2[x*2],".1\0",50);
				strncat(names2[x*2+1],".2\0",50);
			}
			for(xx=0;xx<n_site;xx++) {
				for(x=0;x<nsamuser_eff;x++) {
					if(DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] == 'w') {
						DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] = '1';
						DNA_matr2[(long long)n_site*(unsigned long)(x*2+1) + xx] = '4';
					}
					if(DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] == 'm') {
						DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] = '2';
						DNA_matr2[(long long)n_site*(unsigned long)(x*2+1) + xx] = '4';
					}
					if(DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] == 'r') {
						DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] = '3';
						DNA_matr2[(long long)n_site*(unsigned long)(x*2+1) + xx] = '4';
					}
					if(DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] == 'y') {
						DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] = '1';
						DNA_matr2[(long long)n_site*(unsigned long)(x*2+1) + xx] = '2';
					}
					if(DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] == 'k') {
						DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] = '1';
						DNA_matr2[(long long)n_site*(unsigned long)(x*2+1) + xx] = '3';
					}
					if(DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] == 's') {
						DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] = '2';
						DNA_matr2[(long long)n_site*(unsigned long)(x*2+1) + xx] = '3';
					}
					if(DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] == 't') {
						DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] = '1';
						DNA_matr2[(long long)n_site*(unsigned long)(x*2+1) + xx] = '5';
					}
					if(DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] == 'c') {
						DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] = '2';
						DNA_matr2[(long long)n_site*(unsigned long)(x*2+1) + xx] = '5';
					}
					if(DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] == 'g') {
						DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] = '3';
						DNA_matr2[(long long)n_site*(unsigned long)(x*2+1) + xx] = '5';
					}
					if(DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] == 'a') {
						DNA_matr2[(long long)n_site*(unsigned long)x*2 + xx] = '4';
						DNA_matr2[(long long)n_site*(unsigned long)(x*2+1) + xx] = '5';
					}
				}
			}
			nsamuser_eff *= 2;
		}
		else {
			if((DNA_matr2 = (char *)calloc((long long)n_site*(nsamuser_eff),sizeof(char))) == 0) {
				fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. get_obsdata.23 \n");
				for(x=0;x<n_sam;x++) free(names[x]); free(names);
				free(DNA_matr);
				free(matrix_sizepos);
				free(matrix_segrpos);
				return(0);
			}
			if((names2 = (char **)calloc((nsamuser_eff)*1,sizeof(char *))) == 0) {
				fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. get_obsdata.12 \n");
				return(0);
			}
			for(x=0;x<(nsamuser_eff)*1;x++) {
				if((names2[x] = (char *)calloc(50,sizeof(char))) == 0) {
					fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. get_obsdata.22 \n");
					return(0);
				}
			}
			for(x=0;x<(nsamuser_eff)*1;x++) {
				strncpy(DNA_matr2+(long long)n_site*(long long)x,DNA_matr+(long long)n_site*(long long)x,n_site);				
				strncpy(names2[x],names[x],50);
			}
		}
		if(outgroup_presence == 0) {
			if ((DNA_matr2 = (char *)realloc(DNA_matr2,(long long)n_site*(nsamuser_eff+!outgroup_presence)*sizeof(char))) == 0) {
				fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. get_obsdata.23a \n");
				for(x=0;x<n_sam;x++) free(names[x]); free(names);
				free(DNA_matr);
				free(DNA_matr2);
				free(matrix_sizepos);
				free(matrix_segrpos);
				return(0);
			}
			/*we forced the invented outgroup without gaps or uncertainties, if possible*/
			/*strncpy(DNA_matr2+(unsigned long)n_site*(unsigned long)(n_samp),DNA_matr2+(unsigned long)n_site*(unsigned long)(n_samp-1),n_site);*/
			for(xx=0;xx<n_site;xx++) {
				ns = 0;
				while(ns < nsamuser_eff-1 && DNA_matr2[(long long)n_site*(unsigned long)ns+xx] > '4') ns++;
				DNA_matr2[(unsigned long)n_site*(long long)(nsamuser_eff)+xx] = DNA_matr2[(unsigned long)n_site*(unsigned long)ns+xx];
			}
            nsamuser_eff += 1;
		}
		
		/*define variables for mhits*/
		*nmhits = 0;
		if((mhitbp = (long int *) calloc (n_site, sizeof(long int))) == 0) {
			fzprintf(file_logerr,file_logerr_gz,"Error: memory not reallocated. get_obsstat.6"); 
			return(0);
		}
		/*here include a function to filter positions: to read gff files (if necessary)*/
		if(gfffiles == 1) {
			/*include**name_fileinputgff,*subset_positions,ifgencode,*codename,*genetic_code,*matrix_sizepos,n_samp,n_site,*DNA_matr*/
			/*the function read the gff file and cut the DNA_matr, also gives the number of positions in matrix_sizepos and count the total in n_site*/
			/*only modify values in matrix_sizepos*/
			if(use_gff(name_fileinputgff,subset_positions,genetic_code,matrix_sizepos,nsamuser_eff,n_site,DNA_matr2,matrix_segrpos,
					   file_output,file_output_gz,mainargc,file_logerr,file_logerr_gz,include_unknown,criteria_transcript,output,nmhits,mhitbp,outgroup_presence,nsamuser[npops-1],chr_name,first) == 0) {
				/*if error realloc DNA_matr*/
				for(x=0;x<n_sam;x++) free(names[x]); free(names);
				free(DNA_matr);
				free(DNA_matr2);
				free(matrix_sizepos);
				free(matrix_segrpos);
				free(mhitbp);
				return(0);
			}
		}
        /*function to analyze all data*/
        if(get_obsstats(file_output,file_output_gz,file_mask,file_logerr,file_logerr_gz,nsamuser_eff,n_site,length_al_real,names2,DNA_matr2,matrix_sizepos,matrix_segrpos,
						matrix_pol,matrix_freq,matrix_pos,length_al,length_seg,nsamuser,npops,svratio,missratio,include_unknown,
						sum_sam,tcga,matrix_sv,nmhits,output,ploidy,outgroup_presence,nsites1_pop,nsites1_pop_outg,
						nsites2_pop,nsites2_pop_outg,nsites3_pop,nsites3_pop_outg,anx,bnx,anxo,bnxo,lengthamng,lenghtamng_outg,mhitbp,matrix_pol_tcga,0) == 0) {
			for(x=0;x<n_sam;x++) free(names[x]); free(names);
			free(DNA_matr);
			free(DNA_matr2);
			free(matrix_sizepos);
			free(matrix_segrpos);
			free(mhitbp);
            return(0);
        }
		free(mhitbp);
        for(x=0;x<n_sam;x++) free(names[x]); free(names);
	}
    return(1);
}

int var_char(FILE *file_input,SGZip *input_gz,FILE *file_logerr,SGZip *file_logerr_gz,long int *count,int *c,int *n_sam,long int *n_sit,int *nseq,int *maxsam,char ***names,char **DNA_matr,
	long int *n_site,int excludelines,char *name_excluded,int *n_excl,int includelines,char *name_ingroups,char *name_outgroup,int outgroup,int nsamuser_eff,char *ploidy)
{
    int  aa = 0;
    int  bb = 0;
    long int  dd;
    double ee;
    char *strin;
    /*long long t;*/
    
    aa = assigna(file_input,input_gz,file_logerr,file_logerr_gz,c,nseq,maxsam,names);
	if(aa == 1) {
		if(outgroup > 0) {
			if(((strin = strstr(names[0][*nseq-1],name_outgroup)) == 0)) {
				if(excludelines > 0) {
					if(((strin = strstr(names[0][*nseq-1],name_excluded)) != 0)) {
						*nseq -= 1;
						*n_excl += 1;
						*c = fzgetc(file_input, input_gz);
						while(*c != '*' && *c != -1 && *c != 0 && *c != '>')
							*c = fzgetc(file_input, input_gz);
						aa = 0;
					}
				}
				if(includelines > 0) {
					if(((strin = strstr(names[0][*nseq-1],name_ingroups)) == 0)) {
						*nseq -= 1;
						*n_excl += 1;
						*c = fzgetc(file_input, input_gz);
						while(*c != '*' && *c != -1 && *c != 0 && *c != '>')
							*c = fzgetc(file_input, input_gz);
						aa = 0;
					}
				}
			}
		}
		else {
			if(excludelines > 0) {
				if(((strin = strstr(names[0][*nseq-1],name_excluded)) != 0)) {
					*nseq -= 1;
					*n_excl += 1;
					*c = fzgetc(file_input, input_gz);
					while(*c != '*' && *c != -1 && *c != 0 && *c != '>')
						*c = fzgetc(file_input, input_gz);
					aa = 0;
				}
			}
			if(includelines > 0) {
				if(((strin = strstr(names[0][*nseq-1],name_ingroups)) == 0)) {
					*nseq -= 1;
					*n_excl += 1;
					*c = fzgetc(file_input, input_gz);
					while(*c != '*' && *c != -1 && *c != 0 && *c != '>')
						*c = fzgetc(file_input, input_gz);
					aa = 0;
				}
			}
		}
	}
    if(aa == 1) {
        while(bb == 0) {
            dd = (long int)floor((double)*count/(double)10000);
            ee = (double)*count/(double)10000;
            
            if(dd == ee) {
				/*
                if((t=(((long int)dd+(long int)1)*(long int)10000)) > 2147483647) {
                    fzprintf(file_logerr,file_logerr_gz,"Error: too much positions.\n");
                    return(0);
                }
				*/
                if((*DNA_matr = realloc(*DNA_matr,((long long)dd+(long long)1)*(long long)10000*sizeof(char))) == 0) {
                    fzprintf(file_logerr,file_logerr_gz,"Error: realloc error varchar.1\n");
                    return(0);
                }    
            }
            switch(*c/* = fzgetc(file_input, input_gz)*/) {
                case 'T':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '1';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 't':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '1';
					if(ploidy[0]=='2')
						DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 't';
					*count += 1;
                    *n_sit += 1;
                    break;
                case 'U':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '1';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'u':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '1';
					if(ploidy[0]=='2')
						DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 't';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'C':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '2';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'c':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '2';
					if(ploidy[0]=='2')
						DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 'c';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'G':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '3';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'g':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '3';
					if(ploidy[0]=='2')
						DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 'g';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'A':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '4';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'a':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '4';
					if(ploidy[0]=='2')
						DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 'a';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 0:
                    bb = 1; /*in FASTA case*/
                    break;
                case -1:
                    bb = 1; /*in FASTA case*/
                    break;
                case 10:
                    break;
                case 13:
                    break;
                case 32:
                    break;
                case '*': /* in NBRF case*/
                    bb = 1;
                    break;
                case '>': /* in FASTA case*/
                    bb = 1;
                    break;
				case 'N':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '5';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'n':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '5';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case '?':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '5';
                    *count += 1;
                    *n_sit += 1;
                    break;
				case '-':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '6';
                    *count += 1;
                    *n_sit += 1;
                    break;
                /*degenerated are converted to N*/
                case 'W':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 'w';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'w':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 'w';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'M':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 'm';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'm':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 'm';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'R':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 'r';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'r':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 'r';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'Y':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 'y';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'y':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 'y';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'K':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 'k';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'k':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 'k';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'S':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 's';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 's':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = 's';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'b':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '5';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'B':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '5';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'd':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '5';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'D':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '5';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'h':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '5';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'H':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '5';
                    *count += 1;
                    *n_sit += 1;
                    break;
				case 'v':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '5';
                    *count += 1;
                    *n_sit += 1;
                    break;
                case 'V':
                    DNA_matr[0][(((long long)*n_site*(unsigned long)*n_sam)+(unsigned long)*n_sit)] = '5';
                    *count += 1;
                    *n_sit += 1;
                    break;
				default:
                    fzprintf(file_logerr,file_logerr_gz,"Unexpected value in file");
                    printf("%d",*c);
                    return(0);
                    break;
            }
            if(*c !='>'/*|| *c != 0 || *c != -1*/) *c = fzgetc(file_input, input_gz);
        }
        *n_sam += 1;
        if(*n_site == 0) *n_site = *n_sit;
        else if(*n_site != *n_sit) {
            fzprintf(file_logerr,file_logerr_gz,"The number of sites are not equal in all lines in the alignment.");
            return(0);
        }
    }
    return 1;
}

int assigna(FILE *file_input,SGZip *input_gz,FILE *file_logerr,SGZip *file_logerr_gz,int *c,int *nseq,int *maxsam,char ***names)
{
    int N_VAR = 2;
    char var_file[2][50]  =
    {
        { ">"},
        { ">DL;"},
    };

    int i_;
    int j;
    int nn;
    int x;
    int c0; 
    
    j = 0;
    for(i_=0;i_<N_VAR;i_++) {
        while((var_file[i_][j]) == *c && (var_file[i_][j]) != '\0' && c != '\0') {
            *c = fzgetc(file_input, input_gz);
            j++;
        }
    }
    if(j<4 && j>0) i_= 1;/*FASTA*/
    else if(j==4) i_= 2;/*NBRF*/
        else {
            i_=0;/*NO ACCEPTED*/
            if(*c != 0 && *c != -1) *c = fzgetc(file_input, input_gz);
        }
        
    if((i_ == 2  && *c != 0 && *c != -1)) { /* NBRF files */
        nn = 0;
        while(*c != '\0' && *c != 13 && *c != 10 && nn < 50-2) {
            if(*c != '\t' && *c != 32) {
				names[0][*nseq][nn] = (char)*c;
				nn++;
			}
			*c = fzgetc(file_input, input_gz);
        }
        names[0][*nseq][nn] = '\0';
        *nseq += 1;
        if(*nseq == *maxsam) {
            *maxsam += 32;
            if(*maxsam > 32767) {
                fzprintf(file_logerr,file_logerr_gz,"\n Sorry, no more samples are allowed.");
                    return 0;
            }
            if ((*names = (char **)realloc(*names,*maxsam*sizeof(char *))) == 0) {
                fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. assigna.1 \n");
                return(0);
            }
            for(x=*nseq;x<*maxsam;x++) {
                if ((names[0][x] = (char *)calloc(50,sizeof(char))) == 0) {
                    fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. assigna.2 \n");
                    return(0);
                }
            }
        }
        /*use unix or macos or dos format. begin*/
        while(*c != '\0' && *c != 13 && *c != 10 && *c != -1 && *c != 0)
            *c = fzgetc(file_input, input_gz);

        c0 = *c;
        *c = fzgetc(file_input, input_gz);
        if(c0 == 13 && *c == 10) *c = fzgetc(file_input, input_gz);
        while(*c != 10 && *c != 13 && *c != -1 && c != 0) 
            *c = fzgetc(file_input, input_gz);
        if(*c == -1 || *c == 0) {
            fzprintf(file_logerr,file_logerr_gz,"\n Unexpected end of file");
            return(0);
        }
        c0 = *c;
        *c = fzgetc(file_input, input_gz);
        if(c0 == 13 && *c == 10) *c = fzgetc(file_input, input_gz);
        if(*c == -1 || *c == 0) {
            fzprintf(file_logerr,file_logerr_gz,"\n Unexpected end of file");
            return(0);
        }
        /*use unix or macos or dos format. end*/
        
        return(1);
    }
    else {	
        if(i_ == 1 && *c != 0 && *c != -1) {	/* FASTA files */
            nn = 0;
            while(*c != '\0' && *c != 13 && *c != 10 && nn < 50-2) {
				if(*c != '\t' && *c != 32) {
					names[0][*nseq][nn] = (char)*c;
					nn++;
				}
				*c = fzgetc(file_input, input_gz);
            }
            names[0][*nseq][nn] = '\0';
            *nseq += 1;
            if(*nseq == *maxsam) {
                *maxsam += 32;
                if(*maxsam > 32767) {
                    fzprintf(file_logerr,file_logerr_gz,"\n Sorry, no more samples are allowed.");
                        return 0;
                }
                if ((*names = (char **)realloc(*names,*maxsam*sizeof(char *))) == 0) {
                    fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. assigna.1 \n");
                    return(0);
                }
                for(x=*nseq;x<*maxsam;x++) {
                    if ((names[0][x] = (char *)calloc(50,sizeof(char))) == 0) {
                        fzprintf(file_logerr,file_logerr_gz,"\nError: memory not reallocated. assigna.2 \n");
                        return(0);
                    }
                }
            }
            /*use unix or macos or dos format. begin*/
            while(*c != '\0' && *c != 13 && *c != 10 && *c != -1 && *c != 0)
                *c = fzgetc(file_input, input_gz);
    
            c0 = *c;
            *c = fzgetc(file_input, input_gz);
            if(c0 == 13 && *c == 10) *c = fzgetc(file_input, input_gz);
            if(*c == -1 || *c == 0) {
                fzprintf(file_logerr,file_logerr_gz,"\n Unexpected end of file");
                return 0;
            }
            /*use unix or macos or dos format. end*/
            return(1);
        }
        else return(0);
    }
    return 0;
}


