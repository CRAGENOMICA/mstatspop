/*
 *  get_msdata.c
 *  mstatspop
 *
 *  Created by Sebastian Ramos-Onsins on 2010.
 *  Copyright 2010 CRAG. All rights reserved.
 *
 */

#include "get_msdata.h"
#include "log.h"

#define TOCONSIDER 0

int get_msdata( FILE *file_input,SGZip *input_gz,
        //FILE *file_logerr,SGZip *file_logerr_gz,
		char **matrix_pol, long int **matrix_freq, long int **matrix_pos,
		long int *length_seg, int *nsamuser, int npops, int nsamtot, long int length, long int *nmhits, 
		int *matrix_mask, float *vector_mask, float svratio, double **vector_priors, int *npriors,
		long int **matrix_sv, int outgroup_presence, int force_outgroup, float freq_revert, double *sum_sam,
		double *nsites1_pop, double *nsites1_pop_outg, int formatfile,
		double *nsites2_pop,double *nsites2_pop_outg,double *nsites3_pop,double *nsites3_pop_outg,
		double *anx, double *bnx,double *anxo, double *bnxo,double **lengthamng,double **lengthamng_outg, int include_unknown,
		char *file_mas, double freq_missing_ms, int kind_length,
		double **sum_sam_mask, double *length_mask, long int *length_mask_real, double *missratio, int location_missing_ms,int *sort_nsam)
{
	static char c[1],cc[1];
	long int li,li2,max_biasites = 256;
	char value[20];
	int i,y,y2,k,n2;

#if TOCONSIDER==1
    long int posmhitinit;
	int t_0=0;
	int t_1=0;
	int s1_0=0;
	int s1_1=0;
	int s2_0=0;
	int s2_1=0;
	double r0,r1;
    int physical;
    long int *matrix_phys;
#endif
	int ft;
	long int n,xx,li3;
	int x,z,v,p,z2,v2;
	
	long int *mhitbp;
	long int *discarded;
	double freq_missing_ms2=0.;
    int *matrix_mask2=0;
	
	/*matrix_sv: number of transitions and transversions. We have to give these results to the main function*/
	
	*nmhits = 0;

	/*READ THE MS FILE */
	
	*cc = fzgetc(file_input, input_gz);
	if(*cc==0 || *cc==-1) return 1;
	*c  = fzgetc(file_input, input_gz);
	if(*c==0 || *c==-1) return 1;
	while ((*c != '/' || *cc != '/') && *c != 0 && *c != -1) { /* read initial iteration signal: '//' */
		*cc = *c;
		*c  = fzgetc(file_input, input_gz);
	}
	if(*c==0 || *c==-1) return 1;
	if(*c == '/') {
		*c  = fzgetc(file_input, input_gz);
		if(*c==0 || *c==-1) return 1;
		
		ft=0;
		/*Read priors or other parameters*/
		if(*vector_priors==0) {
			if((*vector_priors = (double *)calloc((long int)1,sizeof(double)))==0) {
				//fprintf(file_logerr,"Error: memory not allocated. get_msdata.01"); 
                log_fatal("Error: memory not allocated, vector_priors get_msdata.01");
				return(1);
			}
			ft = 1;
		}
		*npriors = 0;
		while (*c != 10  && *c != 13 && *c != 0) {
			while(*c == 32 || *c == 9) {*c  = fzgetc(file_input, input_gz);}
			i = 0;
			while(*c!=32 && *c != 9 && *c != 10 && *c != 13) {
				value[i] = *c;
				*c  = fzgetc(file_input, input_gz);
				i++;
			}
			if(*c == 0 || *c == 10 || *c == 13) if(i==0) break;
			value[i] = '\0';
			vector_priors[0][*npriors] = (double)atof(value);
			*npriors += 1;
			if(ft) {
				if((*vector_priors = (double *)realloc(*vector_priors,(long int)(*npriors+1)*sizeof(double)))==0) {
					//fprintf(file_logerr,"Error: memory not allocated. get_msdata.02"); 
                    log_fatal("Error: memory not allocated, vector_priors get_msdata.02");
					return(1);
				}
			}
		}	
		ft = 0;
		/*... go next line*//*
		while (*c != 10) 
			*c  = fzgetc(file_input, input_gz);
		*/
        
        if((matrix_mask2 = (int *)calloc((unsigned int)(nsamtot+1)*length,sizeof(int))) == 0) {
            printf("Error allocating memory");
            exit(1);
        }

        /*BEGIN define mask: define first the ratio of missing values:*/
        /*formatfile=2 IS DEPRECATED*/
		/*if(formatfile == 2 && *npriors >= location_missing_ms && location_missing_ms > 0)
			freq_missing_ms2 = vector_priors[0][location_missing_ms-1];
		else */freq_missing_ms2 = freq_missing_ms;
		
		/*eliminate positions from mask*/ /*it is much faster to include a binomial distribution per position in the calculations!!*/
		if((file_mas[0] == '-'  && file_mas[1] == '1')/* || formatfile == 2*/) { /*all positions are accepted*/
			for(li2=0;li2<length;li2++) vector_mask[li2] = 1.0;
			if(freq_missing_ms2) {
				/* MODIFY TO INCLUDE CONDITIONAL LENGTH (kind_length:0,1,2,3) *//*
				for(li2=0;li2<(nsamtot-nsamuser[npops-1])*length;li2++) {
					if(ran1() <= freq_missing_ms2 && formatfile == 1) 
						matrix_mask[li2] = -1;
					else matrix_mask[li2] = 0;
				}
				*/
				for(li2=0;li2<length;li2++) {
					li3 = -1; /*number of samples with no-Ns from a binomial*/
					/*Even in the case of NGS reads, the distribution of Ns seems a binomial. */
					/*In cases with high depth the variance increase. Another factors also influence. Some regions have clearly lower depth than expected*/
					/*On the other hand, a binomial approach seems justified*/
					while(li3 < kind_length) /*CONDITIONAL on kind_length*/
						li3 = largebinomialdist(1.0-freq_missing_ms2, nsamtot-nsamuser[npops-1]);
					for(n=nsamtot-nsamuser[npops-1]-1;n>=0;n--) { /*define for each sample if N or nt*/
						if(ran1() > (double)li3/(double)(n+1)) {
							matrix_mask[(unsigned int)n*(unsigned int)length+(unsigned int)li2] = -1; /*N*/
						}
						else {
							matrix_mask[(unsigned int)n*(unsigned int)length+(unsigned int)li2] = 0; /*nt*/
							li3 -= 1;
						}
					}
				}
				/*END MODIFICATION*/
				/*the outgroup has no missing values!!*//* <<<<<<<<<<<<-----------------------------------!!!!!!!!*/
				for(li2=(nsamtot-nsamuser[npops-1])*length;li2<(nsamtot)*length;li2++) 
					matrix_mask[li2] = 0;
			}
			else {
				for(li2=0;li2<(nsamtot)*length;li2++) matrix_mask[li2] = 0;
			}

            /*if(include_unknown) {*/
                if((*sum_sam_mask = (double *)calloc(nsamtot,sizeof(double))) == 0) {
                    printf("Error allocating memory");
                    exit(1);
                }
            /*}*/
            *missratio = 0.;
            *length_mask = 0.;
            *length_mask_real = 0;
            li3 = 0;
            for(li2=0;li2<length;li2++) {
                if(vector_mask[li2] > 0.) {
                    if(npops > 1) {
                        x = 0;
                        for(n=0;n<nsamtot-nsamuser[npops-1];n++) {
                            li3 += 1;
                            if(matrix_mask[(unsigned int)n*(unsigned int)length+(unsigned int)li2] == 0) {
                                sum_sam_mask[0][n] += 1.; /*1 sum, 0 no sum*/
                                x += 1;
                            }
                            else *missratio += 1.;
                        }
                    }
                    else x = 1;
                    y = 0;
                    for(n=nsamtot-nsamuser[npops-1];n<nsamtot;n++) {
                        if(npops == 1) li3 += 1;
                        if(matrix_mask[(unsigned int)n*(unsigned int)length+(unsigned int)li2] == 0) {
                            sum_sam_mask[0][n] += 1.; /*1 sum, 0 no sum*/
                            y += 1;
                        }
                        else {
                            if(npops == 1) *missratio += 1.;
                        }
                    }
                    if(x > 0 && y > 0) {
                        *length_mask_real += 1; /*length where the outgroup and one of the rest lines exists*/
                        *length_mask += vector_mask[li2]; /*length where the outgroup and one of the rest lines exists*/
                    }
                    else {
                        if(x==0 || (y==0 && x>0)) {
                            missratio -= ((nsamtot - nsamuser[npops-1]) - x);
                            li3 -= (nsamtot - nsamuser[npops-1]);
                        }
                        if(y==0 && npops == 1) {
                            missratio -= (nsamtot);
                            li3 -= (nsamtot);
                        }
                    }
                }
            }
            if(li3) *missratio = (double)*missratio/(double)li3;
            else *missratio = 1.;
		}
		/*END defining mask in case m option ==-1 and include_unknown ==1*/
		
        while (*c==32 || *c == 9 || *c == 13 || *c == 10) *c = fzgetc(file_input, input_gz);
		/*READ the row with the number of segregating sites*/
		i = 0;
		while(*c!=32 && *c != 9 && *c != 10 && *c != 13) {
			value[i] = *c;
			*c  = fzgetc(file_input, input_gz);
			i++;
		}
		value[i] = '\0';
		if(strcmp(value,"segsites:") == 0) {
			while(*c == 32 || *c == 9) {*c  = fzgetc(file_input, input_gz);}
			i = 0;
			while(*c!=32 && *c != 9 && *c != 10 && *c != 13) {
				value[i] = *c;
				*c  = fzgetc(file_input, input_gz);
				i++;
			}
			value[i] = '\0';
		}
		max_biasites = (long int)atol(value);
		
		/*allocate memory*/
        if((*matrix_pol = (char *) calloc (((long int)nsamtot*(max_biasites+1)),sizeof(char))) == 0) {
            //fprintf(file_logerr,"Error: memory not allocated. get_msdata.3");
            log_fatal("Error: memory not allocated, matrix_pol get_msdata.3");
            return(1);
        }
		if((*matrix_freq = (long int *) calloc (max_biasites+1,sizeof(long int))) == 0) {
			//fprintf(file_logerr,"Error: memory not allocated. get_msdata.5"); 
            log_fatal("Error: memory not allocated, matrix_freq get_msdata.5");
			return(1);
		}
		if((*matrix_pos = (long int *) calloc (max_biasites+1,sizeof(long int))) == 0) {
			//fprintf(file_logerr,"Error: memory not allocated. get_msdata.4"); 
            log_fatal("Error: memory not allocated, matrix_pos get_msdata.4");
			return(1);
		}
		if((*matrix_sv = (long int *) calloc((long int)max_biasites+1,sizeof(long int)))==0) {
			//fprintf(file_logerr,"Error: memory not allocated. get_msdata.00"); 
            log_fatal("Error: memory not allocated, matrix_sv get_msdata.00");
			return(1);
		}
        if(max_biasites == 0)
		{
			*length_seg = 0;
			return 0;
		}
		
		if((mhitbp = (long int *) calloc (length+1, sizeof(long int))) == 0) {
			//fprintf(file_logerr,"Error: memory not reallocated. get_msdata.13"); 
            log_fatal("Error: memory not allocated, mhitbp get_msdata.13");
			return(0);
		}
#if TOCONSIDER == 0
        if((discarded = (long int *) calloc (max_biasites+1, sizeof(long int))) == 0) {
            //fprintf(file_logerr,"Error: memory not reallocated. get_msdata.13");
            log_fatal("Error: memory not allocated, discarded get_msdata.13");
            return(0);
        }
#else
        if((matrix_phys = (long int *) calloc (max_biasites+1,sizeof(long int))) == 0) {
            //fprintf(file_logerr,"Error: memory not allocated. get_msdata.4");
            log_fatal("Error: memory not allocated, matrix_phys get_msdata.4");
            return(1);
        }
		if((discarded = (long int *) calloc (max_biasites+1, sizeof(long int))) == 0) {
			//fprintf(file_logerr,"Error: memory not reallocated. get_msdata.13"); 
            log_fatal("Error: memory not allocated, discarded get_msdata.13");
			return(0);
		}
#endif
		/*READ the row header */
		while (*c != 10) {
			*c  = fzgetc(file_input, input_gz); /*next line*/
		}
		while (*c==32 || *c == 9 || *c == 13 || *c == 10) {*c  = fzgetc(file_input, input_gz);}
		i = 0;
		while(*c!=32 && *c != 9 && *c != 10 && *c != 13) {
			value[i] = *c;
			{*c  = fzgetc(file_input, input_gz);}
			i++;
		}
		value[i] = '\0';

        if(max_biasites) {

#if TOCONSIDER == 1
            physical = 0;
            if(strcmp(value,"physical:") == 0) {
                physical = 1;
                for(li=0;li<max_biasites;li++) {
                    while(*c == 32 || *c == 9 || *c == 13 || *c == 10) {*c  = fzgetc(file_input, input_gz);}
                    i = 0;
                    while(*c!=32 && *c != 9 && *c != 10 && *c != 13) {
                        value[i] = *c;
                        {*c  = fzgetc(file_input, input_gz);}
                        i++;
                    }
                    value[i] = '\0';
                    matrix_phys[li] = (long int)round((double)atof(value));
                }

                /*READ the row with the positions*/
                while (*c != 10) {
                    *c  = fzgetc(file_input, input_gz); /*next line*/
                }
                while (*c==32 || *c == 9 || *c == 13 || *c == 10) {*c  = fzgetc(file_input, input_gz);}
                i = 0;
                while(*c!=32 && *c != 9 && *c != 10 && *c != 13) {
                    value[i] = *c;
                    {*c  = fzgetc(file_input, input_gz);}
                    i++;
                }
                value[i] = '\0';
            }
#endif
            /*Take the "positions" in physical integer positions (starting from 1!!!): */
            /*if mhit, keep the physical position (mhitbp) in order to consider*/
            /*will keep only those positions that reverted to 1/2 variants, otherwise reject*/
            if(strcmp(value,"positions:") == 0) {
#if TOCONSIDER==1
                posmhitinit = -1;
#endif
                while(*c == 32 || *c == 9 || *c == 13 || *c == 10) {*c  = fzgetc(file_input, input_gz);}
                for(li=0;li<max_biasites;li++) {
                    while(*c == 32 || *c == 9) {*c  = fzgetc(file_input, input_gz);}
                    if(*c == 13 || *c == 10) {max_biasites = li+1; break;}
                    i = 0;
                    while(*c!=32 && *c != 9 && *c != 10 && *c != 13) {
                        value[i] = *c;
                        {*c  = fzgetc(file_input, input_gz);}
                        i++;
                    }
                    value[i] = '\0';
                    matrix_pos[0][li] = round(((double)atof(value)*length)+1);/*start from 1!!!*/
                    /*Mandatory to count a negative (missing) value at position 0!!!*/
                    if(matrix_pos[0][li] == (long int)length) matrix_pos[0][li] = (long int)length+1-1;
                    /*in case ms has a position value of 1.000 by rounding (...)*/
    #if TOCONSIDER==0
                    /*TOO COMPLEX: REJECT MHIT POSITIONS DIRECTLY*/
                    /*use this code instead*/
                    if(li) {
                        if(labs(matrix_pos[0][li]) == labs(matrix_pos[0][li])-1) {
                            if(mhitbp[(long int)labs(matrix_pos[0][li])] == 0) *nmhits += 1;
                            mhitbp[(long int)labs(matrix_pos[0][li])] = 1;
                        }
                    }
    #else
                    if(svratio > -1.0) {
                        if(li==0) {/*first position: transition or transversion?*/
                            if((r0=ran1()) < (svratio/(svratio+1.))) {matrix_sv[0][li]=1;t_0 = 1;s1_0 = 0;s2_0 = 0;}
                            else {
                                if(r0 > (svratio/(svratio+1.)) + (1.-svratio/(svratio+1.))/2.) {
                                    matrix_sv[0][li]=2;t_0 = 0;s1_0 = 1;s2_0 = 0;
                                }
                                else {matrix_sv[0][li]=3;t_0 = 0;s1_0 = 0;s2_0 = 1;}
                            }
                        }
                        if(li) {
                            if(labs(matrix_pos[0][li]) == labs(matrix_pos[0][li])-1) {
                                /*mhit: can we recover it? (the same two nucleotides)*/
                                if(posmhitinit == -1) {
                                    posmhitinit = li-1; /*locate first column of mhit position*/
                                    /*define if transition,transversion1 or transversion2 (t,s1,s2) for column0*/
                                }
                                /*if(matrix_pos[0][posmhitinit] >= 0) {*//*now always true*/
                                    /*define if transition,transversion1 or transversion2 (t,s1,s2) for given column*/
                                    if((r1=ran1()) < (svratio/(svratio+1.))) {
                                        t_1 = 1;s1_1 = 0;s2_1 = 0;
                                    }
                                    else {
                                        if(r1 < (svratio/(svratio+1.)) + (1.-svratio/(svratio+1.))/2.) {
                                            t_1 = 0;s1_1 = 1;s2_1 = 0;
                                        }
                                        else {
                                            t_1 = 0;s1_1 = 0;s2_1 = 1;
                                        }
                                    }
                                /*}*/
                                if(!(t_0 == t_1 && s1_0 == s1_1 && s2_0 == s2_1)) { /*mhit*/
                                    /*when it is a mhit */
                                    /* count the position mhit in a matrix*/
                                    mhitbp[(long int)labs(matrix_pos[0][li])] = 1;
                                    if(t_1)  matrix_sv[0][li]=1;
                                    else if(s1_1) matrix_sv[0][li]=2;
                                    else if(s2_1) matrix_sv[0][li]=3;
                                    for(li2=posmhitinit;li2<=li;li2++)
                                        discarded[li2] = 1;/*mhit: discard the positions*/
                                }
                                else {
                                    if(discarded[posmhitinit]==0) {/*if there were more mhits we discard the position*/
                                        /*mhit reverted to the same two nucleotides: no observed mhit*/
                                        /*in fact ms is not using discrete positions. It is only an approach when there is recombination*/
                                        mhitbp[(long int)labs(matrix_pos[0][li])] = 2;/*mhit but not seen*/
                                        matrix_sv[0][li] = matrix_sv[0][posmhitinit];
                                    }
                                    else
                                        discarded[li] = 1;
                                }
                            }
                            else {
                                /*if no mhit calculate if transition or transversion*/
                                if((r0=ran1()) < (svratio/(svratio+1.))) {matrix_sv[0][li]=1;t_0 = 1;s1_0 = 0;s2_0 = 0;}
                                else {
                                    if(r0 < (svratio/(svratio+1.)) + (1.-svratio/(svratio+1.))/2.) {
                                        matrix_sv[0][li]=2;t_0 = 0;s1_0 = 1;s2_0 = 0;
                                    }
                                    else {matrix_sv[0][li]=3;t_0 = 0;s1_0 = 0;s2_0 = 1;}
                                }
                                posmhitinit = -1;
                            }
                        }
                    }
    #endif
                }
            }
                    
#if TOCONSIDER==1
            /*READ the matrix*/
            /*we accept values 0,1 for variants (DEPRECATED!: TO ERASE: and 8 and 9 for gaps and missing values)*/
            for(y=0;y<nsamtot-!outgroup_presence;y++) {
                while(*c == 32 || *c == 9 || *c == 13 || *c == 10) {*c  = fzgetc(file_input, input_gz);}
                li=0; /*li is the counter for max_bialsites.*/
                li2=0;/*li2 only count if the physical position is not equal*/
                k=0;/*index to control hidden mhits (2 mutations in the same position with same variants)*/
                for(li=0;li<max_biasites;li++)/*while(*c != 10 && *c != 13)*/ {
                    if(li==0 || (li > 0 && labs(matrix_pos[0][li]) != labs(matrix_pos[0][li])-1)) k=0; /*{*/
                    if(mhitbp[(long int)labs(matrix_pos[0][li])] == 0) {
                        if(li>0) li2++;
                        if(matrix_mask[y*length+(long int)labs(matrix_pos[0][li])-1] == -1 ||
                           matrix_mask2[y*length+(long int)labs(matrix_pos[0][li])-1] == -1 || *c=='8' || *c == '9') {/*missing*/
                            matrix_pol[0][((li2*nsamtot)+y)] = '-';
                            matrix_mask2[y*length+(long int)labs(matrix_pos[0][li])-1] = -1;
                            matrix_pos[0][li] = -(long int)labs(matrix_pos[0][li]);
                            /* <---- means this position has mising values (Ns)*/
                        }	
                        if(matrix_mask[y*length+(long int)labs(matrix_pos[0][li])-1] == 0 &&
                           matrix_mask2[y*length+(long int)labs(matrix_pos[0][li])-1] == 0 && *c!='8' && *c != '9') {/*no missing*/
                            matrix_pol[0][((li2*nsamtot)+y)] = *c;
                        }
                        if(matrix_mask[y*length+(long int)labs(matrix_pos[0][li])-1] != -1 && *c > '0' && *c != '8' && *c != '9') {
                            matrix_freq[0][li2] += 1;
                        }
                        if(outgroup_presence == 0 && y == nsamtot-2) {
                            /*make a copy of the last line and calculate the total frequency of derived*/
                            if(ran1()<freq_revert) {
                                y2 = y;
                                while(y2 > 0 && matrix_pol[0][((li2*nsamtot)+y2)] == '-') y2--;
                                matrix_pol[0][((li2*nsamtot)+y+1)] = matrix_pol[0][((li2*nsamtot)+y2)];
                                if(matrix_pol[0][((li2*nsamtot)+y2)] != '0' && matrix_pol[0][((li2*nsamtot)+y2)] != '-')
                                    matrix_freq[0][li2] += 1;
                            }
                            else 
                                matrix_pol[0][((li2*nsamtot)+y+1)] = '0';
                        }
                    }
                    /*}*/
                    else {/*mhit*/
                        if(discarded[li] == 0) {
                            /* mhit with two variants (invisible mhit)*/
                            if(matrix_mask[y*length+(long int)labs(matrix_pos[0][li])-1] == -1 ||
                               matrix_mask2[y*length+(long int)labs(matrix_pos[0][li])-1] == -1 || *c=='8' || *c == '9') {
                                /*missing data*/
                                matrix_pol[0][(((li2)*nsamtot)+y)] = '-';
                                matrix_mask2[y*length+(long int)labs(matrix_pos[0][li])-1] = -1;
                                matrix_pos[0][li2] = -(long int)labs(matrix_pos[0][li2]);
                                /* <---- means this position has mising values ('N')*/
                                k+=1;
                            }
                            else {
                                if(matrix_mask[y*length+(int)labs(matrix_pos[0][li])-1] == 0 &&
                                   matrix_mask2[y*length+(int)labs(matrix_pos[0][li])-1] == 0) {
                                    /*no missing position in mask neither missing in *c */
                                    if(matrix_pol[0][(long int)(((li2)*(long int)nsamtot)+y)] == '0' && *c != '0') {
                                        /*mutation in *c, change matrix_pol and count freq*/
                                        matrix_pol[0][(long int)(((li2)*(long int)nsamtot)+y)] = *c;
                                        matrix_freq[0][(long int)li2] += (long int)1;
                                    }
                                    /*if(matrix_pol[0][(((li2)*nsamtot)+y)] != '0' && *c == '0') {*/
                                        /*no mutation in *c, nothing changes*/
                                    /*}*/
                                    if(matrix_pol[0][(long int)(((li2)*(long int)nsamtot)+y)] != '0' && *c != '0') {
                                        /*mutation in both positions: revert to '0'*/
                                        matrix_pol[0][(long int)(((li2)*(long int)nsamtot)+y)] = '0';/*revert the mutation to ancestral*/
                                        matrix_freq[0][(long int)li2] -= (long int)1;
                                    }
                                }
                                k+=1;
                            }
                            if(y==nsamtot-!outgroup_presence-1 && k>1)
                                discarded[li]=1;
                        }
                    }
                    *c  = fzgetc(file_input, input_gz);
                }
            }
#else
            /*READ the matrix*/
            /*we accept values 0,1 for variants (DEPRECATED!: TO ERASE: and 8 and 9 for gaps and missing values)*/
            for(y=0;y<nsamtot-!outgroup_presence;y++) {
                while(*c == 32 || *c == 9 || *c == 13 || *c == 10) {*c  = fzgetc(file_input, input_gz);}
                li=0; /*li is the counter for max_bialsites.*/
                for(li=0;li<max_biasites;li++) {
                    if(matrix_mask[y*length+(long int)labs(matrix_pos[0][li])-1] == -1) {/*missing*/
                        matrix_pol[0][((li*nsamtot)+y)] = '-';
                        matrix_pos[0][li] = -(long int)labs(matrix_pos[0][li]);
                    }
                    if(matrix_mask[y*length+(long int)labs(matrix_pos[0][li])-1] == 0) {/*no missing*/
                        matrix_pol[0][((li*nsamtot)+y)] = *c;
                    }
                    if(matrix_mask[y*length+(long int)labs(matrix_pos[0][li])-1] != -1 && *c > '0') {
                        matrix_freq[0][li] += 1;
                    }
                    if(outgroup_presence == 0 && y == nsamtot-2) {
                        /*make a copy of the last line and calculate the total frequency of derived*/
                        if(ran1()<freq_revert) {
                            y2 = y;
                            while(y2 > 0 && matrix_pol[0][((li*nsamtot)+y2)] == '-') y2--;
                            matrix_pol[0][((li*nsamtot)+y+1)] = matrix_pol[0][((li*nsamtot)+y2)];
                            if(matrix_pol[0][((li*nsamtot)+y2)] != '0' && matrix_pol[0][((li*nsamtot)+y2)] != '-')
                                matrix_freq[0][li] += 1;
                        }
                        else
                            matrix_pol[0][((li*nsamtot)+y+1)] = '0';
                    }
                    *c  = fzgetc(file_input, input_gz);
                }
            }
#endif
        }
        /*eliminate those positions that are partial (syn/nsyn) according to their probability (ex. if it is a Syn with 0.6 it has 0.4 to be NSyn and in that case will not be counted): FOR SIMULATIONS!!! */
        /* (not really correct, it is correlated to svratio, but approx ok): TOO COMPLEX, AVOID IT*/
       /**/
        *length_seg = 0;
        for(li=0;li<max_biasites;li++) {
            if((include_unknown == 0 && matrix_pos[0][li] >= 0) || (include_unknown == 1)) {
                if(li==0 || (li > 0 && (long int)labs(matrix_pos[0][li]) != (long int)labs(matrix_pos[0][li])-1)) {
                    if(vector_mask[(long int)labs(matrix_pos[0][li])-1] == 1.0)
                        *length_seg += 1;
                    else {
                        if(ran1() < vector_mask[(long int)labs(matrix_pos[0][li])-1]) {
                            *length_seg += 1;
                        }
                        else {
                            discarded[li] = 1;
                            /*vector_mask[(long int)labs(matrix_pos[0][li])-1] = 0.;*/
                            /*matrix_pos[0][li] = -999999999*/ /*-(long int)fabs(matrix_pos[0][li])*/;
                        }
                    }
                }
                /**/
                if(outgroup_presence == 1) {
                    n2=0;
                    for(n=nsamtot-nsamuser[npops-1];n<nsamtot;n++) {
                        if(matrix_mask[(unsigned int)n*(unsigned int)length+(unsigned int)labs(matrix_pos[0][li])-1] == -1) {
                            n2 += 1;
                        }
                        else
                            break;
                    }
                    if(n2 == nsamuser[npops-1]) {
                        discarded[li] = 1;
                     }
                }
                /**/
            }
#if TOCONSIDER==0
            if(mhitbp[(long int)labs(matrix_pos[0][li])] == 1) { /*keep the highest frequency*/
                li2 = li + 1;
                while(mhitbp[(long int)labs(matrix_pos[0][li])] == mhitbp[(long int)labs(matrix_pos[0][li2])]) {
                    if(matrix_freq[0][li2] > matrix_freq[0][li]) {
                        discarded[li] = 1;
                        li = li2;
                    }
                    else {
                        discarded[li2] = 1;
                    }
                    li2 += 1;
                }
                li = li2;
            }
#endif
        }
        /**/
        
        /*change '-1' values to '0' in matrix_mask if it is a VARIANT position AND includes '8' or '9' [format ms_e]*/
        /*discarded*/
        /*
        if(formatfile == 2){
            for(li=0;li<max_biasites;li++) {
                if(matrix_pos[0][li] > 0) {
                    for(y=0;y<nsamtot-!outgroup_presence;y++) 
                        matrix_mask2[y*length+matrix_pos[0][li]-1] = 0;
                }
            }
        }
        */
#if TOCONSIDER==1
        /*re-calculate the matrix_sv: 1:transition, 2:transversion*/
		for(li=0;li<max_biasites;li++) {
			if(matrix_sv[0][li] == 3) matrix_sv[0][li] = 2;
		}
        
        *nmhits=0;
        posmhitinit = -1;
        for(li=0;li<length;li++) {
            /*printf("\nmhitbp[%ld]=%ld",li+1,mhitbp[li+1]);*/
            if(mhitbp[li+1] == 1) {
                *nmhits += 1;
            }
        }
#endif
		/*recalculate matrix_pos: eliminate all discarded values. Also count the number of total positions*/
		li2 = 0;
		for(li=0;li<max_biasites;li++) {
            /*printf("\nmatrix_pos[%ld]=%ld\tdiscarded[%ld]=%ld",li,matrix_pos[0][li],li,discarded[li]);*/
            if(discarded[li] == 0) {
 				matrix_pos[0][li2] = matrix_pos[0][li];
                for(y=0;y<nsamtot;y++)
                    matrix_pol[0][li2*nsamtot+y] = matrix_pol[0][li*nsamtot+y];
				li2++;
			}
		}
		for(li=li2;li<max_biasites;li++) {
			matrix_pos[0][li] = -999999999;
            for(y=0;y<nsamtot;y++)
                matrix_pol[0][li*nsamtot+y] = '-';
		}
		
#if TOCONSIDER==1
        /*rename the positions of the matrix_pos by the physical positions, if the information is included:*/
		if(physical) {
			for(li=0;li<max_biasites;li++) {
				if(matrix_pos[0][li] >= 0) {
					matrix_pos[0][li] =  (long int)labs(matrix_phys[li]);
				}
				else {
					matrix_pos[0][li] = -(long int)labs(matrix_phys[li]);
				}
			}
		}
#endif
		/*count number of total positions (excluding gaps and mhits) per line*/
		for(li=0;li<li2;li++) {
			for(y=0;y<nsamtot;y++) {
				if(matrix_pol[0][(li*nsamtot)+y] != '-')
					sum_sam[y] += (double)1;
			}
		}
		
		/*calculate the number of positions per population considering or not the outgroup*/
		/*calculate nsites1_pop and nsites1_pop_outg*/
		for(y=0;y<npops-!outgroup_presence;y++) {
			nsites1_pop[y] = 0.;
			nsites1_pop_outg[y] = 0.;
			nsites2_pop[y] = 0.;
			nsites2_pop_outg[y] = 0.;
			nsites3_pop[y] = 0.;
			nsites3_pop_outg[y] = 0.;
		}
		z = 0;
		for(y=0;y<npops-1;y++) {
			anx[y]  = bnx[y]  = 0.0;
			anxo[y] = bnxo[y] = 0.0;
			for(xx=0;xx<length;xx++) {
				v = 0;
				for(x=z;x<z+nsamuser[y];x++) {
					if(matrix_mask[x*length+xx] == -1 || matrix_mask2[x*length+xx] == -1) v += 1; /*number of Ns*/
				}
/*
				if(discarded[xx] == 0 && ((include_unknown == 1 && vector_mask[xx] > 0) ||
				   (include_unknown==0 && v == 0 && vector_mask[xx] > 0))) {
*/
                if(mhitbp[xx+1] != 1 &&
                   ((include_unknown == 1 && vector_mask[xx] > 0) ||
                    (include_unknown == 0 && vector_mask[xx] > 0  && v == 0))) {

                    if(v < nsamuser[y])
						nsites1_pop[y] += vector_mask[xx];
					if(v < nsamuser[y]-1) 
						nsites2_pop[y] += vector_mask[xx];
					if(v < nsamuser[y]-2) 
						nsites3_pop[y] += vector_mask[xx];
					for(k=1;k<nsamuser[y]-v;k++) {
						anx[y] += 1.0/((double)k);
						bnx[y] += 1.0/((double)k*(double)k);
					}
					p = 0;
					if(outgroup_presence == 1 || force_outgroup == 1) {//modified 20220603
						for(x=nsamtot-1;x>=nsamtot-nsamuser[npops-1];x--)
							if(matrix_mask2[x*length+xx] == -1 || matrix_mask[x*length+xx] == -1 ) p += 1;
						if(p < nsamuser[npops-1] && v < nsamuser[y])
							nsites1_pop_outg[y] += vector_mask[xx];
						if(p < nsamuser[npops-1] && v < nsamuser[y]-1)
							nsites2_pop_outg[y] += vector_mask[xx];
						if(p < nsamuser[npops-1] && v < nsamuser[y]-2)
							nsites3_pop_outg[y] += vector_mask[xx];
						if(p < nsamuser[npops-1] && v < nsamuser[y]) {
							for(k=1;k<nsamuser[y]-v;k++) {
								anxo[y] += 1.0/((double)k);
								bnxo[y] += 1.0/((double)k*(double)k);
							}
						}
					}
				}
			}
			z += nsamuser[y];
			/*nsites1_pop[y] -= *nmhits;
			nsites1_pop_outg[y] -= *nmhits;
			nsites2_pop[y] -= *nmhits;
			nsites2_pop_outg[y] -= *nmhits;*/
			if(nsites2_pop[y]>0) {
				anx[y] /= (double)nsites2_pop[y];
				bnx[y] /= (double)nsites2_pop[y];
			}
			if(nsites2_pop_outg[y]>0) {
				anxo[y] /= (double)nsites2_pop_outg[y];
				bnxo[y] /= (double)nsites2_pop_outg[y];
			}
		}

		
		/*calculate lengthamng*/
        for(y=0;y<npops-1;y++) {
            for(y2=y+1;y2<npops-1;y2++) {
                lengthamng[y][y2] = 0.;
                lengthamng_outg[y][y2] = 0.;
            }
        }
		for(xx=0;xx<length;xx++) {
			z = 0;
			for(y=0;y<npops-1;y++) {
				v = 0;
				for(x=z;x<z+nsamuser[y];x++) {
					if(matrix_mask[x*length+xx] == -1 || matrix_mask2[x*length+xx] == -1) v += 1;
				}
				z2 = 0; x= 0;
				while(x<y+1) {z2 += nsamuser[x]; x++;}
				for(y2=y+1;y2<npops-1;y2++) {
					v2 = 0;
					for(x=z2;x<z2+nsamuser[y2];x++) {
						if(matrix_mask[x*length+xx] == -1 || matrix_mask2[x*length+xx] == -1) v2 += 1;
					}
					p = 0;
					if(mhitbp[xx]!=1 && ((include_unknown == 1 && vector_mask[xx] > 0) ||
					   (include_unknown==0 && v == 0 && vector_mask[xx] > 0))) {
						if(outgroup_presence == 1 || force_outgroup == 1) { //modified 20220603
							for(x=nsamtot-1;x>=nsamtot-nsamuser[npops-1];x--) {
								if(matrix_mask2[x*length+xx] == -1/**/ || matrix_mask[x*length+xx] == -1/**/) p += 1;
							}
							if(p < nsamuser[npops-1] && v < nsamuser[y] && v2 < nsamuser[y2]) 
								lengthamng_outg[y][y2] += vector_mask[xx];
						}
                        if(v < nsamuser[y] && v2 < nsamuser[y2])
                            lengthamng[y][y2] += vector_mask[xx];
					}
					z2 += nsamuser[y2];
				}
				z += nsamuser[y];
			}
		}
		/**/
#if TOCONSIDER==1
        free(matrix_phys);
#endif
        free(discarded);
		free(mhitbp);
        free(matrix_mask2);
		
		return 0;
	}
	else {
		//fprintf(file_logerr,"\nError: ms file finished abruptely.\n");
        log_error("Error: ms file finished abruptely.");
		exit(1);
	}
	
	return 1;
}
