/*
 *  jointfreqdist.c
 *  mstatspop
 *
 *  Created by Sebastian Ramos-Onsins on 25/02/2010.
 *  Copyright 2010 CRAG. All rights reserved.
 *
 */
#include "jointfreqdist.h"
#include "sancestral.h"


int jointfreqdist(	int npops, int *nsam, char *matrix_pol,long int *matrix_pos,
					long int length, struct stats *statistics,long int *sites_matrix,
					double **jfd, int **nfd,int outgroup_presence,int force_outgroup )
{
	long int j;
	int x,sumnsam;
	int a0,a1,a2,pol,polc,anc,ancl,n,popo,pop1,polv;
	int *initsq1,initso,inits;

    initsq1 = (int *)calloc(npops,sizeof(int));

    sumnsam = 0;
    for(x=0;x<npops;x++) {
        initsq1[x] = sumnsam;
        sumnsam += nsam[x];
    }
    
    a0=a1=0;
    for(j=0;j<length;j++) {
        if((ispolmis(matrix_pol,sumnsam,j,length)) > 0)
            continue;
        else {
            polc = n = 0;
            initso = sumnsam-nsam[npops-1];/*initso is the number of sample where the outgroup (or single pop) starts*/
            while((a0 = matrix_pol[j*sumnsam+initso]) == '-' && initso < sumnsam) initso++;
            if(initso >= sumnsam) polc = -1; /*we do not have outgroup (or pop) in this position and we do not count jfd*/
            else {
                n = 1;
                for(popo=initso+1;popo<sumnsam;popo++) {
                    while((a1 = matrix_pol[j*sumnsam+popo]) == '-' && popo < sumnsam) popo++;
                    if(popo >= sumnsam) {
                        break;
                    }
                    else {
                        if(a1 != a0) {
                            polc += 1; /*there is polymorphism within the outgroup (or single pop)*/
                        }
                    }
                    n += 1;
                }
            }
            if(outgroup_presence+force_outgroup == 1 && polc== 0) { /*outgroup*/
                anc = a0;
                for(x=0;x<npops-1;x++) {
                    a0 = a1 = pol = n = 0;
                    inits = initsq1[x];
                    while((a0 = matrix_pol[j*sumnsam+inits]) == '-' && inits < initsq1[x]+nsam[x]) inits += 1;
                    if(inits > initsq1[x]+nsam[x]-1) pol = 0;
                    else {
                        n = 1;
                        for(pop1=inits+1;pop1<initsq1[x]+nsam[x];pop1++) {
                            while((a1 = matrix_pol[j*sumnsam+pop1]) == '-' && pop1 < initsq1[x]+nsam[x]) pop1++;
                            if(pop1 >= initsq1[x]+nsam[x]) break;
                            if(a1 != a0) {
                                pol += 1;
                            }
                            n += 1;
                        }				
                    }
                    if(n) {
                        if(anc != a0) pol = n-pol;
                        jfd[x][j] = (double)pol/(double)n;
                        nfd[x][j] = n;
                        
                        /*count frequency of snps per line*/
                        for(pop1=inits;pop1<initsq1[x]+nsam[x];pop1++) {
                            while((a2 = matrix_pol[j*sumnsam+pop1]) == '-' && pop1 < initsq1[x]+nsam[x]) pop1++;
                            if(pop1 >= initsq1[x]+nsam[x]) break;
                            if(a2 != anc) 
                                statistics[0].linefreq[pop1][pol] += 1;
                            /*else statistics[0].linefreq[pop1][n-pol] += 1;*/
                        }
                    }						
                }
            }
            else { 
                if(polc==0) { /*no outgroup*/
                    anc = a0;
                    for(x=0;x<npops-1;x++) {
                        a0 = a1 = a2 = pol = n = 0;
                        inits = initsq1[x];
                        while((a0 = matrix_pol[j*sumnsam+inits]) == '-' && inits < initsq1[x]+nsam[x]) inits += 1;
                        if(inits > initsq1[x]+nsam[x]-1) pol = 0;
                        else {
                            n = 1;
                            for(pop1=inits+1;pop1<initsq1[x]+nsam[x];pop1++) {
                                while((a1 = matrix_pol[j*sumnsam+pop1]) == '-' && pop1 < initsq1[x]+nsam[x]) pop1++;
                                if(pop1 >= initsq1[x]+nsam[x]) break;
                                if(a1 != a0) {
                                    pol += 1;
                                    a2 = a1;
                                }
                                n += 1;
                            }				
                        }
                        if(n) {
                            if/*(anc != a0)*/(pol > n-pol) pol = n-pol;
                            jfd[x][j] = (double)pol/(double)n;
                            nfd[x][j] = n;
                            /*
                            jfd[x][j] = (pol < n-pol)? (double)pol/(double)n:(double)(n-pol)/(double)n;
                            nfd[x][j] = n;
                            */
                            ancl  = (pol < n-pol)? a0:a2;
                            polv  = (pol < n-pol)? pol:n-pol;

                            /*count frequency of snps per line*/
                            for(pop1=inits;pop1<initsq1[x]+nsam[x];pop1++)
                            {
                                while((a2 = matrix_pol[j*sumnsam+pop1]) == '-' && pop1 < initsq1[x]+nsam[x]) pop1++;
                                if(pop1 >= initsq1[x]+nsam[x]) break;
                                if(a2 != ancl) 
                                    statistics[0].linefreq[pop1][polv] += 1;
                                /*else statistics[0].linefreq[pop1][n-pol] += 1;*/
                            }
                        }						
                    }
                }
            }
        }
    }
    free(initsq1);

	return 1;
}
