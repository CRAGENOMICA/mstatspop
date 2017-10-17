/*
 *  sancestral.c
 *  mstatspop
 *
 *  Created by Sebastian Ramos-Onsins on 16/08/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "sancestral.h"

int calc_sxsfss(int npops, int *nsam, char *matrix_pol,long int *matrix_pos,
		long int length, struct stats *statistics,long int *sites_matrix,int outgroup_presence,int force_outgroup)
{
	int x;
	long int j;
	int sumnsam;
	int ispolmis(char *,int,long int,long int);
	int polc,pop1,pop2,popo,npo;
	int a0,a1;
	int *initsq1,*initsq1a,*initsq2,initso,initsa,*polqa,*polqb;
	int initsb=0;
	long int *y0,*y1,*y2,*y3,*o0,*o1,*o2,*o3;
	long int css_rest;
	
    for(x=0;x<4*npops;x++) {
        statistics[0].Sanc[x] = 0;
    }
    if(length==0) return 1;
    
    initsq1 = (int *)calloc(npops,sizeof(int));
    initsq2 = (int *)calloc(npops,sizeof(int));
    initsq1a= (int *)calloc(npops,sizeof(int));
    polqa   = (int *)calloc(npops,sizeof(int));
    polqb   = (int *)calloc(npops,sizeof(int));
    
    y0   = (long int *)calloc(npops,sizeof(long int));
    y1   = (long int *)calloc(npops,sizeof(long int));
    y2   = (long int *)calloc(npops,sizeof(long int));
    y3   = (long int *)calloc(npops,sizeof(long int));
    o0   = (long int *)calloc(1,sizeof(long int));
    o1   = (long int *)calloc(1,sizeof(long int));
    o2   = (long int *)calloc(1,sizeof(long int));
    o3   = (long int *)calloc(1,sizeof(long int));
    
    sumnsam = 0;
    for(x=0;x<npops;x++) {
        initsq1[x] = sumnsam;
        if(x) initsq2[x] = 0;
        else initsq2[x] = nsam[0];
        sumnsam += nsam[x];
        y0[x]=y1[x]=y2[x]=y3[x]=0;
    }
    o0[0]=o1[0]=o2[0]=0;
    
    a0=a1=0;
    for(j=0;j<length;j++)
    {
        if((ispolmis(matrix_pol,sumnsam,j,length)) > 0) {
            continue;
        }
        else {
            /*outgroup*/
            polc = 0;
            initso = sumnsam-nsam[npops-1];
            while((a0 = matrix_pol[j*sumnsam+initso]) == '-' && initso < sumnsam) initso++;
            if(initso >= sumnsam) {
                polc = -1;
            }
            else {
                for(popo=initso+1;popo<sumnsam;popo++) {
                    while((a1 = matrix_pol[j*sumnsam+popo]) == '-' && popo < sumnsam) popo++;
                    if(popo >= sumnsam) {
                        polc = 0;
                        break;
                    }
                    else {
                        if(a1 != a0) {
                            polc = 1;
                            break;
                        }
                    }
                }				
            }
            /*populations*/
            css_rest = 1;
            for(x=0;x<npops-1;x++) {
                /*check if each group is polymorphic or not*/
                polqa[x] = polqb[x] = 0;
                initsa = initsq1a[x] = initsq1[x];
                while((a0 = matrix_pol[j*sumnsam+initsa]) == '-' && initsa < initsq1[x]+nsam[x]) initsa += 1;
                if(initsa >= initsq1[x]+nsam[x]) {polqa[x] = -1;}
                else {
                    initsq1a[x] = initsa;
                    for(pop1=initsa+1;pop1<initsq1[x]+nsam[x];pop1++) {
                        while((a1 = matrix_pol[j*sumnsam+pop1]) == '-' && pop1 < initsq1[x]+nsam[x]) pop1++;
                        if(pop1 >= initsq1[x]+nsam[x]) break;
                        if(a1 != a0) {
                            polqa[x] = 1;
                            break;
                        }
                    }				
                }
                for(pop1=0;pop1<npops-1;pop1++) {
                    if(pop1 != x) {
                        initsb = initsq2[x];
                        while((a0 = matrix_pol[j*sumnsam+initsb]) == '-' && initsb < sumnsam-nsam[npops-1]) {
                            if(initsb == initsq1[x]) 
                                initsb += nsam[x];
                            else initsb += 1;
                        }
                        if(initsb >= (sumnsam-nsam[npops-1])) polqb[x] = -1;
                        else {
                            /*a0 = nt in inits+0*/
                            for(pop2=initsb+1;pop2<sumnsam-nsam[npops-1];pop2++) {
                                if(pop2 >= initsq1[x] && pop2 < initsq1[x]+nsam[x]) 
                                    pop2 = initsq1[x+1];
                                while((a1 = matrix_pol[j*sumnsam+pop2]) == '-' && pop2 < sumnsam-nsam[npops-1]) {
                                    if(pop2 == initsq1[x]) pop2 += nsam[x];
                                    else pop2++;
                                }
                                if(pop2 >= initsq1[x] && pop2 < initsq1[x]+nsam[x]) 
                                    pop2 = initsq1[x+1];
                                if(pop2 >= (sumnsam-nsam[npops-1])) break;
                                if(a1 != a0) {
                                    polqb[x] = 1;
                                    break;
                                }
                            }				
                        }
                    }
                    if(polqb[x] == 1) break;
                }
                /*define classes*/
                if(polqa[x]==1 && polqb[x]==0 && polc==0) {
                    if((npops == 2) || (matrix_pol[j*sumnsam+initsb] == matrix_pol[j*sumnsam+initso])) {
                        statistics[0].Sanc[x*4+0] += 1;/*Sx1*/
                        sites_matrix[y0[x]*4*npops+4*x+0] = matrix_pos[j];
                        y0[x] += 1;
                        css_rest = 0;
                        continue;
                    }
                    else {
                        if(outgroup_presence+force_outgroup == 1) statistics[0].Sanc[x*4+2] += 1;/*Sx1f2*/
                        if(outgroup_presence+force_outgroup == 0) statistics[0].Sanc[x*4+0] += 1;/*Sx1*/
                        sites_matrix[y2[x]*4*npops+4*x+2] = matrix_pos[j];
                        y2[x] += 1;
                        css_rest = 0;
                        continue;
                    }
                }
                if(polqa[x]==1 && polqb[x]==1 && polc==0) {
                    statistics[0].Sanc[x*4+3] += 1;/*Ssh*/
                    sites_matrix[y3[x]*4*npops+4*x+3] = matrix_pos[j];
                    y3[x] += 1;
                    css_rest = 0;
                    continue;
                }
                if(polqa[x]==0 && polqb[x]==0 && polc==0) {
                    if(outgroup_presence+force_outgroup == 1) {
                        if((npops == 2 && matrix_pol[j*sumnsam+initsa] != matrix_pol[j*sumnsam+initso]) || (matrix_pol[j*sumnsam+initsa] != matrix_pol[j*sumnsam+initso] && matrix_pol[j*sumnsam+initsb] == matrix_pol[j*sumnsam+initso])) {
                            statistics[0].Sanc[x*4+1]+= 1;/*Sf1*/
                            sites_matrix[y1[x]*4*npops+4*x+1] = matrix_pos[j];
                            y1[x] += 1;
                            css_rest = 0;
                            continue;
                        }
                    }
                    if(outgroup_presence+force_outgroup == 0) {
                        if(matrix_pol[j*sumnsam+initsa] != matrix_pol[j*sumnsam+initsb]) {
                            statistics[0].Sanc[x*4+1]+= 1;/*Sf1*/
                            sites_matrix[y1[x]*4*npops+4*x+1] = matrix_pos[j];
                            y1[x] += 1;
                            css_rest = 0;
                            continue;
                        }
                    }
                }
                if(polqa[x]== -1) {
                    /*no information for this population and position, continue*/
                    continue;
                }
                if(polqb[x] == -1) {
                    /*this is like 2 pops with no outgroup: Sx1,Ss,Sf[1 vs outg]: not considered because fixed variants can not be polarized. It is inconsistent with the rest of mutations with 3 pops.*/
                    continue;
                }
            }
            /*for the outgroup: it doesn't matter what is the group of populations(x), we chose pop <- npops-2 if it is not -1*/
            npo=npops-2;
            while(npo > 0 && (polqa[npo] == -1 || polqb[npo] == -1)) npo--;
            if(polc == 1) {
                if(npops-!(outgroup_presence+force_outgroup) > 1) {
                    if(polqa[npo] == 0 && polqb[npo] == 0) {
                        if((npops == 2) || matrix_pol[j*sumnsam+initsb] == matrix_pol[j*sumnsam+initsq1a[npo]]) {
                            statistics[0].Sanc[(npops-1)*4+0] += 1;/*Sxo*/
                            sites_matrix[o0[0]*4*npops+4*(npops-1)+0] = matrix_pos[j];
                            o0[0]++;
                            css_rest = 0;
                            continue;
                        }
                        else {
                           statistics[0].Sanc[(npops-1)*4+2] += 1; /*Sanco*/
                           sites_matrix[o2[0]*4*npops+4*(npops-1)+2] = matrix_pos[j];
                           o2[0]++;
                            css_rest = 0;
                            continue;
                        }
                    }
                    else {
                        if(polqa[npo] != -1 && polqb[npo] != -1) {
                            statistics[0].Sanc[(npops-1)*4+2] += 1; /*Sanco*/
                            sites_matrix[o2[0]*4*npops+4*(npops-1)+2] = matrix_pos[j];
                            o2[0]++;
                            css_rest = 0;
                            continue;
                        }
                        else {
                            if(polqa[npo] == 0 || polqb[npo] == 0) {
                                statistics[0].Sanc[(npops-1)*4+0] += 1;/*Sxo*/
                                sites_matrix[o0[0]*4*npops+4*(npops-1)+0] = matrix_pos[j];
                                o0[0]++;
                                css_rest = 0;
                                continue;
                            }
                            else {
                                statistics[0].Sanc[(npops-1)*4+2] += 1; /*Sanco*/
                                sites_matrix[o2[0]*4*npops+4*(npops-1)+2] = matrix_pos[j];
                                o2[0]++;
                                css_rest = 0;
                                continue;
                            }
                        }
                    }
                }
                else  { 
                    statistics[0].Sanc[(npops-1)*4+0] += 1;/*Sx1*/
                    sites_matrix[y0[0]*4*npops+4*(npops-1)+0] = matrix_pos[j];
                    y0[0] += 1;
                    css_rest = 0;
                    continue;
                }
            }
            if(npops > 1) { 
                if(polqa[npo]==0 && polqb[npo]==0 && polc==0) {
                    if((npops == 2 && matrix_pol[j*sumnsam+initsq1a[npo]] != matrix_pol[j*sumnsam+initso]) || (matrix_pol[j*sumnsam+initsb] == matrix_pol[j*sumnsam+initsq1a[npo]] &&  matrix_pol[j*sumnsam+initsq1a[npo]] != matrix_pol[j*sumnsam+initso])) {
                        statistics[0].Sanc[(npops-1)*4+1] += 1;/*Sfo*/
                        sites_matrix[o1[0]*4*npops+4*(npops-1)+1] = matrix_pos[j];
                        o1[0]++;
                        css_rest = 0;
                        continue;
                    }
                }
            }
            if(polc == -1) {
                /*this is like 2 pops with no outgroup: Sx,Ss,Sf[1 vs 2]: not considered because fixed variants can not be polarized. It is inconsistent with the rest of mutations with 3 pops.*/
                continue;
            }
            if(css_rest == 1) {
                statistics[0].Sanc[(npops-1)*4+3] += 1;/*Ssrest fixed between several pops: NOT COUNTED  BEFORE!*/
                sites_matrix[o3[0]*4*npops+4*(npops-1)+3] = matrix_pos[j];
                o3[0]++;
                css_rest = 0;
            }
        }
    }
    
    free(initsq1);
    free(initsq2);
    free(initsq1a);
    free(polqa);
    free(polqb);

    for(x=0;x<npops-1;x++){
        sites_matrix[y0[x]*4*npops+4*x+0] = 0;
        sites_matrix[y1[x]*4*npops+4*x+1] = 0;
        sites_matrix[y2[x]*4*npops+4*x+2] = 0;
        sites_matrix[y3[x]*4*npops+4*x+3] = 0;
    }
    if(npops >1) {
        sites_matrix[o0[0]*4*npops+4*(npops-1)+0] = 0;
        sites_matrix[o1[0]*4*npops+4*(npops-1)+1] = 0;
        sites_matrix[o2[0]*4*npops+4*(npops-1)+2] = 0;
    }
    free(y0);
    free(y1);
    free(y2);
    free(y3);
    free(o0);
    free(o1);
    free(o2);
    free(o3);

	return 1;
}
		   
int ispolmis(char *matrix_pol,int nsam, long int j,long int length)
{
	int x,y;
	int a0,a1;
	
	y=0;
	while((a0 = matrix_pol[j*nsam+y]) == '-' && y < nsam) {
		y++;
	}
	if(y>=nsam-1) {
		return 2; /*all are missings or monomorphic*/
	}

	for(x=(y+1);x<nsam;x++)
	{
		while( (a1 = matrix_pol[j*nsam+x]) == '-' && x < nsam)
			x++;

		if(a0 != a1) {
			return 0; /*polymorphic*/
		}
	}
	
	return 1; /*monomorphic*/
}

