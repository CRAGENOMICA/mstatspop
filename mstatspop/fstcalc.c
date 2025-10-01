/*
 *  fstcalc.c
 *  mstatspop
 *
 *  Created by Sebastian Ramos-Onsins on 17/08/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "fstcalc.h"


/*FST = 1 - PIWITHIN/PIAMONG*/

int calc_piwpiafst(int flag, int formatfile,int npops, int *nsam, char *matrix_pol, long int length, struct stats *statistics, long int *matrix_sv,int outgroup_presence, int force_outgroup)
{
    int x,y,z;
    long int j;
    int pop1,pop2;
    int **freq,**freq1;
    int sumnsam,*initsq1,inits,max,initso,a0;
    double pia,piw,piT,*pia1all;
    double /*spiw,*/spiwn,spia,spian;
    /*double *mean_nsam,**mean_nsam_amng;*/
    /*long int *n_nsam,**n_nsam_amng;*/
    int ncw,nca;
    double tn93(double,double,double,double, double, double, double,double);
    double gT,gC,gA,gG,P1,P2,Q;
    double total_effective_length,sumg;
    double *piwallno1;
    double *maxlength, maxlength_total;
    
    /*init*/
    z = 0;
    statistics[0].fstALL = -10000;
    for(pop1=0;pop1<npops/*-1*/;pop1++) {
        statistics[0].piw[pop1] = 0.;
        statistics[0].K[pop1] = 0.;
        statistics[0].piwHKY[pop1] = -10000;
        statistics[0].thetaTHKY[pop1] = -10000;
        statistics[0].KHKY[pop1] = -10000;
        statistics[0].fst1all[pop1] = -10000;
        statistics[0].sv[pop1][pop1][0] = 0.;
        statistics[0].sv[pop1][pop1][1] = 0.;
        for(pop2=pop1+1;pop2<npops;pop2++) {
            statistics[0].pia[z] = 0.;
            statistics[0].piT[z] = 0.;
            statistics[0].piaHKY[z] = -10000;
            statistics[0].piTHKY[z] = -10000;
            statistics[0].fst[z] = -10000;
            statistics[0].fstHKY[z] = -10000;
            statistics[0].sv[pop1][pop2][0] = 0.;
            statistics[0].sv[pop1][pop2][1] = 0.;
            z++;
        }
    }
    
    if(length==0) return 1;
    
    initsq1 = (int *)calloc(npops,sizeof(int));
    pia1all = (double *)calloc(npops,sizeof(double));
    piwallno1 = (double *)calloc(npops,sizeof(double));
    freq = (int **)calloc(npops,sizeof(int *));
    for(x=0;x<npops;x++) freq[x] = (int *)calloc(4,sizeof(int));
    freq1 = (int **)calloc(npops,sizeof(int *));
    for(x=0;x<npops;x++) freq1[x] = (int *)calloc(4,sizeof(int));
    /*mean_nsam = (double *)calloc(npops,sizeof(double));
     mean_nsam_amng = (double **)calloc(npops,sizeof(double *));
     for(x=0;x<npops;x++) mean_nsam_amng[x] = (double *)calloc(npops,sizeof(double));
     n_nsam = (long int *)calloc(npops,sizeof(long int));
     n_nsam_amng = (long int **)calloc(npops,sizeof(long int *));
     for(x=0;x<npops;x++) n_nsam_amng[x] = (long int *)calloc(npops,sizeof(long int));*/
    maxlength = (double *)calloc(npops,sizeof(double));
    
    sumnsam = 0;
    for(x=0;x<npops;x++) {
        initsq1[x] = sumnsam;
        sumnsam += nsam[x];
    }
    
    /*NOTE THE HAPLOTYPE VALUES GO FROM 0 to npops-1 BUT FREQUENCY GOES FROM 0 to npops!!!!!!*/
    /*
     printf("\n");
     for(y=0;y<sumnsam;y++) {
     for(j=0;j<length;j++) {
     printf("%c",matrix_pol[j*sumnsam+y]);
     }
     printf("\n");
     }
     printf("\n");
     */
    for(j=0;j<length;j++) {
        /*eliminate if no outgroup*/
        initso = sumnsam-nsam[npops-1];
        while((a0 = matrix_pol[j*sumnsam+initso]) == '-' && initso < sumnsam) initso++;
        if(initso >= sumnsam) {
            continue;
        }
        /*calculate the frequencies of each population for each position*/
        for(pop1=0;pop1<npops-0;pop1++) {
            inits = initsq1[pop1];
            max   = initsq1[pop1]+nsam[pop1];
            freq[pop1][0]=freq[pop1][1]=freq[pop1][2]=freq[pop1][3]=0;
            for(y=inits;y<max;y++) {
                if(matrix_pol[j*sumnsam+y] == '0') {freq[pop1][1] += 1;freq[pop1][0] += 1;}
                if(matrix_pol[j*sumnsam+y] == '1') {freq[pop1][2] += 1;freq[pop1][0] += 1;}
                if(matrix_pol[j*sumnsam+y] == '-') {freq[pop1][3] += 1;}
            }
            if(nsam[pop1] > 1 && freq[pop1][0] > 1) {
                piw = (double)(freq[pop1][1]*freq[pop1][2])/(double)(freq[pop1][0]*(freq[pop1][0]-1.0)/2.0);
                statistics[0].piw[pop1] += piw;
                if(flag == 1) if(matrix_sv[j] == 1) statistics[0].sv[pop1][pop1][0] += piw;
                if(flag == 1) if(matrix_sv[j] == 2) statistics[0].sv[pop1][pop1][1] += piw;
                /*mean_nsam[pop1] += (double)freq[pop1][0];
                 n_nsam[pop1] += 1;*/
            }
            /*including new code for Fst1all:*/
            y =0;
            freq1[pop1][0]=freq1[pop1][1]=freq1[pop1][2]=freq1[pop1][3]=0;
            while(y < initso) {
                if(y<inits || y >= max) {
                    if(matrix_pol[j*sumnsam+y] == '0') {freq1[pop1][1] += 1;freq1[pop1][0] += 1;}
                    if(matrix_pol[j*sumnsam+y] == '1') {freq1[pop1][2] += 1;freq1[pop1][0] += 1;}
                    if(matrix_pol[j*sumnsam+y] == '-') {freq1[pop1][3] += 1;}
                }
                y++;
            }
            if(initso - nsam[pop1] > 1 && freq1[pop1][0] > 1) {
                piw = (double)(freq1[pop1][1]*freq1[pop1][2])/(double)(freq1[pop1][0]*(freq1[pop1][0]-1.0)/2.0);
                piwallno1[pop1] += piw;
            }
            if(initso - nsam[pop1] > 0 && freq[pop1][0] > 0 && freq1[pop1][0] > 0) {
                pia = (double)(freq[pop1][1]*freq1[pop1][2] + freq[pop1][2]*freq1[pop1][1])/(double)(freq[pop1][0]*freq1[pop1][0]);
                pia1all[pop1] += pia;
            }
            /**/
        }
        z=0;
        for(pop1=0;pop1<npops-1;pop1++) {
            for(pop2=pop1+1;pop2<npops-0;pop2++) {
                if((freq[pop1][0]*freq[pop2][0]) > 0) {
                    pia = (double)(freq[pop1][1]*freq[pop2][2] + freq[pop1][2]*freq[pop2][1])/(double)(freq[pop1][0]*freq[pop2][0]);
                    statistics[0].pia[z] += pia;
                    if(flag == 1) if(matrix_sv[j] == 1) statistics[0].sv[pop1][pop2][0] += pia;
                    if(flag == 1) if(matrix_sv[j] == 2) statistics[0].sv[pop1][pop2][1] += pia;
                    if(pop2 == npops-1)
                        statistics[0].K[pop1] += pia;
                    /*we do not show these results. we do not consider those values that are not in Fst*//*
                                                                                                          if((pop2 < npops-1 || npops == 2) && nsam[pop1]>1 && nsam[pop2]>1) {
                                                                                                          pia1all[pop1] += pia;
                                                                                                          pia1all[pop2] += pia;
                                                                                                          }*//*
                                                                                                              mean_nsam_amng[pop1][pop2] += (double)freq[pop1][0] * (double)freq[pop2][0];
                                                                                                              n_nsam_amng[pop1][pop2] += 1;*/
                }
                z++;
            }
        }
        z = 0;
        for(pop1=0;pop1<npops-1;pop1++) {
            for(pop2=pop1+1;pop2<npops-0;pop2++) {
                if((nsam[pop1]+nsam[pop2]) > 1 && (freq[pop1][0]+freq[pop2][0]) > 1) {
                    piT = (double)((freq[pop1][1]+freq[pop2][1])*(freq[pop1][2]+freq[pop2][2]))/(double)((freq[pop1][0]+freq[pop2][0])*((freq[pop1][0]+freq[pop2][0])-1.0)/2.0);
                    statistics[0].piT[z] += piT;
                    if(flag == 1) if(matrix_sv[j] == 1) statistics[0].svT[pop1][pop2][0] += piT;
                    if(flag == 1) if(matrix_sv[j] == 2) statistics[0].svT[pop1][pop2][1] += piT;
                }
                z++;
            }
        }
    }
    
    maxlength_total = 0.;
    for(pop1=0;pop1<npops-1;pop1++) { maxlength[pop1] = 0.;}
    
    for(pop1=0;pop1<npops-1;pop1++) {
        for(pop2=pop1+1;pop2<npops-0;pop2++) {
            if(pop1 != pop2) {
                if(outgroup_presence + force_outgroup == 1) {
                    if(maxlength[pop1] < statistics[0].lengthamng_outg[pop1][pop2])
                        maxlength[pop1] = statistics[0].lengthamng_outg[pop1][pop2];
                    if(maxlength[pop2] < statistics[0].lengthamng_outg[pop1][pop2])
                        maxlength[pop2] = statistics[0].lengthamng_outg[pop1][pop2];
                }
                else {
                    if(maxlength[pop1] < statistics[0].lengthamng[pop1][pop2])
                        maxlength[pop1] = statistics[0].lengthamng[pop1][pop2];
                    if(maxlength[pop2] < statistics[0].lengthamng[pop1][pop2])
                        maxlength[pop2] = statistics[0].lengthamng[pop1][pop2];
                }
            }
        }
        if(maxlength[pop1] > maxlength_total)
            maxlength_total = maxlength[pop1];
    }
    
    z=0;
    for(pop1=0;pop1<npops-1;pop1++)
    {
        for(pop2=pop1+1;pop2<npops-0;pop2++)
        {
            if(outgroup_presence + force_outgroup == 1)
                total_effective_length = (double)statistics[0].lengthamng_outg[pop1][pop2];
            else
                total_effective_length = (double)statistics[0].lengthamng[pop1][pop2];
            statistics[0].piant[z] = statistics[0].pia[z]/(double)(total_effective_length/*statistics[0].total_length*/);/*incorrect...?*/
            statistics[0].piTnt[z] = statistics[0].piT[z]/(double)(/**/total_effective_length/*statistics[0].total_length*/);/*incorrect...?*/
            z++;
        }
    }
    
    /*HKY correction and piant*/
    if(flag == 1 && outgroup_presence + force_outgroup == 1) {
        if(formatfile == 0 || formatfile == 3) {
            z=0;
            for(pop1=0;pop1<npops-1;pop1++)
            {
                if((double)statistics[0].length2[pop1] == 0) {
                    statistics[0].piwHKY[pop1] = -10000;
                    statistics[0].thetaTHKY[pop1] = -10000;
                    for(pop2=pop1+1;pop2<npops-0;pop2++) {
                        statistics[0].piaHKY[z] = -10000;
                        statistics[0].piTHKY[z] = -10000;
                        if(pop2 ==npops-1) 
                            statistics[0].KHKY[pop1] = -10000;
                        /*statistics[0].piant[z] = -10000;*/
                        /*statistics[0].piTnt[z] = -10000;*/
                        z++;
                    }
                    continue;
                }
                /**/total_effective_length = (double)statistics[0].length2[pop1];/**/
                sumg = (double) statistics[0].tcga[pop1][0] + statistics[0].tcga[pop1][1] + statistics[0].tcga[pop1][2] + statistics[0].tcga[pop1][3];
                gT = statistics[0].tcga[pop1][0] /sumg * total_effective_length;
                gC = statistics[0].tcga[pop1][1] /sumg * total_effective_length;
                gG = statistics[0].tcga[pop1][2] /sumg * total_effective_length;
                gA = statistics[0].tcga[pop1][3] /sumg * total_effective_length;
                P1 = statistics[0].sv[pop1][pop1][0] * (gA*gG/(gA*gG + gT*gC));
                P2 = statistics[0].sv[pop1][pop1][0] * (gT*gC/(gA*gG + gT*gC));		
                Q =  statistics[0].sv[pop1][pop1][1];
                statistics[0].piwHKY[pop1] = tn93(gT,gC,gG,gA,P1,P2,Q,/*statistics[0].total_length*/ /**/total_effective_length/**/);
                statistics[0].thetaTHKY[pop1] = statistics[0].piwHKY[pop1];
                for(pop2=pop1+1;pop2<npops-0;pop2++)
                {
                    if((double)statistics[0].lengthamng_outg[pop1][pop2] == 0) {
                        statistics[0].piaHKY[z] = -10000;
                        statistics[0].piTHKY[z] = -10000;
                        if(pop2 ==npops-1) 
                            statistics[0].KHKY[pop1] = -10000;
                        /*statistics[0].piant[z] = -10000;*/
                        /*statistics[0].piTnt[z] = -10000;*/
                        z++;
                        continue;
                    }
                    total_effective_length = (double)statistics[0].lengthamng_outg[pop1][pop2];
                    sumg = (double) statistics[0].tcga[pop1][0] + statistics[0].tcga[pop1][1] + statistics[0].tcga[pop1][2] + statistics[0].tcga[pop1][3] + 
                    statistics[0].tcga[pop2][0] + statistics[0].tcga[pop2][1] + statistics[0].tcga[pop2][2] + statistics[0].tcga[pop2][3];
                    gT = (statistics[0].tcga[pop1][0] + statistics[0].tcga[pop2][0]) /sumg * total_effective_length;
                    gC = (statistics[0].tcga[pop1][1] + statistics[0].tcga[pop2][1]) /sumg * total_effective_length;
                    gG = (statistics[0].tcga[pop1][2] + statistics[0].tcga[pop2][2]) /sumg * total_effective_length;
                    gA = (statistics[0].tcga[pop1][3] + statistics[0].tcga[pop2][3]) /sumg * total_effective_length;		
                    P1 = statistics[0].sv[pop1][pop2][0] * (gA*gG/(gA*gG + gT*gC));
                    P2 = statistics[0].sv[pop1][pop2][0] * (gT*gC/(gA*gG + gT*gC));		
                    Q =  statistics[0].sv[pop1][pop2][1];
                    statistics[0].piaHKY[z] = tn93(gT,gC,gG,gA,P1,P2,Q,/*statistics[0].total_length*/ /**/total_effective_length/**/);/*it is incorrect in case missing values are considered...*/
                    /*piant*/  
                    //statistics[0].piant[z] = statistics[0].pia[z]/(double)(total_effective_length/*statistics[0].total_length*/);/*incorrect*/
                    
                    P1 = statistics[0].svT[pop1][pop2][0] * (gA*gG/(gA*gG + gT*gC));
                    P2 = statistics[0].svT[pop1][pop2][0] * (gT*gC/(gA*gG + gT*gC));		
                    Q =  statistics[0].svT[pop1][pop2][1];		
                    statistics[0].piTHKY[z] = tn93(gT,gC,gG,gA,P1,P2,Q,/*statistics[0].total_length*/ /**/total_effective_length/**/);/*it is incorrect in case missing values are considered...*/
                    /*piTnt*/  
                    //statistics[0].piTnt[z] = statistics[0].piT[z]/(double)(/**/total_effective_length/*statistics[0].total_length*/);/*incorrect*/
                    
                    if(pop2 ==npops-1) 
                        statistics[0].KHKY[pop1] = statistics[0].piaHKY[z];
                    z++;					
                }
            }
        }
        if(formatfile ==1 || formatfile == 2) {
            z=0;
            /*
             for(pop1=0;pop1<npops;pop1++) {
             n1 = ceil(statistics[0].length[pop1]/statistics[0].total_length - 1E-5/);
             gT = statistics[0].tcga[pop1][0];
             gC = statistics[0].tcga[pop1][1];
             gG = statistics[0].tcga[pop1][2];
             gA = statistics[0].tcga[pop1][3];		
             P1 = statistics[0].sv[pop1][pop1][0] * (gA*gG/(gA*gG + gT*gC));
             P2 = statistics[0].sv[pop1][pop1][0] * (gT*gC/(gA*gG + gT*gC));		
             Q = statistics[0].sv[pop1][pop1][1];		
             statistics[0].piwHKY[pop1] = tn93(gT,gC,gG,gA,P1,P2,Q,statistics[0].total_length);
             statistics[0].thetaTHKY[pop1] = statistics[0].piwHKY[pop1];
             for(pop2=pop1+1;pop2<npops;pop2++) {
             n2 = ceil(statistics[0].length[pop2]/statistics[0].total_length- 1E-5);
             gT = statistics[0].tcga[pop1][0] + statistics[0].tcga[pop2][0];
             gC = statistics[0].tcga[pop1][1] + statistics[0].tcga[pop2][1];
             gG = statistics[0].tcga[pop1][2] + statistics[0].tcga[pop2][2];
             gA = statistics[0].tcga[pop1][3] + statistics[0].tcga[pop2][3];		
             P1 = statistics[0].sv[pop1][pop2][0] * (gA*gG/(gA*gG + gT*gC));
             P2 = statistics[0].sv[pop1][pop2][0] * (gT*gC/(gA*gG + gT*gC));		
             Q = statistics[0].sv[pop1][pop2][1];		
             statistics[0].piaHKY[z] = tn93(gT,gC,gG,gA,P1,P2,Q,statistics[0].total_length);
             if(pop2 ==npops-1) 
             statistics[0].KHKY[pop1] = statistics[0].piaHKY[z];
             statistics[0].piant[z] = statistics[0].pia[z]/(double)(statistics[0].total_length);
             z++;					
             }
             }
             */
        }
    }
    else {
        if(formatfile ==1 || formatfile == 2) {
            z=0;
            for(pop1=0;pop1<npops-1;pop1++) {
                statistics[0].piwHKY[pop1] = -10000;
                statistics[0].thetaTHKY[pop1] = -10000;
                for(pop2=pop1+1;pop2<npops-0;pop2++) {
                    statistics[0].piaHKY[z] = -10000;
                    statistics[0].piTHKY[z] = -10000;
                    if(pop2 ==npops-1) 
                        statistics[0].KHKY[pop1] = -10000;
                    statistics[0].piant[z] = -10000;
                    statistics[0].piTnt[z] = -10000;
                    z++;
                }
            }
        }
    }
    
    /*Fsts calculated excluding the outgroup (outgroup always the last pop even in case it is not defined)*/
    /*Fst pair-pair*/
    z=0;
    for(pop1=0;pop1<npops-1;pop1++) {
        for(pop2=pop1+1;pop2<npops-0;pop2++) {
            if(outgroup_presence + force_outgroup == 1) {
                if(statistics[0].pia[z] && nsam[pop1] > 1 && nsam[pop2] > 1 &&
                   (double)statistics[0].length2[pop1] > 0 && 
                   (double)statistics[0].length2[pop2] > 0 &&
                   (double)statistics[0].lengthamng_outg[pop1][pop2] > 0) {
                    statistics[0].fst[z] = 1.0 - ((statistics[0].piw[pop1]/(double)statistics[0].length2[pop1] + 
                                                   statistics[0].piw[pop2]/(double)statistics[0].length2[pop2])/2.0)/
                    (statistics[0].pia[z]/(double)statistics[0].lengthamng_outg[pop1][pop2]);
                } else {
                    statistics[0].fst[z] = -10000;
                }
                if(flag == 1 &&
                   (double)statistics[0].length2[pop1] > 0 && 
                   (double)statistics[0].length2[pop2] > 0 &&
                   (double)statistics[0].lengthamng_outg[pop1][pop2] > 0) {
                    if(statistics[0].piaHKY[z] > 0 && nsam[pop1] > 1 && nsam[pop2] > 1) {
                        statistics[0].fstHKY[z] = 1.0 - ((statistics[0].piwHKY[pop1]/(double)statistics[0].length2[pop1] + 
                                                          statistics[0].piwHKY[pop2]/(double)statistics[0].length2[pop2])/2.0)/
                        (statistics[0].piaHKY[z]/(double)statistics[0].lengthamng_outg[pop1][pop2]);
                    } else {
                        statistics[0].fstHKY[z] = -10000;
                    }
                }
                else if(formatfile ==1 || formatfile == 2) statistics[0].fstHKY[z] = -10000;
            }
            else {
                if(statistics[0].pia[z] && nsam[pop1] > 1 && nsam[pop2] > 1 &&
                   (double)statistics[0].length2[pop1] > 0 &&
                   (double)statistics[0].length2[pop2] > 0 &&
                   (double)statistics[0].lengthamng[pop1][pop2] > 0) {
                    statistics[0].fst[z] = 1.0 - (
                                                  (statistics[0].piw[pop1]/(double)statistics[0].length2[pop1] +
                                                   statistics[0].piw[pop2]/(double)statistics[0].length2[pop2])/2.0
                                                  )/(statistics[0].pia[z]/(double)statistics[0].lengthamng[pop1][pop2]
                                                     );
                } else {
                    statistics[0].fst[z] = -10000;
                }
                statistics[0].fstHKY[z] = -10000;
            }
            z++;
        }
    }
    
    /*Fstall*/
    if(npops > 2) {
        z = 0;
        /*spiw = */spiwn = spia = spian = 0.;
        ncw = nca = 0;
        for(pop1=0;pop1<npops-1;pop1++) {
            if(nsam[pop1] > 1 && (double)statistics[0].length2[pop1] > 0) {
                spiwn += statistics[0].piw[pop1]/(double)statistics[0].length2[pop1];
                /*spiw  += statistics[0].piw[pop1]/(double)statistics[0].length2[pop1] * mean_nsam[pop1] * (mean_nsam[pop1]-1.0) / 2.0;*/
                ncw += 1;
            }
            for(pop2=pop1+1;pop2<npops-0;pop2++) {
                if(nsam[pop1] + nsam[pop2] > 1 && pop2 < npops-1 && (double)statistics[0].lengthamng_outg[pop1][pop2] > 0) {
                    spian += statistics[0].pia[z]/(double)statistics[0].lengthamng_outg[pop1][pop2];
                    /*spia  += statistics[0].pia[z]/(double)statistics[0].lengthamng_outg[pop1][pop2] * mean_nsam_amng[pop1][pop2];*/
                    nca += 1;
                }
                z++;
            }
        }
        for(pop1=0;pop1<npops-1;pop1++) {
            if(nsam[pop1] > 1 && /*spiw+*//*spian*/pia1all[pop1] > 0.) {
                /*statistics[0].fst1all[pop1] = 1.0 - (statistics[0].piw[pop1]/((spiw+spia)/((sumnsam-nsam[npops-1])*(sumnsam-nsam[npops-1]-1.0)/2.0)));*/
                /*statistics[0].fst1all[pop1] = 1.0 - (statistics[0].piw[pop1]/((spia)/((sumnsam-nsam[npops-1])*(sumnsam-nsam[npops-1]-1.0)/2.0)));*/
                /*statistics[0].fst1all[pop1] = 1.0 - ((statistics[0].piw[pop1]/(double)statistics[0].length2[pop1])/((spian)/(double)nca));*/
                statistics[0].fst1all[pop1] =  1.0 - ((statistics[0].piw[pop1]/(double)statistics[0].length2[pop1] +
                                                       piwallno1[pop1]/maxlength[pop1])/2.0) /
                (pia1all[pop1]/maxlength_total);
                /*printf("%.5f\t%.5f\t",(statistics[0].piw[pop1]/(double)statistics[0].length2[pop1]),((spian)/(double)nca));*/
            } else {
                statistics[0].fst1all[pop1] = -10000;
            }
        }
        /*printf("\n");*/
        /*calculate fstALL*/
        if((/*spiw+*/spian) > 0. && nca>0 && ncw>0) {
            /*statistics[0].fstALL = 1.0 - ((spiwn/(double)ncw)/((spiw+spia)/((sumnsam-nsam[npops-1])*(sumnsam-nsam[npops-1]-1.0)/2.0)));*/
            /*statistics[0].fstALL = 1.0 - ((spiwn/(double)ncw)/((spia)/((sumnsam-nsam[npops-1])*(sumnsam-nsam[npops-1]-1.0)/2.0)));*/
            statistics[0].fstALL = 1.0 - ((spiwn/(double)ncw)/((spian)/(double)nca));
        } else {
            statistics[0].fstALL = -10000;
        }
    }
    
    free(initsq1);
    free(pia1all);
    for(x=0;x<npops;x++) free(freq[x]);
    free(freq);
    for(x=0;x<npops;x++) free(freq1[x]);
    free(freq1);
    /*free(mean_nsam);
     for(x=0;x<npops;x++) free(mean_nsam_amng[x]);
     free(mean_nsam_amng);
     free(n_nsam);
     for(x=0;x<npops;x++) free(n_nsam_amng[x]);
     free(n_nsam_amng);*/
    free(piwallno1);
    free(maxlength);
    
    return 1;
}

int calc_hwhafsth(int npops, int *nsam, char *matrix_pol,long int length, struct stats *statistics)
{
    int x,y,z,yy,w;
    long int j,k;
    int flag;
    int pop1,pop2;
    int sumnsam,*initsq1,inits,max;
    double *hapa1all,*hapw1all,*hapT1all,*nsam_allpop1/*,hapw,spiw,spiwn,spia*/;
    char *matrix_hap;
    /*int nn,np;*/
    int **sfreqh;
    char **valh;
    double hs,ht,xi,nm,Js,Hs,Dm,Dkl,Ht;
    int nh,nhp;	
    double shapw,shapa;
    int ncw,nca;
    
    /*init*/
    z = 0;
    statistics[0].fsthALL = -10000;
    statistics[0].GstALL = -10000;
    for(pop1=0;pop1<npops;pop1++) {
        statistics[0].nhpop[pop1] = 0;
        statistics[0].hapw[pop1] = 0.;
        statistics[0].fsth[pop1] = -10000;
        statistics[0].fsth1all[pop1] = -10000;
        for(pop2=pop1+1;pop2<npops;pop2++) {
            statistics[0].hapa[z] = 0.;
            statistics[0].hapT[z] = 0.;
            statistics[0].Gst[z] = -10000;
            z++;
        }
    }
    if(length==0) return 1;
    
    initsq1  = (int *)calloc(npops,sizeof(int));
    
    hapa1all = (double *)calloc(npops,sizeof(double));
    hapT1all = (double *)calloc(npops,sizeof(double));
    hapw1all = (double *)calloc(npops,sizeof(double));
    nsam_allpop1 = (double *)calloc(npops,sizeof(double));
    
    sfreqh   = (int **)calloc(npops,sizeof(int *));
    
    sumnsam = 0;
    for(x=0;x<npops;x++) {
        initsq1[x] = sumnsam;
        sumnsam += nsam[x];
        statistics[0].nhpop[x] = 0;
        statistics[0].hapw[x] = 0.;
    }	
    
    z=0;
    for(pop1=0;pop1<npops-1;pop1++) {
        for(pop2=pop1+1;pop2<npops-0;pop2++) {
            statistics[0].hapa[z] = 0.;
            statistics[0].hapT[z] = 0.;
            z++;
        }
    }
    
    for(x=0;x<npops;x++)
        sfreqh[x] = (int *)calloc(sumnsam,sizeof(int));
    valh = (char **)calloc(sumnsam,sizeof(char *)); /*pointer to each diff haplotype*/	
    
    matrix_hap = (char *)calloc(sumnsam*length,sizeof(char));
    
    /*NOTE THE HAPLOTYPE VALUES GOES FROM 0 to npops-1 BUT FREQUENCY GOES FROM 0 to npops!!!!!!*/
    /*obtain haplotypes in matrix_hap*/
    k = 0;
    for(pop1=0;pop1<npops-0;pop1++) {
        k = 0;
        for(j=0;j<length;j++) {
            inits = initsq1[pop1];
            max   = initsq1[pop1]+nsam[pop1];
            flag = 0;
            for(y=inits;y<max;y++) {
                if(matrix_pol[j*sumnsam+y] != '-') {
                    matrix_hap[y*length+k] = matrix_pol[j*sumnsam+y];
                }
                else {
                    flag = 1; 
                    for(yy=inits;yy<=y;yy++) matrix_hap[yy*length+k] = (char)0;
                    break;
                }
            }
            if(flag == 0) k++; /*number of effective positions in the haplotype sequence*/
        }
    }
    
    /*calculate hapw*/
    nh = 0;
    for(pop1=0;pop1<npops-0;pop1++) {
        statistics[0].hapw[pop1] = 0.;
        inits = initsq1[pop1];
        max   = initsq1[pop1]+nsam[pop1];
        for(y=inits;y<max/*-1*/;y++) {
            /* calculate freq hap*/
            if(y==inits && pop1 == 0) {
                sfreqh[pop1][nh] += 1;
                if(sfreqh[pop1][nh] == 1) statistics[0].nhpop[pop1] += 1;
                valh[nh] = matrix_hap+y*length;
                nh++;
            }
            else {
                w = 0;
                for(z=nh-1;z>=0;z--) {
                    if((memcmp(matrix_hap+y*length,valh[z],k) == 0)) {
                        sfreqh[pop1][z] += 1;
                        if(sfreqh[pop1][z] == 1) statistics[0].nhpop[pop1] += 1;
                        w = 1;
                        break;
                    }
                }
                if(w==0) {
                    sfreqh[pop1][nh] += 1;
                    if(sfreqh[pop1][nh] == 1) statistics[0].nhpop[pop1] += 1;
                    valh[nh] = matrix_hap+y*length;
                    nh++;
                }
            }
            /*end calculation freq hap*/
            for(z=y+1;z<max;z++) {
                if(memcmp(matrix_hap+y*length,matrix_hap+z*length,k) != 0) {
                    statistics[0].hapw[pop1] += 1;
                }
            }
            for(z=0;z<nh;z++) 
                statistics[0].freqh[pop1][z] = sfreqh[pop1][z];
        }
        statistics[0].hapw[pop1] /= (nsam[pop1]*(nsam[pop1]-1.0)/2.0);
        
        /*calculate FstHallpop1*/
        hapw1all[pop1] = 0.;
        nsam_allpop1[pop1] = 0.;
        for(y=0;y<sumnsam-nsam[npops-1];y++) {
            if(y<inits || y >= max) {
                nsam_allpop1[pop1] += 1;
                for(z=y+1;z<sumnsam-nsam[npops-1];z++) {
                    if(z<inits || z >= max) {
                        if(memcmp(matrix_hap+y*length,matrix_hap+z*length,k) != 0)
                            hapw1all[pop1] += 1;
                    }
                }
            }
        }
        hapw1all[pop1] /= (nsam_allpop1[pop1]*(nsam_allpop1[pop1]-1.0)/2.0);
        
        for(y=inits;y<max;y++) {
            for(z=0;z<sumnsam-nsam[npops-1];z++) {
                if(z<inits || z >= max) {
                    if(memcmp(matrix_hap+y*length,matrix_hap+z*length,k) != 0)
                        hapa1all[pop1] += 1;
                }
            }
        }
        hapa1all[pop1] /= (nsam[pop1]*nsam_allpop1[pop1]);
        /**/
    }
    statistics[0].nh = nh;
    
    /*calculate hapa and hapT*/
    z=0;
    for(pop1=0;pop1<npops-1;pop1++) {
        for(pop2=pop1+1;pop2<npops-0;pop2++) {
            /*REDUNDANT?*/
            /*k = 0;
             for(j=0;j<length;j++) {
             inits = initsq1[pop1];
             max   = initsq1[pop1]+nsam[pop1];
             flag = 0;
             for(y=inits;y<max;y++) {
             if(matrix_pol[j*sumnsam+y] != '-') {
             matrix_hap[y*length+k] = matrix_pol[j*sumnsam+y];
             }
             else {
             flag = 1; 
             for(yy=inits;yy<=y;yy++) matrix_hap[yy*length+k] = (char)0;
             break;
             }
             }
             if(flag == 1) continue;
             inits = initsq1[pop2];
             max   = initsq1[pop2]+nsam[pop2];
             for(y=inits;y<max;y++) {
             if(matrix_pol[j*sumnsam+y] != '-') {
             matrix_hap[y*length+k] = matrix_pol[j*sumnsam+y];
             }
             else {
             flag = 1; 
             for(yy=inits;yy<=y;yy++) matrix_hap[yy*length+k] = (char)0;
             break;
             }
             }
             if(flag == 0) k++;
             }*/
            /*calculate hapa*/
            statistics[0].hapa[z] = 0.;
            for(y=initsq1[pop1];y<initsq1[pop1]+nsam[pop1];y++) {
                for(x=initsq1[pop2];x<initsq1[pop2]+nsam[pop2];x++) {
                    if(memcmp(matrix_hap+y*length,matrix_hap+x*length,k) != 0) {
                        statistics[0].hapa[z] += 1;
                        statistics[0].hapT[z] += 1;
                    }
                }
            }
            statistics[0].hapa[z] /= (double)(nsam[pop1]*nsam[pop2]);
            
            /*calculate hapT*/
            for(y=initsq1[pop1];y<initsq1[pop1]+nsam[pop1]-1;y++) {
                for(x=y+1;x<initsq1[pop1]+nsam[pop1];x++) {
                    if(memcmp(matrix_hap+x*length,matrix_hap+y*length,k) != 0) {
                        statistics[0].hapT[z] += 1;
                    }
                    
                }
            }
            for(y=initsq1[pop2];y<initsq1[pop2]+nsam[pop2]-1;y++) {
                for(x=y+1;x<initsq1[pop2]+nsam[pop2];x++) {
                    if(memcmp(matrix_hap+x*length,matrix_hap+y*length,k) != 0) {
                        statistics[0].hapT[z] += 1;
                    }
                    
                }
            }
            statistics[0].hapT[z] /= (double)(nsam[pop1]+nsam[pop2])*(nsam[pop1]+nsam[pop2]-1.0)/2.0;
            z++;
        }
    }
    
    /*calculate statistics related to fsth1all and fsthTall*/
    /*
     if(npops==2) {
     nn = sumnsam;
     np = 0;
     }
     else {
     nn = sumnsam-nsam[npops-1];
     np = 1;
     }*//*
         if(npops > 2) {
         nn = sumnsam-nsam[npops-1];
         np = 1;
         k = 0;
         for(j=0;j<length;j++) {
         flag = 0;
         for(y=0;y<nn;y++) {
         if(matrix_pol[j*sumnsam+y] != '-') {
         matrix_hap[y*length+k] = matrix_pol[j*sumnsam+y];
         }
         else {
         flag = 1; 
         for(yy=0;yy<=y;yy++) matrix_hap[yy*length+k] = (char)0;
         break;
         }
         }
         if(flag == 0) k++;
         }
         for(pop1=0;pop1<npops-np;pop1++) {
         hapw1all[pop1] = 0.;
         inits = initsq1[pop1];
         max   = initsq1[pop1]+nsam[pop1];
         for(y=inits;y<max-1;y++) {
         for(z=y+1;z<max;z++) {
         if(memcmp(matrix_hap+y*length,matrix_hap+z*length,k) != 0) {
         hapw1all[pop1] += 1;
         }
         }
         }
         if(nsam[pop1]>1) hapw1all[pop1] /= (nsam[pop1]*(nsam[pop1]-1.0)/2.0);
         else hapw1all[pop1] = 0.;
         }
         
         for(pop1=0;pop1<npops-np;pop1++) {
         hapa1all[pop1] = 0.;
         hapT1all[pop1] = 0.;
         }
         
         for(pop1=0;pop1<npops-1;pop1++) {				
         for(pop2=pop1+1;pop2<npops-np;pop2++) {
         if(nsam[pop1] > 1 && nsam[pop2] > 1) {
         for(y=initsq1[pop1];y<initsq1[pop1]+nsam[pop1];y++) {
         for(x=initsq1[pop2];x<initsq1[pop2]+nsam[pop2];x++) {
         if(memcmp(matrix_hap+y*length,matrix_hap+x*length,k) != 0) {
         hapa1all[pop1] += 1.0/(double)(nsam[pop1]*nsam[pop2]);
         hapa1all[pop2] += 1.0/(double)(nsam[pop1]*nsam[pop2]);
         hapT1all[pop1] += 1.0;
         hapT1all[pop2] += 1.0;
         }
         }
         }
         }
         if(nsam[pop1] > 1) {
         for(y=initsq1[pop1];y<initsq1[pop1]+nsam[pop1]-1;y++) {
         for(x=y+1;x<initsq1[pop1]+nsam[pop1];x++) {
         if(memcmp(matrix_hap+y*length,matrix_hap+x*length,k) != 0) {
         hapT1all[pop1] += 1.0;
         }
         }
         }
         }
         if(nsam[pop2] > 1) {
         for(y=initsq1[pop2];y<initsq1[pop2]+nsam[pop2]-1;y++) {
         for(x=y+1;x<initsq1[pop2]+nsam[pop2];x++) {
         if(memcmp(matrix_hap+y*length,matrix_hap+x*length,k) != 0) {
         hapT1all[pop2] += 1.0;
         }
         }
         }
         }
         }
         }
         }
         */
    /*calculate fsth*/
    z = 0;
    for(pop1=0;pop1<npops-1;pop1++) {
        for(pop2=pop1+1;pop2<npops-0;pop2++) {
            if(statistics[0].hapa[z] && nsam[pop1] > 1 && nsam[pop2] > 1) 
                statistics[0].fsth[z] = 1.0 - ((statistics[0].hapw[pop1] + statistics[0].hapw[pop2])/2.0)/statistics[0].hapa[z];
            else statistics[0].fsth[z] = -10000;
            z++;
        }
    }
    /*calculate Gst (Nei 1973) for pair-pair comparisons*/
    z = 0;
    for(pop1=0;pop1<npops-1;pop1++) {
        for(pop2=pop1+1;pop2<npops-0;pop2++) {
            if(nsam[pop1] > 1 && nsam[pop2] > 1) {
                Js = Dm = Dkl = 0;
                for(y=0;y<sumnsam;y++) {
                    Js += (double)(sfreqh[pop1][y])/(double)nsam[pop1]*(double)(sfreqh[pop1][y])/(double)nsam[pop1];
                    Js += (double)(sfreqh[pop2][y])/(double)nsam[pop2]*(double)(sfreqh[pop2][y])/(double)nsam[pop2];
                    Dkl += (double)pow(((double)(sfreqh[pop1][y])/(double)nsam[pop1] - (double)(sfreqh[pop2][y])/(double)nsam[pop2]),2.0)/2.0;
                }
                Js /= 2.0;
                Hs = 1.0 - Js;
                Dm += Dkl/2.0;
                Ht = Hs + Dm;
                statistics[0].Gst[z] = Dm / Ht;
                /**/
                nm = 2.0/(1.0/(double)nsam[pop1] + 1.0/(double)nsam[pop2]);
                hs = ht = 0.;
                for(y=0;y<sumnsam;y++) {
                    hs += (double)(sfreqh[pop1][y])/(double)nsam[pop1]*(double)(sfreqh[pop1][y])/(double)nsam[pop1];
                    hs += (double)(sfreqh[pop2][y])/(double)nsam[pop2]*(double)(sfreqh[pop2][y])/(double)nsam[pop2];
                    xi  = (double)(sfreqh[pop1][y])/(double)nsam[pop1]+(double)(sfreqh[pop2][y])/(double)nsam[pop2];
                    xi /= 2.0;
                    ht += xi*xi;
                }
                hs /= 2.0;
                hs = 1.0 - hs;
                hs = nm*hs/(nm - 1.0);
                ht = 1.0 - ht;
                ht = ht + hs/(nm*2.0);
                statistics[0].Gst[z] = 1.0 - hs/ht;
                /**/
            }
            else statistics[0].Gst[z] = -10000;
            z++;
        }
    }
    
    /*calculate fsth1ALL*/
    if(npops > 2) {/*
                    z = 0;
                    spiw = spiwn = spia = 0.;
                    ncw = 0;
                    for(pop1=0;pop1<npops-1;pop1++) {
                    if(nsam[pop1] > 1) {
                    spiwn += statistics[0].hapw[pop1];
                    spiw  += statistics[0].hapw[pop1] * nsam[pop1] * (nsam[pop1]-1.0) / 2.0;
                    ncw += 1;
                    }
                    for(pop2=pop1+1;pop2<npops-1;pop2++) {
                    if(nsam[pop1] + nsam[pop2] > 1) {
                    spia += statistics[0].hapa[z] * nsam[pop1] * nsam[pop2];
                    }
                    z++;
                    }
                    }*/
        for(pop1=0;pop1<npops-1;pop1++) {
            /*if(nsam[pop1] > 1 && (spiw+spia) > 0.)
             statistics[0].fsth1all[pop1] = 1.0 - ((statistics[0].hapw[pop1])/((spiw+spia)/((sumnsam-nsam[npops-1])*(sumnsam-nsam[npops-1]-1.0)/2.0)));*/
            if(nsam[pop1] > 1 && nsam_allpop1[pop1] > 1 && hapa1all[pop1] > 0) {
                statistics[0].fsth1all[pop1] = 1.0 - ((statistics[0].hapw[pop1] + hapw1all[pop1])/2.0) / hapa1all[pop1];
            }
            else statistics[0].fsth1all[pop1] = -10000;
        }
        /*calculate fstALL*/
        /*
         if((spiw+spia) > 0.) 
         statistics[0].fsthALL = 1.0 - ((spiwn/(double)ncw)/((spiw+spia)/((sumnsam-nsam[npops-1])*(sumnsam-nsam[npops-1]-1.0)/2.0)));
         else statistics[0].fsthALL = -10000;
         */
        /*
         if((spia) > 0.) 
         statistics[0].fsthALL = 1.0 - ((spiwn/(double)ncw)/((spia)/((sumnsam-nsam[npops-1])*(sumnsam-nsam[npops-1]-1.0)/2.0)));
         else statistics[0].fsthALL = -10000;
         */
        
        /*
         for(pop1=0;pop1<npops-1;pop1++) {
         hapw = 0.;
         ncw = 0;
         for(pop2=0;pop2<npops-1;pop2++) {
         if(pop2 != pop1 && nsam[pop2] > 1) {
         hapw += hapw1all[pop2];
         ncw += 1;
         }
         }
         hapw /= (double)ncw;
         nca = ncw;
         if(hapa1all[pop1] && nsam[pop1] > 1 && ncw > 0) 
         statistics[0].fsth1all[pop1] = 1.0 - ((hapw)/2.0)/(hapa1all[pop1]/(double)nca);
         else statistics[0].fsth1all[pop1] = -10000;
         }		
         */
        /*calculate fsthALL*/
        /**/
        shapw = shapa = 0.;
        ncw = nca = 0;
        
        z = 0;
        for(pop1=0;pop1<npops-1;pop1++) {
            if(nsam[pop1] > 1) {
                shapw += statistics[0].hapw[pop1];
                ncw += 1;
            }
            for(pop2=pop1+1;pop2<npops-0;pop2++) {
                if(nsam[pop1] > 1 && nsam[pop2] > 1 && pop2 < npops-1) {
                    shapa += statistics[0].hapa[z];
                    nca += 1;
                }
                z++;
            }
        }
        if(shapa) statistics[0].fsthALL = 1.0 - (shapw/(double)ncw)/(shapa/(double)nca);
        else statistics[0].fsthALL = -10000;
        
        /**/
        /*number of pops with nsam>1 nhp*/
        nhp = 0;
        for(pop1=0;pop1<npops-1;pop1++) {
            if(nsam[pop1]>1) nhp++;
        }
        
        /*calculate GstALL (G'st = Dm/H't)*/
        Js = 0.;
        for(pop1=0;pop1<npops-1;pop1++) {
            if(nsam[pop1] > 1) {
                for(z=0;z<sumnsam;z++) {
                    Js += (double)(sfreqh[pop1][z])/(double)nsam[pop1]*(double)(sfreqh[pop1][z])/(double)nsam[pop1];
                }
            }
        }
        Js /= (double)nhp;
        Hs = 1.0 - Js;
        
        Dm = 0.;
        for(pop1=0;pop1<npops-1;pop1++) {
            if(nsam[pop1] > 1) {
                for(pop2=0;pop2<npops-0;pop2++) {
                    if(nsam[pop2] > 1 && pop1 != pop2) {
                        Dkl = 0.;
                        for(z=0;z<sumnsam;z++) {
                            Dkl += (double)pow(((double)(sfreqh[pop1][z])/(double)nsam[pop1] - (double)(sfreqh[pop2][z])/(double)nsam[pop2]),2.0)/2.0;
                        }
                        Dm += Dkl;
                    }
                }
            }
        }
        Dm /= (double)(nhp*(nhp-1));
        Ht = Hs + Dm;
        statistics[0].GstALL = Dm / Ht;
        
        /**/
        nm = 0.;
        for(pop1=0;pop1<npops-1;pop1++) {
            if(nsam[pop1] > 1) {
                nm += 1.0/(double)(nsam[pop1]);
            }
        }
        nm = (double)nhp/nm;		
        
        hs = 0.;
        for(pop1=0;pop1<npops-1;pop1++) {
            if(nsam[pop1] > 1) {
                for(z=0;z<sumnsam;z++) {
                    hs += (double)(sfreqh[pop1][z])/(double)nsam[pop1]*(double)(sfreqh[pop1][z])/(double)nsam[pop1];
                }
            }
        }
        hs /= (double)nhp;
        hs = 1.0 - hs;
        hs = nm*hs/(nm - 1.0);
        
        ht = 0.;
        for(y=0;y<sumnsam;y++) {
            xi = 0.;
            for(pop1=0;pop1<npops-1;pop1++) {
                if(nsam[pop1] > 1) {
                    xi += (double)(sfreqh[pop1][y])/(double)nsam[pop1];
                }
            }
            xi /= (double)nhp;
            ht += xi*xi;
        }
        ht = 1.0 - ht;
        ht = ht + hs/(nm*(double)nhp);
        
        statistics[0].GstALL = 1.0 - hs/ht;
        /**/
        
        
    }	
    
    for(x=0;x<npops;x++)
        free(sfreqh[x]);
    free(sfreqh);
    free(valh);
    free(initsq1);
    /**/
    free(hapa1all);
    free(hapT1all);
    free(hapw1all);
    /**/
    free(matrix_hap);
    free(nsam_allpop1);
    
    return 1;
}

