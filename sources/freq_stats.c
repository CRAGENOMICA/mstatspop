/*
 *  freq_stats.c
 *  mstatspop
 *
 *  Created by Sebastian Ramos-Onsins on 19/08/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "freq_stats.h"

#define DEBUG_EMANUELE 0
#define NITER 1e6

/**
 *
 * @param npops
 * @param nsam
 * @param matrix_pol
 * @param length
 * @param statistics
 * @param outgroup_presence
 * @param force_outgroup
 * @return
 */
int calc_freqstats(int npops, int *nsam, char *matrix_pol,long int length, struct stats *statistics,int outgroup_presence,int force_outgroup, int include_unknown/*,long int *miss_outg*/, int n_ccov, int H1frq)
{
	long int j;
	int pop1;
	long int freq[4],freqo[4],**sfreq,**snfreq,fr;
	int sumnsam,*initsq1,inits,max;
	char ancestral[1];
	int x,y;
	double an 		= 0.;
	double bn 		= 0.;
	double ts       = 0.;
	double thetaS 	= 0.;
	double thetaT 	= 0.;
	double thetaSo 	= 0.;
	double thetaTo 	= 0.;
	double thetaFL 	= 0.;
	double thetaFW 	= 0.;
	double thetaL 	= 0.;
	double thetaSA 	= 0.;
	double thetaTA 	= 0.;
	
	double anx      = 0.;
	double bnx      = 0.;
	double any      = 0.;
	
	/**/double *coef,n;/**/
	
	long int S,Si;
	double Sc,Sco;
	double *w1,*w2;
	long int *sfreqn;
	int maxnsam;

	/*missing values: weights and samples*/
	double **ominx;
	long int **eix;
	int **nx;
	int *no;
	
	struct covar_ij *Covi1i2;
	static long int count_cov = 0;	
	
	double *mean_freqsptr,mean_S;

	/*init*/
	for(pop1=0;pop1<npops-1;pop1++) {
		statistics[0].S[pop1]  = 0;
		statistics[0].thetaS[pop1]  = thetaS;
		statistics[0].thetaT[pop1]  = thetaT;
		statistics[0].So[pop1]  = 0;
		statistics[0].thetaSo[pop1]  = thetaSo;
		statistics[0].thetaTo[pop1]  = thetaTo;
		statistics[0].thetaFL[pop1] = thetaFL;
		statistics[0].thetaFW[pop1] = thetaFW;
		statistics[0].thetaL[pop1]  = thetaL;
		statistics[0].thetaSA[pop1] = thetaSA;
		statistics[0].thetaTA[pop1] = thetaTA;
		
		statistics[0].Dtaj[pop1] = -10000;
		statistics[0].Dfl[pop1] = -10000;
		statistics[0].Ffl[pop1] = -10000;
		statistics[0].Hnfw[pop1] = -10000;
		statistics[0].Ez[pop1] = -10000;
        statistics[0].Yach[pop1] = -10000;
        statistics[0].FH[pop1] = -10000;
        
		statistics[0].To_ii[pop1] = -10000.0;
		statistics[0].To_i0[pop1] = -10000.0;
		statistics[0].ToH0_00[pop1] = -10000.0;
		statistics[0].ToH0_ii[pop1] = -10000.0;
    }
	
    if(length==0) return 1;
    if(count_cov == 0) {
        if((Covi1i2 = (struct covar_ij *)calloc(1,sizeof(struct covar_ij))) == 0) {
            printf("Error calloc memory freqtest1");
            exit(1);
        }
    }

    ancestral[0] 	= 0;
    coef 			= (double *)calloc(24,sizeof(double));
    initsq1 		= (int *)calloc(npops,sizeof(int));
    sfreq   		= (long int **)calloc(npops,sizeof(long int *));
    snfreq   		= (long int **)calloc(npops,sizeof(long int *));
    
    sumnsam = 0;
    maxnsam = 0;

    for(x=0;x<npops;x++)
    {
        initsq1[x] = sumnsam;
        sumnsam += nsam[x];
        if(maxnsam < nsam[x]) maxnsam = nsam[x];
        sfreq[x] = (long int *)calloc(nsam[x]+1,sizeof(long int));
        snfreq[x] = (long int *)calloc(nsam[x]+1,sizeof(long int));
    }
    
    w1 		= (double *) calloc(maxnsam,sizeof(double));
    w2 		= (double *) calloc(maxnsam,sizeof(double));
    sfreqn 	= (long int *)calloc(maxnsam,sizeof(long int));

    ominx	= (double **) calloc(maxnsam,sizeof(double *));	
    eix     = (long int **)calloc(maxnsam,sizeof(long int *));
    nx      = (int **)calloc(maxnsam,sizeof(int *));
    no      = (int *)calloc(length,sizeof(int));

    for(x=0;x<maxnsam;x++) {
        ominx[x]   = (double *) calloc(length,sizeof(double));		
        eix[x]     = (long int *) calloc(length,sizeof(long int));
        nx[x]      = (int *) calloc(length,sizeof(int));
    }

    /*calculate the frequencies of each population for each position*/
    if(outgroup_presence+force_outgroup==1)
    {
        for(pop1=0;pop1<npops-1;pop1++) {
            an = 0;
            bn = 0;
            Si = 0;
            S  = 0;
            Sc = 0.;
            thetaS = 0.;
            thetaT = 0.;
            thetaSo = 0.;
            thetaTo = 0.;
            thetaFL = 0;
            thetaFW = 0.;
            thetaL = 0.;
            thetaSA = 0.;
            thetaTA = 0.;
            
            anx = 0.;
            bnx = 0.;
            any = 0.;

            statistics[0].S[pop1]  = S;
            statistics[0].thetaS[pop1]  = thetaS;
            statistics[0].thetaT[pop1]  = thetaT;
            statistics[0].So[pop1]  = Sc;
            statistics[0].thetaSo[pop1]  = thetaSo;
            statistics[0].thetaTo[pop1]  = thetaTo;
            statistics[0].thetaFL[pop1] = thetaFL;
            statistics[0].thetaFW[pop1] = thetaFW;
            statistics[0].thetaL[pop1]  = thetaL;
            statistics[0].thetaSA[pop1] = thetaSA;
            statistics[0].thetaTA[pop1] = thetaTA;
            
            statistics[0].Dtaj[pop1] = -10000;
            statistics[0].Dfl[pop1] = -10000;
            statistics[0].Ffl[pop1] = -10000;
            statistics[0].Hnfw[pop1] = -10000;
            statistics[0].Ez[pop1] = -10000;
            statistics[0].Yach[pop1] = -10000;
            statistics[0].FH[pop1] = -10000;
            
            for(x=1;x<nsam[pop1];x++) {
                statistics[0].freq[pop1][x] = 0;
                sfreqn[x] = 0;
            }
            if(nsam[pop1] > 1) {
                for(j=0;j<length;j++) {		
                    freqo[0]=freqo[1]=freqo[2]=freqo[3]=0;
                    /*outgroup*/
                    for(y=initsq1[npops-1];y<sumnsam;y++) {
                        if(matrix_pol[j*sumnsam+y] == '0') {freqo[1] += 1;freqo[0] += 1;}
                        if(matrix_pol[j*sumnsam+y] == '1') {freqo[2] += 1;freqo[0] += 1;}
                        if(matrix_pol[j*sumnsam+y] == '-') {freqo[3] += 1;}
                    }
                    if(freqo[0]) {
                        if(freqo[1] != freqo[0] && freqo[1] != 0) {
                            if(pop1==0) sfreq[npops-1][freqo[1]] += (int)1;
                            ancestral[0] = (char)0;/*if the outgroup is polymorphic, we do not consider for neutrality tests with outgroup*/
                        }
                        else {
                            if(freqo[1] == freqo[0]/*>= (double)freqo[0]/2.0*/) {
                                ancestral[0] = '0';
                            }
                            if(freqo[2] == freqo[0]/*> (double)freqo[0]/2.0*/) {
                                ancestral[0] = '1';
                            }	
                        }
                    }
                    else ancestral[0] = (char)0;
                    no[j] = ancestral[0]; /*if the outgroup is valid 1, if not 0*/
                    /*if(ancestral[0]==0) *miss_outg +=1; */
                    /*pop*/
                    inits = initsq1[pop1];
                    max   = initsq1[pop1]+nsam[pop1];
                    freq[0]=freq[1]=freq[2]=freq[3]=0;	
                    for(y=inits;y<max;y++) {eix[y-inits][j] = 0;}
                    for(y=inits;y<max;y++) {
                        if(matrix_pol[j*sumnsam+y] == '0') {freq[1] += 1;freq[0] += 1;nx[y-inits][j] = 1;}
                        if(matrix_pol[j*sumnsam+y] == '1') {freq[2] += 1;freq[0] += 1;nx[y-inits][j] = 1;}
                        if(matrix_pol[j*sumnsam+y] == '-') {freq[3] += 1;nx[y-inits][j] = 0;}
                    }
                    if(freq[0]) {
                        if(freq[2]>0 && freq[2]<freq[0]) {
                            Sc += 1.0;
                            snfreq[pop1][freq[2]] += (int)1; /*no mhits allowed!*/
                        }
                        if(ancestral[0] == '0') {
                            if(freq[2]>0 && freq[2]<freq[0]) {
                                S += 1;
                                sfreq[pop1][freq[2]] += (int)1; /*no mhits allowed!*/
                            }
                            /*else sfreq[pop1][freq[0]+freq[3]] += (int)1; *//*no mhits allowed!*/
                        }
                        if(ancestral[0] == '1') {
                            if(freq[1]>0 && freq[1]<freq[0]) {
                                S += 1;
                                sfreq[pop1][freq[1]] += (int)1;
                            }
                            /*else sfreq[pop1][freq[0]+freq[3]] += (int)1; *//*no mhits allowed!*/
                            fr = freq[1];
                            freq[1] = freq[2];
                            freq[2] = fr;
                        }
                        /*
                        if(ancestral[0] == (char)0) {
                            if(freq[1]>0 && freq[1]<freq[0]) {
                                sfreqn[freq[1]] += (int)1;
                                Sc += 1.0;
                            }
                        }
                        */
                    }
                    /*FOR MISSING VALUES!!!*//**/
                    if(freq[0]) {
                        if(freq[2]>0 && freq[2]<freq[0]) {
                            an = 0.; for(x=1;x<freq[0];x++) an += 1./(double)x;
                            any += an;
                        }
                    }
                    if(freq[2]>0 && freq[2]<freq[0]) {
                        /*an = 0.; for(x=1;x<freq[0];x++) an += 1./(double)x;*/
                        bn = 0.; for(x=1;x<freq[0];x++) bn += 1./((double)x*x);
                        anx += an;
                        bnx += bn;
                        Si += 1;
                        
                        thetaS += (double)1/(double)an;				
                        thetaT += (double)(freq[1]*freq[2])/(double)(freq[0]*(freq[0]-1.0)/2.0);				
                        
                        eix[freq[2]][j] = 1; 
                        
                        if(ancestral[0]) {
                            thetaSo += (double)1/(double)an;				
                            thetaTo += (double)(freq[1]*freq[2])/(double)(freq[0]*(freq[0]-1.0)/2.0);				
                            if(freq[2] == 1) thetaFL += 1.;
                            thetaFW += ((double)freq[2] * (double)freq[2])/(double)(freq[0]*(freq[0]-1.0)/2.0);
                            thetaL  += (double)freq[2]/(double)(freq[0]-1.);					
                            if(freq[2] > 1) {
                                thetaSA += 1./(double)(an-1.0);
                                thetaTA += (double)freq[2] * (double)freq[1]/(double)(freq[0]*(freq[0]-1.0)/2.0) * (double)freq[0]/(double)(freq[0]-2.0);
                            }
                        }
                    }
                    /**/
                }
                /*	
                for(x=1;x<nsam[pop1];x++)
                    an += 1/(double)x;
                thetaS = (double)S/(double)an;
                
                for(x=1;x<nsam[pop1];x++)
                    thetaT += (double)x * (nsam[pop1]-x) * (double)sfreq[pop1][x];
                thetaT *= (double)1/(double)(nsam[pop1]*(nsam[pop1]-1.0)/2.0);
                
                thetaFL = sfreq[pop1][1];
                
                for(x=1;x<nsam[pop1];x++)
                    thetaFW += ((double)x * (double)x) * sfreq[pop1][x];
                thetaFW *= (double)1/(double)(nsam[pop1]*(nsam[pop1]-1.0)/2.0);
                
                for(x=1;x<nsam[pop1];x++)
                    thetaL += (double)x * (double)sfreq[pop1][x];
                thetaL *= (double)1.0/(double)(nsam[pop1]-(double)1);
                
                for(x=2;x<nsam[pop1];x++)
                    thetaSA += sfreq[pop1][x];
                thetaSA *= (double)1/(double)(an-1.0);
                
                for(x=2;x<nsam[pop1];x++)
                    thetaTA += (double)x * (nsam[pop1]-x) * (double)sfreq[pop1][x];
                thetaTA *= (double)1/(double)(nsam[pop1]*(nsam[pop1]-1.0)/2.0);
                pie1 = thetaTA;
                thetaTA *= (double)nsam[pop1]/(double)(nsam[pop1]-2);
                */
                
                statistics[0].So[pop1]      = S;
                statistics[0].thetaSo[pop1] = thetaSo;
                statistics[0].thetaTo[pop1] = thetaTo;
                statistics[0].S[pop1]       = Sc;
                statistics[0].thetaS[pop1]  = thetaS;
                statistics[0].thetaT[pop1]  = thetaT;
                statistics[0].thetaFL[pop1] = thetaFL;
                statistics[0].thetaFW[pop1] = thetaFW;
                statistics[0].thetaL[pop1]  = thetaL;
                statistics[0].thetaSA[pop1] = thetaSA;
                statistics[0].thetaTA[pop1] = thetaTA;
                /*
                statistics[0].anx[pop1]  = anx/(double)Si;
                statistics[0].bnx[pop1]  = bnx/(double)Si;
                */
                if(statistics[0].freq) {
                    for(x=1;x<nsam[pop1];x++) 
                        statistics[0].freq[pop1][x] = sfreq[pop1][x];
                }
                /*MISSING VALUES*/
                /*n = ceil(statistics[0].length[pop1]/statistics[0].total_length- 1E-5*//*round error*//*)*//*nsam[pop1]*//*;*/
                n = nsam[pop1];
                init_coef(coef,n);
                /*Sc = Sc;*/
                /*Sc += S;*//*thetaS * coef[0];*/
                Sco = S;/*thetaSo * coef[0];*/
                /*
                n = ceil(statistics[0].length[pop1]/statistics[0].total_length- 1E-5);
                init_coef(coef,n);
                */
                if(include_unknown == 0) {
                        /*
                        statistics[0].Dtaj[pop1] = tajima_d(thetaT,(long int)Sc,coef);
                        *//*
                        statistics[0].Dfl[pop1]  = fl_d(thetaFL,(long int)Sco,coef );
                        statistics[0].Ffl[pop1]  = fl_f(thetaFL,(long int)Sco,thetaTo,coef );
                        statistics[0].Hnfw[pop1] = fay_wu_normalized2(n,thetaL,thetaSo,(double)Sco,coef,thetaTo);
                        statistics[0].Ez[pop1]   = E_zeng(n,thetaL,thetaSo,(long int)Sco,coef);
                        */
                                    
                        /*check analysis using Achaz approach*/
                        /**/
                        for(x=1;x<=floor(nsam[pop1]/2);x++) {
                            if(x == nsam[pop1]-x) w1[x] = (double)nsam[pop1]/2.0;
                            else w1[x] = (double)nsam[pop1]/1.0;
                            if(x == nsam[pop1]-x) w2[x] = (double)nsam[pop1]/((double)x*(double)(nsam[pop1] - x)*2.0);
                            else w2[x] = (double)nsam[pop1]/((double)x*(double)(nsam[pop1] - x)*1.0);
                            if(x == nsam[pop1]-x) sfreqn[x] = (int)((snfreq[pop1][x] + snfreq[pop1][nsam[pop1]-x])/2.0);
                            else sfreqn[x] = snfreq[pop1][x] + snfreq[pop1][nsam[pop1]-x];
                        }
                        statistics[0].Dtaj[pop1] = freqtestn_achaz(nsam[pop1],sfreqn,1,w1,w2);
                        /*
                        for(x=1;x<=floor(nsam[pop1]/2);x++) {
                            if(x == nsam[pop1]-x) w1[x] = (double)nsam[pop1]/2.0;
                            else w1[x] = (double)nsam[pop1]/1.0;
                            if(x == nsam[pop1]-x) w2[x] = (double)nsam[pop1]/((double)x*(double)(nsam[pop1] - x)*2.0);
                            else w2[x] = (double)nsam[pop1]/((double)x*(double)(nsam[pop1] - x)*1.0);
                            
                            if(x == nsam[pop1]-x) sfreqn[x] = (int)((sfreq[pop1][x] + sfreq[pop1][nsam[pop1]-x] + sfreqn[x] + sfreqn[nsam[pop1]-x])/2.0);
                            else sfreqn[x] = sfreq[pop1][x] + sfreq[pop1][nsam[pop1]-x] + sfreqn[x] + sfreqn[nsam[pop1]-x];
                        }
                        statistics[0].Dtaj[pop1] = freqtestn_achaz(nsam[pop1],sfreqn,1,w1,w2);
                        */
                        /*
                        for(x=1;x<nsam[pop1];x++) {
                            w1[x] = (double)(nsam[pop1] - x);
                            w2[x] = 1.0/(double)x;
                        }
                        statistics[0].Dtaj[pop1] = freqtesto_achaz(nsam[pop1],sfreq[pop1],1,w1,w2);
                        */
                        for(x=1;x<nsam[pop1];x++) {
                            w1[x] = 1.0/(double)x;
                            if(x==1) w2[x] = 1.0;
                            else w2[x] = 0.0;
                        }
                        statistics[0].Dfl[pop1] = freqtesto_achaz(nsam[pop1],sfreq[pop1],1,w1,w2);
                        /*fu li F*/
                        for(x=1;x<nsam[pop1];x++) {
                            w1[x] = (double)(nsam[pop1] - x);
                            if(x==1) w2[x] = 1.0;
                            else w2[x] = 0.0;
                        }
                        statistics[0].Ffl[pop1] = freqtesto_achaz(nsam[pop1],sfreq[pop1],1,w1,w2);
                        /*fay wu H*/
                        for(x=1;x<nsam[pop1];x++) {
                            w1[x] = (double)(nsam[pop1] - x);
                            w2[x] = (double)x;
                        }
                        statistics[0].Hnfw[pop1] = freqtesto_achaz(nsam[pop1],sfreq[pop1],1,w1,w2);
                        /*zeng Z*/
                        for(x=1;x<nsam[pop1];x++) {
                            w1[x] = (double)1;
                            w2[x] = 1.0/(double)x;
                        }
                        statistics[0].Ez[pop1] = freqtesto_achaz(nsam[pop1],sfreq[pop1],1,w1,w2);
                        /*achaz Y*/
                        for(x=1;x<nsam[pop1];x++) {
                            if(x==1) w1[x] = 0.0;
                            else w1[x] = (double)(nsam[pop1] - x);
                            if(x==1) w2[x] = 0.0;
                            else w2[x] = 1.0/(double)x;
                        }
                        statistics[0].Yach[pop1] = freqtesto_achaz(nsam[pop1],sfreq[pop1],0,w1,w2);
                        /*statistics[0].Yach[pop1] = compute_Omega_Y(nsam[pop1],sfreq[pop1],1);*/
                        /*faywu H watterson*/
                        for(x=1;x<nsam[pop1];x++) {
                            w1[x] = 1.0/(double)x;
                            w2[x] = (double)x;
                       }
                       statistics[0].FH[pop1] = freqtesto_achaz(nsam[pop1],sfreq[pop1],1,w1,w2);
                }
                else {					
                    if(nsam[pop1] <= SAMPLE_LARGE){
                        if(ominx_tajD(ominx,nx,no,any,nsam[pop1],length,Sco))
                            statistics[0].Dtaj[pop1] = freqtest_outg_missing(ominx,eix,nx,no,thetaSo,nsam[pop1],length,n_ccov,&Covi1i2,&count_cov,1,&bn,&ts);
                        if(ominx_Dfl(ominx,nx,no,any,nsam[pop1],length,Sco))
                            statistics[0].Dfl[pop1]  = freqtest_outg_missing(ominx,eix,nx,no,thetaSo,nsam[pop1],length,n_ccov,&Covi1i2,&count_cov,1,&bn,&ts);
                        if(ominx_Ffl(ominx,nx,no,nsam[pop1],length,Sco))
                            statistics[0].Ffl[pop1]  = freqtest_outg_missing(ominx,eix,nx,no,thetaSo,nsam[pop1],length,n_ccov,&Covi1i2,&count_cov,1,&bn,&ts);
                        if(ominx_Hnfw(ominx,nx,no,nsam[pop1],length,Sco))
                            statistics[0].Hnfw[pop1] = freqtest_outg_missing(ominx,eix,nx,no,thetaSo,nsam[pop1],length,n_ccov,&Covi1i2,&count_cov,1,&bn,&ts);
                        if(ominx_Ezeng(ominx,nx,no,any,nsam[pop1],length,Sco))
                            statistics[0].Ez[pop1]   = freqtest_outg_missing(ominx,eix,nx,no,thetaSo,nsam[pop1],length,n_ccov,&Covi1i2,&count_cov,1,&bn,&ts);
                        if(ominx_Yachaz(ominx,nx,no,any,nsam[pop1],length,Sco))
                            statistics[0].Yach[pop1] = freqtest_outg_missing(ominx,eix,nx,no,thetaSA,nsam[pop1],length,n_ccov,&Covi1i2,&count_cov,0,&bn,&ts);
                        if(ominx_FH(ominx,nx,no,any,nsam[pop1],length,Sco))
                            statistics[0].FH[pop1] = freqtest_outg_missing(ominx,eix,nx,no,thetaSo,nsam[pop1],length,n_ccov,&Covi1i2,&count_cov,1,&bn,&ts);
                        if(H1frq) {
                            mean_freqsptr = (double *)calloc(nsam[pop1],sizeof(double));
                            for(y=0;y<nsam[pop1]-1;y++)  mean_freqsptr[y] = (double)statistics[0].H1freq[pop1][y+1];
                            mean_S = 0.0; for(y=0;y<nsam[pop1]-1;y++) mean_S += (double)mean_freqsptr[y];
                            for(y=0;y<nsam[pop1]-1;y++) mean_freqsptr[y] = mean_freqsptr[y]/mean_S;

                            statistics[0].To_ii[pop1] = -10000.0;
                            statistics[0].To_i0[pop1] = -10000.0;
                            statistics[0].ToH0_00[pop1] = -10000.0;
                            statistics[0].ToH0_ii[pop1] = -10000.0;
                            if(ominx_Optimal_unfolded(ominx,nx,no,any,nsam[pop1],length,Sco,mean_freqsptr))
                                statistics[0].To_00[pop1] = freqtest_outg_missing(ominx,eix,nx,no,thetaSA,nsam[pop1],length,n_ccov,&Covi1i2,&count_cov,1,&bn,&ts);
                            free(mean_freqsptr);
                        }
                    }
                }
            }
        }
    }
    if(outgroup_presence+force_outgroup==0) {
        for(pop1=0;pop1<npops-1;pop1++) {
            an = 0;
            bn = 0;
            S = 0;
            Sc = 0.;
            thetaS = 0;
            thetaT = 0.;
            thetaFL = 0.;
            thetaSA = 0.;
            thetaTA = 0.;
            
            anx = 0.;
            bnx = 0.;
            any = 0.;

            statistics[0].S[pop1]  = S;
            statistics[0].thetaS[pop1]  = thetaS;
            statistics[0].thetaT[pop1]  = thetaT;
            statistics[0].thetaFL[pop1] = thetaFL;
            statistics[0].thetaSA[pop1] = thetaSA;
            statistics[0].thetaTA[pop1] = thetaTA;
            
            statistics[0].Dtaj[pop1] = -10000;
            statistics[0].Dfl[pop1] = -10000;
            statistics[0].Ffl[pop1] = -10000;
            statistics[0].Yach[pop1] = -10000;

            if(nsam[pop1] > 1) {
                for(j=0;j<length;j++) {	
                    inits = initsq1[pop1];
                    max   = initsq1[pop1]+nsam[pop1];
                    freq[0]=freq[1]=freq[2]=freq[3]=0;
                    for(y=inits;y<max;y++) {eix[y-inits][j] = 0;}
                    for(y=inits;y<max;y++) {
                        if(matrix_pol[j*sumnsam+y] == '0') {freq[1] += 1;freq[0] += 1;nx[y-inits][j] = 1;}
                        if(matrix_pol[j*sumnsam+y] == '1') {freq[2] += 1;freq[0] += 1;nx[y-inits][j] = 1;}
                        if(matrix_pol[j*sumnsam+y] == '-') {freq[3] += 1;nx[y-inits][j] = 0;}
                    }
                    if(freq[0] > 1) {
                        if(freq[2] && freq[1]) {
                            S += 1;
                            sfreq[pop1][freq[2]] += (int)1; /*no mhits allowed!*/
                            //if(matrix_pol[j*sumnsam+max]== '-')
                            //    printf("Error!");
                        }
                        /*else sfreq[pop1][freq[0]+freq[3]] += (int)1; *//*no mhits allowed!*/
                    }
                    /*FOR MISSING VALUES!!!*//**/
                    if(freq[0] > 1) {
                        if(freq[2] && freq[1]) {
                            an = 0.; for(x=1;x<freq[0];x++) an += 1./(double)x;
                            any += an;
                        }
                    }
                    if(freq[2]>0 && freq[1]>0) {
                        /*an = 0.; for(x=1;x<freq[0];x++) an += 1./(double)x;*/
                        bn = 0.; for(x=1;x<freq[0];x++) bn += 1./((double)x*x);
                        anx += an;
                        bnx += bn;
                        thetaS += (double)1/(double)an;				
                        thetaT += (double)(freq[1]*freq[2])/(double)(freq[0]*(freq[0]-1.0)/2.0);
                        if((freq[2] == 1 || freq[1] == 1) && freq[0]==2) thetaFL += 1.0 * (freq[0]-1.)/(double)freq[0] * 2.0;
                        if((freq[2] == 1 || freq[1] == 1) && freq[0]!=2) thetaFL += 1.0 * (freq[0]-1.)/(double)freq[0];
                        if(freq[2] > 1 && freq[1] > 1) {
                            thetaSA += 1.0/(double)(an - (double)freq[0]/(double)(freq[0]-1.0));
                            thetaTA += (double)freq[2] * (double)freq[1]/(double)(freq[0]*(freq[0]-1.0)/2.0) * (double)(freq[0]-1.0)/(double)(freq[0]-3.0);
                        }
                        eix[freq[2]][j] = 1;
                    }
                    /**/
                }
                /*
                for(x=1;x<nsam[pop1];x++)
                    an += 1/(double)x;
                thetaS = (double)S/((double)an);
                
                for(x=1;x<nsam[pop1];x++)
                    thetaT += (double)x * (nsam[pop1]-x) * (double)sfreq[pop1][x];
                thetaT *= (double)1/(double)(nsam[pop1]*(nsam[pop1]-1.0)/2.0);

                thetaFL = (sfreq[pop1][1] + sfreq[pop1][nsam[pop1]-1])*(nsam[pop1]-1.0)/(nsam[pop1]);
                
                for(x=2;x<nsam[pop1]-1;x++)
                    thetaSA += sfreq[pop1][x];
                thetaSA *= (double)1/(double)(an - (double)nsam[pop1]/(double)(nsam[pop1]-1));
                
                for(x=2;x<nsam[pop1]-1;x++)
                    thetaTA += (double)x * (nsam[pop1]-x) * (double)sfreq[pop1][x];
                thetaTA *= (double)1/(double)(nsam[pop1]*(nsam[pop1]-1.0)/2.0);
                pie1 = thetaTA;
                thetaTA *= (double)(nsam[pop1] - 1.0)/(double)(nsam[pop1] - 3.0);
                */
                
                statistics[0].S[pop1]  = S;
                statistics[0].thetaS[pop1]  = thetaS;
                statistics[0].thetaT[pop1]  = thetaT;
                statistics[0].thetaFL[pop1] = thetaFL;
                statistics[0].thetaSA[pop1] = thetaSA;
                statistics[0].thetaTA[pop1] = thetaTA;
                /*
                statistics[0].anx[pop1]  = anx/(double)S;
                statistics[0].bnx[pop1]  = bnx/(double)S;
                */
                if(statistics[0].freq) {
                    for(x=1;x<nsam[pop1];x++) {
                        if(x != nsam[pop1]-x) statistics[0].freq[pop1][x] = sfreq[pop1][x] + sfreq[pop1][nsam[pop1]-x];
                        else statistics[0].freq[pop1][x] = sfreq[pop1][x];
                    }
                }
                /**/
                /*n = ceil(statistics[0].length[pop1]/statistics[0].total_length- 1E-5*//*round error*//*)*//*nsam[pop1]*//*;*/
                n = nsam[pop1];
                init_coef(coef,n);
                Sc = S;/*thetaS * coef[0];*/
                /*
                n = ceil(statistics[0].length[pop1]/statistics[0].total_length- 1E-5);
                init_coef(coef,n);
                */
                if(include_unknown == 0) {
                        /*
                        statistics[0].Dtaj[pop1] = tajima_d(thetaT,(long int)Sc,coef);
                        *//*
                        statistics[0].Dfl[pop1]  = fl_d2(n,thetaFL,(long int)Sc,coef);
                        statistics[0].Ffl[pop1]  = fl_f2(n,thetaFL,(long int)Sc,thetaT,coef);
                        */
                        /*check analysis using Achaz approach*/
                        /**/
                        for(x=1;x<=floor(nsam[pop1]/2);x++) {
                            if(x == nsam[pop1]-x) w1[x] = (double)nsam[pop1]/2.0;
                            else w1[x] = (double)nsam[pop1]/1.0;
                            if(x == nsam[pop1]-x) w2[x] = (double)nsam[pop1]/((double)x*(double)(nsam[pop1] - x)*2.0);
                            else w2[x] = (double)nsam[pop1]/((double)x*(double)(nsam[pop1] - x)*1.0);
                            if(x == nsam[pop1]-x) sfreqn[x] = (int)((sfreq[pop1][x] + sfreq[pop1][nsam[pop1]-x])/2.0);
                            else sfreqn[x] = sfreq[pop1][x] + sfreq[pop1][nsam[pop1]-x];
                        }
                        statistics[0].Dtaj[pop1] = freqtestn_achaz(nsam[pop1],sfreqn,1,w1,w2);
                        /**/
                        /*
                        for(x=1;x<=floor(nsam[pop1]/2);x++) {
                            if(x == nsam[pop1]-x) w1[x-1] = (double)nsam[pop1]/2.0;
                            else w1[x-1] = (double)nsam[pop1]/1.0;
                            if(x == nsam[pop1]-x) w2[x-1] = (double)nsam[pop1]/((double)x*(double)(nsam[pop1] - x)*2.0);
                            else w2[x-1] = (double)nsam[pop1]/((double)x*(double)(nsam[pop1] - x)*1.0);
                            if(x == nsam[pop1]-x) sfreqn[x-1] = (sfreq[pop1][x] + sfreq[pop1][nsam[pop1]-x])/2.0;
                            else sfreqn[x-1] = sfreq[pop1][x] + sfreq[pop1][nsam[pop1]-x];
                        }
                        TD = compute_T_star(nsam[pop1],sfreqn,w1,w2,1);
                        */
                        /**/
                        for(x=1;x<=floor(nsam[pop1]/2);x++) {
                            if(x == nsam[pop1]-x) w1[x] = (double)nsam[pop1]/((double)x*(double)(nsam[pop1] - x)*2.0);
                            else w1[x] = (double)nsam[pop1]/((double)x*(double)(nsam[pop1] - x)*1.0);
                            if(x==1) w2[x] = nsam[pop1];
                            else w2[x] = 0.0;
                            if(x == nsam[pop1]-x) sfreqn[x] = (int)((sfreq[pop1][x] + sfreq[pop1][nsam[pop1]-x])/2.0);
                            else sfreqn[x] = sfreq[pop1][x] + sfreq[pop1][nsam[pop1]-x];
                        }
                        statistics[0].Dfl[pop1] = freqtestn_achaz(nsam[pop1],sfreqn,1,w1,w2);
                        for(x=1;x<=floor(nsam[pop1]/2);x++) {
                            if(x == nsam[pop1]-x) w1[x] = (double)nsam[pop1]/2.0;
                            else w1[x] = (double)nsam[pop1]/1.0;
                            if(x==1) w2[x] = nsam[pop1];
                            else w2[x] = 0.0;
                            if(x == nsam[pop1]-x) sfreqn[x] = (int)((sfreq[pop1][x] + sfreq[pop1][nsam[pop1]-x])/2.0);
                            else sfreqn[x] = sfreq[pop1][x] + sfreq[pop1][nsam[pop1]-x];
                        }
                        statistics[0].Ffl[pop1] = freqtestn_achaz(nsam[pop1],sfreqn,1,w1,w2);
                        /**/
                        for(x=1;x<=floor(nsam[pop1]/2);x++) {
                            if(x == 1) w1[x] = 0.0;
                            else {
                                if(x == nsam[pop1]-x) w1[x] = (double)nsam[pop1]/2.0;
                                else w1[x] = (double)nsam[pop1]/1.0;
                            }
                            if(x == 1) w2[x] = 0.0;
                            else {
                                if(x == nsam[pop1]-x) w2[x] = (double)nsam[pop1]/((double)x*(double)(nsam[pop1] - x)*2.0);
                                else w2[x] = (double)nsam[pop1]/((double)x*(double)(nsam[pop1] - x)*1.0);
                                    }
                            if(x == nsam[pop1]-x) sfreqn[x] = (int)((sfreq[pop1][x] + sfreq[pop1][nsam[pop1]-x])/2.0);
                            else sfreqn[x] = sfreq[pop1][x] + sfreq[pop1][nsam[pop1]-x];
                        }
                        statistics[0].Yach[pop1] = freqtestn_achaz(nsam[pop1],sfreqn,0,w1,w2);
                 }
                else {			
                    if(nsam[pop1] <= SAMPLE_LARGE) {
                        if(ominx_tajDnooutg(ominx,nx,any,nsam[pop1],length,Sc))
                            statistics[0].Dtaj[pop1] = freqtest_noutg_missing(ominx,eix,nx,thetaS,nsam[pop1],length,n_ccov,&Covi1i2,&count_cov,1,&bn,&ts);
                        if(ominx_Dfl2(ominx,nx,any,nsam[pop1],length,Sc))
                            statistics[0].Dfl[pop1]  = freqtest_noutg_missing(ominx,eix,nx,thetaS,nsam[pop1],length,n_ccov,&Covi1i2,&count_cov,1,&bn,&ts);
                        if(ominx_Ffl2(ominx,nx,nsam[pop1],length,Sc))
                            statistics[0].Ffl[pop1]  = freqtest_noutg_missing(ominx,eix,nx,thetaS,nsam[pop1],length,n_ccov,&Covi1i2,&count_cov,1,&bn,&ts);
                        if(ominx_Yachaz2(ominx,nx,any,nsam[pop1],length,Sc))
                            statistics[0].Yach[pop1] = freqtest_noutg_missing(ominx,eix,nx,thetaSA,nsam[pop1],length,n_ccov,&Covi1i2,&count_cov,0,&bn,&ts);
                        if(H1frq) {
                            mean_freqsptr = (double *)calloc(nsam[pop1],sizeof(double));
                            for(y=0;y<nsam[pop1]/2;y++)  mean_freqsptr[y] = (double)statistics[0].H1freq[pop1][y+1];
                            mean_S = 0.0; for(y=0;y<nsam[pop1]-1;y++) mean_S += (double)mean_freqsptr[y];
                            for(y=0;y<nsam[pop1]/2;y++) mean_freqsptr[y] = mean_freqsptr[y]/mean_S;
                            
                            statistics[0].To_ii[pop1] = -10000.0;
                            statistics[0].To_i0[pop1] = -10000.0;
                            statistics[0].ToH0_00[pop1] = -10000.0;
                            statistics[0].ToH0_ii[pop1] = -10000.0;
                            if(ominx_Optimal_folded(ominx,nx,any,nsam[pop1],length,Sc,mean_freqsptr))
                                statistics[0].To_00[pop1] = freqtest_noutg_missing(ominx,eix,nx,thetaSA,nsam[pop1],length,n_ccov,&Covi1i2,&count_cov,1,&bn,&ts);
                            free(mean_freqsptr);
                        }
                    }
                }
            }
            else {
                statistics[0].thetaS[npops-1]  = -10000;
                statistics[0].thetaT[npops-1]  = -10000;
                statistics[0].thetaFL[npops-1] = -10000;
                statistics[0].thetaSA[npops-1] = -10000;
                statistics[0].thetaTA[npops-1] = -10000;
                
                statistics[0].Dtaj[npops-1] = -10000;
                statistics[0].Dfl[npops-1]  = -10000;
                statistics[0].Ffl[npops-1]  = -10000;
                statistics[0].Yach[npops-1] = -10000;
            }
        }
    }
    /*free pointers*/
    free(Covi1i2);
    for(x=0;x<maxnsam;x++) {
        free(ominx[x]);
        free(eix[x]);
        free(nx[x]);
    }
    free(ominx);
    free(eix);
    free(nx);
    free(no);
    
    for(x=0;x<npops;x++) {
        free(sfreq[x]);
        free(snfreq[x]);
    }
    free(sfreq);
    free(snfreq);
    free(initsq1);
    free(w1);
    free(w2);
    free(sfreqn);
    free(coef);

	return 1;
}

void init_coef(double *p,int sample_size)
{
	/*Tajima and Fu coefficients... and Achaz coefficients. Total = 23 Most from Simonsen 1995*/
	int x,n;
	double an,bn,cn;
	
	an = bn = 0.0;
	n = sample_size;
	
	for(x=1;x<n;x++) {
		an += 1.0/(double)x;
		bn += 1.0/((double)x*(double)x);
	}
	
	p[0] = an;
	p[1] = bn;
	
	/* vt */
	p[2] = (2.0*((double)n*(double)n + (double)n + 3.0)/(9.0*(double)n*((double)n-1.0))
	        -((double)n+2.0)/(an * (double)n) + bn/(an * an)) / (an*an + bn);	
	/* ut */
	p[3] = (( ((double)n+1.0)/(3.0*((double)n-1.0)) - 1.0/an)/an) - p[2];
	/* vd* */
	p[4] = (bn/(an*an) - 2.0/(double)n *(1.0 + 1.0/an - an + an/(double)n) - 1.0/((double)n*(double)n))
	/ (an*an + bn);
	/* ud* */
	p[5] = ((((double)n-1.0)/(double)n - 1.0/an) / an) - p[4];
	/* vf* */
	p[6] = (((2.0*(double)n*(double)n*(double)n + 110.0 * (double)n*(double)n - 255.0 * (double)n + 153.0)
			 / (9.0 * (double)n * (double)n * ((double)n-1.0)) + (2.0*((double)n-1.0) *an)/ ((double)n*(double)n) 
			 - (8.0 * bn)/(double)n)) / (an*an + bn);
	/* uf* */
	p[7] = (((4.0*(double)n*(double)n + 19.0*(double)n + 3.0 - 12.0 * ((double)n +1.0) * (an + 1.0/(double)n))
			 / (3.0*(double)n *((double)n-1.0))) / an) - p[6];
	/* cn */
	cn = 2.0 * ((double)n*an-2.0*((double)n-1.0)) / (((double)n-1.0)*((double)n-2.0));
	/* vd */
	p[8] = 1.0 + (an*an/(bn+an*an)) * (cn - (((double)n+1.0)/((double)n-1.0)));
	/* ud* */
	p[9] = an -1.0 - p[8];
	/* vf */
	p[10] = (cn + (2.0*((double)n*(double)n+(double)n+3.0))/(9.0*(double)n*((double)n-1.0)) - 2.0
			 /((double)n-1.0)) / (an*an + bn);
	/* uf */
	p[11] = (1.0 + ((double)n+1.0)/(3.0*((double)n-1.0)) - 4*((double)n+1.0)
			 /(((double)n-1.0)*((double)n-1.0)) * (an+(1.0/(double)n) - 2.0*(double)n/((double)n+1.0)))
	/ an  -  p[10];
	/* Achaz f */
	p[12] = ((double)n - 2.0)/((double)n * (an - 1.0));
	/* Achaz f* */
	p[13] = ((double)n - 3.0)/(an * ((double)n - 1.0) - (double)n);
	/* alfa_n */
	p[14] = p[12] * p[12] * (an - 1.0) + 
	p[12] * (an * (4.0*((double)n+1.0))/(((double)n-1.0)*((double)n-1.0)) - (2.0 *((double)n + 1.0)*((double)n + 2.0))/((double)n * ((double)n-1.0))) -
	an * (8.0 * ((double)n+1.0)/((double)n*((double)n-1.0)*((double)n-1.0))) +
	((double)n*(double)n*(double)n + (double)n*(double)n + 60.0*(double)n + 12.0)/(3.0*(double)n*(double)n*((double)n-1.0));
	
	/* beta_n */
	p[15] = p[12] * p[12] * (bn + an * (4.0/(((double)n - 1.0) * ((double)n - 2.0))) - 4.0/((double)n-2.0)) +
	p[12] * (-an * (4.0*((double)n+2.0))/((double)n*((double)n-1.0)*((double)n-2.0)) - ((double)n*(double)n*(double)n - 3.0*(double)n*(double)n - 16.0*(double)n + 20.)/((double)n*((double)n-1.0)*((double)n-2.0))) +
	an * (8.0/((double)n*((double)n-1.0)*((double)n-2.0))) +
	2.0 * ((double)n*(double)n*(double)n*(double)n - (double)n*(double)n*(double)n - 17.0*(double)n*(double)n - 42.0*(double)n + 72.0)/(9.0*(double)n*(double)n*((double)n-1.0)*((double)n-2.0));
	
	/* alfa*_n */
	p[16] = p[13] * p[13] * (an - (double)n/((double)n-1.0)) + 
	p[13] * (an * 4.0 * ((double)n+1.0)/(((double)n-1.0)*((double)n-1.0)) - 2.0 * ((double)n+3.0)/((double)n-1.0)) -
	an * (8.0*((double)n+1.0))/((double)n*((double)n-1.0)*((double)n-1.0)) +
	((double)n*(double)n + (double)n + 60.0)/(3.0*(double)n*((double)n-1.0));
	
	/* beta*_n */
	p[17] = p[13] * p[13] * (bn - (2.0*(double)n-1.0)/(((double)n-1.0)*((double)n-1.0))) +
	p[13] * (bn * 8.0/((double)n-1.0) - an * 4.0/((double)n*((double)n-1.0)) - ((double)n*(double)n*(double)n + 12.0*(double)n*(double)n - 35.0*(double)n +18.0)/((double)n*((double)n-1.0)*((double)n-1.0))) -
	bn * 16.0/((double)n*((double)n-1.0)) +
	an * 8.0/((double)n*(double)n*((double)n-1.0)) +
	2.0 * ((double)n*(double)n*(double)n*(double)n + 110.0*(double)n*(double)n - 255.0*(double)n + 126.0)/(9.0*(double)n*(double)n*((double)n-1.0)*((double)n-1.0));
	
	/*b1*/
	p[18] =(double)(n+1.0)/(3.0*((double)n-1.0));
	/*b2*/
	p[19] =2.0*((double)n*(double)n + (double)n + 3)/(9*(double)n*((double)n-1));
	/*c1*/
	p[20] = p[18] - 1.0/p[0];
	/*c2*/
	p[21] =p[19] - ((double)n+2.0)/(p[0]*(double)n) + (p[1]/(p[0]*p[0]));
	/*e1*/
	p[22] =p[20]/p[0];
	/*e2*/
	p[23] =p[21]/(p[0]*p[0] + p[1]);
}

/*Tajima's D*/
double tajima_d(double k_, long int S_, double *coef_taj)
{
	double an/*,ut,vt*/;
	double S_D = -10000;
	
	if(S_ == 0 || *(coef_taj+0) < 1.51) return(-10000); 
	
	an = *(coef_taj+0);
	/*
	ut = *(coef_taj+3);
	vt = *(coef_taj+2);
	
	S_D = (k_ - ((double)S_/an)) / (sqrt((ut*(double)S_) + (vt*(double)S_*((double)S_))));
	*/
	S_D = (k_ - ((double)S_/an)) / (sqrt((coef_taj[22]*(double)S_) + ((coef_taj[23])*(double)S_*((double)S_-1.0))));
	
	if (fabs(S_D) < 1.0E-15)
		S_D = 0.0;
	
	return S_D;
}
/*Fu and Li's D*/
/**
 *
 * @param fr1
 * @param S
 * @param coef
 * @return
 */
double fl_d( long int fr1,long int S, double *coef ) /* amb outgroup */
{
	double an;
	double ud,vd;
	long int re;
	double D 	= -10000;

	if( S == 0 || *(coef+0) < 1.5) return(-10000); /**< todo: Como se debe leer? ambiguo */
	
	re = fr1;	
	an = *(coef+0);
	vd = *(coef+8);
	ud = *(coef+9);
	D  = ((double)S - an*(double)re) / 
		sqrt(ud*(double)S + vd*(double)S*(double)S);
	
	return D;
}

/*Fu and Li's F*/
double fl_f( long int fr1, long int S, double pi, double *coef) /* amb outgroup */
{
	double 		uf,vf;
	long int 	re;
	double 		F;

	if(S == 0 || *(coef+0) < 1.5) return(-10000);
	
	re = fr1;		
	vf = *(coef+10);
	uf = *(coef+11);
	
	F  = (pi - (double)re) / sqrt(uf*(double)S + vf*(double)S*(double)S);
	
	return F;
}

double fay_wu_normalized2(int n,double thetaL,double thetaw,double S,double *coef,double pi) /* eq 11 and 12 from Fay and Wu H NORMALIZED (Zeng et al. Genetics 2006 174: 1431-9)*/
{
    double H,varpiTL,an,bn;
    
    if(pi == (double)0.0 || n < 2) return(-10000);
    an = coef[0];
	bn = coef[1];
	
	varpiTL = thetaw * ((double)(n-(double)2))/((double)6*((double)(n-(double)1))) + 
	 (((S*(S-(double)1)/(an*an+bn)) / ((double)9*((double)n*((double)n-(double)1)*((double)n-(double)1)))) *
	 ((double)18*(double)n*(double)n*((double)3*(double)n+(double)2)*(bn+(double)1.0/((double)n*(double)n)) - 
	 ((double)88*(double)n*(double)n*(double)n + (double)9*(double)n*(double)n - 13*(double)n + (double)6)));
	
	H = (pi - thetaL)/(double)sqrt(varpiTL);
	
    return H;
}

double E_zeng(int n,double thetaL,double thetaw,double S,double *coef) /* (eq 13 and 14 from Zeng et al. Genetics 2006 174: 1431-9)*/
{
    double E,varLW,an,bn;
    
    if(thetaw == (double)0.0 || n < 2) return(-10000);
    an = coef[0];
	bn = coef[1];
	varLW = thetaw * ((double)n/((double)2*(double)(n-1)) - (double)1/an) +
	S*(S-(double)1)/(an*an+bn) * 
	(bn/(an*an) + (double)2*bn*((double)n/(double)(n-1))*((double)n/(double)(n-1)) - 
	 (double)2*((double)n*bn-(double)n+(double)1)/((double)(n-1)*an) - 
	 ((double)3*(double)n+(double)1)/((double)(n-1)));
	
	E = (thetaL - thetaw)/(double)sqrt(varLW);
	
    return E;
}
/* Fu and Li's D* */
double fl_d2(int sample_size,long int fr1w,long int S, double *coef) /* NO outgroup */
{
	double an;
	int n;
	double ud2,vd2;
	long int rs;
	double D2 = -10000;
	
	if(S == 0 || *(coef+0) < 1.51) return(-10000);
	
	rs = fr1w;
	
	n = sample_size;
	an = *(coef+0);
	
	vd2 = *(coef+4);
	ud2 = *(coef+5);
	D2  = ((double)S/an - (double)rs*(((double)n-(double)1)/(double)n)) /
	sqrt(ud2*(double)S + vd2*(double)S*(double)S);
	
	return D2;
}
/* Fu and Li's F* */
double fl_f2(int sample_size,long int fr1w, long int S, double pi, double *coef) /* NO outgroup */
{
	int n;
	double uf2,vf2;
	long int rs;
	double F2;
	
	if(S == 0 || *(coef+0) < 1.51) return(-10000);
	
	rs = fr1w;
	
	n   = sample_size;
	vf2 = *(coef+6);
	uf2 = *(coef+7);
	
	F2  = (pi - ((((double)n-1.0)/(double)n)*(double)rs)) / 
	sqrt(uf2*(double)S + vf2*(double)S*(double)S);
	
	return F2;
}

/*Frequency tests from Achaz Genetics 2009: outgroup*/
double freqtesto_achaz(int sample_size,long int *fr,int singleton,double *w1,double *w2) /* nomes outgroup */
{
    int i,j,ss;
    double Th1,Th2,Test,Thw,Thw2,*ww;
	double sumw1,sumw2,sumww,omi;
	double alfan,betan,alfat,betat;
    double betan2,betat2;
    long int it,nit;
   
    if(sample_size < 2) return(-10000);
	
	ww = (double *)calloc(sample_size,sizeof(double));
    i=an(sample_size);/*just calculate all an*/
    
    Th1 = 0.;
	sumw1 = 0.;
    for(i=1;i<sample_size;i++) {
		Th1 += ((double)*(fr+i))*((double)i*(double)w1[i]);
		sumw1 += w1[i];
	}
    Th1 /= sumw1;

	Th2 = 0.;
	sumw2 = 0.;
    for(i=1;i<sample_size;i++) {
		Th2 += ((double)*(fr+i))*((double)i*(double)w2[i]);
		sumw2 += w2[i];
	}
    Th2 /= sumw2;
	
    if(Th1 == 0. && Th2 == 0.) {
        free(ww);
        return(-10000);
    }

	Thw = 0.;
	sumww = 0.;
    for(i=1;i<sample_size;i++) {
		if(i==1) ss = singleton;
		else ss = 1;
		ww[i] = (double)1/(double)i * (double)ss;
		Thw += (double)*(fr+i)*(double)i*ww[i];
		sumww += ww[i];
	}
    Thw /= sumww;
	
	alfan = 0.;
	for(i=1;i<sample_size;i++) {
		omi = omegai(sample_size,i,w1,w2,sumw1,sumw2);
		alfan += i*(omi*omi);
	}
	
    if(sample_size > SAMPLE_LARGE){
        betan = 0.;
        betan2 = 0;
        nit = maxd(1,(double)NITER/((double)sample_size-1));
        for(i=1;i<sample_size;i++) {
            omi = omegai(sample_size,i,w1,w2,sumw1,sumw2);
            betan += i*i * (omi*omi) * sigmaii(sample_size,i);
            for(it=0;it<nit;it++) {
                do{j=(int)floor(ran1() * (sample_size-1))+1;}while(i==j);
                betan2 += 2.0 * i*j * omegai(sample_size,i,w1,w2,sumw1,sumw2) * omegai(sample_size,j,w1,w2,sumw1,sumw2) * sigmaij(sample_size,j,i);
            }
        }
        betan2 /= nit*(sample_size-1);
        betan2 *= ((sample_size-1)*(sample_size-2)/2.0)-sample_size-1;
        betan  += betan2;
    }
    else {
        betan = 0.;
        for(i=1;i<sample_size;i++) {
            omi = omegai(sample_size,i,w1,w2,sumw1,sumw2);
            betan += (double)i*(double)i * (omi*omi) * sigmaii(sample_size,i);
            /*printf("\nbetan=%f\tsigmaii[n=%d,i=%d]=%f",betan,sample_size,i,sigmaii(sample_size,i));*/
            for(j=i+1;j<sample_size;j++) {
                betan += 2.0 * (double)i*(double)j * omegai(sample_size,i,w1,w2,sumw1,sumw2) * omegai(sample_size,j,w1,w2,sumw1,sumw2) * sigmaij(sample_size,j,i);
                /*printf("\nbetan=%f\tsigmaij[n=%d,i=%d,j=%d]=%f",betan,sample_size,i,j,sigmaij(sample_size,j,i));*/
            }
        }
    }
	/*Theta2*/
	alfat = 0.;
	for(i=1;i<sample_size;i++) {
		alfat += (ww[i]/sumww * ww[i]/sumww)*i;
	}	
    if(sample_size > SAMPLE_LARGE){
        betat = 0.;
        betat2 = 0.;
        nit = maxd(1,(double)NITER/((double)sample_size-1));
        for(i=1;i<sample_size;i++) {
            betat += i*i * (ww[i]/sumww * ww[i]/sumww) * sigmaii(sample_size,i);
            for(it=0;it<nit;it++) {
                do{j=(int)floor(ran1() * (sample_size-1))+1;}while(i==j);
                betat2 += 2.0 * i*j * ww[i]/sumww * ww[j]/sumww * sigmaij(sample_size,j,i);
            }
        }
        betat2 /= nit*(sample_size-1);
        betat2 *= ((sample_size-1)*(sample_size-2)/2.0)-sample_size-1;
        betat  += betat2;
    }
    else {
        betat = 0.;
        for(i=1;i<sample_size;i++) {
            betat += i*i * (ww[i]/sumww * ww[i]/sumww) * sigmaii(sample_size,i);
            for(j=i+1;j<sample_size;j++) {
                betat += 2.0 * i*j * ww[i]/sumww * ww[j]/sumww * sigmaij(sample_size,j,i);
            }
        }
    }
	Thw2 = (Thw*Thw - alfat*Thw)/(1.0 + betat);
	
	/*Test*/
    free(ww);
    if((sqrt(alfan*Thw + betan*Thw2)) == 0)
        return -10000;
    Test = (Th1 - Th2)/(sqrt(alfan*Thw + betan*Thw2));
	
	if (fabs(Test) < 1.0E-15)
		return 0.0;

    return Test;
}
/*Frequency tests from Achaz Genetics 2009: NO outgroup*/
double freqtestn_achaz(int sample_size,long int *fr,int singleton,double *w1,double *w2) /* NO outgroup */
{
    int i,j,ss;
    double Th1,Th2,Test,Thw,Thw2,*ww;
	double sumw1,sumw2,sumww,omi,omj,psi,psj;
	double alfan,betan,alfat,betat;
    double betan2,betat2;
    long int it,nit;
   
    if(sample_size < 2) return(-10000);
	
	ww = (double *)calloc(sample_size,sizeof(double));
    i=an(sample_size);/*just calculate all an*/

    Th1 = 0.;
	sumw1 = 0.;
    for(i=1;i<=floor(sample_size/2);i++) {
		Th1 += ((double)*(fr+i))*(double)w1[i]/((double)psii(sample_size,i));
		sumw1 += w1[i];
	}
    Th1 /= sumw1;
	
	Th2 = 0.;
	sumw2 = 0.;
    for(i=1;i<=floor(sample_size/2);i++) {
		Th2 += ((double)*(fr+i))*(double)w2[i]/((double)psii(sample_size,i));
		sumw2 += w2[i];
	}
    Th2 /= sumw2;
	
    if(Th1 == 0. && Th2 == 0.) {
        free(ww);
        return(-10000);
    }
	
	Thw = 0.;
	sumww = 0.;
	for(i=1;i<=floor(sample_size/2);i++) {
		if(i==1) ss = singleton;
		else ss = 1;
		if(i == sample_size-i) ww[i] = (double)sample_size/((double)i*(double)(sample_size - i)*2.0)*(double)ss;
		else ww[i] = (double)sample_size/((double)i*(double)(sample_size - i)*1.0)*(double)ss;
		Thw += ((double)*(fr+i))*ww[i]/((double)psii(sample_size,i));
		sumww += ww[i];
	}
    Thw /= sumww;
	
	alfan = 0.;
	for(i=1;i<=floor(sample_size/2);i++) {
		omi = omegain(sample_size,i,w1,w2,sumw1,sumw2);
		psi = psii(sample_size,i);
		alfan += (omi*omi)/psi;
	}
	
    if(floor(sample_size/*/2.*/) > SAMPLE_LARGE) {
        betan = 0.;
        betan2 = 0.;
        nit = maxd(1,(double)NITER/((double)sample_size/2.));
        for(i=1;i<=floor(sample_size/2);i++) {
            omi = omegain(sample_size,i,w1,w2,sumw1,sumw2);
            psi = psii(sample_size,i);
            betan += omi/psi * omi/psi * rhoii(sample_size,i);
            for(it=0;it<nit;it++) {
                do{j=(int)floor(ran1() * (sample_size/2.))+1;}while(i==j);
                omj = omegain(sample_size,j,w1,w2,sumw1,sumw2);
                psj = psii(sample_size,j);
                betan2 += 2.0 * omi/psi * omj/psj * rhoij(sample_size,j,i);
            }
        }
        betan2 /= nit*(sample_size/2.);
        betan2 *= ((sample_size/2.)*(sample_size/2.-1)/2.0)-sample_size/2.;
        betan  += betan2;
    }
    else {
        betan = 0.;
        for(i=1;i<=floor(sample_size/2);i++) {
            omi = omegain(sample_size,i,w1,w2,sumw1,sumw2);
            psi = psii(sample_size,i);
            betan += omi/psi * omi/psi * rhoii(sample_size,i);
            for(j=i+1;j<=floor(sample_size/2);j++) {
                omj = omegain(sample_size,j,w1,w2,sumw1,sumw2);
                psj = psii(sample_size,j);
                betan += 2.0 * omi/psi * omj/psj * rhoij(sample_size,j,i);
                /*printf("\nrhoij[n=%d,j=%d,i=%d] = %f",sample_size,j,i,rhoij(sample_size,j,i));*/
            }
        }
    }
	/*Theta2*/
	alfat = 0.;
	for(i=1;i<=floor(sample_size/2);i++) {
		psi = psii(sample_size,i);
		alfat += (ww[i]/sumww * ww[i]/sumww)/psi;
	}	
    if(floor(sample_size/*/2.*/) > SAMPLE_LARGE) {
        betat = 0.;
        betat2 = 0.;
        nit = maxd(1,(double)NITER/((double)sample_size/2.));
        for(i=1;i<=floor(sample_size/2);i++) {
            psi = psii(sample_size,i);
            betat += (ww[i]/sumww)/psi * (ww[i]/sumww)/psi * rhoii(sample_size,i);
            for(it=0;it<nit;it++) {
                do{j=(int)floor(ran1() * (sample_size/2.))+1;}while(i==j);
                psj = psii(sample_size,j);
                betat2 += 2.0 * (ww[i]/sumww)/psi * (ww[j]/sumww)/psj * rhoij(sample_size,j,i);
            }
        }
        betat2 /= nit*(sample_size/2.);
        betat2 *= ((sample_size/2.)*(sample_size/2.-1)/2.0)-sample_size/2.;
        betat  += betat2;
    }
    else {
        betat = 0.;
        for(i=1;i<=floor(sample_size/2);i++) {
            psi = psii(sample_size,i);
            betat += (ww[i]/sumww)/psi * (ww[i]/sumww)/psi * rhoii(sample_size,i);
            for(j=i+1;j<=floor(sample_size/2);j++) {
                psj = psii(sample_size,j);
                betat += 2.0 * (ww[i]/sumww)/psi * (ww[j]/sumww)/psj * rhoij(sample_size,j,i);
            }
        }
    }
	Thw2 = (Thw*Thw - alfat*Thw)/(1.0 + betat);

    free(ww);

    /*Test*/
    if((sqrt(alfan*Thw + betan*Thw2)) == 0)
        return -10000;
	Test = (Th1 - Th2)/(sqrt(alfan*Thw + betan*Thw2));
						
	if (fabs(Test) < 1.0E-15)
		return 0.0;

    return Test;
}
/*Achaz code for debugging*//*************************************************************************************//*
double compute_Omega_Y( int n, long int *Xi_array, char norm ){
	double *w1,
	*w2;
	int n_xi=n-1;
	int i;
	double T=0;
	double Estim_theta; 
	double Estim_theta2;
	w1 = (double *) malloc( (size_t) n_xi*sizeof(double) );
	w2 = (double *) malloc( (size_t) n_xi*sizeof(double) );
	if( !w1 || !w2 )
		fprintf(stderr, "compute_Omega_Y: cannot allocate memory for w1 or w2, bye\n"), exit(3);
	w1[0]=w2[0]=0;
	for(i=2;i<=n_xi; i++){
		w1[i-1] = (n-i+0.0);
		w2[i-1] = 1.0/i;
	}
	compute_theta2estimator_omega(  n,  Xi_array, w2, &Estim_theta, &Estim_theta2  );
	T = compute_T( n, Xi_array, w1, w2, norm, Estim_theta, Estim_theta2 );
    free(w1);
    free(w2);
    return T;
}
double compute_T( int n, long int *Xi_array, double *w1, double *w2, char norm, double theta, double theta2 ){
	double Sum_w1=0.0;
	double Sum_w2=0.0;
	int i,j;
	double alpha=0,
	beta=0;
	double T=0;
	long int S=0;
	double *HarmonicSums=NULL;
	HarmonicSums = compute_HarmonicSums( n );
	for(i=0; i<n-1; i++)
		S += Xi_array[i+1];
	if ( S == 0 )
		return 0;
	for(i=0;i<n-1;i++){
		Sum_w1 += w1[i];
		Sum_w2 += w2[i];
	}
	if(Sum_w1 == 0 || Sum_w2 == 0)return 10000;
	for(i=0; i < n-1 ; i++)
		T += (double)(i+1.0) * Xi_array[i+1] * ( (w1[i]/Sum_w1) - (w2[i]/Sum_w2) );
	if(norm==0)
		return T;
	alpha=0;
	beta=0;
	for(i=0; i < n-1 ; i++){	
		double Omega = ( (w1[i]/Sum_w1) - (w2[i]/Sum_w2) );	
		alpha += Omega*Omega*(i+1.0);
	}
	beta=0;
	for(i=1; i < n ; i++){
		for( j=i; j< n ; j++){
			double Omega_i = ( (w1[i-1]/Sum_w1) - (w2[i-1]/Sum_w2) );
			double Omega_j = ( (w1[j-1]/Sum_w1) - (w2[j-1]/Sum_w2) );
			if( i == j)
				beta  +=  i*i* Omega_i*Omega_i * sigma_ii_FU1995( i,  HarmonicSums , n   );
			else
				beta  +=  2*i*j* Omega_i*Omega_j * sigma_ij_FU1995( i, j,  HarmonicSums , n   );		
		}
	}
	T /= sqrt( alpha*theta + beta*theta2   );
	free(HarmonicSums);
	return T;
}
void compute_theta2estimator_omega(  int n,  long int *Xi_array, double *w, double *theta, double *theta2  ){
	double *HarmonicSums=NULL; 
	int i,j;
	double Sum_w=0.0;
	double alpha=0,
	beta=0;
	for(i=0;i<n-1;i++)
		Sum_w += w[i];
	HarmonicSums = compute_HarmonicSums( n );
	*theta=0;
	for(i=0; i < n-1 ; i++){
		*theta +=  Xi_array[i+1] * (i+1.0)* w[i]/Sum_w;
	}
	alpha=0;
	for(i=0; i < n-1 ; i++){
		double omega = w[i]/Sum_w;
		alpha += omega*omega*(i+1.0);
	}
	beta=0;
	for(i=1; i < n ; i++){
		for( j=i; j< n ; j++){
			double omega_i = (w[i-1]/Sum_w) ;
			double omega_j = (w[j-1]/Sum_w) ;
			if( i == j)
				beta  +=  i*i* omega_i*omega_i * sigma_ii_FU1995( i,  HarmonicSums , n   );
			else
				beta  +=  2*i*j* omega_i*omega_j * sigma_ij_FU1995( i, j,  HarmonicSums , n   );
		}
	}
	free(HarmonicSums);
	*theta2 = ( pow( *theta, 2) - alpha* (*theta) )/ (1.0 + beta);
}
double * compute_HarmonicSums( int n ){
	double *HarmonicSums;
	int i;
	HarmonicSums = (double *)malloc( sizeof(double)*n );
	if(!HarmonicSums)
		fprintf(stderr, "compute_HarmonicSums: malloc error for HarmonicSums\n"), exit(3);
	HarmonicSums[0]=0;
	i=1;
	while(i<n){
		HarmonicSums[i] = HarmonicSums[i-1] + 1.0/i;
		i++;
	}
	return HarmonicSums;
}
double beta_FU1995( int i, double *HarmonicSums, int n  ){
	double ai = HarmonicSums[i-1],
	an = HarmonicSums[n-1];
	double beta=0;
	double nloci= (double)n;
	beta = 2.0 * nloci * ( an + (1.0/nloci) - ai )/(  (nloci-i+1.0 )*(nloci-i) ) - 2.0/(nloci - i);
	return beta;
}
double sigma_ii_FU1995( int i, double *HarmonicSums, int n ){
	double nloci= (double)n;
	double sigma_ii=0;
	double ai = HarmonicSums[i-1],
	an = HarmonicSums[n-1];
	if( 2*i < n )
	{
		sigma_ii = beta_FU1995( i+1, HarmonicSums, n  );
	}
	else
	{
		if( 2*i == n  )
		{
			sigma_ii = 2.0*(an-ai)/(nloci - i) - 1.0/(i*i);
		}
		else
		{
			sigma_ii = beta_FU1995( i , HarmonicSums, n  ) - 1.0/(i*i);
		}
	}
	return sigma_ii;
}
double sigma_ij_FU1995( int i, int j, double *HarmonicSums, int n ){
	double nloci= (double)n;
	double sigma_ij=0;
	if(i==j){
		return sigma_ii_FU1995( i, HarmonicSums, n );
	}
 	if(i<j){
		int tmp=i;
		i=j;
		j=tmp;
	}
	double  ai=HarmonicSums[i-1],
	aj=HarmonicSums[j-1],
	an=HarmonicSums[n-1];
	if( i+j < n )
	{
		sigma_ij = ( beta_FU1995( i+1, HarmonicSums, n  ) -  beta_FU1995( i, HarmonicSums, n  ) ) / 2.0;
	}
	else
	{
		if( i+j == n  )
		{
			sigma_ij  =  ((an - ai)/(nloci - i) + (an - aj)/(nloci - j)) 
			- ( ( beta_FU1995( i, HarmonicSums, n   ) +  beta_FU1995( j+1 , HarmonicSums, n ) )  / 2.0 ) 
			- (1.0/(i*j));
		}
		else
		{
			sigma_ij = (( beta_FU1995( j, HarmonicSums, n  ) -  beta_FU1995( j+1, HarmonicSums, n ) )/2.0) - (1.0/(i*j));
		}
	}
	return sigma_ij;
}
*//***********************************************************************************************/
double an(int n)
{
    static double *an_m=0;
    static long int nsam=0;
    long int i;
    
    if(nsam == 0) {
        an_m = (double *)calloc(n+2, sizeof(double));
        for(i=1;i<n+1;i++)
            an_m[i+1] = an_m[i] + 1.0/(double)i;
        nsam = n;
    } else {
        if(nsam < n) {
            an_m = (double *)realloc(an_m,(n+2)*sizeof(double));
            for(i=nsam+1;i<n+1;i++)
                an_m[i+1] = an_m[i] + 1.0/(double)i;
            nsam = n;
        }
    }
    return an_m[n];
}
double a2n(int n)
{
	double a2ni = 0.;
	int i;
	for(i=1;i<n;i++)
		a2ni += (1.0/((double)i*(double)i));
	return a2ni;
}
				   
double bn(int n,int i)
{
	double an(int);
	double bni;
	bni = (2.0*(double)n)/((double)(n-i+1)*(double)(n-i)) * (an(n+1) - an(i)) - 2.0/((double)(n-i));
	return bni;
}
double sigmaii(int n,int i)
{
	double an(int),bn(int,int);
	double sii=0;
	if(2*i<n)
		sii = bn(n,i+1);
	if(2*i==n)
		sii = 2.0 * (an(n)-an(i))/((double)(n-i)) - 1.0/((double)i*(double)i);
	if(2*i>n)
		sii = bn(n,i) - 1.0/((double)i*(double)i);
	return sii;
}
double sigmaij(int n,int i,int j)
{
	double an(int),bn(int,int);
	double sigmaii(int,int);
	double sij=0.0;
	int ii,jj;
	if(i < j) {
		ii = j;
		jj = i;
	}
	else {
		if(i==j) {
			return(sigmaii(n,i));
		}
		else {
			ii = i;
			jj = j;
		}
	}
	if((ii+jj)<n) 
		sij = (bn(n,ii+1) - bn(n,ii))/2.0;
	if((ii+jj)==n) 
		sij = (an(n)-an(ii))/((double)(n-ii)) + (an(n)-an(jj))/((double)(n-jj)) - (bn(n,ii) + bn(n,jj+1))/2.0 - 1.0/((double)ii*(double)jj);
	if((ii+jj)>n) 
		sij = (bn(n,jj) - bn(n,jj+1))/2.0 - 1.0/((double)ii*(double)jj);
	return sij;
}
double omegai(int n,int i,double *w1,double *w2,double sumw1,double sumw2)
{
    double omi;
    /*
     long int x;
     double sumw1,sumw2;
     
     sumw1=0.0;
     for(x=1;x<n;x++) sumw1 += w1[x];
     sumw2=0.0;
     for(x=1;x<n;x++) sumw2 += w2[x];
    */
    omi = w1[i]/sumw1 - w2[i]/sumw2;
    return omi;
}
/*
double omegai(int n,int i,double *w1,double *w2)
{
	double omi;
	int x;
	double sumw1,sumw2;
	
	sumw1=0.0;
	for(x=1;x<n;x++) sumw1 += w1[x];
	sumw2=0.0;
	for(x=1;x<n;x++) sumw2 += w2[x];
	omi = w1[i]/sumw1 - w2[i]/sumw2;
	return omi;
}
*/
double psii(int n,int i)
{
	double psi;
	int krond;
	
	if((int)i==(int)(n-i)) krond = 1;
	else krond = 0;
	psi= (double)n/((double)(1.+krond)*(double)i*(double)(n-i));
	return psi;
}
double rhoii(int n,int i)
{
	double sigmaii(int,int);
	double sigmaij(int,int,int);
	double rhoi;
	int krond;
	
	if((int)i==(int)(n-i)) 
		krond = 1;
	else 
		krond = 0;
	
	rhoi  = (sigmaii(n,i)+sigmaii(n,n-i)+2.0*sigmaij(n,i,n-i));
	rhoi /= ((1.0 + krond) * (1.0 + krond));
	return rhoi;
}
double rhoij(int n,int i,int j)
{
	double sigmaii(int,int);
	double sigmaij(int,int,int);
	double rhoj;
	int krondi,krondj;
	
	if(i==(n-i)) krondi = 1;
	else krondi = 0;
	if(j==(n-j)) krondj = 1;
	else krondj = 0;
	
	rhoj  = (sigmaij(n,i,j)+sigmaij(n,i,n-j)+sigmaij(n,n-i,j)+sigmaij(n,n-i,n-j));
	rhoj /= ((1.0 + krondi) * (1.0 + krondj));
	return rhoj;
}
double omegain(int n,int i,double *w1,double *w2,double sumw1,double sumw2)
{
    double omi;
    /*
     long int x;
     double sumw1,sumw2;
     
     sumw1=0.0;
     for(x=1;x<=floor(n/2);x++) sumw1 += w1[x];
     sumw2=0.0;
     for(x=1;x<=floor(n/2);x++) sumw2 += w2[x];
     */
    omi = w1[i]/sumw1 - w2[i]/sumw2;
    return omi;
}
/*
double omegain(int n,int i,double *w1,double *w2)
{
	double omi;
	int x;
	double sumw1,sumw2;
	
	sumw1=0.0;
	for(x=1;x<=floor(n/2);x++) sumw1 += w1[x];
	sumw2=0.0;
	for(x=1;x<=floor(n/2);x++) sumw2 += w2[x];
	omi = w1[i]/sumw1 - w2[i]/sumw2;
	return omi;
}
*/
/*weigths for each test suming 0 per position (each wi,nx sum 1 per position)*/
int ominx_tajD(double **ominx, int **nx, int *no, double any, int nsam, long int length, long int S)
{
	long int j;
	int i,n;
	/**/double sum1;/**/
	
	for(j=0;j<length;j++) {
		if(no[j]) {
            n=0;
            for(i=0;i<nsam;i++) {n += nx[i][j];} /*number of samples*/
            /**/sum1 = 0.;/**/
            /**/for(i=1;i<n;i++) sum1 += (double)1./(double)i;/**/
            for(i=1;i<n;i++) {/*frequencies*/
                ominx[i][j] = 2.*((double)n-(double)i)/((double)n*((double)n-1.)) - (1./((double)i))/sum1;
                /*ominx[i][j] = 2.*((double)n-(double)i)/((double)n*((double)n-1.)) - 1./((double)i*any/(double)S);*/
            }
        }
	}
	return(1);
}
int ominx_Dfl(double **ominx, int **nx, int *no, double any, int nsam, long int length, long int S) 
{
	long int j;
	int i,n;
	/**/double sum1;/**/

	for(j=0;j<length;j++) {
		if(no[j]) {
			n=0;
			for(i=0;i<nsam;i++) {n += nx[i][j];}
			/**/sum1=0.;/**/
			/**/for(i=1;i<n;i++) sum1 += (double)1./(double)i;/**/

			if(n>1) {
				/*ominx[1][j] = 1./((double)1.*any/(double)S) - 1.;*/
				/**/ominx[1][j] = (1./((double)1))/sum1 - 1.;/**/
			}
			for(i=2;i<n;i++) {
				/*ominx[i][j] = 1./((double)i*any/(double)S) - 0.;*/
				/**/ominx[i][j] = (1./((double)i))/sum1 - 0.;/**/
			}
		}
	}
	return(1);
}
int ominx_Ffl(double **ominx, int **nx, int *no, int nsam, long int length, long int S) 
{
	long int j;
	int i,n;
	
	for(j=0;j<length;j++) {
		if(no[j]) {
			n=0;
			for(i=0;i<nsam;i++) {n += nx[i][j];}
			for(i=1;i<n;i++) {
				if(i==1) {
					ominx[i][j] = 2.*((double)n-(double)i)/((double)n*((double)n-1.)) - 1.;
				}
				else {
					ominx[i][j] = 2.*((double)n-(double)i)/((double)n*((double)n-1.)) - 0.;
				}
			}
		}
	}
	return(1);
}
int ominx_Hnfw(double **ominx, int **nx, int *no, int nsam, long int length, long int S) 
{
	long int j;
	int i,n;
	for(j=0;j<length;j++) {
		if(no[j]) {
			n=0;
			for(i=0;i<nsam;i++) {n += nx[i][j];}
			for(i=1;i<n;i++) {
				ominx[i][j] = 2.*((double)n-(double)i)/((double)n*((double)n-1.)) - (2.*(double)i)/((double)n*((double)n-1.));
			}
		}
	}
	return(1);
}
int ominx_Ezeng(double **ominx, int **nx, int *no, double any, int nsam, long int length, long int S) 
{
	long int j;
	int i,n;
	/**/double sum1;/**/
	
	for(j=0;j<length;j++) {
		if(no[j]) {
			n=0;
			for(i=0;i<nsam;i++) {n += nx[i][j];}
			/**/sum1=0.;/**/
			/**/for(i=1;i<n;i++) sum1 += (double)1./(double)i;/**/
			for(i=1;i<n;i++) {
				/*ominx[i][j] = 1./(double)(n-1.) - 1./((double)i*any/(double)S);*/
				/**/ominx[i][j] = 1./(double)(n-1.) - (1./((double)i))/sum1;/**/
			}
		}
	}
	return(1);
}
int ominx_Yachaz(double **ominx, int **nx, int *no, double any, int nsam, long int length, long int S) 
{
	long int j;
	int i,n;
	double sum1;
	for(j=0;j<length;j++) {
		if(no[j]) {
			n=0;
			for(i=0;i<nsam;i++) {n += nx[i][j];} /*number of samples*/
			sum1=0.;
			for(i=2;i<n;i++) {sum1 += (double)1./(double)i;}
			for(i=1;i<n;i++) {/*frequencies*/
				if(i!=1) {
					/*ominx[i][j] = 2.*(double)(n-i)/((double)(n-1.)*(double)(n-2.)) - 1./((double)i*(any/(double)S-1.));*/
					/**/ominx[i][j] = 2.*(double)(n-i)/((double)(n-1.)*(double)(n-2.)) - (1./((double)i))/sum1;/**/
				}
				else {
					ominx[i][j] = 0.;
				}
			}
		}
	}
	return(1);
}
int ominx_FH(double **ominx, int **nx, int *no, double any, int nsam, long int length, long int Sco)
{
    long int j;
    int i,n;
    double sum1;
    
    for(j=0;j<length;j++) {
        if(no[j]) {
            n=0;
            for(i=0;i<nsam;i++) {n += nx[i][j];}
            sum1=0.;
            for(i=1;i<n;i++) sum1 += (double)1./(double)i;
            for(i=1;i<n;i++) {
                ominx[i][j] = (1./((double)i))/sum1 - (2.*(double)i)/((double)n*((double)n-1.));
            }
        }
    }
    return(1);
}
/*Optimal*/
int ominx_Optimal_unfolded(double **ominx, int **nx, int *no, double any, int nsam, long int length, long int S, double *mean_freqsptr)
{
	long int j;
	int i,n;
	/**/double sum1;/**/
	
	for(j=0;j<length;j++) {
		if(no[j]) {
			n=0;
			for(i=0;i<nsam;i++) {n += nx[i][j];} /*number of samples*/
			/**/sum1 = 0.;/**/
			/**/for(i=1;i<n;i++) sum1 += (double)1./(double)i;/**/
			for(i=1;i<n;i++) {/*frequencies*/
				ominx[i][j] = mean_freqsptr[i] - (1./((double)i))/sum1;
			}
		}
	}
	return(1);
}
/*NO OUTGROUP*/
int ominx_Optimal_folded(double **ominx, int **nx, double any, int nsam, long int length, long int S, double *mean_freqsptr)
{
	long int j;
	int i,n;
	double sum1;
	
	for(j=0;j<length;j++) {
		n=0;
		for(i=0;i<nsam;i++) {n += nx[i][j];} /*number of samples*/
		sum1 = 0.;
		for(i=2;i<=floor((double)n/2.);i++) {
			if(i==n-i) {
				sum1 += (double)n/((double)i*(double)(n-i)*2.);
			}
			else {
				sum1 += (double)n/((double)i*(double)(n-i));
			}
		}
		for(i=1;i<=floor((double)n/2.);i++) {
			if(i==n-i) {
				ominx[i][j] = mean_freqsptr[i] - ((double)n/(i*(double)(n-i)*2.))/sum1;
			}
			else {
				ominx[i][j] = mean_freqsptr[i] - ((double)n/(i*(double)(n-i)*1.))/sum1;
			}
		}
	}
	return(1);
}
int ominx_tajDnooutg(double **ominx, int **nx, double any, int nsam, long int length, long int S)
{
	long int j;
	int i,n;
	/**/double sum1;/**/
	//double wt=0.;
	for(j=0;j<length;j++) {
		n=0;
		for(i=0;i<nsam;i++) {n += nx[i][j];} /*number of samples*/
		/**/sum1 = 0.;/**/
		/**/for(i=1;i<n;i++) sum1 += (double)1./(double)i;/**/
		for(i=1;i<n;i++) {/*frequencies*/
			/**/ominx[i][j] = 2.*((double)n-(double)i)/((double)n*((double)n-1.)) - (1./((double)i))/sum1;/**/
			/*  ominx[i][j] = 2.*((double)n-(double)i)/((double)n*((double)n-1.)) -  1./((double)i*any/(double)S);*/
		}
		for(i=1;i<=floor((double)n/2.);i++) {/*from Achaz 2009*/
			if(i==n-i) ominx[i][j] *= (double)n/((double)(n-i)*2.);
			else ominx[i][j] *= (double)n/((double)(n-i)*1.);
			//wt+=ominx[i][j];
		}
	}
	return(1);
}
int ominx_Dfl2(double **ominx, int **nx, double any, int nsam, long int length, long int S) 
{
	long int j;
	int i,n;
	double sum1;
	
	for(j=0;j<length;j++) {
		n=0;
		for(i=0;i<nsam;i++) {n += nx[i][j];} /*number of samples*/
		sum1 = 0.;
		for(i=1;i<=floor((double)n/2.);i++) {
			if(i==n-i) {sum1 += (double)n/((double)i*(double)(n-i)*2.);}
			else {sum1 += (double)n/((double)i*(double)(n-i));}
		}
		for(i=1;i<=floor((double)n/2.);i++) {
			if(i==1) {
				if(i==n-i) {
					ominx[i][j] = ((double)n/((double)i*(double)(n-i)*2.))/sum1 - 1.;
				}
				else {
					ominx[i][j] = ((double)n/((double)i*(double)(n-i)*1.))/sum1 - 1.;
				}
			}
			else {
				if(i==n-i) {
					ominx[i][j] = ((double)n/((double)i*(double)(n-i)*2.))/sum1 - 0.;
				}
				else {
					ominx[i][j] = ((double)n/((double)i*(double)(n-i)*1.))/sum1 - 0.;
				}
			}
		}
	}
	return(1);
}
int ominx_Ffl2(double **ominx, int **nx, int nsam, long int length, long int S) 
{
	long int j;
	int i,n;
	double sum1;
	
	for(j=0;j<length;j++) {
		n=0;
		for(i=0;i<nsam;i++) {n += nx[i][j];} /*number of samples*/
		sum1 = 0.;
		for(i=1;i<=floor((double)n/2.);i++) {
			if(i==n-i) {sum1 += (double)n/2.;}
			else {sum1 += (double)n/1.;}
		}
		for(i=1;i<=floor((double)n/2.);i++) {
			if(i==1) {
				if(i==n-i) {
					ominx[i][j] = ((double)n/2.)/sum1 - 1.;
				}
				else {
					ominx[i][j] = ((double)n/1.)/sum1 - 1.;
				}
			}
			else {
				if(i==n-i) {
					ominx[i][j] = ((double)n/2.)/sum1 - 0.;
				}
				else {
					ominx[i][j] = ((double)n/1.)/sum1 - 0.;
				}
			}
		}
	}
	return(1);
}
int ominx_Yachaz2(double **ominx, int **nx, double any, int nsam, long int length, long int S) 
{
	long int j;
	int i,n;
	double sum1,sum2;
	
	for(j=0;j<length;j++) {
		n=0;
		for(i=0;i<nsam;i++) {n += nx[i][j];} /*number of samples*/
		sum1 = 0.;
		sum2 = 0.;
		for(i=2;i<=floor((double)n/2.);i++) {
			if(i==n-i) {
				sum1 += (double)n/2.;
				sum2 += (double)n/((double)i*(double)(n-i)*2.);
			}
			else {
				sum1 += (double)n;
				sum2 += (double)n/((double)i*(double)(n-i));
			}
		}
		for(i=1;i<=floor((double)n/2.);i++) {
			if(i!=1) {
				if(i==n-i) {
					ominx[i][j] = ((double)n/2.)/sum1 - ((double)n/(i*(double)(n-i)*2.))/sum2;
				}
				else {
					ominx[i][j] = ((double)n/1.)/sum1 - ((double)n/(i*(double)(n-i)*1.))/sum2;
				}
			}
			else ominx[i][j] = 0.;
		}
	}
	return(1);
}

double freqtest_outg_missing(double **ominx,long int **eix,int **nx, int *no, double thetaw, int nsam, long int length,int ncomp, struct covar_ij **Covi1i2, long int *count_cov, int singleton, double *bnxp, double *theta_square) 
{
	double Num,Num0,Den,Den2,Den12;
	long int j,k,ll,*llv,h,x,ss;
	int i0,i,ii,i1,i2,i12,n1,n2,n12;
	struct covar_m {
		int n1;
		int n2;
		int n12;
		double Den12;
		double Cov12;
	} *Denm12;
	long int count;
	double anx,est_theta2;
	double betan,omi,omj;
	double SumCov,Cov12,CovMt;
	double alfat,betat,*ww,sumww;
	
    CovMt = 0.;
    if(nsam < 2) return(-10000);

	if((llv = (long int *)calloc(length,sizeof(long int))) == 0) {
		printf("Error calloc memory freqtest0");
		exit(1);
	}
	if((Denm12 = (struct covar_m *)calloc(1,sizeof(struct covar_m))) == 0) {
		printf("Error calloc memory freqtest1");
		exit(1);
	}
	count = 0;
	
	Num = 0.;
	ll = 0;
	n1 = 0;
	for(j=0;j<length;j++) {
		if(no[j]) {
			n1=0;
			for(i=0;i<nsam;i++) {n1 += nx[i][j];} /*number of samples*/
			Num0 = 0;
			for(i=1;i<n1;i++) {
				if(eix[i][j]) {
					Num0 += (double)i*ominx[i][j]*eix[i][j];
					break;
				}
			}
			Num += Num0;
			if(Num0) {llv[ll] = j;ll += 1;}
		}
	}
    if(Num == 0.) {
        free(llv);
        free(Denm12);
        return(-10000);
    }
	
	/*calculate the first term of the variance*/
	Den = 0.;
	for(j=0;j<ll;j++) {
		if(no[j]) {
			n1=0;
			for(i=0;i<nsam;i++) {n1 += nx[i][llv[j]];} /*number of samples*/
			for(i=1;i<n1;i++) 
				Den  += (double)i*ominx[i][llv[j]]*ominx[i][llv[j]]; /*first part of the variance*/
		}
	}
	Den *= thetaw/(double)ll;
	
	/*calculate the covariance: ncomp0=fast_approach(incorrect), ncomp_large=total OR a ncomp_small=random set of comparisons*/
	Den2 = 0.;alfat = 0.;betat = betan = 0.;
	if(ncomp == 0) {
		ww = (double *)calloc(nsam,sizeof(double));
		for(j=0;j<ll;j++) {
			n1 = 0;
			for(i=0;i<nsam;i++) {n1 += nx[i][llv[j]];} /*number of samples*/
			alfat = 0.;
			sumww = 0.;
			for(i=1;i<n1;i++) {
				if(i==1) ss = singleton;
				else ss = 1;
				ww[i] = (double)1/(double)i * (double)ss;
				sumww += ww[i];
			}	
			alfat = 0.;
			betat = 0.;
			betan = 0.;
			for(i1=1;i1<n1;i1++) {
				alfat += (ww[i1]/sumww * ww[i1]/sumww)*i1;
				omi = ominx[i1][llv[j]];
				betan += i1*i1 * (omi*omi) * sigmaii(n1,i1);
				betat += i1*i1 * (ww[i1]/sumww * ww[i1]/sumww) * sigmaii(n1,i1);
				for(i2=i1+1;i2<n1;i2++) {
					omj = ominx[i2][llv[j]];
					betan += 2.0 * i1*i2 * omi * omj * sigmaij(n1,i1,i2);
					betat += 2.0 * i1*i2 * ww[i1]/sumww * ww[i2]/sumww * sigmaij(n1,i1,i2);
				}
			}
			Den2 += betan;
		}
		*bnxp  = Den2;
		
		est_theta2 = (thetaw*thetaw - alfat*thetaw)/(1.0 + betat);	
		Den2 *= est_theta2/ll;
		free(ww);

		*theta_square  = est_theta2;
	}
	else {
		anx = 0.;
		for(j=0;j<ll;j++) {
			n1=0; 
			for(i1=0;i1<nsam;i1++) {
				n1 += (double)nx[i1][llv[j]];
			} 
			for(i1=1;i1<n1;i1++) {
				anx+= 1.0/(double)i1;
			} 
			if(singleton == 0) anx -= n1/(double)(n1-1.0);
		}
		SumCov = 0.;
		if(ncomp >= ll*(ll-1)/2.) {
			for(j=0;j<ll;j++) {
				if(no[j]) {
					n1=0;
					for(i1=0;i1<nsam;i1++) {n1 += nx[i1][llv[j]];} /*number of samples*/
					for(k=j+1;k<ll;k++) {
						if(no[k]) {
							n2=0;
							for(i2=0;i2<nsam;i2++) {n2 += nx[i2][llv[k]];} 
							n12=0;
							for(i=0;i<nsam;i++) {n12 += nx[i][llv[j]] & nx[i][llv[k]];} 
							
							for(i=0;i<count;i++) {
								if(((Denm12[i].n1 == n1 && Denm12[i].n2 == n2) || (Denm12[i].n1 == n2 && Denm12[i].n2 == n1)) && Denm12[i].n12 == n12) {
									Den2 += Denm12[i].Den12;
									SumCov += Denm12[i].Cov12;
									break;
								}
							}
							if(i>=count) {
								Den12 = 0.;
								Cov12 = 0.;
								for(i1=1+!singleton;i1<n1;i1++) {
									for(i2=1+!singleton;i2<n2;i2++) {
										for(x=0;x<*count_cov;x++) {
											if((Covi1i2[0][x].i1==i1 && Covi1i2[0][x].n1==n1 && Covi1i2[0][x].i2==i2 && Covi1i2[0][x].n2==n2 && Covi1i2[0][x].n12==n12) ||
											   (Covi1i2[0][x].i1==i2 && Covi1i2[0][x].n1==n2 && Covi1i2[0][x].i2==i1 && Covi1i2[0][x].n2==n1 && Covi1i2[0][x].n12==n12)) {
												CovMt = Covi1i2[0][x].Covij;
												break;
											}
										}
										if(x>=*count_cov) {
											CovMt = (double)cov_missing(i1,i2,n1,n2,n12,(float)1);
											Covi1i2[0][*count_cov].i1 = i1;
											Covi1i2[0][*count_cov].i2 = i2;
											Covi1i2[0][*count_cov].n1 = n1;
											Covi1i2[0][*count_cov].n2 = n2;
											Covi1i2[0][*count_cov].n12 = n12;
											Covi1i2[0][*count_cov].Covij = CovMt;
											*count_cov += 1;
											if((Covi1i2[0] = (struct covar_ij *)realloc(Covi1i2[0],(count_cov[0]+1)*sizeof(struct covar_ij))) == 0) {
												printf("Error calloc memory freqtestl");
												exit(1);
											}
										}
										Cov12 += CovMt;
										Den12 += (double)i1 * (double)i2 * ominx[i1][llv[j]] * ominx[i2][llv[k]] * CovMt;
									}
								}
								Den12 /=(double)(ll-1.0)*ll;
								Denm12[count].n1 = n1;
								Denm12[count].n2 = n2;
								Denm12[count].n12 = n12;
								Denm12[count].Den12 = Den12;
								Denm12[count].Cov12 = Cov12;
								count += 1;
								if((Denm12 = (struct covar_m *)realloc(Denm12,(count+1)*sizeof(struct covar_m))) == 0) {
									printf("Error calloc memory freqtest2");
									exit(1);
								}
								Den2 += Den12;
								SumCov += Cov12;
							}
						}
					}
				}
			}
			Den2 *= 2.;
			SumCov *= 2.;
			SumCov /= ll*(ll-1.0);
			est_theta2 = ((double)ll*(double)ll-(double)ll)/(anx/(double)ll*anx/(double)ll+SumCov/(1.0*1.0));
			
			*bnxp  = Den2;
			
			Den2 *= est_theta2/(1.0*1.0);
			
			*theta_square  = est_theta2;
		}
		else {
			ii=0;i0=0;
			while(ii<ncomp && i0 < ncomp*10) {
				j=0;k=0;
				while(j==k) {
					j = (long int)floor(ran1()*ll);
					k = (long int)floor(ran1()*ll);
				}
				if(j>k) {h=k;k=j;j=h;}
				if(no[j]) {
					n1=0;
					for(i1=0;i1<nsam;i1++) {n1 += nx[i1][llv[j]];}
					n2=0;
					for(i2=0;i2<nsam;i2++) {n2 += nx[i2][llv[k]];} 
					n12=0;
					for(i12=0;i12<nsam;i12++) {n12 += nx[i12][llv[j]] & nx[i12][llv[k]];} 
					
					for(i=0;i<count;i++) {
						if(((Denm12[i].n1 == n1 && Denm12[i].n2 == n2) || (Denm12[i].n1 == n2 && Denm12[i].n2 == n1)) && Denm12[i].n12 == n12) {
							Den2 += Denm12[i].Den12;
							SumCov += Denm12[i].Cov12;
							break;
						}
					}
					if(i>=count) {
						Den12 = 0;
						Cov12 = 0.;
						for(i1=1+!singleton;i1<n1;i1++) {
							for(i2=1+!singleton;i2<n2;i2++) {
								for(x=0;x<*count_cov;x++) {
									if((Covi1i2[0][x].i1==i1 && Covi1i2[0][x].n1==n1 && Covi1i2[0][x].i2==i2 && Covi1i2[0][x].n2==n2 && Covi1i2[0][x].n12==n12) ||
									   (Covi1i2[0][x].i1==i2 && Covi1i2[0][x].n1==n2 && Covi1i2[0][x].i2==i1 && Covi1i2[0][x].n2==n1 && Covi1i2[0][x].n12==n12)) {
										CovMt = Covi1i2[0][x].Covij;
										break;
									}
								}
								if(x>=*count_cov) {
									CovMt = (double)cov_missing(i1,i2,n1,n2,n12,(float)1);
									Covi1i2[0][*count_cov].i1 = i1;
									Covi1i2[0][*count_cov].i2 = i2;
									Covi1i2[0][*count_cov].n1 = n1;
									Covi1i2[0][*count_cov].n2 = n2;
									Covi1i2[0][*count_cov].n12 = n12;
									Covi1i2[0][*count_cov].Covij = CovMt;
									*count_cov += 1;
									if((Covi1i2[0] = (struct covar_ij *)realloc(Covi1i2[0],(count_cov[0]+1)*sizeof(struct covar_ij))) == 0) {
										printf("Error calloc memory freqtestl");
										exit(1);
									}
								}
								Cov12 += CovMt;
								Den12 += (double)i1 * (double)i2 * ominx[i1][llv[j]] * ominx[i2][llv[k]] * CovMt;
							}
						}
						Den12 /=(double)(ll-1.0)*ll;
						Denm12[count].n1 = n1;
						Denm12[count].n2 = n2;
						Denm12[count].n12 = n12;
						Denm12[count].Den12 = Den12;
						Denm12[count].Cov12 = Cov12;
						count += 1;
						if((Denm12 = (struct covar_m *)realloc(Denm12,(count+1)*sizeof(struct covar_m))) == 0) {
							printf("Error calloc memory freqtest2");
							exit(1);
						}
						
						Den2 += Den12;
						SumCov += Cov12;
					}
					ii += 1;
				}
				i0+=1;
			} 
			Den2 = Den2 * ll*(ll-1.0)/ncomp;
			SumCov = SumCov * ll*(ll-1.0)/ncomp;
			SumCov /= ll*(ll-1.0);
			est_theta2 = ((double)ll*(double)ll-(double)ll)/(anx/(double)ll*anx/(double)ll+SumCov/(1.0*1.0));
			
			*bnxp  = Den2;
			
			Den2 *= est_theta2/(1.0*1.0);
			
			*theta_square  = est_theta2;
		}
	}
	Den += Den2;
	free(Denm12);
	free(llv);
	
	return(Num/sqrt(Den));
}
double freqtest_noutg_missing(double **ominx,long int **eix,int **nx, double thetaw, int nsam, long int length, int ncomp, struct covar_ij **Covi1i2, long int *count_cov, int singleton,double *bnxp, double *theta_square) 
{
	double Num,Num0,Den,Den2,Den12;
	long int j,k,ll,*llv,h,x,ss;
	int i,ii,i1,i2,i12,n1,n2,n12;
	struct covar_m {
		int n1;
		int n2;
		int n12;
		double Den12;
		double Cov12;
	} *Denm12;
	long int count;
	double anx,est_theta2;
	double betan,omi,omj,psi,psj;
	double SumCov,Cov12,CovMt,CovMtS;
	double alfat,betat,*ww,sumww;
	
#if DEBUG_EMANUELE	
	int *n1x,*n2y,*n12xy;
	float wv,var_d,var_h;
#endif	
    CovMtS = 0.;
    if(nsam < 2) return(-10000);

	if((llv = (long int *)calloc(length,sizeof(long int))) == 0) {
		printf("Error calloc memory freqtest0");
		exit(1);
	}
	if((Denm12 = (struct covar_m *)calloc(1,sizeof(struct covar_m))) == 0) {
		printf("Error calloc memory freqtest1");
		exit(1);
	}
	count = 0;

	Num = 0.;
	ll = 0;
	n1 = 0;
	i = 0;
	for(j=0;j<length;j++) {
		n1=0;
		for(i=0;i<nsam;i++) {n1 += nx[i][j];} /*number of samples*/
		Num0 = 0.;
		for(i=1;i<=floor((double)n1/2.);i++) {
			if(i == (n1/2.)/*n1-i*/) {
				Num0 += 2.*(double)i*((double)n1-i)/(double)n1 * ominx[i][j] * (eix[i][j] + eix[n1-i][j])/2.0;
			}
			else {
				Num0 += 1.*(double)i*((double)n1-i)/(double)n1 * ominx[i][j] * (eix[i][j] + eix[n1-i][j])/1.0;
			} 
		}
		Num += Num0;
		if(Num0) {llv[ll] = j;ll += 1;}
	}
    if(Num == 0.) {
        free(llv);
        free(Denm12);
        return(-10000);
    }

#if DEBUG_EMANUELE	
	/*calculation of Emanuele*/
	if((n1x = (int *)calloc(length,sizeof(long int))) == 0) {
		printf("Error calloc memory freqtest0");
		exit(1);
	}
	if((n2y = (int *)calloc(length,sizeof(long int))) == 0) {
		printf("Error calloc memory freqtest0");
		exit(1);
	}
	if((n12xy = (int *)calloc(length,sizeof(long int))) == 0) {
		printf("Error calloc memory freqtest0");
		exit(1);
	}
	for(j=0;j<length;j++) {
		for(i=0;i<nsam;i++) {n1x[j] += nx[i][j];}
		for(i=0;i<nsam;i++) {n2y[length-1-j] += nx[i][length-1-j];}
		for(i=0;i<nsam;i++) {n12xy[j] += nx[i][j] & nx[i][length-1-j];}
	}
	
	wv = watterson_variance(thetaw,ll,n1x,n2y,n12xy,,&var_h);

	free(n1x);
	free(n2y);
	free(n12xy);
#endif
	
	/*calculate the first term of the variance*/
	Den = 0.;
	for(j=0;j<ll;j++) {
		n1=0;
		for(i=0;i<nsam;i++) {n1 += nx[i][llv[j]];} /*number of samples*/
		for(i=1;i<=floor(n1/2.);i++) {
			if(i == (n1/2.)/*n1-i*/) {
				Den  += 2.*(double)i*(double)(n1-i)/(double)n1 * ominx[i][llv[j]] * ominx[i][llv[j]]; /*first part of the variance*/
			}
			else {
				Den  += 1.*(double)i*(double)(n1-i)/(double)n1 * ominx[i][llv[j]] * ominx[i][llv[j]]; /*first part of the variance*/
			} 
		}
	}
	Den *= thetaw/(double)ll;
	
	/*calculate the covariance: ncomp=0(fast_approach except for large S), ncomp=large(total) OR a ncomp=small(random set of comparisons)*/
	Den2 = 0.;alfat = 0.;betat = betan = 0.;
	if(ncomp==0) {
		ww = (double *)calloc(nsam,sizeof(double));
		for(j=0;j<ll;j++) {
			n1 = 0;
			for(i=0;i<nsam;i++) {n1 += nx[i][llv[j]];} /*number of samples*/
			sumww = 0.;
			for(i1=1;i1<=floor(n1/2);i1++) {
				if(i1==1) ss = singleton;
				else ss = 1;
				if(i1 == n1-i1) ww[i1] = (double)n1/((double)i1*(double)(n1 - i1)*2.0)*(double)ss;
				else ww[i1] = (double)n1/((double)i1*(double)(n1 - i1)*1.0)*(double)ss;
				sumww += ww[i1];
			}
			betan = 0.;
			alfat = 0.;
			betat = 0.;
			for(i1=1;i1<=floor(n1/2);i1++) {
				omi = ominx[i1][llv[j]];
				psi = psii(n1,i1);
				alfat += (ww[i1]/sumww * ww[i1]/sumww)/psi;
				betat += (ww[i1]/sumww)/psi * (ww[i1]/sumww)/psi * rhoii(n1,i1);
				betan += omi/psi * omi/psi * rhoii(n1,i1);
				for(i2=i1+1;i2<=floor(n1/2);i2++) {
					omj = ominx[i2][llv[j]];
					psj = psii(n1,i2);
					betat += 2.0 * (ww[i1]/sumww)/psi * (ww[i2]/sumww)/psj * rhoij(n1,i1,i2);
					betan += 2.0 * omi/psi * omj/psj * rhoij(n1,i1,i2);
				}
			}
			Den2 += betan;
		}
		*bnxp  = Den2;

		est_theta2 = (thetaw*thetaw - alfat*thetaw)/(1.0 + betat);	
		Den2 *= est_theta2/ll;
		free(ww);
		
		*theta_square  = est_theta2;
	}
	else {
		anx = 0.;
		for(j=0;j<ll;j++) {
			n1=0; 
			for(i1=0;i1<nsam;i1++) {
				n1 += (double)nx[i1][llv[j]];
			} 
			for(i1=1;i1<n1;i1++) {
				anx+= 1.0/(double)i1;
			} 
			if(singleton == 0) anx -= n1/(double)(n1-1.0);
		}
		SumCov = 0.;
		if(ncomp >= ll*(ll-1)/2.) {
			for(j=0;j<ll;j++) {
				n1=0;
				for(i1=0;i1<nsam;i1++) {n1 += nx[i1][llv[j]];} /*number of samples*/
				for(k=j+1;k<ll;k++) {
					n2=0;
					for(i2=0;i2<nsam;i2++) {n2 += nx[i2][llv[k]];} 
					n12=0;
					for(i=0;i<nsam;i++) {n12 += nx[i][llv[j]] & nx[i][llv[k]];} 
					for(i=0;i<count;i++) {
						if(((Denm12[i].n1 == n1 && Denm12[i].n2 == n2) || (Denm12[i].n1 == n2 && Denm12[i].n2 == n1)) && Denm12[i].n12 == n12) {
							Den2 += Denm12[i].Den12;
							SumCov += Denm12[i].Cov12;
							break;
						}
					}
					if(i>=count) {
						Den12 = 0.;
						Cov12 = 0.;
						for(i1=1+!singleton;i1<=(int)floor(n1/2.);i1++) {
							for(i2=1+!singleton;i2<=(int)floor(n2/2.);i2++) {
								for(x=0;x<*count_cov;x++) {
									if((Covi1i2[0][x].i1==i1 && Covi1i2[0][x].n1==n1 && Covi1i2[0][x].i2==i2 && Covi1i2[0][x].n2==n2 && Covi1i2[0][x].n12==n12) ||
									   (Covi1i2[0][x].i1==i2 && Covi1i2[0][x].n1==n2 && Covi1i2[0][x].i2==i1 && Covi1i2[0][x].n2==n1 && Covi1i2[0][x].n12==n12)) {
										CovMtS = Covi1i2[0][x].Covij;
										break;
									}
								}
								if(x>=*count_cov) {
									CovMtS = (double)cov_missing(i1,i2,n1,n2,n12,(float)1);
									Covi1i2[0][*count_cov].i1 = i1;
									Covi1i2[0][*count_cov].i2 = i2;
									Covi1i2[0][*count_cov].n1 = n1;
									Covi1i2[0][*count_cov].n2 = n2;
									Covi1i2[0][*count_cov].n12 = n12;
									Covi1i2[0][*count_cov].Covij = CovMtS;
									*count_cov += 1;
									if((Covi1i2[0] = (struct covar_ij *)realloc(Covi1i2[0],(count_cov[0]+1)*sizeof(struct covar_ij))) == 0) {
										printf("Error calloc memory freqtestl");
										exit(1);
									}
								}
								CovMt = CovMtS;
								for(x=0;x<*count_cov;x++) {
									if((Covi1i2[0][x].i1==n1-i1 && Covi1i2[0][x].n1==n1 && Covi1i2[0][x].i2==i2 && Covi1i2[0][x].n2==n2 && Covi1i2[0][x].n12==n12) ||
									   (Covi1i2[0][x].i1==i2 && Covi1i2[0][x].n1==n2 && Covi1i2[0][x].i2==n1-i1 && Covi1i2[0][x].n2==n1 && Covi1i2[0][x].n12==n12)) {
										CovMtS = Covi1i2[0][x].Covij;
										break;
									}
								}
								if(x>=*count_cov) {
									CovMtS = (double)cov_missing(n1-i1,i2,n1,n2,n12,(float)1);
									Covi1i2[0][*count_cov].i1 = n1-i1;
									Covi1i2[0][*count_cov].i2 = i2;
									Covi1i2[0][*count_cov].n1 = n1;
									Covi1i2[0][*count_cov].n2 = n2;
									Covi1i2[0][*count_cov].n12 = n12;
									Covi1i2[0][*count_cov].Covij = CovMtS;
									*count_cov += 1;
									if((Covi1i2[0] = (struct covar_ij *)realloc(Covi1i2[0],(count_cov[0]+1)*sizeof(struct covar_ij))) == 0) {
										printf("Error calloc memory freqtestl");
										exit(1);
									}
								}
								CovMt += CovMtS;
								for(x=0;x<*count_cov;x++) {
									if((Covi1i2[0][x].i1==i1 && Covi1i2[0][x].n1==n1 && Covi1i2[0][x].i2==n2-i2 && Covi1i2[0][x].n2==n2 && Covi1i2[0][x].n12==n12) ||
									   (Covi1i2[0][x].i1==n2-i2 && Covi1i2[0][x].n1==n2 && Covi1i2[0][x].i2==i1 && Covi1i2[0][x].n2==n1 && Covi1i2[0][x].n12==n12)) {
										CovMtS = Covi1i2[0][x].Covij;
										break;
									}
								}
								if(x>=*count_cov) {
									CovMtS = (double)cov_missing(i1,n2-i2,n1,n2,n12,(float)1);
									Covi1i2[0][*count_cov].i1 = i1;
									Covi1i2[0][*count_cov].i2 = n2-i2;
									Covi1i2[0][*count_cov].n1 = n1;
									Covi1i2[0][*count_cov].n2 = n2;
									Covi1i2[0][*count_cov].n12 = n12;
									Covi1i2[0][*count_cov].Covij = CovMtS;
									*count_cov += 1;
									if((Covi1i2[0] = (struct covar_ij *)realloc(Covi1i2[0],(count_cov[0]+1)*sizeof(struct covar_ij))) == 0) {
										printf("Error calloc memory freqtestl");
										exit(1);
									}
								}
								CovMt += CovMtS;
								for(x=0;x<*count_cov;x++) {
									if((Covi1i2[0][x].i1==n1-i1 && Covi1i2[0][x].n1==n1 && Covi1i2[0][x].i2==n2-i2 && Covi1i2[0][x].n2==n2 && Covi1i2[0][x].n12==n12) ||
									   (Covi1i2[0][x].i1==n2-i2 && Covi1i2[0][x].n1==n2 && Covi1i2[0][x].i2==n1-i1 && Covi1i2[0][x].n2==n1 && Covi1i2[0][x].n12==n12)) {
										CovMtS = Covi1i2[0][x].Covij;
										break;
									}
								}
								if(x>=*count_cov) {
									CovMtS = (double)cov_missing(n1-i1,n2-i2,n1,n2,n12,(float)1);
									Covi1i2[0][*count_cov].i1 = n1-i1;
									Covi1i2[0][*count_cov].i2 = n2-i2;
									Covi1i2[0][*count_cov].n1 = n1;
									Covi1i2[0][*count_cov].n2 = n2;
									Covi1i2[0][*count_cov].n12 = n12;
									Covi1i2[0][*count_cov].Covij = CovMtS;
									*count_cov += 1;
									if((Covi1i2[0] = (struct covar_ij *)realloc(Covi1i2[0],(count_cov[0]+1)*sizeof(struct covar_ij))) == 0) {
										printf("Error calloc memory freqtestl");
										exit(1);
									}
								}
								CovMt += CovMtS;
								/*
								 CovMt = ((double)cov_missing(i1,i2,n1,n2,n12,(float)1) +
								 (double)cov_missing(n1-i1,i2,n1,n2,n12,(float)1) + 
								 (double)cov_missing(i1,n2-i2,n1,n2,n12,(float)1) +
								 (double)cov_missing(n1-i1,n2-i2,n1,n2,n12,(float)1));
								 */
								if(i1 != (n1/2.) && i2 != (n2/2.)/*i1 != n1-i1 && i2 != n2-i2*/) {
									CovMt /= 1.0;
									Cov12 += CovMt;
									Den12 += 1.*(double)i1*((double)n1-i1)/(double)n1 * 
									1.*(double)i2*((double)n2-i2)/(double)n2 * ominx[i1][llv[j]] * ominx[i2][llv[k]] * CovMt;
								}
								else {
									if(i1 == (n1/2.) && i2 == (n2/2.)/*i1 == n1-i1 && i2 == n2-i2*/) {
										CovMt /= 4.0;
										Cov12 += CovMt;
										Den12 += 2.*(double)i1*((double)n1-i1)/(double)n1 * 
										2.*(double)i2*((double)n2-i2)/(double)n2 * ominx[i1][llv[j]] * ominx[i2][llv[k]] * CovMt;
									}
									else {
										CovMt /= 2.0;
										Cov12 += CovMt;
										Den12 += 2.*(double)i1*((double)n1-i1)/(double)n1 * 
										1.*(double)i2*((double)n2-i2)/(double)n2 * ominx[i1][llv[j]] * ominx[i2][llv[k]] * CovMt;
									}
								}
							}
						}
						Den12 /=(double)(ll-1.0)*ll;
						Denm12[count].n1 = n1;
						Denm12[count].n2 = n2;
						Denm12[count].n12 = n12;
						Denm12[count].Den12 = Den12 /*include some factor here*/;
						Denm12[count].Cov12 = Cov12 /*include some factor here*/;
						count += 1;
						if((Denm12 = (struct covar_m *)realloc(Denm12,(count+1)*sizeof(struct covar_m))) == 0) {
							printf("Error calloc memory freqtest2");
							exit(1);
						}
						
						Den2 += Den12;
						SumCov += Cov12;
					}
				}
			}
			Den2 *= 2.;
			SumCov *= 2.;
			SumCov /= ll*(ll-1.0);
			est_theta2 = ((double)ll*(double)ll-(double)ll)/(anx/(double)ll*anx/(double)ll+SumCov/(1.0*1.0));

			*bnxp  = Den2;

			Den2 *= est_theta2/(1.0*1.0);
			
			*theta_square  = est_theta2;
		}
		else {
			ii=0;
			while(ii<ncomp) {
				j=0;k=0;
				while(j==k) {
					j = (long int)floor(ran1()*ll);
					k = (long int)floor(ran1()*ll);
				}
				if(j>k) {h=k;k=j;j=h;}
				n1=0;
				for(i1=0;i1<nsam;i1++) {n1 += nx[i1][llv[j]];}
				n2=0;
				for(i2=0;i2<nsam;i2++) {n2 += nx[i2][llv[k]];} 
				n12=0;
				for(i12=0;i12<nsam;i12++) {n12 += nx[i12][llv[j]] & nx[i12][llv[k]];} 
				for(i=0;i<count;i++) {
					if(((Denm12[i].n1 == n1 && Denm12[i].n2 == n2) || (Denm12[i].n1 == n2 && Denm12[i].n2 == n1)) && Denm12[i].n12 == n12) {
						Den2 += Denm12[i].Den12;
						SumCov += Denm12[i].Cov12;
						break;
					}
				}
				if(i>=count) {
					Den12 = 0.;
					Cov12 = 0.;
					for(i1=1+!singleton;i1<=(int)floor(n1/2.);i1++) {
						for(i2=1+!singleton;i2<=(int)floor(n2/2.);i2++) {
							for(x=0;x<*count_cov;x++) {
								if((Covi1i2[0][x].i1==i1 && Covi1i2[0][x].n1==n1 && Covi1i2[0][x].i2==i2 && Covi1i2[0][x].n2==n2 && Covi1i2[0][x].n12==n12) ||
								   (Covi1i2[0][x].i1==i2 && Covi1i2[0][x].n1==n2 && Covi1i2[0][x].i2==i1 && Covi1i2[0][x].n2==n1 && Covi1i2[0][x].n12==n12)) {
									CovMtS = Covi1i2[0][x].Covij;
									break;
								}
							}
							if(x>=*count_cov) {
								CovMtS = (double)cov_missing(i1,i2,n1,n2,n12,(float)1);
								Covi1i2[0][*count_cov].i1 = i1;
								Covi1i2[0][*count_cov].i2 = i2;
								Covi1i2[0][*count_cov].n1 = n1;
								Covi1i2[0][*count_cov].n2 = n2;
								Covi1i2[0][*count_cov].n12 = n12;
								Covi1i2[0][*count_cov].Covij = CovMtS;
								*count_cov += 1;
								if((Covi1i2[0] = (struct covar_ij *)realloc(Covi1i2[0],(count_cov[0]+1)*sizeof(struct covar_ij))) == 0) {
									printf("Error calloc memory freqtestl");
									exit(1);
								}
							}
							CovMt = CovMtS;
							for(x=0;x<*count_cov;x++) {
								if((Covi1i2[0][x].i1==n1-i1 && Covi1i2[0][x].n1==n1 && Covi1i2[0][x].i2==i2 && Covi1i2[0][x].n2==n2 && Covi1i2[0][x].n12==n12) ||
								   (Covi1i2[0][x].i1==i2 && Covi1i2[0][x].n1==n2 && Covi1i2[0][x].i2==n1-i1 && Covi1i2[0][x].n2==n1 && Covi1i2[0][x].n12==n12)) {
									CovMtS = Covi1i2[0][x].Covij;
									break;
								}
							}
							if(x>=*count_cov) {
								CovMtS = (double)cov_missing(n1-i1,i2,n1,n2,n12,(float)1);
								Covi1i2[0][*count_cov].i1 = n1-i1;
								Covi1i2[0][*count_cov].i2 = i2;
								Covi1i2[0][*count_cov].n1 = n1;
								Covi1i2[0][*count_cov].n2 = n2;
								Covi1i2[0][*count_cov].n12 = n12;
								Covi1i2[0][*count_cov].Covij = CovMtS;
								*count_cov += 1;
								if((Covi1i2[0] = (struct covar_ij *)realloc(Covi1i2[0],(count_cov[0]+1)*sizeof(struct covar_ij))) == 0) {
									printf("Error calloc memory freqtestl");
									exit(1);
								}
							}
							CovMt += CovMtS;
							for(x=0;x<*count_cov;x++) {
								if((Covi1i2[0][x].i1==i1 && Covi1i2[0][x].n1==n1 && Covi1i2[0][x].i2==n2-i2 && Covi1i2[0][x].n2==n2 && Covi1i2[0][x].n12==n12) ||
								   (Covi1i2[0][x].i1==n2-i2 && Covi1i2[0][x].n1==n2 && Covi1i2[0][x].i2==i1 && Covi1i2[0][x].n2==n1 && Covi1i2[0][x].n12==n12)) {
									CovMtS = Covi1i2[0][x].Covij;
									break;
								}
							}
							if(x>=*count_cov) {
								CovMtS = (double)cov_missing(i1,n2-i2,n1,n2,n12,(float)1);
								Covi1i2[0][*count_cov].i1 = i1;
								Covi1i2[0][*count_cov].i2 = n2-i2;
								Covi1i2[0][*count_cov].n1 = n1;
								Covi1i2[0][*count_cov].n2 = n2;
								Covi1i2[0][*count_cov].n12 = n12;
								Covi1i2[0][*count_cov].Covij = CovMtS;
								*count_cov += 1;
								if((Covi1i2[0] = (struct covar_ij *)realloc(Covi1i2[0],(count_cov[0]+1)*sizeof(struct covar_ij))) == 0) {
									printf("Error calloc memory freqtestl");
									exit(1);
								}
							}
							CovMt += CovMtS;
							for(x=0;x<*count_cov;x++) {
								if((Covi1i2[0][x].i1==n1-i1 && Covi1i2[0][x].n1==n1 && Covi1i2[0][x].i2==n2-i2 && Covi1i2[0][x].n2==n2 && Covi1i2[0][x].n12==n12) ||
								   (Covi1i2[0][x].i1==n2-i2 && Covi1i2[0][x].n1==n2 && Covi1i2[0][x].i2==n1-i1 && Covi1i2[0][x].n2==n1 && Covi1i2[0][x].n12==n12)) {
									CovMtS = Covi1i2[0][x].Covij;
									break;
								}
							}
							if(x>=*count_cov) {
								CovMtS = (double)cov_missing(n1-i1,n2-i2,n1,n2,n12,(float)1);
								Covi1i2[0][*count_cov].i1 = n1-i1;
								Covi1i2[0][*count_cov].i2 = n2-i2;
								Covi1i2[0][*count_cov].n1 = n1;
								Covi1i2[0][*count_cov].n2 = n2;
								Covi1i2[0][*count_cov].n12 = n12;
								Covi1i2[0][*count_cov].Covij = CovMtS;
								*count_cov += 1;
								if((Covi1i2[0] = (struct covar_ij *)realloc(Covi1i2[0],(count_cov[0]+1)*sizeof(struct covar_ij))) == 0) {
									printf("Error calloc memory freqtestl");
									exit(1);
								}
							}
							CovMt += CovMtS;
							/*
							 CovMt = ((double)cov_missing(i1,i2,n1,n2,n12,(float)1) +
							 (double)cov_missing(n1-i1,i2,n1,n2,n12,(float)1) + 
							 (double)cov_missing(i1,n2-i2,n1,n2,n12,(float)1) +
							 (double)cov_missing(n1-i1,n2-i2,n1,n2,n12,(float)1));
							 */
							if(i1 != (n1/2.) && i2 != (n2/2.)/*i1 != n1-i1 && i2 != n2-i2*/) {
								CovMt /= 1.0;
								Cov12 += CovMt;
								Den12 += 1.*(double)i1*((double)n1-i1)/(double)n1 * 
								1.*(double)i2*((double)n2-i2)/(double)n2 * ominx[i1][llv[j]] * ominx[i2][llv[k]] * CovMt;
							}
							else {
								if(i1 == (n1/2.) && i2 == (n2/2.)/*i1 == n1-i1 && i2 == n2-i2*/) {
									CovMt /= 4.0;
									Cov12 += CovMt;
									Den12 += 2.*(double)i1*((double)n1-i1)/(double)n1 * 
									2.*(double)i2*((double)n2-i2)/(double)n2 * ominx[i1][llv[j]] * ominx[i2][llv[k]] * CovMt;
								}
								else {
									CovMt /= 2.0;
									Cov12 += CovMt;
									Den12 += 2.*(double)i1*((double)n1-i1)/(double)n1 * 
									1.*(double)i2*((double)n2-i2)/(double)n2 * ominx[i1][llv[j]] * ominx[i2][llv[k]] * CovMt;
								}
							}
						}
					}
					Den12 /=(double)(ll-1.0)*ll;
					Denm12[count].n1 = n1;
					Denm12[count].n2 = n2;
					Denm12[count].n12 = n12;
					Denm12[count].Den12 = Den12 /*include some factor here*/;
					Denm12[count].Cov12 = Cov12;
					count += 1;
					if((Denm12 = (struct covar_m *)realloc(Denm12,(count+1)*sizeof(struct covar_m))) == 0) {
						printf("Error calloc memory freqtest2");
						exit(1);
					}
					
					Den2 += Den12;
					SumCov += Cov12;
				}
				ii += 1;
			}
			Den2 = Den2 * ll*(ll-1.0)/ncomp;
			SumCov = SumCov * ll*(ll-1.0)/ncomp;
			SumCov /= ll*(ll-1.0);
			est_theta2 = ((double)ll*(double)ll-(double)ll)/(anx/(double)ll*anx/(double)ll+SumCov/(1.0*1.0));
			
			*bnxp  = Den2;
			
			Den2 *= est_theta2/(1.0*1.0);
			
			*theta_square  = est_theta2;
		}
	}
	Den += Den2;
	free(Denm12);
	free(llv);
	
    if(Den==0.) return(-10000);
	return(Num/sqrt(Den));
}
long int maxd(double x,double y) {
    if(x > y) return x;
    else return y;
}

