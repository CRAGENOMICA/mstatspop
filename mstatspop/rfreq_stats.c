//
//  rfreq_stats.c
//  mstatspop
//
//  Created by Sebastian Ramos-Onsins on 11/6/25.
//

#include "rfreq_stats.h"

/*
 * @param npops
 * @param nsam
 * @param matrix_pol
 * @param length
 * @param statistics
 * @param outgroup_presence
 * @param force_outgroup
 * @param include_unknown
 * @return
 */

/*What it does:
 - pre-calculate the weights
 - first calculates the frequency and nsam for each position considering the total population samples. length -> nsam/pos -> freq/pos
 - second counts for each pop (individual) if there is a mutation (derived if outgroup, or the major/minor frequency if no-outgroup) and assign a value=total frequency from the whole population. length -> subnsam -> S? -> freq/pos
 - third calculates the pi for each subgroup (individual). heterozigosity -> theta
 - fourth calculates the relative variability for the different theta statatistics. theta
 - fifth calculates non-normalized neutrality tests. test.
 - Optimize
 */
/* this function works instead the calculation of SFS */
int calc_rfreqstats(int npops, int *nsam, char *matrix_pol, long int length, struct stats *statistics, struct rweights *rw, int outgroup_presence, int force_outgroup, int include_unknown)
{
    int x,y,z,pop1,fr,frn,frc;//,fr0,f0;
    long int j;
    double pia;
    
    //freqs per position for total samples (freq) and for the outgroup (freqo)
    int freq[4],freqo[4];
    //freqs per position/population and relative frequency in relation to total
    int **freqp,**freqr;
    
    //int *no;//valid outgroup or not per position
    
    char ancestral[1]; //determine if there is ancestral
    
    int *initsq1;//sample location per population
    int inits,max;//initial and final sample location per pop
    int sumnsam;//,maxnsam;//sum all samples and larger sample pop
    
    /*weights*/
    /*
    double **ww;
    double **wt;
    double **wfw;
    double **wfl;
    double **wl;
    double **wwA;
    double **wtA;
    double **wfwA;
    double **wflA;
    double **wlA;
    double **wphi_i;
    double ***wpsi_ij;
    */
    
    /*init*/
    for(pop1=0;pop1<npops-1;pop1++) {
        statistics[0].piw[pop1]     = 0;
        /* we will use this definitions to include all rSFS statistics */
        statistics[0].S[pop1]       = (double)0;
        statistics[0].rthetaS[pop1]  = (double)0;
        statistics[0].rthetaT[pop1]  = (double)0;
        statistics[0].So[pop1]      = (double)0;
        statistics[0].rthetaSo[pop1] = (double)0;
        statistics[0].rthetaTo[pop1] = (double)0;
        statistics[0].rthetaFL[pop1] = (double)0;
        statistics[0].rthetaFW[pop1] = (double)0;
        statistics[0].rthetaL[pop1]  = (double)0;
        statistics[0].rthetaSA[pop1] = (double)0;
        statistics[0].rthetaTA[pop1] = (double)0;
        
        statistics[0].rDtaj[pop1]    = (double)0;//rthetaT - rthetaS;
        statistics[0].rDfl[pop1]     = (double)0;//rthetaS - rthetaFL;
        statistics[0].rFfl[pop1]     = (double)0;//rthetaT - rthetaFL;
        statistics[0].rHnfw[pop1]    = (double)0;//rthetaT - rthetaFW;
        statistics[0].rEz[pop1]      = (double)0;//rthetaL - rthetaS;
        statistics[0].rYach[pop1]    = -10000.0;
        statistics[0].rFH[pop1]      = -10000.0;
        
        statistics[0].K[pop1]       = (double)0;
        
        for(x=1;x<nsam[pop1];x++) {
            statistics[0].freq[pop1][x] = 0;
        }
    }
    
    if(length==0) return 1;
    
    ancestral[0] = 0;
    initsq1      = (int *)calloc(npops,sizeof(int));
    
    sumnsam = 0;
    //maxnsam = 0;
    for(x=0;x<npops;x++) {
        initsq1[x] = sumnsam;
        sumnsam += nsam[x];
        //if(maxnsam < nsam[x]) maxnsam = nsam[x];
    }
    
    //for each position, count the fequency (and sample size)
    //no = (int *)calloc(length,sizeof(int));/*if outgroup*/
    freqp   = (int **)calloc(npops,sizeof(int *));/*4 values for each pop*/
    freqr   = (int **)calloc(npops,sizeof(int *));/*4 values for each pop*/
    for(x=0;x<npops;x++) {
        freqp[x]   = (int *)calloc(4,sizeof(int));
        freqr[x]   = (int *)calloc(4,sizeof(int));
    }
    
    /*precalculation of weights:
     for each possible nsam define an array of size (1..nsam-1, here nsam+1 to include 0).
     For hypergeometrical, include the subsample
    
    ww   = (double **) calloc(sumnsam,sizeof(double *));
    wt   = (double **) calloc(sumnsam,sizeof(double *));
    wfl  = (double **) calloc(sumnsam,sizeof(double *));
    wfw  = (double **) calloc(sumnsam,sizeof(double *));
    wl   = (double **) calloc(sumnsam,sizeof(double *));
    wwA  = (double **) calloc(sumnsam,sizeof(double *));
    wtA  = (double **) calloc(sumnsam,sizeof(double *));
    wflA = (double **) calloc(sumnsam,sizeof(double *));
    wfwA = (double **) calloc(sumnsam,sizeof(double *));
    wlA  = (double **) calloc(sumnsam,sizeof(double *));
    
    wphi_i  = (double **) calloc(sumnsam,sizeof(double *));
    wpsi_ij = (double ***) calloc(sumnsam,sizeof(double **));
    
    for(x=0;x<sumnsam;x++) {
        ww[x]   = (double *) calloc(x+1,sizeof(double));
        wt[x]   = (double *) calloc(x+1,sizeof(double));
        wfl[x]  = (double *) calloc(x+1,sizeof(double));
        wfw[x]  = (double *) calloc(x+1,sizeof(double));
        wl[x]   = (double *) calloc(x+1,sizeof(double));
        wwA[x]  = (double *) calloc(x+1,sizeof(double));
        wtA[x]  = (double *) calloc(x+1,sizeof(double));
        wflA[x] = (double *) calloc(x+1,sizeof(double));
        wfwA[x] = (double *) calloc(x+1,sizeof(double));
        wlA[x]  = (double *) calloc(x+1,sizeof(double));
        wphi_i[x]  = (double *) calloc(x+1,sizeof(double));
        
        if(x>1) {
            if(outgroup_presence+force_outgroup==1) {
                //alculate the weights for unfolded
                weights_unfolded(ww[x], wt[x],  wfw[x], wfl[x], wl[x], x, 0.0);
                //calculate the weights for unfolded excluding singletons
                weights_unfolded(wwA[x], wtA[x], wfwA[x], wflA[x], wlA[x], x, 1.0/(double)x);
            }
            else {
                //calculate the weights for folded
                weights_folded(ww[x], wt[x],  wfl[x], wl[x], wphi_i[x], x, 0.0);
                //calculate the weights for folded excluding singletons
                weights_folded(wwA[x], wtA[x],  wflA[x],wlA[x], wphi_i[x], x, 1.0/(double)x);
            }
        }
        
        wpsi_ij[x] = (double **) calloc(x+1,sizeof(double *));
        for(z=0;z<x+1;z++) {
            wpsi_ij[x][z] = (double *) calloc(maxnsam+1,sizeof(double));
            if(x>1) {
                if(outgroup_presence+force_outgroup==1) weights_unfolded_wpsi(wpsi_ij[x][z], x, z, maxnsam);
                else weights_folded_wpsi(wpsi_ij[x][z], x, z, maxnsam);
            }
        }
    }
    */
    //ww[5][2] means the watt weight for n=5 and freq=2
    //wpsi_ij[5][2][2] is the hyperg. prob. for n=5,freq=2 and subsample=2
    
    //calculate the freqs
    //f0=fr0=0;
    for(j=0;j<length;j++) {
        /*calculate the relative frequencies of each population for each position*/
        if(outgroup_presence+force_outgroup==1)
        {
            pia=0.;
            freqo[0]=freqo[1]=freqo[2]=freqo[3]=0;
            /*outgroup*/
            for(y=initsq1[npops-1];y<sumnsam;y++) {
                if(matrix_pol[j*sumnsam+y] == '0') {freqo[1] += 1;freqo[0] += 1;}
                if(matrix_pol[j*sumnsam+y] == '1') {freqo[2] += 1;freqo[0] += 1;}
                if(matrix_pol[j*sumnsam+y] == '-') {freqo[3] += 1;}
            }
            if(freqo[0]) {
                if(freqo[1] != freqo[0] && freqo[1] != 0) {
                    ancestral[0] = (char)0;/*if the outgroup is polymorphic, we do not consider for neutrality tests with outgroup*/
                }
                else {
                    if(freqo[1] == freqo[0]) {
                        ancestral[0] = '0';
                    }
                    if(freqo[2] == freqo[0]) {
                        ancestral[0] = '1';
                    }
                }
            }
            else
                ancestral[0] = (char)0;
            //no[j] = ancestral[0]; /*if the outgroup is valid 1, if not 0*/
            
            /* calculate freqs for each pop*/
            //init
            for(pop1=0;pop1<npops-1;pop1++) {
                freqp[pop1][0]=freqp[pop1][1]=freqp[pop1][2]=freqp[pop1][3]=0;
            }
            
            //start counting freqs for entire pop and for each individually
            pop1 = 0;
            inits = initsq1[pop1];
            max   = initsq1[pop1]+nsam[pop1];
            freq[0]=freq[1]=freq[2]=freq[3]=0; //freq for the whole sample
            
            for(y=0;y<sumnsam-nsam[npops-1];y++) {
                if(y==max) {/*additionally calculate the freq per pop*/
                    if(ancestral[0] == '1') {/*count freq once the pop freq is calculated (considering outgroup) and revert 0->1 and 1->0 if necessary*/
                        fr = freqp[pop1][1];
                        freqp[pop1][1] = freqp[pop1][2];
                        freqp[pop1][2] = fr;
                    }
                    pop1 += 1;
                    inits = initsq1[pop1];
                    max   = initsq1[pop1]+nsam[pop1];
                }
                if(matrix_pol[j*sumnsam+y] == '0') {
                    freq[1]  += 1;freq[0]  += 1;
                    freqp[pop1][1] += 1;freqp[pop1][0] += 1;
                }
                if(matrix_pol[j*sumnsam+y] == '1') {
                    freq[2]  += 1;freq[0]  += 1;
                    freqp[pop1][2] += 1;freqp[pop1][0] += 1;
                }
                if(matrix_pol[j*sumnsam+y] == '-') {
                    freq[3]  += 1;
                    freqp[pop1][3] += 1;
                }
            }
            if(ancestral[0] == '1') {
                fr = freq[1];
                freq[1] = freq[2];
                freq[2] = fr;
            }
            //calculate statistics
            if(freqo[0] && freq[0]) {
                statistics[0].K[0] += (double)(freq[2]*freqo[1] + freq[1]*freqo[2])/(double)(freq[0]*freqo[0]);
            }
            for(pop1=0;pop1<npops;pop1++) {
                /*calculate standard pi and S per population*/
                statistics[0].piw[pop1] += Calc_Theta_unfolded(freqp[pop1][2], rw->wt[freqp[pop1][0]][freqp[pop1][2]],
                                                               freqp[pop1][0], rw->wt[freqp[pop1][0]][0]);
                if(freqp[pop1][2]>0 && freqp[pop1][2]<freqp[pop1][0]) {//variant no fixed, count freq
                    statistics[0].S[pop1] += 1.0; //no consider outgroup
                    if(ancestral[0] != (char)0) {//consider outgroup
                        statistics[0].So[pop1]  += 1.0;
                     }
                }
                
                /*calculate relative variability*/
                fr=frn=0;
                if(freqp[pop1][2]>0 && freq[2]<freq[0]) {//variant exist in pop1 and is no fixed in global
                    //count global freq
                    frn = freq[2];
                    if(ancestral[0])
                        fr = freq[2];
                    
                    statistics[0].rS[pop1] += (frn>0);
                    statistics[0].rthetaS[pop1] += Calc_rTheta_unfolded(frn,rw->ww[freq[0]][freq[2]],
                                                                        rw->wpsi_ij[freq[0]][freq[2]][freqp[pop1][0]],
                                                                        freq[0], rw->ww[freq[0]][0]);
                    statistics[0].rthetaT[pop1] += Calc_rTheta_unfolded(frn,rw->wt[freq[0]][freq[2]],
                                                                        rw->wpsi_ij[freq[0]][freq[2]][freqp[pop1][0]],
                                                                        freq[0], rw->wt[freq[0]][0]);
                    if(ancestral[0]) {
                        //if(pop1==0)printf("%d\t%d\n",pop1,fr);
                        statistics[0].rSo[pop1] += (fr>0);
                        statistics[0].rthetaSo[pop1] += Calc_rTheta_unfolded(fr,rw->ww[freq[0]][freq[2]],
                                                                             rw->wpsi_ij[freq[0]][freq[2]][freqp[pop1][0]],
                                                                             freq[0], rw->ww[freq[0]][0]);
                        statistics[0].rthetaTo[pop1] += Calc_rTheta_unfolded(fr,rw->wt[freq[0]][freq[2]],
                                                                             rw->wpsi_ij[freq[0]][freq[2]][freqp[pop1][0]],
                                                                             freq[0], rw->wt[freq[0]][0]);
                        statistics[0].rthetaFL[pop1]  += Calc_rTheta_unfolded(fr,rw->wfl[freq[0]][freq[2]],
                                                                              rw->wpsi_ij[freq[0]][freq[2]][freqp[pop1][0]],
                                                                              freq[0], rw->wfl[freq[0]][0]);
                        statistics[0].rthetaFW[pop1] += Calc_rTheta_unfolded(fr,rw->wfw[freq[0]][freq[2]],
                                                                             rw->wpsi_ij[freq[0]][freq[2]][freqp[pop1][0]],
                                                                             freq[0], rw->wfw[freq[0]][0]);
                        statistics[0].rthetaL[pop1] += Calc_rTheta_unfolded(fr,rw->wl[freq[0]][freq[2]],
                                                                            rw->wpsi_ij[freq[0]][freq[2]][freqp[pop1][0]],
                                                                            freq[0], rw->wl[freq[0]][0]);
                        statistics[0].rthetaSA[pop1] += Calc_rTheta_unfolded(fr,rw->wwA[freq[0]][freq[2]],
                                                                             rw->wpsi_ij[freq[0]][freq[2]][freqp[pop1][0]],
                                                                             freq[0], rw->wwA[freq[0]][0]);
                        statistics[0].rthetaTA[pop1] += Calc_rTheta_unfolded(fr,rw->wtA[freq[0]][freq[2]],
                                                                             rw->wpsi_ij[freq[0]][freq[2]][freqp[pop1][0]],
                                                                             freq[0], rw->wtA[freq[0]][0]);
                    }
                }
            }
        }
        if(outgroup_presence+force_outgroup==0) {
            //init
            for(pop1=0;pop1<npops-1;pop1++) {
                freqp[pop1][0]=freqp[pop1][1]=freqp[pop1][2]=freqp[pop1][3]=0;
            }
            
            //start counting freqs for entire pop and for each individually
            pop1 = 0;
            inits = initsq1[pop1];
            max   = initsq1[pop1]+nsam[pop1];
            freq[0]=freq[1]=freq[2]=freq[3]=0;
            
            for(y=0;y<sumnsam-nsam[npops-1];y++) {
                if(y==max) {/*additionally calculate the freq per pop. Not cahnge to minimum freq*/
                    pop1 += 1;
                    inits = initsq1[pop1];
                    max   = initsq1[pop1]+nsam[pop1];
                }
                if(matrix_pol[j*sumnsam+y] == '0') {
                    freq[1]  += 1;freq[0]  += 1;
                    freqp[pop1][1] += 1;freqp[pop1][0] += 1;
                }
                if(matrix_pol[j*sumnsam+y] == '1') {
                    freq[2]  += 1;freq[0]  += 1;
                    freqp[pop1][2] += 1;freqp[pop1][0] += 1;
                }
                if(matrix_pol[j*sumnsam+y] == '-') {
                    freq[3]  += 1;
                    freqp[pop1][3] += 1;
                }
            }
            //minimum freq on freq[2]
            if(freq[2]>freq[1]) {
                frc = freq[1];
                freq[1] = freq[2];
                freq[2] = frc;
                for(pop1=0;pop1<npops-1;pop1++) {
                    frc = freqp[pop1][1];
                    freqp[pop1][1] = freqp[pop1][2];
                    freqp[pop1][2] = frc;
                }
                //fr0+=1;
            }
            //f0+=1;
            //calculate statistics
            for(pop1=0;pop1<npops-1;pop1++) {
                /*calculate standard pi and S per population*/
                if(freqp[pop1][2]>0 && freqp[pop1][2]<freqp[pop1][0]) {//variant no fixed, count freq
                    statistics[0].piw[pop1] += Calc_Theta_folded(freqp[pop1][2], rw->wt[freqp[pop1][0]][freqp[pop1][2]],
                                                             rw->wphi_i[freqp[pop1][0]][freqp[pop1][2]], freqp[pop1][0],
                                                             rw->wt[freqp[pop1][0]][0]);
                    statistics[0].S[pop1] += 1.0; //no  outgroup
                }
                /*calculate relative variability*/
                fr=frn=0;
                if(freqp[pop1][2]>0 && freq[2]<freq[0]) {//variant exist in pop1 and is no fixed in global
                    //count global freq
                    frn = freq[2];
                    
                    statistics[0].rS[pop1] += (frn>0);
                    statistics[0].rthetaS[pop1] += Calc_rTheta_folded(frn,rw->ww[freq[0]][freq[2]],
                                                                      rw->wpsi_ij[freq[0]][freq[2]][freqp[pop1][0]],
                                                                      rw->wphi_i[freq[0]][freq[2]], freq[0],
                                                                      rw->ww[freq[0]][0]);
                    statistics[0].rthetaT[pop1] += Calc_rTheta_folded(frn,rw->wt[freq[0]][freq[2]],
                                                                      rw->wpsi_ij[freq[0]][freq[2]][freqp[pop1][0]],
                                                                      rw->wphi_i[freq[0]][freq[2]], freq[0],
                                                                      rw->wt[freq[0]][0]);
                    statistics[0].rthetaFL[pop1]  += Calc_rTheta_folded(frn,rw->wfl[freq[0]][freq[2]],
                                                                        rw->wpsi_ij[freq[0]][freq[2]][freqp[pop1][0]],
                                                                        rw->wphi_i[freq[0]][freq[2]], freq[0],
                                                                        rw->wfl[freq[0]][0]);
                    statistics[0].rthetaSA[pop1] += Calc_rTheta_folded(frn,rw->wwA[freq[0]][freq[2]],
                                                                       rw->wpsi_ij[freq[0]][freq[2]][freqp[pop1][0]],
                                                                       rw->wphi_i[freq[0]][freq[2]], freq[0],
                                                                       rw->wwA[freq[0]][0]);
                    statistics[0].rthetaTA[pop1] += Calc_rTheta_folded(frn,rw->wtA[freq[0]][freq[2]],
                                                                       rw->wpsi_ij[freq[0]][freq[2]][freqp[pop1][0]],
                                                                       rw->wphi_i[freq[0]][freq[2]], freq[0],
                                                                       rw->wtA[freq[0]][0]);
                }
            }
        }
    }
    /*Calculate the UNNORMALIZED TESTS: differences between thetas*/
    for(pop1=0;pop1<npops;pop1++) {
        statistics[0].rDtaj[pop1]    = statistics[0].rthetaT[pop1] - statistics[0].rthetaS[pop1];//rthetaT - rthetaS;
        statistics[0].rDfl[pop1]     = statistics[0].rthetaS[pop1] - statistics[0].rthetaFL[pop1];//rthetaS - rthetaFL;
        statistics[0].rFfl[pop1]     = statistics[0].rthetaT[pop1] - statistics[0].rthetaFL[pop1];//rthetaT - rthetaFL;
        statistics[0].rYach[pop1]    = statistics[0].rthetaTA[pop1] - statistics[0].rthetaSA[pop1];
        if(outgroup_presence+force_outgroup) {
            statistics[0].rHnfw[pop1]    = statistics[0].rthetaT[pop1] - statistics[0].rthetaFW[pop1];//rthetaT - rthetaFW;
            statistics[0].rEz[pop1]      = statistics[0].rthetaL[pop1] - statistics[0].rthetaS[pop1];//rthetaL - rthetaS;
            statistics[0].rFH[pop1]      = statistics[0].rthetaS[pop1] - statistics[0].rthetaFW[pop1];//rthetaS - rthetaFW
        }
    }
    
    /*free pointers: NOT optimal for small scaffolds*/
    for(x=0;x<npops;x++) {
        free(freqp[x]);
        free(freqr[x]);
    }
    free(freqp);
    free(freqr);
    free(initsq1);
    
    /*
    for(x=0;x<sumnsam;x++) {
        free(ww[x]);
        free(wt[x]);
        free(wfw[x]);
        free(wfl[x]);
        free(wl[x]);
        free(wwA[x]);
        free(wtA[x]);
        free(wfwA[x]);
        free(wflA[x]);
        free(wlA[x]);
        free(wphi_i[x]);
        for(z=0;z<x+1;z++)
            free(wpsi_ij[x][z]);
        free(wpsi_ij[x]);
    }
    free(ww);
    free(wt);
    free(wfw);
    free(wfl);
    free(wl);
    free(wwA);
    free(wtA);
    free(wfwA);
    free(wflA);
    free(wlA);
    free(wphi_i);
    free(wpsi_ij);
    */
    //printf("f0=%d fr0=%d",f0,fr0);
    return 1;
}

/* UNFOLDED */
int weights_unfolded(double *ww, double *wt, double *wfw, double *wfl, double *wl,/*double **wpsi_ij,*/ int nsam, /*int subnsam, */double freq_cut) {
    int i;
    for(i=1;i<nsam;i++) {
        //wpsi_ij[i] = (1.0-gsl_ran_hypergeometric_pdf(0, i, nsam-i, subnsam));
        //w[0] is the total sum of weights
        if(i > nsam*freq_cut && i < nsam*(1.0-freq_cut)) {
            ww[i]  = 1.0/(double)i; ww[0] +=  ww[i];
            wt[i]  = (double)(nsam-i); wt[0] +=  wt[i];
            wfw[i] = (double)i; wfw[0] +=  wfw[i];
            wl[i]  = (double)1; wl[0] +=  wl[i];
            if(i==1) {
                wfl[i] = (double)1; wfl[0] +=  wfl[i];
            } else wfl[i] = 0;
        }
        else {
            ww[i]  = 0;
            wt[i]  = 0;
            wfw[i] = 0;
            wl[i]  = 0;
            wfl[i] = 0;
        }
    }
    return(1);
}
int weights_unfolded_wpsi(double *wpsi_ij, int nsam, int fr, int subnsam) {
    int i;
    for(i=1;i<=subnsam;i++) {
        wpsi_ij[i] = (1.0-gsl_ran_hypergeometric_pdf(0, fr, nsam-fr, i));
    }
    return(1);
}
double Calc_rTheta_unfolded(double rfr, double w, double psi, int nsam, double sumw) {
    double th = 0.0;
    th = th + w * rfr * 1.0/psi;
    if(sumw)
        th = th/sumw;
    return(th);
}
double Calc_Theta_unfolded(double fr, double w, int nsam, double sumw) {
    //int i;
    double th = 0.0;
    th = th + w * fr;
    if(sumw)
        th = th/sumw;
    return(th);
}

/* FOLDED */
double deltak(int i, int j) {
    if(i==j) return(1.0);
    return(0.0);
}
int weights_folded(double *ww, double *wt, double *wfl, double *wl, double *wphi_i,/* double **wpsi_ij,*/ int nsam,/* int subnsam,*/ double freq_cut) {
    int i;
    for(i=1;i<=(int)floor(nsam/2);i++) {
        wphi_i[i] = ((double)nsam/((double)i*(double)(nsam-i))) / (1.0+deltak(i,nsam-i));
        if(i > nsam*freq_cut) {
            ww[i]  = (double)nsam/((double)i*(double)(nsam-i)*(1.0+deltak(i,nsam-i))); ww[0] +=  ww[i];
            wt[i]  = (double)nsam/(1+deltak(i,nsam-i)); wt[0] +=  wt[i];
            if(i==1) {wfl[i] = (double)nsam; wfl[0] +=  wfl[i];} else {wfl[i] = 0;}
            wl[i]  = (double)nsam/((double)(nsam-i)*(1.0+deltak(i,nsam-i))); wl[0] +=  wl[i];
        }
        else {
            ww[i]  = 0;
            wt[i]  = 0;
            wfl[i] = 0;
            wl[i]  = 0;
        }
    }
    return(1);
}
int weights_folded_wpsi(double *wpsi_ij, int nsam, int fr, int subnsam) {
    int i;
    for(i=1;i<=subnsam;i++) {
        wpsi_ij[i] = (1.0-gsl_ran_hypergeometric_pdf(0, fr, nsam-fr, i));
    }
    return(1);
}
double Calc_rTheta_folded(double rfr, double w, double phi, double psi, int nsam, double sumw) {
    double th = 0.0;
    th = th + w * 1.0/(1.0+deltak(rfr,nsam-rfr)) * 1.0/phi * 1.0/psi;
    if(sumw)
        th = th/sumw;
    return(th);
}
double Calc_Theta_folded(double fr, double w, double phi, int nsam, double sumw) {
    double th = 0.0;
    th = th + w * 1.0/(1.0+deltak(fr,nsam-fr)) * 1.0/phi;
    if(sumw)
        th = th/sumw;
    return(th);
}
