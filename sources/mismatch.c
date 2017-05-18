/*
 *  mismatch.c
 *  mstatspop
 *
 *  Created by Sebastian Ramos-Onsins on 26/05/2010.
 *  Copyright 2010 CRAG. All rights reserved.
 *
 */
#include "mismatch.h"
#include "sancestral.h"

int calc_mismatch(int npops, int *nsam, char *matrix_pol,long int length, struct stats *statistics, char *ploidy)
{
	/*
		Error: using missing values the weight (comb) is not correctly calculated. Do it for each position independently...
	*/
	int x,y;
	long int z,zz,ncomb;
	long int j;
	int pop1;
    double thetaT=0.0;
    double w,moment;
	double *freq1,*freq2,*freq3,*freq4,s,s2,g1,g2;
	long int sumnsam,maxnsam,initso;
	int *initsq1,inits,max,a0;
	
	for(pop1=0;pop1<npops-1;pop1++){
		statistics[0].mdg1[pop1] = -10000.;
		statistics[0].mdg2[pop1] = -10000.;
	}
    
    initsq1 = (int *)calloc(npops,sizeof(int));
	sumnsam = 0;
	maxnsam = 0;
	
	for(x=0;x<npops;x++) {
		initsq1[x] = (int)sumnsam;
		sumnsam += nsam[x];
		if(maxnsam < nsam[x]) maxnsam = nsam[x];
	}
	    
    freq1 = (double *)calloc(maxnsam*(maxnsam-1)/2,sizeof(double));
	freq2 = (double *)calloc(maxnsam*(maxnsam-1)/2,sizeof(double));
	freq3 = (double *)calloc(maxnsam*(maxnsam-1)/2,sizeof(double));
	freq4 = (double *)calloc(maxnsam*(maxnsam-1)/2,sizeof(double));
	
	for(pop1=0;pop1<npops-1;pop1++)
	{
		if(ploidy[0] == '1') 
			thetaT = statistics[0].thetaTo[pop1];
		if(ploidy[0] == '2') {
			thetaT = 0.;
			w = 0.;
			inits = initsq1[pop1];
			max   = initsq1[pop1]+nsam[pop1];
			for(x=inits;x<max-1;x+=2)  {
				for(y=x+1;y<max;y+=2)  {
					zz = 0;
					for(j=0;j<length;j++) {		
                        /*outgroup*/
                        initso = sumnsam-nsam[npops-1];
                        while((a0 = matrix_pol[j*sumnsam+initso]) == '-' && initso < sumnsam) initso++;
                        if(initso >= sumnsam-nsam[npops-1]) {
                            continue;
                        }
                        /*polymorphism*/
						if((ispolmis(matrix_pol,(int)sumnsam,j,length)) > 0)
							continue;
						else {
							if(matrix_pol[j*sumnsam+x] != '-' && matrix_pol[j*sumnsam+y] != '-') {
								zz++;
								if(matrix_pol[j*sumnsam+x] != matrix_pol[j*sumnsam+y]) {
									thetaT += 1;
								}
							}
						}
					}
					w+=(double)zz/(double)length;
				}
			}
			thetaT /= (double)w;
		}
		
		if(ploidy[0] == '1') {
			ncomb = nsam[pop1]*(nsam[pop1]-1)/2;
			for(z=0;z<ncomb;z++)  {
				freq1[z] = 0.;
				freq2[z] = 0.;
				freq3[z] = 0.;
				freq4[z] = 0.;
			}
			
			z = 0;
			inits = initsq1[pop1];
			max   = initsq1[pop1]+nsam[pop1];
			for(x=inits;x<max-1;x++)  {
				for(y=x+1;y<max;y++)  {
					for(j=0;j<length;j++) {		
                        /*outgroup*/
                        initso = sumnsam-nsam[npops-1];
                        while((a0 = matrix_pol[j*sumnsam+initso]) == '-' && initso < sumnsam) initso++;
                        if(initso >= sumnsam-nsam[npops-1]) {
                            continue;
                        }
						if((ispolmis(matrix_pol,(int)sumnsam,j,length)) > 0)
							continue;
						else {
							if(matrix_pol[j*sumnsam+x] != matrix_pol[j*sumnsam+y] && matrix_pol[j*sumnsam+x] != '-' && matrix_pol[j*sumnsam+y] != '-') {
								if(ploidy[0] == '1') freq1[z] += 1;
							}
						}
					}
					moment = (freq1[z] - thetaT);
					freq2[z] = moment * moment;
					freq3[z] = moment * moment * moment;
					freq4[z] = moment * moment * moment * moment;
					statistics[0].mdw[pop1][z] = freq1[z]; /*mismatch distribution*/
					z++;
				}
			}
			s2 = g1 = g2 = 0.;				
			for(z=0;z<ncomb;z++)  {
				s2 += freq2[z];
				g1 += freq3[z];
				g2 += freq4[z];
			}
            if(ncomb > 1) {
                s = sqrt(s2/((double)ncomb-1));
                statistics[0].mdsd[pop1] = s;
            }
            else {
                statistics[0].mdsd[pop1] = -10000.;
                s = -1;
            }
            if(s>=0) {
				statistics[0].mdg1[pop1] = g1 * (double)ncomb/(((double)ncomb-2)*((double)ncomb-1)*s*s*s);
				statistics[0].mdg2[pop1] = g2 * ((double)ncomb-1) * ncomb /(((double)ncomb-3)*((double)ncomb-2)*((double)ncomb-1)*s*s*s*s)
										  - (3. * ((double)ncomb-1)*((double)ncomb-1))/(((double)ncomb-3)*((double)ncomb-2));
			}
			else {
				statistics[0].mdg1[pop1] = -10000.;
				statistics[0].mdg2[pop1] = -10000.;
			}
		}
		if(ploidy[0] == '2') {
			ncomb = (nsam[pop1])*((nsam[pop1])-1)/2;
			for(z=0;z<ncomb;z++)  {
				freq1[z] = 0.;
				freq2[z] = 0.;
				freq3[z] = 0.;
				freq4[z] = 0.;
			}
			
			z = 0;
			inits = initsq1[pop1];
			max   = initsq1[pop1]+nsam[pop1];
			w = 0.;
			for(x=inits;x<max-1;x+=1)  {
				for(y=x+1;y<max;y+=1)  {
					if(x/2 == (double)x/2 && y==x+1) 
						continue;
					zz = 0;
					for(j=0;j<length;j++) {		
                        /*outgroup*/
                        initso = sumnsam-nsam[npops-1];
                        while((a0 = matrix_pol[j*sumnsam+initso]) == '-' && initso < sumnsam) initso++;
                        if(initso >= sumnsam-nsam[npops-1]) {
                            continue;
                        }
						if((ispolmis(matrix_pol,(int)sumnsam,j,length)) > 0)
							continue;
						else {
							if(matrix_pol[j*sumnsam+x] != '-' && matrix_pol[j*sumnsam+y] != '-') {
								zz++;
								if(matrix_pol[j*sumnsam+x] != matrix_pol[j*sumnsam+y]) {
									freq1[z] += 1;
								}
							}
						}
					}
					w+=(double)zz/(double)length;
					moment = (freq1[z] - thetaT);
					freq2[z] = moment * moment;
					freq3[z] = moment * moment * moment;
					freq4[z] = moment * moment * moment * moment;
					statistics[0].mdw[pop1][z] = freq1[z]; /*mismatch distribution*/
					z++;
				}
			}
			s2 = g1 = g2 = 0.;				
			for(z=0;z<ncomb;z++)  {
				s2 += freq2[z];
				g1 += freq3[z];
				g2 += freq4[z];
			}
			/*ncomb = (long int)w; CHANGE ncomb by w: IMPORTANT!! ???? <<<<---------------*/
            if(w > 1) {
                s = sqrt(s2/((double)w-1));
                statistics[0].mdsd[pop1] = s;
            }
            else {
                statistics[0].mdsd[pop1] = -10000.;
                s = -1;
            }
            if(s>=0) {
				statistics[0].mdg1[pop1] = g1 * w/((w-2)*(w-1)*s*s*s);
				statistics[0].mdg2[pop1] = g2 * (w-1) * w /((w-3)*(w-2)*(w-1)*s*s*s*s) - (3. * (w-1)*(w-1))/((w-3)*(w-2));
			}
			else {
				statistics[0].mdg1[pop1] = -10000.;
				statistics[0].mdg2[pop1] = -10000.;
			}
		}
	}
	
	free(initsq1);
	free(freq1);
	free(freq2);
	free(freq3);
	free(freq4);
	
	return 1;
}
