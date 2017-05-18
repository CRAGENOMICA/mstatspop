/*
 *  calcR2.c
 *  mstatspop
 *
 *  Created by Sebastian Ramos-Onsins on 14/09/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "calcR2.h"

int calcR2(int npops, int *nsam, char *matrix_pol, long int length, struct stats *statistics, char *ploidy)
{
	
    int *initsq1;
    long int *unic;
    int x,y,pop1,inits,max,sumnsam;
    int freq[4];
    long int j;

    /*init*/
    for( pop1=0;pop1<npops;pop1++ ) statistics[0].R2[pop1] = -10000;
    
    if(length==0) return 1;
    
    initsq1 = (int *)calloc(npops,sizeof(int));

    sumnsam = 0;
    for(x=0;x<npops;x++) {
            initsq1[x] = sumnsam;
            sumnsam += nsam[x];
    }

    unic = (long int *)calloc(sumnsam,sizeof(long int));
	inits = 0;
	max = sumnsam;

    for( pop1=0;pop1<npops;pop1++ )
    {
        if( nsam[pop1] > 1 && length > 0)
        {
            for(j=0;j<length;j++)
            {
                inits   = initsq1[pop1];
                max     = initsq1[pop1]+nsam[pop1];
                freq[0] = freq[1] = freq[2] = freq[3] = 0;
                for( y=inits;y<max;y++)
                {
                    if(matrix_pol[j*sumnsam+y] == '0') {freq[1] += 1;freq[0] += 1;}
                    if(matrix_pol[j*sumnsam+y] == '1') {freq[2] += 1;freq[0] += 1;}
                    if(matrix_pol[j*sumnsam+y] == '-') {freq[3] += 1;}
                }
                if(freq[0])
                {
                    if(freq[1]==1)
                    {
                        /**< @todo Revisar este trozo de codigo */
                        for(y=inits;y<max;y++) if(matrix_pol[j*sumnsam+y] == '0') break;
                        unic[y] += 1;
                    }
                    else {
                            if(freq[2]==1) {
                                    for(y=inits;y<max;y++) if(matrix_pol[j*sumnsam+y] == '1') break;
                                    unic[y] += 1;
                            }
                    }
                }
            }
			if(ploidy[0] == '1') {
				statistics[0].R2[pop1] = 
						R2( unic+initsq1[pop1],
							statistics[0].thetaT[pop1],
							nsam[pop1],
							(int)statistics[0].S[pop1] );
			}
			if(ploidy[0] == '2') {
				for(y=inits;y<max;y+=2) unic[y] += unic[y+1];
				statistics[0].R2[pop1] = 
					R2( unic+initsq1[pop1],
					   2.0*statistics[0].thetaT[pop1],
					   nsam[pop1]/2,
					   (int)statistics[0].S[pop1] );
			}
        }
		else {
			statistics[0].R2[pop1] = -10000;
		}
    }

    free(initsq1);
    free(unic);

    return 1;
}


double R2(long int *unic,double pi,int sample_size,long int S)
{
    double sm2 = 0.0;
    int i;
    
    if(S == 0 || sample_size == 0) return(-10000);
    for( i=0;i<sample_size;i++ ) {
        sm2 += ((double)unic[i] - pi/2.0) * ((double)unic[i] - pi/2.0);
    }
    
    sm2 = sqrt(sm2/((double)sample_size))/(double)S;
	
    if( fabs(sm2) < 1.0E-15 )
		sm2 = 0.0;
	
    return (double)sm2;
}


int calcR2p(int npops, int *nsam, char *matrix_pol, long int length, struct stats *statistics, double *sum_sam, int *r2i_ploidies)
{
	
    int *initsq1;
    long int *unic;
    int x,y,pop1,inits,max,sumnsam;
    int freq[4];
    int freqo[4];
    long int j;
	double efflength;
	
    /*init*/
    for( pop1=0;pop1<npops;pop1++ ) {
        statistics[0].R2[pop1] = -10000;
        for(x=0;x<r2i_ploidies[0];x++) {
            statistics[0].R2p[x][pop1] = -10000;
        }
    }
    if(length==0) return 1;
    
    initsq1 = (int *)calloc(npops,sizeof(int));
	
    sumnsam = 0;
    for(x=0;x<npops;x++) {
		initsq1[x] = sumnsam;
		sumnsam += nsam[x];
    }
	
    unic = (long int *)calloc(sumnsam,sizeof(long int));
	inits = 0;
	max = sumnsam;
	
    for( pop1=0;pop1<npops;pop1++ )
    {
        if( nsam[pop1] > 1 && length > 0)
        {
			inits   = initsq1[pop1];
			max     = initsq1[pop1]+nsam[pop1];
            
			for(j=0;j<length;j++)
            {
                /*eliminate those positions with no outgroup*/
                freqo[0]=freqo[1]=freqo[2]=freqo[3]=0;
                for(y=initsq1[npops-1];y<sumnsam;y++) {
                    if(matrix_pol[j*sumnsam+y] == '0') {freqo[1] += 1;freqo[0] += 1;}
                    if(matrix_pol[j*sumnsam+y] == '1') {freqo[2] += 1;freqo[0] += 1;}
                    if(matrix_pol[j*sumnsam+y] == '-') {freqo[3] += 1;}
                }
                if(freqo[0]) {
                    if(freqo[1] != freqo[0] && freqo[1] != 0) {
                        continue;/*if the outgroup is polymorphic, we do not consider*/
                    }
                }
                else
                    continue;
                /*end outgroup*/

                freq[0] = freq[1] = freq[2] = freq[3] = 0;
                for( y=inits;y<max;y++)
                {
                    if(matrix_pol[j*sumnsam+y] == '0') {freq[1] += 1;freq[0] += 1;}
                    if(matrix_pol[j*sumnsam+y] == '1') {freq[2] += 1;freq[0] += 1;}
                    if(matrix_pol[j*sumnsam+y] == '-') {freq[3] += 1;}
                }
                if(freq[0])
                {
                    if(freq[1]==1)
                    {
                        /**< @todo Revisar este trozo de codigo */
                        for(y=inits;y<max;y++) if(matrix_pol[j*sumnsam+y] == '0') break;
                        unic[y] += 1;
                    }
                    else {
						if(freq[2]==1) {
							for(y=inits;y<max;y++) if(matrix_pol[j*sumnsam+y] == '1') break;
							unic[y] += 1;
						}
                    }
                }
            }
			
			efflength = 0.;
			for( y=inits;y<max;y++)
				efflength += sum_sam[y];
            
            /**calculate sum_sam considering the otgroup!!!!!*/
			
			for(x=0;x<r2i_ploidies[0];x++) {
				statistics[0].R2p[x][pop1] = R2p(unic+initsq1[pop1],
												statistics[0].thetaTo[pop1],
												floor(nsam[pop1]/(double)r2i_ploidies[x+1]),
												(int)statistics[0].So[pop1],sum_sam+initsq1[pop1],efflength,r2i_ploidies[x+1]);
				/*
				if(x==1 && ploidy[0] == '1') statistics[0].R2[pop1] = statistics[0].R2p[x][pop1];
				if(x==2 && ploidy[0] == '2') statistics[0].R2[pop1] = statistics[0].R2p[x][pop1];
				*/
			}
        }
		else {
			statistics[0].R2[pop1] = -10000;
			/*
			if(x==1 && ploidy[0] == '1') statistics[0].R2[pop1] = -10000;
			if(x==2 && ploidy[0] == '2') statistics[0].R2[pop1] = -10000;
			*/
		}
    }
	
    free(initsq1);
    free(unic);
	
    return 1;
}


double R2p(long int *unic,double pi,int sample_size,long int S, double *efflength_sam, double efflength, int ploidy)
{
    double sm2 = 0.0;
    int i,j;
    long int unicp;
	double effposi;
	
    if(S == 0 || sample_size == 0) return(-10000);
    for( i=0;i<sample_size*ploidy;i+=ploidy ) {
		unicp=0;
		effposi=0.0;
		for(j=i;j<i+ploidy;j++) {
			unicp += unic[j];
			effposi += (double)efflength_sam[j];
        }
		sm2 += ((double)unicp - effposi/((double)efflength) * pi/2.0*(double)ploidy) * 
			   ((double)unicp - effposi/((double)efflength) * pi/2.0*(double)ploidy);
    }
    sm2 = sqrt(sm2) / (double)S;
	
    if( fabs(sm2) < 1.0E-15 )
		sm2 = 0.0;
	
    return (double)sm2;
}
