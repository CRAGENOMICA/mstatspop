/*
 *  permute.c
 *  mstatspop
 *
 *  Created by Sebastian Ramos-Onsins on 15/08/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "permute.h"

int permute( char *matrix_pol,long int length,int nsamtot,char *matrix_perm,int *nsam2,int *psam2,int npops/*,int *nsam*/,int nsam_outg, int psam_outg)
{
	/*length is the number of variants in matrix_pol
	 nsamtot is the number of samples
	 matrix_pol is a vector with data: (length x nsamtot)
	 matrix_perm is the new vector with the data (samples permuted)
	 psam2 and nsam2 are included to choose 2 populations in order to do Fst
	*/
	/*nsam2 is the samplesize of each of the 2 pops. psam2 is the position (cummulative value from nsamtot) of each of the 2 pops. 
	 ex: pop1=10 pop2 =15. nsam2[0]=10 nsam2[1]=15 psam2[0]=0 psam2[1]=10.*/
	double ran1(void);
	
	int *count;
	int sumnsam2;
	int x,y,c;
	long int j;
	
	count = (int *)calloc(nsamtot,sizeof(int));
	
	/*samples not considered are 0, considered are 1 in count vector*/
	for(y=0;y<nsamtot;y++) count[y] = 0;
	for(x=0;x<2;x++) {
		for(y=psam2[x];y<psam2[x]+nsam2[x];y++) {
			count[y] = 1;
		}
	}
	
	/*permute*/
    /*printf("\n");*/
	sumnsam2 = nsam2[0]+nsam2[1];
	while(sumnsam2) {
		c = (int)floor((ran1()*(double)sumnsam2));
        for(y=0;y<nsamtot;y++) {
			if(count[y] == 1) c--;
			if(c<0)
                break;
		}
		for(j=0;j<length;j++) {
			/*strncpy(matrix_perm+x*nsamtot+((nsam2[0]+nsam2[1])-sumnsam2),matrix_pol+x*nsamtot+y,1);*/
			matrix_perm[j*(nsam2[0]+nsam2[1]+nsam_outg)+((nsam2[0]+nsam2[1])-sumnsam2)] = matrix_pol[j*nsamtot+y];
            /*printf("%c",matrix_perm[j*(nsam2[0]+nsam2[1])+((nsam2[0]+nsam2[1])-sumnsam2)]);*/
		}
        /*printf("\n");*/
		count[y] = 0;
		sumnsam2--;
	}
	/*
	c=0;
	for(y=nsamtot-nsam[npops-1];y<nsamtot;y++,c++)
		for(x=0;x<length;x++)
			strncpy(matrix_perm+x*nsamtot+(nsam2[0]+nsam2[1])+c,matrix_pol+x*nsamtot+y,1);
	*/
	free(count);

    for(y=psam_outg,x=0;y<psam_outg+nsam_outg;y++,x++) {
        for(j=0;j<length;j++) {
            matrix_perm[j*(nsam2[0]+nsam2[1]+nsam_outg)+sumnsam2+x] = matrix_pol[j*nsamtot+y];
        }
    }

    return 1;
}
