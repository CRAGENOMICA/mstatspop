/*
 * fstcalc.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gvera
 */
 
#ifndef FSTCALC_H_
#define FSTCALC_H_

#include "common.h"

#ifdef	__cplusplus
extern "C" {
#endif


#ifndef inMISSINGH
	#define MISSINGH 0
#else
	#define MISSINGH 2./32.
#endif


int calc_piwpiafst(int flag, int formatfile,int npops, int *nsam, char *matrix_pol,
                       long int length, struct stats *statistics, long int *matrix_sv,
                       int outgroup_presence, int force_outgroup);
int calc_hwhafsth(int npops, int *nsam, char *matrix_pol,long int length, struct stats *statistics);


#ifdef	__cplusplus
}
#endif

#endif /* FSTCALC_H_ */
