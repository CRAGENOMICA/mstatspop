/*
 * permute.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gvera
 */

#ifndef PERMUTE_H_
#define PERMUTE_H_

#include "common.h"

#ifdef	__cplusplus
extern "C" {
#endif

int permute( char *matrix_pol,long int length,int nsamtot,char *matrix_perm,int *nsam2,int *psam2,int npops/*,int *nsam*/,int nsam_outg, int psam_outg);

#ifdef	__cplusplus
}
#endif



#endif /* PERMUTE_H_ */
