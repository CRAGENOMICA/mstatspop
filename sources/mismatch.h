/*
 * mismatch.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gvera
 */

#ifndef MISMATCH_H_
#define MISMATCH_H_

#include "common.h"

#ifdef	__cplusplus
extern "C" {
#endif

int calc_mismatch(int npops, int *nsam, char *matrix_pol,long int length, struct stats *statistics, char *ploidy);

#ifdef	__cplusplus
}
#endif



#endif /* MISMATCH_H_ */
