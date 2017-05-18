/*
 * sancestral.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gvera
 */

#ifndef SANCESTRAL_H_
#define SANCESTRAL_H_

#include "common.h"

#ifdef	__cplusplus
extern "C" {
#endif


int calc_sxsfss(int npops, int *nsam, char *matrix_pol,long int *matrix_pos,
		long int length, struct stats *statistics,long int *sites_matrix,int outgroup_presence,int force_outgroup);
int ispolmis(char *matrix_pol,int nsam, long int j,long int length);


#ifdef	__cplusplus
}
#endif

#endif /* SANCESTRAL_H_ */
