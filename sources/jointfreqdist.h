/*
 * jointfreqdist.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gvera
 */

#ifndef JOINTFREQDIST_H_
#define JOINTFREQDIST_H_

#include "common.h"

#ifdef	__cplusplus
extern "C" {
#endif


int jointfreqdist(int npops, int *nsam, char *matrix_pol,long int *matrix_pos, long int length, struct stats *statistics,long int *sites_matrix, double **jfd, int **nfd,int outgroup_presence,int force_outgroup);

#ifdef	__cplusplus
}
#endif



#endif /* JOINTFREQDIST_H_ */
