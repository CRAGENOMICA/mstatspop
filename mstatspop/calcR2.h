/* 
 * File:   calcR2.h
 * Author: gvera
 *
 * Created on April 3, 2012, 2:46 PM
 */

#ifndef CALCR2_H
#define	CALCR2_H

#include "common.h"

#ifdef	__cplusplus
extern "C" {
#endif

/**
 * @brief calcR2
 * @details bla cla
 * @param npops
 * @param nsam
 * @param matrix_pol
 * @param length
 * @param statistics
 * @return
 * @todo Review documentation
 */
int calcR2(int npops, int *nsam, char *matrix_pol, long int length, struct stats *statistics, char *);

/**
 * @brief
 * @details R2 Ramos & Rozas: "*unic" is the number of singletons in each sequence (in comparison to the sample studied)
 * @param unic
 * @param pi
 * @param sample_size
 * @param S
 * @return
 * @todo Review documentation
 */
double R2(long int *unic,double pi,int sample_size,long int S);

int calcR2p(int npops, int *nsam, char *matrix_pol, long int length, struct stats *statistics, double *sum_sam, int *r2i_ploidies,int outgroup);
double R2p(long int *unic,double pi,int sample_size,long int S, double *efflength_sam, double efflength, int ploidy);

#ifdef	__cplusplus
}
#endif

#endif	/* CALCR2_H */

