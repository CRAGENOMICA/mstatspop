/* 
 * File:   calcFs.h
 * Author: gvera
 *
 * Created on April 3, 2012, 1:49 PM
 */

#ifndef CALCFS_H
#define	CALCFS_H

#include "common.h"

#ifdef	__cplusplus
extern "C" {
#endif

/**
 * @brief calcFS function brief description HERE!
 * @param npops : Number of Populations
 * @param nsam : Number of Samples
 * @param statistics
 * @return ALWAYS 1, so why an int and not void?
 * @todo Check documentation
 */
int calcFs(int npops, int *nsam, struct stats *statistics);

/**
 * @brief Fs from Fu
 * @param Nsample
 * @param pi
 * @param NumAlelos
 * @return
 */
double Fs(int Nsample, double pi, int NumAlelos);

/**
 *
 * @param N
 * @param i
 * @param theta
 * @param qew_
 * @return
 */
double FunEq23Ewens(int N,int i,double theta, double *qew_);

#ifdef	__cplusplus
}
#endif

#endif	/* CALCFS_H */

