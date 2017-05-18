/*
 * util.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gvera
 *
 *  Auxiliary functions. Independent of mstatspop and common.h!
 */

#ifndef UTIL_H_
#define UTIL_H_

#include "common.h"
#include "ran1.h"

#ifdef	__cplusplus
extern "C" {
#endif

#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#define MAX(a,b)  ((a) > (b) ? (a) : (b))

/*
double min(int a, int b);
double max(int a, int b);
*/
double gammln(double zz);
double largebinomialdist(double pp, double n);

#ifdef	__cplusplus
}
#endif


#endif /* UTIL_H_ */
