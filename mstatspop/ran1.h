/*
 * ran1.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gvera
 */

#ifndef RAN1_H_
#define RAN1_H_

#include "common.h"

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <math.h>

#define IA 16807
#define IM 2147483647
#define IQ 127773
#define IR 2836
#define NTAB 32
#define EPS (1.2E-07)

void init_seed1(long int seed);
double ran1(void);

#ifdef	__cplusplus
}
#endif

#endif /* RAN1_H_ */
