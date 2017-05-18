/*
 * optimal-RinfR0.h
 *
 *  Created on: Apr 18, 2012
 *      Author: gvera
 */

#ifndef OPTIMAL_RINFR0_H_
#define OPTIMAL_RINFR0_H_

#include "common.h"

#ifdef	__cplusplus
extern "C" {
#endif


/**< todo: inGSL tendra que definirse en tiempo de compilacion, en el Makefile, no aqui, no??? */
#ifndef inGSL
	#define inGSL 1
#endif

/* Include libraries */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#if inGSL
	#include <gsl/gsl_rng.h>
	#include <gsl/gsl_randist.h>
	#include <gsl/gsl_sf.h>
	#include <gsl/gsl_math.h>

	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_block.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
#endif

#define DEBUG_DEN 0

int generate_optimal_test_rinf_gennull(struct test_rinf *opt,int n,double *xi_b, double *xi_null);
int generate_optimal_test_rinf(struct test_rinf *opt,int n,double *xi_b) ;
double optimal_test_rinf(struct test_rinf *t_rinf, double *xi);
double optimal_test_rinf(struct test_rinf *t_rinf, double *xi);
int generate_optimal_test_r0(struct test_r0 *opt,int n, double *xi_b,/* temporary: double theta*/ long sw, long l) ;
int generate_optimal_test_r0_gennull(struct test_r0 *opt,int n, double *xi_b,double *xi_null, /* temporary: double theta*/ long sw, long l);
double optimal_test_r0(struct test_r0 *t_r0, double *xi);
int generate_optimal_test_rinf0 (struct test_rinf *opt,int n, double *xi_b, long sw, long l);
double optimal_test_rinf_rvar0 (struct test_rinf *t_rinf, double *xi);
int generate_optimal_test_quad (struct test_quad *opt,int n, double *xi_b,double theta, long l);
double optimal_test_quad(struct test_quad *t_quad, double *xi);
int generate_optimal_test_lc (struct test_lc *opt, int n, double *xi_b,double theta, long l);
double optimal_test_lc(struct test_lc *t_lc, double *xi);
int generate_optimal_test_quad_wc (struct test_quad_wc *opt, int n, double *xi_b,double theta, long l);
double optimal_test_quad_wc(struct test_quad_wc *t_quad_wc, double *xi);


#ifdef	__cplusplus
}
#endif

#endif /* OPTIMAL_RINFR0_H_ */
