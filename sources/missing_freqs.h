#ifndef MISSING_FREQS_H
#define	MISSING_FREQS_H

#include "common.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>


#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) > (Y) ? (X) : (Y))

#define INF 10000
#define CONFSIZE 10000 //size of configuration array for Hellman's probability

#ifdef	__cplusplus
extern "C" {
#endif
	
int min3(int a, int b, int c);
int max3(int a, int b, int c);
int myround(float r);
float lnfact(int n);
float lnbinomial(int n, int k);
float lnmultinomial(int n, int k1, int k2);
float a(int n);
float beta(int n,int i);
float ps(int k, int l, int n, float theta);
float sigma(int i, int j,int n);
float pe(int k, int l, int n, float theta);
float cs(int i, int j, int k, int l, int nx, int ny, int nxy,float theta);
float ce(int i, int j, int k, int l, int nx, int ny, int nxy,float theta);
float p(int i, int j, int nx, int ny, int nxy, float theta);
float cov_missing(int i, int j, int nx, int ny, int nxy, float theta);

float hellman_sum(int *conf,int nconf, int noind, int *r);
int check_conf(int j, int n,int* conf);
float pc(int j,int *r,int n);

float hn(int n);
void print_conf_matrix(int *conf,int nconf, int n);
float watterson(int* a, int n);
float corrected_watterson(int*a, int* r, int n);
void store_configurations(int n,int s ,int* conf,int *row);
float pc(int j,int *r,int n);
int check_conf(int j, int n,int* conf);
float hellman_prod(int noind, int* r,int* i);
float hellman_sum(int *conf,int nconf, int noind, int *r);
float ck(int r);
int ak(int m,int r);
float fhat(int* a, float*c, int n,float* var);
void _log(char *s);
float tajima_pi(int *a,float *c, int n);
float watterson_variance(float theta, int l, int* nx, int* ny,int* nxy, float* var_d, float* var_h);	

#ifdef	__cplusplus
}
#endif

#endif	/* MISSING_FREQS_H */

