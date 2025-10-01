//
//  rfreq_stats.h
//  mstatspop
//
//  Created by Sebastian Ramos-Onsins on 11/6/25.
//

#ifndef RFREQ_STATS_H
#define RFREQ_STATS_H

#include "common.h"
#include "ran1.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#ifdef    __cplusplus
extern "C" {
#endif


int calc_rfreqstats(int npops, int *nsam, char *matrix_pol,long int length, struct stats *statistics,
                   struct rweights *rw, int outgroup_presence,int force_outgroup,int include_unknown);
int weights_unfolded(double *ww, double *wt, double *wfw, double *wfl, double *wl, /*double **wpsi_ij, */int nsam, /*int subnsam,*/ double freq_cut);
int weights_unfolded_wpsi(double *wpsi_ij, int nsam, int fr, int subnsam);

double Calc_rTheta_unfolded(double rfr, double w, double psi, int nsam, double sumw);

double Calc_Theta_unfolded(double fr, double w, int nsam, double sumw);

double deltak(int i, int j);

int weights_folded(double *ww, double *wt, double *wfl, double *wl, double *wphi_i, /*double **wpsi_ij, */ int nsam, /*int subnsam, */double freq_cut);
int weights_folded_wpsi(double *wpsi_ij, int nsam, int fr, int subnsam);

double Calc_rTheta_folded(double rfr, double w, double phi, double psi, int nsam, double sumw);

double Calc_Theta_folded(double fr, double w, double phi, int nsam, double sumw);



#ifdef    __cplusplus
}
#endif

#endif    /* RFREQ_STATS_H */

