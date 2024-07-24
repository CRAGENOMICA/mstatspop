/*
 *  calc_Toptimal_tests.c 
 *  mstatspop
 *  Sebastian E. Ramos-Onsins
 *
 *  Algorithms of optimal tests developed by Luca Ferretti and Giacomo Marmorini (2010)
 *
 */

#include "calc_Toptimal_tests.h"
#include "optimal-RinfR0.h"

int	calc_Toptimal_tests(int npops,int *nsam,struct stats *statistics)
{
	int x,y;
	double 	*mean_freqsptr,
			*data_freqsptr,
			*H0freqsptr,
			mean_S,
			data_S,
			H0_S,
			an,
			thetaw;
	
	struct test_rinf *omega_test_rinf,*omega_test_rinfH0,*omega_test_rinfvar0;
	struct test_r0 *omega_test_r0,*omega_test_r0H0;
	struct test_quad *omega_test_quad;
	struct test_lc *omega_test_line;
	struct test_quad_wc *omega_test_quad_wc;

	for(x=0;x<npops-1;x++)
	{
		omega_test_rinf = (struct test_rinf *)calloc(1,sizeof(struct test_rinf));	
		omega_test_rinf->w=(double *)malloc((nsam[x]-1)*sizeof(double));
		
		omega_test_rinfH0 = (struct test_rinf *)calloc(1,sizeof(struct test_rinf));	
		omega_test_rinfH0->w=(double *)malloc((nsam[x]-1)*sizeof(double));
		
		omega_test_rinfvar0 = (struct test_rinf *)calloc(1,sizeof(struct test_rinf));	
		omega_test_rinfvar0->w=(double *)malloc((nsam[x]-1)*sizeof(double));

		omega_test_r0H0=(struct test_r0 *)malloc(sizeof(struct test_r0));
		omega_test_r0H0->w=(double *)malloc((nsam[x]-1)*sizeof(double));
		omega_test_r0H0->dw=(double *)malloc((nsam[x]-1)*sizeof(double));
		omega_test_r0H0->d2w=(double *)malloc((nsam[x]-1)*sizeof(double));

		omega_test_r0=(struct test_r0 *)malloc(sizeof(struct test_r0));
		omega_test_r0->w=(double *)malloc((nsam[x]-1)*sizeof(double));
		omega_test_r0->dw=(double *)malloc((nsam[x]-1)*sizeof(double));
		omega_test_r0->d2w=(double *)malloc((nsam[x]-1)*sizeof(double));

		omega_test_quad = (struct test_quad *)calloc(1,sizeof(struct test_quad));	
		omega_test_quad->w=(double *)malloc((nsam[x]-1)*sizeof(double));
		omega_test_quad->w2=(double **)malloc((nsam[x]-1)*sizeof(double *));
		for(y=0;y<nsam[x]-1;y++) 
			omega_test_quad->w2[y]=(double *)malloc((nsam[x]-1)*sizeof(double));
		
		omega_test_quad_wc = (struct test_quad_wc *)calloc(1,sizeof(struct test_quad_wc));	
		omega_test_quad_wc->w=(double *)malloc((nsam[x]-1)*sizeof(double));
		omega_test_quad_wc->w2=(double **)malloc((nsam[x]-1)*sizeof(double *));
		for(y=0;y<nsam[x]-1;y++) 
			omega_test_quad_wc->w2[y]=(double *)malloc((nsam[x]-1)*sizeof(double));
		
		omega_test_line = (struct test_lc   *)calloc(1,sizeof(struct test_lc  ));	
		omega_test_line->w=(double *)malloc((nsam[x]-1)*sizeof(double));
				
		/*frequency spectrum for H1, H0 and data*/
		mean_freqsptr = (double *)calloc(nsam[x],sizeof(double));
		for(y=0;y<nsam[x]-1;y++)  mean_freqsptr[y] = (double)statistics[0].H1freq[x][y+1];
		
		data_freqsptr = (double *)calloc(nsam[x],sizeof(double));
		for(y=0;y<nsam[x]-1;y++)  data_freqsptr[y] = (double)statistics[0].freq[x][y+1];
		
		H0freqsptr = (double *)calloc(nsam[x],sizeof(double));
		for(y=0;y<nsam[x]-1;y++)  H0freqsptr[y] = (double)statistics[0].H0freq[x][y+1];
		
		mean_S = 0; for(y=0;y<nsam[x]-1;y++) mean_S += (double)mean_freqsptr[y];
		data_S = 0; for(y=0;y<nsam[x]-1;y++) data_S += data_freqsptr[y];
		H0_S   = 0; for(y=0;y<nsam[x]-1;y++) H0_S   += H0freqsptr[y];
		
		/*NO MISSING VALUES*/
		an     = 0.;for(y=1;y<nsam[x];y++) an += 1./(double)y;
		thetaw = data_S/(an * statistics[0].total_length);/*theta(data)*/
		
		/*SUM TO 1*/
		for(y=0;y<nsam[x]-1;y++) mean_freqsptr[y] = mean_freqsptr[y]/mean_S;
		for(y=0;y<nsam[x]-1;y++) H0freqsptr[y] = H0freqsptr[y]/H0_S;
		
		/*rinfvarinf_H0nonnull*/
		if(generate_optimal_test_rinf_gennull(omega_test_rinfH0,nsam[x],mean_freqsptr,H0freqsptr) == (double)-10000) 
			statistics[0].ToH0_ii[x] = (double)-10000;
		else
			statistics[0].ToH0_ii[x] = optimal_test_rinf(omega_test_rinfH0,data_freqsptr);

		/*rinfvarinf*/
		if(generate_optimal_test_rinf(omega_test_rinf,nsam[x],mean_freqsptr) == (double)-10000) 
			statistics[0].To_ii[x] = (double)-10000;
		else
			statistics[0].To_ii[x] = optimal_test_rinf(omega_test_rinf,data_freqsptr);
		
		/*rinfvar0*/
		if(generate_optimal_test_rinf0(omega_test_rinfvar0,nsam[x],mean_freqsptr,data_S,statistics[0].total_length) == (double)-10000) 
			statistics[0].To_i0[x] = (double)-10000;
		else 
			statistics[0].To_i0[x] = optimal_test_rinf_rvar0(omega_test_rinfvar0,data_freqsptr);
		
		/*r0var0*/
		if(generate_optimal_test_r0(omega_test_r0,nsam[x],mean_freqsptr,data_S,statistics[0].total_length) == (double)-10000) 
			statistics[0].To_00[x] = (double)-10000;
		else 
			statistics[0].To_00[x] = optimal_test_r0(omega_test_r0,data_freqsptr);
		
		/*r0var0H0nonnull*/
		if(generate_optimal_test_r0_gennull(omega_test_r0H0,nsam[x],mean_freqsptr,H0freqsptr,data_S,statistics[0].total_length) == (double)-10000) 
			statistics[0].ToH0_00[x] = (double)-10000;
		else 
			statistics[0].ToH0_00[x] = optimal_test_r0(omega_test_r0H0,data_freqsptr);					

		/*ALTERNATIVE FREQUENCY SPECTRUM CORRECTION: SUM TO S(data)/theta(H1)*/
		for(y=0;y<nsam[x]-1;y++) 
			mean_freqsptr[y] = mean_freqsptr[y] * mean_S/(statistics[0].thetaH1[x]*statistics[0].total_length);
		
		/*rquad*/
		if(generate_optimal_test_quad(omega_test_quad,nsam[x],mean_freqsptr,thetaw,statistics[0].total_length) == (double)-10000) 
			statistics[0].To_Qc_ii[x] = (double)-10000;
		else
			statistics[0].To_Qc_ii[x] = optimal_test_quad(omega_test_quad,data_freqsptr);					

		/*rquad_wc*/
		if(generate_optimal_test_quad_wc(omega_test_quad_wc,nsam[x],mean_freqsptr,thetaw,statistics[0].total_length) == (double)-10000) 
			statistics[0].To_Qw_ii[x] = (double)-10000;
		else
			statistics[0].To_Qw_ii[x] = optimal_test_quad_wc(omega_test_quad_wc,data_freqsptr);					

		/*rline*/
		if(generate_optimal_test_lc(omega_test_line,nsam[x],mean_freqsptr,thetaw,statistics[0].total_length) == (double)-10000) 
			statistics[0].To_Lc_ii[x] = (double)-10000;
		else 
			statistics[0].To_Lc_ii[x] = optimal_test_lc(omega_test_line,data_freqsptr);	
		
		free(mean_freqsptr);
		free(data_freqsptr);
		
		free(omega_test_rinf->w);
		free(omega_test_rinf);
		free(omega_test_rinfH0->w);
		free(omega_test_rinfH0);
		free(omega_test_rinfvar0->w);
		free(omega_test_rinfvar0);

		free(omega_test_r0->w);
		free(omega_test_r0->dw);
		free(omega_test_r0->d2w);
		free(omega_test_r0);
		
		free(omega_test_r0H0->w);
		free(omega_test_r0H0->dw);
		free(omega_test_r0H0->d2w);
		free(omega_test_r0H0);
		
		free(omega_test_quad->w);
		for(y=0;y<nsam[x]-1;y++) free(omega_test_quad->w2[y]);
		free(omega_test_quad->w2);
		free(omega_test_quad);
		
		free(omega_test_quad_wc->w);
		for(y=0;y<nsam[x]-1;y++) free(omega_test_quad_wc->w2[y]);
		free(omega_test_quad_wc->w2);
		free(omega_test_quad_wc);
		
		free(omega_test_line->w);
		free(omega_test_line);
        
        free(H0freqsptr);
	}
	
	return 1;
}

