/*
 *  print_output.c
 *  mstatspop
 *
 *  Created by Sebastian Ramos-Onsins on 10/08/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include "print_output.h"

static char tripletsU[64][3] =
{
	{"UUU"},/*0*/	{"UUC"},/*1*/	{"UUA"},/*2*/	{"UUG"},/*3*/	{"UCU"},/*4*/	{"UCC"},/*5*/	{"UCA"},/*6*/	{"UCG"},/*7*/
	{"UAU"},/*8*/	{"UAC"},/*9*/	{"UAA"},/*10*/	{"UAG"},/*11*/	{"UGU"},/*12*/	{"UGC"},/*13*/	{"UGA"},/*14*/	{"UGG"},/*15*/
	{"CUU"},/*16*/	{"CUC"},/*17*/	{"CUA"},/*18*/	{"CUG"},/*19*/	{"UCU"},/*20*/	{"CCC"},/*21*/	{"CCA"},/*22*/	{"CCG"},/*23*/
	{"CAU"},/*24*/	{"CAC"},/*25*/	{"CAA"},/*26*/	{"CAG"},/*27*/	{"CGU"},/*28*/	{"CGC"},/*29*/	{"CGA"},/*30*/	{"CGG"},/*31*/
	{"AUU"},/*32*/	{"AUC"},/*33*/	{"AUA"},/*34*/	{"AUG"},/*35*/	{"ACU"},/*36*/	{"ACC"},/*37*/	{"ACA"},/*38*/	{"ACG"},/*39*/
	{"AAU"},/*40*/	{"AAC"},/*41*/	{"AAA"},/*42*/	{"AAG"},/*43*/	{"AGU"},/*44*/	{"AGC"},/*45*/	{"AGA"},/*46*/	{"AGG"},/*47*/
	{"GUU"},/*48*/	{"GUC"},/*49*/	{"GUA"},/*50*/	{"GUG"},/*51*/	{"GCU"},/*52*/	{"GCC"},/*53*/	{"GCA"},/*54*/	{"GCG"},/*55*/
	{"GAU"},/*56*/	{"GAC"},/*57*/	{"GAA"},/*58*/	{"GAG"},/*59*/	{"GGU"},/*60*/	{"GGC"},/*61*/	{"GGA"},/*62*/	{"GGG"},/*63*/
};

/* prints out ALL the results... */
int print_output( int mainargc,int npops,int *nsam,
					FILE *file_out, SGZip *file_out_gz,
                    char *file_input, char *file_output,
					int gfffiles, char *file_GFF, char *subset_positions,
					char *code_name, char *genetic_code,
					long int length, long int length_seg,
					double length_al,long int length_al_real,
					struct stats *statistics, struct probs *piter,
					long int niter, long int *sites_matrix, char *ploidy,
					double svratio, double missratio, int include_unknown,
					long int *matrix_pos,
					double **jfd, int **nfd, int output, int H1frq,int H0frq,
					long int nseed, char *file_H1f, char *file_H0f,
					double *vector_priors, int npriors, int formatfile,
					int outgroup_presence,int force_outgroup, double freq_missing_ms,
					double *nsites1_pop, double *nsites1_pop_outg,
				    double *nsites2_pop, double *nsites2_pop_outg,
					double *nsites3_pop, double *nsites3_pop_outg,
					long int niterdata, char *matrix_pol, int *r2i_ploidies,
                    char *matrix_pol_tcga,char *chr_name)
{
	int x=0;
	int y=0;
	int z,z1,z2,np,oo,ss,xx;
	long int zz;
	int sumnsam;
	int *initsq1;
	double mean_z1,mean_z2;

	double cv;
	double delta,der;
	double sf;
	/*
	double an;
	double bn;
	int j;
	*/
	/*long int freqo[4];*/
	/* char ancestral[1]; */
	long int fabsmatsize;
    char nt1[1];
    char nt2[1];
    
    nt1[0] = 0;
    nt2[0] = 0;
	
/**< todo: WTF is that LUCA_CRs???? */
#if LUCA_CRs == 1
	double luca_cvo2, luca_cvo1, luca_cv0, luca_cv1;
	double luca_dero2, luca_dero1, luca_der0, luca_der1;

	double luca_mean_z1, luca_mean_z2;
	double luca_cso2, luca_cso1, luca_cs0, luca_cs1;
	double luca_cro2, luca_cro1, luca_cr0, luca_cr1;

	double luca_cvo2d, luca_cvo1d, luca_cv0d, luca_cv1d;
	double luca_dero2d, luca_dero1d, luca_der0d, luca_der1d;
	
	double luca_cso2d, luca_cso1d, luca_cs0d, luca_cs1d;
	double luca_cro2d, luca_cro1d, luca_cr0d, luca_cr1d;
#endif
	
	initsq1 = (int *)calloc(npops,sizeof(int));
	sumnsam = 0;
	for(x=0;x<npops;x++) {
		initsq1[x] = sumnsam;
		sumnsam += nsam[x];
	}

	
	if(output == 0 || output == 10) {
		if(file_out == 0) file_out = stdout;
		fprintf(file_out,"\nInput file: %s",file_input);
		if(ploidy[0] == '2') fprintf(file_out,"\nEach sample represents two chromosomes ('.1' and '.2') IUPAC symbols must be used. Haplotype based analysis NOT SHOWN.");
		/*fprintf(file_out,"\nOutput file: %s\n",file_output);*/
		if(H1frq) {
			fprintf(file_out,"\nH1 model frequency file: %s",file_H1f);
			if(H0frq) {
				fprintf(file_out,"\nH0 model frequency file: %s",file_H0f);
			}
		}
		/**/
		if(npriors) {
            if(formatfile != 3) {
                fprintf(file_out,"\nPriors:\n");
                for(x=0;x<npriors;x++) {
                    fprintf(file_out,"prior%02d: %f\t",x,vector_priors[x]);
                }
            }
            if(formatfile == 3) { /*tfasta*/
                if(npriors == 2) {
                    fprintf(file_out,"scaffold_name: %s\t",chr_name);
                    fprintf(file_out,"start_window: %.0f\t",vector_priors[0]);
                    fprintf(file_out,"end_window: %.0f\t",vector_priors[1]);
                }
                else {
                    for(x=0;x<npriors;x++) {
                        fprintf(file_out,"prior%02d: %f\t",x,vector_priors[x]);
                    }
                }
            }
		}
		/**/
		if(gfffiles == 1) {
			fprintf(file_out,"\nGFF file: %s",file_GFF);
			fprintf(file_out,"\nPositions studied: %s",subset_positions);
			if(strcmp(subset_positions,"synonymous")==0 || strcmp(subset_positions,"nonsynonymous")==0 || strcmp(subset_positions,"silent")==0)
				fprintf(file_out,"\nGenetic Code: %s",code_name);
			if(strcmp(code_name,"Other")==0) {
				fputs("\nGenetic Code assigned to each triplet:\n",file_out);
				for(x=0;x<64;x++) {
					fprintf(file_out,"[%c%c%c] = %c\n",tripletsU[x][0],tripletsU[x][1],tripletsU[x][2],genetic_code[x]);
				}
			}
		}
		fprintf(file_out,"\n\nInclude Unknown positions [1/0]: %d",include_unknown);
		fprintf(file_out,"\nLenght of the Total alignment (including gaps): %ld",length);
		fprintf(file_out,"\nLenght of the Total alignment (excluding fixed gaps (or missing outgroup, if considered) but counting all positions as 1, eg. Syn positions < 1 count as 1): %ld",length_al_real);
		fprintf(file_out,"\nmultiple hits: %ld",statistics[0].nmhits);
		fprintf(file_out,"\nLenght of the Selected alignment (including unknown bp if selected, excluding bp with no outgroup, if considerd): %.2f",length_al);
		fprintf(file_out,"\nNumber of Variant (biallelic) sites (including unknown bp if selected and excluding codons with more than two variants, if considered, excluding bp with no outgroup, if considered): %ld",length_seg);
		
		if(svratio > -10000) fprintf(file_out,"\nRatio_S/V: %.3f\n",svratio);
		else fprintf(file_out,"\nRatio_S/V: NA\n");

		if(include_unknown ==1) fprintf(file_out,"Ratio_Missing (if outgroup, counting only positions where outgroup is present): %.6f",missratio);
        else fprintf(file_out,"Ratio_Missing (if outgroup, counting only positions where outgroup is present): NA");

		fprintf(file_out,"\n\nNumber of populations: %d",npops-!outgroup_presence);
		if(ploidy[0] == '2') fprintf(file_out,"\nNumber of samples (duplicated) for each population:");

		fprintf(file_out,"\nNumber of samples for each population:");
		for(x=0;x<npops-!outgroup_presence;x++) {
			fprintf(file_out,"\nnsam[%d]: %d",x,nsam[x]);
		}
		/*if(missratio > 0.) {*/
			/*
			fprintf(file_out,"\nEffective Number of samples for each population:");
			for(x=0;x<npops-!outgroup_presence;x++) {
				fprintf(file_out,"\nEff_nsam[%d]: %d",x,(int)ceil(statistics[0].length[x]/(double)length_al- 1E-5));
			}
			*/
			fprintf(file_out,"\n\nEffective length for each population (with at least one sequence, two or three sequences per pop) and excluding or including outgroup (if defined):");
			for(x=0;x<npops-1;x++) {
                /*if(outgroup_presence==1) {fprintf(file_out,"\nEff_length1_pop[%d]: %.2f\tEff_length2_pop[%d]: %.2f\tEff_length3_pop[%d]: %.2f\tEff_length1_pop_outg[%d]: %.2f\tEff_length2_pop_outg[%d]: %.2f\tEff_length3_pop_outg[%d]: %.2f",x,(double)nsites1_pop[x],x,(double)nsites2_pop[x],x,(double)nsites3_pop[x],x,(double)nsites1_pop_outg[x],x,(double)nsites2_pop_outg[x],x,(double)nsites3_pop_outg[x]);}*/
                if(outgroup_presence + force_outgroup) {fprintf(file_out,"\nEff_length1_pop_outg[%d]: %.2f\tEff_length2_pop_outg[%d]: %.2f\tEff_length3_pop_outg[%d]: %.2f",x,(double)nsites1_pop_outg[x],x,(double)nsites2_pop_outg[x],x,(double)nsites3_pop_outg[x]);}
				if(outgroup_presence + force_outgroup==0) {fprintf(file_out,"\nEff_length1_pop[%d]: %.2f\tEff_length2_pop[%d]: %.2f\tEff_length3_pop[%d]: %.2f",x,(double)nsites1_pop[x],x,(double)nsites2_pop[x],x,(double)nsites3_pop[x]);}
			}
		/*}*/
				
		fprintf(file_out,"\n\nSTATISTICS:\n");
		
		fprintf(file_out,"\nEstimates of variability for each population (an and bn for the variant positions):\n");
		if(outgroup_presence + force_outgroup) {
			np = npops-1;
			for(x=0;x<np;x++) {
				if(nsam[x] > 1) {
					/*
					fprintf(file_out,"S[%d]: %d\t",x,(int)statistics[0].S[x]);
					fprintf(file_out,"Thetaw(Wat)[%d]: %f\t",x,statistics[0].thetaS[x]);
					fprintf(file_out,"Thetaw(Taj)[%d]: %f\t",x,statistics[0].thetaT[x]);
					*/
					fprintf(file_out,"S[%d]: %d\t",x,(int)statistics[0].So[x]);
					fprintf(file_out,"Theta(Wat)[%d]: %f\t",x,statistics[0].thetaSo[x]);
					fprintf(file_out,"Theta(Taj)[%d]: %f\t",x,statistics[0].thetaTo[x]);
					fprintf(file_out,"Theta(Fu&Li)[%d]: %f\t",x,statistics[0].thetaFL[x]);
					fprintf(file_out,"Theta(Fay&Wu)[%d]: %f\t",x,statistics[0].thetaFW[x]);
					fprintf(file_out,"Theta(Zeng)[%d]: %f\t",x,statistics[0].thetaL[x]);
					fprintf(file_out,"Theta(Achaz,Wat)[%d]: %f\t",x,statistics[0].thetaSA[x]);
					fprintf(file_out,"Theta(Achaz,Taj)[%d]: %f\t",x,statistics[0].thetaTA[x]);
					fprintf(file_out,"Divergence[%d]: %f\t",x,statistics[0].K[x]);
					/*if(statistics[0].S[x]>0) {*/
                        /*fprintf(file_out,"an_x[%d]: %f\t",x,statistics[0].anx[x]);
                        fprintf(file_out,"bn_x[%d]: %f\t",x,statistics[0].bnx[x]);*/
                        fprintf(file_out,"an_xo[%d]: %f\t",x,statistics[0].anxo[x]);
                        fprintf(file_out,"bn_xo[%d]: %f\t",x,statistics[0].bnxo[x]);
					/*}
					else {
                        fprintf(file_out,"an_x[%d]: NA\t",x);
                        fprintf(file_out,"bn_x[%d]: NA\t",x);
                        fprintf(file_out,"an_xo[%d]: NA\t",x);
                        fprintf(file_out,"bn_xo[%d]: NA\t",x);
					}*/
					fprintf(file_out,"\n");
				}
				else {
					fprintf(file_out,"S[%d]: NA\t",x);
					fprintf(file_out,"Theta(Wat)[%d]: NA\t",x);
					fprintf(file_out,"Theta(Taj)[%d]: NA\t",x);
					fprintf(file_out,"Theta(Fu&Li)[%d]: NA\t",x);
					fprintf(file_out,"Theta(Fay&Wu)[%d]: NA\t",x);
					fprintf(file_out,"Theta(Zeng)[%d]: NA\t",x);
					fprintf(file_out,"Theta(Achaz,Wat)[%d]: NA\t",x);
					fprintf(file_out,"Theta(Achaz,Taj)[%d]: NA\t",x);
					fprintf(file_out,"Divergence[%d]: %f\t",x,statistics[0].K[x]);
                    /*fprintf(file_out,"an_x[%d]: NA\t",x);
                    fprintf(file_out,"bn_x[%d]: NA\t",x);*/
                    fprintf(file_out,"an_xo[%d]: NA\t",x);
                    fprintf(file_out,"bn_xo[%d]: NA\t",x);
					fprintf(file_out,"\n");
				}
			}
		}
		else {
			np = npops-1;
			for(x=0;x<np;x++) {
				if(nsam[x] > 1) {
					fprintf(file_out,"S[%d]: %d\t",x,(int)statistics[0].S[x]);
					fprintf(file_out,"Theta(Wat)[%d]: %f\t",x,statistics[0].thetaS[x]);
					fprintf(file_out,"Theta(Taj)[%d]: %f\t",x,statistics[0].thetaT[x]);
					fprintf(file_out,"Theta(Fu&Li)[%d]: %f\t",x,statistics[0].thetaFL[x]);
					fprintf(file_out,"Theta(Achaz,Wat)[%d]: %f\t",x,statistics[0].thetaSA[x]);
					fprintf(file_out,"Theta(Achaz,Taj)[%d]: %f\t",x,statistics[0].thetaTA[x]);
					/*if(statistics[0].S[x]>0) {*/
                        fprintf(file_out,"an_x[%d]: %f\t",x,statistics[0].anx[x]/*(double)statistics[0].S[x]/(double)statistics[0].thetaS[x]*/);
                        fprintf(file_out,"bn_x[%d]: %f\t",x,statistics[0].bnx[x]);
                        /*fprintf(file_out,"an_xo[%d]: %f\t",x,statistics[0].anxo[x]);
                        fprintf(file_out,"bn_xo[%d]: %f\t",x,statistics[0].bnxo[x]);*/
 					/*}
					else {
                        fprintf(file_out,"an_x[%d]: NA\t",x);
                        fprintf(file_out,"bn_x[%d]: NA\t",x);
                        fprintf(file_out,"an_xo[%d]: NA\t",x);
                        fprintf(file_out,"bn_xo[%d]: NA\t",x);
					}*/
					fprintf(file_out,"\n");
				}
				else {
					fprintf(file_out,"S[%d]: NA\t",x);
					fprintf(file_out,"Theta(Wat)[%d]: NA\t",x);
					fprintf(file_out,"Theta(Taj)[%d]: NA\t",x);
					fprintf(file_out,"Theta(Fu&Li)[%d]: NA\t",x);
					fprintf(file_out,"Theta(Achaz,Wat)[%d]: NA\t",x);
					fprintf(file_out,"Theta(Achaz,Taj)[%d]: NA\t",x);
                    fprintf(file_out,"an_x[%d]: NA\t",x);
                    fprintf(file_out,"bn_x[%d]: NA\t",x);
                    /*fprintf(file_out,"an_xo[%d]: NA\t",x);
                    fprintf(file_out,"bn_xo[%d]: NA\t",x);*/
                    /*fprintf(file_out,"theta_square[%d]: NA\n",x);*/
					fprintf(file_out,"\n");
				}
			}
		}
		if(formatfile != 2) {
			fprintf(file_out,"\nEstimates of NUCLEOTIDE variability for each population (if missing, corrected by the averaged effective positions):\n");
			if(outgroup_presence == 1) {
				np = npops-1;
				for(x=0;x<np;x++) {
					if(nsam[x] > 1) {
						/*
						fprintf(file_out,"S[%d]: %d\t",x,(int)statistics[0].S[x]);
						fprintf(file_out,"Thetaw/nt(Wat)[%d]: %f\t",x,statistics[0].thetaS[x]/(double)nsites2_pop[x]);
						fprintf(file_out,"Thetaw/nt(Taj)[%d]: %f\t",x,statistics[0].thetaT[x]/(double)nsites2_pop[x]);
						*/
						fprintf(file_out,"S[%d]: %d\t",x,(int)statistics[0].So[x]);
						if(nsites2_pop_outg[x] > 0) fprintf(file_out,"Theta/nt(Wat)[%d]: %f\t",x,statistics[0].thetaSo[x]/(double)nsites2_pop_outg[x]);
                        else fprintf(file_out,"Theta/nt(Wat)[%d]: NA\t",x);
						if(nsites2_pop_outg[x] > 0) fprintf(file_out,"Theta/nt(Taj)[%d]: %f\t",x,statistics[0].thetaTo[x]/(double)nsites2_pop_outg[x]);
                        else fprintf(file_out,"Theta/nt(Taj)[%d]: NA\t",x);
						if(nsites2_pop_outg[x] > 0) fprintf(file_out,"Theta/nt(Fu&Li)[%d]: %f\t",x,statistics[0].thetaFL[x]/(double)nsites2_pop_outg[x]);
                        else fprintf(file_out,"Theta/nt(Fu&Li)[%d]: NA\t",x);
						if(nsites2_pop_outg[x] > 0) fprintf(file_out,"Theta/nt(Fay&Wu)[%d]: %f\t",x,statistics[0].thetaFW[x]/(double)nsites2_pop_outg[x]);
                        else fprintf(file_out,"Theta/nt(Fay&Wu)[%d]: NA\t",x);
						if(nsites2_pop_outg[x] > 0) fprintf(file_out,"Theta/nt(Zeng)[%d]: %f\t",x,statistics[0].thetaL[x]/(double)nsites2_pop_outg[x]);
                        else fprintf(file_out,"Theta/nt(Zeng)[%d]: NA\t",x);
						if(nsites3_pop_outg[x] > 0) fprintf(file_out,"Theta/nt(Achaz,Wat)[%d]: %f\t",x,statistics[0].thetaSA[x]/(double)nsites3_pop_outg[x]);
                        else fprintf(file_out,"Theta/nt(Achaz,Wat)[%d]: NA\t",x);
						if(nsites3_pop_outg[x] > 0) fprintf(file_out,"Theta/nt(Achaz,Taj)[%d]: %f\t",x,statistics[0].thetaTA[x]/(double)nsites3_pop_outg[x]);
                        else fprintf(file_out,"Theta/nt(Achaz,Taj)[%d]: NA\t",x);
						if(nsites1_pop_outg[x] > 0) fprintf(file_out,"Divergence[%d]: %f\t",x,statistics[0].K[x]/(double)nsites1_pop_outg[x]);
                        else fprintf(file_out,"Divergence/nt[%d]: NA\t",x);
						if(statistics[0].thetaTHKY[x] > -10000/* && missratio == 0.*/) fprintf(file_out,"Theta/nt(Taj)HKY[%d]: %f\t",x,statistics[0].thetaTHKY[x]);
						else fprintf(file_out,"Theta/nt(Taj)HKY[%d]: NA\t",x);
						if(statistics[0].KHKY[x] > -10000/* && missratio == 0.*/) fprintf(file_out,"Divergence/nt_HKY[%d]: %f\n",x,statistics[0].KHKY[x]);
						else fprintf(file_out,"Divergence/nt_HKY[%d]: NA\n",x);
					}
					else {
						fprintf(file_out,"S[%d]: NA\t",x);
						fprintf(file_out,"Theta/nt(Wat)[%d]: NA\t",x);
						fprintf(file_out,"Theta/nt(Taj)[%d]: NA\t",x);
						fprintf(file_out,"Theta/nt(Fu&Li)[%d]: NA\t",x);
						fprintf(file_out,"Theta/nt(Fay&Wu)[%d]: NA\t",x);
						fprintf(file_out,"Theta/nt(Zeng)[%d]: NA\t",x);
						fprintf(file_out,"Theta/nt(Achaz,Wat)[%d]: NA\t",x);
						fprintf(file_out,"Theta/nt(Achaz,Taj)[%d]: NA\t",x);
						fprintf(file_out,"Divergence/nt[%d]: NA\t",x);
						fprintf(file_out,"Theta/nt(Taj)HKY[%d]: NA\t",x);
						fprintf(file_out,"Divergence/nt_HKY[%d]: NA\n",x);
					}
				}
			}
			else {
				np = npops-1;
				for(x=0;x<np;x++) {
					if(nsam[x] > 1) {
						fprintf(file_out,"S[%d]: %d\t",x,(int)statistics[0].S[x]);
						if(nsites2_pop[x] > 0) fprintf(file_out,"Theta/nt(Wat)[%d]: %f\t",x,statistics[0].thetaS[x]/(double)nsites2_pop[x]);
                        else fprintf(file_out,"Theta/nt(Wat)[%d]: NA\t",x);
						if(nsites2_pop[x] > 0) fprintf(file_out,"Theta/nt(Taj)[%d]: %f\t",x,statistics[0].thetaT[x]/(double)nsites2_pop[x]);
                        else fprintf(file_out,"Theta/nt(Taj)[%d]: NA\t",x);
						if(nsites2_pop[x] > 0) fprintf(file_out,"Theta/nt(Fu&Li)[%d]: %f\t",x,statistics[0].thetaFL[x]/(double)nsites2_pop[x]);
                        else fprintf(file_out,"Theta/nt(Fu&Li)[%d]: NA\t",x);
						if(nsites3_pop[x] > 0) fprintf(file_out,"Theta/nt(Achaz,Wat)[%d]: %f\t",x,statistics[0].thetaSA[x]/(double)nsites3_pop[x]);
                        else fprintf(file_out,"Theta/nt(Achaz,Wat)[%d]: NA\t",x);
						if(nsites3_pop[x] > 0) fprintf(file_out,"Theta/nt(Achaz,Taj)[%d]: %f\t",x,statistics[0].thetaTA[x]/(double)nsites3_pop[x]);
                        else fprintf(file_out,"Theta/nt(Achaz,Taj)[%d]: NA\n",x);
						if(statistics[0].thetaTHKY[x] > -10000/* && missratio == 0.*/) fprintf(file_out,"Theta/nt(Taj)HKY[%d]: %f\t",x,statistics[0].thetaTHKY[x]);
						else fprintf(file_out,"Theta/nt(Taj)HKY[%d]: NA\t",x);
						fprintf(file_out,"\n");
					}
					else {
						fprintf(file_out,"S[%d]: NA\t",x);
						fprintf(file_out,"Theta/nt(Wat)[%d]: NA\t",x);
						fprintf(file_out,"Theta/nt(Taj)[%d]: NA\t",x);
						fprintf(file_out,"Theta/nt(Fu&Li)[%d]: NA\t",x);
						fprintf(file_out,"Theta/nt(Achaz,Wat)[%d]: NA\t",x);
						fprintf(file_out,"Theta/nt(Achaz,Taj)[%d]: NA\t",x);
						fprintf(file_out,"Theta/nt(Taj)HKY[%d]: NA\n",x);
						fprintf(file_out,"\n");
					}
				}
			}		
		}
		if(include_unknown == 0 && ploidy[0] == '1') {
			fprintf(file_out,"\nHaplotype diversity and number of haplotypes for each population:\n");
			np = npops-1;
			for(x=0;x<np;x++) {
				if(nsam[x]>1) {
					fprintf(file_out,"HapW[%d]: %f\tnHap[%d]: %d",x,statistics[0].hapw[x],x,statistics[0].nhpop[x]);
					fprintf(file_out,"\n");
				}
				else {
					fprintf(file_out,"HapW[%d]: NA\tnHap[%d]: 1",x,x);
					fprintf(file_out,"\n");
				}
			}
		}
		fprintf(file_out,"\nNeutrality tests for each population:\n");
		if(outgroup_presence == 1) {
			np = npops-1;
			for(x=0;x<np;x++) {
				if(statistics[0].Dtaj[x] > -10000)
					fprintf(file_out,"Tajima D[%d]: %f\t",x,statistics[0].Dtaj[x]);
				else fprintf(file_out,"Tajima D[%d]: NA\t",x);
				if(statistics[0].Dfl[x] > -10000)
					fprintf(file_out,"Fu&Li D[%d]: %f\t",x,statistics[0].Dfl[x]);
				else  fprintf(file_out,"Fu&Li D[%d]: NA\t",x);
				if(statistics[0].Ffl[x] > -10000)
					fprintf(file_out,"Fu&Li F[%d]: %f\t",x,statistics[0].Ffl[x]);
				else fprintf(file_out,"Fu&Li F[%d]: NA\t",x);
				if(statistics[0].Hnfw[x] > -10000)
					fprintf(file_out,"Fay&Wu norm H[%d]: %f\t",x,statistics[0].Hnfw[x]);
				else fprintf(file_out,"Fay&Wu norm H[%d]: NA\t",x);
				/*if(statistics[0].thetaT[x]-statistics[0].thetaFW[x] != -10000)
					fprintf(file_out,"Fay&WuH[%d]:\t%f\t",x,statistics[0].thetaT[x]-statistics[0].thetaFW[x]);
				else fprintf(file_out,"Fay&WuH[%d]: NA\t",x);*/
				if(statistics[0].Ez[x] > -10000)
					fprintf(file_out,"Zeng E[%d]: %f\t",x,statistics[0].Ez[x]);
				else fprintf(file_out,"Zeng E[%d]: NA\t",x);
                if(statistics[0].Yach[x] > -10000)
                    fprintf(file_out,"Achaz Y[%d]: %f\t",x,statistics[0].Yach[x]);
                else fprintf(file_out,"Achaz Y[%d]: NA\t",x);
/**/
                if(statistics[0].FH[x] > -10000)
                    fprintf(file_out,"Ferretti L[%d]: %f\t",x,statistics[0].FH[x]);
                else fprintf(file_out,"Ferretti L[%d]: NA\t",x);
/**/
				if(ploidy[0] == '1') {
					/*
					if(statistics[0].R2[x] > -10000)
						fprintf(file_out,"R2[%d]: %f\t",x,statistics[0].R2[x]);
					else fprintf(file_out,"R2[%d]: NA\t",x);
					*/
					if(include_unknown == 0) {
						if(statistics[0].Fs[x] > -10000 && missratio < 1e-6)
							fprintf(file_out,"Fs[%d]: %f\t",x,statistics[0].Fs[x]);
						else fprintf(file_out,"Fs[%d]: NA\t",x);
					}
				}/*
				else {
					if(statistics[0].R2[x] > -10000)
						fprintf(file_out,"R2d[%d]: %f\t",x,statistics[0].R2[x]);
					else fprintf(file_out,"R2d[%d]: NA\t",x);	
				}
				*/
				for(xx=0;xx<r2i_ploidies[0];xx++) {
					if(statistics[0].R2p[xx][x] > -10000) {
						if(r2i_ploidies[0]==1 && ploidy[0] == '1') fprintf(file_out,"R2[%d]: %f\t",x,statistics[0].R2p[xx][x]);
						if(r2i_ploidies[0]==1 && ploidy[0] == '2') fprintf(file_out,"R2d[%d]: %f\t",x,statistics[0].R2p[xx][x]);
						if(r2i_ploidies[0]!=1) fprintf(file_out,"R2p%d[%d]: %f\t",r2i_ploidies[xx+1],x,statistics[0].R2p[xx][x]);
					} else {
						if(r2i_ploidies[0]==1 && ploidy[0] == '1') fprintf(file_out,"R2[%d]: NA\t",x);
						if(r2i_ploidies[0]==1 && ploidy[0] == '2') fprintf(file_out,"R2d[%d]: NA\t",x);
						if(r2i_ploidies[0]!=1) fprintf(file_out,"R2p%d[%d]: NA\t",r2i_ploidies[xx+1],x);
						
					}
				}
				
				/*fprintf(file_out,"Fay&WuH[%d]:\t%f\t",x,statistics[0].thetaT[x]-statistics[0].thetaFW[x]);*/
				/*fprintf(file_out,"Fay&WuH2[%d]:\t%f\t",x,2.*(statistics[0].thetaT[x]-statistics[0].thetaL[x]));*/
				fprintf(file_out,"\n");
			}
		}
		else {
			np = npops-1;
			for(x=0;x<np;x++) {
				if(statistics[0].Dtaj[x] > -10000)
					fprintf(file_out,"Tajima D[%d]: %f\t",x,statistics[0].Dtaj[x]);
				else fprintf(file_out,"Tajima D[%d]: NA\t",x);
				if(statistics[0].Dfl[x] > -10000)
					fprintf(file_out,"Fu&Li D*[%d]: %f\t",x,statistics[0].Dfl[x]);
				else  fprintf(file_out,"Fu&Li D*[%d]: NA\t",x);
				if(statistics[0].Ffl[x] > -10000)
					fprintf(file_out,"Fu&Li F*[%d]: %f\t",x,statistics[0].Ffl[x]);
				else fprintf(file_out,"Fu&Li F*[%d]: NA\t",x);
				if(statistics[0].Yach[x] > -10000)
					fprintf(file_out,"Achaz Y*[%d]: %f\t",x,statistics[0].Yach[x]);
				else fprintf(file_out,"Achaz Y*[%d]: NA\t",x);
				if(ploidy[0] == '1') {
					/*
					if(statistics[0].R2[x] > -10000)
						fprintf(file_out,"R2[%d]: %f\t",x,statistics[0].R2[x]);
					else fprintf(file_out,"R2[%d]: NA\t",x);
					*/
					if(include_unknown == 0) {
						if(statistics[0].Fs[x] > -10000 && missratio < 1e-6)
							fprintf(file_out,"Fs[%d]: %f\t",x,statistics[0].Fs[x]);
						else fprintf(file_out,"Fs[%d]: NA\t",x);
					}
				}/*
				else {
					if(statistics[0].R2[x] > -10000)
						fprintf(file_out,"R2d[%d]: %f\t",x,statistics[0].R2[x]);
					else fprintf(file_out,"R2d[%d]: NA\t",x);					
				}*/
				for(xx=0;xx<r2i_ploidies[0];xx++) {
					if(statistics[0].R2p[xx][x] > -10000) {
						if(r2i_ploidies[0]==1 && ploidy[0] == '1') fprintf(file_out,"R2[%d]: %f\t",x,statistics[0].R2p[xx][x]);
						if(r2i_ploidies[0]==1 && ploidy[0] == '2') fprintf(file_out,"R2d[%d]: %f\t",x,statistics[0].R2p[xx][x]);
						if(r2i_ploidies[0]!=1) fprintf(file_out,"R2p%d[%d]: %f\t",r2i_ploidies[xx+1],x,statistics[0].R2p[xx][x]);
					} else {
						if(r2i_ploidies[0]==1 && ploidy[0] == '1') fprintf(file_out,"R2[%d]: NA\t",x);
						if(r2i_ploidies[0]==1 && ploidy[0] == '2') fprintf(file_out,"R2d[%d]: NA\t",x);
						if(r2i_ploidies[0]!=1) fprintf(file_out,"R2p%d[%d]: NA\t",r2i_ploidies[xx+1],x);
						
					}
				}

				fprintf(file_out,"\n");
			}
		}
		
		if(H1frq/* && missratio == 0.*/) {
			if(outgroup_presence == 1) {		
				fprintf(file_out,"\nAlternative Expected Frequency Spectrum of variants for each population:");
				for(x=0;x<npops-1;x++) {
					fprintf(file_out,"\n");
					for(y=1;y<nsam[x];y++) {
						fprintf(file_out,"frH1[%d,%d]: %.5f\t",x,y,statistics[0].H1freq[x][y]);
					}
				}
				fprintf(file_out,"\nOptimal tests given SNM as null model: To[inf,inf], To[0,0], To[inf,0]");/*, ToQc[inf,inf], ToQw[inf,inf], ToLc[inf,inf].");*/
				for(x=0;x<npops-1;x++) {
					fprintf(file_out,"\n");
					if(statistics[0].To_ii[x] > -10000)
						fprintf(file_out,"To_ii[%d]: %f\t",x,statistics[0].To_ii[x]);
					else  fprintf(file_out,"To_ii[%d]: NA\t",x);
					if(statistics[0].To_00[x] > -10000)
						fprintf(file_out,"To_00[%d]: %f\t",x,statistics[0].To_00[x]);
					else  fprintf(file_out,"To_00[%d]: NA\t",x);
					if(statistics[0].To_i0[x] > -10000)
						fprintf(file_out,"To_i0[%d]: %f\t",x,statistics[0].To_i0[x]);
					else  fprintf(file_out,"To_i0[%d]: NA\t",x);
				#if TO_NEW
					if(statistics[0].To_Qc_ii[x] > -10000)
						fprintf(file_out,"To_Qc_ii[%d]: %f\t",x,statistics[0].To_Qc_ii[x]);
					else  fprintf(file_out,"To_Qc_ii[%d]: NA\t",x);
					if(statistics[0].To_Qw_ii[x] > -10000)
						fprintf(file_out,"To_Qw_ii[%d]: %f\t",x,statistics[0].To_Qw_ii[x]);
					else  fprintf(file_out,"To_Qw_ii[%d]: NA\t",x);
					if(statistics[0].To_Lc_ii[x] > -10000)
						fprintf(file_out,"To_Lc_ii[%d]: %f\t",x,statistics[0].To_Lc_ii[x]);
					else  fprintf(file_out,"To_Lc_ii[%d]: NA\t",x);
				#endif
				}
				fprintf(file_out,"\n");
				if(H0frq) {
					fprintf(file_out,"\nNULL Expected Frequency Spectrum of variants for each population:");
					for(x=0;x<npops-1;x++) {
						fprintf(file_out,"\n");
						for(y=1;y<nsam[x];y++) {
							fprintf(file_out,"frH1[%d,%d]: %.5f\t",x,y,statistics[0].H0freq[x][y]);
						}
					}
					fprintf(file_out,"\nOptimal tests given the null expected frequency Spectrum as null model: ToH0[inf,inf], ToH0[inf,0].");
					for(x=0;x<npops-1;x++) {
						fprintf(file_out,"\n");
						if(statistics[0].ToH0_ii[x] > -10000)
							fprintf(file_out,"To_H0_ii[%d]: %f\t",x,statistics[0].ToH0_ii[x]);
						else  fprintf(file_out,"To_H0_ii[%d]: NA\t",x);
						if(statistics[0].ToH0_00[x] > -10000)
							fprintf(file_out,"To_H0_i0[%d]: %f\t",x,statistics[0].ToH0_00[x]);
						else  fprintf(file_out,"To_H0_i0[%d]: NA\t",x);
					}
					fprintf(file_out,"\n");
				}
			}
		}

		fprintf(file_out,"\nVariants assigned to exclusive, fixed, polymorphic but fixed in rest of pops, and shared.\n");
		fprintf(file_out,"Ss[rest] are shared variants between populations but fixed within:\n");
		if(npops>2 && outgroup_presence == 1) {
			for(x=0;x<npops-1;x++) {
				fprintf(file_out,"Sx[%d]: %ld\t",x,statistics[0].Sanc[x*4+0]);
				fprintf(file_out,"Sf[%d]: %ld\t",x,statistics[0].Sanc[x*4+1]);
				fprintf(file_out,"Sxf[%d,rest]: %ld\t",x,statistics[0].Sanc[x*4+2]);
				fprintf(file_out,"Ss[%d]: %ld\n",x,statistics[0].Sanc[x*4+3]);
			}
			fprintf(file_out,"Sx[outg]: %ld\t",statistics[0].Sanc[(npops-1)*4+0]);
			fprintf(file_out,"Sf[outg]: %ld\t",statistics[0].Sanc[(npops-1)*4+1]);
			fprintf(file_out,"Ss[outg,rest]: %ld\t",statistics[0].Sanc[(npops-1)*4+2]);
			if(npops > 3) fprintf(file_out,"Ss[rest]:\t%ld\n",statistics[0].Sanc[(npops-1)*4+3]);
			else fprintf(file_out,"Ss[rest]:\tNA\n");
		}
		else {
			if(npops>2 && outgroup_presence == 0) {
				for(x=0;x<npops-1;x++) {
					fprintf(file_out,"Sx[%d]: %ld\t",x,statistics[0].Sanc[x*4+0]/* + statistics[0].Sanc[x*4+2]*/);
					fprintf(file_out,"Sf[%d]: %ld\t",x,statistics[0].Sanc[x*4+1]);
					fprintf(file_out,"Ss[%d]: %ld\n",x,statistics[0].Sanc[x*4+3]);
				}
				fprintf(file_out,"Ss[rest]: %ld\n",statistics[0].Sanc[(npops-1)*4+3]);
			}
			else {
				if(npops > 1) {
					x = 0;
					fprintf(file_out,"Sx[%d]: %ld\t",x,statistics[0].Sanc[x*4+0]);
					x = 1;
					fprintf(file_out,"Sx[%d]: %ld\t",x,statistics[0].Sanc[x*4+0]);
					x=0;
					fprintf(file_out,"Sf: %ld\t",statistics[0].Sanc[x*4+1]);
					fprintf(file_out,"Ss: %ld\n",statistics[0].Sanc[x*4+3]);
				}
				else
					fprintf(file_out,"Sx[%d]: %ld\n",x,statistics[0].Sanc[0*4+0]);
			}
		}
		/*if(ploidy[0] == '1' && include_unknown == 0) { */
			fprintf(file_out,"\nmismatch distribution statistics:\n");
			for(x=0;x<npops-1/*!outgroup_presence*/;x++) {
				if(nsam[x] > 1) {
                    if(statistics[0].mdsd[x] > -10000.) fprintf(file_out,"SDev[%d]: %f\t",x,statistics[0].mdsd[x]);
                    else fprintf(file_out,"SDev[%d]: NA\t",x);
                }
                else
                    fprintf(file_out,"SDev[%d]: NA\t",x);
                if(nsam[x] > 2) {
					if(statistics[0].mdg1[x] > -10000. && statistics[0].mdsd[x] > -10000.) fprintf(file_out,"Skewness[%d]: %f\t",x,statistics[0].mdg1[x]);
					else fprintf(file_out,"Skewness[%d]: NA\t",x);
				}
				else {
					fprintf(file_out,"Skewness[%d]: NA\t",x);
				}
				if(nsam[x] > 3) {
					if(statistics[0].mdg2[x] > -10000. && statistics[0].mdsd[x] > -10000.) fprintf(file_out,"Kurtosis[%d]: %f\t",x,statistics[0].mdg2[x]);
					else fprintf(file_out,"Kurtosis[%d]: NA\t",x);
				}
				else
					fprintf(file_out,"Kurtosis[%d]: NA\t",x);
					
				fprintf(file_out,"\n");
			}
		/*}*/
		if(npops-1/*!outgroup_presence*/ > 1) {						
			fprintf(file_out,"\nFst (calculated as 1-piw/pia):\n"); /*ANALYSIS using estimator from Hudson, Slatkin and Maddison (Genetics, 1992) and Gst' from Nei (1987)*/
			fprintf(file_out,"\nseed: %ld",nseed);
			if(npops-1/*!outgroup_presence*/ > 2) {
				fprintf(file_out,"\nFst(nucleotide) and FstH(haplotype)for the whole populations (except outgroup)");/* and Gst' */
				if(niter) fprintf(file_out," and P-values using %ld iterations with permutation test:\n",niter);
				else fprintf(file_out,"\n");
				if(statistics[0].fstALL > -10000) {
					if(niter && piter[0].niteriall > 0 && include_unknown == 0)
						fprintf(file_out,"Fst: %f\t P-value: %f\t",statistics[0].fstALL,(double)piter[0].iall/(double)piter[0].niteriall);
					else 
						fprintf(file_out,"Fst: %f\t P-value: NA\t",statistics[0].fstALL);
				}	
				else {
					if(niter && piter[0].niteriall > 0 && include_unknown == 0)
						fprintf(file_out,"Fst: NA\t P-value: NA\t");
					else
						fprintf(file_out,"Fst: NA\t P-value: NA\t");
				}
				if(ploidy[0] == '1' && include_unknown == 0) {
					if(statistics[0].fsthALL > -10000) {
						if(niter && piter[0].niterihall > 0 && include_unknown == 0)
							fprintf(file_out,"FstH: %f\t P-value: %f\t",statistics[0].fsthALL,(double)piter[0].ihall/(double)piter[0].niterihall);
						else
							fprintf(file_out,"FstH: %f\t P-value: NA\t",statistics[0].fsthALL);
					}
					else {
						if(niter && piter[0].niterihall > 0 && include_unknown == 0)
							fprintf(file_out,"FstH: NA\t P-value: NA\t");
						else
							fprintf(file_out,"FstH: NA\t P-value: NA\t");
					}
					/*
					if(statistics[0].GstALL > -10000) {
						if(niter && piter[0].niterighall > 0)
							fprintf(file_out,"Gst': %f\t P-value: %f\t",statistics[0].GstALL,(double)piter[0].ighall/(double)piter[0].niterighall);
						else
							fprintf(file_out,"Gst': %f\t P-value: NA\t",statistics[0].GstALL);
					}
					else {
						if(niter && piter[0].niterighall > 0)
							fprintf(file_out,"Gst': NA\t P-value: NA\t");
						else
							fprintf(file_out,"Gst': NA\t P-value: NA\t");
					}
					 */
				}
				fprintf(file_out,"\n");
			
				fprintf(file_out,"\nFst(nucleotide) and FstH(haplotype) of each population ACROSS ALL the rest (except outgroup)");
				if(niter) fprintf(file_out," and P-values using %ld iterations with permutation test:\n",niter);
				else fprintf(file_out,"\n");
				for(x=0;x<npops-1;x++) {
					if(statistics[0].fst1all[x] > -10000) {
						if(niter && piter[0].niteri1[x] > 0 && include_unknown == 0)
							fprintf(file_out,"Fst1[%d,rest]: %f\t P-value: %f\t",x,statistics[0].fst1all[x],(double)piter[0].i1[x]/(double)piter[0].niteri1[x]);
						else 
							fprintf(file_out,"Fst1[%d,rest]: %f\t P-value: NA\t",x,statistics[0].fst1all[x]);
					}	
					else {
						if(niter && piter[0].niteri1[x] > 0)  
							fprintf(file_out,"Fst1[%d,rest]: NA\t P-value: NA\t",x);
						else
							fprintf(file_out,"Fst1[%d,rest]: NA\t P-value: NA\t",x);
					}
					if(ploidy[0] == '1' && include_unknown == 0) {
						if(statistics[0].fsth1all[x] > -10000) {
							if(niter && piter[0].niterih1[x] > 0 && include_unknown == 0)
								fprintf(file_out,"Fst1H[%d,rest]: %f\t P-value: %f\t",x,statistics[0].fsth1all[x],(double)piter[0].ih1[x]/(double)piter[0].niterih1[x]);
							else
								fprintf(file_out,"Fst1H[%d,rest]: %f\t P-value: NA\t",x,statistics[0].fsth1all[x]);
						}
						else {
							if(niter && piter[0].niterih1[x] > 0)
								fprintf(file_out,"Fst1H[%d,rest]: NA\t P-value: NA\t",x);
							else
								fprintf(file_out,"Fst1H[%d,rest]: NA\t P-value: NA\t",x);
						}
					}
					fprintf(file_out,"\n");
				}
			}
			
			fprintf(file_out,"\nFst(nucleotide), FstH(haplotype)"); /* and Gst' BETWEEN populations*/
			if(niter) fprintf(file_out," and P-values using %ld iterations with permutation test:\n",niter);
			else fprintf(file_out,"\n");
			z=0;
			for(x=0;x<npops-1;x++) {
				for(y=x+1;y<npops-0;y++) {
					if(y==npops-1) {z++;continue;}
					if(statistics[0].fst[z] > -10000) {
						if(niter && piter[0].niteri[z] > 0 && include_unknown == 0)
							fprintf(file_out,"Fst[%d,%d]: %f\t P-value: %f\t",x,y,statistics[0].fst[z],(double)piter[0].i[z]/(double)piter[0].niteri[z]);
						else
							fprintf(file_out,"Fst[%d,%d]: %f\t P-value: NA\t",x,y,statistics[0].fst[z]);
					}
					else {
						if(niter && piter[0].niteri[z] > 0 && include_unknown == 0)
							fprintf(file_out,"Fst[%d,%d]: NA\t P-value: NA\t",x,y);
						else
							fprintf(file_out,"Fst[%d,%d]: NA\t P-value: NA\t",x,y);
					}
					fprintf(file_out,"\n");
					z++;
				}
			}
			fprintf(file_out,"\n");
			z=0;
			for(x=0;x<npops-1;x++) { /*NOTE THE HAPLOTYPE VALUES GOES FROM 0 to npops-1 BUT FREQUENCY GOES FROM 0 to npops!!!!!!*/
				for(y=x+1;y<npops-0;y++) {
					if(y==npops-1) {z++;continue;}
					if(ploidy[0] == '1' && include_unknown == 0) {
						if(statistics[0].fsth[z] > -10000) {
							if(niter && piter[0].niterih[z] > 0 && include_unknown == 0)
								fprintf(file_out,"FstH[%d,%d]: %f\t P-value: %f\t",x,y,statistics[0].fsth[z],(double)piter[0].ih[z]/(double)piter[0].niterih[z]);
							else
								fprintf(file_out,"FstH[%d,%d]: %f\t P-value: NA\t",x,y,statistics[0].fsth[z]);
						}
						else {
							if(niter && piter[0].niterih[z] > 0 && include_unknown == 0)
								fprintf(file_out,"FstH[%d,%d]: NA\t P-value: NA\t",x,y);
							else
								fprintf(file_out,"FstH[%d,%d]: NA\t P-value: NA\t",x,y);
						}
						/*
						if(statistics[0].Gst[z] > -10000) {
							if(niter && piter[0].niterigh[z] > 0)
								fprintf(file_out,"Gst'[%d,%d]: %f\t P-value: %f\t",x,y,statistics[0].Gst[z],(double)piter[0].igh[z]/(double)piter[0].niterigh[z]);
							else
								fprintf(file_out,"Gst'[%d,%d]: %f\t P-value: NA\t",x,y,statistics[0].Gst[z]);
						}
						else {
							if(niter && piter[0].niterigh[z] > 0)  
								fprintf(file_out,"Gst'[%d,%d]: NA\t P-value: NA\t",x,y);
							else
								fprintf(file_out,"Gst'[%d,%d]: NA\t P-value: NA\t",x,y);
						}
						 */
						fprintf(file_out,"\n");
					}
					z++;
				}
			}

			if(include_unknown == 0 && outgroup_presence==1) {
				fprintf(file_out,"\nFst corrected with HKY (Hasegawa, Kishino, and Yano. 1985): \n"); /*ANALYSIS using estimator from Hudson, Slatkin and Maddison (Genetics, 1992) */
				if(npops > 2) {
					z=0;
					for(x=0;x<npops-1;x++) {
						for(y=x+1;y<npops-0;y++) {
							if(y==npops-1) {z++;continue;}
							if(statistics[0].fstHKY[z] > -10000) {
								fprintf(file_out,"FstHKY[%d,%d]: %f\tPiWHKY[%d]: %f\tPiWHKY[%d]: %f\tPiAHKY[%d,%d]: %f\tPiTHKY[%d,%d]: %f\t",
										x,y,statistics[0].fstHKY[z],x,statistics[0].piwHKY[x],y,statistics[0].piwHKY[y],x,y,statistics[0].piaHKY[z],x,y,statistics[0].piTHKY[z]);
							}
							else {
								fprintf(file_out,"FstHKY[%d,%d]: NA\t",x,y);
								if(statistics[0].piwHKY[x] > -10000) {
									fprintf(file_out,"PiWHKY[%d]: %f\t",x,statistics[0].piwHKY[x]);
								}
								else fprintf(file_out,"PiWHKY[%d]: NA\t",x);
								if(statistics[0].piwHKY[y] > -10000) {
									fprintf(file_out,"PiWHKY[%d]: %f\t",y,statistics[0].piwHKY[y]);
								}
								else fprintf(file_out,"PiWHKY[%d]: NA\t",y);
								if(statistics[0].piaHKY[z] > -10000) {
									fprintf(file_out,"PiAHKY[%d,%d]: %f\tPiTHKY[%d,%d]: %f\t",x,y,statistics[0].piaHKY[z],x,y,statistics[0].piTHKY[z]);
								}
								else fprintf(file_out,"PiAHKY[%d,%d]: NA\tPiTHKY[%d,%d]: NA\t",x,y,x,y);
							}
							fprintf(file_out,"\n");
							z++;
						}
					}
				}
			}
		}

		if(include_unknown == 0) {
			fprintf(file_out,"\nNucleotide and haplotype diversity WITHIN populations:\n");
			for(x=0;x<npops-1/*!outgroup_presence*/;x++) {
				if(nsam[x]>1) {
					fprintf(file_out,"PiW[%d]: %f\t",x,statistics[0].piw[x]);
					if(ploidy[0] == '1' && include_unknown == 0) fprintf(file_out,"HapW[%d]: %f\t",x,statistics[0].hapw[x]);
					fprintf(file_out,"\n");
				}
				else {
					fprintf(file_out,"PiW[%d]: NA\t",x);
					if(ploidy[0] == '1' && include_unknown == 0) fprintf(file_out,"HapW[%d]: NA\t",x);
					fprintf(file_out,"\n");
				}
			}
			if(npops-1/*!outgroup_presence*/ > 1) {
				fprintf(file_out,"\nNucleotide and haplotype diversity BETWEEN populations and TOTAL:\n");
                z=0;
                for(x=0;x<npops-1;x++) {
                    for(y=x+1;y<npops-0;y++) {
                        if(y==npops-1) {z++;continue;}
                        if(statistics[0].piant[z] > -10000) {
                            fprintf(file_out,"PiA[%d,%d]: %f\tPiA/nt[%d,%d]: %f\tPiT[%d,%d]: %f\tPiT/nt[%d,%d]: %f\t",x,y,statistics[0].pia[z],x,y,statistics[0].piant[z],x,y,statistics[0].piT[z],x,y,statistics[0].piTnt[z]);
                        }
                        else {
                            fprintf(file_out,"PiA[%d,%d]: %f\tPiA/nt[%d,%d]: NA\tPiT[%d,%d]: %f\tPiT/nt[%d,%d]: NA\t",x,y,statistics[0].pia[z],x,y,x,y,statistics[0].piT[z],x,y);
                        }
                        fprintf(file_out,"\n");
                        z++;
                    }
                }
                fprintf(file_out,"\n");
				z=0;
				for(x=0;x<npops-1;x++) {
					for(y=x+1;y<npops-0;y++) {
						if(y==npops-1) {z++;continue;}
						if(ploidy[0] == '1' && include_unknown == 0) {
							fprintf(file_out,"HapA[%d,%d]: %f\tHapT[%d,%d]: %f\t",x,y,statistics[0].hapa[z],x,y,statistics[0].hapT[z]);
							fprintf(file_out,"\n");
						}
						z++;
					}
				}
			}
		}
		else {
			fprintf(file_out,"\nNucleotide and haplotype diversity WITHIN populations:\n");
			for(x=0;x<npops-1/*!outgroup_presence*/;x++) {
				if(nsam[x]>1) {
					if(statistics[0].length2[x] > 0.)
                        fprintf(file_out,"PiW[%d]: %f\tPiW/nt[%d]: %f\t",x,statistics[0].piw[x],x,statistics[0].piw[x]/(double)statistics[0].length2[x]);
                    else
                        fprintf(file_out,"PiW[%d]: %f\tPiW/nt[%d]: NA\t",x,statistics[0].piw[x],x);
					fprintf(file_out,"\n");
				}
				else {
					fprintf(file_out,"PiW[%d]: NA\tPiW/nt[%d]: NA\t",x,x);
					fprintf(file_out,"\n");
				}
			}
			if(npops-1/*!outgroup_presence*/ > 1) {
				fprintf(file_out,"\nNucleotide and haplotype diversity BETWEEN populations:\n");
				z=0;
				for(x=0;x<npops-1;x++) {
					for(y=x+1;y<npops-0;y++) {
						if(y==npops-1) {z++;continue;}
                        if(outgroup_presence+force_outgroup) {
                            if(statistics[0].piant[z] > -10000) {
                                fprintf(file_out,"PiA[%d,%d]: %f\tPiA/nt[%d,%d]: %f\tPiT[%d,%d]: %f\tPiT/nt[%d,%d]: %f\tlen[%d,%d]:",x,y,statistics[0].pia[z],x,y,statistics[0].piant[z],x,y,statistics[0].piT[z],x,y,statistics[0].piTnt[z],x,y);
                                if(outgroup_presence+force_outgroup)
                                    fprintf(file_out," %f\t",statistics[0].lengthamng_outg[x][y]);
                                else
                                    fprintf(file_out," %f\t",statistics[0].lengthamng[x][y]);
                            }
                            else {
                                fprintf(file_out,"PiA[%d,%d]: %f\tPiA/nt[%d,%d]: NA\tPiT[%d,%d]: %f\tPiT/nt[%d,%d]: NA\tlen[%d,%d]:",x,y,statistics[0].pia[z],x,y,x,y,statistics[0].piT[z],x,y,x,y);
                                if(outgroup_presence+force_outgroup)
                                    fprintf(file_out," %f\t",statistics[0].lengthamng_outg[x][y]);
                                else
                                    fprintf(file_out," %f\t",statistics[0].lengthamng[x][y]);
                            }
                        }
                        else {
                            if(statistics[0].piant[z] > -10000) {
                                fprintf(file_out,"PiA[%d,%d]: %f\tPiA/nt[%d,%d]: %f\tPiT[%d,%d]: %f\tPiT/nt[%d,%d]: %f\tlen[%d,%d]:",x,y,statistics[0].pia[z],x,y,statistics[0].piant[z],x,y,statistics[0].piT[z],x,y,statistics[0].piTnt[z],x,y);
                                if(outgroup_presence+force_outgroup)
                                    fprintf(file_out," %f\t",statistics[0].lengthamng_outg[x][y]);
                                else
                                    fprintf(file_out," %f\t",statistics[0].lengthamng[x][y]);
                            }
                            else {
                                fprintf(file_out,"PiA[%d,%d]: %f\tPiA/nt[%d,%d]: NA\tPiT[%d,%d]: %f\tPiT/nt[%d,%d]: NA\tlen[%d,%d]:",x,y,statistics[0].pia[z],x,y,x,y,statistics[0].piT[z],x,y,x,y);
                                if(outgroup_presence+force_outgroup)
                                    fprintf(file_out," %f\t",statistics[0].lengthamng_outg[x][y]);
                                else
                                    fprintf(file_out," %f\t",statistics[0].lengthamng[x][y]);
                            }
                        }
						fprintf(file_out,"\n");
						z++;
					}
				}
				fprintf(file_out,"\n");
			}
		}
		
		fprintf(file_out,"\nFrequency of variants for each population:");
		if(npops-!outgroup_presence > 1) {		
			for(x=0;x<npops-1/*!outgroup_presence*/;x++) {
				fprintf(file_out,"\n");
				if(outgroup_presence == 1) {
					for(y=1;y<nsam[x];y++) {
						fprintf(file_out,"fr[%d,%d]: %ld\t",x,y,statistics[0].freq[x][y]);
					}
				}
				else {
					for(y=1;y<=floor(nsam[x]/2);y++) {
						fprintf(file_out,"fr[%d,%d]: %ld\t",x,y,statistics[0].freq[x][y]);
					}
				}
			}
		}
		else {
			fprintf(file_out,"\n");
			for(y=1;y<=floor(nsam[0]/2);y++) {
				fprintf(file_out,"fr[0,%d]: %ld\t",y,statistics[0].freq[0][y]);
			}
		}
		fprintf(file_out,"\n");
		if(ploidy[0] == '1' && include_unknown == 0) {
			fprintf(file_out,"\nFrequency of each haplotype in the populations:");
			for(x=0;x<npops-1/*!outgroup_presence*/;x++) {
				fprintf(file_out,"\n");
				for(y=0;y<statistics[0].nh;y++) { ///******* problem!: in case force an outgroup we may have one more haplotype that is not real *********/
					fprintf(file_out,"frH[%d,hap%02d]: %ld\t",x,y,statistics[0].freqh[x][y]);
				}
			}
		}
		fprintf(file_out,"\n\nPositions of each variant (negative indicates it contains any missing 'N' values):");
		if(npops-!outgroup_presence>1) {
			for(x=0;x<npops-1;x++) {
				fprintf(file_out,"\nSx[%d]: ",x);				
				y=0;
				zz = sites_matrix[y*4*npops+4*x+0];
				while(zz!=0) {
					if(y >= length_seg) break;
                    if(formatfile == 3) {
                        if(zz<0) fprintf(file_out,"%ld ",zz-(long int)vector_priors[0]+1);
                        else fprintf(file_out,"%ld ",zz+(long int)vector_priors[0]-1);
                    } else fprintf(file_out,"%ld ",zz);
					y++;
					zz = sites_matrix[y*4*npops+4*x+0];
				}	
				if(outgroup_presence == 0) {
					y=0;
					zz = sites_matrix[y*4*npops+4*x+2];
					while(zz!=0) {
						if(y >= length_seg) break;
                        if(formatfile == 3) {
                            if(zz<0) fprintf(file_out,"%ld ",zz-(long int)vector_priors[0]+1);
                            else fprintf(file_out,"%ld ",zz+(long int)(long int)vector_priors[0]-1);
                        } else fprintf(file_out,"%ld ",zz);
						y++;
						zz = sites_matrix[y*4*npops+4*x+2];
					}
				}
				fprintf(file_out,"\nSf[%d]: ",x);
				y=0;
				zz = sites_matrix[y*4*npops+4*x+1];
				while(zz!=0) {
					if(y >= length_seg) break;
                    if(formatfile == 3) {
                        if(zz<0) fprintf(file_out,"%ld ",zz-(long int)vector_priors[0]+1);
                        else fprintf(file_out,"%ld ",zz+(long int)vector_priors[0]-1);
                    } else fprintf(file_out,"%ld ",zz);
					y++;
					zz = sites_matrix[y*4*npops+4*x+1];
				}
				if(outgroup_presence == 1) {
					fprintf(file_out,"\nSxf[%d]: ",x);
					y=0;
					zz = sites_matrix[y*4*npops+4*x+2];
					while(zz!=0) {
						if(y >= length_seg) break;
                        if(formatfile == 3) {
                            if(zz<0) fprintf(file_out,"%ld ",zz-(long int)vector_priors[0]+1);
                            else fprintf(file_out,"%ld ",zz+(long int)vector_priors[0]-1);
                        } else fprintf(file_out,"%ld ",zz);
						y++;
						zz = sites_matrix[y*4*npops+4*x+2];
					}
				}
				fprintf(file_out,"\nSs[%d]: ",x);
				y=0;
				zz = sites_matrix[y*4*npops+4*x+3];
				while(zz!=0) {
					if(y >= length_seg) break;
                    if(formatfile == 3) {
                        if(zz<0) fprintf(file_out,"%ld ",zz-(long int)vector_priors[0]+1);
                        else fprintf(file_out,"%ld ",zz+(long int)vector_priors[0]-1);
                    } else fprintf(file_out,"%ld ",zz);
					y++;
					zz = sites_matrix[y*4*npops+4*x+3];
				}
			}
			x = npops-1;
			fprintf(file_out,"\nSs[rest]: ");				
			y=0;
			zz = sites_matrix[y*4*npops+4*x+3];
			while(zz!=0) {
				if(y >= length_seg) break;
                if(formatfile == 3) {
                    if(zz<0) fprintf(file_out,"%ld ",zz-(long int)vector_priors[0]+1);
                    else fprintf(file_out,"%ld ",zz+(long int)vector_priors[0]-1);
                } else fprintf(file_out,"%ld ",zz);
				y++;
				zz = sites_matrix[y*4*npops+4*x+3];
			}				
			if(outgroup_presence == 1) {
				x = npops-1;
				fprintf(file_out,"\nSx[outg]: ");				
				y=0;
				zz = sites_matrix[y*4*npops+4*x+0];
				while(zz!=0) {
					if(y >= length_seg) break;
                    if(formatfile == 3) {
                        if(zz<0) fprintf(file_out,"%ld ",zz-(long int)vector_priors[0]+1);
                        else fprintf(file_out,"%ld ",zz+(long int)vector_priors[0]-1);
                    } else fprintf(file_out,"%ld ",zz);
					y++;
					zz = sites_matrix[y*4*npops+4*x+0];
				}				
				fprintf(file_out,"\nSf[outg]: ");
				y=0;
				zz = sites_matrix[y*4*npops+4*x+1];
				while(zz!=0) {
					if(y >= length_seg) break;
                    if(formatfile == 3) {
                        if(zz<0) fprintf(file_out,"%ld ",zz-(long int)vector_priors[0]+1);
                        else fprintf(file_out,"%ld ",zz+(long int)vector_priors[0]-1);
                    } else fprintf(file_out,"%ld ",zz);
					y++;
					zz = sites_matrix[y*4*npops+4*x+1];
				}
				fprintf(file_out,"\nSf[outg,rest]: ");
				y=0;
				zz = sites_matrix[y*4*npops+4*x+2];
				while(zz!=0) {
					if(y >= length_seg) break;
                    if(formatfile == 3) {
                        if(zz<0) fprintf(file_out,"%ld ",zz-(long int)vector_priors[0]+1);
                        else fprintf(file_out,"%ld ",zz+(long int)vector_priors[0]-1);
                    } else fprintf(file_out,"%ld ",zz);
					y++;
					zz = sites_matrix[y*4*npops+4*x+2];
				}
			}
		}
		else {
			x = 0;
			fprintf(file_out,"\nSx[%d]: ",x);				
			y=0;
			zz = sites_matrix[y*4*npops+4*x+0];
			while(zz!=0) {
				if(y >= length_seg) break;
                if(formatfile == 3) {
                    if(zz<0) fprintf(file_out,"%ld ",zz-(long int)vector_priors[0]+1);
                    else fprintf(file_out,"%ld ",zz+(long int)vector_priors[0]-1);
                } else fprintf(file_out,"%ld ",zz);
				y++;
				zz = sites_matrix[y*4*npops+4*x+0];
			}				
		}
		
		oo = 1;
		fprintf(file_out,"\n\nJoint frequency distribution for each variant and population (No included variants that are missing or polymorphic in the outgroup):");
		fprintf(file_out,"\n");
		for(x=0;x<npops-oo;x++) fprintf(file_out,"\tpop[%d]",x);
		for(zz=0;zz<length_seg;zz++) {
			ss=0; for(x=0;x<npops-oo;x++) {ss += nfd[x][zz];}
			if(ss) {
                if(formatfile ==3) {
                    if(matrix_pos[zz] < 0) fprintf(file_out,"\nSNP[%ld]",matrix_pos[zz]-(long int)vector_priors[0]+1);
                    else fprintf(file_out,"\nSNP[%ld]",matrix_pos[zz]+(long int)vector_priors[0]-1);
                }
                else fprintf(file_out,"\nSNP[%ld]",matrix_pos[zz]);
				for(x=0;x<npops-oo;x++) {
					fprintf(file_out,"\t%.3f",jfd[x][zz]);
				}
			}
		}
        fprintf(file_out,"\n");
		if(missratio > 0.) {
			fprintf(file_out,"\n\nNumber of samples for each variant and population (No included variants are missing or polymorphic in the outgroup):");
			fprintf(file_out,"\n");
			for(x=0;x<npops-oo;x++) fprintf(file_out,"\tpop[%d]",x);
			for(zz=0;zz<length_seg;zz++) {
				ss=0; for(x=0;x<npops-oo;x++) {ss += nfd[x][zz];}
				if(ss) {
                    if(formatfile ==3) {
                        if(matrix_pos[zz] < 0) fprintf(file_out,"\nSNP[%ld]",matrix_pos[zz]-(long int)vector_priors[0]+1);
                        else fprintf(file_out,"\nSNP[%ld]",matrix_pos[zz]+(long int)vector_priors[0]-1);
                    }
                    else fprintf(file_out,"\nSNP[%ld]",matrix_pos[zz]);
					for(x=0;x<npops-oo;x++) {
						fprintf(file_out,"\t%d",nfd[x][zz]);
					}
				}
			}
		}
        fprintf(file_out,"\n");
		if(include_unknown == 0) {
			if(ploidy[0] == '1') {
				fprintf(file_out,"\n\nFrequency of variants for each line and population:");
				if(outgroup_presence==0) {
					for(x=0;x<npops-oo;x++) {
						if(nsam[x] > 1) {
							fprintf(file_out,"\nPop[%d]\n\t",x);
							for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) fprintf(file_out,"line[%d]\t",ss-initsq1[x]);
							fprintf(file_out,"\n");
							for(z1=1;z1<=(long int)floor(nsam[x]/2);z1++) {
								fprintf(file_out,"freq[%d]:\t",z1);
								for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
									fprintf(file_out,"%ld\t",(long int)statistics[0].linefreq[ss][z1]);
								}
								fprintf(file_out,"\n");
							}
						}
					}
				}
				else {
					for(x=0;x<npops-oo;x++) {
						if(nsam[x] > 1) {
							fprintf(file_out,"\nPop[%d]\n\t",x);
							for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) fprintf(file_out,"line[%d]\t",ss-initsq1[x]);
							fprintf(file_out,"\n");
							for(z1=1;z1<nsam[x];z1++) {
								fprintf(file_out,"freq[%d]:\t",z1);
								for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
									fprintf(file_out,"%ld\t",(long int)statistics[0].linefreq[ss][z1]);
								}
								fprintf(file_out,"\n");
							}
						}
					}
				}
			}
		}
		
		if(output == 10) {
			if(npops-!outgroup_presence==1) oo = 0;
			else oo = 1;
			/*header*/
			fprintf(file_out,"\n#dadi format for joint frequency spectrum.\n");
			/*fprintf(file_out,"#Note: All alleles are defined arbitrarily to A and T.\n");*/
			fprintf(file_out,"Ref_int\tRef_out\tAllele1\t");
			for(x=0;x<npops-oo;x++) fprintf(file_out,"Pop_%03d\t",x);
			fprintf(file_out,"Allele2\t");
			for(x=0;x<npops-oo;x++) fprintf(file_out,"Pop_%03d\t",x);
			fprintf(file_out,"Position");
			fprintf(file_out, "\n");
			/*table*/
			for(zz=0;zz<length_seg;zz++) {
				ss=0; sf = 0.0; for(x=0;x<npops-oo;x++) {ss += nfd[x][zz]; sf += jfd[x][zz];}
				if(ss && sf > 0.0 && sf < 1.0) {
					/*Define nucleotide*/
                    for(y=0;y<sumnsam;y++) {
                        if(matrix_pol[zz*sumnsam+y] == '0') {
                            if(matrix_pol_tcga[zz*sumnsam+y] == '1')
                                nt1[0] = 'T';
                            if(matrix_pol_tcga[zz*sumnsam+y] == '2')
                                nt1[0] = 'C';
                            if(matrix_pol_tcga[zz*sumnsam+y] == '3')
                                nt1[0] = 'G';
                            if(matrix_pol_tcga[zz*sumnsam+y] == '4')
                                nt1[0] = 'A';
                            break;
                        }
                    }
                    for(y=0;y<sumnsam;y++) {
                        if(matrix_pol[zz*sumnsam+y] == '1') {
                            if(matrix_pol_tcga[zz*sumnsam+y] == '1')
                                nt2[0] = 'T';
                            if(matrix_pol_tcga[zz*sumnsam+y] == '2')
                                nt2[0] = 'C';
                            if(matrix_pol_tcga[zz*sumnsam+y] == '3')
                                nt2[0] = 'G';
                            if(matrix_pol_tcga[zz*sumnsam+y] == '4')
                                nt2[0] = 'A';
                            break;
                        }
                    }
					fprintf(file_out, "-%c-\t",nt1[0]);
					if(outgroup_presence) fprintf(file_out, "-%c-\t",nt2[0]);
					else fprintf(file_out, "---\t");
					/*Allele1*/
					fprintf(file_out, "%c\t",nt1[0]);
					for(x=0;x<npops-oo;x++) fprintf(file_out,"%.0f\t",(jfd[x][zz])*nfd[x][zz]);
					/*Allele2*/
					fprintf(file_out, "%c\t",nt2[0]);
					for(x=0;x<npops-oo;x++) fprintf(file_out,"%.0f\t",(1.0-jfd[x][zz])*nfd[x][zz]);
					/*position*/
					if(formatfile==3) fprintf(file_out,"%ld",labs(matrix_pos[zz])+(long int)vector_priors[0]-1);
                    else fprintf(file_out,"%ld",labs(matrix_pos[zz]));
					fprintf(file_out,"\n");
				}
			}
			/*fprintf(file_out,"\n");*/

			fprintf(file_out,"\n\nAll pairwise comparisons (mismatch distribution) per population:");
			for(x=0;x<npops-oo;x++) {
				fprintf(file_out,"\npop[%d:%ld combinations]:",x,(long int)(nsam[x])*(nsam[x]-1)/2);
				for(zz=0;zz<nsam[x]*(nsam[x]-1)/2;zz++) {
					if(ploidy[0] == '1') {fprintf(file_out,"\t%.2f",statistics[0].mdw[x][zz]);}
					if(ploidy[0] == '2') {fprintf(file_out,"\tNA");}					
				}
			}
			if(ploidy[0] == '1') {
				fprintf(file_out,"\n\nFrequency of SNPs of each line :");
				if(outgroup_presence==0) {
					for(x=0;x<npops-oo;x++) {
						fprintf(file_out,"\n\npop[%d]",x);
						for(ss=1;ss<=(int)floor(nsam[x]/2);ss++) fprintf(file_out,"\tfreq[%d]",ss);
						for(z=initsq1[x];z<initsq1[x]+nsam[x];z++) {
							fprintf(file_out,"\nline[%d]:",z);
							for(ss=1;ss<=(long int)floor(nsam[x]/2);ss++) {
								fprintf(file_out,"\t%.3f",statistics[0].linefreq[z][ss]);
							}
						}
					}
				}
				else {
					for(x=0;x<npops-oo;x++) {
						fprintf(file_out,"\n\npop[%d]",x);
						for(ss=1;ss<nsam[x];ss++) fprintf(file_out,"\tfreq[%d]",ss);
						for(z=initsq1[x];z<initsq1[x]+nsam[x];z++) {
							fprintf(file_out,"\nline[%d]:",z);
							for(ss=1;ss<nsam[x];ss++) {
								fprintf(file_out,"\t%.3f",statistics[0].linefreq[z][ss]);
							}
						}
					}
				}
				/*
				fprintf(file_out,"\n\nCovariance matrix of SNPs frequency per line based on SNM:");
				if(outgroup_presence==0) {
					for(x=0;x<npops-oo;x++) {
						fprintf(file_out,"\n\npop[%d]:",x);
						for(ss=1;ss<=(long int)floor(nsam[x]/2);ss++) fprintf(file_out,"\tfreq[%d]",ss);
						for(z1=1;z1<=(long int)floor(nsam[x]/2);z1++) {
							fprintf(file_out,"\nfreq[%d]",z1);
							for(z2=1;z2<=(long int)floor(nsam[x]/2);z2++) {
								cv = 0.;
								for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
									mean_z1 = mean_z2 = statistics[0].thetaS[x]/nsam[x];
									if(!((float)z1 == floor(nsam[x]/2) && (float)z1 == (float)nsam[x]/2.)) 
										mean_z1 *= 2.;
									if(!((float)z2 == floor(nsam[x]/2) && (float)z2 == (float)nsam[x]/2.)) 
										mean_z2 *= 2.;
									
									cv += (statistics[0].linefreq[ss][z1] - mean_z1) *
										  (statistics[0].linefreq[ss][z2] - mean_z2);
								}
								fprintf(file_out,"\t%.3f",cv/nsam[x]);
							}
						}
					}
				}
				else {
					for(x=0;x<npops-oo;x++) {
						fprintf(file_out,"\n\npop[%d]:",x);
						for(ss=1;ss<nsam[x];ss++) fprintf(file_out,"\tfreq[%d]",ss);
						for(z1=1;z1<nsam[x];z1++) {
							fprintf(file_out,"\nfreq[%d]",z1);
							for(z2=1;z2<nsam[x];z2++) {
								cv = 0.;
								for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
									cv += (statistics[0].linefreq[ss][z1] - (statistics[0].thetaS[x]/nsam[x])) *
										  (statistics[0].linefreq[ss][z2] - (statistics[0].thetaS[x]/nsam[x]));
								}
								fprintf(file_out,"\t%.3f",cv/nsam[x]);
							}
						}
					}
				}
				fprintf(file_out,"\n\nCovariance matrix of SNPs frequency per line based on empirical data:");
				if(outgroup_presence==0) {
					for(x=0;x<npops-oo;x++) {
						fprintf(file_out,"\n\npop[%d]:",x);
						for(ss=1;ss<=(long int)floor(nsam[x]/2);ss++) fprintf(file_out,"\tfreq[%d]",ss);
						for(z1=1;z1<=(long int)floor(nsam[x]/2);z1++) {
							fprintf(file_out,"\nfreq[%d]",z1);
							for(z2=1;z2<=(long int)floor(nsam[x]/2);z2++) {
								cv = 0.;
								for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
									mean_z1 = statistics[0].freq[x][z1] * z1/nsam[x];
									if(!((float)z1 == floor(nsam[x]/2) && (float)z1 == (float)nsam[x]/2.)) 
										mean_z1 += statistics[0].freq[x][z1] * z1/nsam[x];
									mean_z2 = statistics[0].freq[x][z2] * z2/nsam[x];
									if(!((float)z2 == floor(nsam[x]/2) && (float)z2 == (float)nsam[x]/2.)) 
										mean_z2 += statistics[0].freq[x][z2] * z2/nsam[x];
									
									cv += (statistics[0].linefreq[ss][z1] - mean_z1) *
										  (statistics[0].linefreq[ss][z2] - mean_z2);
								}
								fprintf(file_out,"\t%.3f",cv/nsam[x]);
							}
						}
					}
				}
				else {
					for(x=0;x<npops-oo;x++) {
						fprintf(file_out,"\n\npop[%d]:",x);
						for(ss=1;ss<nsam[x];ss++) fprintf(file_out,"\tfreq[%d]",ss);
						for(z1=1;z1<nsam[x];z1++) {
							fprintf(file_out,"\nfreq[%d]",z1);
							for(z2=1;z2<nsam[x];z2++) {
								cv = 0.;
								for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
									cv += (statistics[0].linefreq[ss][z1] - (statistics[0].freq[x][z1] * z1/nsam[x])) *
										  (statistics[0].linefreq[ss][z2] - (statistics[0].freq[x][z2] * z2/nsam[x]));
								}
								fprintf(file_out,"\t%.3f",cv/nsam[x]);
							}
						}
					}
				}
				*/
			}
		}
		fprintf(file_out,"\n");
		/*
		if (file_output[0]) 
			fprintf(file_out,"\nResults sent to file %s\nProgram finished succesfully.\n",file_output);
		*/fflush(file_out);
	}
	else { /*single stdout line*/		
		if(output == 6) {
			/*Matrix SNPs*/
			fprintf(file_out,"#Run: %ld",niterdata);
            fprintf(file_out,"\nscaffold_name: %s",chr_name);
			fprintf(file_out,"\n#SNP chr");
			for(x=0;x<npops-!outgroup_presence;x++) { for(y=0;y<nsam[x];y++) { fprintf(file_out," ind%03d%03d",x,y);}}
			fprintf(file_out,"\n");
            if(force_outgroup==1 || outgroup_presence == 0) oo = 1; else oo = 0;
			/*table*/
			for(zz=0;zz<length_seg;zz++) {
				/* ancestral[0] = 0; */
				if(matrix_pos[zz]<=0) fabsmatsize = -matrix_pos[zz];
				else fabsmatsize = matrix_pos[zz];
				fprintf(file_out,"%-9ld %-9c",fabsmatsize,'0');
				
				/*calculate the ancestral sequence if outgroup*/
				/*
                if(outgroup_presence+force_outgroup==1) {
					freqo[0]=freqo[1]=freqo[2]=freqo[3]=0;
					for(y=sumnsam-nsam[npops-1];y<sumnsam;y++) {
                        if(matrix_pol[zz*sumnsam+y] == '0') {
                            freqo[1] += 1;freqo[0] += 1;
                        }
						if(matrix_pol[zz*sumnsam+y] == '1') {
                            freqo[2] += 1;freqo[0] += 1;
                        }
						if(matrix_pol[zz*sumnsam+y] == '-') {
                            freqo[3] += 1;}
					}
					if(freqo[0]) {
						if(freqo[1] != freqo[0] && freqo[1] != 0) {
							ancestral[0] = (char)0;
						}
						else {
							if(freqo[1] == freqo[0]) {
								ancestral[0] = '0';
							}
							if(freqo[2] == freqo[0]) {
								ancestral[0] = '1';
							}	
						}
					}
					else ancestral[0] = (char)0;
				}
                */
				/*print genotype, including the outgroup if exist*/
				for(y=0;y<sumnsam/*-oo*/;y++) {
                    if(matrix_pol[zz*sumnsam+y] == '0') {
                        if(matrix_pol_tcga[zz*sumnsam+y] == '1')
                            nt1[0] = 'T';
                        if(matrix_pol_tcga[zz*sumnsam+y] == '2')
                            nt1[0] = 'C';
                        if(matrix_pol_tcga[zz*sumnsam+y] == '3')
                            nt1[0] = 'G';
                        if(matrix_pol_tcga[zz*sumnsam+y] == '4')
                            nt1[0] = 'A';
                    }
                    if(matrix_pol[zz*sumnsam+y] == '1') {
                        if(matrix_pol_tcga[zz*sumnsam+y] == '1')
                            nt2[0] = 'T';
                        if(matrix_pol_tcga[zz*sumnsam+y] == '2')
                            nt2[0] = 'C';
                        if(matrix_pol_tcga[zz*sumnsam+y] == '3')
                            nt2[0] = 'G';
                        if(matrix_pol_tcga[zz*sumnsam+y] == '4')
                            nt2[0] = 'A';
                    }
					if(matrix_pol[zz*sumnsam+y] == '0') fprintf(file_out," %c",nt1[0]);
                    if(matrix_pol[zz*sumnsam+y] == '1') fprintf(file_out," %c",nt2[0]);
                    if(matrix_pol[zz*sumnsam+y] == '-') fprintf(file_out," N");
                    /*
                    if(ancestral[0] == matrix_pol[zz*sumnsam+y]) fprintf(file_out," 0");
					else {
						if(matrix_pol[zz*sumnsam+y] == '-') fprintf(file_out," 9");
						else {
							if(ancestral[0] == 0) fprintf(file_out," %c",matrix_pol[zz*sumnsam+y]);
							else {
								if(ancestral[0] != matrix_pol[zz*sumnsam+y]) fprintf(file_out," 1");
							}
						}
					}
                    */
				}
				fprintf(file_out,"\n");
			}
			/**/
		}
		else {
			if(output == 3) {
				if(npops-!outgroup_presence==1) oo = 0;
				else oo = 1;
				/*header*/
				fprintf(file_out,"#mstatspop output for dadi format of joint frequency spectrum.\n");
				fprintf(file_out,"#Input file: %s\n",file_input);
                fprintf(file_out,"#scaffold_name: %s\n",chr_name);
				fprintf(file_out,"#Note: All alleles are defined arbitrarily to A and T.\n");
				fprintf(file_out,"Ref_int\tRef_out\tAllele1\t");
				for(x=0;x<npops-oo;x++) fprintf(file_out,"Pop_%03d\t",x);
				fprintf(file_out,"Allele2\t");
				for(x=0;x<npops-oo;x++) fprintf(file_out,"Pop_%03d\t",x);
				fprintf(file_out,"Position\t");
				if(*file_input) fprintf(file_out,"Location");
				/*table*/
				for(zz=0;zz<length_seg;zz++) {
					ss=0; sf = 0.0; for(x=0;x<npops-oo;x++) {ss += nfd[x][zz]; sf += jfd[x][zz];}
					if(ss && sf > 0.0 && sf < 1.0*(npops-oo)) {
						fprintf(file_out,"\n");
						/*Define nucleotide*/
                        for(y=0;y<sumnsam;y++) {
                            if(matrix_pol[zz*sumnsam+y] == '0') {
                                if(matrix_pol_tcga[zz*sumnsam+y] == '1')
                                    nt1[0] = 'T';
                                if(matrix_pol_tcga[zz*sumnsam+y] == '2')
                                    nt1[0] = 'C';
                                if(matrix_pol_tcga[zz*sumnsam+y] == '3')
                                    nt1[0] = 'G';
                                if(matrix_pol_tcga[zz*sumnsam+y] == '4')
                                    nt1[0] = 'A';
                                break;
                            }
                        }
                        for(y=0;y<sumnsam;y++) {
                            if(matrix_pol[zz*sumnsam+y] == '1') {
                                if(matrix_pol_tcga[zz*sumnsam+y] == '1')
                                    nt2[0] = 'T';
                                if(matrix_pol_tcga[zz*sumnsam+y] == '2')
                                    nt2[0] = 'C';
                                if(matrix_pol_tcga[zz*sumnsam+y] == '3')
                                    nt2[0] = 'G';
                                if(matrix_pol_tcga[zz*sumnsam+y] == '4')
                                    nt2[0] = 'A';
                                break;
                            }
                        }
                        fprintf(file_out, "-%c-\t",nt1[0]);
                        if(outgroup_presence) fprintf(file_out, "-%c-\t",nt2[0]);
                        else fprintf(file_out, "---\t");
                        /*Allele1*/
                        fprintf(file_out, "%c\t",nt1[0]);
                        for(x=0;x<npops-oo;x++) fprintf(file_out,"%.0f\t",(jfd[x][zz])*nfd[x][zz]);
                        /*Allele2*/
                        fprintf(file_out, "%c\t",nt2[0]);
                        for(x=0;x<npops-oo;x++) fprintf(file_out,"%.0f\t",(1.0-jfd[x][zz])*nfd[x][zz]);
                        /*position*/
                        fprintf(file_out,"%ld",labs(matrix_pos[zz]));
						if(*file_input) fprintf(file_out,"%s",file_input);
					}
				}
				fprintf(file_out,"\n");
            }
			else {
				if(file_out == 0) file_out = stdout;
				if(*file_input)
					fprintf(file_out,"infile:\t%s\t",file_input);
				else
					fprintf(file_out,"infile:\tNA\t");
				if(gfffiles == 1) {
					fprintf(file_out,"GFFfile:\t%s\t",file_GFF);
					fprintf(file_out,"subset_pos:\t%s\t",subset_positions);
					if(strcmp(subset_positions,"synonymous")==0 || strcmp(subset_positions,"nonsynonymous")==0 || strcmp(subset_positions,"silent")==0)
						fprintf(file_out,"genetic_code:\t%s\t",code_name);
				}
				/**/
				if(npriors) {
                    if(formatfile != 3) {
                        for(x=0;x<npriors;x++) {
                            fprintf(file_out,"prior%02d:\t%f\t",x,vector_priors[x]);
                        }
                    }
                    if(formatfile == 3) {
                        if(npriors == 2) {
                            fprintf(file_out,"scaffold_name:\t%s\t",chr_name);
                            fprintf(file_out,"start_window:\t%.0f\t",vector_priors[0]);
                            fprintf(file_out,"end_window:\t%.0f\t",vector_priors[1]);
                        }
                        else {
                            for(x=0;x<npriors;x++) {
                                fprintf(file_out,"prior%02d:\t%f\t",x,vector_priors[x]);
                            }
                        }
                    }
                }
				/**/
				if(H1frq) {
					fprintf(file_out,"H1file:\t%s\t",file_H1f);
					if(H0frq) {
						fprintf(file_out,"H0file:\t%sºt",file_H0f);
					}
				}
				fprintf(file_out,"missing:\t%d\t",include_unknown);
				fprintf(file_out,"iteration:\t%ld\t",niterdata);
				fprintf(file_out,"npermutations:\t%ld\t",niter);
				if(formatfile==0 || formatfile==3) fprintf(file_out,"seed:\t%ld\t",nseed);
				else fprintf(file_out,"seed:\tna\t");
				if(formatfile==0 || formatfile==3) fprintf(file_out,"Length:\t%.2f\t",length_al);
				else fprintf(file_out,"Length:\t%.2f\t",length_al);
				fprintf(file_out,"Lengtht:\t%ld\t",length_al_real);
				fprintf(file_out,"mh:\t%ld\t",statistics[0].nmhits);
				if(svratio > -10000) fprintf(file_out,"Ratio_S/V:\t%.3f\t",svratio);
				else fprintf(file_out,"Ratio_S/V:\tNA\t");
				if(include_unknown == 1)fprintf(file_out,"Ratio_Missing:\t%f\t",missratio);
                else fprintf(file_out,"Ratio_Missing:\tNA\t");
				fprintf(file_out,"Variants:\t%ld\t",length_seg);
				
				fprintf(file_out,"npops:\t%d\t",npops-!outgroup_presence);
				for(x=0;x<npops-!outgroup_presence;x++) {
					fprintf(file_out,"nsam[%d]:\t%d\t",x,nsam[x]);
				}
				for(x=0;x<npops-1;x++) {
                    /*if(outgroup_presence==1) {fprintf(file_out,"Eff_length1_pop[%d]:\t%.2f\tEff_length2_pop[%d]:\t%.2f\tEff_length3_pop[%d]:\t%.2f\tEff_length1_pop_outg[%d]:\t%.2f\tEff_length2_pop_outg[%d]:\t%.2f\tEff_length3_pop_outg[%d]:\t%.2f\t",x,(double)nsites1_pop[x],x,(double)nsites2_pop[x],x,(double)nsites3_pop[x],x,(double)nsites1_pop_outg[x],x,(double)nsites2_pop_outg[x],x,(double)nsites3_pop_outg[x]);}*/
                    if(outgroup_presence==1) {fprintf(file_out,"Eff_length1_pop_outg[%d]:\t%.2f\tEff_length2_pop_outg[%d]:\t%.2f\tEff_length3_pop_outg[%d]:\t%.2f\t",x,(double)nsites1_pop_outg[x],x,(double)nsites2_pop_outg[x],x,(double)nsites3_pop_outg[x]);}
					if(outgroup_presence==0) {fprintf(file_out,"Eff_length1_pop[%d]:\t%.2f\tEff_length2_pop[%d]:\t%.2f\tEff_length3_pop[%d]:\t%.2f\t",x,(double)nsites1_pop[x],x,(double)nsites2_pop[x],x,(double)nsites3_pop[x]);}
				}
				oo = 1;
				
				if(output == 11) {
					if(ploidy[0] == '1') {
						if(outgroup_presence==0) {
							for(x=0;x<npops-oo;x++) {
								fprintf(file_out,"Cov_Star_pop[%d]:\t",x);
								for(z1=1;z1<=(long int)floor(nsam[x]/2);z1++) {
									for(z2=1;z2<=(long int)floor(nsam[x]/2);z2++) {
										cv = 0.;
										for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
											mean_z1 = mean_z2 = 0.;
											if(z1 == 1) mean_z1 = statistics[0].thetaT[x]/2.;
											if(z2 == 1) mean_z2 = statistics[0].thetaT[x]/2.;
											
											cv += (statistics[0].linefreq[ss][z1] - mean_z1) *
											(statistics[0].linefreq[ss][z2] - mean_z2);
										}
										cv = cv/(double)nsam[x];
		#if DENCOV_CORRECTION == 1
										if(statistics[0].freq[x][z2] * statistics[0].freq[x][z1] > 0.) {
											if(z1==z2) delta =1.; 
											else delta=0.;
											cv = cv - delta * (double)z1/(double)nsam[x] * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1];
											if((der = (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta))) != 0.)
												cv = cv/der;
										}
										else cv = 0.;
		#endif
										fprintf(file_out,"%.6f\t",cv);
									}
								}
							}
						}
						else {
							for(x=0;x<npops-oo;x++) {
								fprintf(file_out,"Cov_Star_pop[%d]:\t",x);
								for(z1=1;z1<nsam[x];z1++) {
									for(z2=1;z2<nsam[x];z2++) {
										cv = 0.;
										for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
											mean_z1 = mean_z2 = 0.;
											if(z1 == 1) mean_z1 = statistics[0].thetaT[x]/2.;
											if(z2 == 1) mean_z2 = statistics[0].thetaT[x]/2.;
											cv += (statistics[0].linefreq[ss][z1] - mean_z1) *
											(statistics[0].linefreq[ss][z2] - mean_z2);
										}
										cv = cv/(double)nsam[x];
		#if DENCOV_CORRECTION == 1
										if(statistics[0].freq[x][z2] * statistics[0].freq[x][z1] > 0.) {
											if((der = (statistics[0].freq[x][z1]*statistics[0].freq[x][z2]))!=0.)
												cv = cv/der;
										}
										else cv = 0.;
		#endif
										fprintf(file_out,"%.6f\t",cv);
									}
								}
							}
						}
						fprintf(file_out,"\n");
					}
				}
				if(output == 9) {
					if(ploidy[0] == '1') {
						if(outgroup_presence==0) {
							for(x=0;x<npops-oo;x++) {
								fprintf(file_out,"linefreq_pop[%d]:\t",x);
								for(z1=1;z1<=(long int)floor(nsam[x]/2);z1++) {
									for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
										fprintf(file_out,"%ld\t",(long int)statistics[0].linefreq[ss][z1]);
									}
								}
							}
						}
						else {
							for(x=0;x<npops-oo;x++) {
								fprintf(file_out,"linefreq_pop[%d]:\t",x);
								for(z1=1;z1<nsam[x];z1++) {
									for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
										fprintf(file_out,"%ld\t",(long int)statistics[0].linefreq[ss][z1]);
									}
								}
							}
						}
						fprintf(file_out,"\n");
					}
				}
				if(output == 8) {
					if(ploidy[0] == '1' && include_unknown == 0) {
						for(y=0;y<statistics[0].nh-!outgroup_presence;y++) {
							for(x=0;x<npops-1;x++) {
								fprintf(file_out,"frHap[%d,hap%02d]:\t%ld\t",x,y,statistics[0].freqh[x][y]);
							}
						}
						for(y=statistics[0].nh-!outgroup_presence;y<sumnsam;y++) {
							for(x=0;x<npops-1;x++) {
								fprintf(file_out,"frHap[%d,hap%02d]:\t0\t",x,y);
							}
						}
					}
					fprintf(file_out,"\n");
				}
				if(output == 70) {
					if(ploidy[0] == '1') {/*haploid*/
						if(outgroup_presence==0) {
							for(x=0;x<npops-oo;x++) {
								#if MATRIXCOV == 1
								fprintf(file_out,"Cov_EMP_pop[%d]:\t",x);
								for(z1=1;z1<=(long int)floor(nsam[x]/2);z1++) {
									for(z2=1;z2<=(long int)floor(nsam[x]/2);z2++) {
										cv = 0.;
										for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
											mean_z1 = statistics[0].freq[x][z1] * (float)z1/(float)nsam[x];
											if(!((float)z1 == floor(nsam[x]/2) && (float)z1 == (float)nsam[x]/2.)) 
												mean_z1 += statistics[0].freq[x][z1] * (float)z1/(float)nsam[x];
											mean_z2 = statistics[0].freq[x][z2] * (float)z2/(float)nsam[x];
											if(!((float)z2 == floor(nsam[x]/2) && (float)z2 == (float)nsam[x]/2.)) 
												mean_z2 += statistics[0].freq[x][z2] * (float)z2/(float)nsam[x];
											
											cv += (statistics[0].linefreq[ss][z1] - mean_z1) *
												  (statistics[0].linefreq[ss][z2] - mean_z2);
										}
										cv = cv/(double)nsam[x];
										#if DENCOV_CORRECTION == 1
										if(statistics[0].freq[x][z2] * statistics[0].freq[x][z1] > 0.) {
											if(z1==z2) delta =1.; 
											else delta=0.;
											cv = cv - delta * (double)z1/(double)nsam[x] * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1];
											if((der = (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta))) != 0.)
											cv = cv/der;
										}
										else cv = 0.;
										#endif
										fprintf(file_out,"%.6f\t",cv);
									}
								}
								#endif
								#if LUCA_CRs == 1
								/*LUCA: here I calculate the statistics for folded tests with different weights 1/i^2,1/i,1,i*/
								luca_cvo2=0.; luca_cvo1=0.; luca_cv0=0.; luca_cv1=0.;
								luca_dero2=0.; luca_dero1=0.; luca_der0=0.; luca_der1=0.;
								luca_cvo2d=0.; luca_cvo1d=0.; luca_cv0d=0.; luca_cv1d=0.;
								luca_dero2d=0.; luca_dero1d=0.; luca_der0d=0.; luca_der1d=0.;
								for(z1=1;z1<=(long int)floor(nsam[x]/2);z1++) {
									for(z2=1;z2<=(long int)floor(nsam[x]/2);z2++) {
										for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
											mean_z1 = statistics[0].freq[x][z1] * (float)z1/(float)nsam[x];
											if(!((float)z1 == floor(nsam[x]/2) && (float)z1 == (float)nsam[x]/2.)) 
												mean_z1 += statistics[0].freq[x][z1] * (float)z1/(float)nsam[x];
											mean_z2 = statistics[0].freq[x][z2] * (float)z2/(float)nsam[x];
											if(!((float)z2 == floor(nsam[x]/2) && (float)z2 == (float)nsam[x]/2.)) 
												mean_z2 += statistics[0].freq[x][z2] * (float)z2/(float)nsam[x];
											if(z1==z2) delta =1.; 
											else delta=0.;
											
											luca_cvo2 += 1/(double)(z1*z1*z2*z2) * ((statistics[0].linefreq[ss][z1] - mean_z1) *
												  (statistics[0].linefreq[ss][z2] - mean_z2)
												  - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
											luca_cvo1 += 1/(double)(z1*z2) * ((statistics[0].linefreq[ss][z1] - mean_z1) *
												  (statistics[0].linefreq[ss][z2] - mean_z2)
												  - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
											luca_cv0 += ((statistics[0].linefreq[ss][z1] - mean_z1) *
												  (statistics[0].linefreq[ss][z2] - mean_z2)
												  - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
											luca_cv1 += (double)(z1*z2) * ((statistics[0].linefreq[ss][z1] - mean_z1) *
												  (statistics[0].linefreq[ss][z2] - mean_z2)
												  - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
											luca_dero2 += 1/(double)(z1*z1*z2*z2) * (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
											luca_dero1 += 1/(double)(z1*z2) * (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
											luca_der0 += (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
											luca_der1 += (double)(z1*z2) * (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
											
											if(z1==z2) {
												luca_cvo2d += 1/(double)(z1*z1*z2*z2) * ((statistics[0].linefreq[ss][z1] - mean_z1) *
																						(statistics[0].linefreq[ss][z2] - mean_z2)
																						- delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
												luca_cvo1d += 1/(double)(z1*z2) * ((statistics[0].linefreq[ss][z1] - mean_z1) *
																				  (statistics[0].linefreq[ss][z2] - mean_z2)
																				  - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
												luca_cv0d += ((statistics[0].linefreq[ss][z1] - mean_z1) *
															 (statistics[0].linefreq[ss][z2] - mean_z2)
															 - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
												luca_cv1d += (double)(z1*z2) * ((statistics[0].linefreq[ss][z1] - mean_z1) *
																			   (statistics[0].linefreq[ss][z2] - mean_z2)
																			   - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
												luca_dero2d += 1/(double)(z1*z1*z2*z2) * (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
												luca_dero1d += 1/(double)(z1*z2) * (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
												luca_der0d += (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
												luca_der1d += (double)(z1*z2) * (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
											}
										}
									}
								}
								
								luca_cvo2 = luca_cvo2/luca_dero2/(double)nsam[x];
								luca_cvo1 = luca_cvo1/luca_dero1/(double)nsam[x];
								luca_cv0 = luca_cv0/luca_der0/(double)nsam[x];
								luca_cv1 = luca_cv1/luca_der1/(double)nsam[x];
								
								luca_cvo2d = luca_cvo2d/luca_dero2d/(double)nsam[x];
								luca_cvo1d = luca_cvo1d/luca_dero1d/(double)nsam[x];
								luca_cv0d = luca_cv0d/luca_der0d/(double)nsam[x];
								luca_cv1d = luca_cv1d/luca_der1d/(double)nsam[x];

								fprintf(file_out,"CVO2[%d]:\t%.6f\t",x,luca_cvo2);
								fprintf(file_out,"CVO1[%d]:\t%.6f\t",x,luca_cvo1);
								fprintf(file_out,"CV0[%d]:\t%.6f\t",x,luca_cv0);
								fprintf(file_out,"CV1[%d]:\t%.6f\t",x,luca_cv1);
								fprintf(file_out,"CVO2d[%d]:\t%.6f\t",x,luca_cvo2d);
								fprintf(file_out,"CVO1d[%d]:\t%.6f\t",x,luca_cvo1d);
								fprintf(file_out,"CV0d[%d]:\t%.6f\t",x,luca_cv0d);
								fprintf(file_out,"CV1d[%d]:\t%.6f\t",x,luca_cv1d);
								fprintf(file_out,"R2[%d]:\t%.6f\t",x,statistics[0].R2[x]);
								/*end LUCA*/
								#endif
							}
						}
						else {
							for(x=0;x<npops-oo;x++) {
								#if MATRIXCOV == 1
								fprintf(file_out,"Cov_EMP_pop[%d]:\t",x);
								for(z1=1;z1<nsam[x];z1++) {
									for(z2=1;z2<nsam[x];z2++) {
										cv = 0.;
										for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
											cv += (statistics[0].linefreq[ss][z1] - (statistics[0].freq[x][z1] * (float)z1/(float)nsam[x])) *
												  (statistics[0].linefreq[ss][z2] - (statistics[0].freq[x][z2] * (float)z2/(float)nsam[x]));
										}
										cv = cv/(double)nsam[x];
										#if DENCOV_CORRECTION == 1
										if(statistics[0].freq[x][z2] * statistics[0].freq[x][z1] > 0.) {
											if((der = (statistics[0].freq[x][z1]*statistics[0].freq[x][z2]))!=0.)
												cv = cv/der;
										}
										else cv = 0.;
										#endif
										fprintf(file_out,"%.6f\t",cv);
									}
								}
								#endif
								#if LUCA_CRs == 1
								/*LUCA: here I calculate the statistics for unfolded tests with different weights 1/i^2,1/i,1,i*/
								luca_cvo2=0.; luca_cvo1=0.; luca_cv0=0.; luca_cv1=0.;
								luca_cro2=0.; luca_cro1=0.; luca_cr0=0.; luca_cr1=0.;
								luca_cso2=0.; luca_cso1=0.; luca_cs0=0.; luca_cs1=0.;
								luca_dero2=0.; luca_dero1=0.; luca_der0=0.; luca_der1=0.;
								luca_cvo2d=0.; luca_cvo1d=0.; luca_cv0d=0.; luca_cv1d=0.;
								luca_cro2d=0.; luca_cro1d=0.; luca_cr0d=0.; luca_cr1d=0.;
								luca_cso2d=0.; luca_cso1d=0.; luca_cs0d=0.; luca_cs1d=0.;
								luca_dero2d=0.; luca_dero1d=0.; luca_der0d=0.; luca_der1d=0.;
								for(z1=1;z1<nsam[x];z1++) {
									for(z2=1;z2<=nsam[x]-z1;z2++) {
										for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
											mean_z1 = (statistics[0].freq[x][z1] * (float)z1/(float)nsam[x]);
											mean_z2 = (statistics[0].freq[x][z2] * (float)z2/(float)nsam[x]);
											luca_mean_z1 = (statistics[0].freq[x][nsam[x]-z1] * (float)(nsam[x]-z1)/(float)nsam[x]);
											luca_mean_z2 = (statistics[0].freq[x][nsam[x]-z2] * (float)(nsam[x]-z2)/(float)nsam[x]);
											if(z1==z2) delta =1.; 
											else delta=0.;
											luca_cvo2 += 1/(double)(z1*z1*z2*z2) * ((statistics[0].linefreq[ss][z1] - mean_z1) *
												  (statistics[0].linefreq[ss][z2] - mean_z2)
												  - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
											luca_cvo1 += 1/(double)(z1*z2) * ((statistics[0].linefreq[ss][z1] - mean_z1) *
												  (statistics[0].linefreq[ss][z2] - mean_z2)
												  - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
											luca_cv0 += ((statistics[0].linefreq[ss][z1] - mean_z1) *
												  (statistics[0].linefreq[ss][z2] - mean_z2)
												  - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
											luca_cv1 += (double)(z1*z2) * ((statistics[0].linefreq[ss][z1] - mean_z1) *
												  (statistics[0].linefreq[ss][z2] - mean_z2)
												  - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
											
											luca_dero2 += 1/(double)(z1*z1*z2*z2) * (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
											luca_dero1 += 1/(double)(z1*z2) * (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
											luca_der0 += (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
											luca_der1 += (double)(z1*z2) * (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
											
											luca_cso2 += 1/(double)(z1*z1*z2*z2) * ((statistics[0].linefreq[ss][nsam[x]-z1] - luca_mean_z1) *
												  (statistics[0].linefreq[ss][nsam[x]-z2] - luca_mean_z2)
												  - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][nsam[x]-z1]);
											luca_cso1 += 1/(double)(z1*z2) * ((statistics[0].linefreq[ss][nsam[x]-z1] - luca_mean_z1) *
												  (statistics[0].linefreq[ss][nsam[x]-z2] - luca_mean_z2)
												  - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][nsam[x]-z1]);
											luca_cs0 += ((statistics[0].linefreq[ss][nsam[x]-z1] - luca_mean_z1) *
												  (statistics[0].linefreq[ss][nsam[x]-z2] - luca_mean_z2)
												  - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][nsam[x]-z1]);
											luca_cs1 += (double)(z1*z2) * ((statistics[0].linefreq[ss][nsam[x]-z1] - luca_mean_z1) *
												  (statistics[0].linefreq[ss][nsam[x]-z2] - luca_mean_z2)
												  - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][nsam[x]-z1]);
											
											luca_cro2 += ((double)min(z1,z2)/(double)nsam[x]-(double)(z1*z2)/(double)(nsam[x]*nsam[x])) * 
											1/(double)(z1*z1*z2*z2) * (statistics[0].freq[x][nsam[x]-z1]*(statistics[0].freq[x][nsam[x]-z2]-delta));
											luca_cro1 += ((double)min(z1,z2)/(double)nsam[x]-(double)(z1*z2)/(double)(nsam[x]*nsam[x])) * 
											1/(double)(z1*z2) * (statistics[0].freq[x][nsam[x]-z1]*(statistics[0].freq[x][nsam[x]-z2]-delta));
											luca_cr0 += ((double)min(z1,z2)/(double)nsam[x]-(double)(z1*z2)/(double)(nsam[x]*nsam[x])) * 
											(statistics[0].freq[x][nsam[x]-z1]*(statistics[0].freq[x][nsam[x]-z2]-delta));
											luca_cr1 += ((double)min(z1,z2)/(double)nsam[x]-(double)(z1*z2)/(double)(nsam[x]*nsam[x])) * 
											(double)(z1*z2) * (statistics[0].freq[x][nsam[x]-z1]*(statistics[0].freq[x][nsam[x]-z2]-delta));
											
											if(z1==z2) {
												luca_cvo2d += 1/(double)(z1*z1*z2*z2) * ((statistics[0].linefreq[ss][z1] - mean_z1) *
																						(statistics[0].linefreq[ss][z2] - mean_z2)
																						- delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
												luca_cvo1d += 1/(double)(z1*z2) * ((statistics[0].linefreq[ss][z1] - mean_z1) *
																				  (statistics[0].linefreq[ss][z2] - mean_z2)
																				  - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
												luca_cv0d += ((statistics[0].linefreq[ss][z1] - mean_z1) *
															 (statistics[0].linefreq[ss][z2] - mean_z2)
															 - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
												luca_cv1d += (double)(z1*z2) * ((statistics[0].linefreq[ss][z1] - mean_z1) *
																			   (statistics[0].linefreq[ss][z2] - mean_z2)
																			   - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1]);
												
												luca_dero2d += 1/(double)(z1*z1*z2*z2) * (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
												luca_dero1d += 1/(double)(z1*z2) * (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
												luca_der0d += (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
												luca_der1d += (double)(z1*z2) * (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta));
												
												luca_cso2d += 1/(double)(z1*z1*z2*z2) * ((statistics[0].linefreq[ss][nsam[x]-z1] - luca_mean_z1) *
																						(statistics[0].linefreq[ss][nsam[x]-z2] - luca_mean_z2)
																						- delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][nsam[x]-z1]);
												luca_cso1d += 1/(double)(z1*z2) * ((statistics[0].linefreq[ss][nsam[x]-z1] - luca_mean_z1) *
																				  (statistics[0].linefreq[ss][nsam[x]-z2] - luca_mean_z2)
																				  - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][nsam[x]-z1]);
												luca_cs0d += ((statistics[0].linefreq[ss][nsam[x]-z1] - luca_mean_z1) *
															 (statistics[0].linefreq[ss][nsam[x]-z2] - luca_mean_z2)
															 - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][nsam[x]-z1]);
												luca_cs1d += (double)(z1*z2) * ((statistics[0].linefreq[ss][nsam[x]-z1] - luca_mean_z1) *
																			   (statistics[0].linefreq[ss][nsam[x]-z2] - luca_mean_z2)
																			   - delta * (double)z1 * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][nsam[x]-z1]);
												
												luca_cro2d += ((double)min(z1,z2)/(double)nsam[x]-(double)(z1*z2)/(double)(nsam[x]*nsam[x])) * 
												1/(double)(z1*z1*z2*z2) * (statistics[0].freq[x][nsam[x]-z1]*(statistics[0].freq[x][nsam[x]-z2]-delta));
												luca_cro1d += ((double)min(z1,z2)/(double)nsam[x]-(double)(z1*z2)/(double)(nsam[x]*nsam[x])) * 
												1/(double)(z1*z2) * (statistics[0].freq[x][nsam[x]-z1]*(statistics[0].freq[x][nsam[x]-z2]-delta));
												luca_cr0d += ((double)min(z1,z2)/(double)nsam[x]-(double)(z1*z2)/(double)(nsam[x]*nsam[x])) * 
												(statistics[0].freq[x][nsam[x]-z1]*(statistics[0].freq[x][nsam[x]-z2]-delta));
												luca_cr1d += ((double)min(z1,z2)/(double)nsam[x]-(double)(z1*z2)/(double)(nsam[x]*nsam[x])) * 
												(double)(z1*z2) * (statistics[0].freq[x][nsam[x]-z1]*(statistics[0].freq[x][nsam[x]-z2]-delta));
											}
										}
									}
								}
								
								luca_cvo2 = luca_cvo2*luca_cro2/luca_cso2/luca_dero2/(double)nsam[x];
								luca_cvo1 = luca_cvo1*luca_cro1/luca_cso1/luca_dero1/(double)nsam[x];
								luca_cv0 = luca_cv0*luca_cr0/luca_cs0/luca_der0/(double)nsam[x];
								luca_cv1 = luca_cv1*luca_cr1/luca_cs1/luca_der1/(double)nsam[x];
								
								luca_cvo2d = luca_cvo2d*luca_cro2d/luca_cso2d/luca_dero2d/(double)nsam[x];
								luca_cvo1d = luca_cvo1d*luca_cro1d/luca_cso1d/luca_dero1d/(double)nsam[x];
								luca_cv0d = luca_cv0d*luca_cr0d/luca_cs0d/luca_der0d/(double)nsam[x];
								luca_cv1d = luca_cv1d*luca_cr1d/luca_cs1d/luca_der1d/(double)nsam[x];
								
								fprintf(file_out,"CVO2[%d]:\t%.6f\t",x,luca_cvo2);
								fprintf(file_out,"CVO1[%d]:\t%.6f\t",x,luca_cvo1);
								fprintf(file_out,"CV0[%d]:\t%.6f\t",x,luca_cv0);
								fprintf(file_out,"CV1[%d]:\t%.6f\t",x,luca_cv1);
								
								fprintf(file_out,"CVO2d[%d]:\t%.6f\t",x,luca_cvo2d);
								fprintf(file_out,"CVO1d[%d]:\t%.6f\t",x,luca_cvo1d);
								fprintf(file_out,"CV0d[%d]:\t%.6f\t",x,luca_cv0d);
								fprintf(file_out,"CV1d[%d]:\t%.6f\t",x,luca_cv1d);
								fprintf(file_out,"R2[%d]:\t%.6f\t",x,statistics[0].R2[x]);
								/*end LUCA*/
								#endif
							}
						}
						fprintf(file_out,"\n");
					}
				}
				if(output == 60) {
					if(ploidy[0] == '1') {
						if(outgroup_presence==0) {
							for(x=0;x<npops-oo;x++) {
								fprintf(file_out,"Cov_SNM_pop[%d]:\t",x);
								for(z1=1;z1<=(long int)floor(nsam[x]/2);z1++) {
									for(z2=1;z2<=(long int)floor(nsam[x]/2);z2++) {
										cv = 0.;
										for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
											mean_z1 = mean_z2 = statistics[0].thetaS[x]/(float)nsam[x];
											if(!((float)z1 == floor(nsam[x]/2) && (float)z1 == (float)nsam[x]/2.)) 
												mean_z1 *= 2.;
											if(!((float)z2 == floor(nsam[x]/2) && (float)z2 == (float)nsam[x]/2.)) 
												mean_z2 *= 2.;
											
											cv += (statistics[0].linefreq[ss][z1] - mean_z1) *
												  (statistics[0].linefreq[ss][z2] - mean_z2);
										}
										cv = cv/(double)nsam[x];
										#if DENCOV_CORRECTION == 1
										if(statistics[0].freq[x][z2] * statistics[0].freq[x][z1] > 0.) {
											if(z1==z2) delta =1.; 
											else delta=0.;
											cv = cv - delta * (double)z1/(double)nsam[x] * (1.-(double)z1/(double)nsam[x])*statistics[0].freq[x][z1];
											if((der = (statistics[0].freq[x][z1]*(statistics[0].freq[x][z2]-delta))) != 0.)
												cv = cv/der;
										}
										else cv = 0.;
										#endif
										fprintf(file_out,"%.6f\t",cv);
									}
								}
							}
						}
						else {
							for(x=0;x<npops-oo;x++) {
								fprintf(file_out,"Cov_SNM_pop[%d]:\t",x);
								for(z1=1;z1<nsam[x];z1++) {
									for(z2=1;z2<nsam[x];z2++) {
										cv = 0.;
										for(ss=initsq1[x];ss<initsq1[x]+nsam[x];ss++) {
											cv += (statistics[0].linefreq[ss][z1] - (statistics[0].thetaS[x]/(float)nsam[x])) *
												  (statistics[0].linefreq[ss][z2] - (statistics[0].thetaS[x]/(float)nsam[x]));
										}
										cv = cv/(double)nsam[x];
										#if DENCOV_CORRECTION == 1
										if(statistics[0].freq[x][z2] * statistics[0].freq[x][z1] > 0.) {
											if((der = (statistics[0].freq[x][z1]*statistics[0].freq[x][z2]))!=0.)
												cv = cv/der;
										}
										else cv = 0.;
										#endif
										fprintf(file_out,"%.6f\t",cv);
									}
								}
							}
						}
						fprintf(file_out,"\n");
					}
				}
				if(output == 5) {
					if(ploidy[0] == '1') {
						/*fprintf(file_out,"\n\nFrequency of SNPs in a population that contain each line:");*/
						if(outgroup_presence==0) {
							for(x=0;x<npops-oo;x++) {
								for(z=initsq1[x];z<initsq1[x]+nsam[x];z++) {
									for(ss=1;ss<(long int)floor(nsam[x]/2);ss++) {
										fprintf(file_out,"line_freq[%d,%d]:\t%.3f\t",z,ss,statistics[0].linefreq[z][ss]);
									}
								}
							}
						}
						else {
							for(x=0;x<npops-oo;x++) {
								for(z=initsq1[x];z<initsq1[x]+nsam[x];z++) {
									for(ss=1;ss<nsam[x];ss++) {
										fprintf(file_out,"line_freq[%d,%d]:\t%.3f\t",z,ss,statistics[0].linefreq[z][ss]);
									}
								}
							}
						}
					}
					fprintf(file_out,"\n");
				}
				if(output == 4) {
					if(npops-!outgroup_presence==1) oo = 0;
					else oo = 1;
					/*fprintf(file_out,"\nAll pairwise comparisons (mismatch distribution) per population:");*/
					for(x=0;x<npops-oo;x++) {
						fprintf(file_out,"Pairwise_Comparisons_npop[%d:%ld comb]:\t",x,(long int)(nsam[x])*(nsam[x]-1)/2);
						for(zz=0;zz<nsam[x]*(nsam[x]-1)/2;zz++) {
							if(ploidy[0] == '1') {fprintf(file_out,"%.1f\t",statistics[0].mdw[x][zz]);}
							if(ploidy[0] == '2') {fprintf(file_out,"NA\t");}
						}
					}
					fprintf(file_out,"\n");
				}
				if(output == 30) {
					if(npops-!outgroup_presence==1) oo = 0;
					else oo = 1;
					/*fprintf(file_out,"\n");*/
					/*for(x=0;x<npops-oo;x++) fprintf(file_out,"\tFREQpop[%d]",x);*/
					/*if(missratio > 0.) {for(x=0;x<npops-oo;x++) fprintf(file_out,"\tSAMPpop[%d]",x);}*/
					for(zz=0;zz<length_seg;zz++) {
						ss=0; for(x=0;x<npops-oo;x++) {ss += nfd[x][zz];}
						if(ss) {
							fprintf(file_out,"\tSNP[%ld]",matrix_pos[zz]);
							for(x=0;x<npops-oo;x++) fprintf(file_out,"\t%.3f",jfd[x][zz]);
							if(missratio > 0.) {for(x=0;x<npops-oo;x++) fprintf(file_out,"\t%d",nfd[x][zz]);}
						}
					}
					fprintf(file_out,"\n");
				}
				if(output == 2) {
					/*ffprintf(file_out,file_out,"\nFrequency of variants for each population:");*/
					if(outgroup_presence == 1) {		
						for(x=0;x<npops-1;x++) {
							fprintf(file_out,"\tnsam[%d]:\t%d\tS[%d]:\t%ld",x,nsam[x],x,(long int)statistics[0].S[x]);
							for(y=1;y<nsam[x];y++) {
								fprintf(file_out,"\tfr[%d,%d]:\t%ld",x,y,statistics[0].freq[x][y]);
							}
						}
					}
					else {
						for(x=0;x<npops-1;x++) {
							fprintf(file_out,"\tnsam[%d]:\t%d\tS[%d]:\t%ld",x,nsam[x],x,(long int)statistics[0].S[x]);
							for(y=1;y<=floor(nsam[x]/2);y++) {
								fprintf(file_out,"\tfr[%d,%d]:\t%ld",x,y,statistics[0].freq[x][y]);
							}
						}
					}
					/*fprintf(file_out,"\n");
					 
					 if(npops>1) {
					 for(x=0;x<npops-oo;x++) {
					 fprintf(file_out,"%d\t%ld\t",nsam[x],(long int)statistics[0].S[x]);
					 for(y=1;y<nsam[x];y++) {
					 fprintf(file_out,"%ld\t",statistics[0].freq[x][y]);
					 }
					 }
					 }
					 else {
					 x=0;
					 fprintf(file_out,"%d\t%ld\t",nsam[x],(long int)statistics[0].S[x]);
					 for(y=1;y<=floor(nsam[0]/2);y++) {
					 fprintf(file_out,"%ld\t",statistics[0].freq[0][y]);
					 }
					 }
					 */
					fprintf(file_out,"\n");
				}
				if(output == 1) {
					 /*
					 if(outgroup_presence == 1) {
						np = npops-1;
						for(x=0;x<np;x++) {
							if(nsam[x] > 1) {
								fprintf(file_out,"S[%d]:\t%d\t",x,(int)statistics[0].S[x]);
								fprintf(file_out,"Theta(Wat)[%d]:\t%f\t",x,statistics[0].thetaS[x]);
								fprintf(file_out,"Theta(Taj)[%d]:\t%f\t",x,statistics[0].thetaT[x]);
								fprintf(file_out,"Theta(Fu&Li)[%d]:\t%f\t",x,statistics[0].thetaFL[x]);
								fprintf(file_out,"Theta(Fay&Wu)[%d]:\t%f\t",x,statistics[0].thetaFW[x]);
								fprintf(file_out,"Theta(Zeng)[%d]:\t%f\t",x,statistics[0].thetaL[x]);
								fprintf(file_out,"Theta(Achaz,Wat)[%d]:\t%f\t",x,statistics[0].thetaSA[x]);
								fprintf(file_out,"Theta(Achaz,Taj)[%d]:\t%f\t",x,statistics[0].thetaTA[x]);
								if(statistics[0].thetaTHKY[x] > -10000 && missratio == 0.) fprintf(file_out,"Theta(Taj)HKY[%d]:\t%f\t",x,statistics[0].thetaTHKY[x]*length_al);
								else fprintf(file_out,"Theta(Taj)HKY[%d]:\tNA\t",x);
							}
							else {
								fprintf(file_out,"S[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Wat)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Taj)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Fu&Li)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Fay&Wu)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Zeng)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Achaz,Wat)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Achaz,Taj)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Taj)HKY[%d]:\tNA\t",x);
							}
							if(force_outgroup==0) fprintf(file_out,"Divergence[%d]:\t%f\t",x,statistics[0].K[x]);
							else fprintf(file_out,"Divergence[%d]:\tNA\t",x);
							if(statistics[0].KHKY[x] > -10000 && missratio == 0.) fprintf(file_out,"DivergenceHKY[%d]:\t%f\t",x,statistics[0].KHKY[x]*length_al);
							else fprintf(file_out,"DivergenceHKY[%d]:\tNA\t",x);
						}
					}
					else {
						np = npops-1;
						for(x=0;x<np;x++) {
							if(nsam[x] > 1) {
								fprintf(file_out,"S[%d]:\t%d\t",x,(int)statistics[0].S[x]);
								fprintf(file_out,"Theta(Wat)[%d]:\t%f\t",x,statistics[0].thetaS[x]);
								fprintf(file_out,"Theta(Taj)[%d]:\t%f\t",x,statistics[0].thetaT[x]);
								fprintf(file_out,"Theta(Fu&Li)[%d]:\t%f\t",x,statistics[0].thetaFL[x]);
								fprintf(file_out,"Theta(Achaz,Wat)[%d]:\t%f\t",x,statistics[0].thetaSA[x]);
								fprintf(file_out,"Theta(Achaz,Taj)[%d]:\t%f\t",x,statistics[0].thetaTA[x]);
								if(statistics[0].thetaTHKY[x] > -10000 && missratio == 0.) fprintf(file_out,"Theta(Taj)HKY[%d]:\t%f\t",x,statistics[0].thetaT[x]*length_al);
								else fprintf(file_out,"Theta(Taj)HKY[%d]:\tNA\t",x);
							}
							else {
								fprintf(file_out,"S[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Wat)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Taj)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Fu&Li)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Achaz,Wat)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Achaz,Taj)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Taj)HKY[%d]:\tNA\t",x);
							}
						}
					}
					*//**/
					if(outgroup_presence + force_outgroup) {
						np = npops-1;
						for(x=0;x<np;x++) {
							if(nsam[x] > 1) {
								/*
								fprintf(file_out,"S[%d]:\t%d\t",x,(int)statistics[0].S[x]);
								fprintf(file_out,"Thetaw(Wat)[%d]:\t%f\t",x,statistics[0].thetaS[x]);
								fprintf(file_out,"Thetaw(Taj)[%d]:\t%f\t",x,statistics[0].thetaT[x]);
								*/
								fprintf(file_out,"S[%d]:\t%d\t",x,(int)statistics[0].So[x]);
								fprintf(file_out,"Theta(Wat)[%d]:\t%f\t",x,statistics[0].thetaSo[x]);
								fprintf(file_out,"Theta(Taj)[%d]:\t%f\t",x,statistics[0].thetaTo[x]);
								fprintf(file_out,"Theta(Fu&Li)[%d]:\t%f\t",x,statistics[0].thetaFL[x]);
								fprintf(file_out,"Theta(Fay&Wu)[%d]:\t%f\t",x,statistics[0].thetaFW[x]);
								fprintf(file_out,"Theta(Zeng)[%d]:\t%f\t",x,statistics[0].thetaL[x]);
								fprintf(file_out,"Theta(Achaz,Wat)[%d]:\t%f\t",x,statistics[0].thetaSA[x]);
								fprintf(file_out,"Theta(Achaz,Taj)[%d]:\t%f\t",x,statistics[0].thetaTA[x]);
								if(statistics[0].thetaTHKY[x] > -10000/* && missratio == 0.*/) fprintf(file_out,"Theta/nt(Taj)HKY[%d]:\t%f\t",x,statistics[0].thetaTHKY[x]);
								else fprintf(file_out,"Theta/nt(Taj)HKY[%d]:\tNA\t",x);
								/*if(statistics[0].S[x] > 0) {*/
                                    /*fprintf(file_out,"an_x[%d]:\t%f\t",x,statistics[0].anx[x]);
                                    fprintf(file_out,"bn_x[%d]:\t%f\t",x,statistics[0].bnx[x]);*/
                                    fprintf(file_out,"an_xo[%d]:\t%f\t",x,statistics[0].anxo[x]);
                                    fprintf(file_out,"bn_xo[%d]:\t%f\t",x,statistics[0].bnxo[x]);
								/*}
								else {
                                    fprintf(file_out,"an_x[%d]:\tNA\t",x);
                                    fprintf(file_out,"bn_x[%d]:\tNA\t",x);
                                    fprintf(file_out,"an_xo[%d]:\tNA\t",x);
                                    fprintf(file_out,"bn_xo[%d]:\tNA\t",x);
								}*/
							}
							else {
								fprintf(file_out,"S[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Wat)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Taj)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Fu&Li)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Fay&Wu)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Zeng)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Achaz,Wat)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Achaz,Taj)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta/nt(Taj)HKY[%d]:\tNA\t",x);
								/*fprintf(file_out,"an_x[%d]:\tNA\t",x);
								fprintf(file_out,"bn_x[%d]:\tNA\t",x);*/
								fprintf(file_out,"an_xo[%d]:\tNA\t",x);
								fprintf(file_out,"bn_xo[%d]:\tNA\t",x);
							}
							if(force_outgroup==0) fprintf(file_out,"Divergence[%d]:\t%f\t",x,statistics[0].K[x]);
							else fprintf(file_out,"Divergence[%d]:\tNA\t",x);
							if(statistics[0].KHKY[x] > -10000/* && missratio == 0.*/) fprintf(file_out,"Divergence/nt_HKY[%d]:\t%f\t",x,statistics[0].KHKY[x]);
							else fprintf(file_out,"Divergence/nt_HKY[%d]:\tNA\t",x);
						}
					}
					else {
						np = npops-1;
						for(x=0;x<np;x++) {
							if(nsam[x] > 1) {
								fprintf(file_out,"S[%d]:\t%d\t",x,(int)statistics[0].S[x]);
								fprintf(file_out,"Theta(Wat)[%d]:\t%f\t",x,statistics[0].thetaS[x]);
								fprintf(file_out,"Theta(Taj)[%d]:\t%f\t",x,statistics[0].thetaT[x]);
								fprintf(file_out,"Theta(Fu&Li)[%d]:\t%f\t",x,statistics[0].thetaFL[x]);
								fprintf(file_out,"Theta(Achaz,Wat)[%d]:\t%f\t",x,statistics[0].thetaSA[x]);
								fprintf(file_out,"Theta(Achaz,Taj)[%d]:\t%f\t",x,statistics[0].thetaTA[x]);
								if(statistics[0].thetaTHKY[x] > -10000/* && missratio == 0.*/) fprintf(file_out,"Theta/nt(Taj)HKY[%d]:\t%f\t",x,statistics[0].thetaTHKY[x]);
								else fprintf(file_out,"Theta/nt(Taj)HKY[%d]:\tNA\t",x);
								/*if(statistics[0].S[x] > 0) {*/
									fprintf(file_out,"an_x[%d]:\t%f\t",x,statistics[0].anx[x]/*(double)statistics[0].S[x]/(double)statistics[0].thetaS[x]*/);
									fprintf(file_out,"bn_x[%d]:\t%f\t",x,statistics[0].bnx[x]);
                                    /*fprintf(file_out,"an_xo[%d]:\t%f\t",x,statistics[0].anxo[x]);
                                    fprintf(file_out,"bn_xo[%d]:\t%f\t",x,statistics[0].bnxo[x]);*/
								/*}
								else {
									fprintf(file_out,"an_x[%d]:\tNA\t",x);
									fprintf(file_out,"bn_x[%d]:\tNA\t",x);
									fprintf(file_out,"an_xo[%d]:\tNA\t",x);
									fprintf(file_out,"bn_xo[%d]:\tNA\t",x);
								}*/
							}
							else {
								fprintf(file_out,"S[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Wat)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Taj)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Fu&Li)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Achaz,Wat)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta(Achaz,Taj)[%d]:\tNA\t",x);
								fprintf(file_out,"Theta/nt(Taj)HKY[%d]:\tNA\t",x);
								fprintf(file_out,"an_x[%d]:\tNA\t",x);
								fprintf(file_out,"bn_x[%d]:\tNA\t",x);
								/*fprintf(file_out,"an_xo[%d]:\tNA\t",x);
								fprintf(file_out,"bn_xo[%d]:\tNA\t",x);*/
							}
						}
					}	
					/**/
					np = npops-1;
					for(x=0;x<np;x++) {
						if(include_unknown == 0) {
							if(nsam[x]>1) {
								if(ploidy[0] == '1' && include_unknown == 0) 
									fprintf(file_out,"HapW[%d]:\t%f\tnHap[%d]:\t%d\t",x,statistics[0].hapw[x],x,statistics[0].nhpop[x]);
							}
							else {
								if(ploidy[0] == '1' && include_unknown == 0) 
									fprintf(file_out,"HapW[%d]:\tNA\tnHap[%d]:\t1\t",x,x);
							}
						}
					}
					if(outgroup_presence == 1) {
						np = npops-1;
						for(x=0;x<np;x++) {
							if(statistics[0].Dtaj[x] > -10000)
								fprintf(file_out,"TajimaD[%d]:\t%f\t",x,statistics[0].Dtaj[x]);
							else fprintf(file_out,"TajimaD[%d]:\tNA\t",x);
							if(statistics[0].Dfl[x] > -10000)
								fprintf(file_out,"Fu&LiD[%d]:\t%f\t",x,statistics[0].Dfl[x]);
							else  fprintf(file_out,"Fu&LiD[%d]:\tNA\t",x);
							if(statistics[0].Ffl[x] > -10000)
								fprintf(file_out,"Fu&LiF[%d]:\t%f\t",x,statistics[0].Ffl[x]);
							else fprintf(file_out,"Fu&LiF[%d]:\tNA\t",x);
							if(statistics[0].Hnfw[x] > -10000)
								fprintf(file_out,"Fay&WunormH[%d]:\t%f\t",x,statistics[0].Hnfw[x]);
							else fprintf(file_out,"Fay&WunormH[%d]:\tNA\t",x);
							/*if(statistics[0].thetaT[x]-statistics[0].thetaFW[x] != -10000)
								fprintf(file_out,"Fay&WuH[%d]:\t%f\t",x,statistics[0].thetaT[x]-statistics[0].thetaFW[x]);
							else fprintf(file_out,"Fay&WuH[%d]: NA\t",x);*/
							if(statistics[0].Ez[x] > -10000)
								fprintf(file_out,"ZengE[%d]:\t%f\t",x,statistics[0].Ez[x]);
							else fprintf(file_out,"ZengE[%d]:\tNA\t",x);
                            if(statistics[0].Yach[x] > -10000)
                                fprintf(file_out,"AchazY[%d]:\t%f\t",x,statistics[0].Yach[x]);
                            else fprintf(file_out,"AchazY[%d]:\tNA\t",x);
/**/
                            if(statistics[0].FH[x] > -10000)
                                fprintf(file_out,"FerrettiL[%d]:\t%f\t",x,statistics[0].FH[x]);
                            else fprintf(file_out,"FerrettiL[%d]:\tNA\t",x);
/**/
                            if(ploidy[0] == '1') {
								/*
								if(statistics[0].R2[x] > -10000)
									fprintf(file_out,"R2[%d]:\t%f\t",x,statistics[0].R2[x]);
								else fprintf(file_out,"R2[%d]:\tNA\t",x);
								*/
								if(include_unknown == 0) {
									if(statistics[0].Fs[x] > -10000 && missratio < 1e-6)
										fprintf(file_out,"Fs[%d]:\t%f\t",x,statistics[0].Fs[x]);
									else fprintf(file_out,"Fs[%d]:\tNA\t",x);
								}
							}/*
							else {
								if(statistics[0].R2[x] > -10000
									fprintf(file_out,"R2d[%d]:\t%f\t",x,statistics[0].R2[x]);
								else fprintf(file_out,"R2d[%d]:\tNA\t",x);					
							}*/
							for(xx=0;xx<r2i_ploidies[0];xx++) {
								if(statistics[0].R2p[xx][x] > -10000) {
									if(r2i_ploidies[0]==1 && ploidy[0] == '1') fprintf(file_out,"R2[%d]:\t%f\t",x,statistics[0].R2p[xx][x]);
									if(r2i_ploidies[0]==1 && ploidy[0] == '2') fprintf(file_out,"R2d[%d]:\t%f\t",x,statistics[0].R2p[xx][x]);
									if(r2i_ploidies[0]!=1) fprintf(file_out,"R2p%d[%d]:\t%f\t",r2i_ploidies[xx+1],x,statistics[0].R2p[xx][x]);
								} else {
									if(r2i_ploidies[0]==1 && ploidy[0] == '1') fprintf(file_out,"R2[%d]:\tNA\t",x);
									if(r2i_ploidies[0]==1 && ploidy[0] == '2') fprintf(file_out,"R2d[%d]:\tNA\t",x);
									if(r2i_ploidies[0]!=1) fprintf(file_out,"R2p%d[%d]:\tNA\t",r2i_ploidies[xx+1],x);
									
								}
							}
							
							/*fprintf(file_out,"Fay&WuH[%d]:\t%f\t",x,statistics[0].thetaT[x]-statistics[0].thetaFW[x]);*/
							/*fprintf(file_out,"Fay&WuH2[%d]:\t%f\t",x,2.*(statistics[0].thetaT[x]-statistics[0].thetaL[x]));*/
						}
					}
					else {
						np = npops-1;
						for(x=0;x<np;x++) {
							if(statistics[0].Dtaj[x] > -10000)
								fprintf(file_out,"TajimaD[%d]:\t%f\t",x,statistics[0].Dtaj[x]);
							else fprintf(file_out,"TajimaD[%d]:\tNA\t",x);
							if(statistics[0].Dfl[x] > -10000)
								fprintf(file_out,"Fu&LiD*[%d]:\t%f\t",x,statistics[0].Dfl[x]);
							else  fprintf(file_out,"Fu&LiD*[%d]:\tNA\t",x);
							if(statistics[0].Ffl[x] > -10000)
								fprintf(file_out,"Fu&LiF*[%d]:\t%f\t",x,statistics[0].Ffl[x]);
							else fprintf(file_out,"Fu&LiF*[%d]:\tNA\t",x);
							if(statistics[0].Yach[x] > -10000)
								fprintf(file_out,"AchazY*[%d]:\t%f\t",x,statistics[0].Yach[x]);
							else fprintf(file_out,"AchazY*[%d]:\tNA\t",x);
							if(ploidy[0] == '1') {
								/*
								if(statistics[0].R2[x] > -10000)
									fprintf(file_out,"R2[%d]:\t%f\t",x,statistics[0].R2[x]);
								else fprintf(file_out,"R2[%d]:\tNA\t",x);
								*/
								if(include_unknown == 0) {
									if(statistics[0].Fs[x] > -10000 && missratio < 1e-6)
										fprintf(file_out,"Fs[%d]:\t%f\t",x,statistics[0].Fs[x]);
									else fprintf(file_out,"Fs[%d]:\tNA\t",x);
								}
							}/*
							else {
								if(statistics[0].R2[x] > -10000)
									fprintf(file_out,"R2d[%d]:\t%f\t",x,statistics[0].R2[x]);
								else fprintf(file_out,"R2d[%d]:\tNA\t",x);					
							}*/
							for(xx=0;xx<r2i_ploidies[0];xx++) {
								if(statistics[0].R2p[xx][x] > -10000) {
									if(r2i_ploidies[0]==1 && ploidy[0] == '1') fprintf(file_out,"R2[%d]:\t%f\t",x,statistics[0].R2p[xx][x]);
									if(r2i_ploidies[0]==1 && ploidy[0] == '2') fprintf(file_out,"R2d[%d]:\t%f\t",x,statistics[0].R2p[xx][x]);
									if(r2i_ploidies[0]!=1) fprintf(file_out,"R2p%d[%d]:\t%f\t",r2i_ploidies[xx+1],x,statistics[0].R2p[xx][x]);
								} else {
									if(r2i_ploidies[0]==1 && ploidy[0] == '1') fprintf(file_out,"R2[%d]:\tNA\t",x);
									if(r2i_ploidies[0]==1 && ploidy[0] == '2') fprintf(file_out,"R2d[%d]:\tNA\t",x);
									if(r2i_ploidies[0]!=1) fprintf(file_out,"R2p%d[%d]:\tNA\t",r2i_ploidies[xx+1],x);
									
								}
							}
						}
					}
					
					if(H1frq/* && missratio == 0.*/) {
						if(outgroup_presence == 1) {		
							for(x=0;x<npops-1;x++) {
								if(statistics[0].To_ii[x] > -10000)
									fprintf(file_out,"To_ii[%d]:\t%f\t",x,statistics[0].To_ii[x]);
								else  fprintf(file_out,"To_ii[%d]:\tNA\t",x);
								if(statistics[0].To_00[x] > -10000)
									fprintf(file_out,"To_00[%d]:\t%f\t",x,statistics[0].To_00[x]);
								else  fprintf(file_out,"To_00[%d]:\tNA\t",x);
								if(statistics[0].To_i0[x] > -10000)
									fprintf(file_out,"To_i0[%d]:\t%f\t",x,statistics[0].To_i0[x]);
								else  fprintf(file_out,"To_i0[%d]:\tNA\t",x);
							#if TO_NEW
								if(statistics[0].To_Qc_ii[x] > -10000)
									fprintf(file_out,"To_Qc_ii[%d]:\t%f\t",x,statistics[0].To_Qc_ii[x]);
								else  fprintf(file_out,"To_Qc_ii[%d]:\tNA\t",x);
								if(statistics[0].To_Qw_ii[x] > -10000)
									fprintf(file_out,"To_Qw_ii[%d]:\t%f\t",x,statistics[0].To_Qw_ii[x]);
								else  fprintf(file_out,"To_Qw_ii[%d]:\tNA\t",x);
								if(statistics[0].To_Lc_ii[x] > -10000)
									fprintf(file_out,"To_Lc_ii[%d]:\t%f\t",x,statistics[0].To_Lc_ii[x]);
								else  fprintf(file_out,"To_Lc_ii[%d]:\tNA\t",x);
							#endif
								if(H0frq) {
									if(statistics[0].ToH0_ii[x] > -10000)
										fprintf(file_out,"To_H0_ii[%d]:\t%f\t",x,statistics[0].ToH0_ii[x]);
									else  fprintf(file_out,"To_H0_ii[%d]:\tNA\t",x);
									if(statistics[0].ToH0_00[x] > -10000)
										fprintf(file_out,"To_H0_i0[%d]:\t%f\t",x,statistics[0].ToH0_00[x]);
									else  fprintf(file_out,"To_H0_i0[%d]:\tNA\t",x);
								}
							}
						}
					}
					if(npops>2 && outgroup_presence == 1) {
						for(x=0;x<npops-1;x++) {
							fprintf(file_out,"Sx[%d]:\t%ld\t",x,statistics[0].Sanc[x*4+0]);
							fprintf(file_out,"Sf[%d]:\t%ld\t",x,statistics[0].Sanc[x*4+1]);
							fprintf(file_out,"Sxf[%d,rest]:\t%ld\t",x,statistics[0].Sanc[x*4+2]);
							fprintf(file_out,"Ss[%d]:\t%ld\t",x,statistics[0].Sanc[x*4+3]);
						}
						fprintf(file_out,"Sx[outg]:\t%ld\t",statistics[0].Sanc[(npops-1)*4+0]);
						fprintf(file_out,"Sf[outg]:\t%ld\t",statistics[0].Sanc[(npops-1)*4+1]);
						fprintf(file_out,"Ss[outg,rest]:\t%ld\t",statistics[0].Sanc[(npops-1)*4+2]);
						if(npops > 3) fprintf(file_out,"Ss[rest]:\t%ld\t",statistics[0].Sanc[(npops-1)*4+3]);
						else fprintf(file_out,"Ss[rest]:\tNA\t");

					}
					else {
						if(npops>2 && outgroup_presence == 0) {
							for(x=0;x<npops-1;x++) {
								fprintf(file_out,"Sx[%d]:\t%ld\t",x,statistics[0].Sanc[x*4+0] /*+ statistics[0].Sanc[x*4+2]*/);
								fprintf(file_out,"Sf[%d]:\t%ld\t",x,statistics[0].Sanc[x*4+1]);
								fprintf(file_out,"Ss[%d]:\t%ld\t",x,statistics[0].Sanc[x*4+3]);
							}
							fprintf(file_out,"Ss[rest]:\t%ld\t",statistics[0].Sanc[(npops-1)*4+3]);
						}
						else {
							if(npops > 1) {
								x=0;
								fprintf(file_out,"Sx[%d]:\t%ld\t",x,statistics[0].Sanc[x*4+0]);
								x=1;
								fprintf(file_out,"Sx[%d]:\t%ld\t",x,statistics[0].Sanc[x*4+0]);
								x=0;
								fprintf(file_out,"Sf:\t%ld\t",statistics[0].Sanc[x*4+1]);
								fprintf(file_out,"Ss:\t%ld\t",statistics[0].Sanc[x*4+3]);
							}
							else
								fprintf(file_out,"Sx[%d]:\t%ld\t",x,statistics[0].Sanc[0*4+0]);
						}
					}
					
					/*if(ploidy[0] == '1' && include_unknown==0) {*/
						for(x=0;x<npops-1/*!outgroup_presence*/;x++) {
							if(nsam[x] > 1) {
                                if(statistics[0].mdsd[x] > -10000.)
                                    fprintf(file_out,"MD_SDev[%d]:\t%f\t",x,statistics[0].mdsd[x]);
                                else fprintf(file_out,"MD_SDev[%d]:\tNA\t",x);
                            }
                            else
                                fprintf(file_out,"MD_SDev[%d]:\tNA\t",x);
                            if(nsam[x] > 2) {
								if(statistics[0].mdg1[x] > -10000. && statistics[0].mdsd[x] > -10000.) fprintf(file_out,"MD_Skewness[%d]:\t%f\t",x,statistics[0].mdg1[x]);
								else fprintf(file_out,"MD_Skewness[%d]:\tNA\t",x);
							}
							else {
								fprintf(file_out,"MD_Skewness[%d]:\tNA\t",x);
							}
							if(nsam[x] > 3) {
								if(statistics[0].mdg2[x] > -10000. && statistics[0].mdsd[x] > -10000.) fprintf(file_out,"MD_Kurtosis[%d]:\t%f\t",x,statistics[0].mdg2[x]);
								else fprintf(file_out,"MD_Kurtosis[%d]:\tNA\t",x);
							}
							else
								fprintf(file_out,"MD_Kurtosis[%d]:\tNA\t",x);
						}
					/*}*/
					
					if(npops-1/*!outgroup_presence*/ > 1) {
						if(npops-1/*!outgroup_presence*/ > 2) {
							if(statistics[0].fstALL > -10000) {
								if(niter && piter[0].niteriall > 0 && include_unknown == 0) fprintf(file_out,"Fst:\t%f\t P-value:\t%f\t",statistics[0].fstALL,(double)piter[0].iall/(double)piter[0].niteriall);
								else fprintf(file_out,"Fst:\t%f\t P-value:\tNA\t",statistics[0].fstALL);
							}	
							else {
								if(niter && piter[0].niteriall > 0 && include_unknown == 0) fprintf(file_out,"Fst:\tNA\t P-value:\tNA\t");
								else fprintf(file_out,"Fst:\tNA\t P-value:\tNA\t");
							}
										
							for(x=0;x<npops-1;x++) {
								if(statistics[0].fst1all[x] > -10000) {
									if(niter && piter[0].niteri1[x] > 0 && include_unknown == 0) fprintf(file_out,"Fst1[%d,rest]:\t%f\t P-value:\t%f\t",x,statistics[0].fst1all[x],(double)piter[0].i1[x]/(double)piter[0].niteri1[x]);
									else fprintf(file_out,"Fst1[%d,rest]:\t%f\t P-value:\tNA\t",x,statistics[0].fst1all[x]);
								}
								else {
									if(niter && piter[0].niteri1[x] > 0 && include_unknown == 0)
                                        fprintf(file_out,"Fst1[%d,rest]:\tNA\t P-value:\tNA\t",x);
									else fprintf(file_out,"Fst1[%d,rest]:\tNA\t P-value:\tNA\t",x);
								}
							}
						}
						
						z = 0;
						for(x=0;x<npops-1;x++) {
							for(y=x+1;y<npops-0;y++) {
								if(y==npops-1) {z++;continue;}
								if(statistics[0].fst[z] > -10000) {
									if(niter && piter[0].niteri[z] > 0 && include_unknown == 0) fprintf(file_out,"Fst[%d,%d]:\t%f\t P-value:\t%f\t",x,y,statistics[0].fst[z],(double)piter[0].i[z]/(double)piter[0].niteri[z]);
									else fprintf(file_out,"Fst[%d,%d]:\t%f\t P-value:\tNA\t",x,y,statistics[0].fst[z]);
								}
								else {
									if(niter && piter[0].niteri[z] > 0 && include_unknown == 0) fprintf(file_out,"Fst[%d,%d]:\tNA\t P-value:\tNA\t",x,y);
									else fprintf(file_out,"Fst[%d,%d]:\tNA\t P-value:\tNA\t",x,y);
								}
								z++;
							}
						}
					}
					if(include_unknown == 0) {
						if(npops-1/*!outgroup_presence*/ > 1) {
							if(npops-1/*!outgroup_presence*/ > 2) {
								z=0;
								for(x=0;x<npops-1;x++) {
									for(y=x+1;y<npops-0;y++) {
										if(y==npops-1) {z++;continue;}
										if(statistics[0].fstHKY[z] > -10000) {
											fprintf(file_out,"FstHKY[%d,%d]:\t%f\tPiWHKY[%d]:\t%f\tPiWHKY[%d]:\t%f\tPiAHKY[%d,%d]:\t%f\tPiTHKY[%d,%d]:\t%f\t",
													x,y,statistics[0].fstHKY[z],x,statistics[0].piwHKY[x],y,statistics[0].piwHKY[y],x,y,statistics[0].piaHKY[z],x,y,statistics[0].piTHKY[z]);
										}
										else {
											fprintf(file_out,"FstHKY[%d,%d]:\tNA\t",x,y);
											if(statistics[0].piwHKY[x] > -10000) {
												fprintf(file_out,"PiWHKY[%d]:\t%f\t",x,statistics[0].piwHKY[x]);
											}
											else fprintf(file_out,"PiWHKY[%d]:\tNA\t",x);
											if(statistics[0].piwHKY[y] > -10000) {
												fprintf(file_out,"PiWHKY[%d]:\t%f\t",y,statistics[0].piwHKY[y]);
											}
											else fprintf(file_out,"PiWHKY[%d]:\tNA\t",y);
											if(statistics[0].piaHKY[z] > -10000) {
												fprintf(file_out,"PiAHKY[%d,%d]:\t%f\tPiTHKY[%d,%d]:\t%f\t",x,y,statistics[0].piaHKY[z],x,y,statistics[0].piTHKY[z]);
											}
											else fprintf(file_out,"PiAHKY[%d,%d]:\tNA\tPiTHKY[%d,%d]:\tNA\t",x,y,x,y);
										}
										z++;
									}
								}
							}
						}
					}
					if(ploidy[0] == '1' && include_unknown == 0) {
						for(x=0;x<npops-1/*!outgroup_presence*/;x++) {
							if(nsam[x]>1) fprintf(file_out,"HapW[%d]:\t%f\t",x,statistics[0].hapw[x]);
							else fprintf(file_out,"HapW[%d]:\tNA\t",x);
						}
						if(npops-1/*!outgroup_presence*/ > 1) {
							z=0;
							for(x=0;x<npops-1;x++) {
								for(y=x+1;y<npops-0;y++) {
									if(y==npops-1) {z++;continue;}
									fprintf(file_out,"HapA[%d,%d]:\t%f\tHapT[%d,%d]:\t%f\t",x,y,statistics[0].hapa[z],x,y,statistics[0].hapT[z]);
									z++;
								}
							}
							if(npops-1/*!outgroup_presence*/ > 2) {
								if(statistics[0].fsthALL > -10000) {
									if(niter && piter[0].niterihall > 0 && include_unknown == 0) fprintf(file_out,"FstH:\t%f\t P-value:\t%f\t",statistics[0].fsthALL,(double)piter[0].ihall/(double)piter[0].niterihall);
									else fprintf(file_out,"FstH:\t%f\t P-value:\tNA\t",statistics[0].fsthALL);
								}
								else {
									if(niter && piter[0].niterihall > 0) fprintf(file_out,"FstH:\tNA\t P-value:\tNA\t");
									else fprintf(file_out,"FstH:\tNA\t P-value:\tNA\t");
								}

								for(x=0;x<npops-1;x++) {
									if(statistics[0].fsth1all[x] > -10000 && include_unknown == 0) {
										if(niter && piter[0].niterih1[x]) fprintf(file_out,"FstH1[%d,rest]:\t%f\t P-value:\t%f\t",x,statistics[0].fsth1all[x],(double)piter[0].ih1[x]/(double)piter[0].niterih1[x]);
										else fprintf(file_out,"FstH1[%d,rest]:\t%f\t P-value:\tNA\t",x,statistics[0].fsth1all[x]);
									}
									else {
										if(niter && piter[0].niterih1[x] > 0 && include_unknown == 0) fprintf(file_out,"FstH1[%d,rest]:\tNA\t P-value:\tNA\t",x);
										else fprintf(file_out,"FstH1[%d,rest]:\tNA\t P-value:\tNA\t",x);
									}
								}
								/*
								if(statistics[0].GstALL > -10000) {
									if(niter && piter[0].niterighall > 0) fprintf(file_out,"Gst':\t%f\t P-value:\t%f\t",statistics[0].GstALL,(double)piter[0].ighall/(double)piter[0].niterighall);
									else fprintf(file_out,"Gst':\t%f\t P-value:\tNA\t",statistics[0].GstALL);
								}
								else {
									if(niter && piter[0].niterighall > 0) fprintf(file_out,"Gst':\tNA\t P-value:\tNA\t");
									else fprintf(file_out,"Gst':\tNA\t P-value:\tNA\t");
								}
								 */
							}
							z = 0;
							for(x=0;x<npops-1;x++) {
								for(y=x+1;y<npops-0;y++) {
									if(y==npops-1 && !outgroup_presence) {z++;continue;}
									if(statistics[0].fsth[z] > -10000) {
										if(niter && piter[0].niterih[z] > 0 && include_unknown == 0) fprintf(file_out,"FstH[%d,%d]:\t%f\t P-value:\t%f\t",x,y,statistics[0].fsth[z],(double)piter[0].ih[z]/(double)piter[0].niterih[z]);
										else fprintf(file_out,"FstH[%d,%d]:\t%f\t P-value:\tNA\t",x,y,statistics[0].fsth[z]);
									}
									else {
										if(niter && piter[0].niterih[z] > 0 && include_unknown == 0) fprintf(file_out,"FstH[%d,%d]:\tNA\t P-value:\tNA\t",x,y);
										else fprintf(file_out,"FstH[%d,%d]:\tNA\t P-value:\tNA\t",x,y);
									}
									z++;
								}
							}
							/*
							z = 0;
							for(x=0;x<npops-1;x++) {
								for(y=x+1;y<npops;y++) {
									if(y==npops-1 && !outgroup_presence) {z++;continue;}
									if(statistics[0].Gst[z] > -10000) {
										if(niter && piter[0].niterigh[z] > 0) fprintf(file_out,"Gst'[%d,%d]:\t%f\t P-value:\t%f\t",x,y,statistics[0].Gst[z],(double)piter[0].igh[z]/(double)piter[0].niterigh[z]);
										else fprintf(file_out,"Gst'[%d,%d]:\t%f\t P-value:\tNA\t",x,y,statistics[0].Gst[z]);
									}
									else {
										if(niter && piter[0].niterigh[z] > 0) fprintf(file_out,"Gst'[%d,%d]:\tNA\t P-value:\tNA\t",x,y);
										else fprintf(file_out,"Gst'[%d,%d]:\tNA\t P-value:\tNA\t",x,y);
									}
									z++;
								}
							}
							 */
						}
					}
					if(include_unknown == 0) {/*
						for(x=0;x<npops-1;x++) {
							if(nsam[x]>1) fprintf(file_out,"PiW[%d]:\t%f\t",x,statistics[0].piw[x]);
							else fprintf(file_out,"PiW[%d]:\tNA\t",x);
						}*/
						if(npops-1/*!outgroup_presence*/ > 1) {
							z=0;
							for(x=0;x<npops-1;x++) {
								for(y=x+1;y<npops-0;y++) {
									if(y==npops-1) {z++;continue;}
									fprintf(file_out,"PiA[%d,%d]:\t%f\tPiT[%d,%d]:\t%f\t",x,y,statistics[0].pia[z],x,y,statistics[0].piT[z]);
									z++;
								}
							}
						}
					}
					else {/*
						for(x=0;x<npops-1;x++) {
							if(nsam[x]>1) fprintf(file_out,"PiW[%d]:\t%f\t",x,statistics[0].piw[x]);
							else fprintf(file_out,"PiW[%d]:\tNA\t",x);
						}*/
						if(npops-1/*!outgroup_presence*/ > 1) {
							z=0;
							for(x=0;x<npops-1;x++) {
								for(y=x+1;y<npops-0;y++) {
									if(y==npops-1) {z++;continue;}
									fprintf(file_out,"PiA[%d,%d]:\t%f\tPiT[%d,%d]:\t%f\t",x,y,statistics[0].pia[z],x,y,statistics[0].piT[z]);
									z++;
								}
							}
                            for(x=0;x<npops-1;x++) {
                                for(y=x+1;y<npops-0;y++) {
                                    if(y==npops-1) {continue;}
                                    if(outgroup_presence+force_outgroup)
                                        fprintf(file_out,"len[%d,%d]:\t%f\t",x,y,statistics[0].lengthamng_outg[x][y]);
                                    else
                                        fprintf(file_out,"len[%d,%d]:\t%f\t",x,y,statistics[0].lengthamng[x][y]);
                                }
                            }
						}
					}
					if(outgroup_presence==1) {		
						for(x=0;x<npops-1/*!outgroup_presence*/;x++) {
							for(y=1;y<nsam[x];y++) {
								fprintf(file_out,"fr[%d,%d]:\t%ld\t",x,y,statistics[0].freq[x][y]);
							}
						}
					}
					else {
						for(x=0;x<npops-1/*!outgroup_presence*/;x++) {
							for(y=1;y<=floor(nsam[x]/2);y++) {
								fprintf(file_out,"fr[%d,%d]:\t%ld\t",x,y,statistics[0].freq[0][y]);
							}
						}
					}
					/*no include because the number of columns is very long*//*
					if(ploidy[0] == '1' && include_unknown == 0) {
						for(y=0;y<statistics[0].nh-!outgroup_presence;y++) {
							for(x=0;x<npops-1;x++) {
								fprintf(file_out,"frHap[%d,hap%02d]:\t%ld\t",x,y,statistics[0].freqh[x][y]);
							}
						}
						for(y=statistics[0].nh-!outgroup_presence;y<sumnsam;y++) {
							for(x=0;x<npops-1;x++) {
								fprintf(file_out,"frHap[%d,hap%02d]:\t0\t",x,y);
							}
						}
					}
					*/
					fprintf(file_out,"\n");
				}
			}
		}
	}
	/*fflush(file_out);*/
    free(initsq1);
	return (1);
}

