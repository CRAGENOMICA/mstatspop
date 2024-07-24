/*
 *  tn93.c
 *  mstatspop
 *
 *  Created by Sebastian Ramos-Onsins on 26/02/2010.
 *  Copyright 2010 CRAG. All rights reserved.
 *
 */


#include "tn93.h"

/*P1 is the A/G, P2 is T/C and Q is transversions.*/
double tn93(double gT,double gC,double gG,double gA, double P1, double P2, double Q,double length)
{
	double d1,d2,d3,gR,gY,sum,d;
	
	sum = gT+gC+gG+gA;
	gT 	= gT/sum;
	gC 	= gC/sum;
	gG 	= gG/sum;
	gA 	= gA/sum;
	P1 	= P1/length;
	P2 	= P2/length;
	Q  	= Q/length;
	
	gR = gA + gG;
	gY = gT + gC;
	
	if(P1+P2+Q == 0) return 0.;
	if(gR/(2.*gA*gG) * P1 - 1./(2.*gR) * Q >= 1.) return -10000;
	if(gY/(2.*gT*gC) * P2 - 1./(2.*gY) * Q >= 1.) return -10000;
	if(1./(2.*gR*gY) * Q >= 1.) return -10000;
	
	d1 = -2.*gA*gG/gR * log(1. - gR/(2.*gA*gG) * P1 - 1./(2.*gR) * Q);
	d2 = -2.*gT*gC/gY * log(1. - gY/(2.*gT*gC) * P2 - 1./(2.*gY) * Q);
	d3 = -2.*(gR*gY - gA*gG*gY/gR - gT*gC*gR/gY) * log(1. - 1./(2.*gR*gY) * Q);
	
	d = d1 + d2 + d3;
	
	if(d < 0.) return -10000;
	return d;
}
