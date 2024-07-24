

#include "ran1.h"
#include "util.h"

/**< todo: Que hace una variable global??? */
static long int idum = 0;

void init_seed1(long int seed)
{
    idum = -seed;
}

/* Random number generator ran1 from Computers in Physics */
/* Volume 6 No. 5, 1992, 522-524, Press and Teukolsky */
/* To generate real random numbers 0.0-1.0 */
/* Should be seeded with a negative integer */



double ran1(void)
{
	long int j,k;
	static long int iv[NTAB],iy=0;
	static double NDIV = (double)1.0/((double)1.0+((double)IM-(double)1.0)/(double)NTAB);
	static double RNMX = ((double)1.0-(double)EPS);
	static double AM = ((double)1.0/(double)IM);
	double ran;

	/*double precisionv(double,int);*/

	if( (idum <= 0) || (iy == 0) )
	{
		idum = MAX(-idum,idum);
		for( j=NTAB+7; j>=0; j-- )
		{
			k 		= (long int)((double)idum/(double)IQ);
			idum 	= IA*(idum-k*IQ)-IR*k;
			if(idum < 0)
				idum += IM;
			if(j < NTAB)
				iv[j] = idum;
		}
		iy = iv[0];
	}
	k 		= (long int)((double)idum/(double)IQ);
	idum 	= IA*(idum-k*IQ)-IR*k;

	if(idum<0)
		idum += IM;

	j 		= (long int)((double)iy*NDIV);
	iy 		= iv[j];
	iv[j] 	= idum;
	ran 	= MIN((double)AM*(double)iy,(double)RNMX);
	
	/*precision of 1e-7*/
	/*
	ran = precisionv(ran,(int)7);
	*/
	return (double)ran;
}
