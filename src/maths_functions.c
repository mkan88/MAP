/*********************
* Date of creation 09/10/2017
* Author: Oliver Hinds, Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
**************************************
* History
*
**************************************/

#include <stdlib.h>
#include <math.h>

#include "maths_functions.h"

//
// Performs basic kronecker delta
//

int kronecker_delta(int i, int j)
{
    if(i == j)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

//
// Performs the levi_civita_density permutation
//

int levi_civita_density(int i, int j)
{
    if(i == j)
    {
        return 0;
    }
    else if( (i == 0 && j == 1) || (i == 1 && j == 2) || (i == 2 && j == 0))
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

int linear_index_from_coordinates(int x_Max,int y_Max, int x, int y, int z)
{
    return x+x_Max*y+x_Max*y_Max*z;
}

void coordinates_from_linear_index(int location,int x_Max,int y_Max, int *x, int *y, int *z)
{
     *x = location%x_Max;
     location=(location - *x)/x_Max;
     *y = location%y_Max;
     location=(location-*y)/y_Max;
     *z = location;
}




double randSign(long int tSeed) // Randomly generate -1 or +1
{
	if (ran1(&tSeed)>0.5)
	{
		return 1.0;
	} else
	{
		return -1.0;
	}
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX



float gaussdev(long *idum)
{
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if (*idum < 0)
    iset=0;
	if (iset == 0)
	{
		do
		{
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1*v2*v2;
		}
		while (rsq>=1.0 || rsq == 0.0);
		fac = sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	}
	else
	{
		iset=0;
		return gset;
	}
}
