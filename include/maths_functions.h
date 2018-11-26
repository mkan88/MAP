/*********************
* Date of creation 09/10/2017
* Authors:Oliver Hinds, Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
**************************************
* History
*
**************************************/

#ifndef _MATHS_FUNCTIONS_H
#define _MATHS_FUNCTIONS_H

#include <stdlib.h>


int kronecker_delta(int i, int k);

int levi_civita_density(int i, int j);

int linear_index_from_coordinates(int x_Max,int y_Max, int x, int y, int z);

void coordinates_from_linear_index(int location,int x_Max,int y_Max, int *x, int *y, int *z);


double randSign(long int tSeed);

float ran1(long *idum);

float gaussdev(long *idum);


#endif  // _MATHS_FUNCTIONS_H
