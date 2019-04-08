/*********************
* Date of creation 17/10/2017
* Author: Oliver Hinds
* Contact:
**************************************
* History
*
**************************************/

#ifndef _STOCHASTIC_FORCE_H
#define _STOCHASTIC_FORCE_H

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


#ifdef __cplusplus
  extern "C" {
#endif

void stochastic_displacement_creation(int numberOfParticles, double *stochasticWeighting, double *stochasticDisplacement, gsl_rng *rndarray[],double *rndNumArray, double timestep);

#ifdef __cplusplus
  }
#endif

#endif // _STOCHASTIC_FORCE_H
