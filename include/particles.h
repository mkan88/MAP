/*********************
* Date of creation 09/10/2017
* Author: Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
* Other Authors: N/A
**************************************
* History
*
* 09/10/2017 -Created particl_read_in prototype. Also created the struct
*             particleVariables to store the varibles in use.
**************************************/

#ifndef _PARTICLES_H
#define _PARTICLES_H

#include <stdio.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "initial_finalisation.h"

#ifdef __cplusplus
  extern "C" {
#endif

typedef struct
{
    double x;
    double y;
    double z;
    double alpha;
    double beta;
    double gamma;
}particleVariables;

int  particle_read_in(particleVariables **particles);

int generate_particle_data(int numberOfParticles, particleVariables **particles, gsl_rng *tSeed, double radius, double xMax, double yMax, double zMax);

double *generalised_coordinate_creation(int numberOfParticles, particleVariables *particles);

#ifdef __cplusplus
  }
#endif

#endif // _PARTICLES_H
