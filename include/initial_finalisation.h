#ifndef _INITIAL_FINALISATION_H
#define _INITIAL_FINALISATION_H

#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <stdarg.h>
#include "particles.h"

#ifdef __cplusplus
  extern "C" {
#endif

typedef struct
{
	double temperature; // K
    double viscosity; //N m^-2 s
    double radius; // m
    double currentTime;
    double deltaTime; // Seconds
    double endTime; // Seconds
	double mass; // kg
    double xMax; //m
    double yMax; //m
    double zMax; //m
    int numberOfParticles;
	int fileNum;
	double drivingForceMagnitude; // N
}environmentVariables;

int cmd_line_read_in(int argc, char *argv[], environmentVariables *conditions);

void boilerplate_variables(environmentVariables *conditions);

gsl_rng** rand_array_allocation();

double* generalised_coordinate_initialisation(environmentVariables *conditions, gsl_rng *rndarray[]);

void free_memory(int listSize, ...);

void initialise_gpu();

void finalise_gpu();

#ifdef __cplusplus
  }
#endif

#endif //_INITIAL_FINALISATION_H
