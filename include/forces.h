#ifndef _FORCES_H
#define _FORCES_H

#include <stdlib.h>

#include "particles.h"

typedef struct
{
	double mag; // Magnitude
	double alpha; // Alpha Euler angle
	double beta; // Beta Euler angle
}field_t;

void force_torque_summation(double *additionalForces,double *generalisedCoordinates, int numberOfCells, int *forceList, int numberOfForces,	environmentVariables conditions, field_t drivingField);

enum forces_available
{
    NONE ,
    GRAVITY ,
    VAN_DER_WAALS ,
    EXP_REPULSION ,
	ALIGN_TORQUE ,
	DRIVING_FIELD,
	POLAR_DRIVING_FORCE,
	VISECK_ALIGN_TORQUE
};


#endif //_FORCES_H
