#ifndef _DIFFUSIONMATRIX_H
#define _DIFFUSIONMATRIX_H

#include <stdio.h>

#include "particles.h"

void diffusion_matrix_creation(int numberOfParticles, double *diffusionMatrix, double *stochasticWeighting, double *generalisedCoordinates, environmentVariables *conditions);

void translational_tensor_creation(double *tempMatrix, double *generalisedCoordinates, environmentVariables *conditions, int i, int j);

void rotational_tensor_creation(double *tempMatrix, double *generalisedCoordinates, environmentVariables *conditions, int i, int j);

void translation_rotation_coupling_tensor_creation(double *tempMatrix, double *generalisedCoordinates, environmentVariables *conditions, int i, int j);

#endif //_DIFFUSIONMATRIX_H
