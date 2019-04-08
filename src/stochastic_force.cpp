
/*********************
* Date of creation 17/10/2017
* Author: Oliver Hinds
* Contact:
**************************************
* History
**************************************/

#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include "stochastic_force.h"
#include "maths_functions.h"

#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "cusolverDn.h"
#include "helper_cuda.h"
#include "helper_cusolver.h"


extern "C" void stochastic_displacement_creation(int numberOfParticles, double *stochasticWeighting, double *stochasticDisplacement, gsl_rng *rndarray[], 
double *rndNumArray, double timestep);


extern int gDebug;
extern int gNumOfthreads;

void stochastic_displacement_creation(int numberOfParticles, double *stochasticWeighting, double *stochasticDisplacement, gsl_rng *rndarray[],
double *rndNumArray, double timestep){


	/*// HOME-MADE METHOD
	int N = 6 * numberOfParticles; // array is a linearized (6*N) by (6*N) array
    #pragma omp parallel
    {
        gsl_rng *tSeed = rndarray[omp_get_thread_num()];

		#pragma omp single
    	for (int k = 0; k < N; k++) // iterates over diagonals
    	{

    		stochasticWeighting[k*(N+1)] = sqrt(stochasticWeighting[k*(N+1)]); // square roots diagonals

    		for (int i = k + 1; i < N; i++) // iterates over the elements in the column below the diagonal
    		{
    			stochasticWeighting[i*N + k] = stochasticWeighting[i*N + k] / stochasticWeighting[k*(N+1)];
    		}
            //
    		// iterate over lower triangle subtended by (k,k) element
            //
    		for (int j = k + 1; j < N; j++) // iterates over the columns j>k
    		{
    			for (int i = j; i < N; i++) // iterates over the rows i>k
    			{
    				stochasticWeighting[i*N + j] = stochasticWeighting[i*N + j] - (stochasticWeighting[i*N + k] * stochasticWeighting[j*N + k]);
    			}
    		}
        }
        //
        // iterate over upper triangle and set all values to zero
        //
        #pragma omp for
    	for (int k = 0; k < N; k++)
    	{
    		for (int j = k + 1; j < N; j++)
    		{
    			stochasticWeighting[k*N+j] = 0;
    		}
    	}

        #pragma omp for
    	for (int i = 0; i < N; i++)
    	{
            double ran_num;
    		stochasticDisplacement[i] = 0;
    		for (int j = 0; j < N; j++)
    		{
    			ran_num = gsl_ran_gaussian(tSeed, 1);
    			stochasticDisplacement[i] += stochasticWeighting[i*N + j] * ran_num * sqrt(2*timestep);
    		}
    		//if (gDebug == 1) printf("%+1.5e\n\n", stochasticDisplacement[i]);
    	}
    }
/*
*/

	// LITERATURE METHOD (Cholesky-Crout)
	int N = 6 * numberOfParticles; // array is a linearized (6*N) by (6*N) array
	double sum;
	int i, j, k;
	int cutoff = N/2;



	// copy contents from input matrix to output matrix
	//  this is done to simplify the code and potential testing
	//  and should not be counted as algorithm time
	//for (i = 0; i < dimensionSize * dimensionSize; i++)
		//L[i] = A[i];

	// for (j = 0; j < N; j++)
	// {
	// 	sum = 0;
	// 	for (k = 0; k < j; k++)
	// 	{
	// 		sum += stochasticWeighting[j * N + k] * stochasticWeighting[j * N + k];
	// 	}
	// 	stochasticWeighting[j * N + j] = sqrt(stochasticWeighting[j * N + j] - sum);

	// 	#pragma omp parallel for private(i,k,sum) shared (stochasticWeighting,j) schedule(static) if (j < N - cutoff)
	// 	for (i = j + 1; i < N; i++)
	// 	{
	// 		sum = 0;
	// 		for (k = 0; k < j; k++)
	// 		{
	// 			sum += stochasticWeighting[i * N + k] * stochasticWeighting[j * N + k];
	// 		}
	// 		stochasticWeighting[i * N + j] = (1.0 / stochasticWeighting[j * N + j] * (stochasticWeighting[i * N + j] - sum));
	// 	}
	// }

	// CUDA method
	double * stochasticWeightingDevice = NULL;	// Copy of stochasticWeighting
	cusolverDnHandle_t handle = NULL;
    cudaStream_t stream = NULL;
    int bufferSize = 0;
    int *info = NULL;
    double *buffer = NULL;
    double *stochasticWeightingWorkspace = NULL;	//Copy of stochasticWeighting being decomposed
    int h_info = 0;
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_UPPER;

	checkCudaErrors(cusolverDnCreate(&handle));
    checkCudaErrors(cudaStreamCreate(&stream));

    checkCudaErrors(cusolverDnSetStream(handle, stream));

    checkCudaErrors(cudaMalloc((void **)&stochasticWeightingDevice, sizeof(double)*N*N));

    checkCudaErrors(cudaMemcpy(stochasticWeightingDevice, stochasticWeighting, sizeof(double)*N*N, cudaMemcpyHostToDevice));

    //cuSolverDn_LinearSolver->linearSolverCHOL
    checkCudaErrors(cusolverDnDpotrf_bufferSize(handle, uplo, N, (double*)stochasticWeightingDevice, N, &bufferSize));

    checkCudaErrors(cudaMalloc(&info, sizeof(int)));
    checkCudaErrors(cudaMalloc(&buffer, sizeof(double)*bufferSize));
    checkCudaErrors(cudaMalloc(&stochasticWeightingWorkspace, sizeof(double)*N*N));

    checkCudaErrors(cudaMemcpy(stochasticWeightingWorkspace, stochasticWeightingDevice, sizeof(double)*N*N, cudaMemcpyDeviceToDevice));
    checkCudaErrors(cudaMemset(info, 0, sizeof(int)));

    checkCudaErrors(cusolverDnDpotrf(handle, uplo, N, stochasticWeightingWorkspace, N, buffer, bufferSize, info));

    checkCudaErrors(cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost));

    if ( 0 != h_info ){
        fprintf(stderr, "Error: Cholesky factorization failed\n");
    }

    checkCudaErrors(cudaMemcpy(stochasticWeighting, stochasticWeightingWorkspace, sizeof(double)*N*N, cudaMemcpyDeviceToHost));

    if (info  ) { checkCudaErrors(cudaFree(info)); }
    if (buffer) { checkCudaErrors(cudaFree(buffer)); }
    if (stochasticWeightingWorkspace     ) { checkCudaErrors(cudaFree(stochasticWeightingWorkspace)); }

    if (handle) { checkCudaErrors(cusolverDnDestroy(handle)); }
    // if (cublasHandle) { checkCudaErrors(cublasDestroy(cublasHandle)); }
    if (stream) { checkCudaErrors(cudaStreamDestroy(stream)); }

    if (stochasticWeightingDevice) { checkCudaErrors(cudaFree(stochasticWeightingDevice)); }

    // cudaDeviceReset();


	// reset upper triangle
	//  this is done to simplify the code and potential testing
	//  and should not be counted as algorithm time
	#pragma omp parallel for
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < i; j++) {
			stochasticWeighting[j * N + i] = 0;
		}
	}
	
	#pragma omp parallel
	{
		gsl_rng *tSeed = rndarray[omp_get_thread_num()];
		#pragma omp for
		for (i = 0; i < N; i++)
		{
			rndNumArray[i] = gsl_ran_gaussian(tSeed, 1);
		}
	}


	#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		stochasticDisplacement[i] = 0;
		for (int j = 0; j < N; j++)
		{
			stochasticDisplacement[i] += stochasticWeighting[i*N + j] * rndNumArray[j] * sqrt(2*timestep);
		}
		//if (gDebug == 1) printf("%+1.5e\n\n", stochasticDisplacement[i]);
	}


/*
	// GSL METHOD
	int N = 6 * numberOfParticles; // array is a linearized (6*N) by (6*N) array

	gsl_matrix *test = gsl_matrix_alloc(N, N);

	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			gsl_matrix_set(test, i, j, stochasticWeighting[i*N + j]);

	gsl_linalg_cholesky_decomp1(test);

	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			stochasticWeighting[i*N + j] = gsl_matrix_get(test, i, j);

	#pragma omp for
	for (int i = 0; i < N; i++)
	{
		gsl_rng *tSeed = rndarray[omp_get_thread_num()];
		double ran_num;
		stochasticDisplacement[i] = 0;
		for (int j = 0; j < N; j++)
		{
			ran_num = gsl_ran_gaussian(tSeed, 1);
			stochasticDisplacement[i] += stochasticWeighting[i*N + j] * ran_num * sqrt(2*timestep);
		}
		//if (gDebug == 1) printf("%+1.5e\n\n", stochasticDisplacement[i]);
	}*/
}
