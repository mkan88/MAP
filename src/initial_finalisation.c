
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <stdarg.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#include "initial_finalisation.h"

extern int gDebug;
extern int gSerial;
extern int gNumOfthreads;
extern double gPi;


//
// Read in cmd line arguments and handle them
//
int cmd_line_read_in(int argc, char *argv[], environmentVariables *conditions)
{
	if (argc > 1)
	{
		for (int i=1; i<argc; i++)
		{
			if(strstr(argv[i],"-debug") != NULL) gDebug = 1;
			else if(strstr(argv[i],"-serial") != NULL)	gSerial = 1;

			else if (strstr(argv[i],"-num") != NULL || strstr(argv[i],"-n") != NULL)
			{
				if (sscanf(argv[i+1],"%d", &conditions->numberOfParticles) != 1)
				{
					printf("Invalid number of particles\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-temp") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &conditions->temperature) != 1)
				{
					printf("Invalid temperature\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-drivemag") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &conditions->drivingForceMagnitude) != 1)
				{
					printf("Invalid polar driving force magnitude\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-filenum") != NULL)

			{
				if (sscanf(argv[i+1],"%d", &conditions->fileNum) != 1)
				{
					printf("Invalid file number\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-threads") != NULL)
			{
				if (sscanf(argv[i+1],"%d", &gNumOfthreads) != 1)
				{
					printf("Invalid number of threads\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-cube") != NULL)
			{
				double temp_num;
				if (strstr(argv[i+1],"R") != NULL)
				{
					if (sscanf(argv[i+1],"%lfR", &temp_num) != 2)
					{
						conditions->xMax = temp_num * conditions->radius;
						conditions->yMax = temp_num * conditions->radius;
						conditions->zMax = temp_num * conditions->radius;
					}
					else
					{
						printf("Invalid maximum dimension values\n");
						return -1;
					}
				}
				else if (sscanf(argv[i+1],"%lf", &temp_num) == 1)
				{
					conditions->xMax = temp_num;
					conditions->yMax = temp_num;
					conditions->zMax = temp_num;
				}
				else
				{
					printf("Invalid maximum dimension values\n");
					return -1;
				}

				printf("Set dimensions to %1.1em\n", conditions->xMax);

			}
			else if (strstr(argv[i],"-x") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &conditions->xMax) != 1)
				{
					printf("Invalid maximum x-dimension value\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-y") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &conditions->yMax) != 1)
				{
					printf("Invalid maximum y-dimension value\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-z") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &conditions->zMax) != 1)
				{
					printf("Invalid maximum z-dimension value\n");
					return -1;
				}
			}
			else if (strstr(argv[i],"-t") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &conditions->endTime) != 1)
				{
					printf("Invalid duration\n");
					return -1;
				}
				else printf("Simulation duration set to %1.1e seconds\n", conditions->endTime);
			}
			else if (strstr(argv[i],"-dt") != NULL)
			{
				if (sscanf(argv[i+1],"%lf", &conditions->deltaTime) != 1)
				{
					printf("Invalid timestep\n");
					return -1;
				}
				else printf("Simulation timestep set to %1.1e seconds\n", conditions->deltaTime);
			}
		}
		if (gDebug == 1 && gSerial == 1) printf("Debug & serial modes active\n");
		else if (gDebug == 1 && gSerial == 0) printf("Debug mode active\n");
		else if (gDebug == 0 && gSerial == 1) printf("Serial mode active\n");

	}
	return 1;
}

//
// Set default values of environmental variables
//
void boilerplate_variables(environmentVariables *conditions)
{
    conditions->temperature = 298; // K
	conditions->viscosity = 8.9E-4; //N m^-2 s
	conditions->radius = 50E-9; // m
	conditions->currentTime = 0; // Seconds
	conditions->deltaTime = 1E-7; // Seconds
	conditions->endTime = 1E-3; // Seconds
	conditions->mass = (4/3) * gPi * pow(conditions->radius,3)*19320; // kg - density of gold
	conditions->xMax = 1E-6;
	conditions->yMax = 1E-6;
	conditions->zMax = 1E-6;
	conditions->numberOfParticles = 0;
	conditions->fileNum = 0;
	conditions->drivingForceMagnitude=1E-13; //N,  set to 0.1pN, 1pN too large


	gNumOfthreads = omp_get_max_threads();
}

//
// Allocate an array of random number generators
//
gsl_rng** rand_array_allocation()
{
	unsigned int seed;
	gsl_rng **rndarray = calloc(gNumOfthreads,sizeof(*rndarray));
	for(int i = 0; i<gNumOfthreads; i++)
	{
		seed =(i+1) * (unsigned int)time(NULL);
		rndarray[i] = gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(rndarray[i], seed);
	}
	return rndarray;
}

//
// Generate generalised coordinate data from either file or randomly
//
double* generalised_coordinate_initialisation(environmentVariables *conditions, gsl_rng *rndarray[])
{
	particleVariables* particles = NULL;

	if (conditions->numberOfParticles == 0)
	{
	    //
	    // Call function to read in particle data
	    //
	    if( (conditions->numberOfParticles = particle_read_in(&particles)) <= 0)
	    {
	        return NULL;
	    }

	    printf("Data read in success\n" );
		fflush(stdout);
	}
	else
	{
		if(generate_particle_data(conditions->numberOfParticles, &particles, rndarray[0], conditions->radius, conditions->xMax, conditions->yMax, conditions->zMax) != 0)
		{
			free_memory(1,particles);
			particles = NULL ;
			return NULL;
		}
		printf("Particle Generation Success\n" );
		fflush(stdout);
	}

    double *generalisedCoordinates = generalised_coordinate_creation( conditions->numberOfParticles, particles);

	free_memory(1,particles);
	particles = NULL ;

	return generalisedCoordinates;

}

void free_memory(int listSize, ...)
{
	va_list args;
	va_start(args, listSize);
	for(int i = 0; i<listSize; i++)
	{
		void *ptr = va_arg(args,void*);

		if(ptr!=NULL)
		free(ptr);
	}
	va_end(args);
}
