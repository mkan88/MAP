/*********************
* Date of creation 09/10/2017
* Authors: Oliver Hinds, Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
**************************************
* Change History
**************************************/
//
// GNU and external libraries
//
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <omp.h>
#include "mpi.h"
//
// Bespoke libraries
//
#include "particles.h"
#include "diffusionmatrix.h"
#include "stochastic_force.h"
#include "moving_on.h"
#include "forces.h"
#include "initial_finalisation.h"

#define MASTER 0
#define TOP_SLAVE 1


double gBoltzmannConst = 1.38064852E-23; // m^2 kg s^-2 K^-1
double gPi = 3.14159265359;
double gGrav = 9.80665; // m s^-2

int gDebug = 0;
int gSerial = 0;
int gNumOfthreads;
int gNumOfNodes = 1;

enum MESSAGE_TAGS
{
    NUM_OF_PARTICLES,
    NUM_OF_FORCES ,         // Number of forces to expect
    FORCE_LIST,             // Which forces to use
    COORDINATES ,           // Corordinates to use
    ADDITIONAL_FORCES       // Total forces
};

int main(int argc, char *argv[])
{
    //
    // MPI Initilisation
    //
    int taskid;
    int MPI_error = 0;
    MPI_Status status;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &gNumOfNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);


	//
	// Create conditions variable
	//
	environmentVariables conditions;
	//
	// Initilise conditions with boilerplate variables
	//
	boilerplate_variables(&conditions);
	//
	// Read in cmd line arguments and adjust conditions as neccessary
	//
	cmd_line_read_in(argc, argv, &conditions);

    omp_set_num_threads(gNumOfthreads);
    // Create driving field
    field_t drivingField;
    drivingField.mag = 1E-10;
    drivingField.alpha = 0;
    drivingField.beta = gPi/2;


    //
    //  Choose forces to be included
    //
    int numberOfForces = 2; // must be at least 1, with the force none chosen
    //
    // Copy of enum to understand force forceList
    //
    //enum forces_available
    //{
    //    NONE ,
    //    GRAVITY ,
    //    VAN_DER_WAALS ,
    //    EXP_REPULSION ,
    //	  ALIGN_TORQUE ,
    //	  DRIVING_FIELD,
    //	  POLAR_DRIVING_FORCE,
    //    VISECK_ALIGN_TORQUE
    //};

    int forceList[4] = {VAN_DER_WAALS,EXP_REPULSION,POLAR_DRIVING_FORCE,VISECK_ALIGN_TORQUE};
    //
    // Master process environment
    //
    if( MASTER == taskid )
    {
        //
    	//Open files for output
    	//
        char filename1[sizeof "../results/output1000000.csv"];
        char filename2[sizeof "../results/angle_output1000000.csv"];
        char filename3[sizeof "../results/forces_output1000000.csv"];
        char filename4[sizeof "../results/rms_output1000000.csv"];
        char filename5[sizeof "../results/diffusion_output1000000.csv"];
        char filename6[sizeof "../results/details1000000.csv"];

        sprintf(filename1, "../results/output%07d.csv", conditions.fileNum);
        sprintf(filename2, "../results/angle_output%07d.csv", conditions.fileNum);
        sprintf(filename3, "../results/forces_output%07d.csv", conditions.fileNum);
        sprintf(filename4, "../results/rms_output%07d.csv", conditions.fileNum);
        sprintf(filename5, "../results/diffusion_output%07d.csv", conditions.fileNum);
        sprintf(filename6, "../results/details%07d.csv", conditions.fileNum);


        FILE *output = fopen(filename1,"w");
        FILE *angle_output = fopen(filename2,"w");
    	FILE *forces_output = fopen(filename3,"w");
    	FILE *rms_output = fopen(filename4,"w");
        FILE *details = fopen(filename6,"w");
        if(output == NULL || angle_output == NULL || forces_output == NULL || rms_output == NULL  || details==NULL)
        {
            printf("-Error %d : %s\n : File %s : Line : %d", errno, strerror( errno ), __FILE__, __LINE__);
            MPI_Abort(MPI_COMM_WORLD, MPI_error);
            return -errno;
        }

    	//
    	// Initilise random variables
    	//
    	gsl_rng **rndarray=rand_array_allocation();

    	if(rndarray == NULL)
    	{
            printf("-Problem with rndarray");
            MPI_Abort(MPI_COMM_WORLD, MPI_error);
    		return -1;
    	}

        //
        // Initilise generalised coordinates
        //

    	double *generalisedCoordinates = generalised_coordinate_initialisation(&conditions,rndarray);

    	if(generalisedCoordinates == NULL)
    	{
            printf("Problem with generalisedCoordinates");
            MPI_Abort(MPI_COMM_WORLD, MPI_error);
    		return -1;
    	}

        int vectorSize = 6 * conditions.numberOfParticles;
        //
        // Send number of particles to slaves
        //
        if(gNumOfNodes > 1)
        {
            MPI_Send(&conditions.numberOfParticles, 1, MPI_INT, TOP_SLAVE, NUM_OF_PARTICLES, MPI_COMM_WORLD);
        }

        //---------------------------- DEBUG ------------------------------//
        //
        // Prints the generalisedCoordinates to a file for inspection
        //
        if( gDebug == 1 && generalisedCoordinates != NULL)
        {
            FILE *genCoordOutput = fopen("../bin/genCoord_output.txt","w");

            for(int i = 0; i < vectorSize; i++)
            {
                fprintf(genCoordOutput, "%e\n", generalisedCoordinates[i]);
            }
            fclose (genCoordOutput);
        }
        //--------------------------- END ---------------------------------//


        //
        // Allocate memory required for the program.
        // Requires: Diffusion matrix, stochastic displacement,
        //          additional forces, stochasticWeighting,
        //
        //


        double *diffusionMatrix = NULL ;
        double *stochasticWeighting = NULL;
        double *additionalForces = NULL;
        double *stochasticDisplacement = NULL;
		double *rndNumArray = NULL;
        diffusionMatrix = calloc( pow(vectorSize, 2), sizeof *diffusionMatrix) ;
        stochasticWeighting = calloc( pow( vectorSize, 2), sizeof *stochasticWeighting);
        stochasticDisplacement = calloc( vectorSize, sizeof *stochasticDisplacement);
        additionalForces = calloc( vectorSize, sizeof *additionalForces);
		rndNumArray = calloc( vectorSize, sizeof *rndNumArray);

        if(  diffusionMatrix==NULL  || stochasticWeighting==NULL || stochasticDisplacement==NULL || additionalForces==NULL || rndNumArray==NULL)
        {
    		free_memory(6,diffusionMatrix, generalisedCoordinates, stochasticWeighting, stochasticDisplacement,additionalForces,rndNumArray);
    		diffusionMatrix = generalisedCoordinates = stochasticWeighting = stochasticDisplacement = additionalForces = rndNumArray=NULL ;
            printf("-Error %d : %s\n : File %s : Line : %d", errno, strerror( errno ), __FILE__, __LINE__);
            MPI_Abort(MPI_COMM_WORLD, MPI_error);
            return -errno;
        }

        fprintf(details,"Particle Num %d\n",conditions.numberOfParticles);
        fprintf(details,"Temp %g K\n",conditions.temperature);
        fprintf(details,"Polar driving magnitude %g N\n",conditions.drivingForceMagnitude);
        fprintf(details,"Viscosity %g Nm^-2 s\n",conditions.viscosity);
        fprintf(details,"Radius %g m\n",conditions.radius);
        fflush(details);

        //
        // GPU initialisation
        //
        initialise_gpu();

        //
        // Loop through time, output each time step to a file.
        //
        int loop = 0;
        int count = 0;
        int maxLoop = conditions.endTime/(double)conditions.deltaTime;

        double progTime = omp_get_wtime();

        while(conditions.currentTime<conditions.endTime)
        {
			//
			// Save data
			//

			if(loop%100 == 0)
            {
				/* Count number of particles within the volume in which the particles were initially generated,
					centred on the average position*/

				// Calculate the mean position
				double totalX = 0;
				double totalY = 0;
				double totalZ = 0;
				double meanX,meanY,meanZ;

				for (int i=0; i < 3*conditions.numberOfParticles; i+= 3)
				{
					// Sum each particle's position
					totalX += generalisedCoordinates[i + 0];
					totalY += generalisedCoordinates[i + 1];
					totalZ += generalisedCoordinates[i + 2];
				}

				meanX = totalX/conditions.numberOfParticles;
				meanY = totalY/conditions.numberOfParticles;
				meanZ = totalZ/conditions.numberOfParticles;

				// Calculate the root mean square displacement relative to the average position
				double totalSquare = 0;
				double rootMeanSquare;
				for (int i=0; i < 3*conditions.numberOfParticles; i+= 3)
				{
					/* Sum each particle's total relative displacement (i.e. its position relative to the average position now,
						minus its initial position relative to the initial average position)*/
					totalSquare +=	  pow((generalisedCoordinates[i + 0] - meanX) , 2)
									+ pow((generalisedCoordinates[i + 1] - meanY) , 2)
									+ pow((generalisedCoordinates[i + 2] - meanZ) , 2);

					/* Count the total number of particles within the volume in which the particles were initially generated,
						centred on the average position*/
				}
				rootMeanSquare = sqrt(totalSquare/conditions.numberOfParticles);



    			int angle_offset = 3*conditions.numberOfParticles;
                fprintf(output, "%e, ", conditions.currentTime);
                fprintf(angle_output, "%e, ", conditions.currentTime);
    			fprintf(forces_output, "%e, ",conditions.currentTime);
    			fprintf(rms_output, "%e, %e, %e\n",conditions.currentTime, rootMeanSquare, sqrt(2 * diffusionMatrix[0] * conditions.currentTime));
                fflush(output);
                fflush(angle_output);
                fflush(forces_output);
                fflush(rms_output);
                for(int i = 0; i < 3 * conditions.numberOfParticles; i++)
                {
                    fprintf(output, "%e", generalisedCoordinates[i]);
    				fprintf(forces_output, "%e", additionalForces[i]);
                    double temp =fmod(generalisedCoordinates[angle_offset + i],2*gPi);
                    if(temp < 0)
                    {
                        temp = 2*gPi + temp;
                    }
                    fprintf(angle_output, "%e", temp);
                    fflush(output);
                    fflush(angle_output);
                    fflush(forces_output);
    				if (i < 3*conditions.numberOfParticles - 1)
                    {
    	                fprintf(output, ", ");
    	                fprintf(angle_output, ", ");
    					fprintf(forces_output, ", ");
                        fflush(output);
                        fflush(angle_output);
                        fflush(forces_output);
                    }
                }
                fprintf(output, "\n");
                fprintf(angle_output, "\n");
    			fprintf(forces_output, "\n");
                fflush(output);
                fflush(angle_output);
                fflush(forces_output);
            }
            conditions.currentTime+=conditions.deltaTime; // time step
    		if((maxLoop/10)*count == loop)
    		{
    			//printf("%d%%\n",count*10);
    			count++;
    		}
            loop ++;

            //
            // Create diffusion matrix
            //
            if(gNumOfNodes > 1)
            {
                MPI_Send(&generalisedCoordinates[0], vectorSize, MPI_DOUBLE, TOP_SLAVE, COORDINATES, MPI_COMM_WORLD);
            }

    	    diffusion_matrix_creation( conditions.numberOfParticles, diffusionMatrix, stochasticWeighting, generalisedCoordinates, &conditions);

    	    //---------------------------- DEBUG------------------------------//
    	    //
    	    // Prints the diffusionMatrix to a file for inspection
    	    //
    	    if( gDebug == 1 && diffusionMatrix != NULL)
    	    {
    	        //conditions.currentTime = conditions.endTime+1;
    	        FILE *matrixOutput = fopen("../bin/matrix_output.txt","w");

    	        for(int i = 0; i < 6 * conditions.numberOfParticles; i++)
    	        {
    	            for(int j = 0; j < 6 * conditions.numberOfParticles; j++)
    	            {
    	                fprintf(matrixOutput, "%e\t", diffusionMatrix[i * 6 * conditions.numberOfParticles + j]);
    	            }
    	            fprintf(matrixOutput, "\n");

    	        }
    	        fclose (matrixOutput);
    	    }
    	    //---------------------------END---------------------------------//


    	    //
    	    // Create the stochastic displacement
    	    //

    	    stochastic_displacement_creation( conditions.numberOfParticles, stochasticWeighting, stochasticDisplacement, rndarray,rndNumArray, conditions.deltaTime);

    		if( gDebug == 1 && stochasticWeighting != NULL)
    	    {
    	    	//conditions.currentTime = conditions.endTime+1;
    	        FILE *stochasticOutput = fopen("../bin/stochastic_matrix_output.txt","w");

    	        for(int i = 0; i < 6 * conditions.numberOfParticles; i++)
    	        {
    	            for(int j = 0; j < 6 * conditions.numberOfParticles; j++)
    	            {
    	                fprintf(stochasticOutput, "%e\t", stochasticWeighting[i * 6 * conditions.numberOfParticles + j]);
    	            }
    	            fprintf(stochasticOutput, "\n");
    	        }
    	        fclose (stochasticOutput);
    		}


    		//
    		// Include additional forces (only calculates if its the only node)
    		//
            if(gNumOfNodes <= 1)
            {
                force_torque_summation(additionalForces, generalisedCoordinates, 6 * conditions.numberOfParticles, forceList, numberOfForces, conditions, drivingField);
            }
            else
            {
                MPI_Recv(&additionalForces[0], vectorSize, MPI_DOUBLE, TOP_SLAVE, ADDITIONAL_FORCES, MPI_COMM_WORLD, &status);
            }

            //
            // Calculate time step.
            //
            moving_on_routine(conditions.numberOfParticles, &conditions, diffusionMatrix, additionalForces, stochasticDisplacement, generalisedCoordinates, NULL);
        }
    	progTime = omp_get_wtime() - progTime;

    	fprintf(details,"Run time %lf s\n",progTime);

        //
        // Free memory
        //

        fclose(output);
    	fclose(angle_output);
    	fclose(forces_output);
        fclose(details);
        fclose(rms_output);

    	free_memory(5,diffusionMatrix, generalisedCoordinates, stochasticWeighting, stochasticDisplacement, additionalForces);
    	diffusionMatrix = generalisedCoordinates = stochasticWeighting = stochasticDisplacement = additionalForces = NULL ;

        finalise_gpu();
    }

    if( TOP_SLAVE == taskid )
    {
        //
        // Recieve number of particles and update conditions
        //
        MPI_Recv(&conditions.numberOfParticles, 1, MPI_INT, MASTER, NUM_OF_PARTICLES, MPI_COMM_WORLD, &status);

        int vectorSize = 6 * conditions.numberOfParticles;

        //
        // Allocate memory required for the program.
        // Requires: additional forces, generalised coordinates
        //
        double *additionalForces = NULL;
        double *generalisedCoordinates = NULL;

        generalisedCoordinates = calloc( vectorSize, sizeof *generalisedCoordinates);
        additionalForces = calloc( vectorSize, sizeof *additionalForces);

        if( generalisedCoordinates==NULL || additionalForces==NULL)
        {
            free_memory(2, generalisedCoordinates,additionalForces);
            generalisedCoordinates = additionalForces = NULL ;
            printf("-Error %d : %s\n : File %s : Line : %d", errno, strerror( errno ), __FILE__, __LINE__);
            MPI_Abort(MPI_COMM_WORLD, MPI_error);
            return -errno;
        }

        while(conditions.currentTime<=conditions.endTime)
        {
            //
            // Recieve coordinates
            //
            MPI_Recv(&generalisedCoordinates[0], vectorSize, MPI_DOUBLE, MASTER, COORDINATES, MPI_COMM_WORLD, &status);
    		//
    		// Calculate additional forces
    		//
            force_torque_summation(additionalForces, generalisedCoordinates, 6 * conditions.numberOfParticles, forceList, numberOfForces, conditions, drivingField);
            //
            // Send new forces
            //
            MPI_Send(&additionalForces[0], vectorSize, MPI_DOUBLE, MASTER, ADDITIONAL_FORCES, MPI_COMM_WORLD);

            conditions.currentTime+=conditions.deltaTime; // time step
    	}
        free_memory(2, generalisedCoordinates,additionalForces);
        generalisedCoordinates = additionalForces = NULL ;

    }
    MPI_Finalize();
    return 0;
}
