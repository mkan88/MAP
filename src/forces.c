/*********************
* Date of creation 19/02/2018
* Authors: Oliver Hinds, Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
**************************************
* History
*
**************************************/

#include <stdlib.h>
#include <omp.h>
#include <math.h>

#include "forces.h"
#include "initial_finalisation.h"

static void force_gravity(double *additionalForces, int numberOfCells, double mass);

static void force_van_der_waals(double *additionalForces, double *generalisedCoordinates, int numberOfCells, double radius);

static void force_exp_repulsion(double *additionalForces, double *generalisedCoordinates, int numberOfCells);

static void alignment_torque(double *additionalForces, double *generalisedCoordinates, int numberOfCells, environmentVariables conditions);

static void viseck_alignment_torque(double *additionalForces, double *generalisedCoordinates, environmentVariables conditions, double rCutoff);

static void driving_force(double *additionalForces, double *generalisedCoordinates, int numberOfCells, field_t drivingField);

static void polar_driving_force(double *additionalForces, double *generalisedCoordinates, int numberOfCells, double drivingForceMagnitude);

extern int gDebug;

extern double gGrav;
extern double gPi;
//
// Tuning varibles for the exponential
// F = Aexp(-r/lamda)
// r= r2 -r1
//
static double gHamaker =2*3.7E-19;
static double gExpTuningA = 7.4E4;
static double gExpTuningB = 4.8E8;
static double gTuningC = -550E-9;
static double gTuningD = 7;


//
// Calculate the sum of the forces and torques from the forces and torques wanted
//
void force_torque_summation(double *additionalForces,double *generalisedCoordinates, int numberOfCells, int *forceList, int numberOfForces,	environmentVariables conditions, field_t drivingField)
{

    double rCutoff = 5E-6;
    //
    // Set entire array to zero
    //
    //#pragma omp parallel for
    for(int i = 0; i < numberOfCells; i++)
    {
        additionalForces[i] = 0;
    }


    for(int i = 0; i < numberOfForces; i++)
    {
        switch(forceList[i])
        {
            case NONE : break;
            case GRAVITY : force_gravity(additionalForces, numberOfCells, conditions.mass); break;
            case VAN_DER_WAALS : force_van_der_waals(additionalForces, generalisedCoordinates, numberOfCells, conditions.radius); break;
            case EXP_REPULSION : force_exp_repulsion(additionalForces, generalisedCoordinates, numberOfCells); break;
			case ALIGN_TORQUE : alignment_torque(additionalForces, generalisedCoordinates, numberOfCells, conditions); break;
			case DRIVING_FIELD : driving_force(additionalForces, generalisedCoordinates, numberOfCells, drivingField); break;
            case POLAR_DRIVING_FORCE : polar_driving_force(additionalForces, generalisedCoordinates, numberOfCells, conditions.drivingForceMagnitude); break;
            case VISECK_ALIGN_TORQUE : viseck_alignment_torque(additionalForces, generalisedCoordinates,  conditions, rCutoff);break;
            default : break;
        }
    }
}

static void force_gravity(double *additionalForces, int numberOfCells, double mass)
{
    //
    // Half the number of cells to ignore torques.
    // Step 3 each time to only hit the z axis
    // Add the forces onto the preexisting values
    //
    #pragma omp parallel for
    for(int i = 0; i < numberOfCells/2; i+=3)
    	additionalForces[i] -= mass * gGrav;
}

static void force_van_der_waals(double *additionalForces, double *generalisedCoordinates, int numberOfCells, double radius)
{
    //
    // get the number of particles.
    //
    int numberOfParticles = numberOfCells/6;
    // F = -32ArR^6/|r|^4(|r|^2 - 4R^2)^2
    // r = r2-r1
    // Add the forces onto the preexisting values
    //


    #pragma omp parallel for
    for(int i = 0; i < numberOfParticles; i++)
    {
        for(int j = 0; j < numberOfParticles; j++)
        {
            double x,y,z,r,forceConst;
            //
            // Ignore self interaction
            //
            if(i == j)
            continue ;
            //
            // Calculate the vector components of r
            //
            x = generalisedCoordinates[j * 3] - generalisedCoordinates[i * 3];
            y = generalisedCoordinates[j * 3 + 1] - generalisedCoordinates[i * 3 + 1];
            z = generalisedCoordinates[j * 3 + 2] - generalisedCoordinates[i * 3 + 2];
            //
            // Calculate |r|
            //
            r = sqrt(x * x + y * y + z * z);
            //
            // Calculate  -32AR^6 / 3|r|^3(|r|^2 - 4R^2)^2
            //
            forceConst = 32*gHamaker*pow(radius,6)*pow(gTuningD,6) / ( 3 * pow(gTuningC-r,3)*pow( ( pow(gTuningC -r,2) - 4*pow(radius*gTuningD,2)) , 2) );
            //
            // Vectorise
            //
            additionalForces[j * 3] += forceConst * x/r;
            additionalForces[j * 3 + 1] += forceConst * y/r;
            additionalForces[j * 3 + 2] += forceConst * z/r;
        }
    }
}

static void force_exp_repulsion(double *additionalForces, double *generalisedCoordinates, int numberOfCells)
{
    //
    // get the number of particles.
    //
    int numberOfParticles = numberOfCells/6;
    // F = Aexp(-r/lamda)
    // r = r2-r1
    // Add the forces onto the preexisting values
    //

    #pragma omp parallel for
    for(int i = 0; i < numberOfParticles; i++)
    {
        for(int j = 0; j < numberOfParticles; j++)
        {
            double x,y,z,r,forceConst;
            //
            // Ignore self interaction
            //
            if(i == j)
            continue ;
            //
            // Calculate the vector components of r
            //
            x = generalisedCoordinates[j * 3] - generalisedCoordinates[i * 3];
            y = generalisedCoordinates[j * 3 + 1] - generalisedCoordinates[i * 3 + 1];
            z = generalisedCoordinates[j * 3 + 2] - generalisedCoordinates[i * 3 + 2];
            //
            // Calculate |r|
            //
            r = sqrt(x * x + y * y + z * z);
            //
            // Calculate 1/R * Aexp(-|r|/lamda)
            //
            forceConst =  gExpTuningA * gExpTuningB * ( 1/gTuningD)*exp(gExpTuningB * (gTuningC - r)/gTuningD);
            //
            // Vectorise
            //
            additionalForces[j * 3] += forceConst * x/r;
            additionalForces[j * 3 + 1] += forceConst * y/r;
            additionalForces[j * 3 + 2] += forceConst * z/r;
        }
    }

}

static void alignment_torque(double *additionalForces, double *generalisedCoordinates, int numberOfCells, environmentVariables conditions)
{
	int numberOfParticles = numberOfCells/6;
	int rotOffset = numberOfCells/2;

	// Calculated from rotational version of Langevin equation, substituting that the average angular displacement per timestep is pi/2
	// Torque T= pi*I/dt^2 for change in angle ~pi/2 per timestep
	// Moment of intertia I = (2/5)MR^2 for solid sphere of radius R and mass M
	double forceConst = gPi * 0.4 * conditions.mass * pow(conditions.radius, 2) * conditions.radius / pow(conditions.deltaTime, 2);
	/*															  |						  |
															Kept separate simply for clarity, might
														  be more comptuationally efficient to combine.
		 The second radius multiplier comes from using the inverse distance mutliplier in the force equation
		 	- the result is that the force is normalised based on particles separated by a distance of the same order as their radius
	*/
	double totalX = 0;
	double totalY = 0;
	double totalZ = 0;
	double totalAlpha = 0;
	double totalBeta = 0;

	double meanX, meanY, meanZ;
	double meanAlpha, meanBeta, difAlpha, difBeta;




    #pragma omp parallel
    {
	// Sum position and angles
        #pragma omp for reduction(+:totalX, totalY,totalZ,totalAlpha,totalBeta)
    	for (int i=0; i<numberOfParticles; i++)
    	{
    		totalX += generalisedCoordinates[3*i + 0];
    		totalY += generalisedCoordinates[3*i + 1];
    		totalZ += generalisedCoordinates[3*i + 2];


    		// Summing only alpha and beta angles
    		totalAlpha += generalisedCoordinates[rotOffset + 3*i + 0];
    		totalBeta += generalisedCoordinates[rotOffset + 3*i + 1];
    	}

    	// Calculate average positions
        #pragma omp single
        {
        	meanX = totalX/numberOfParticles;
        	meanY = totalY/numberOfParticles;
        	meanZ = totalZ/numberOfParticles;

        	// Calculate average angles

        	meanAlpha = totalAlpha/numberOfParticles;
        	meanBeta = totalBeta/numberOfParticles;
        }
    	// Calculate torques on each particle according to their alignment with average angle
        #pragma omp for
    	for (int i=0; i<numberOfParticles; i++)
    	{
            double dist,distMul;
    		// Calculate the distance between the particle and the average position
    		dist = sqrt(pow(meanX - generalisedCoordinates[3*i + 0],2) + pow(meanY - generalisedCoordinates[3*i + 1],2) + pow(meanZ - generalisedCoordinates[3*i + 2],2));
    		distMul = 1/(dist + conditions.radius); // Distance multiplier
    		// Calculate torques in alpha and beta directions such that maximum torque is when the particle's axes are maximally separated from the average angle.
    		additionalForces[rotOffset + 3*i + 0] += difAlpha = distMul * forceConst * sin(meanAlpha - generalisedCoordinates[rotOffset + 3*i + 0]);
    		additionalForces[rotOffset + 3*i + 1] += difBeta = distMul * forceConst * sin(meanBeta - generalisedCoordinates[rotOffset + 3*i + 1]);
    		/* Applies inverse distance multiplier as previously mentioned, which normalises the force so that it should be producing the average ~pi/2 angular change per timestep when the
    			particle is in the order of a radius of the average position*/


    		// Print stuff for debugging
    		if (i<3 && gDebug==1) // Avoid some spam
    		{
    			printf("\ts:%e\tc:%e\tTheta:%e\tPhi:%e\n", dist, forceConst, generalisedCoordinates[rotOffset + 3*i + 0], generalisedCoordinates[rotOffset + 3*i + 1]);
    			printf("Theta:\t%e\t%e\t%e\t%e\t%e\n", fmod(generalisedCoordinates[rotOffset + 3*i + 0], 2*gPi), distMul, forceConst, sin(meanAlpha - generalisedCoordinates[rotOffset + 3*i + 0]), additionalForces[rotOffset + 3*i + 0]);
    			printf("Phi:\t%e\t%e\t%e\t%e\t%e\n", fmod(generalisedCoordinates[rotOffset + 3*i + 1], 2*gPi), distMul, forceConst, sin(meanBeta - generalisedCoordinates[rotOffset + 3*i + 1]), additionalForces[rotOffset + 3*i + 1]);
    		}
    	}
    }

}

static void driving_force(double *additionalForces, double *generalisedCoordinates, int numberOfCells, field_t drivingField)
{
    #pragma omp parallel for
	for (int i=0; i<numberOfCells/6; i++)
	{
		additionalForces[3*i + 0] +=  drivingField.mag * abs(cos((drivingField.alpha - generalisedCoordinates[numberOfCells/2 + 3*i + 0])/2)) * abs(cos((drivingField.beta - generalisedCoordinates[numberOfCells/2 + 3*i + 1])/2));
		//additionalForces[3*i + 1] +=  drivingField.mag * cos(drivingField.beta - generalisedCoordinates[numberOfCells/2 + 3*i + 1]);
	}

}
static void polar_driving_force(double *additionalForces, double *generalisedCoordinates, int numberOfCells, double drivingForceMagnitude)
{
    #pragma omp parallel for
	for (int i=0; i<numberOfCells/6; i++)
	{
		additionalForces[3*i + 0] +=  drivingForceMagnitude * cos(generalisedCoordinates[numberOfCells/2 + 3*i + 0]) * sin (generalisedCoordinates[numberOfCells/2 + 3*i + 1]) ;
		additionalForces[3*i + 1] +=  drivingForceMagnitude * sin (generalisedCoordinates[numberOfCells/2 + 3*i + 0]) * sin(generalisedCoordinates[numberOfCells/2 + 3*i + 1]) ;
        additionalForces[3*i + 2] +=  drivingForceMagnitude * cos(generalisedCoordinates[numberOfCells/2 + 3*i + 1]) ;
	}

}

static void viseck_alignment_torque(double *additionalForces, double *generalisedCoordinates, environmentVariables conditions, double rCutoff)
{

    //***ALGORITHM***//
    // - Pick a particle
    // - Calculate average angle within r
    // - adjust angle of particle
    // - move to next particle

	int numberOfCells = conditions.numberOfParticles/6;
	int rotOffset = numberOfCells/2;

	// Calculated from rotational version of Langevin equation, substituting that the average angular displacement per timestep is pi/2
	// Torque T= pi*I/dt^2 for change in angle ~pi/2 per timestep
	// Moment of intertia I = (2/5)MR^2 for solid sphere of radius R and mass M
	double forceConst = gPi * 0.4 * conditions.mass * pow(conditions.radius, 2) * conditions.radius / pow(conditions.deltaTime, 2);
	/*															  |						  |
															Kept separate simply for clarity, might
														  be more comptuationally efficient to combine.
		 The second radius multiplier comes from using the inverse distance mutliplier in the force equation
		 	- the result is that the force is normalised based on particles separated by a distance of the same order as their radius
	*/

    #pragma omp parallel for
    for(int i = 0; i < conditions.numberOfParticles; i++)
    {
        int meanAlpha, meanBeta, difAlpha,difBeta,totalAlpha=0,totalBeta=0;
        int numNeighbours = 0;
        for(int j = 0; j < conditions.numberOfParticles; j++)
        {
            double x,y,z,r2;
            //
            // Calculate the vector components of r
            //
            x = generalisedCoordinates[j * 3] - generalisedCoordinates[i * 3];
            y = generalisedCoordinates[j * 3 + 1] - generalisedCoordinates[i * 3 + 1];
            z = generalisedCoordinates[j * 3 + 2] - generalisedCoordinates[i * 3 + 2];
            //
            // Calculate |r|^2
            //
            r2 = x * x + y * y + z * z;

            if(sqrt(r2) <= rCutoff)
            {
            	// Summing only alpha and beta angles
            	totalAlpha += generalisedCoordinates[rotOffset + 3*j + 0];
            	totalBeta += generalisedCoordinates[rotOffset + 3*j + 1];
                numNeighbours++;
            }

        }

        meanAlpha = totalAlpha/numNeighbours;
        meanBeta = totalBeta/numNeighbours;
        // Calculate torques in alpha and beta directions such that maximum torque is when the particle's axes are maximally separated from the average angle.

        additionalForces[rotOffset + 3*i + 0] += difAlpha = forceConst * sin(meanAlpha - generalisedCoordinates[rotOffset + 3*i + 0]);
        additionalForces[rotOffset + 3*i + 1] += difBeta = forceConst * sin(meanBeta - generalisedCoordinates[rotOffset + 3*i + 1]);
    }

}
