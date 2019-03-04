# Modelling Active Particles

A program for modeling active nanoparticle spheres to simulation or induce swarm like behaviour. This is of interest to the scientific community as it has many potential applications in medical and other fields.

##Dependencies
* MPI (Message parsing interface) - [OpenMPI](https://www.open-mpi.org/), [MS-MPI](https://msdn.microsoft.com/en-us/library/bb524831(v=vs.85).aspx)
* GSL (GNU Scientific Library) - [GSL](https://www.gnu.org/software/gsl/)
* OpenMP - Comes with GCC compilers 

##Input/Output
File set up

File must be called particleInput.txt

----------------------------------

numberOfParticles \n

x y z alpha beta gamma\n

... \n

... \n

---------------------------------------------------

This is the order the input values are read.
The first number is the number of expected particles to allocate memory. If
negative or zero it will fail with error 0
Following rows must have 6 values to be read in per row in the order given above
with a newline (\n) at the end.  


Output file outputs as:

time x1 y1 z1 x2 y2 z2 ... xn yn zn alpha1 beta1 gamma1 ... alphan betan gamman


### Input 

If not generating particles using -num [VALUE], an input file is expected to be provided. The file must be as follows:

* Name - particleInput.txt
* Format :

 NumberOfParticles\n

 x0  y0  z0 alpha0 beta0 gamma0\n   

 x1  y1  z1 alpha1 beta1 gamma1\n

 Alpha is the angle in the x-y plane. Beta is the angle in the z-y plane. Gamma is the rotation of the particle along its N-S axis 



### Output

The output comes as 5 files. Each file has a standard name component followed by a number and is in csv format. The number is default to be 0000000 and will overwrite every time the files are written to. TO prevent overwrite you can change the number of the file by using the flag -filenum [VALUE] where VALUE will be appended to the file name padded with zeros. For example:

$ ./program -filenum 1 , this will return a file filename0000001.csv
$ ./program -filenum 5000 , this will return a file filename0005000.csv

Note the maximum value of VALUE is 9999999. Larger values will cause undefined behaviour.

The 5 files outputted are:
* output0000000.csv, this contains the coordinate positions of each particle
* forces_output0000000.csv this contains the coordinate forces acting on each particle
* angle_output0000000.csv this contains the coordinate angles of each particle
* rms_output0000000.csv this contains the root mean square of the displacements
* details0000000.csv this contains the details of the simulation


The files, other than details and rms_output, come in the format:

time0, x0 , y0, z0 , x1, y1, z1, ...

time1, x0 , y0, z0 , x1, y1, z1, ...

Where you replace x y and z with Fx Fy Fz for forces and alpha beta gamma for angles. 

Alpha is the angle in the x-y plane. Beta is the angle in the z-y plane. Gamma is the rotation of the particle along its N-S axis 


#Run Flags
Run flags are used to control the flow of a program by giving a set of conditions the program must follow. Note the space is required.

They usually appear as -{FLAG} [VALUE], for example:

* If you wish to give the maximum number of threads as 2, use, -threads 2

These flags must come after the program run command, for example:

* $ ./program -threads 2  In terminal this tells a program to run with a maximum of 2 threads.
* program.exe -threads 2  In CMD prompt this tells a program to run with a maximum of 2 threads.

Multiple flags can be used on a program for more complex control, for example:

* $ ./program -threads 2 -num 4  In terminal this tells a program to run with a maximum of 2 threadsand 4 particles.
* program.exe -threads 2 -num 4  In CMD prompt this tells a program to run with a maximum of 2 threads and 4 particles.

If using MPI, the number of nodes needs to be specified. This leads to two levels of flags:

* $ mpirun -n 2 ./program -threads 2 -num 4  
* mpiexec -n 2 program.exe -threads 2 -num 4  
Note how the order is important here. The MPI flags must come BEFORE the program name, then the program flags come AFTER the program name. 

### Available Flags
* -debug, tells the program to write all values and variable to file for debugging purposes.
* -num [VALUE] or -n [VALUE], the number of particles to be randomly generated and simulated (ignore if you want to use input in particleInput.txt). VALUE must be positive integer
* -threads [VALUE], the maximum number of threads per node
* -filenum [VALUE], number appended to file name. See [Input Output]
* -cube [VALUE], defines a cube of edge length VALUE where particles will be placed if spawned in program
* -x [VALUE 1], -y [VALUE 2] , -z [VALUE 3] , defines the dimensions of a box in which particles if spawned in program
* -t [VALUE], the maximum simulation time the program should run for
* -dt [VALUE], the time step used to increment each loop
* -temp [VALUE], the temperature the system is to run at
* -drivemag [VALUE], the magnitude of the driving force

