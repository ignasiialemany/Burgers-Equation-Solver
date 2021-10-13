This code solves the burgers equation using MPI paralellisation. It was part of an HPC project in Imperial College London. 

The code can be run through the makefile using the following orders:

`Make diffp` - Compiles if needed and runs diffusion with Nx=101, Ny=101, Nt=1000 and 2 Processors

`Make advxp` - Compiles if needed and runs advection in the x axis with Nx=101, Ny=101, Nt=1000 and 2 Processors

`Make advyp` - Compiles if needed and runs advection in the y axis with Nx=101, Ny=101, Nt=1000 and 2 Processors

`Make burgp` - Compiles if needed and runs the full burgers equation with Nx=101, Ny=101, Nt=1000 and 2 Processors

`Make clean` - Cleans all files

The code can also be run withouth the makefile using boost libraries. The following terminal command can be called.

`mpiexec -n X .\myprog --options (where X is number of procs)`

The command help shows all the available options and requirements.

`--help` 



