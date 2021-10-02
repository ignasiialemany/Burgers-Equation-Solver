This code solves the burgers equation using MPI paralellisation. It was part of an HPC project in Imperial College London. 

The code can be run through the makefile using the following orders:

`Make diffp` - Compiles if needed and runs diffusion with Nx=101, Ny=101, Nt=1000 and 2 Processors

`Make advxp` - Compiles if needed and runs advection in the x axis with Nx=101, Ny=101, Nt=1000 and 2 Processors

`Make advyp` - Compiles if needed and runs advection in the y axis with Nx=101, Ny=101, Nt=1000 and 2 Processors

`Make burgp` - Compiles if needed and runs the full burgers equation with Nx=101, Ny=101, Nt=1000 and 2 Processors

`Make clean` - Cleans all files

The code can also be run withouth the makefile using boost libraries. The following terminal command can be called.

`mpiexec -n X .\myprog --options (where X is number of procs)`

options are:

`ax` (default value is 1) "Requirement : Parameter ax, ax>0"

`ay` (default value is 0.5) "Requirement : Parameter ay, ay>0"

`b` (default value is 1) "Requirement : Parameter b, b>0"

`c` (default value is 0.02) "Requirement : Parameter c, c>0"

`Nx`(default value is 101) "Requirement : Nx, Nx>0"

`Ny` (default value is 101) "Requirement : Ny, Ny>0"

`Nt` (default value is 100) "Requirement : Nt, Nt>0"

`Lx` (default value is 10) "Requirement : Lx, Lx>0"

`Ly` (default value is 10) "Requirement : Ly, Ly>0"

`T` (default value is 1) "Requirement : T, T>0"

`Px` (default value is 1) "Requirement : Px >= 1"

`Py` (default value is 1)"Requirement : Py >= 1"

The command help shows all the available options and requirements.

`--help` 



