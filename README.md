MAKEFILE:

Make --> compiles the files needed to run the code
Make diffp --> (and compile if needed) runs diffusion with Nx=101, Ny=101, Nt=1000 and 2 Processors
Make advxp --> (and compile if needed) runs advx with Nx=101, Ny=101, Nt=1000 and 2 Processors
Make advyp --> (and compile if needed) runs advy with Nx=101, Ny=101, Nt=1000 and 2 Processors
Make burgp --> (and compile if needed) runs burg with Nx=101, Ny=101, Nt=1000 and 2 Processors
Make clean --> cleans files

RUN THE CODE WITHOUT MAKEFILE (WITH BOOST LIBRARY):

mpiexec -n X .\myprog --options (where X is number of procs)

--options are:
--ax (default value is 1) "Requirement : Parameter ax, ax>0"
--ay (default value is 0.5) "Requirement : Parameter ay, ay>0"
--b (default value is 1) "Requirement : Parameter b, b>0"
--c (default value is 0.02) "Requirement : Parameter c, c>0"
--Nx (default value is 101) "Requirement : Nx, Nx>0"
--Ny (default value is 101) "Requirement : Ny, Ny>0"
--Nt (default value is 100) "Requirement : Nt, Nt>0"
--Lx (default value is 10) "Requirement : Lx, Lx>0"
--Ly (default value is 10) "Requirement : Ly, Ly>0"
--T (default value is 1) "Requirement : T, T>0"
--Px (default value is 1) "Requirement : Px >= 1"
--Py (default value is 1)"Requirement : Py >= 1"
--help (help message, shows all the available options and requirements). 
       


