CXX=mpicxx 
LDLIBS=-lboost_system -lboost_program_options

.PHONY: clean all default

all: compile

default: compile

compile: main.o model.o burgers.o 
	$(CXX) -o myprog main.o model.o burgers.o $(LDLIBS)
main.o: main.cpp 
	$(CXX) -o main.o -c main.cpp
model.o: Model.cpp
	$(CXX) -o model.o -c Model.cpp
burgers.o: BurgersNova.cpp
	$(CXX) -o burgers.o -c BurgersNova.cpp

diffp: compile
	mpiexec -n 2 ./myprog --ax=0 --ay=0  --b=0 --c=1 --Nx=101 --Ny=101 --Nt=1000 --Px=1 --Py=2

advxp: compile
	mpiexec -n 2 ./myprog --ax=1 --ay=0  --b=0 --c=0 --Nx 101 --Ny=101 --Nt=1000 --Px=1 --Py=2

advyp: compile
	mpiexec -n 2 ./myprog --ax=0 --ay=1  --b=0 --c=0 --Nx=101 --Ny=101 --Nt=1000 --Px=1 --Py=2

burgp: compile 
	mpiexec -n 2 ./myprog --ax=1 --ay=0.5  --b=1 --c=0.02 --Nx=101 --Ny=101 --Nt=1000 --Px=1 --Py=2

diff: compile
	mpiexec -n 1 ./myprog --ax=0 --ay=0  --b=0 --c=1 --Nx=101 --Ny=101 --Nt=1000 --Px=1 --Py=1

advx: compile
	mpiexec -n 1 ./myprog --ax=1 --ay=0  --b=0 --c=0 --Nx 101 --Ny=101 --Nt=1000 --Px=1 --Py=1

advy: compile
	mpiexec -n 1 ./myprog --ax=0 --ay=1  --b=0 --c=0 --Nx=101 --Ny=101 --Nt=1000 --Px=1 --Py=1

burg: compile 
	mpiexec -n 1 ./myprog --ax=1 --ay=0.5  --b=1 --c=0.02 --Nx=101 --Ny=101 --Nt=1000 --Px=1 --Py=1

clean:
	rm *.o myprog


