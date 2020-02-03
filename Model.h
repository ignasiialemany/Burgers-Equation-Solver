#ifndef MODEL_H
#define MODEL_H
#include <iostream>
#include "Model.h"
#include "BurgersNova.h"
#include <math.h>
#include <fstream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <mpi.h>

using std::cerr;
using std::cout;
using std::endl;


namespace opt = boost::program_options;

class Model {
    
public:
		
        //CONSTRUCTORS AND DESTRUCTOR
        
                            Model(int argc, char *argv[]);
                            //Ask why default constructor is called when I separate
                            //Function print.parameters
                            Model ();
                            ~Model();

        // Getters
                          
                            double GetX0()     const { return x0; }
                            double GetY0()     const { return y0; }
                            double GetLx()     const { return Lx; }
                            double GetLy()     const { return Ly; }
                            int    GetNx()     const { return Nx; }
                            int    GetNy()     const { return Ny; }
                            int    GetNt()     const { return Nt; }
                            double GetDx()     const { return dx; }
                            double GetDy()     const { return dy; }
                            double GetDt()     const { return dt; }
                            double GetAx()     const { return ax; }
                            double GetAy()     const { return ay; }
                            double GetB()      const { return b; }
                            double GetC()      const { return c; }
                            

                            int Getpy()      const { return py; }
                            int Getpx()      const { return px; }
                            
                            int GetCx()      const { return C_x; }
                            int GetCy()      const { return C_y; }
                           
                            int Getprocessorx()      const { return processor_x; }
                            int Getprocessory()      const { return processor_y; }

private:
       
       //Private Functions
                           
                           //Program options (using boost library) 
                           void programoptions(int argc, char *argv[]);
                           
                           //Prints parameters in terminal after checking if they are valid 
                           void printparameters();
                           
                           //Validats parameters and help message
                           void ValidateParameters();
                           
                           //Calculates other parameters in function of program options and
                           //distributes nodes in the processors. 
                           void CalculateParameters();



        //Private Variables
        
                            //Bool to validate parameters
                            bool verbose;
                                   
                            // Numerics
                            double x0;
                            double y0;
                            double Lx;
                            double T;
                            double Ly;
                            int    Nx;
                            int    Ny;
                            double Nt;
                            double dx;
                            double dy;
                            double dt;
                            
                            //MPI NUMERICS
                            
                            int pNx;
                            int pNy;
                            int sizex;
                            int sizey;
                            int C_x;
                            int C_y;
                            
                            // Physics
                            double ax;
                            double ay;
                            double b;
                            double c;
                            
                            //MPI
                            int px;
                            int py;
                            int processor_x;
                            int processor_y;
                            
                             //Rank and size
                            int rank; 
                            int size;

       
};

#endif // MODEL_H
