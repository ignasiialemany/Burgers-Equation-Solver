#include <iostream>
#include <chrono>
#include "Model.h"
#include "BurgersNova.h"
#include <math.h>
#include <cstdlib>
#include <iomanip>
#include <boost/program_options.hpp>
#include <mpi.h>
#define UNDERLINE "\033[4m"
#define CLOSEUNDERLINE "\033[0m"
#include <fstream>

using namespace std;
namespace opt = boost::program_options;

int main (int argc, char*argv[])
{  
    //Init MPI
    MPI_Init(&argc,&argv);
    
    //Call Model and Burgers class. 
    Model h(argc,argv);
    BurgersNova b(h);
    
    //Call functions from object BurgersNova class (b).

    b.setvalue();
    
    //Gather initial velocity field and compute energy. 
    b.gatherparameters();
    if(b.rank==0)
    {
    cout << "Initial Condition written in file Initial.txt" << endl;
    ofstream file("Initial.txt");
    file << b;
    file.close();    
    cout << UNDERLINE << "Initial Energy"<< CLOSEUNDERLINE << endl;
    }
    
    b.energy();
    
    //Initialize chrono
    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    
    //Iterative process. 
    hrc::time_point start = hrc::now();
    b.iterativeprocess();
    hrc::time_point end = hrc::now();
    
    //Gather and energy after iterative process.
    b.gatherparameters();
    
    if(b.rank==0)
    {
    cout << endl;
    cout << "Final velocity written in file V_final.txt" << endl;    
    ofstream file("V_final.txt");
    file << b;
    file.close();    
    cout << UNDERLINE << "Final Energy"<< CLOSEUNDERLINE << endl;
    }
    
    b.energy();    
    
    MPI_Barrier(b.com2d);
    
    //Cout time. 
    std::cout << "rank " <<  b.rank <<  "Iteration time is "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << "milliseconds" << endl;
    MPI_Finalize();
}


   
  
    
    
    
