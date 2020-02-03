#include <iostream>
#include "Model.h"
#include "BurgersNova.h"
#include <math.h>
#include <fstream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <mpi.h>
#define UNDERLINE "\033[4m"
#define CLOSEUNDERLINE "\033[0m"


using std::cerr;
using std::cout;
using std::endl;

namespace opt = boost::program_options;

//Default constructor
Model::Model ()
{
    
}
//Defined constructor
Model::Model(int argc, char *argv[])
{
    //Declare rank for each processor and the num of processors.
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    /*
       PROGRAM OPTIONS(Parameters)
        * ax    Nx      Px  
        * ay    Ny      Py
        * b     CFL     help
        * c     
    */
    //ax,ay,c,b are required
    //Px,Py =3 and Nx,Ny=11 (for default)
    
    //Program options. Defines,validates and calculates other parameters
    programoptions(argc,argv); 
    
    MPI_Barrier(MPI_COMM_WORLD);
    
}
    
Model::~Model()
{
}
void Model::programoptions(int argc, char *argv[])
{
    //PROGRAM OPTIONS;
        opt::options_description desc ("All Options");
        desc.add_options()
        //Parameters ax,ay,b,c
        ("ax",opt::value<double>()->default_value(1),"Requirement : Parameter ax, ax>0")
        ("ay",opt::value<double>()->default_value(0.5),"Requirement : Parameter ay, ay>0")
        ("b",opt::value<double>()->default_value(1),"Requirement : Parameter b, b>0")
        ("c",opt::value<double>()->default_value(0.02),"Requirement : Parameter c, c>0")
        //Nx,Ny,CFL
        ("Nx",opt::value<double>()->default_value(101),"Requirement : Nx, Nx>0")
        ("Ny",opt::value<double>()->default_value(101),"Requirement : Ny, Ny>0")
        ("Nt",opt::value<double>()->default_value(1000),"Requirement : Nt, Nt>0")
        ("Lx",opt::value<double>()->default_value(10),"Requirement : Lx, Lx>0")
        ("Ly",opt::value<double>()->default_value(10),"Requirement : Ly, Ly>0")
        ("T",opt::value<double>()->default_value(1),"Requirement : T, T>0")
        //Px,Py (subdomains in the x and y directions)
        ("Px",opt::value<int>()->default_value(1),"Requirement : Px >= 1")
        ("Py",opt::value<int>()->default_value(1),"Requirement : Py >= 1")
        ("help","Help message");
        
        
        //Saves variables into vm and stores it
        opt::variables_map vm;
        opt::store(opt::parse_command_line(argc, argv, desc),vm);
        opt::notify(vm);
        
        //All processors store the data inside the class.
         
        ay = vm["ay"].as<double>();
        ax = vm["ax"].as<double>();
        b  = vm["b"].as<double>();
        c  = vm["c"].as<double>();
        Nx  = vm["Nx"].as<double>();
        Ny  = vm["Ny"].as<double>();
        Nt = vm["Nt"].as<double>();
        px = vm["Px"].as<int>();
        py = vm["Py"].as<int>();
        Lx = vm["Lx"].as<double>();
        Ly = vm["Ly"].as<double>();
        T = vm["T"].as<double>();
        
        //Calculate parameters  
        CalculateParameters();
        
        
        //Processor 0 checks parameters, if any conditions is not satisfied aborts all other 
        // processors. 
      
      MPI_Barrier(MPI_COMM_WORLD);
      
    if(rank==0)
        {
                                ValidateParameters();
            
                                if (vm.count("help"))   {
                                cout << UNDERLINE << "VALUES FOR BURGERS EQUATION" << CLOSEUNDERLINE << endl;
                                cout << desc << endl;
                                cout << "Aborting all processors" << endl;
                                MPI_Abort(MPI_COMM_WORLD, 0);
                                exit(0);
                                                        }
                                                      
                                if (verbose == false)   {
                                cout << "Some value does not follow the requirements. " << endl;
                                cout << desc << endl;
                                cout << "Aborting all processors" << endl;
                                MPI_Abort(MPI_COMM_WORLD,  0);
                                exit(0);
                                                        }
                                if ((size) != (px*py))
                                                        {
                                 cout << "The number of processors doesn't match the partitions in x and y" << endl;
                                 cout << "Size of processors = Px · Py " << endl;
                                 cout << "Aborting all processors" << endl;
                                 cout << desc << endl;
                                 MPI_Abort(MPI_COMM_WORLD, 0);
                                 exit(0);
                                                        }
                                                        
            printparameters();                                           
        
        }
        
}
void Model::CalculateParameters()
{
    dx=Lx/(Nx-1);
    dy=Ly/(Ny-1);
    dt=T/(Nt-1);
    Nt=Nt;
    x0=Lx/2;
    y0=Ly/2;
    //Nodes without boundaries
    pNx=Nx-2; 
    pNy=Ny-2; 
    
    //Process size of each processor.
    //Processor 0 will change its size to fit into partitions. 
    
    //Divide the number of total nodes for the partitions in x and y
    sizex=pNx/px;
    sizey=pNy/py;
    
    //processor_x are the nodes in the x direction in the processor
    //processor_y are the nodes in the y direction in the processor
 
    //Processors in coordinates x = 0 store nodes needed to fill all pNx
    if(sizex*px != pNx && rank%px==0)
    {
        processor_x=pNx-(sizex*(px-1))+2;
        C_x=pNx-(sizex*(px-1));
    }
    else{
        processor_x= sizex + 2;
        C_x=pNx-(sizex*(px-1));
    }
    //Processors in coordinates y=0 store nodes needed to fill all pNy
    if(sizey*py != pNy && rank<px)
    {
        processor_y=pNy-(sizey*(py-1))+2;
        C_y=pNy-(sizey*(py-1));
    }
    else{
    processor_y= sizey +2;
    C_y=pNy-(sizey*(py-1));
    }
  
}
void Model::ValidateParameters()
{
    //Validate parameters Nx,Ny,Nt,ax,ay,b,c.
    if (Nx>0 && Ny>0 && Nt>0 && ax>=0 && b>=0 && c>=0 && ay>=0 && px >=0 && py >=0 && Lx >=0 
    && Ly>=0 && T>=0)
    {
        verbose = true;
    }
    else {
        verbose = false;
    }
}
void Model::printparameters()
{
    cout << endl;
    cout << UNDERLINE << "Parameters Burger's Equation" << CLOSEUNDERLINE << endl;
    cout << "ax = " << ax << endl;
    cout << "ay = " << ay << endl;
    cout << "b = " << b << endl;
    cout << "c = " << c << endl;
    cout << UNDERLINE  << "MPI Parameters " << CLOSEUNDERLINE << endl;
    cout << "Px = " << px << endl;
    cout << "Py = " << py << endl;
    cout << "Number of processors = " << size << endl;
    cout << UNDERLINE  << "Grid Parameters" CLOSEUNDERLINE << endl;
    cout << "Nx (Including Boundaries) = " << Nx << endl;
    cout << "Ny (Including Boundaries) = " << Ny << endl;
    cout << "Number of steps (Nt) = " << Nt << endl;
    cout << "dx =  " << dx << " m" <<  endl;
    cout << "dy =  " << dy << " m" <<  endl;
    cout << "dt =  " << dt << " s" <<  endl;
    cout << UNDERLINE  << "Problem Parameters" CLOSEUNDERLINE << endl;
    cout << "Lx =  " << Lx << " m" <<  endl;
    cout << "Ly =  " << Ly << " m" <<  endl;
    cout << "T  =  " << T  << " s" << endl;
}


