#ifndef BURGERSNOVA_H
#define BURGERSNOVA_H
#include <iostream>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include "Model.h"
#include <mpi.h>
#include <iomanip>
#define underline "\033[4m"

using namespace std;

class Model;
class BurgersNova
{
public:

                        //CONSTRUCTOR, calls Model.h
                         BurgersNova(Model &c);
                         
                        //DECONSTRUCTOR
                         ~BurgersNova();

                        //Public Functions to be called in main.cpp
                        void setvalue();
                        void iterativeprocess();
                        void energy();
                        void gatherparameters();
                        double * v;

                        //Getters
                        int getX()       const  {return X;}
                        int getY()       const  {return Y;}
                        
                        //MPI rank size and com2d
                        MPI_Comm com2d;
                        int rank;
                        int size;
    
private:

                        //Variables stored from Model.h
                        Model c();
                        void setmodel (Model &c);
                        int X;
                        int Y;
                        int Z;
                        double pax;
                        double pay;
                        double pb;
                        double pc;
                        double deltaX;
                        double deltaY;
                        double deltaT;
                        double offsetx;
                        double offsety;
                        int node_x;
                        int px;
                        int node_y;
                        int py;
                        int C_x;
                        int C_y;

                        //Matrix Allocations for Iterative Process and Initial Condition
                        double * u;
                        double * temp;
                        void set0(double *v);
                        double * matrixalloc();
                        double initilfunction(int posi, int posj);
                        
                        //Set boundaries in final velocity
                        void boundaries();
                        
                        //MPI GRID VARIABLES and Function
                        void Gridprocessors();
                        int dim[2];
                        int period[2];
                        int reorder;
                        int coord[2];
                        int S=0, E=1, N=2, W=3;
                        int neighbor[4];
                        int indx;
                        int indy;
                        int recv;
                        
                         //Derivatives
                        double dery(double * m,int lp, int j);
                        double derx(double * m,int lp, int j);
                        double derxx(double * m,int lp, int j);
                        double deryy(double * m,int lp, int j);
                        
                        
                        //Iterative Process functions and MPI_DATATYPE columns
                        MPI_Datatype columns;
                        void exchangeprocessors();
                        void passvalues();
                        void calculatevalues();
};

//Overwrite operator << to print into file. 
ofstream& operator <<(ofstream& ofs, const BurgersNova &b);



#endif // BURGERSNOVA_H
