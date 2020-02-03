#include "BurgersNova.h"
#include <iostream>
#include "Model.h"
#include <math.h>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <cstdlib>
#define UNDERLINE "\033[4m"
#define CLOSEUNDERLINE "\033[0m"

using namespace std;

BurgersNova::BurgersNova(Model &c)
{
            //Set again rank and size for parallel computing 
            MPI_Comm_rank(MPI_COMM_WORLD,&rank);
            MPI_Comm_size(MPI_COMM_WORLD,&size);
            
            if(rank==0) cout << endl;

            //Get values from Model class    
            setmodel(c);
            //Allocate matrix u (velocity field)
            u=matrixalloc();
            //Set all values to 0
            set0(u);
            //Set Cartesian Grid for processors, also computes neighbour, coordinates.
            Gridprocessors();
}

BurgersNova::~BurgersNova()
{
    //Other arrays are freed in the functions were also called.  
    free(u);
    free(temp);
    free(v);
}

//Store values from Model class.
void BurgersNova::setmodel(Model &c)
{   
    //Gets nodes Nx, Ny and Nt and stores it as X, Y , Z
    X=c.GetNx();
    Y=c.GetNy();
    Z=c.GetNt();
    
    //Gets deltax,deltay,deltaz and the offsetx/offsety
    offsetx=c.GetX0();
    offsety=c.GetY0();
    deltaX=c.GetDx();
    deltaY=c.GetDy();
    deltaT=c.GetDt();
    
    //Declare offset without boundaries
    offsetx=offsetx-deltaX;
    offsety=offsety-deltaY;
    
    //Parameters ax,ay,b,c
    pc=c.GetC();
    pax=c.GetAx();
    pay=c.GetAy();
    pb=c.GetB();
    
    //MPI parameters for the processors
    
    //py and px are the partitions in the x and y domain (how many processors)
    py=c.Getpy(); 
    px=c.Getpx();
    
    //the nodes for each subdomain including halo exchange columns and rows 
    node_x=c.Getprocessorx();
    node_y=c.Getprocessory();
    C_x = c.GetCx();
    C_y = c.GetCy();
}

//Matrix allocation for any 1darry from node_y and node_x
double * BurgersNova::matrixalloc()
{
    //Allocated on the heap. ("new").
    double *m = new double [node_y*node_x];
    return m;
}

 //1d array, row major (node_x is in rows,node_y is in col).
void BurgersNova::set0(double * v)
{
   
    for(int i=0;i<node_y;i++)
    {
        for(int j=0;j<node_x;j++)
        {
            v[i*node_x+j]=0;
            //cout << v[i*node_x+j] << " ";
        }
        //cout << endl;
    }
}

//Declare Cartesian Grid, neighbours and coordinates.
//Declares indices (row,col) for the first node of each processor relative to final velocity field
//This indices will be useful for the gathering.
void BurgersNova::Gridprocessors()
{
    //Declare neighbors null at the beginning
    //Declare dimensions in x,y, period and order 
    dim[0]=py; dim[1]=px;
    period[0]=0; period[1]=0;
    reorder=0;
    
    
    //Creates grid
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Cart_create(MPI_COMM_WORLD,2,dim,period,reorder,&com2d);
    
    //Fills neighbor array for each processor (north,south,east,west).
    MPI_Cart_shift(com2d, 0, 1, &neighbor[N], &neighbor[S]);
    MPI_Cart_shift(com2d, 1, 1, &neighbor[W], &neighbor[E]);
    
    //cout << rank << " " << neighbor[N] << " " << neighbor[S] <<  neighbor[E] << " " << neighbor[W] << endl;
    
    //Get cords for each rank(processor) and stores it in array coord
    MPI_Cart_coords(com2d,rank,2,coord);
    
    MPI_Barrier(com2d);

    MPI_Comm_rank(com2d,&rank);
    
    //Compute indx,indy (first node position relative to final velocity field)(useful for gathering)
    if(coord[1]!=0) indx = 1+(C_x)+(node_x-2)*(coord[1]-1);
    if(coord[1]==0) indx = 1;
    if(coord[0]!=0) indy = 1+(C_y)+(node_y-2)*(coord[0]-1);
    if(coord[0]==0) indy = 1;
    
}
//Initl Function computes radius and the value of the initial condition
double BurgersNova::initilfunction(int posi, int posj)
{
   //Declare variables for this initilfunction
   double f;
   double posx;
   double posy;
   double r;
   
    //Compute posx,posy from posj and posi. 
    posx = (-offsetx)+(indx-1)*deltaX+(posj-1)*deltaX;
    posy = (offsety)-(indy-1)*deltaY-(posi-1)*deltaY;
    
  //Define radius from the positions in x and y 
     r=sqrt(posx*posx+posy*posy);
     
  //Fill field with function f. 
   if (r<=1)
    {
        f=2*pow((1-r),4)*(4*r+1);
    }
    else
    {
        f=0;
    }
    
    return f;
}
//Set initial condition , calls function initilfunction
void BurgersNova::setvalue()
{
    double func;
    for (int lp = 1 ; lp<(node_y-1) ;lp++)
        {
            for (int j=1; j<(node_x-1); j++)
            {
                 func = initilfunction(lp,j);
                 u[lp*(node_x)+j]=func;
                 //cout << u[lp*node_x+j] << " ";
            }
            //cout << endl;
        }
        
}

//Derivatives
double BurgersNova::derx(double * m,int lp, int j)
{
    return ((m[lp*node_x+j]- m[lp*node_x+(j-1)])/(deltaX));
}
double BurgersNova::dery(double * m,int lp, int j)
{
    return ((m[lp*node_x+j]- m[(lp-1)*node_x+j])/(deltaY));
}   
double BurgersNova::deryy(double * m,int lp, int j)
{
    return ((m[((lp+1)*node_x)+j]-2*(m[lp*node_x+j])+m[((lp-1)*node_x)+j])/(deltaY*deltaY));
}
double BurgersNova::derxx(double * m,int lp, int j)
{
    return ((m[(lp*node_x)+(j+1)]-2*(m[(lp*node_x)+j])+m[(lp*node_x)+(j-1)])/(deltaX*deltaX));
}

//Iterativeprocess calls three functions. (exchangeprocessors,calculatevalues and passvalues)
void BurgersNova::iterativeprocess()
{
    temp = matrixalloc();
     int k;
        for (k=0;k<Z;k++)
        {
            exchangeprocessors();
            calculatevalues();
            passvalues();
        }
}
//Exchange rows and columns between neighbour processors.
void BurgersNova::exchangeprocessors()
{    
    //Define Mpi datatype (col)
    //Stripe of node_x, block 1, count node_y-2
    MPI_Type_vector((node_y-2), 1 , node_x , MPI_DOUBLE, &columns);
    MPI_Type_commit(&columns);

    //Sendrecv for North and South. 
    MPI_Sendrecv(&u[1*node_x+1], (node_x-2), MPI_DOUBLE, neighbor[N], 0, 
    &u[(node_y-1)*node_x+1], (node_x-2) , MPI_DOUBLE, neighbor[S], 0, com2d, NULL);
    
    MPI_Sendrecv(&u[(node_y-2)*node_x+1], (node_x-2), MPI_DOUBLE, neighbor[S], 0, 
    &u[0*node_x+1], (node_x-2) , MPI_DOUBLE, neighbor[N], 0, com2d, NULL);
    
    //Sendrecv for West East
    MPI_Sendrecv(&u[1*node_x+(node_x-2)],1,columns, neighbor[E], 0, &u[1*node_x+0], 1,columns,
                neighbor[W], 0, com2d, NULL);
                
    MPI_Sendrecv(&u[1*node_x+1],1, columns, neighbor[W], 0, &u[1*node_x+(node_x-1)],1, columns, 
                neighbor[E], 0, com2d, NULL);
}
//Compute viscous terms, advec1, advec 2
//Only 1 velocity field (initial condition is the same). 
void BurgersNova::calculatevalues()
{
    double visc;
    double advec1;
    double advec2;
    int lp;
    int j;
    
    for (lp=1;lp<=(node_y-2);lp++)
    {  
        for (j=1;j<=(node_x-2);j++)
        {
            visc = (deltaT*pc)*(derxx(u,lp,j)+deryy(u,lp,j));
            advec1 = ((pax+pb*u[lp*node_x+j])*deltaT)*(derx(u,lp,j));
            advec2 = ((pay+pb*u[lp*node_x+j])*deltaT)*(dery(u,lp,j));
            temp[lp*node_x+j]=visc-advec1-advec2+u[(lp*node_x)+j];
        }
    }
    
}
//Pass values from temp vector to u
void BurgersNova::passvalues()
{   
    int lp;
    int j;
    for ( lp=1;lp<=(node_y-2);lp++)
                    {  
                            for (j=1;j<=(node_x-2);j++)
                        {
                            u[lp*node_x+j]=temp[lp*node_x+j];    
                        }
                    } 
}
//Gather all information into rank=0
void BurgersNova::gatherparameters()
{
    //For Gather_v (offsets and elements)
    int recv = (node_x-2)*(node_y-2);
    int*  recvcounts = new int [size];
    int*  displs = new int [size];
    
    //Gather indices of each processor and nodes_x and nodes_y
    int * indxgather = new int [size];
    int * indygather = new int [size];
    int * nodexgather = new int [size];
    int * nodeygather = new int [size];
    
    //Declare vectors useful to write into final velocity field. 
    int gathersize=(X-2)*(Y-2);
    double * gatherv = new double [gathersize];
    double * gatherindividual = new double [recv];
    
    //Allocate final velocity field and set boundaries (finalsize (X)*(Y))
    int finalsize = (X)*(Y);
    v = new double [finalsize];
    boundaries();
    
    //Gather all recv,indx,indy,node_x,node_y to rank 0. 
    MPI_Gather(&recv,1,MPI_INT,recvcounts,1,MPI_INT,0,com2d);
    MPI_Gather(&indx,1,MPI_INT,indxgather,1,MPI_INT,0,com2d);
    MPI_Gather(&indy,1,MPI_INT,indygather,1,MPI_INT,0,com2d);
    MPI_Gather(&node_x,1,MPI_INT,nodexgather,1,MPI_INT,0,com2d);
    MPI_Gather(&node_y,1,MPI_INT,nodeygather,1,MPI_INT,0,com2d);
    
    //Compute offsets to do Gather_v from recvcounts. 
    if(rank==0)
    {   
        displs[0]=0;
        for(int i=1;i<size;i++) displs[i]=displs[i-1]+(recvcounts[i-1]);
    }
    
    //New matrix to not store ghost cells into Gather_v 
    for (int i=1;i<(node_y-1);i++)
    {
        for(int j=1;j<(node_x-1);j++)
        {
            gatherindividual[(i-1)*(node_x-2)+(j-1)]= u[(i)*(node_x)+(j)];
        }
    }
    
    //Gatherv all the processors nodes (without ghost cells) into gatherv
    MPI_Gatherv(&gatherindividual[(0*node_x)+0],recv,MPI_DOUBLE,gatherv,recvcounts,displs,MPI_DOUBLE,0,com2d);

    //Write this into final velocity field using the gathering indices, and nodes for each processor
    if(rank==0)
    {   
        int c=0;
        for (int s=0;s<size;s++)
            {
                for( int i=indygather[s]; i<(indygather[s]+(nodeygather[s]-2));i++)
                {
                   for( int j=indxgather[s]; j<indxgather[s]+(nodexgather[s]-2);j++)
                            {
                                v[i*(X)+j] = gatherv[c];
                                v[i*X+j] = sqrt(2*pow(v[i*X+j],2));
                                c++;
                            } 
                }
            }  
    }
    
    //Free memory allocated on heap
    free(indxgather);
    free(indygather);
    free(nodexgather);
    free(nodeygather);
    free(recvcounts);
    free(displs);
    free(gatherindividual);
    free(gatherv);
}
//Computes boundaries, since are 0, sets all field to 0.
//This is called before writting into matrix v
void BurgersNova::boundaries()
{
    for (int i=0;i<Y;i++)
    {
        for (int j=0;j<X;j++)
        {
            v[i*X+j]=0;
        }
    }
}
//Compute energy and cout in terminal the value. 
void BurgersNova::energy()
{
    if(rank==0)
    {
    double sum=0;
    double e=0;
    for(int i=0;i<Y;i++)
    {
        for(int j=0;j<X;j++)
        {
            e=v[i*X+j]*v[i*X+j]*deltaX*deltaY;
            sum=sum+e;
        }
    }
    
    cout << "Energy of the velocity field is  " << sum/2 << endl;
    }
}
//Outside class (Print Final Velocity into file.txt)
//Overwrite operator. 
ofstream& operator <<(ofstream& ofs, const BurgersNova& b)
{
        for(int lp=0; lp < b.getY() ; lp++ )
        {
            for (int j=0; j<b.getX() ; j++)
            {
                
                ofs << setprecision(5) <<  b.v[(lp*b.getX())+j] << " ";
                
            }
           
           ofs << endl;
           
        }
    return ofs;
}
    
    
    
    
    
    
    
    
    
    


