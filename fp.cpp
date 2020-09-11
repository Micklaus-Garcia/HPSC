
#include "mpi.h"
#include "fp.h"
#include "particles.h"
#include "mpiInfo.h"


class Mesh
{

public:

  double x0, x1, y0, y1;
  VD x,y;
  VD Qval;
  int nRealx    , nRealy   , nField, nReal;
  double lengthx, lengthy, dx, dy;
  int myPE;

  //  ==
  //  ||
  //  ||  Constructor:  Initialize values
  //  ||
  //  ==

  Mesh(double _x0  , double _x1,  double _y0, double _y1 , int ncell_x , int ncell_y, mpiInfo &myMPI )
  {
    // Copy incoming values

    x0   = _x0;          x1 = _x1;
    y0   = _y0;          y1 = _y1;
    myPE = myMPI.myPE;

    // Compute number of real (physical) nodes, and their spacing
    
    nRealx      = ncell_x + 1; 
    nRealy      = ncell_y + 1;
    nReal       = nRealx*nRealy;
    dx          = (x1-x0) / ncell_x;
    dy          = (y1-y0) / ncell_y; 

    // Compute the size of the field variables, which also lie on ghost nodes
    
    nField      = (nRealx+2)*(nRealy+2);

    x.resize(nField+1); y.resize(nField+1);  Qval.resize(nField+1);

    for ( int i = 1 ; i <= nRealx ; ++i )
      for ( int j = 1 ; j <= nRealy ; ++j )
	{
	  int p = pid(i,j);
	  x[p] = x0 + (i-1)*dx;
	  y[p] = y0 + (j-1)*dy;
	}
  }
    
  //  ==
  //  ||
  //  ||  ParticlesOnMesh
  //  ||
  //  ||
  //  ==

  void ParticlesOnMesh(particles &PTCL, mpiInfo &myMPI)
  {
    double hx, hy;
    double w[5];
    int    p[5];
    int    iL, iR, jB, jT;
    
    // -
    // |
    // | Determine which particles are still on this mesh and which have left
    // |
    // -

    for ( int k = 1 ; k <= PTCL.n ; ++k )
      {
	if ( PTCL.active[k] == 1 )
	  {
	    if ( PTCL.x[k] < x0     ) PTCL.active[k] = -1; 
	    if ( PTCL.x[k] > x1     ) PTCL.active[k] = -1; 
	    if ( PTCL.y[k] < y0     ) PTCL.active[k] = -1; 
	    if ( PTCL.y[k] > y1     ) PTCL.active[k] = -1; 
	  }

      }
    
    // -
    // |
    // | Accumulate particles to the nodes
    // |
    // -

    for ( int k = 1 ; k <= nField ; ++k ) Qval[k] = 0.;

    for ( int k = 1 ; k <= PTCL.n ; ++k )
      {
	if ( PTCL.active[k] == 1 )
	  {
	    iL = int ( (PTCL.x[k]-x0) / dx ) + 1;    // point to the left
	    jB = int ( (PTCL.y[k]-y0) / dy ) + 1;    // point below 
	    iR = iL + 1;                             // point to the right
	    jT = jB + 1;                             // point above
	    
	    // 1-2-3-4 are temporary numbers that refer to the lower left, lower right,
	    // upper right, upper left, in that order.
	    
     	    p[1] = pid(iL,jB);
     	    p[2] = pid(iR,jB);
     	    p[3] = pid(iR,jT);
     	    p[4] = pid(iL,jT);
	    
     	    // Compute weights for spreading particle k to the four surrounding nodes.
     	    // Here, hx is the fractional x-distance from particle k to node 1.  Similar for hy.
     	    // See [1].
	    
     	    hx = (PTCL.x[k] - x[ p[1] ])/dx;
     	    hy = (PTCL.y[k] - y[ p[1] ])/dy;
	    
     	    w[1] = ( 1. - hx ) * ( 1. - hy );
     	    w[2] =        hx   * ( 1. - hy );
     	    w[3] =        hx   *        hy  ;
     	    w[4] = ( 1. - hx ) *        hy  ;
	    
     	    // Spread particle k to the four surrounding points
	
	    for ( int i = 1 ; i <= 4 ; ++i ) Qval[ p[i] ] +=  w[i] * PTCL.Qp;
     	  }
       }

    // -
    // |
    // | Sum Qval on processor boundaries
    // |
    // -

    //    myMPI.PEsum(MESH.Qval);
    //    MESH.plotter("Qval", 0, myMPI);
                                          // <--- In Lab:  Write this routine
        
  }

  
#include "mesh_plotter.h"
  
  int pid(int i,int j) { return (i+1) + (j)*(nRealx+2); }
  
};





//  ==
//  ||
//  ||
//  ||  Main Program
//  ||
//  ||
//  ==
//

int main(int argc, char *argv[])
{
  
   mpiInfo myMPI;
   MPI_Init     (&argc         , &argv       );
   MPI_Comm_size(MPI_COMM_WORLD, &myMPI.numPE);
   MPI_Comm_rank(MPI_COMM_WORLD,&myMPI.myPE  );

   int nPEx, nPEy, nCellx, nCelly;
   double tEnd, dt;
   double flux;

   for (int count =  0 ; count < argc; ++count)
     {
       if ( !strcmp(argv[count],"-nPEx"  ) ) nPEx   = atoi(argv[count+1]);
       if ( !strcmp(argv[count],"-nPEy"  ) ) nPEy   = atoi(argv[count+1]);
       if ( !strcmp(argv[count],"-nCellx") ) nCellx = atoi(argv[count+1]);
       if ( !strcmp(argv[count],"-nCelly") ) nCelly = atoi(argv[count+1]);
       if ( !strcmp(argv[count],"-flux"  ) ) flux   = atof(argv[count+1]);
     }

   // -
   // |
   // | MPI / Processor ID
   // |
   // -
   
   myMPI.GridDecomposition(nPEx,nPEy,nCellx,nCelly);

   // -
   // |
   // | Parallel Grid Generation
   // |
   // -
   
   double totalLength = 1.;
   double eachPElength_x = totalLength / nPEx;
   double eachPElength_y = totalLength / nPEy;

   double x0 = eachPElength_x * myMPI.iPE;   double x1 = x0 + eachPElength_x;
   double y0 = eachPElength_y * myMPI.jPE;   double y1 = y0 + eachPElength_y;
   
    Mesh MESH( x0  , x1 , y0 , y1 , nCellx , nCelly , myMPI );

   // -
   // |
   // | Set up Particles
   // |
   // -
   
//   particles PTCL(500); 
     myMPI.PEsum(MESH.Qval);
     MESH.plot("Qval", 0, myMPI);

   // -
   // |
   // | Time Marching Loop
   // |
   // -
   
//   for ( double t = 0. ; t <= tEnd ; t += dt )
//    {
//      //     cout << endl;
//    cout << "myPE: " << myMPI.myPE << " Time = " << t << endl;

       // Inject particles
       
       // Move particles
       
       // Map between particles and the mesh

       // Plot

    // }


   // -
   // |
   // | Wrap-Up
   // |
   // -

   
   MPI_Finalize();
   
   return 0;

}
