//  ================================================================================
//  ||                                                                            ||
//  ||              fp                                                            ||
//  ||              -------------------------------------------                   ||
//  ||              F R E E   P A R T I C L E                                     ||
//  ||                                                                            ||
//  ||              D E M O N S T R A T I O N   C O D E                           ||
//  ||              -------------------------------------------                   ||
//  ||                                                                            ||
//  ||       Developed by: Scott R. Runnels, Ph.D.                                ||
//  ||                     University of Colorado Boulder                         ||
//  ||                                                                            ||
//  ||                For: CU Boulder CSCI 4576/5576 and associated labs          ||
//  ||                                                                            ||
//  ||           Copyright 2020 Scott Runnels                                     ||
//  ||                                                                            ||
//  ||                     Not for distribution or use outside of the             ||
//  ||                     this course.                                           ||
//  ||                                                                            ||
//  ================================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include "stdio.h"
#include "math.h"
#include "string.h"
using std :: string;
using std :: vector;
using std :: stringstream;
using std :: cout;
using std :: endl;

#define rLOOPf for ( int r = 1 ; r <= nrowsf ; ++r )
#define rLOOP  for ( int r = 1 ; r <= nrows  ; ++r )
#define cLOOPf for ( int c = 1 ; c <= nrowsf ; ++c )
#define cLOOP  for ( int c = 1 ; c <= nrows  ; ++c )
#define iLOOP  for ( int i = 1 ; i <= nRealx ; ++i )
#define iLOOPi for ( int i = 2 ; i <= nx-1   ; ++i )
#define jLOOP  for ( int j = 1 ; j <= nRealy ; ++j )
#define jLOOPi for ( int j = 2 ; j <= ny-1   ; ++j )

#define sLOOP  for ( int s = 0 ; s <= nRealx+1 ; ++s )
#define tLOOP  for ( int t = 0 ; t <= nRealy+1 ; ++t )

#define pLOOP  for ( int p = 1 ; p < npts*npts ; ++p )

#define D_ double*
#define _D_ double
#define _I_ int
#define VD  vector<double>
#define VDD vector<vector<double> >
#define VI  vector<int>

#define Here   I.pid (iC   , jC   ) 
#define West   I.pid (iC-1 , jC   ) 
#define East   I.pid (iC+1 , jC   ) 
#define North  I.pid (iC   , jC+1 ) 
#define South  I.pid (iC   , jC-1 ) 
#define NW     I.pid (iC-1 , jC+1 ) 
#define NE     I.pid (iC+1 , jC+1 ) 
#define SE     I.pid (iC-1 , jC-1 ) 
#define SW     I.pid (iC+1 , jC-1 ) 


void Exit()
{
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(0);
}


void FatalError(string msg)
{

  cout << " " << endl;
  cout << " " << endl;
  cout << "Fatal Error on this PE: " << msg << endl;
  cout << " " << endl;
  Exit();


}


