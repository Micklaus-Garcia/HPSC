//  ==
//  ||
//  ||   C L A S S:    m p i I n f o
//  ||
//  ==
class mpiInfo {
  public:

  int myPE;
  int numPE;
  int nRealx,  nRealy;
  int nx,  ny;
  int nPEx, nPEy;
  int iPE , jPE;
  int iMin, iMax, jMin, jMax ; // The global i-j numbers on this processor
  int nei_n, nei_s, nei_e, nei_w;
  
  int countx, county;

  double *phiL, *phiR;
  double *phiT, *phiB;
    
  double *phiSend_n,  *phiSend_s;
  double *phiSend_e,  *phiSend_w;
  double *phiRecv_n,  *phiRecv_s, *phiRecv_e,  *phiRecv_w;
  
  MPI_Status  status;
  int         err;
  int         tag;
  MPI_Request request;

  //  -
  //  |
  //  |   GridDecomposition: Set up PE numbering system in figure below and
  //  |                      establish communication arrays.
  //  |
  //  |                      nPEx -- number of PEs in the x-direction
  //  |                      nPEy -- number of PEs in the y-direction
  //  |                      numPE = total number of PEs
  //  |
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       |       |       |         |       |
  //  |                       |       |       |         | numPE |
  //  |                       |       |       |         |       |
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       .       .       .         .       .
  //  |                       .       .       .         .       .
  //  |                       .       .       .         .       .
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       |       |       |         |       |
  //  |                       | nPEx  | nPEx+1|         |       |
  //  |                       |       |       |         |       |
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       |       |       |         |       |
  //  |                       |   0   |   1   |         | nPEx-1|
  //  |                       |       |       |         |       |
  //  |                       +-------+-------+ . . .   +-------+
  //  |
  //  |
  //  -
  void GridDecomposition(int _nPEx, int _nPEy, int nCellx, int nCelly) {

    nRealx = nCellx + 1;
    nRealy = nCelly + 1;

    // Store and check incoming processor counts 
    nPEx = _nPEx;
    nPEy = _nPEy;
    
    if (nPEx*nPEy != numPE) {
      if ( myPE == 0 ) cout << "Fatal Error:  Number of PEs in x-y directions do not add up to numPE" << endl;
    	MPI_Barrier(MPI_COMM_WORLD);
    	MPI_Finalize();
    	exit(0);
    }
    
    // Get the i-j location of this processor, given its number.  See figure above: 
    jPE = int(myPE/nPEx);
    iPE = myPE - jPE*nPEx;

    // Set neighbor values
    nei_n = nei_s = nei_e = nei_w = -1;

    if ( iPE > 0      ) { nei_w = myPE - 1    ; }
    if ( jPE > 0      ) { nei_s = myPE - nPEx ; }
    if ( iPE < nPEx-1 ) { nei_e = myPE + 1    ; }
    if ( jPE < nPEy-1 ) { nei_n = myPE + nPEx ; }

    // Acquire memory for the communication between adjacent processors:
    countx = nRealx + 2;
    county = nRealy + 2;
    
    phiL = new double [ county ];	
    phiR = new double [ county ];	
    phiT = new double [ countx ];	
    phiB = new double [ countx ];
    
    phiSend_n = new double [ countx ];
    phiSend_s = new double [ countx ];
    phiSend_e = new double [ county ];
    phiSend_w = new double [ county ];
    phiRecv_n = new double [ countx ];
    phiRecv_s = new double [ countx ];
    phiRecv_e = new double [ county ];
    phiRecv_w = new double [ county ];

    tag = 0;
  }
  void PEsum(VD &field) {

      sLOOP phiSend_n[s] = 0.;
	sLOOP phiSend_s[s] = 0.;
	tLOOP phiSend_e[t] = 0.;
	tLOOP phiSend_w[t] = 0.;
 
      tLOOP phiL[t] = field[pid(1     , t     )];
      tLOOP phiR[t] = field[pid(nRealx, t     )];
      sLOOP phiT[s] = field[pid(s     , nRealy)];
      sLOOP phiB[s] = field[pid(s     , 1     )];

      sLOOP phiSend_n[s] = phiT[s];
      sLOOP phiSend_s[s] = phiB[s];
      tLOOP phiSend_e[t] = phiR[t];
      tLOOP phiSend_w[t] = phiL[t];
       
      if (nei_n>=0) err = MPI_Isend(phiSend_n, countx, MPI_DOUBLE, nei_n, tag, MPI_COMM_WORLD, &request);
      if (nei_s>=0) err = MPI_Isend(phiSend_s, countx, MPI_DOUBLE, nei_s, tag, MPI_COMM_WORLD, &request);
      if (nei_e>=0) err = MPI_Isend(phiSend_e, county, MPI_DOUBLE, nei_e, tag, MPI_COMM_WORLD, &request);
      if (nei_w>=0) err = MPI_Isend(phiSend_w, county, MPI_DOUBLE, nei_w, tag, MPI_COMM_WORLD, &request);

      if (nei_n>=0) {
      err = MPI_Irecv(phiRecv_n, countx, MPI_DOUBLE, nei_n, tag, MPI_COMM_WORLD, &request);   
            MPI_Wait(&request, &status); 
      }
      if (nei_s>=0) {
      err = MPI_Irecv(phiRecv_s, countx, MPI_DOUBLE, nei_s, tag, MPI_COMM_WORLD, &request);   
            MPI_Wait(&request, &status); 
      }
      if (nei_e>=0) {
      err = MPI_Irecv(phiRecv_e, county, MPI_DOUBLE, nei_e, tag, MPI_COMM_WORLD, &request);   
            MPI_Wait(&request, &status); 
      }
      if (nei_w>=0) {
      err = MPI_Irecv(phiRecv_w, county, MPI_DOUBLE, nei_w, tag, MPI_COMM_WORLD, &request);   
            MPI_Wait(&request, &status); 
      }

      if (nei_n>=0) sLOOP phiT[s] = phiRecv_n[s];
      if (nei_s>=0) sLOOP phiB[s] = phiRecv_s[s];
      if (nei_e>=0) tLOOP phiR[t] = phiRecv_e[t];
      if (nei_w>=0) tLOOP phiL[t] = phiRecv_w[t];

      if (nei_n>=0) sLOOP field[pid( s        , nRealy )] += phiT[s];
      if (nei_s>=0) sLOOP field[pid( s        , 1      )] += phiB[s];
      if (nei_e>=0) tLOOP field[pid( nRealx   , t      )] += phiR[t];
      if (nei_w>=0) tLOOP field[pid( 1        , t      )] += phiL[t];
   }

   int pid(int i,int j) {return (i+1) + (j)*(nRealx+2);} 
};

