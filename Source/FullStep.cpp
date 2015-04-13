/*
Copyright (c) 2013, Christoph Vogel, ETH Zurich

The code may be used free of charge for non-commercial and
educational purposes, the only requirement is that this text is
preserved within the derivative work. For any other purpose you
must contact the authors for permission. This code may not be
redistributed without written permission from the authors.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES 
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF 
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE 
FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY 
DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, 
WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, 
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

//mex CXX=g++-4.2 CXX=g++-4.2 LD=g++-4.2 -I/home/christop/CPP -lm -output mexTest IRLS_mex.cpp

// less memory consumption but slower
//#define _noeigen_

#include <vector>
#include <algorithm>
#include <omp.h>
#include "mex.h"

#include "BicubicIp.h"
#include "CensusImgMask.h"
#include "CSAD_Parameter.h"
#ifdef _noeigen_
  #include "TGV_SparseMatrix_own.h"
#else
  #include "TGV_SparseMatrix.h"
#endif
#include "MedSimple.h"

/// matlab calling
void fullStep ( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{

  typedef double Scalar;

  Scalar *output1;

  Scalar *input1 = (Scalar*)   mxGetPr(prhs[0]);//image1
  Scalar *input2 = (Scalar*)   mxGetPr(prhs[1]);//image2
  int ring       = (int)  (*mxGetPr(prhs[2]));//ringsize
  Scalar *u      = (Scalar*)   mxGetPr(prhs[3]);//u
  Scalar *v      = (Scalar*)   mxGetPr(prhs[4]);//v
  Scalar  lambda = (Scalar ) (*mxGetPr(prhs[5]));//lambda

  // matrix
  Scalar* ids_i     = (Scalar*) mxGetPr(prhs[6]); // 
  Scalar* ids_j     = (Scalar*) mxGetPr(prhs[7]); // 
  Scalar* ids_v     = (Scalar*) mxGetPr(prhs[8]); // 

  Scalar  n_rows    = (Scalar) (*mxGetPr(prhs[9]));
  Scalar  n_cols    = (Scalar) (*mxGetPr(prhs[10]));

  Scalar* uvw       = (Scalar*) mxGetPr(prhs[11]); // uvw
  Scalar* pq        = (Scalar*) mxGetPr(prhs[12]); // pq

  int warps         = (int)  (*mxGetPr(prhs[13]));// outer iterations
  int its           = (int)  (*mxGetPr(prhs[14]));// inner iterations

  Scalar cEps = 0.15/255.;
  int doStereo = 0;
  int doCensus = 0;
  int doTV = 0;

  if (nrhs > 15)
    cEps         = (Scalar)  (*mxGetPr(prhs[15]));// for census
  if (nrhs > 16)
    doStereo     = (int)     (*mxGetPr(prhs[16]));// do stereo == ignore y direction
  if (nrhs > 17)
    doCensus     = (int)     (*mxGetPr(prhs[17]));// do census=1, csad=2, sad=0
  if (nrhs > 18)
    doTV         = (int)     (*mxGetPr(prhs[18]));// do total variation instead

  // constraints are treated in rhs of matrix and weighted already - therefore l1 constraints
  // sum w_i u_i = u_constraint; with sum_i w_i = 1* weight of constraint
  Scalar* constraints = NULL;
  int nConstraints = 0;
  if (nrhs > 19)
  {
    constraints  = (Scalar *) (mxGetPr(prhs[19]));// positional constraints on flow field; simply already in matrix
    nConstraints=  (int)mxGetNumberOfElements(prhs[19]);
  }

  int doPf = 1;
  if (nrhs > 20) // 0: off no median filter, 1: filter peaks only 2: filter all no matter what
    doPf          = (int)     (*mxGetPr(prhs[20]));// do median filtering

  if(doCensus==0)
    ring =0;

  int elements=(int)mxGetNumberOfElements(prhs[7]);

  ///////////////////////
  int nsubs=mxGetNumberOfDimensions(prhs[0]);
  const mwSize* dims = mxGetDimensions(prhs[0]);
  int pixel  = dims[0]*dims[1];
  int dx = dims[0];
  int dy = dims[1];

  censusMaskImaging<Scalar> cmib(dx, dy, ring, lambda, cEps);
  peakFilter<Scalar> pf( dx, dy, 1 );// 3x3 filter

  cmib.setEps( cEps );//0.15/255 );
  if (doCensus==1) // CENSUS
    cmib.computeEvents( input1 );
  
  if (doCensus==2) // CSAD
    cmib.compute( input1 );

//  printf("cEps: %.2f, doSt: %d, doCen:%d  ", 255*cEps, doStereo, doCensus);

  Scalar* tau_ = (Scalar*) malloc( sizeof(Scalar)*pixel );
#ifdef _noeigen_
  ownMatrixVecor<Scalar> solver( n_rows, n_cols, uvw, pq, elements, ids_i, ids_j, ids_v, doTV, constraints, nConstraints );
#else
  eigenMatrixVecor<Scalar> solver( n_rows, n_cols, uvw, pq, elements, ids_i, ids_j, ids_v, doTV, constraints, nConstraints );
#endif
  solver.store_uvw();
  solver.copy_tau( tau_ );

  CsadParameter<Scalar> csadP (dx, dy, ring, lambda, tau_, input2, u, v );
  int ringsize = (2*ring+1)*(2*ring+1)-1;
  csadP.setbImg( cmib.getbImg() );

  for (int www =0; www < warps; www++ )
  {
    std::vector<Scalar> temp_u; 
    std::vector<Scalar> temp_v; 
    solver.extract_uv (temp_u, temp_v);

    csadP.compute_Indices_mask( temp_u, temp_v );
    if (doCensus==1) // CENSUS
    {
      csadP.interpolateCensus( cmib.getEps(), doStereo );
      csadP.preSortDataCensus( cmib.getEvents() );
    }
    if (doCensus==2) // CSAD
    {
      csadP.interpolate( doStereo );
      csadP.preSortData();
    }
    if (doCensus==0) // SAD
      csadP.interpolateSAD( input1, doStereo );


    for (int iii =0; iii < its; iii++ )
    {
      solver.step1( doTV );
      solver.extract_uv ( temp_u, temp_v );


      if (doCensus==1) // TCENSUS
//        csadP.dataStepCensus_new( temp_u, temp_v ); // why worse ?
        csadP.dataStepCensus( temp_u, temp_v, doStereo );
      if (doCensus==2)// CSAD
        csadP.dataStep( temp_u, temp_v );
      if (doCensus==0) // SAD
        csadP.dataStepSAD( temp_u, temp_v );

      if (iii == its-1 && doPf)
      {
        pf.compute( temp_u, doPf>1 );
        pf.compute( temp_v, doPf>1 );
      }
      solver.set_uv( temp_u, temp_v );
      solver.store2_uvw();
    }
  }

  plhs[0] = mxCreateDoubleMatrix( dx, dy, mxREAL);
  output1  = mxGetPr(plhs[0]);
  solver.copy_u ( output1 );

  plhs[1] = mxCreateDoubleMatrix( dx, dy, mxREAL);
  output1  = mxGetPr(plhs[1]);
  solver.copy_v ( output1 );

  mwSize dimsOut[3];
  dimsOut[0] = dx;dimsOut[1] = dy;dimsOut[2] = 4;
  plhs[2]    =  mxCreateNumericArray(3, dimsOut, mxDOUBLE_CLASS, mxREAL);
  output1  = mxGetPr(plhs[2]);
  if (!doTV)
    solver.copy_w ( output1 );
  else
    memset( output1, 0, sizeof(double)*dimsOut[0] * dimsOut[1] * dimsOut[2] );
//  dimsOut[0] = dx;dimsOut[1] = dy;dimsOut[1] = 4;
  plhs[3]    =  mxCreateNumericArray(3, dimsOut, mxDOUBLE_CLASS, mxREAL);
  output1  = mxGetPr(plhs[3]);
  solver.copy_p ( output1 );

  dimsOut[2] = 8;
  plhs[4]    =  mxCreateNumericArray(3, dimsOut, mxDOUBLE_CLASS, mxREAL);
  output1  = mxGetPr(plhs[4]);
  if (!doTV)
    solver.copy_q ( output1 );
  else
    memset( output1, 0, sizeof(double)*dimsOut[0] * dimsOut[1] * dimsOut[2] );
  free(tau_);
}