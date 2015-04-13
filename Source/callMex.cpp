/*
 * MexInterface.cpp
 *
 *  Created on: Feb 12, 2010
 *      Author: christoph vogel
 */
//#define EIGEN_USE_MKL_ALL
//C:\Program Files (x86)\Intel\Composer XE\mkl\include
//C:\Program Files %28x86%29\Intel\Composer XE\mkl\lib\intel64
//mkl_intel_lp64.lib
//mkl_core.lib
//mkl_blas95_ilp64.lib
//mkl_sequential.lib

#include <ctime>
//#define _NO_OPENMP
#include <iostream>
#include "mex.h"
#include <mat.h>
//#include "engine.h"
#include "FullStep.cpp"
#include "nrutil.h"
#include "nrutil.c"

#define zero 1

typedef void (*mexFunction_t) (int nargout, mxArray *pargout[], int nargin, const mxArray *pargin[]);

int main ( int argc, const char *argv[])
{

  char buf[1024];

  MATFile * pmat;

//  pmat = matOpen("testFlow.mat","r");
//  pmat = matOpen("testFlowCSAD.mat","r");
  pmat = matOpen("of_censPic.mat","r");

  const char **dir;
  const char *name;
  int	  ndir;
  mxArray *pa;
  
  dir = (const char **)matGetDir(pmat, &ndir);

  mxFree(dir);

  /* In order to use matGetNextXXX correctly, reopen file to read in headers. */
  if (matClose(pmat) != 0) {
    return(1);
  }

//  pmat = matOpen("testFlow.mat","r");
//  pmat = matOpen("testFlowCSAD.mat","r");
  pmat = matOpen("of_censPic.mat","r");

  if (pmat == NULL) {
    return(1);
  }

  for (int i=0; i < ndir; i++) 
  {
    pa = matGetNextVariableInfo(pmat, &name);
    if (pa == NULL) 
    {
	    return(1);
    }
  }

  if (matClose(pmat) != 0) {
    return(1);
  }

//  pmat = matOpen("testFlow.mat","r");
//  pmat = matOpen("testFlowCSAD.mat","r");
  pmat = matOpen("of_censPic.mat","r");
  

  if (pmat == NULL) {
    return(1);
  }

  const mxArray *argin[20];
  /*
  for (int i=0; i < ndir; i++) 
  {
    argin[i] = matGetNextVariable(pmat, &name);
    if (argin[i] == NULL) 
    {
	    return(1);
    }
  }
*/
  for (int i=0; i < ndir; i++) 
  {
    mxArray* temp = matGetNextVariable(pmat, &name);
    int number = atoi( &name[1] );
    argin[number-1] = temp;
    if (argin[number-1] == NULL) 
	    return(1);
  }


  mxArray *pargout[2] = {0};//,0};

  int nlhs = 1;
  int nrhs = ndir;

  std::clock_t fullStart(std::clock());
  fullStep( nlhs, pargout, nrhs, argin );
  std::clock_t fullStop(std::clock());
  printf("Full data took %f\n", double(fullStop-fullStart)/CLOCKS_PER_SEC );

}