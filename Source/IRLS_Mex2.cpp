/*
 * Sort for Census Transform.cpp
 *
 *  Created on: Dec 1, 2009
 *      Author: Christoph Vogel
 */

//#define _NO_OPENMP

#include <iostream>
#include <string.h>

#include "mex.h"
#include "nrutil.h"
#include "nrutil.c"

#include "FullStep.cpp"


// mex CXX=g++-4.2 CXX=g++-4.2 LD=g++-4.2 -I/home/christop/CPP -lm -output mexTest IRLS_mex.cpp
// mex CC=g++-4.2 CXX=g++-4.2 CFLAGS='-fPIC -fopenmp' LD=g++-4.2 -O -I/home/christop/CPP -lm -lgomp -output IRLS_mex IRLS_mex.cpp
void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
    fullStep( nlhs, plhs, nrhs, prhs );
}
