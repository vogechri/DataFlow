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

#ifndef _CSAD_PAR_H
#define _CSAD_PAR_H

#include <vector>
#include <algorithm>

#include "BicubicIp.h"
#include "CensusImgMask.h"
#include "quickSelect.cpp"

//function warpParameter = computeWarpParametersCSAD (b1Img, I2, u,v, par, lambda)
//warpParameter.Nids     = Nids; % per pixel a number
//warpParameter.u_i      = u_i;  % 
//warpParameter.dx       = dx; // gradient dx
//warpParameter.dy       = dy; // gradient dy
//warpParameter.lambdanormG = lambda.*normG;
//warpParameter.u0 = u; // old val
//warpParameter.v0 = v;

// needs censusmaskImg
//ring = par.ring;
//lambda = lambda *8 /((2*(ring+1)-1)^2-1);

// 1. gradients from warped positions, id+flow(1,2, ..) y grad =0 if stereo
// 2. 
//  u_i  = bsxfun( @rdivide,  bImg-b1Img, normG );
//  Nids = sum(~bMask, 3);% now can compute elements of first part
//  u_i(bMask) = 0; % those are ignored in the mex file,

// 3. u_i =preSortMed_quick(reshape( u_i, [numel(Nids), size(u_i, 3)])', Nids(:));

// ideally keep in memory (as a class)


//b1Img, I2, u,v, par, lambda

template <long i> class Tag{};

template<long i, long N, long P, typename T>
inline void partial_insertionsort_for(T a[], Tag<N>, Tag<i>)
{   long j = i <= P+1 ? i : P+1;  // partial sort
    T temp = a[i];
    a[i] = a[j];       // compiler should optimize this away where possible
    while(temp < a[j - 1] && j > 0)
    { a[j] = a[j - 1];
      j--;}
    a[j] = temp;
    partial_insertionsort_for<i+1,N,P,T>(a,Tag<N>(),Tag<i+1>());
};

template<long i, long N, long P, typename T>
inline void partial_insertionsort_for(T a[], Tag<N>, Tag<N>){};

template <long N, long P, typename T>
inline void partial_insertionsort(T a[])
 {partial_insertionsort_for<0,N,P,T>(a, Tag<N>(), Tag<0>());};


template <typename S_>
class CsadParameter
{
  public:

  typedef S_ Scalar;

  CsadParameter(int _M, int _N, int _ring, Scalar _lambda, Scalar*& _tau, Scalar*& _Img2, Scalar*& _u, Scalar*& _v ) 
    :Img2(_Img2), N(_N), M(_M), ring(_ring), lambda(_lambda), tau(_tau), bcIP(_N, _M, _Img2), cmi( _M, _N, _ring ), 
    ringsize((2*_ring+1)*(2*_ring+1)-1), bestSpotsInit(false)
  {
    eventScore = lambda/ringsize;
    generateMem(_u,_v);
  };

  ~CsadParameter()
  {
    freeMem();    
//    deleteMatrices();
  };

  void generateMem( Scalar*& _u0, Scalar*& _v0 )
  {
    u0.resize(N*M);
    std::copy (_u0, &_u0[N*M], u0.begin());
    v0.resize(N*M);
    std::copy (_v0, &_v0[N*M], v0.begin());

//    ringsize = (2*ring+1)*(2*ring+1)-1);
    lambda = lambda * Scalar(8.) / std::max(Scalar(8.), Scalar(ringsize));//((2*ring+1)*(2*ring+1)-1);

    lambdanormG.resize(N*M,0.);
    nIds.resize(N*M,0.);
    u_i.resize(N*M * ringsize, 0. );

    bcIP.approximateImageGradients();
//        bcIP.interpolate( elementsX, X, Y, outputF, outputdFX, outputdFY );
    // index memory: 
    compute_Indices_mask( u0, v0 );
  };

  void freeMem(){};

  /// the mask is defined by the indices not within the image coordinates
//  void compute_Indices_mask( Scalar*& _u0, Scalar*& _v0 )
  void compute_Indices_mask( std::vector<Scalar>& _u0, std::vector<Scalar>& _v0 )
  {
    u0.resize(N*M);
    std::copy (_u0.begin(), _u0.end(), u0.begin());
    v0.resize(N*M);
    std::copy (_v0.begin(), _v0.end(), v0.begin());

    Iu.resize(N*M);Iv.resize(N*M);mask.resize(N*M);maskGrad.resize(N*M);

    int pos = 0;
    for(int j=0; j<N;j++)
      for(int i=0; i<M;i++)
      {
        //twist (see below):
        Iu[pos] = _v0[pos] + (i+1);// 1,2,3,4,5,6,7,8,1,2,...
        Iv[pos] = _u0[pos] + (j+1);// 1111,..,2222,..,333, ....
        mask[pos]     = (Iu[pos]<1) || (Iu[pos]>M)   || (Iv[pos]<1) || (Iv[pos]>N);
        maskGrad[pos] = (Iu[pos]<2) || (Iu[pos]>M-1) || (Iv[pos]<2) || (Iv[pos]>N-1);
        pos++;
      }
  }

  /// interpolates image and computes gradients, stores dx,dy, normG, etc - but not? all 
  void interpolate( int doStereo =0 )
  {
    int patchSize = (2*ring+1)*(2*ring+1)-1;
    std::vector<Scalar> val(N*M);
    dx.resize(N*M);dy.resize(N*M);
    bcIP.interpolate( Iv, Iu, val, dx, dy );//twist (see above)

    if (doStereo)
      dy.assign(N*M, 0.);

    // normal mask!
    cmi.compute( val, mask);

    // done by class censusMask:
    //[bImg, bMask] = computeCensusMask( val, ring, mask, 1 );
    //bMask = bsxfun(@or, bMask, mask );
    const std::vector<Scalar>    bImg2 = cmi.getbImg();//(N*M*patchSize);
    const std::vector<int>       bmask = cmi.getbMask();//(N*M*patchSize);
    // now compute csad stuff:

    //dx(maskGrad) = 0;
    //dy(maskGrad) = 0;
    //normG = max( 0.00001, sqrt(I2_y_warped.^2 + I2_x_warped.^2) );
    //dx = I2_x_warped./normG;
    //dy = I2_y_warped./normG;

    //u_i  = bsxfun( @rdivide,  bImg-b1Img, normG );
    //Nids = sum(~bMask, 3);% now can compute elements of first part
    //u_i(bMask) = 0; % those are ignored in the mex file,

    u_i.assign(N*M*patchSize, 0.0);
    nIds.assign(N*M, 0.0);
    lambdanormG.assign(N*M, 0.0);
//    dy.assign(N*M, 0.0);
//    dx.assign(N*M, 0.0);
    // apply masks on gradients, compute bimg2
    int pos = 0;
    for(int j=0; j<N;j++)
      for (int i=0;i<M;i++)
      {
        if(!maskGrad[pos])
        {
          Scalar normG = sqrt(dx[pos]*dx[pos] + dy[pos]*dy[pos]);
          if (normG > 0.00001)
          {
            dx[pos] = dx[pos] / normG;

            if (!doStereo)
              dy[pos] = dy[pos] / normG;

            lambdanormG[pos] = normG*lambda*tau[pos];
            //            nIds[pos] = 0;
            for ( int k=0;k<patchSize;k++ )
            {
              if( !bmask[patchSize*pos+k] )
              {
                u_i[patchSize*pos+k] = (bImg2[patchSize*pos+k]-b1Img[patchSize*pos+k]) / normG;
                nIds[pos]++;
              }
              //            else
              //              u_i[patchSize*pos+k] = 0.;
            }
            // apply mask
            //lambdaNormG[pos] = lambda * normG;
          }
        }
        else
        {
          dx[pos]=0.0;
          dy[pos]=0.0;
        }
        pos++;
      }
  }


  /// interpolates image and .. ?
  void interpolateSAD( Scalar *& img1, int doStereo =0 )
  {
    std::vector<Scalar> val(N*M);
    dx.resize(N*M);dy.resize(N*M);
    bcIP.interpolate( Iv, Iu, val, dx, dy );//twist (see above)

    if (doStereo)
      dy.assign(N*M, 0.);

    u_i.assign(N*M, 0.0);
    lambdanormG.assign(N*M, 0.0);

    // apply masks on gradients, compute bimg2
    int pos = 0;
    for(int j=0; j<N;j++)
      for (int i=0;i<M;i++)
      {
        if(!mask[pos])
        {
          Scalar normG = sqrt(dx[pos]*dx[pos] + dy[pos]*dy[pos]);
          if (normG > 0.00001)
          {
            dx[pos] = dx[pos] / normG;

            if (!doStereo)
              dy[pos] = dy[pos] / normG;

            lambdanormG[pos] = normG*lambda*tau[pos];
            u_i[pos] = (val[pos]-img1[pos]) / normG;
          }
        }
        else
        {
          dx[pos]=0.0;
          dy[pos]=0.0;
        }
        pos++;
      }
  }

//events( cat(3, bMask, bMask ) ) = 0;
//
//%there are points still indside but whose cells are outside
//normG = max( 0.00001, sqrt(I2_y_warped.^2 + I2_x_warped.^2) );
//dx = I2_x_warped./normG;
//dy = I2_y_warped./normG;
//
//% intersection points: deltaM < deltaP
//deltaP = bsxfun(@rdivide,  cEps - bImg, normG );
//deltaM = bsxfun(@rdivide, -cEps - bImg, normG );
//
//deltaP = bsxfun( @times, deltaP, 1-mGrad );
//deltaM = bsxfun( @times, deltaM, 1-mGrad );
//deltas = cat(3, deltaP, deltaM );

  /// interpolates image and computes gradients, stores dx,dy, normG, etc - also deltas of events
  void interpolateCensus( Scalar cEps, int doStereo = 0 )
  {
    int patchSize = (2*ring+1)*(2*ring+1)-1;
    std::vector<Scalar> val(N*M);
    dx.resize(N*M);dy.resize(N*M, 0.);
    bcIP.interpolate( Iv, Iu, val, dx, dy );//twist (see below)

    // normal mask!
    cmi.compute( val, mask);

    const std::vector<Scalar>    bImg2 = cmi.getbImg();
    const std::vector<int>       bmask = cmi.getbMask();
    // now compute csad stuff:

    deltas.assign(N*M*patchSize*2, 0.0);
    nIds.assign(N*M, 0.0);
    lambdanormG.assign(N*M, 0.0);

    if (doStereo)
      dy.assign(N*M, 0.);

    // apply masks on gradients, compute bimg2
    int pos = 0;
    for(int j=0; j<N;j++)
      for (int i=0;i<M;i++)
      {
        if(!maskGrad[pos])
        {
          Scalar normG = sqrt(dx[pos]*dx[pos] + dy[pos]*dy[pos]);
          if (normG > 0.00001)
          {
            dx[pos] = dx[pos] / normG;

            if (!doStereo)
              dy[pos] = dy[pos] / normG;

            // tau later
            lambdanormG[pos] = normG;//*lambda*tau[pos];
            //            nIds[pos] = 0;
            for ( int k=0;k<patchSize;k++ )
            {
              if( !bmask[patchSize*pos+k] )
              {
                Scalar temp = (bImg2[patchSize*pos+k]);//-b1Img[patchSize*pos+k]);
                // 1. intesity diff becomes > eps if delta*normG + temp > cEps <-> delta = 
                deltas [2*patchSize*pos+k          ] = ( cEps - temp) / normG;
                // 2. intesity diff becomes <-eps if delta*normG + temp <-cEps <-> delta = 
                deltas [2*patchSize*pos+k+patchSize] = (-cEps - temp) / normG;
//                u_i[patchSize*pos+k] = (bImg2[patchSize*pos+k]-b1Img[patchSize*pos+k]) / normG;
//                nIds[pos]++;
              }
            }
          }
        }
        else
        {
          dx[pos]=0.0;
          dy[pos]=0.0;
        }
        pos++;
      }
  }

  void preSortData()
  {
//    #pragma omp for schedule(static)//num_threads(omp_get_num_procs())
    for(int n=0;n<N*M;n++)// each pixel
    {
      int N = nIds[n]; // number of independent summands

      if (N<1) continue;

      // insert those not equal to 0 !
      int i = 0;int pos = n*ringsize;

      // swap elements s.t. elements != 0 are first
      for( int it=0; it < ringsize; it++ )
      {
        Scalar diff = u_i[ pos + it ];
        if ( diff != 0.0 && i < N )
          u_i[pos+i++] = diff;
//        events[i++] = diff;
      }
      quick_select((Scalar*) (&u_i[pos]), N);
//    std::sort( (Scalar*) (&output1[pos]), (Scalar*) (&output1[pos+N]) );
    }
  }

  //[sDeltas, sEvents, sElems] = SortAndCumsum_noZeros(deltas, events);
//
//warpParameter.dx       = dx;
//warpParameter.dy       = dy;
//warpParameter.normG    = normG;
//warpParameter.nElems   = sElems;
//
//warpParameter.deltas_  = sDeltas;
//warpParameter.events_  = sEvents;

  /// sort_cumsum_zero
  void preSortDataCensus( const std::vector<Scalar>& events )
  {
    bestSpotsInit = false;
    bestSpot.assign( N*M, 0);
    //     const std::vector<Scalar> events = cmi.getEvents();
    sEvents.assign( N*M*ringsize*2, 0.0 );
    const std::vector<int>       bMask = cmi.getbMask();

    // storing in:
    //  sDeltas,  = deltas
    //  sEvents : new
    //  sElems: nIds

#pragma omp for schedule(static)//num_threads(omp_get_num_procs())
    for(int idv=0;idv<N*M;idv++)//M
    {
      int nEntries = 0;
      int pos = idv*ringsize*2;
//      int breakhere = 0;
//      if( idv == 19 )
//        breakhere = 1;

      if(maskGrad[idv])
        continue;

      std::vector< std::pair<Scalar,Scalar> > mapSort(ringsize*2);

      for(int it=0; it<ringsize*(1-maskGrad[idv]); it++)
      {
        if ( !bMask[ it+idv*ringsize ] )
        {
          if ( events[ it + pos ] != Scalar(0) )
            mapSort[nEntries++] = std::pair<Scalar,Scalar> (deltas[ it + pos ], events[ it + pos ]);
          if ( events[ it + pos + ringsize ] != Scalar(0) )
            mapSort[nEntries++] = std::pair<Scalar,Scalar> (deltas[ it + pos + ringsize ], events[ it + pos + ringsize ]);
        }
      }
      // uses already insertionsort for less than 32 entries:
      std::sort( mapSort.begin(), mapSort.begin()+nEntries );

      int nonZeroEntries = 0;

      // VERIFY THIS (later) or still use above: OK
      // what about this and skip 0 events right away
      for(int it=1;it<nEntries; it++)
      {
        if ((mapSort[nonZeroEntries].first == mapSort[it].first))
        {
          mapSort[nonZeroEntries].second += mapSort[it].second;//mapSort[it].second = 0;
        }
        else
        {
          //          nonZeroEntries++;
          mapSort[++nonZeroEntries] = mapSort[it];
        }
      }
      nonZeroEntries++;

      for(int it=0;it<nonZeroEntries; it++)
      {
        deltas [it + pos] = mapSort[it].first;
        sEvents[it + pos] = mapSort[it].second;
      }
      nIds[idv] = nonZeroEntries;
    }
  }

  /// csad data step
  void dataStep( std::vector<Scalar>& _u, std::vector<Scalar>& _v )
  {
    // first step compute bAdd:
    //    std::vector<Scalar> bAddVec(N*M, 0);

    // the fnuction is sum_i | x-b_i | + 0.5*(x-0)^2
    // input: first is the b's, second is the number of entries in that pixel (e.g. 1 ring: 8 in general, less at iamge border, corners)
    // 3rd input is the weights, assumed to be constant per i

    // this is not the case for zero mean stuff here there is a weight for each (only 2 different actually) i
    //    typedef double Scalar;

//    Scalar *output1, *output2;
    //mwIndex id, idv;

//    std::vector<Scalar> lambdanormG;
//    std::vector<Scalar> nIds;
//    std::vector<Scalar> u_i;

    int vM =  ringsize; // vM components per pixel
    int vN =  M*N; // vN pixel

    int elements=N*M*ringsize;
//    int vM = M*N;

#pragma omp for schedule(static)//num_threads(omp_get_num_procs())
    for(int n=0;n<vN;n++)
    {
      int N_ids = nIds[n]; // number of independent summands
      if (N_ids ==0) continue;
      //    bAdd = (u_-warpParam.u0) .* warpParam.dx + (v_-warpParam.v0) .* warpParam.dy;
      const Scalar bAdd = (_u[n]-u0[n]) * dx[n] + (_v[n]-v0[n]) * dy[n];

      const int pivot = (N_ids-1)/2;
      const Scalar lambda = lambdanormG[n]; // the weight for all of the summands

      //find the i's which bracket the median value
      Scalar x = u_i[ vM*n + pivot ] + bAdd;
      Scalar x_ = (N_ids-x/lambda)/2;
      //    int period = std::max(-1, std::min( int(x_), N+1) );

      int nmid = int((N_ids-1)/2)+1;
      int M_r = N_ids-nmid;
      int M_l = N_ids-M_r-1;

      int k = std::max(-1, std::min( int(x_+1)-1, N_ids) )+1;// int(x_+1)-1 behave like floor for x smaller 0 (greater -1 but after that not important)

      Scalar val(0);

      if (M_l+N_ids-k+1 > k+M_r) // if equal x is median
      {
        std::vector< Scalar > events( 2*M_l+1, 0 );
        int place = 0;
        int add = std::min(k,M_l-k);
        for (int id= add+M_l; id >= k; id--, place++)  // at most k+1 largest elements
//        for (int id= k+M_l; id >= k; id--, place++)  // at most k+1 largest elements
          events[place] = lambda*(N_ids-2*id);      

        for (int id= 0; id < pivot; id++, place++) // M_l elements all smaller than x, if those are super small, than only the k largest of events have to be considered
          events[place] = u_i[vM*n + id] + bAdd;

//        val = quick_select_pick((Scalar*) (&events[0]), events.size(), M_l+k );
        val = quick_select_pick((Scalar*) (&events[0]), place, M_l+add );
        //      std::nth_element( events.begin(), events.begin()+M_l+k, events.end() );// totally dominated by this ?! - superslow
        //      val = events[M_l+k];//% M_l+k are smaller, so one more
      }
      if (M_l+N_ids-k+1 < k+M_r) // if equal x is median
      {
        //% larger: keep M_r
        //% 1+M_l+N-k+1 are discarded, so instead of Nth smallest of 2N+1 we now need
        //% N-(1+M_l+N-k+1) smallest of 2N+1 - (1+M_l+N-k+1)
        //% assuming that the unsorted M_r elements are super large,
        //% we need to keep all elements from entries which can become this
        //% N-(1+M_l+N-k+1)=k-M_l-2 smallest element:
        //% pick the k-M_l-2 smallest of the entries, discard the rest, 
        //% starting at k-1 down to ?:
        //% k - (k-M_l-2 ) - 1 = M_l+1 is the index
        //% restVek = cat(1, lambda*(N-2*(M_l+1:-1:k-1))', vi( M_l+2:N, i) ); %old, fine

        //% should remove the largest of the remainding entries, so go from k-1

        // we start with 2n+1 u_i and 2n lambda*(2n-2i) i= 0..2n
        // now x <0 so the part of elements > x is larger, due to symmetry of lambda*(n-2i)
        // so the median is one of the n elements larger x and 0, .. all lambda*(n-2i) > x
        // assume there are k of these, k>n by definition
        // then we search in the remaining elements the median of 4n+1 points before
        // after discarding n+1 and 2n-k we therefore search for the element no. a with n+1+2n-k + a == 2n+1
        // so a == k-n. but of these k we can directly conclude that only the k-n smallest are of interest
        // the k largest can never become the k-n largest anyway

        //      restVek = cat(1, lambda*(N-2*(k-1:-1:M_l+1))', vi( M_l+2:N, i) );
        std::vector< Scalar > events( 2*M_r+1, 0 );

        int place = 0;
        for (int id= k-1; id > M_l ; id--, place++)
          events[place] = lambda*(N_ids-2*id);

        for (int id= pivot+1; id < N_ids; id++, place++)
          events[place] = u_i[vM*n + id] + bAdd;

        val = quick_select_pick((Scalar*) (&events[0]), place, k-M_l-2 );
        //      std::nth_element( events.begin(), events.begin()+(k-M_l-2), events.begin()+place );// totally dominated by this ?! - superslow
        //      val = events[k-M_l-2];//% M_l+k are smaller, so one more    
      }

      if (M_l+N_ids-k+1 == k+M_r) // if equal x is median
        val = x;

      // delta[n] = val;
      // now push input/whateverin direction of 
      _u[n] -= dx[n]*val;
      _v[n] -= dy[n]*val;
    }
  }


    /// csad data step
  void dataStepSAD( std::vector<Scalar>& _u, std::vector<Scalar>& _v )
  {
    // first step compute bAdd:
    //    std::vector<Scalar> bAddVec(N*M, 0);

    // the fnuction is | x-b_i | + 0.5*(x-0)^2, i=1
    // input: first is the b's, second is the number of entries in that pixel (e.g. 1 ring: 8 in general, less at iamge border, corners)
    // 3rd input is the weights, assumed to be constant per i

    int vN =  M*N; // vN pixel
    int elements=N*M;

#pragma omp for schedule(static)//num_threads(omp_get_num_procs())
    for(int n=0;n<vN;n++)
    {
      if ( mask[n] ) continue;
      //    bAdd = (u_-warpParam.u0) .* warpParam.dx + (v_-warpParam.v0) .* warpParam.dy;
      const Scalar bAdd = (_u[n]-u0[n]) * dx[n] + (_v[n]-v0[n]) * dy[n];

      const Scalar lambdaTauNormG = lambdanormG[n]; // the weight for all of the summands

      //find the i's which bracket the median value
      Scalar G_hat = u_i[ n ] + bAdd;// G^

      Scalar val(0);

      // now standard thresh:
      if ( lambdanormG[n] > G_hat)
      {
        if (G_hat > -lambdanormG[n])
          val = G_hat;
        else // g_hat smallest, so take smaller
        {
//          if (lambdanormG[n] > -lambdanormG[n])// always since its a norm
            val = -lambdanormG[n];
//          else
//            val =  lambdanormG[n];
        }
      }
      else //lambdanormG[n] <= G_hat
      {
        if (G_hat < -lambdanormG[n])
          val = G_hat;
        else  // G_hat largest, so take larger
        {
//          if (lambdanormG[n] > -lambdanormG[n])
            val =  lambdanormG[n];
//          else
//            val = -lambdanormG[n];
        }      
      }
      // delta[n] = val;
      // now push input/whateverin direction of 
      _u[n] -= dx[n]*val;
      _v[n] -= dy[n]*val;
    }
  }

  /// census data step
  void dataStepCensus( std::vector<Scalar>& _u, std::vector<Scalar>& _v, int doStereo )
  {
//  Scalar *input1 = (Scalar*)  mxGetPr(prhs[0]);//deltas
//  Scalar *input2 = (Scalar*)  mxGetPr(prhs[1]);//events
//  Scalar *thetas = (Scalar*)  mxGetPr(prhs[2]);//tau
//  Scalar *nElems = (Scalar*)  mxGetPr(prhs[3]);//nIds

    // well what i could do is:
    // find the minimum cumsum closest to the middle/0 motion
    // then from that spot i just have to run in 1 direction instead of 2
    // namely towards 0 - also knowing the lowest spot
    // i only need to go as far to the other side of 0
    // as the current spot is located !

#pragma omp for schedule(static)//num_threads(omp_get_num_procs())
    for(int idv=0;idv<N*M;idv++)
    {
      const int   n_elem = nIds[idv];
      if (n_elem==0) continue;

      const Scalar bAdd = (_u[idv]-u0[idv]) * dx[idv] + (_v[idv]-v0[idv]) * dy[idv];
      const Scalar theta = tau[idv];
      const int   pos    = idv*2*ringsize;
      // now find smallest index: as  mapSort[it].second + delta^2/theta and  (0, mapSort[it-1].second + delta^2/theta
      // also add the 0 event at delta=0 with effect 0:

      Scalar zerosEventScore(0);
      bool zeroSpot(false); 

      Scalar sum( sEvents[pos] );
      Scalar currentDelta = (deltas[pos]-bAdd)*(deltas[pos]-bAdd) / theta;

      // inconsistency: do i start at o for the eventscore or at the first score in events?
      // case A: start at 0: APPEARS CORRECT
      // case B: start at event[0], then best has no minimum here
      Scalar best = std::min( sEvents[pos] + currentDelta, currentDelta);
      Scalar bestDelta = (deltas[pos]-bAdd);

      if( (deltas[pos]-bAdd) > 0 )
      {
        // case A: this should be 0 then
        zerosEventScore = 0;//events[idv];
        zeroSpot = true;
      }

      // if many inner its:
      // could work with a bound; and currentBest; like best cumsum to current bin
      // then after zeroSpot currentDelta+bound > currentBest -> break
      // if currentDelta will increase only check if 
      // bound == (n_elem-it) * (-evScore) + currentDelta > currentBest -> break
      // if (zeroSpot && ( (n_elem-it) * (-lambda/Scalar(ringsize)) + currentDelta > currentBest) break;) 
      for(int it=1; it<n_elem; it++)
      {
        // really? VERIFY THIS: NO its wrong, since 0's can be assigned see (X)
        if (sEvents[pos+it]==0) continue; // it does not matter

        currentDelta = (deltas[pos+it]-bAdd) * (deltas[pos+it]-bAdd) / theta;

        if (!zeroSpot && (deltas[pos+it]-bAdd) > 0)
        {
          // (X): wrong, should be NOT min(sum, sum+ ..) but only sum unless 0 is an event by itself
//          zerosEventScore = std::min( sum, sum + events[idv+it*step] ); // performs best ?!
          zerosEventScore = sum; // only the score within the bracket
          zeroSpot = true;
        }

        // if one could start at 0 .. or at least from the better side could also use a predictor (min(n..end) (cumsum) )
        if (zeroSpot && ( (n_elem-it) * (-eventScore) + currentDelta > best-sum)) 
          break;

        if (sum + currentDelta < best)
        {
          bestDelta = (deltas[pos+it]-bAdd);
          best = sum + currentDelta;
        }

        sum += sEvents[pos+it];

        if (sum + currentDelta < best)
        {
          bestDelta = (deltas[pos+it]-bAdd);
          best = sum + currentDelta;
        }
      }

      if(!zeroSpot)
      {
        zerosEventScore = sum;
      }
      // now compare with 0 event:

      /////////////////////////////////

      if ( best > zerosEventScore)
      {
        bestDelta = 0;
        best = zerosEventScore; 
      }

      // basically 'set free'
      if (doStereo && _u[idv] + dx[idv]*bestDelta > 0) {bestDelta=0;}

      _u[idv] += dx[idv]*bestDelta;
      _v[idv] += dy[idv]*bestDelta;
    }
  }


  // till a small bug ?!
  void dataStepCensus_init( std::vector<Scalar>& _u, std::vector<Scalar>& _v )
  {
#pragma omp for schedule(static)//num_threads(omp_get_num_procs())
    for(int idv=0;idv<N*M;idv++)
    {
      const int   n_elem = nIds[idv];
      if (n_elem==0) continue;

      const Scalar bAdd = (_u[idv]-u0[idv]) * dx[idv] + (_v[idv]-v0[idv]) * dy[idv];
      const Scalar theta = tau[idv];
      const int   pos    = idv*2*ringsize;
      // now find smallest index: as  mapSort[it].second + delta^2/theta and  (0, mapSort[it-1].second + delta^2/theta
      // also add the 0 event at delta=0 with effect 0:

      Scalar zerosEventScore(0);
      bool zeroSpot(false); 

      Scalar sum( 0 );//sEvents[pos] );
      Scalar currentDelta = 0;//(deltas[pos]-bAdd)*(deltas[pos]-bAdd) / theta;

      /// keep track of
      Scalar bestScore = std::min(0.,sum);
      int startSpot = bestSpot[idv]-1;
      // must look in both directions to get the minimum
      // so spot 

      Scalar best(1000);
      Scalar bestDelta (0);
      Scalar startDelta (deltas[pos+std::max(0,startSpot)]-bAdd);
      // if < 0 go to right as usual

#ifdef _DEBUG
      std::vector< Scalar > ds ( n_elem,0 );
      std::vector< Scalar > es ( n_elem,0 );
      std::copy( deltas.begin()+pos, deltas.begin()+pos+n_elem, ds.begin() );
      std::copy( sEvents.begin()+pos, sEvents.begin()+pos+n_elem, es.begin() );
#endif

      // go right:
      if (startDelta <0)
      {
        Scalar sum(0);
        // best at delta=0 here
        bestDelta = startDelta;
        best = startDelta*startDelta/theta;

//        Scalar spotDelta = deltas[pos+startSpot+1]-bAdd;

        if ( deltas[pos+startSpot+1]-bAdd > 0)// found best: delta=0
        {
          bestDelta = 0;// since left <0 right >0 and best cumsum
          // end search
        }
        else // run along
        {
          zeroSpot = false;
//          best = spotDelta*spotDelta/theta;// best is this since the same cumsum and closer to 0
//          bestDelta = spotDelta;
          // now iterate as usual:
          // add cumsum, compare both delta posistions if smaller
          // remember 0 crossing score
          // break if currentDelta-bAdd > -spotDelta (since cumsum score >=0 by default)

          for(int it=startSpot+1; it<n_elem; it++)
          {
            // really? VERIFY THIS: NO its wrong, since 0's can be assigned see (X)
            if (sEvents[pos+it]==0) continue; // it does not matter

            Scalar delta = deltas[pos+it]-bAdd;
            currentDelta = delta * delta / theta;

            if (!zeroSpot && delta > 0)
            {
              // (X): wrong, should be NOT min(sum, sum+ ..) but only sum unless 0 is an event by itself
              //      zerosEventScore = std::min( sum, sum + events[idv+it*step] ); // performs best ?!
              zerosEventScore = sum; // only the score within the bracket
              zeroSpot = true;
            }

            // if one could start at 0 .. or at least from the better side could also use a predictor (min(n..end) (cumsum) )
            if ( zeroSpot && ( currentDelta > best) )
              break;

            if (sum + currentDelta < best)
            {
              bestDelta = (deltas[pos+it]-bAdd);
              best = sum + currentDelta;
            }

            sum += sEvents[pos+it];

            int notGood=0;
            if (sum<-0.000001)
              notGood =1;

            if (sum + currentDelta < best)
            {
              bestDelta = (deltas[pos+it]-bAdd);
              best = sum + currentDelta;
            }
          }
          // end for
          if(!zeroSpot)
            zerosEventScore = sum;

          // now compare with 0 event:
          /////////////////////////////////
          if ( best > zerosEventScore)
          {
            bestDelta = 0;
            best = zerosEventScore; 
          }
        }
        // new algorithm from here giong to the left
      }
      else // go left as startDelta >= 0
      {
        Scalar sum(0);
        // best at delta=0 here
        bestDelta = startDelta;
        best = startDelta*startDelta/theta;
        // unfortunately we checked out the left part already:   delta checked -> | cumsum==0  |
        // check the same again
        zeroSpot = false;
        for(int it=startSpot; it>=0; it--)
          {
            // really? VERIFY THIS: NO its wrong, since 0's can be assigned see (X)
            if (sEvents[pos+it]==0) continue; // it does not matter

            Scalar delta = deltas[pos+it]-bAdd;
            currentDelta = delta * delta / theta;

            if (!zeroSpot && delta < 0)
            {
              zerosEventScore = sum; // only the score within the bracket
              zeroSpot = true;
            }

            // if one could start at 0 .. or at least from the better side could also use a predictor (min(n..end) (cumsum) )
            if ( zeroSpot && ( currentDelta > best ) )
              break;

            if (sum + currentDelta < best)
            {
              bestDelta = delta;
              best = sum + currentDelta;
            }

            sum -= sEvents[pos+it];//newver < 0 : impossible as the best cumcum was picked so .. if : something is fishy

            int notGood=0;
            if (sum<-0.000001)
              notGood =1;

            if (sum + currentDelta < best)
            {
              bestDelta = delta;
              best = sum + currentDelta;
            }
          }
          // end for
          if(!zeroSpot)
            zerosEventScore = sum;

          // now compare with 0 event:
          /////////////////////////////////
          if ( best > zerosEventScore)
          {
            bestDelta = 0;
            best = zerosEventScore; 
          }    
      }

      _u[idv] += dx[idv]*bestDelta;
      _v[idv] += dy[idv]*bestDelta;
    }
  
  }

  /// census data step with preparation for predictive estimation dataStepCensus_init
  void dataStepCensus_new( std::vector<Scalar>& _u, std::vector<Scalar>& _v )
  {
//  Scalar *input1 = (Scalar*)  mxGetPr(prhs[0]);//deltas
//  Scalar *input2 = (Scalar*)  mxGetPr(prhs[1]);//events
//  Scalar *thetas = (Scalar*)  mxGetPr(prhs[2]);//tau
//  Scalar *nElems = (Scalar*)  mxGetPr(prhs[3]);//nIds

    // well what i could do is:
    // find the minimum cumsum closest to the middle/0 motion
    // then from that spot i just have to run in 1 direction instead of 2
    // namely towards 0 - also knowing the lowest spot
    // i only need to go as far to the other side of 0
    // as the current spot is located !

    if ( bestSpotsInit )
    {
      dataStepCensus_init( _u, _v );// with know spots
      return;
    }

    bestSpotsInit = true;

#pragma omp for schedule(static)//num_threads(omp_get_num_procs())
    for(int idv=0;idv<N*M;idv++)
    {
      const int   n_elem = nIds[idv];
      if (n_elem==0) continue;

      const Scalar bAdd = (_u[idv]-u0[idv]) * dx[idv] + (_v[idv]-v0[idv]) * dy[idv];
      const Scalar theta = tau[idv];
      const int   pos    = idv*2*ringsize;
      // now find smallest index: as  mapSort[it].second + delta^2/theta and  (0, mapSort[it-1].second + delta^2/theta
      // also add the 0 event at delta=0 with effect 0:

      Scalar zerosEventScore(0);
      bool zeroSpot(false); 

      Scalar sum( sEvents[pos] );
      Scalar currentDelta = (deltas[pos]-bAdd)*(deltas[pos]-bAdd) / theta;

      // inconsistency: do i start at o for the eventscore or at the first score in events?
      // case A: start at 0: APPEARS CORRECT
      // case B: start at event[0], then best has no minimum here
      Scalar best = std::min( sEvents[pos] + currentDelta, currentDelta);
      Scalar bestDelta = (deltas[pos]-bAdd);

      /// keep track of
      Scalar bestScore = std::min(0.,sum);
      int bestSpots = bestScore<0 ? 1:0;// spot between delta[bestSpot-1] and delta[bestSpot]
      Scalar bestSpotDelta = currentDelta;

      if( (deltas[pos]-bAdd) > 0 )
      {
        // case A: this should be 0 then
        zerosEventScore = 0;//events[idv];
        zeroSpot = true;
      }

      // if many inner its:
      // could work with a bound; and currentBest; like best cumsum to current bin
      // then after zeroSpot currentDelta+bound > currentBest -> break
      // if currentDelta will increase only check if 
      // bound == (n_elem-it) * (-evScore) + currentDelta > currentBest -> break
      // if (zeroSpot && ( (n_elem-it) * (-lambda/Scalar(ringsize)) + currentDelta > currentBest) break;) 
      for(int it=1; it<n_elem; it++)
      {
        // really? VERIFY THIS: NO its wrong, since 0's can be assigned see (X)
        if (sEvents[pos+it]==0) continue; // it does not matter captured as 0 spot

        currentDelta = (deltas[pos+it]-bAdd) * (deltas[pos+it]-bAdd) / theta;

        if (!zeroSpot && (deltas[pos+it]-bAdd) > 0)
        {
          // (X): wrong, should be NOT min(sum, sum+ ..) but only sum unless 0 is an event by itself
//          zerosEventScore = std::min( sum, sum + events[idv+it*step] ); // performs best ?!
          zerosEventScore = sum; // only the score within the bracket
          zeroSpot = true;
        }

        // try without maybe influence on bestSpot picked ? YES
        // if one could start at 0 .. or at least from the better side could also use a predictor (min(n..end) (cumsum) )
//        if (zeroSpot && ( (n_elem-it) * (-eventScore) + currentDelta > best-sum)) 
//          break;

        if (sum + currentDelta < best)
        {
          bestDelta = (deltas[pos+it]-bAdd);
          best = sum + currentDelta;
        }

        sum += sEvents[pos+it];

        // tie breaker : close to 0/middle OR tie breaker smaller |delta|?
//        if ( bestScore > sum || (bestScore == sum && abs(bestSpots-n_elem/2) > abs(it+1-n_elem/2) ))
        // tiebreaker small |delta| ? :
        if ( bestScore > sum || (bestScore == sum && ( currentDelta < bestSpotDelta ) ) )
        {
          bestSpotDelta = currentDelta;
          bestSpots = it+1;
          bestScore = sum;
        }

        if (sum + currentDelta < best)
        {
          bestDelta = (deltas[pos+it]-bAdd);
          best = sum + currentDelta;
        }
      }

      if(!zeroSpot)
      {
        zerosEventScore = sum;
      }
      // now compare with 0 event:

      /////////////////////////////////

      if ( best > zerosEventScore)
      {
        bestDelta = 0;
        best = zerosEventScore; 
      }

      _u[idv] += dx[idv]*bestDelta;
      _v[idv] += dy[idv]*bestDelta;
      bestSpot[idv] = bestSpots;
    }
  }

  // like this from outside, or computed from inside at start
  void setbImg( Scalar*& _bImg_in, int size )
  {
    b1Img.resize( size );
    std::copy(_bImg_in, &(_bImg_in[size-1]), b1Img.begin());
  }

  void setbImg( std::vector<Scalar>& _bImg_in )
  {
    b1Img.resize( _bImg_in.size() );
    std::copy(_bImg_in.begin(), _bImg_in.end(), b1Img.begin());
  }

//  void setLambda(Scalar _lambda) { lambda = _lambda;}

  std::vector<Scalar>& getLambdaNormG() {return lambdanormG;};
  std::vector<Scalar>& getNIds()        {return nIds;};
  std::vector<Scalar>& get_ui()         {return u_i;};
  std::vector<Scalar>& get_u0()         {return u0;};
  std::vector<Scalar>& get_v0()         {return v0;};
  std::vector<Scalar>& get_dx()         {return dx;};
  std::vector<Scalar>& get_dy()         {return dy;};

    /// csad data step
  void dataStep_ins( std::vector<Scalar>& _u, std::vector<Scalar>& _v )
  {
    // the fnuction is sum_i | x-b_i | + 0.5*(x-0)^2
    // input: first is the b's, second is the number of entries in that pixel (e.g. 1 ring: 8 in general, less at iamge border, corners)
    // 3rd input is the weights, assumed to be constant per i

    // this is not the case for zero mean stuff here there is a weight for each (only 2 different actually) i
    typedef double Scalar;

    int vM =  ringsize; // vM components per pixel
    int vN =  M*N; // vN pixel
    int elements=N*M*ringsize;


#pragma omp for schedule(static)//num_threads(omp_get_num_procs())
    for(int n=0;n<vN;n++)
    {
      int N_ids = nIds[n]; // number of independent summands
      if (N_ids ==0) continue;
      //    bAdd = (u_-warpParam.u0) .* warpParam.dx + (v_-warpParam.v0) .* warpParam.dy;
      const Scalar bAdd = (_u[n]-u0[n]) * dx[n] + (_v[n]-v0[n]) * dy[n];

      int pivot = (N_ids-1)/2;
      Scalar lambda = lambdanormG[n]; // the weight for all of the summands

      //find the i's which bracket the median value
      Scalar x = u_i[ vM*n + pivot ] + bAdd;
      Scalar x_ = (N_ids-x/lambda)/2;
      //    int period = std::max(-1, std::min( int(x_), N+1) );

      int nmid = int((N_ids-1)/2)+1;
      int M_r = N_ids-nmid;
      int M_l = N_ids-M_r-1;

      int k = std::max(-1, std::min( int(x_+1)-1, N_ids) )+1;// int(x_+1)-1 behave like floor for x smaller 0 (greater -1 but after that not important)

      Scalar val(0);

      if (M_l+N_ids-k+1 > k+M_r) // if equal x is median
      {
        std::vector< Scalar > events( 2*M_l+1, 0 );
        int place = 0;

        for (int id= k+M_l; id >= k; id--, place++)
          events[place] = lambda*(N_ids-2*id);      

        for (int id= 0; id < pivot; id++, place++)
          events[place] = u_i[vM*n + id] + bAdd;

//        partial_insertionsort<events.size(), M_l+k, Scalar>(&(events[0])); val = events [M_l+k];
        val = insertionMedian( events, events.size(), M_l+k );
//        val = quick_select_pick((Scalar*) (&events[0]), events.size(), M_l+k );

        //      std::nth_element( events.begin(), events.begin()+M_l+k, events.end() );// totally dominated by this ?! - superslow
        //      val = events[M_l+k];//% M_l+k are smaller, so one more
      }
      if (M_l+N_ids-k+1 < k+M_r) // if equal x is median
      {
        //% larger: keep M_r
        //% 1+M_l+N-k+1 are discarded, so instead of Nth smallest of 2N+1 we now need
        //% N-(1+M_l+N-k+1) smallest of 2N+1 - (1+M_l+N-k+1)
        //% assuming that the unsorted M_r elements are super large,
        //% we need to keep all elements from entries which can become this
        //% N-(1+M_l+N-k+1)=k-M_l-2 smallest element:
        //% pick the k-M_l-2 smallest of the entries, discard the rest, 
        //% starting at k-1 down to ?:
        //% k - (k-M_l-2 ) - 1 = M_l+1 is the index
        //% restVek = cat(1, lambda*(N-2*(M_l+1:-1:k-1))', vi( M_l+2:N, i) ); %old, fine

        //% should remove the largest of the remainding entries, so go from k-1

        //      restVek = cat(1, lambda*(N-2*(k-1:-1:M_l+1))', vi( M_l+2:N, i) );
        std::vector< Scalar > events( 2*M_r+1, 0 );

        int place = 0;
        for (int id= k-1; id > M_l ; id--, place++)
          events[place] = lambda*(N_ids-2*id);

        for (int id= pivot+1; id < N_ids; id++, place++)
          events[place] = u_i[vM*n + id] + bAdd;

        val = insertionMedian( events, place, k-M_l-2 );
//        val = quick_select_pick((Scalar*) (&events[0]), place, k-M_l-2 );

        //      std::nth_element( events.begin(), events.begin()+(k-M_l-2), events.begin()+place );// totally dominated by this ?! - superslow
        //      val = events[k-M_l-2];//% M_l+k are smaller, so one more    
      }

      if (M_l+N_ids-k+1 == k+M_r) // if equal x is median
        val = x;

      // delta[n] = val;
      // now push input/whateverin direction of 
      _u[n] -= dx[n]*val;
      _v[n] -= dy[n]*val;
    }
  }

  Scalar insertionMedian( std::vector<Scalar>& events, int number, int pos)
  {
    switch( number )
    {
     case 1: return events[0];
     case 2:
      switch( pos )
      {
      case 0: partial_insertionsort<2, 0, Scalar>(&(events[0])); return events [0];
      case 1: partial_insertionsort<2, 1, Scalar>(&(events[0])); return events [1];
      }
     case 3: 
      switch( pos )
      {
      case 0: partial_insertionsort<3, 0, Scalar>(&(events[0])); return events [0];
      case 1: partial_insertionsort<3, 1, Scalar>(&(events[0])); return events [1];
      case 2: partial_insertionsort<3, 2, Scalar>(&(events[0])); return events [2];
      }
     case 4: 
      switch( pos )
      {
      case 0: partial_insertionsort<4, 0, Scalar>(&(events[0])); return events [0];
      case 1: partial_insertionsort<4, 1, Scalar>(&(events[0])); return events [1];
      case 2: partial_insertionsort<4, 2, Scalar>(&(events[0])); return events [2];
      case 3: partial_insertionsort<4, 3, Scalar>(&(events[0])); return events [3];
      }
     case 5: 
      switch( pos )
      {
      case 0: partial_insertionsort<5, 0, Scalar>(&(events[0])); return events [0];
      case 1: partial_insertionsort<5, 1, Scalar>(&(events[0])); return events [1];
      case 2: partial_insertionsort<5, 2, Scalar>(&(events[0])); return events [2];
      case 3: partial_insertionsort<5, 3, Scalar>(&(events[0])); return events [3];
      case 4: partial_insertionsort<5, 4, Scalar>(&(events[0])); return events [4];
      }
     case 6:
      switch( pos )
      {
      case 0: partial_insertionsort<6, 0, Scalar>(&(events[0])); return events [0];
      case 1: partial_insertionsort<6, 1, Scalar>(&(events[0])); return events [1];
      case 2: partial_insertionsort<6, 2, Scalar>(&(events[0])); return events [2];
      case 3: partial_insertionsort<6, 3, Scalar>(&(events[0])); return events [3];
      case 4: partial_insertionsort<6, 4, Scalar>(&(events[0])); return events [4];
      case 5: partial_insertionsort<6, 5, Scalar>(&(events[0])); return events [5];
      }
     case 7:
      switch( pos )
      {
      case 0: partial_insertionsort<7, 0, Scalar>(&(events[0])); return events [0];
      case 1: partial_insertionsort<7, 1, Scalar>(&(events[0])); return events [1];
      case 2: partial_insertionsort<7, 2, Scalar>(&(events[0])); return events [2];
      case 3: partial_insertionsort<7, 3, Scalar>(&(events[0])); return events [3];
      case 4: partial_insertionsort<7, 4, Scalar>(&(events[0])); return events [4];
      case 5: partial_insertionsort<7, 5, Scalar>(&(events[0])); return events [5];
      case 6: partial_insertionsort<7, 6, Scalar>(&(events[0])); return events [6];
      }
     case 8:
      switch( pos )
      {
      case 0: partial_insertionsort<8, 0, Scalar>(&(events[0])); return events [0];
      case 1: partial_insertionsort<8, 1, Scalar>(&(events[0])); return events [1];
      case 2: partial_insertionsort<8, 2, Scalar>(&(events[0])); return events [2];
      case 3: partial_insertionsort<8, 3, Scalar>(&(events[0])); return events [3];
      case 4: partial_insertionsort<8, 4, Scalar>(&(events[0])); return events [4];
      case 5: partial_insertionsort<8, 5, Scalar>(&(events[0])); return events [5];
      case 6: partial_insertionsort<8, 6, Scalar>(&(events[0])); return events [6];
      case 7: partial_insertionsort<8, 7, Scalar>(&(events[0])); return events [7];
      }
     case 9:
      switch( pos )
      {
      case 0: partial_insertionsort<9, 0, Scalar>(&(events[0])); return events [0];
      case 1: partial_insertionsort<9, 1, Scalar>(&(events[0])); return events [1];
      case 2: partial_insertionsort<9, 2, Scalar>(&(events[0])); return events [2];
      case 3: partial_insertionsort<9, 3, Scalar>(&(events[0])); return events [3];
      case 4: partial_insertionsort<9, 4, Scalar>(&(events[0])); return events [4];
      case 5: partial_insertionsort<9, 5, Scalar>(&(events[0])); return events [5];
      case 6: partial_insertionsort<9, 6, Scalar>(&(events[0])); return events [6];
      case 7: partial_insertionsort<9, 7, Scalar>(&(events[0])); return events [7];
      case 8: partial_insertionsort<9, 8, Scalar>(&(events[0])); return events [8];
      }
    }
  }
  /*
  // return the median value in a vector of 27 floats pointed to by a
  Scalar heapMedian3( Scalar *a )
  {
   Scalar left[14], right[14], median, *p;
   unsigned char nLeft, nRight;

   // pick first value as median candidate
   p = a;
   median = *p++;
   nLeft = nRight = 1;

   for(;;)
   {
       // get next value
       float val = *p++;

       // if value is smaller than median, append to left heap
       if( val < median )
       {
           // move biggest value to the heap top
           unsigned char child = nLeft++, parent = (child - 1) / 2;
           while( parent && val > left[parent] )
           {
               left[child] = left[parent];
               child = parent;
               parent = (parent - 1) / 2;
           }
           left[child] = val;

           // if left heap is full
           if( nLeft == 14 )
           {
               // for each remaining value
               for( unsigned char nVal = 27 - (p - a); nVal; --nVal )
               {
                   // get next value
                   val = *p++;

                   // if value is to be inserted in the left heap
                   if( val < median )
                   {
                       child = left[2] > left[1] ? 2 : 1;
                       if( val >= left[child] )
                           median = val;
                       else
                       {
                           median = left[child];
                           parent = child;
                           child = parent*2 + 1;
                           while( child < 14 )
                           {
                               if( child < 13 && left[child+1] > left[child] )
                                   ++child;
                               if( val >= left[child] )
                                   break;
                               left[parent] = left[child];
                               parent = child;
                               child = parent*2 + 1;
                           }
                           left[parent] = val;
                       }
                   }
               }
               return median;
           }
       }

       // else append to right heap
       else
       {
           // move smallest value to the heap top
           unsigned char child = nRight++, parent = (child - 1) / 2;
           while( parent && val < right[parent] )
           {
               right[child] = right[parent];
               child = parent;
               parent = (parent - 1) / 2;
           }
           right[child] = val;

           // if right heap is full
           if( nRight == 14 )
           {
               // for each remaining value
               for( unsigned char nVal = 27 - (p - a); nVal; --nVal )
               {
                   // get next value
                   val = *p++;

                   // if value is to be inserted in the right heap
                   if( val > median )
                   {
                       child = right[2] < right[1] ? 2 : 1;
                       if( val <= right[child] )
                           median = val;
                       else
                       {
                           median = right[child];
                           parent = child;
                           child = parent*2 + 1;
                           while( child < 14 )
                           {
                               if( child < 13 && right[child+1] < right[child] )
                                   ++child;
                               if( val <= right[child] )
                                   break;
                               right[parent] = right[child];
                               parent = child;
                               child = parent*2 + 1;
                           }
                           right[parent] = val;
                       }
                   }
               }
               return median;
           }
       }
   }
}
*/
  ////////
  private:

  int N;
  int M;

  Scalar* Img2;
  Scalar* tau;

  std::vector<Scalar> dx;
  std::vector<Scalar> dy;

  std::vector<Scalar> u0;
  std::vector<Scalar> v0;

  /// image coordinates matlab style, use in interpoaltion, mask for bimg computation
  std::vector<Scalar> Iu;
  std::vector<Scalar> Iv;
  /// mask defining good points
  std::vector<int>    mask;
  /// mask defining bad gradients: one more than good points!
  std::vector<int>    maskGrad;

  std::vector<Scalar> lambdanormG;
  std::vector<Scalar> nIds;
  std::vector<Scalar> u_i;

  std::vector<Scalar> b1Img;

  std::vector<Scalar> deltas;
  std::vector<Scalar> sEvents;

  int    ringsize;
  int    ring;
  Scalar lambda;
  Scalar cEps;
  Scalar eventScore;

  BiCubicIp<Scalar> bcIP;
  // also possible external?
  censusMaskImaging<Scalar> cmi;

  bool bestSpotsInit;
  std::vector<int> bestSpot;
};

#endif // #define _CSAD_PAR_H