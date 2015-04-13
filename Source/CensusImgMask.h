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

#ifndef _CENSUS_MASK_IMG_H
#define _CENSUS_MASK_IMG_H

#include <vector>
#include <algorithm>
#include <cmath>

template <typename S_>
class censusMaskImaging
{
public:

  typedef S_ Scalar;

  censusMaskImaging(int _M, int _N, int _ring, Scalar _lambda = 12.33, Scalar _cEps = 0) 
    :N(_N), M(_M), dx(_M), dy(_N), ring(_ring),
    ringsize ((2*_ring+1)*(2*_ring+1)-1), pixel(N*M), cEps(_cEps), lambda(_lambda)
  {generateMem();eventScore = lambda/ringsize;};

  ~censusMaskImaging() { freeMem(); };

  void freeMem(){};

  void generateMem()
  {
    int pos = 0;
    for(int j=-ring;j<=ring;j++)
      for(int i=-ring;i<=ring;i++)
      {
        if (i==0 && j==0)
          continue;
        if((j<i) || (j<=i && i>0) ) //so j>0
        {
          displacementsY.push_back(j);
          displacementsX.push_back(i);
          displacementsP.push_back(pos);
        }
        pos++;
      }
  }

  void compute(Scalar *img, Scalar *bimg, int* mask=NULL, int* bmask = NULL)
  {
#pragma omp for schedule(static)//num_threads(omp_get_num_procs())
    for( int n=0; n<pixel; n++ )//N*M
    {
      int idv  = n;
      int px   = n%dx;// next column
      int py   = n/dx;// next row
      int pos0 = n*ringsize; 

      if( mask && mask[n] ) // all invalid, bimg: not defined
        continue;

      // pixel not masked
      for (int i=0;i<displacementsX.size();i++)
      {
        int inpos  = displacementsP[i];
        int inposi = ringsize-1 - displacementsP[i];
        int ppx = px+displacementsX[i];
        int ppy = py+displacementsY[i];
        // other position:
        int pos1 = ppx + ppy*dx;
        if( ppx>=0 && ppx<dx && ppy>=0 && ppy<dy && !(mask && mask[pos1]) )
        {
          Scalar val1 = img[n];
          Scalar val2 = img[pos1];

          bimg[pos0 +inpos] = val1-val2;
          bmask[pos0 +inpos] = 0;
          // neighbor as well(if exists):
          bimg[pos1*ringsize + inposi] = val2-val1;
          bmask[pos1*ringsize + inposi] = 0;
        }
      }
    }
  }

  void compute(std::vector<Scalar>& img, std::vector<int>&  mask)
  {
    bimg.resize(pixel*ringsize, 0);
    bimg.assign(pixel*ringsize, 0.0);
    bmask.resize(pixel*ringsize, 1);
    bmask.assign(pixel*ringsize, 1);

    int boolMask = (mask.size() == pixel);

#pragma omp for schedule(static)//num_threads(omp_get_num_procs())
    for( int n=0; n<pixel; n++ )//N*M
    {
      int idv  = n;
      int px   = n%dx;// next column
      int py   = n/dx;// next row
      int pos0 = n*ringsize; 

      if( boolMask && mask[n] ) // all invalid, bimg: not defined
        continue;

      // pixel not masked
      for (int i=0;i<displacementsX.size();i++)
      {
        int inpos  = displacementsP[i];
        int inposi = ringsize-1 - displacementsP[i];
        int ppx = px+displacementsX[i];
        int ppy = py+displacementsY[i];
        // other position:
        int pos1 = ppx + ppy*dx;
        if( ppx>=0 && ppx<dx && ppy>=0 && ppy<dy && !(boolMask && mask[pos1]) )
        {
          Scalar val1 = img[n];
          Scalar val2 = img[pos1];

          bimg[pos0 +inpos] = val1-val2;
          bmask[pos0 +inpos] = 0;
          // neighbor as well(if exists):
          bimg[pos1*ringsize + inposi] = val2-val1;
          bmask[pos1*ringsize + inposi] = 0;
        }
      }
    }
  }

  void compute(Scalar*& img)//, std::vector<int>&  mask)
  {
    bimg.resize(pixel*ringsize, 0);
    bimg.assign(pixel*ringsize, 0.0);
    bmask.resize(pixel*ringsize, 1);
    bmask.assign(pixel*ringsize, 1);

//    int boolMask = (mask.size() == pixel);

#pragma omp for schedule(static)//num_threads(omp_get_num_procs())
    for( int n=0; n<pixel; n++ )//N*M
    {
      int idv  = n;
      int px   = n%dx;// next column
      int py   = n/dx;// next row
      int pos0 = n*ringsize; 

//      if( boolMask && mask[n] ) // all invalid, bimg: not defined
//        continue;

      // pixel not masked
      for (int i=0;i<displacementsX.size();i++)
      {
        int inpos  = displacementsP[i];
        int inposi = ringsize-1 - displacementsP[i];
        int ppx = px+displacementsX[i];
        int ppy = py+displacementsY[i];
        // other position:
        int pos1 = ppx + ppy*dx;
        if(ppx>=0 && ppx<dx && ppy>=0 && ppy<dy )//&& !(boolMask && mask[pos1]) )
        {
          Scalar val1 = img[n];
          Scalar val2 = img[pos1];

          bimg[pos0 +inpos] = val1-val2;
          bmask[pos0 +inpos] = 0;
          // neighbor as well(if exists):
          bimg[pos1*ringsize + inposi] = val2-val1;
          bmask[pos1*ringsize + inposi] = 0;
        }
      }
    }
  }

//eventP =  (lambda/cells) * ones( M, N, (2*(radius+1)-1)^2-1 );
//eventM = (-lambda/cells) * ones( M, N, (2*(radius+1)-1)^2-1 );
//
//if numel(cEps)==1                          %   eM   eP
//  eventP(fullImg >  cEps) = -lambda/cells; % ~~~~~~~~|_____
//  eventP(fullImg < -cEps) = 0;             % ___|~~~~~~~~~~
//
//  eventM(fullImg >  cEps) = 0;             % ~~~~~~~~|_____
//  eventM(fullImg < -cEps) = lambda/cells;  % ___|~~~~~~~~~~
//else  
//  eventP(bsxfun(@gt ,fullImg,  cEps)) = -lambda/cells;
//  eventP(bsxfun(@lt ,fullImg, -cEps)) = 0;
//  
//  eventM(bsxfun(@gt ,fullImg,  cEps)) = 0;
//  eventM(bsxfun(@lt ,fullImg, -cEps)) = lambda/cells;
//end
//
//% those havbe no partner in the image and therefore are not to be looked at
//eventM(fullMask) = 0;
//eventP(fullMask) = 0;
//events           = cat ( 3, eventP, eventM );

  void   setEps(Scalar cEps_) {cEps = cEps_;};
  Scalar getEps( )            {return cEps;};
  void setLambda(Scalar lambda_) {lambda = lambda_;};

  void computeEvents( Scalar*& img )//, std::vector<int>&  mask)
  {
    events.resize(pixel*ringsize*2, 0);
    events.assign(pixel*ringsize*2, 0.0);
    bmask.resize(pixel*ringsize, 1);
    bmask.assign(pixel*ringsize, 1);
//    bimg.resize(pixel*ringsize, 0);
//    bimg.assign(pixel*ringsize, 0.0);

    Scalar value = lambda/Scalar(ringsize);
    eventScore   = value;

#pragma omp for schedule(static)//num_threads(omp_get_num_procs())
    for( int n=0; n<pixel; n++ )//N*M
    {
      int idv  = n;
      int px   = n%dx;// next column
      int py   = n/dx;// next row
      int pos0 = n*ringsize*2;
      int pos0b= n*ringsize;

      // pixel not masked
      for (int i=0;i<displacementsX.size();i++)
      {
        int inpos  = displacementsP[i];
        int inposi = ringsize-1 - displacementsP[i];
        int ppx = px+displacementsX[i];
        int ppy = py+displacementsY[i];
        // other position:
        int pos1 = ppx + ppy*dx;
        if(ppx>=0 && ppx<dx && ppy>=0 && ppy<dy )//&& !(boolMask && mask[pos1]) )
        {
          Scalar val1 = img[n];
          Scalar val2 = img[pos1];

          //  eventP(bsxfun(@gt ,fullImg,  cEps)) = -lambda/cells;
          //  eventP(bsxfun(@lt ,fullImg, -cEps)) = 0;
          //  eventM(bsxfun(@gt ,fullImg,  cEps)) = 0;
          //  eventM(bsxfun(@lt ,fullImg, -cEps)) = lambda/cells;
#define _original_version_
#ifdef _original_version_
          if (val1-val2 > cEps)
          {
            events[pos0            + inpos            ] = -value;
            events[pos0            + inpos + ringsize ] = Scalar(0);
            // if val1-val2 > cEps then val2-val1 < -cEps
            events[pos1*ringsize*2 + inposi           ] = Scalar(0);
            events[pos1*ringsize*2 + inposi + ringsize] = value;
          }
          else if (val1-val2 < -cEps)
          {
            events[pos0            + inpos            ] = Scalar(0);
            events[pos0            + inpos + ringsize ] = value;
            // neighbor as well(if exists):
            events[pos1*ringsize*2 + inposi           ] = -value;
            events[pos1*ringsize*2 + inposi + ringsize] = Scalar(0);
          }
          else 
          {
            events[pos0            + inpos            ] =  value;
            events[pos0            + inpos + ringsize ] = -value;// +
            // neighbor as well(if exists):
            events[pos1*ringsize*2 + inposi           ] =  value;// -
            events[pos1*ringsize*2 + inposi + ringsize] = -value;
          }

#else
          Scalar localDiff      = fabs(val1-val2) - cEps;
//          Scalar localWeight    = Scalar(2.)-std::pow(2., Scalar(0.6-localDiff*32.));
          Scalar localWeight    = Scalar(1.5)-pow(2., Scalar(-localDiff*128.));// super sharp
          if (val1-val2 > cEps)
          {
            events[pos0            + inpos            ] = -value*localWeight;
            events[pos0            + inpos + ringsize ] = Scalar(0);
            // if val1-val2 > cEps then val2-val1 < -cEps
            events[pos1*ringsize*2 + inposi           ] = Scalar(0);
            events[pos1*ringsize*2 + inposi + ringsize] = value*localWeight;
          }
          else if (val1-val2 < -cEps)
          {
            events[pos0            + inpos            ] = Scalar(0);
            events[pos0            + inpos + ringsize ] = value*localWeight;
            // neighbor as well(if exists):
            events[pos1*ringsize*2 + inposi           ] = -value*localWeight;
            events[pos1*ringsize*2 + inposi + ringsize] = Scalar(0);
          }
          else
          {
            events[pos0            + inpos            ] =  value*localWeight;
            events[pos0            + inpos + ringsize ] = -value*localWeight;// +
            // neighbor as well(if exists):
            events[pos1*ringsize*2 + inposi           ] =  value*localWeight;// -
            events[pos1*ringsize*2 + inposi + ringsize] = -value*localWeight;
          }

#endif

          bmask[pos0b+inpos] = 0;
          // neighbor as well(if exists):
          bmask[pos1*ringsize + inposi] = 0;
        }
      }
    }
  }


  std::vector<Scalar>& getbImg() {return bimg;};

  std::vector<int>&    getbMask() {return bmask;}

  std::vector<Scalar>& getEvents() {return events;};

private:

  std::vector<Scalar> bimg;
  std::vector<int>    bmask;

  std::vector<int> displacementsX;
  std::vector<int> displacementsY;
  std::vector<int> displacementsP; // storage position
  
  int N;
  int M;
  int dx;
  int dy;

  int ring;
  int pixel;
  int ringsize;

  // census stuff
  Scalar cEps;
  Scalar lambda;
  Scalar eventScore;
  std::vector<Scalar> events;
};

#endif // _CENSUS_MASK_IMG_H