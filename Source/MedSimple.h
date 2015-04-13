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

#ifndef _MED_SIMPLE_H
#define _MED_SIMPLE_H

#include <vector>
#include <algorithm>
#include "quickSelect.cpp"

template <typename S_>
class peakFilter
{
public:

  typedef S_ Scalar;

  peakFilter( int _M, int _N, int _ring ) 
    :N(_N), M(_M), ring(_ring), ringW(2*_ring+1)
  { generateMem(); };

  ~peakFilter() {};

  void freeMem(){};
  void generateMem(){};

  void compute_noRepeat( std::vector<Scalar>& img )
  {
    std::vector<Scalar> medImg(img.size(), 0);
    std::vector<Scalar> diff  (img.size(), 0);

//    printf("M : %d, N :%d\n", M, N);

    //u_ = medfilt2(u,[3 3],'symmetric');
    //diff = abs(u-u_);
    //v = mean(abs(u_(:)));
    Scalar meanU(0);
    int medianSize = (2*ring+1)*(2*ring+1);
    Scalar nElements(0);
//#pragma omp for schedule(static)//num_threads(omp_get_num_procs())
    for( int m=ring; m<M-ring; m++ )
    {
//      if (M-ringW < m) continue;
      // init roi
      std::vector<Scalar> mfs( medianSize, Scalar(0) );
      std::vector<Scalar> medField( medianSize, Scalar(0) );

      for (int rr = 0; rr < ringW; rr++)
        std::copy( img.begin()+M*rr+m-ring, img.begin()+M*rr+m+ring+1, mfs.begin()+rr*ringW );

      // to be replaced first , then (rr+1) % ringW
      int rr = ringW-1;
      for( int n=ring; n<N-ring; n++ )
      {
        int idv  = n*M+m; // middle position
        // first replace old entry
        std::copy( img.begin()+M*(n+ring)+m-ring, img.begin()+M*(n+ring)+m+ring+1, mfs.begin()+rr*ringW );
        rr = (rr+1) % ringW;

        std::copy(mfs.begin(), mfs.end(), medField.begin());
        Scalar median = quick_select( (Scalar*) (&(medField[0])), medianSize );
        medImg[idv] = median;
        diff[idv]   = fabs( median-img[idv] );
//        Scalar tet1 = fabs(median)/(nElements+Scalar(1));
//        Scalar test2 = meanU * nElements / (nElements+Scalar(1));

        meanU = fabs(median)/(nElements+Scalar(1)) + meanU * nElements / (nElements+Scalar(1));
        nElements+=Scalar(1);
      }
    }

    // u(diff > v) = u_(diff > v);
    for( int m=ring; m<M-ring; m++ )
      for( int n=ring; n<N-ring; n++ )
      {
        int idv  = n*M+m;
        if( diff[idv] > meanU )
          img[idv] = medImg[idv];
      }
  }

  void compute( std::vector<Scalar>& img, bool force = false )
  {
    std::vector<Scalar> medImg(img.size(), 0);
    std::vector<Scalar> diff  (img.size(), 0);

//    printf("M : %d, N :%d\n", M, N);

    //u_ = medfilt2(u,[3 3],'symmetric');
    //diff = abs(u-u_);
    //v = mean(abs(u_(:)));
    Scalar meanU(0);
    int medianSize = (2*ring+1)*(2*ring+1);
    Scalar nElements(0);
//#pragma omp for schedule(static)//num_threads(omp_get_num_procs())
    for( int m=0; m<M; m++ )
    {
//      if (M-ringW < m) continue;
      // init roi
      std::vector<Scalar> mfs( medianSize, Scalar(0) );
      std::vector<Scalar> medField( medianSize, Scalar(0) );
      int n=0;

      if( m>=ring && m < M-ring )
      {
      for (int rr = -ring; rr <= ring; rr++)
        std::copy( img.begin()+M*std::max(rr,0)+m-ring, img.begin()+M*std::max(rr,0)+m+ring+1, mfs.begin()+(rr+ring)*ringW );
      }
      else // have to repeat borders -but also everywhere else :(
      {
        int pos =0;
        for (int i=-ring;i<=ring;i++)
          for (int j=-ring;j<=ring;j++)
            mfs[pos++] = img[ std::min(N-1, std::max(n+i, 0))*M + std::min( M-1, std::max(0, m+j ) ) ];
      }

      // to be replaced first , then (rr+1) % ringW
      int rr = ringW-1;
      for( ; n<N; n++ )
      {
        int idv  = n*M+m; // middle position
        // first replace old entry
        if( m>=ring && m < M-ring )
           std::copy( img.begin()+M*std::min(N-1,n+ring)+m-ring, img.begin()+M*std::min(N-1,n+ring)+m+ring+1, mfs.begin()+rr*ringW );
        else
        {
          int pos =0;
          for (int i=-ring;i<=ring;i++)
            for (int j=-ring;j<=ring;j++)
              mfs[pos++] = img[ std::min(N-1, std::max(n+i, 0))*M + std::min( M-1, std::max(0, m+j ) ) ];
        }

        rr = (rr+1) % ringW;

        std::copy(mfs.begin(), mfs.end(), medField.begin());
        Scalar median = quick_select( (Scalar*) (&(medField[0])), medianSize );
        medImg[idv] = median;
        diff[idv]   = fabs( median-img[idv] );

        meanU = fabs(median)/(nElements+Scalar(1)) + meanU * nElements / (nElements+Scalar(1));
        nElements+=Scalar(1);
      }
    }

    // u(diff > v) = u_(diff > v);
    for( int m=0; m<M; m++ )
      for( int n=0; n<N; n++ )
      {
        int idv  = n*M+m;
        if( diff[idv] > meanU || force)
          img[idv] = medImg[idv];
      }
  }

private:
  
  int N;
  int M;

  int ringW;
  int ring;
};

#endif // _MED_SIMPLE_H