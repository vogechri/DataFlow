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

#ifndef _TGV_SPARSE_MATRIX_H
#define _TGV_SPARSE_MATRIX_H

#include <vector>
#include <algorithm>

#ifndef _noeigen_
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#else
#include "OwnSparseMatrix.h"
#endif

//nr = full(sum(abs(A),2));
//nc = full(sum(abs(A),1))';
//normalizerR = 1./max(1, nr);
//normalizerC = 1./max(1, nc);
//
//tt     = spdiags(normalizerR , 0, size(A,1), size(A,1));
//tau    = reshape(normalizerC(1:N*M) , M,N);
//ss     = spdiags(normalizerC, 0, size(A,2), size(A,2));
//
//ttA = tt*A;
//ssA = ss*A';

//    y = y + ttA*(vv_);
//    y = reprojectY(y);
//    
//    vv_=vv;
//    vv=vv-ssA*y;
//
//    % extract u,v:
//    u = reshape(vv(1:N*M), M,N);
//    v = reshape(vv(1+N*M:2*N*M),M,N);    
//
//  if par.dataTerm == 1    
//    [u, v] = solveCsadData(  u,v, warpParam );
//  end
//  if par.dataTerm == 0
//    %    [u, v] = solveNCCInner( u,v, lambda.*tau, false, warpParam );
//    [u, v] = solveCensusData( warpParam, u,v, tau );
//  end
//  
//  vv(1:N*M)       = u(:);
//  vv(1+N*M:2*N*M) = v(:);
//    
//  vv_ = 2*vv-vv_;

// usage: step1();  
//        extract_uv(std::vector<Scalar>& vv, std::vector<Scalar>& uu)
//        data step
//        set_uv(std::vector<Scalar> vv, std::vector<Scalar>& uu)
//        store2_uvw()
// peakfilter
template <typename S_>
class ownMatrixVecor
{
public:

  typedef S_ Scalar;
#ifndef _noeigen_
  typedef Eigen::Triplet<Scalar> T;
#endif

  ownMatrixVecor( int _rows, int _cols, Scalar*& _uvw, Scalar*& _pq, int _elements, Scalar*& ids_i, Scalar*& ids_j, Scalar*& ids_v, 
    int doTV=0, Scalar* _cons = NULL, int _nCons =0)
    : rows(_rows), cols(_cols), mat( _rows, _cols ), matT( _cols, _rows ), uvw( _cols ), uvw_( _cols ), pq( _rows ), elements(_elements), cons(_cons), nCons(_nCons)
  {
    if (!doTV)
      pixel = ((_rows-nCons)/12);
    else
      pixel = ((_rows-nCons)/4);
    // with constraints, n_rows is higher, also pq do not fit well and pixel also not correct
    // first idea: put into pq on rhs, later restore all but the rest
    // pixel must be adjusted
    createVectors( _uvw, _pq );
    createMatrices( ids_i, ids_j, ids_v );
  };

  ~ownMatrixVecor(){};

  void extract_uv(std::vector<Scalar>& uu, std::vector<Scalar>& vv)
  {
    vv.resize(pixel, Scalar(0.));
    uu.resize(pixel, Scalar(0.));

    for (int i=0;i<pixel;i++)
    {
      uu[i] = uvw[i];
      vv[i] = uvw[i+pixel];
    }
  }

  void set_uv(std::vector<Scalar>& uu, std::vector<Scalar>& vv)
  {
    for (int i=0;i<pixel;i++)
    {
      uvw[i] = uu[i];
      uvw[i+pixel] = vv[i];
    }
  }

  /// choices can be total variation or total generalized variation (order 1) == default
  void step1(int doTV = 0)
  {
    // ok: store2_uvw called from outside after set_uv
#ifndef _noeigen_
    pq = pq + mat * uvw_;//y+sigma*D* (2*v^k+1 - v^k)
#else
    mat.multplus ( uvw_, pq );
    Scalar test = pq[2];
#endif
    if (!doTV)
      reproject_PQ_2_4();
    else
      reproject_PQ_TV();

    store_uvw(); // saves v^k for a moment (needed, see above)
#ifndef _noeigen_
    uvw = uvw - matT * pq;// v^k-tau D^T y^k -- proxmap of G from outside
#else
    mat.multTminus ( pq, uvw );
#endif
  }

  void store_uvw() {uvw_=uvw;}

#ifndef _noeigen_
  void store2_uvw() {uvw_=2*uvw-uvw_;}
#else
  void store2_uvw() {uvw_.addTwiceSub(uvw);}
#endif

  void createVectors(Scalar*& _uvw, Scalar*& _pq)
  {
    for (int i=0; i< uvw.size(); i++)
      uvw[i] = _uvw[i];
    // add extra constraints here ! 
    for (int i=0; i< pq.size()-nCons; i++)
      pq[i] = _pq[i];

    for (int i=pq.size()-nCons;i<pq.size(); i++)
      pq[i] = 0;

//    int restCons = pq.size()-nCons;
//    for (int i=0; i < nCons; i++)
//      pq[i+restCons] = cons[i];

    store_uvw();
  }

  void createMatrices(Scalar*& ids_i, Scalar*& ids_j, Scalar*& ids_v)
  {
    // first compute normalizers
    std::vector<Scalar> normalizerR(rows,0);
    std::vector<Scalar> normalizerC(cols,0);

    for(int i=0; i < elements; i++)
    {
      normalizerR[ ids_i[i]-1 ] += fabs(ids_v[i]);
      normalizerC[ ids_j[i]-1 ] += fabs(ids_v[i]);
    }

    tau.resize( pixel );
    for (int i =0; i< pixel ; i++)
      if (normalizerC[i] > 0.00001)
        tau[i] = Scalar(1) / normalizerC[i];
      else
        tau[i] = Scalar(0);

#ifndef _noeigen_
    std::vector<T> tripletList;
    tripletList.reserve( elements );
    std::vector<T> tripletListT;
    tripletListT.reserve( elements );
    for(int i=0; i < elements; i++)
    {
      //int col_    = ids_i[i];
      //int row_    = ids_j[i];
      //Scalar val_ = ids_v[i];

      if ( normalizerR[ ids_i[i]-1 ] > 0 )
        tripletList.push_back ( T( ids_i[i]-1, ids_j[i]-1, ids_v[i]/normalizerR[ ids_i[i]-1 ] ) );

      if ( normalizerC[ ids_j[i]-1 ] > 0 )
        tripletListT.push_back( T( ids_j[i]-1, ids_i[i]-1, ids_v[i]/normalizerC[ ids_j[i]-1 ] ) );
    }

    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    matT.setFromTriplets(tripletListT.begin(), tripletListT.end());
#else

    mat.setValues ( elements, ids_i, ids_j, ids_v, normalizerR, normalizerC );

#endif



    if (nCons > 0)
      cons_rhs.resize(nCons);

    for (int i=0;i<nCons;i++)
      if (normalizerR[rows-nCons+i]>0)
       cons_rhs[i] = cons[i] / normalizerR[rows-nCons+i];

  }

  void reproject_PQ()
  {
//    int pixel = rows/12;
    for (int i =0;i < pixel; i++) // x dir, here jump by vM
    {
      Scalar p0  = pq[ i           ];
      Scalar p1  = pq[ i+    pixel ];
      Scalar p2  = pq[ i+  2*pixel ];
      Scalar p3  = pq[ i+  3*pixel ];
      Scalar p4  = pq[ i+  4*pixel ];
      Scalar p5  = pq[ i+  5*pixel ];
      Scalar p6  = pq[ i+  6*pixel ];
      Scalar p7  = pq[ i+  7*pixel ];
      Scalar p8  = pq[ i+  8*pixel ];
      Scalar p9  = pq[ i+  9*pixel ];
      Scalar p10 = pq[ i+ 10*pixel ];
      Scalar p11 = pq[ i+ 11*pixel ];

      Scalar rp1 = Scalar(1) / std::max( Scalar(1), sqrt( p1*p1 + p2*p2 + p3*p3 + p0*p0 ) );
      Scalar rp2 = Scalar(1) / std::max( Scalar(1), sqrt( p5*p5 + p6*p6 + p7*p7 + p4*p4 ) );
      Scalar rp3 = Scalar(1) / std::max( Scalar(1), sqrt( p9*p9 + p10*p10 + p11*p11 + p8*p8 ) );

      // out:
      pq[ i           ] = p0 * rp1;
      pq[ i+    pixel ] = p1 * rp1;
      pq[ i+  2*pixel ] = p2 * rp1;
      pq[ i+  3*pixel ] = p3 * rp1;
      pq[ i+  4*pixel ] = p4 * rp2;
      pq[ i+  5*pixel ] = p5 * rp2;
      pq[ i+  6*pixel ] = p6 * rp2;
      pq[ i+  7*pixel ] = p7 * rp2;
      pq[ i+  8*pixel ] = p8 * rp3;
      pq[ i+  9*pixel ] = p9 * rp3;
      pq[ i+ 10*pixel ] = p10 * rp3;
      pq[ i+ 11*pixel ] = p11 * rp3;
    }
    // if constraints reproject these here: use xtra function since called multiple times
    reproject_Constraint( 12*pixel );
  }


  void reproject_PQ_TV()
  {
    for (int i =0;i < pixel; i++) // x dir, here jump by vM
    {
      Scalar p0  = pq[ i           ];
      Scalar p1  = pq[ i+    pixel ];
      Scalar p2  = pq[ i+  2*pixel ];
      Scalar p3  = pq[ i+  3*pixel ];

      Scalar rp0 = Scalar(1) / std::max( Scalar(1), sqrt( p1*p1 + p0*p0 + p2*p2 + p3*p3 ) );

      // out:
      pq[ i           ] = p0 * rp0;
      pq[ i+    pixel ] = p1 * rp0;
      pq[ i+  2*pixel ] = p2 * rp0;
      pq[ i+  3*pixel ] = p3 * rp0;
    }
    // if constraints reproject here
    reproject_Constraint( 4*pixel );
  }

  void reproject_PQ_2_4()
  {
//    int pixel = rows/12;
    for (int i =0;i < pixel; i++) // x dir, here jump by vM
    {
      Scalar p0  = pq[ i           ];
      Scalar p1  = pq[ i+    pixel ];
      Scalar p2  = pq[ i+  2*pixel ];
      Scalar p3  = pq[ i+  3*pixel ];
      Scalar p4  = pq[ i+  4*pixel ];
      Scalar p5  = pq[ i+  5*pixel ];
      Scalar p6  = pq[ i+  6*pixel ];
      Scalar p7  = pq[ i+  7*pixel ];
      Scalar p8  = pq[ i+  8*pixel ];
      Scalar p9  = pq[ i+  9*pixel ];
      Scalar p10 = pq[ i+ 10*pixel ];
      Scalar p11 = pq[ i+ 11*pixel ];

      Scalar rp0 = Scalar(1) / std::max( Scalar(1), sqrt( p1*p1 + p0*p0 ) );
      Scalar rp1 = Scalar(1) / std::max( Scalar(1), sqrt( p2*p2 + p3*p3 ) );
      Scalar rp2 = Scalar(1) / std::max( Scalar(1), sqrt( p5*p5 + p6*p6 + p7*p7 + p4*p4 ) );
      Scalar rp3 = Scalar(1) / std::max( Scalar(1), sqrt( p9*p9 + p10*p10 + p11*p11 + p8*p8 ) );

      // out:
      pq[ i           ] = p0 * rp0;
      pq[ i+    pixel ] = p1 * rp0;
      pq[ i+  2*pixel ] = p2 * rp1;
      pq[ i+  3*pixel ] = p3 * rp1;
      pq[ i+  4*pixel ] = p4 * rp2;
      pq[ i+  5*pixel ] = p5 * rp2;
      pq[ i+  6*pixel ] = p6 * rp2;
      pq[ i+  7*pixel ] = p7 * rp2;
      pq[ i+  8*pixel ] = p8 * rp3;
      pq[ i+  9*pixel ] = p9 * rp3;
      pq[ i+ 10*pixel ] = p10 * rp3;
      pq[ i+ 11*pixel ] = p11 * rp3;
    }
    reproject_Constraint( 12*pixel );
  }

  void reproject_PQ_4_8()
  {
//    int pixel = rows/12;
    for (int i =0;i < pixel; i++) // x dir, here jump by vM
    {
      Scalar p0  = pq[ i           ];
      Scalar p1  = pq[ i+    pixel ];
      Scalar p2  = pq[ i+  2*pixel ];
      Scalar p3  = pq[ i+  3*pixel ];
      Scalar p4  = pq[ i+  4*pixel ];
      Scalar p5  = pq[ i+  5*pixel ];
      Scalar p6  = pq[ i+  6*pixel ];
      Scalar p7  = pq[ i+  7*pixel ];
      Scalar p8  = pq[ i+  8*pixel ];
      Scalar p9  = pq[ i+  9*pixel ];
      Scalar p10 = pq[ i+ 10*pixel ];
      Scalar p11 = pq[ i+ 11*pixel ];

      Scalar rp1 = Scalar(1) / std::max( Scalar(1), sqrt( p1*p1 + p2*p2 + p3*p3 + p0*p0 ) );
      Scalar rp2 = Scalar(1) / std::max( Scalar(1), sqrt( p5*p5 + p6*p6 + p7*p7 + p4*p4 + p9*p9 + p10*p10 + p11*p11 + p8*p8) );

      // out:
      pq[ i           ] = p0  * rp1;
      pq[ i+    pixel ] = p1  * rp1;
      pq[ i+  2*pixel ] = p2  * rp1;
      pq[ i+  3*pixel ] = p3  * rp1;
      pq[ i+  4*pixel ] = p4  * rp2;
      pq[ i+  5*pixel ] = p5  * rp2;
      pq[ i+  6*pixel ] = p6  * rp2;
      pq[ i+  7*pixel ] = p7  * rp2;
      pq[ i+  8*pixel ] = p8  * rp2;
      pq[ i+  9*pixel ] = p9  * rp2;
      pq[ i+ 10*pixel ] = p10 * rp2;
      pq[ i+ 11*pixel ] = p11 * rp2;
    }
    reproject_Constraint( 12*pixel );
  }

  void reproject_Constraint( int startPos )
  {
    // constraints have to be in succesive order - which means
    // horizontal flow at pixel i 
    // vertical  flow at pixel i 
    // horizontal flow at pixel j 
    // vertical  flow at pixel j , ... etc.
    for (int i=0;i< nCons;i+=2)
    {

      Scalar d1 = pq [ startPos+i   ] - cons_rhs[i];
      Scalar d2 = pq [ startPos+i+1 ] - cons_rhs[i+1];

      Scalar rp = Scalar(1) / std::max( Scalar(1), sqrt(d1*d1 + d2*d2));
      pq [ startPos+i ]   =  d1*rp;
      pq [ startPos+i+1 ] =  d2*rp;
    }
  }


  void copy_u(Scalar*& _out)
  {
    for (int i=0;i<pixel;i++)
      _out[i] = uvw[i];
  };

  void copy_v(Scalar*& _out)
  {
    for (int i=0;i<pixel;i++)
      _out[i] = uvw[i+pixel];
  };
    
  void copy_w(Scalar*& _out)
  {
    for (int i=0;i<4*pixel;i++)
      _out[i] = uvw[i+2*pixel];
  };

  void copy_p(Scalar*& _out)
  {
    for (int i=0;i<4*pixel;i++)
      _out[i] = pq[i];
  };
    
  void copy_q(Scalar*& _out)
  {
    for (int i=0;i<8*pixel;i++)
      _out[i] = pq[i+4*pixel];
  };

  void copy_tau(Scalar*& _out)
  {
    for (int i=0;i<pixel;i++)
      _out[i] = tau[i];
  };


  ////////////////////////
private:

#ifndef _noeigen_
  Eigen::SparseMatrix<Scalar> mat;
  Eigen::SparseMatrix<Scalar> matT;
  Eigen::VectorXd uvw;
  Eigen::VectorXd uvw_;
  Eigen::VectorXd pq;
#else
  Own::SparseMatrix<Scalar> mat;
  Own::SparseMatrix<Scalar> matT;//dummy
  Own::Vector<Scalar> uvw;
  Own::Vector<Scalar> uvw_;
  Own::Vector<Scalar> pq;
#endif
  /// to apply on proximal map 
  std::vector<Scalar> tau;

  int rows;
  int cols;
  int elements;
  int pixel;

  /// potential constraints - flow = uv - in the matrix
  Scalar* cons;
  /// after preconditioning
  std::vector<Scalar> cons_rhs;
  /// number of constraints
  int nCons;
};

#endif // _TGV_SPARSE_MATRIX_H