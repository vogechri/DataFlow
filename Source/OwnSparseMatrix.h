#ifndef _OWN_SPARSE_MATRIX_T
#define _OWN_SPARSE_MATRIX_T

//== INCLUDES =======================================================

//== NAMESPACE ================================================================

namespace Own{

  //== CLASS DEFINITION =========================================================

  /// silly class of a vector
  template <class Scalar>
  class Vector
  {
  public:
    typedef Scalar            value_type;
    typedef Vector<Scalar> type;
    typedef Vector<Scalar> Self;

    Vector(): init(false){ };

    Vector( int nRows ): init(true), x(nRows,0) { };
    
    ~Vector() {};

    Scalar& operator[] (int i) {return x[i];}

    int size() {return x.size();};

//    OwnVector<Scalar>& 
    void addTwiceSub ( Vector<Scalar>& y ) 
    {
      //std::vector<Scalar> z(x.size(),0);
      for (int i=0;i<x.size();i++)
        x[i] = Scalar(2.)*y[i] - x[i];
//      return z;
    }

    //template< typename Scalar >    inline VectorT<Scalar,N> operator*(Scalar _s, const VectorT<Scalar,N>& _v) return VectorT<Scalar,N>(_v)

    void multS(Scalar _s) 
    {
      for (int i=0;i<x.size();i++)
        x[i] *= _s;
    }

    std::vector<Scalar> x;
    bool init;
  };

  /** /class OwnSparseMatrix
  
  Stores a sparse symmetric! matrix efficiently.
  Expects Arrays as input to construct
  Holds diagonal preconditioners as well -- for column and rows.
  **/

  template <class Scalar>
  class SparseMatrix
  {
  public:
    typedef Scalar                  value_type;
    typedef SparseMatrix<Scalar> type;
    typedef SparseMatrix<Scalar> Self;


    SparseMatrix(): init(false), pcg(false)  { };

    SparseMatrix( int nRows, int nCols ): init(false), pcg(false), nrows(nRows), ncols(nCols)  { };
/*
    SparseMatrix( Scalar*& ids_i, Scalar*& ids_j, Scalar*& ids_v, int nRows, int nCols ): init(true), pcg(false), values(ids_v), columns(ids_i), rows(ids_j), nrows(nRows), ncols(nCols)  { };

    SparseMatrix( Scalar*& ids_i, Scalar*& ids_j, Scalar*& ids_v, std::vector<Scalar>& normalizerR, std::vector<Scalar>& normalizerC ): _init(true), _pcg(false), values(ids_v), columns(ids_i), rows(ids_j)  
    {
      pcg_row = normalizerR;
      pcg_col = normalizerC;
    };
*/
    ~SparseMatrix() { };

    /////////////////////////////////

    void setValues ( int _elements, Scalar*& ids_i, Scalar*& ids_j, Scalar*& ids_v, std::vector<Scalar>& normalizerR, std::vector<Scalar>& normalizerC  )
    {
      elements = _elements;
      values = ids_v; columns= ids_i; rows = ids_j;init=true;
      pcg_row = normalizerR;
      pcg_col = normalizerC;

//#pragma omp parallel
      for(int i=0;i < pcg_row.size();i++)
        if ( pcg_row[ i ] > 0 ) 
          pcg_row[ i ] = Scalar(1.) / pcg_row[ i ];
        else
          pcg_row[ i ] = Scalar(0);

// #pragma omp parallel
      for(int i=0;i < pcg_col.size();i++)
        if ( pcg_col[ i ] > 0 ) 
          pcg_col[ i ] = Scalar(1.) / pcg_col[ i ];
        else
          pcg_col[ i ] = Scalar(0);

      pcg = true;
    }

    void setValues ( int _elements, Scalar*& ids_i, Scalar*& ids_j, Scalar*& ids_v )
    {
      elements = _elements;
      values = ids_v; columns= ids_i; rows = ids_j;init=true;
      pcg = false;
    }

    bool initialized() { return init; }

#ifdef _matrixcolMajor_
    /// z = M*x+z
    void multplus( std::vector<Scalar>&x, std::vector<Scalar>&z )
    {
      if (!pcg)
      for (int i=0; i<ncols; i++) // Loop through columns
      {
         int stop = columns[i+1]-1;
         Scalar rhs  = x[i];
         for (int k=columns[i]-1; k<stop; k++) // Loop through non-zeros in ith column
         {
            z[ rows[k]-1 ] += values[k] * rhs;
         }
      }
      else
      for (int i=0; i<ncols; i++) // Loop through columns
      {
         int stop = columns[i+1]-1;
         Scalar rhs  = x[i];
         for (int k=columns[i]-1; k<stop; k++) // Loop through non-zeros in ith column
         {
            z[ rows[k]-1 ] += pcg_row[rows[k]-1] * values[k] * rhs;
         }
      }
    }


    /// z = z-M^tx
    void multTminus( std::vector<Scalar>&x, std::vector<Scalar>&z )
    {
      if (!pcg)
#pragma omp parallel
      for (int i=0; i<ncols; i++) // Loop through columns
      {
         int stop = columns[i+1]-1;
         for (int k=columns[i]-1; k<stop; k++) // Loop through non-zeros in ith column
         {
            z[ i ] -= values[k] * x[ rows[k]-1 ];
         }
      }
      else
#pragma omp parallel
      for (int i=0; i<ncols; i++) // Loop through columns
      {
         int stop = columns[i+1]-1;
         for (int k=columns[i]-1; k<stop; k++) // Loop through non-zeros in ith column
         {
            z[ i ] -= pcg_col[i] * values[k] * x [ rows[k]-1 ];
         }
      }    
    }

#else // just indices -- of course SLOW:

    /// z = M*x+z
    void multplus( std::vector<Scalar>&x, std::vector<Scalar>&z )
    {
      if (!pcg)
        for (int i=0;i<elements;i++)
          z[ columns[i]-1 ] += values[i] * x[ rows[i]-1 ];
      else
        for (int i=0;i<elements;i++)
          z[ columns[i]-1 ] += pcg_row[ columns[i]-1 ] * values[i] * x[ rows[i]-1 ];
    }


    /// z = z-M^tx
    void multTminus( std::vector<Scalar>&x, std::vector<Scalar>&z )
    {
      if (!pcg)
        for (int i=0;i<elements;i++)
          z[ rows[i]-1 ] -= values[i] * x[ columns[i]-1 ];
      else
        for (int i=0;i<elements;i++)
          z[ rows[i]-1 ] -= pcg_col[ rows[i]-1 ] * values[i] * x[ columns[i]-1 ];
    }
#endif


    /// z = M*x+z
    void multplus( Vector<Scalar>&x, Vector<Scalar>&z )
    {
      multplus(x.x, z.x);
    }

    /// z = z-M^tx
    void multTminus( Vector<Scalar>&x, Vector<Scalar>&z )
    {
      multTminus(x.x, z.x);
    }

  private:

    bool init;
    bool pcg; 

    int elements;
    int nrows;
    int ncols;

    Scalar*        values;
    Scalar*        columns;
    Scalar*        rows;

    std::vector<Scalar> pcg_row;
    std::vector<Scalar> pcg_col;

  };
};
#endif
