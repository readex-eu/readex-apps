/*!
 * @file    FullMatrix.h
 * @author  Michal Merta 
 * @date    July 5, 2013
 * @brief   Header file for the FullMatrix class
 * 
 */

#ifndef FULLMATRIX_H
#define FULLMATRIX_H

#include <cstring>
#include <vector>

#include "Matrix.h"
#include "SparseMatrix.h"
#include "Macros.h"

#if N_MIC > 0
#pragma offload_attribute( push, target( mic ) )
#endif
#include "BLAS_LAPACK_wrapper.h"
#if N_MIC > 0
#pragma offload_attribute( pop )
#endif


namespace bem4i {

  /*! 
   * Class representing a full matrix
   * 
   * The matrix is stored in a column-major order.
   * Matrix operations are performed by calling the BLAS/LAPACK routines.
   * 
   */
  template<class LO, class SC>
  class FullMatrix : public Matrix<LO, SC>, Offloadable {
    typedef typename GetType<LO, SC>::SCVT SCVT;

    friend class FullMatrix<LO, SCVT>;
    friend class FullMatrix<LO, std::complex<SCVT> >;

  public:

    //! default constructor
    FullMatrix( );

    //! copy constructor
    FullMatrix(
      const FullMatrix& orig
      );

    /*!
     * Constructor allocating a full matrix
     * 
     * @param[in]   nRows number of rows
     * @param[in]   nCols number of columns
     * @param[in]   zeroOut (optional) whether to set all matrix elements to zero
     * @param[in]   allocWork (optional) whether to allocate working space for BLAS/LAPACK
     */
    FullMatrix(
      LO nRows,
      LO nCols,
      bool zeroOut = true,
      bool allocWork = true
      );

    /*!
     * Constructor taking an user preallocated array with matrix data
     * 
     * @param[in]   nRows number of rows
     * @param[in]   nCols number of columns
     * @param[in]   data pointer to an array representing matrix entries
     * @param[in]   allocWork (optional) whether to allocate working space for BLAS/LAPACK
     */
    FullMatrix(
      LO nRows,
      LO nCols,
      SC * data,
      bool allocWork = true
      );

    //! destructor
    virtual ~FullMatrix( );

    //! deletes matrix values and allocates new resized matrix
    void resize(
      LO nRows,
      LO nCols,
      bool zeroOut = true
      );

    void copy(
      FullMatrix<LO, SC> & copy
      ) const;

    void copyToComplex(
      FullMatrix< LO, std::complex< SCVT > > & copy
      ) const;

    //! returns an (m,n) element of a matrix

    inline SC get(
      LO m,
      LO n
      ) const {
      return data[n * this->nRows + m];
    }

    //! sets the (m,n) element to the value of val

    inline void set(
      LO m,
      LO n,
      SC val
      ) {
      data[n * this->nRows + m] = val;
    }

    //! adds val to the (m,n) element

    inline void add(
      LO m,
      LO n,
      SC val
      ) {
      data[n * this->nRows + m] += val;
    }

    inline void addAtomic(
      LO m,
      LO n,
      SC val
      ) {

#pragma omp atomic update
      data[n * this->nRows + m] += val;
    }

    //! sets all elements to the value of val

    inline void setAll(
      SC val
      ) {

      SC *p = this->data, *last = this->data + this->nRows * this->nCols;
      while ( p != last ) *( p++ ) = val;
    }

    /*!
     * Sums values specified by input array to the subset of  matrix 
     * 
     * @param[in]   rows  indices of rows
     * @param[in]   cols  indices of cols
     * @param[in]   mat   input matrix
     */
    void addToPositions(
      std::vector<LO> const & rows,
      std::vector<LO> const & cols,
      FullMatrix<LO, SC> const & mat
      );

    void addToPositionsAtomic(
      std::vector<LO> const & rows,
      std::vector<LO> const & cols,
      FullMatrix<LO, SC> const & mat
      );

    /*!
     * Sums values specified by input array to the subset of  matrix to positions 
     * (i, j, v) - i specified by rows, j by Cols
     * 
     * @param[in]   rows    indices of rows
     * @param[in]   cols    indices of cols
     * @param[in]   values  input values
     */
    void addToPositions(
      std::vector<LO> const & rows,
      std::vector<LO> const & cols,
      std::vector<SC> const & values
      );

    //! copies a column of a matrix into a <em>user-preallocated</em> array
    void getCol(
      LO idx,
      SC * outCol
      ) const;

    //! copies a row of a matrix into a <em>user-preallocated</em> array
    void getRow(
      LO idx,
      SC * outCol
      ) const;

    //! scales all values of matrix by alpha
    void scale(
      SC alpha
      );

    //! performs in-place conjugation
    void conjugate( );

    //! computes 1-norm of a matrix
    SC norm1( ) const;

    //! computes Frobenius norm of a matrix
    SC normFro( );

    //! computes Inf-norm of a matrix
    SC normI( );

    /*!
     * @brief Matrix-matrix addition
     * 
     * Computes the sum this = this + alpha*A
     * @param A 
     * @param alpha
     */
    void add(
      FullMatrix<LO, SC> &A,
      SC alpha = 1.0 );

    /*!
     * @brief Matrix-matrix addition
     * 
     * Computes the sum this = this + alpha*A
     * @param A 
     * @param alpha
     */
    void add(
      SparseMatrix<LO, SC, Eigen::ColMajor> &A,
      SC alpha = 1.0 );

    /*!
     * @brief Matrix-matrix multiplication
     * 
     * Computes a sum this = beta*this + alpha*A*B
     * @param A 
     * @param B
     * @param alpha
     * @param beta
     */
    void multiply(
      FullMatrix<LO, SC> &A,
      FullMatrix<LO, SC> &B,
      bool transA = false,
      bool transB = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

    /*!
     * @brief Matrix-matrix multiplication
     * 
     * Computes a sum this = beta*this + alpha*A*B
     * @param A 
     * @param B
     * @param alpha
     * @param beta
     */
    void multiply(
      FullMatrix<LO, SC> &A,
      SparseMatrix<LO, SC, Eigen::ColMajor> &B,
      bool transA = false,
      bool transB = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

    /*!
     * @brief Matrix-matrix multiplication
     * 
     * Computes a sum this = beta*this + alpha*A*B
     * @param A 
     * @param B
     * @param alpha
     * @param beta
     */
    void multiply(
      SparseMatrix<LO, SC, Eigen::ColMajor> &A,
      FullMatrix<LO, SC> &B,
      bool transA = false,
      bool transB = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

    /*!
     * @brief Matrix-matrix multiplication
     * 
     * Computes a sum this = beta*this + alpha*A(1:nRows, 1:nCols)*B(1:nCols, :)
     * @param A 
     * @param B
     * @param alpha
     * @param beta
     */
    void multiply(
      FullMatrix<LO, SC> &A,
      FullMatrix<LO, SC> &B,
      LO ARows,
      LO ACols,
      LO BRows,
      LO BCols,
      bool transA,
      bool transB,
      SC alpha = 1.0,
      SC beta = 0.0 );

    /*!
     * @brief Performs a matrix-vector multiplication
     * 
     * Computes y = beta*y + alpha*this*x
     * @param A 
     * @param x
     * @param y
     * @param alpha
     * @param beta
     */
    virtual void apply(
      Vector<LO, SC> const & x,
      Vector<LO, SC> & y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

    void applyMIC(
      Vector<LO, SC> const & x,
      Vector<LO, SC> & y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0,
      int device = 0
      );

    /*!
     * @brief Performs a matrix-vector multiplication using submatrix
     * 
     * Computes y = beta*y + alpha*this(1:ARows, 1:ACols)*x
     * @param A 
     * @param x
     * @param y
     * @param ARows
     * @param ACols
     * @param alpha
     * @param beta
     */
    void applySubmatrix(
      Vector<LO, SC> const &x,
      Vector<LO, SC> &y,
      LO ARows,
      LO ACols,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0 );

    /*!
     * @brief Performs a matrix-vector multiplication on a subvector
     * 
     * Computes y(indicesY) = this*x(indicesX)
     * @param A 
     * @param x
     * @param y
     */
    void apply(
      Vector<LO, SC> const &x,
      std::vector<LO>& indicesX,
      Vector<LO, SC> &y,
      std::vector<LO> &indicesY,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

    /*!
     * @brief Solves a system of linear equations using LU decomposition
     * 
     * Solves the system this*x = rhs. WARNING: The original matrix is 
     * overwritten by its factors!
     * @param x on input: right-hand side, on output: result
     */
    void LUSolve(
      Vector<LO, SC> & x,
      LO nRhs = 1
      );

    /*!
     * @brief Solves a system of linear equations using Choleski decomposition
     * 
     * Solves the system this*x = rhs. WARNING: The original matrix is 
     * overwritten by its factors!
     * @param x on input: right-hand side, on output: result
     */
    void CholeskiSolve(
      Vector<LO, SC> & x,
      LO nRhs = 1
      );

    /*!
     * @brief Solves a system of linear equations using LU decomposition
     * 
     * Solves the system this*x = rhs. WARNING: The original matrix is 
     * overwritten by its factors!
     * @param x on input: right-hand side, rhs given by columns of x
     */
    void LUSolve(
      FullMatrix<LO, SC> & x,
      LO nRhs = 0
      );

    /*!
     * @brief Solves a system of linear equations using Choleski decomposition
     * 
     * Solves the system this*x = rhs. WARNING: The original matrix is 
     * overwritten by its factors!
     * @param x on input: right-hand side, rhs given by columns of x
     */
    void CholeskiSolve(
      FullMatrix<LO, SC> & x,
      LO nRhs = 0
      );


    /*!
     * @brief Solves a system with upper triangular matrix by backward substitution
     * 
     * Solves the system this*x = rhs. WARNING: The original matrix is 
     * overwritten by its factors!
     * @param x     on input: right-hand side, on output: resul
     * @param nRhs  number of right-hand sides
     * @param n     dimension of a system matrix this(1:n, 1:n), if 0 => n = this.nRows
     */
    void backward(
      Vector<LO, SC> & x,
      LO nRhs = 1,
      LO n = 0
      );

    /*!
     * @brief Computes eigenvectors and eigenvalues of symmetric real matrix
     * 
     * todo: Not implemented for complex matrices!
     * @param[in,out]   eigenvectors  array of eigenvectors
     * @param[in,out]   eigenvalues   array of eigenvalues
     */
    void eigs(
      SC * eigenvectors,
      SC * eigenvalues
      );

    /*!
     * @brief Computes the Schur decompositin of a Hessenberg matrix
     * 
     * Computes decomposition of this matrxi into H = Z T Z**H. The input matrix 
     * must be of Hessenberg form.
     * 
     * WARNING: May rewrite the original matrix
     * 
     * @param[in, out]  Z  reference to the full matrix Z from Schur decomposition
     * @param[out]      eigs pointer to user preallocated array to store 
     *                  eigenvalues 
     */
    void schurOfHessenberg(
      FullMatrix<LO, SC> & Z,
      std::complex<SCVT> * eigs
      );

    //! prints the matrix
    void print(
      std::ostream & stream = std::cout
      ) const;

    /*!
     * @brief Returns a pointer to the array of data. Use with caution!
     * 
     * 
     * @param[in, out]  data
     */
    inline SC * getData( ) const {
      return data;
    }

    //  void setFactReuseOn( ) {
    //    this->reuseFact = true;
    //  }

    void xferToMIC(
      int device = 0
      ) const {

#if N_MIC > 0

      const SCVT * data = reinterpret_cast < SCVT * > ( this->data );
      const SCVT * work = reinterpret_cast < SCVT * > ( this->work );
      int mult = 1;
      if ( !std::is_same< SC, SCVT >::value ) mult = 2;

      LO nRows = this->nRows;
      LO nCols = this->nCols;

#pragma omp target enter data device( device ) \
map( to : data[ 0 : mult * nRows * nCols ] ) \
map( to : work[ 0 : mult * nRows ] )

#endif
    }

    void xferToHost(
      int device = 0
      ) const {

#if N_MIC > 0

      const SCVT * data = reinterpret_cast < SCVT * > ( this->data );
      int mult = 1;
      if ( !std::is_same< SC, SCVT >::value ) mult = 2;

      LO nRows = this->nRows;
      LO nCols = this->nCols;

#pragma omp target update device( device ) \
from( data[ 0 : mult * nRows * nCols ] )

#endif
    }

    void updateMIC(
      int device = 0
      ) const {

#if N_MIC > 0

      const SCVT * data = reinterpret_cast < SCVT * > ( this->data );
      const SCVT * work = reinterpret_cast < SCVT * > ( this->work );
      int mult = 1;
      if ( !std::is_same< SC, SCVT >::value ) mult = 2;

      LO nRows = this->nRows;
      LO nCols = this->nCols;

#pragma omp target update device( device ) \
to( data[ 0 : mult * nRows * nCols ] )

#endif
    }

    void deleteMIC(
      int device = 0
      ) const {

#if N_MIC > 0

      const SCVT * data = reinterpret_cast < SCVT * > ( this->data );
      const SCVT * work = reinterpret_cast < SCVT * > ( this->work );

#pragma omp target exit data device( device ) \
map( delete : data ) \
map( delete : work )

#endif

    }

  private:

    //! 1D array with matrix data in column-major order
    SC * data;

    /*!
     * @brief auxiliary workspace for some LAPACK routines
     * 
     * Note, that some LAPACK routines may need more allocated space.
     */
    mutable SC * work;

    //! whether the destructor should delete array with data
    bool deleteData;

    //! whether the matrix has been factorized for LAPACK LU solve already
    //  bool factorized;

    //! pivoting info from LU factorization
    //  int * ipiv;

    //  bool reuseFact;

  };

}
// include .cpp file to overcome linking problems due to templates
#include "FullMatrix.cpp"


#endif /* FULLMATRIX_H */

