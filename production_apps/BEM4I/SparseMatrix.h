/*!
 * @file    SparseMatrix.h
 * @author  Michal Merta 
 * @date    November 22, 2013
 * @brief   Header file for the SparseMatrix class
 * 
 */

#ifndef SPARSEMATRIX_H
#define	SPARSEMATRIX_H

#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include "Eigen/Sparse"
#include "Eigen/Core"
#include "Eigen/SparseLU"
#include "Eigen/SparseQR"
#include "Eigen/IterativeLinearSolvers"

#include "Matrix.h"
#include "Macros.h"



// optional packages
#ifdef HAS_PARDISO
#include "Eigen/PardisoSupport"
#endif

namespace bem4i {

/*! 
 * Class representing a sparse matrix
 * 
 * Class serves as a wrapper on Eigen SparseMatrix<_Scalar, _Options, _Index>
 * 
 */
template<class LO, class SC, int storage = Eigen::ColMajor>
class SparseMatrix : public Matrix<LO, SC> {

public:

  //! default constructor
  SparseMatrix( );

  //! copy constructor
  SparseMatrix(
      const SparseMatrix & orig
      );

  /*!
   * Constructor allocating a sparse matrix
   * 
   * @param[in]   nRows number of rows
   * @param[in]   nCols number of columns
   */
  SparseMatrix(
      LO nRows,
      LO nCols
      );

  /*!
   * Constructor allocating a sparse matrix with preallocated elements
   * 
   * @param[in]   nRows number of rows
   * @param[in]   nCols number of columns
   * @param[in]   nnzPerRow number of non zero elements per matrix row
   */
  SparseMatrix(
      LO nRows,
      LO nCols,
      LO reserveSize
      );

  /*!
   * Constructor assembling matrix from triplets (i, j, v)
   * 
   * @param[in]   nRows number of rows
   * @param[in]   nCols number of columns
   * @param[in]   rowInd indices of rows
   * @param[in]   colInd indices of cols
   * @param[in]   values values to insert to given indices
   */
  SparseMatrix(
      LO nRows,
      LO nCols,
      std::vector<LO> & rowInd,
      std::vector<LO> & colInd,
      std::vector<SC> & values );

  //! destructor
  virtual ~SparseMatrix( );

  //! returns pointer to the underlying Eigen::SparseMatrix

  inline Eigen::SparseMatrix<SC, storage, LO>* getEigenSparseMatrix( ) const {
    return this->matrix;
  }

  inline LO nnz( ) const {
    return this->matrix->nonZeros( );
  }

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
      const Vector<LO, SC> & x,
      Vector<LO, SC> & y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );


  //! method compresses the matrix into CSR or CCS format

  inline void makeCompressed( ) {
    if ( !this->matrix->isCompressed( ) ) {
      this->matrix->makeCompressed( );
    }
  }

  //! sets the (m,n) element to the value of val

  inline void set( LO m, LO n, SC val ) {
    if ( std::abs( val ) > 0.0 ) {
      this->matrix->coeffRef( m, n ) = val;
    }
  }

  //! adds val to the (m,n) element 

  inline void add( LO m, LO n, SC val ) {
    if ( std::abs( val ) > 0.0 ) {
      this->matrix->coeffRef( m, n ) += val;
    }
  }

  //! returns an (m,n) element of a matrix

  inline SC get( LO m, LO n ) const {
    return this->matrix->coeffRef( m, n );
  }

  //! scales all values of matrix by alpha

  inline void scale( SC alpha ) {
    *( this->matrix ) *= alpha;
  }

  //! computes Frobenius norm of a matrix

  inline SC normFro( ) const {
    return this->matrix->norm( );
  }

  //! resize matrix to size (m, n) and sets all values to zero

  inline void resize( LO m, LO n ) {
    this->nRows = m;
    this->nCols = n;
    this->matrix->resize( m, n );
  }

  /*!
   * @brief Matrix-matrix addition
   * 
   * Computes the sum this = this + alpha*A
   * @param A 
   * @param alpha
   */
  inline void add( SparseMatrix<LO, SC, storage> &A, SC alpha = 1.0 ) {
    *( this->matrix ) += alpha * ( *A.getEigenSparseMatrix( ) );
  }

  void addToPositions( std::vector<LO> const & rows, std::vector<LO> const & cols, FullMatrix<LO, SC> const & mat );

  /*!
   * @brief Set from triplets
   * 
   * @param nRows 
   * @param nCols
   * @param rowInd
   * @param colInd
   * @param values
   */
  void setFromTriplets(
      LO nRows,
      LO nCols,
      std::vector<LO> &rowInd,
      std::vector<LO> &colInd,
      std::vector<SC> &values
      );

  /*!
   * @brief Set from triplets
   * 
   * @param nRows 
   * @param nCols
   * @param rowInd
   * @param colInd
   * @param values
   */
  void setFromTriplets(
      LO nRows,
      LO nCols,
      LO nnz,
      LO *rowInd,
      LO *colInd,
      SC *values
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
  void multiply( SparseMatrix<LO, SC> &A, SparseMatrix<LO, SC> &B, bool transA = false,
      bool transB = false, SC alpha = 1.0, SC beta = 0.0 );

  void print( std::ostream &stream = std::cout ) const {
    for ( int k = 0; k < matrix->outerSize( ); ++k )
      for ( typename Eigen::SparseMatrix<SC, storage, LO>::InnerIterator it( *matrix, k ); it; ++it ) {
        std::cout << it.row( ) << " "; // row index
        std::cout << it.col( ) << " "; // col index (here it is equal to k)
        std::cout << it.value( ) << std::endl;
      }
  }

  //! saves matrix triplets in binary file
  void saveTripletsBin( const std::string &fileName ) const;

  //! loads matrix from binary file containing triplets
  void loadTripletsBin( const std::string &fileName ) const;

  //! loads matrix from a text file containing triplets
  void loadTriplets(
      const std::string &fileName
      );

  //! method solves the system of linear equations using LU factorization
  void LUSolve( const Vector<LO, SC>& rhs, Vector<LO, SC>& x );

  //! method solves the system of linear equations using QR factorization
  void QRSolve( const Vector<LO, SC>& rhs, Vector<LO, SC>& x );

  //! method solves the system of linear equations using BiCGStab solver
  void BiCGStabSolve( Vector<LO, SC>& rhs, Vector<LO, SC>& x );

  inline int getStorageType( ) {
    return storage;
  }

#ifdef HAS_PARDISO
  //! method solves the system of linear equations using PARDISO
  void PARDISOSolve( Vector<LO, SC>& rhs, Vector<LO, SC>& x );
#endif



private:

  //! pointer to Eigen::SparseMatrix<SC, storage, LO>
  Eigen::SparseMatrix<SC, storage, LO> *matrix;

  //! was matrix factorized by Eigen?
  bool factorized;

  //! was matrix QR factorized by Eigen?
  bool QRFactorized;

  //! Eigen sparse lu solver
  Eigen::SparseLU<Eigen::SparseMatrix<SC, storage, LO>, Eigen::COLAMDOrdering<LO> > *solver;

  //! Eigen sparse qr solver
  Eigen::SparseQR<Eigen::SparseMatrix<SC, storage, LO>, Eigen::COLAMDOrdering<LO> > *QRSolver;

};

}
// include .cpp file to overcome linking problems due to templates
#include "SparseMatrix.cpp"


#endif	/* SPARSEMATRIX_H */
