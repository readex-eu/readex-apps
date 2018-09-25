/*!
 * @file    MPIBlockMatrix.h
 * @author  Michal Merta 
 * @date    March 10, 2014
 * @brief   Header file for the MPIBlockMatrix class
 * 
 */

#ifndef MPIBLOCKMATRIX_H
#define	MPIBLOCKMATRIX_H

#include "BlockMatrix.h"
//#include <mpi.h>

namespace bem4i {

/*! 
 * Class representing a distributed block matrix
 * 
 * 
 */
template<class LO, class SC>
class MPIBlockMatrix : public BlockMatrix<LO, SC> {

public:

  //! default constructor
  MPIBlockMatrix( );

  /*! Constructor taking number of blocks their sizes and rank assignments as arguments
   * 
   * @param[in]   nBlockRows  number of block row
   * @param[in]   nBlockCols  number of block columns
   * @param[in]   numbersOfRows array of legth nBlockRows containing numbers of rows of each block row
   * @param[in]   numbersOfCols array of legth nBlockCols containing numbers of columns of each block column
   * @param[in]   ranks   array of length nBlockRows*nBlockCols containing the assignment of blocks to MPI ranks
   */
  MPIBlockMatrix(
      int nBlockRows,
      int nBlockCols,
      LO * numbersOfRows,
      LO * numbersOfCols,
      int * ranks,
      MPI_Comm communicator
      );

  //! destructor
  virtual ~MPIBlockMatrix( );

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
      Vector<LO, SC> const &x,
      Vector<LO, SC> &y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

  //! returns owner of block on position (i,j)

  inline int getOwner(
      LO i,
      LO j
      ) {
    return this->ranks[j * this->nBlockRows + i];
  }

  //! returns a communicator with which the matrix is associated

  inline MPI_Comm getCommunicator( ) {
    return communicator;
  }

  /*! Method resizes a matrix according to the given parameters. 
   * 
   *  If no rank specifying array is given, all blocks will be located on rank 0!
   * 
   * @param[in]   nBlockRows  number of block row
   * @param[in]   nBlockCols  number of block columns
   * @param[in]   numbersOfRows array of legth nBlockRows containing numbers of rows of each block row
   * @param[in]   numbersOfCols array of legth nBlockCols containing numbers of columns of each block column
   */
  virtual void resize(
      int nBlockRows,
      int nBlockCols,
      LO * numbersOfRows,
      LO * numbersOfCols
      );


  /*! Method resizes a matrix according to the given parameters. 
   * 
   * 
   * @param[in]   nBlockRows  number of block row
   * @param[in]   nBlockCols  number of block columns
   * @param[in]   numbersOfRows array of legth nBlockRows containing numbers of rows of each block row
   * @param[in]   numbersOfCols array of legth nBlockCols containing numbers of columns of each block column
   */
  virtual void resize(
      int nBlockRows,
      int nBlockCols,
      LO * numbersOfRows,
      LO * numbersOfCols,
      int * ranks
      );

protected:

  //! ranks of owners of individual blocks
  int * ranks;

  //! communicator of the matrix
  MPI_Comm communicator;

  //! returns whether the current process owns the block (i,j)

  inline bool amIOwner( LO i, LO j ) {
    int rank;
    MPI_Comm_rank( communicator, &rank );
    return (this->ranks[j * this->nBlockRows + i] == rank );
  }

};

} // end of namespace bem4i

// include .cpp file to overcome linking problems due to templates
#include "MPIBlockMatrix.cpp"
#endif