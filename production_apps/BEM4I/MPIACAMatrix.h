/*!
 * @file    MPIACAMarix.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    March 13, 2015
 * @brief   Header file for class MPIACAMatrix
 * 
 */

#ifndef MPIACAMATRIX_H
#define	MPIACAMATRIX_H

#include "ACAMatrix.h"

namespace bem4i {

template<class LO, class SC>
class MPIACAMatrix : public ACAMatrix<LO, SC> {

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! default constructor
  MPIACAMatrix( );

  /*!
   * Constructor allocating a full matrix
   * 
   * @param[in]   nRows number of rows
   * @param[in]   nCols number of columns
   */
  MPIACAMatrix(
      LO nRows,
      LO nCols,
      MPI_Comm communicator
      );


  //! destructor
  virtual ~MPIACAMatrix( );


  //! applies matrix to a vector
  virtual void apply(
      Vector<LO, SC> const &x,
      Vector<LO, SC> &y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

  //! returns a communicator with which the matrix is associated

  inline MPI_Comm getCommunicator( ) {
    return communicator;
  }

private:

  //! copy constructor
  MPIACAMatrix( const MPIACAMatrix& orig );

  //! communicator of the matrix
  MPI_Comm communicator;


};

}

// include .cpp file to overcome linking problems due to templates
#include "MPIACAMatrix.cpp"

#endif	/* MPIACAMATRIX_H */

