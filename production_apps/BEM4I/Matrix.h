/*!
 * @file    Matrix.h
 * @author  Michal Merta 
 * @date    July 5, 2013
 * @brief   Header file for abstract class Matrix
 * 
 */


#ifndef MATRIX_H
#define	MATRIX_H

#include "LinearOperator.h"
#include "IterativeSolver.h"

namespace bem4i {

/*! 
 * Abstract class representing matrix-type operator
 * 
 */
template<class LO, class SC>
class Matrix : public LinearOperator<LO, SC> {

public:

  typedef typename GetType<LO, SC>::SCVT SCVT;

  //! default constructor
  Matrix( );

  //! destructor
  virtual ~Matrix( );

  //! returns a number of rows

  inline LO getNRows( ) const {
    return nRows;
  }

  //! returns a number of cols

  inline LO getNCols( ) const {
    return nCols;
  }

  //! returns a number of rows

  inline void setNRows(
      LO nRows
      ) {
    this->nRows = nRows;
  }

  //! returns a number of cols

  inline void setNCols(
      LO nCols
      ) {
    this->nCols = nCols;
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
      Vector<LO, SC> const &x,
      Vector<LO, SC> &y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      ) = 0;

  virtual void print( std::ostream &stream = std::cout ) const = 0;

protected:
  
  // {} instead of =, workaround because of bug in gcc
  //! number of rows
  LO & nRows{LinearOperator< LO, SC >::dimRange}; //this->dimRange;//

  //! number of cols
  LO & nCols{LinearOperator< LO, SC >::dimDomain};//this->dimDomain;

private:

};

}

// include .cpp file to overcome linking problems due to templates
#include "Matrix.cpp"

#endif	/* MATRIX_H */

