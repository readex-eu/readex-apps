/*!
 * @file    IdentityOperator.h
 * @author  Michal Merta 
 * @date    July 18, 2013
 * @brief   Header file for class IdentityOperator
 * 
 */


#ifndef IDENTITYOPERATOR_H
#define	IDENTITYOPERATOR_H

#include <vector>

#include "LinearOperator.h"
#include "Vector.h"
#include "BESpace.h"
#include "SparseMatrix.h"
#include "IterativeSolver.h"

namespace bem4i {

/*!
 * Class representing a weak identity operator on given discretized space
 */
template<class LO, class SC>
class IdentityOperator : public LinearOperator<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! constructor taking a boundary element space as argument
  IdentityOperator(
      BESpace<LO, SC> * space
      );

  /*!
   * @brief Applies operator on a vector
   * 
   * Computes y = beta*y + alpha*this*x
   * @param A 
   * @param x
   * @param y
   * @param alpha
   * @param beta
   */
  void apply(
      Vector<LO, SC> const &x,
      Vector<LO, SC> &y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

  void assemble(
      SparseMatrix< LO, SC > & M
      ) const;

  //! destructor
  virtual ~IdentityOperator( );

private:

  BESpace<LO, SC> * space;
  
    //! default constructor
  IdentityOperator( ) {};

  //! copy constructor
  IdentityOperator(
      const IdentityOperator & orig
      ) {};

  //! applies operator for combination p0, p0
  void applyP0P0(
      Vector<LO, SC> const &x,
      Vector<LO, SC> &y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

  //! applies operator for combination p0, p1
  void applyP0P1(
      Vector<LO, SC> const &x,
      Vector<LO, SC> &y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

  //! applies operator for combination p1, p0
  void applyP1P0(
      Vector<LO, SC> const &x,
      Vector<LO, SC> &y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

  //! applies operator for combination p1, p1
  void applyP1P1(
      Vector<LO, SC> const &x,
      Vector<LO, SC> &y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

  void assembleP0P1(
      SparseMatrix< LO, SC > & M
      ) const;
  
  void assembleP1P0(
      SparseMatrix< LO, SC > & M
      ) const;
  
   void assembleP1P1(
      SparseMatrix< LO, SC > & M
      ) const;
   
   void assembleP0P0(
      SparseMatrix< LO, SC > & M
      ) const;

};

}
// include .cpp file to overcome linking problems due to templates
#include "IdentityOperator.cpp"

#endif	/* IDENTITYOPERATOR_H */

