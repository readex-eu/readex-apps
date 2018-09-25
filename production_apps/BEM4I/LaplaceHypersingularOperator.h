/*!
 * @file    LaplaceHypersingularOperator.h
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    November 21, 2014
 * @brief   Header file for class LaplaceHypersingularOperator
 */

#ifndef LAPLACEHYPERSINGULAROPERATOR
#define	LAPLACEHYPERSINGULAROPERATOR

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
class LaplaceHypersingularOperator : public LinearOperator<LO, SC> {

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! todo: assembly without V
  /*!
   * @brief Constructor with BESpace and user supplied matrix V (p1p1)
   * 
   * @param space  BESpace specifying surface mesh and test and ansatz spaces  
   * @param V      user supplied single-layer operator matrix
   */
  LaplaceHypersingularOperator(
      BESpace< LO, SC > * space,
      Matrix< LO, SC > * V
      );

  /*!
   * @brief Returns true if regularization by a'*a is applied
   */
  bool isRegularized( ) const {
    return regularized;
  }

  /*!
   * @brief Regularization of D by a'*a
   * 
   * @param regularized  set to true if regularization by a'*a should be applied
   */
  void setRegularized(
      bool regularized = true
      ) {
    this->regularized = regularized;
  }

//  /*!
//   * @brief Get hypersingular operator matrix D (+ a'*a if regularized) 
//   * 
//   * @param empty matrix
//   */
//  void getMatrix(
//      FullMatrix< LO, SC > & D
//      );

  /*!
   * @brief Assembles transformation matrix T and stabilizing vector a
   */
  void assemble( ) {
    this->assembleT( );
    if ( this->regularized ) {
      this->assembleA( );
    }
  }

  /*!
   * @brief Applies operator on a vector ( y = beta * y + alpha * this * x )
   *  
   * @param x      operator argument
   * @param y      result (user pre-allocated)
   * @param alpha  y = beta * y + alpha * this * x
   * @param beta   y = beta * y + alpha * this * x
   */
  void apply(
      Vector<LO, SC> const & x,
      Vector<LO, SC> & y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

  //! Destructor
  virtual ~LaplaceHypersingularOperator( );

private:

  //! BESpace specifying mesh and ansatz and test spaces
  BESpace<LO, SC> * space;

  //! Single-layer operator matrix
  Matrix< LO, SC > * V;

  // todo: should be SCVT!
  //! Transformation matrix
  SparseMatrix< LO, SC > * T;

  // todo: should be SCVT!
  //! Stabilization vector
  Vector< LO, SC > * a;

  //! Regularization by a'*a, a[i] = \int \phi_i 
  bool regularized;

  //! Default constructor (empty)
  LaplaceHypersingularOperator( );

  //! Copy constructor (empty)
  LaplaceHypersingularOperator(
      const LaplaceHypersingularOperator& orig
      );

  //! Assembles transformation matrix T
  void assembleT( );

  //! Assembles regularization vector a
  void assembleA( );

  //! Applies operator for p1p1 setting
  void applyP1P1(
      Vector<LO, SC> const & x,
      Vector<LO, SC> & y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

};

}
// include .cpp file to overcome linking problems due to templates
#include "LaplaceHypersingularOperator.cpp"

#endif	/* LAPLACEHYPERSINGULAROPERATOR */
