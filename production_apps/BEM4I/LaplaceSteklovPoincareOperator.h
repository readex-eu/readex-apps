/*!
 * @file    LaplaceSteklovPoincareOperator.h
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    November 24, 2014
 * @brief   Header file for class LaplaceSteklovPoincareOperator (DtN map)
 */

#ifndef LAPLACESTEKLOVPOINCAREOPERATOR_H
#define	LAPLACESTEKLOVPOINCAREOPERATOR_H

#include <vector>

#include "LinearOperator.h"
#include "Vector.h"
#include "BESpace.h"
#include "SparseMatrix.h"
#include "IterativeSolver.h"

namespace bem4i {

/*!
 * Class representing a weak SP operator on given discretized space
 */
template<class LO, class SC>
class LaplaceSteklovPoincareOperator : public LinearOperator<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! todo: assembly without V, K, D, M
  /*!
   * @brief Constructor with BESpace and user supplied matrix V (p1p1)
   * 
   * @param space  BESpace specifying surface mesh and test and ansatz spaces  
   * @param V      user supplied single-layer operator
   * @param K      user supplied double-layer operator
   * @param D      user supplied REGULARIZED hypersingular operator
   */
  LaplaceSteklovPoincareOperator(
      BESpace< LO, SC > * space,
      LinearOperator< LO, SC > * V,
      LinearOperator< LO, SC > * K,
      LinearOperator< LO, SC > * D,
      LinearOperator< LO, SC > * M
      );

  //! Returns epsilon for CG solver for inv(V)

  SCVT getEpsSingleLayer( ) const {
    return epsSingleLayer;
  }

  /*!
   * @brief Sets precision for the application of inv(V)
   * 
   * @param epsSingleLayer epsilon for CG solver for inv(V)
   */
  void setEpsSingleLayer(
      SCVT epsSingleLayer
      ) {
    this->epsSingleLayer = epsSingleLayer;
  }

  //! Returns maximumnumber of iterations for CG solver for inv(V)

  LO getMaxItSingleLayer( ) const {
    return maxItSingleLayer;
  }

  /*!
   * @brief Sets maximum number of iterations for the application of inv(V)
   * 
   * @param maxItSingleLayer maximum number of iterations for inv(V)
   */
  void setMaxItSingleLayer(
      LO maxItSingleLayer
      ) {
    this->maxItSingleLayer = maxItSingleLayer;
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
  virtual ~LaplaceSteklovPoincareOperator( );


private:

  //! BESpace specifying mesh and ansatz and test spaces
  BESpace<LO, SC> * space;

  //! Single-layer operator
  LinearOperator< LO, SC > * V;

  //! Double-layer operator
  LinearOperator< LO, SC > * K;

  //! Hypersingular operator
  LinearOperator< LO, SC > * D;

  //! Identity operator
  LinearOperator< LO, SC > * M;

  //! Precision for inv(V)
  SCVT epsSingleLayer;

  //! Maximum number of iterations for inv(V)
  LO maxItSingleLayer;

  //! Default constructor (empty)
  LaplaceSteklovPoincareOperator( );

  //! Copy constructor (empty)
  LaplaceSteklovPoincareOperator(
      const LaplaceSteklovPoincareOperator& orig
      );

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
#include "LaplaceSteklovPoincareOperator.cpp"

#endif	/* LAPLACESTEKLOVPOINCAREOPERATOR_H */
