/*!
 * @file    HelmholtzRegularizedExteriorNeumannOperator.h
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    April 7, 2015
 * @brief   Header file for class HelmholtzRegularizedExteriorNeumannOperator
 */

#ifndef HELMHOLTZREGULARIZEDEXTERIORNEUMANNOPERATOR_H
#define	HELMHOLTZREGULARIZEDEXTERIORNEUMANNOPERATOR_H

#include <vector>

#include "LinearOperator.h"
#include "Vector.h"
#include "BESpace.h"
#include "SparseMatrix.h"
#include "IterativeSolver.h"

namespace bem4i {

/*!
 * Class representing a weak operator on the given discretized space
 */
template<class LO, class SC>
class HelmholtzRegularizedExteriorNeumannOperator :
public LinearOperator<LO, SC> {

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! todo: assembly without V, K, D, M
  /*!
   * @brief Constructor with BESpace and user supplied matrix V (p1p1)
   * 
   * @param space  BESpace specifying surface mesh and test and ansatz spaces  
   * @param Vlap      user supplied Laplace single-layer operator
   * @param K      user supplied Helmholtz double-layer operator
   * @param D      user supplied Helmholtz hypersingular operator
   */
  HelmholtzRegularizedExteriorNeumannOperator(
      BESpace< LO, SC > * space,
      LinearOperator< LO, SC > * Vlap,
      LinearOperator< LO, SC > * K,
      LinearOperator< LO, SC > * D,
      LinearOperator< LO, SC > * M
      );

  //! Returns epsilon for CG solver for inv(V)

  SCVT getEpsSingleLayer( ) const {
    return this->epsSingleLayer;
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
    return this->maxItSingleLayer;
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

  //! Returns epsilon for coupling

  SCVT getEta( ) const {
    return this->eta;
  }

  /*!
   * @brief Sets coupling parameter
   * 
   * @param eta coupling parameter
   */
  void setEta(
      SCVT eta
      ) {
    this->eta = eta;
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

  /*!
   * @brief Get Vlap^-1(-1/2I+K_-k)density
   *  
   * @param density   operator argument
   * @param w         result
   */
  void getW(
      const Vector< LO, SC > & density,
      Vector< LO, SC > & w
      ) const;

  //! Destructor
  virtual ~HelmholtzRegularizedExteriorNeumannOperator( );


private:

  //! BESpace specifying mesh and ansatz and test spaces
  BESpace<LO, SC> * space;

  //! wave number
  SC kappa;

  //! eta for linear combination of operators
  SCVT eta;

  //! Single-layer operator
  LinearOperator< LO, SC > * Vlap;

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
  HelmholtzRegularizedExteriorNeumannOperator( );

  //! Copy constructor (empty)
  HelmholtzRegularizedExteriorNeumannOperator(
      const HelmholtzRegularizedExteriorNeumannOperator& orig
      );

  //! Applies operator for p1p1 setting
  void applyP1P1(
      Vector<LO, SC> const & x,
      Vector<LO, SC> & y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

  //! computes w for p1p1 setting
  void getWP1P1(
      const Vector< LO, SC > & density,
      Vector< LO, SC > & w
      ) const;

};

}
// include .cpp file to overcome linking problems due to templates
#include "HelmholtzRegularizedExteriorNeumannOperator.cpp"

#endif	/* HELMHOLTZREGULARIZEDEXTERIORNEUMANNOPERATOR_H */
