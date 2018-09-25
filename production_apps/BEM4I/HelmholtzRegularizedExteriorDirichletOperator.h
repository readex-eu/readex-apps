/*!
 * @file    HelmholtzRegularizedExteriorDirichletOperator.h
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    May 15, 2015
 * @brief   Header file for class HelmholtzRegularizedExteriorDirichletOperator
 */

#ifndef HELMHOLTZREGULARIZEDEXTERIORDIRICHLETOPERATOR_H
#define	HELMHOLTZREGULARIZEDEXTERIORDIRICHLETOPERATOR_H

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
class HelmholtzRegularizedExteriorDirichletOperator :
public LinearOperator<LO, SC> {

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! todo: assembly without V, K, D, M
  /*!
   * @brief Constructor with BESpace and user supplied matrix V (p1p1)
   * 
   * @param space  BESpace specifying surface mesh and test and ansatz spaces  
   * @param V      user supplied Helmholtz  single-layer operator
   * @param K      user supplied Helmholtz double-layer operator
   * @param Dlap   user supplied REGULARIZED Laplace hypersingular operator
   */
  HelmholtzRegularizedExteriorDirichletOperator(
      BESpace< LO, SC > * space,
      LinearOperator< LO, SC > * V,
      LinearOperator< LO, SC > * K,
      LinearOperator< LO, SC > * Dlap,
      LinearOperator< LO, SC > * M
      );

  //! Returns epsilon for CG solver for inv(D)

  SCVT getEpsSingleLayer( ) const {
    return this->epsHypersingular;
  }

  /*!
   * @brief Sets precision for the application of inv(D)
   * 
   * @param epsHypersingular epsilon for CG solver for inv(D)
   */
  void setEpsHypersingular(
      SCVT epsHypersingular
      ) {
    this->epsHypersingular = epsHypersingular;
  }

  //! Returns maximumnumber of iterations for CG solver for inv(D)

  LO getMaxItHypersingular( ) const {
    return this->maxItHypersingular;
  }

  /*!
   * @brief Sets maximum number of iterations for the application of inv(D)
   * 
   * @param maxItHypersingular maximum number of iterations for inv(D)
   */
  void setMaxItHypersingular(
      LO maxItHypersingular
      ) {
    this->maxItHypersingular = maxItHypersingular;
  }

  //! Returns eta for coupling

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
   * @brief Get Dlap^-1(1/2I+K*_-k)density
   *  
   * @param density   operator argument
   * @param v         result
   */
  void getV(
      const Vector< LO, SC > & density,
      Vector< LO, SC > & v
      ) const;

  //! Destructor
  virtual ~HelmholtzRegularizedExteriorDirichletOperator( );


private:

  //! BESpace specifying mesh and ansatz and test spaces
  BESpace<LO, SC> * space;

  //! wave number
  SC kappa;

  //! eta for linear combination of operators
  SCVT eta;

  //! Single-layer operator
  LinearOperator< LO, SC > * V;

  //! Double-layer operator
  LinearOperator< LO, SC > * K;

  //! Hypersingular operator
  LinearOperator< LO, SC > * Dlap;

  //! Identity operator
  LinearOperator< LO, SC > * M;

  //! Precision for inv(V)
  SCVT epsHypersingular;

  //! Maximum number of iterations for inv(V)
  LO maxItHypersingular;

  //! Default constructor (empty)
  HelmholtzRegularizedExteriorDirichletOperator( );

  //! Copy constructor (empty)
  HelmholtzRegularizedExteriorDirichletOperator(
      const HelmholtzRegularizedExteriorDirichletOperator & orig
      );

  //! Applies operator for p1p1 setting
  void applyP0P0(
      Vector<LO, SC> const & x,
      Vector<LO, SC> & y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

  //! computes w for p1p1 setting
  void getVP0P0(
      const Vector< LO, SC > & density,
      Vector< LO, SC > & w
      ) const;

};

}
// include .cpp file to overcome linking problems due to templates
#include "HelmholtzRegularizedExteriorDirichletOperator.cpp"

#endif	/* HELMHOLTZREGULARIZEDEXTERIORDIRICHLETOPERATOR_H */
