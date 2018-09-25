/*!
 * @file    HelmholtzNeumannRobinOperator.h
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    July 30, 2015
 * @brief   Header file for class HelmholtzNeumannRobinOperator
 */

#ifndef HELMHOLTZNEUMANNROBINOPERATOR_H
#define	HELMHOLTZNEUMANNROBINOPERATOR_H

#include <vector>

#include "LinearOperator.h"
#include "Vector.h"
#include "BESpace.h"
#include "SparseMatrix.h"
#include "IterativeSolver.h"
#include "RepresentationFormulaHelmholtz.h"

namespace bem4i {

/*!
 * Class representing a weak operator on the given discretized space
 */
template<class LO, class SC>
class HelmholtzNeumannRobinOperator :
public LinearOperator<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! todo: assembly without V, K, D, M
  /*!
   * @brief Constructor with SurfaceMesh3D
   *
   * @param mesh   surface mesh
   * @param inE    list of input elems
   * @param outE   list of output elems
   * @param kappa  kappa
   * @param V      user supplied Helmholtz single-layer operator
   * @param K      user supplied Helmholtz double-layer operator
   * @param D      user supplied Helmholtz hypersingular operator
   * @param M      user supplied identity operator
   */
  HelmholtzNeumannRobinOperator(
      SurfaceMesh3D< LO, SC > * mesh,
      std::vector< LO > & inE,
      std::vector< LO > & outE,
      SC kappa,
      LinearOperator< LO, SC > * V,
      LinearOperator< LO, SC > * K,
      LinearOperator< LO, SC > * D,
      LinearOperator< LO, SC > * M,
      basisType type
      );

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

  bool LUSolve(
      Vector<LO, SC> & rhs
      );

  void getRHS(
      Vector< LO, SC > & rhs
      ) const;

  void solveDirichletProblem(
      const Vector< LO, SC > & u,
      Vector< LO, SC > & dudn
      ) const;
  
  // todo: grid, u and dudn should be const
  void evalInside(
      SurfaceMesh3D< LO, SC > & grid, 
      Vector< LO, SC > & u,
      Vector< LO, SC > & dudn,
      Vector< LO, SC > & res
      ) const;

  //! Destructor
  virtual ~HelmholtzNeumannRobinOperator( );


private:

  basisType type;
  
  //! robin condition parameter
  SC ubar;
  
  //! Surface mesh
  SurfaceMesh3D< LO, SC > * mesh;

  //! interior wave number
  SC kappa;
  
  //! input elems
  std::vector< LO > inE;
  
  //! output elems
  std::vector< LO > outE;

  //! single-layer operator
  LinearOperator< LO, SC > * V;

  //! double-layer operator
  LinearOperator< LO, SC > * K;

  //! hypersingular operator
  LinearOperator< LO, SC > * D;

  //! Identity operator
  LinearOperator< LO, SC > * M;

  SparseMatrix< LO, SC > * MTilde;
  
  void assembleMTilde( );
  
  void assembleMTildeP1P1( );
  
  void assembleMTildeP0P0( );

  //! Default constructor (empty)
  HelmholtzNeumannRobinOperator( );

  //! Copy constructor (empty)
  HelmholtzNeumannRobinOperator(
      const HelmholtzNeumannRobinOperator & orig
      );

};

}
// include .cpp file to overcome linking problems due to templates
#include "HelmholtzNeumannRobinOperator.cpp"

#endif	/* HELMHOLTZNEUMANNROBINOPERATOR_H */
