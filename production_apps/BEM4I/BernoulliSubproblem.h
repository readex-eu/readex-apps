/*!
 * @file    BernoulliSubproblem.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    May 23, 2014
 * @brief   Header file for class BernoulliSubproblem
 * 
 */

#ifndef BERNOULLISUBPROBLEM_H
#define	BERNOULLISUBPROBLEM_H

#include "SurfaceMesh3D.h"
#include "Vector.h"
#include "OptimizationSubproblem.h"
#include "BEBilinearFormLaplace1Layer.h"
#include "BEBilinearFormLaplace2Layer.h"
#include "IdentityOperator.h"

namespace bem4i {

enum costType {
  L2Norm,
  VNorm
};

template< class LO, class SC >
class BernoulliSubproblem : public OptimizationSubproblem< LO, SC > {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  /*!
   * Constructs BernoulliSubproblem from the fixed and free meshes,
   * given Dirichlet data, and the fitting constant Q.
   * 
   * @param[in] dirichletDataFree Dirichlet data for the free component
   * @param[in] dirichletDataFixed Dirichlet data for the fixed component
   * @param[in] Q fitting constant
   */
  BernoulliSubproblem(
      SC dirichletDataFree,
      const Vector< LO, SC > & dirichletDataFixed,
      SC Q,
      basisType dirBasis,
      basisType neuBasis
      );

  //! Destructor
  virtual ~BernoulliSubproblem( );

  //! Prints basic info about this.
  virtual void printInfo( ) const;

  //! Prints available data to Paraview vtu file.
  virtual void printVtu(
      const string & filename
      ) const;

  /*!
   * Solves the boundary value subproblem and fills in data necessary to 
   * compute the shape gradient.
   */
  virtual bool solve( );

  /*!
   * Returns the value of the cost functional (need to solve first).
   */
  virtual void getCost(
      SCVT & cost
      ) const;

  /*!
   * Returns the shape gradient corresponding to the perturbation given.
   * 
   * @param[in] perturbation perturbation vector on the current mesh
   * @param[in,out] dx1 shape gradient for perturbation in x1
   * @param[in,out] dx2 shape gradient for perturbation in x2
   * @param[in,out] dx3 shape gradient for perturbation in x3
   */
  virtual void getShapeGradient(
      const Vector< LO, SCVT > & perturbation,
      SCVT & dx1,
      SCVT & dx2,
      SCVT & dx3
      ) const;

  /*!
   * Returns the shape gradient corresponding to the perturbation given.
   * 
   * @param[in] perturbationX1 perturbation vector on the current mesh in x1
   * @param[in] perturbationX2 perturbation vector on the current mesh in x2
   * @param[in] perturbationX3 perturbation vector on the current mesh in x3
   * @param[in,out] dx shape gradient for perturbation in x3
   */
  virtual void getShapeGradient(
      const Vector< LO, SCVT > & perturbationX1,
      const Vector< LO, SCVT > & perturbationX2,
      const Vector< LO, SCVT > & perturbationX3,
      SCVT & dx
      ) const;

  /*!
   * Returns the shape gradient kernel in each node of the free mesh
   *
   * @param[in,out] grad shape gradient kernel
   */
  virtual void getShapeGradient(
      Vector< LO, SC > & grad
      ) const;

  /*!
   * Updates current mesh and mesh data
   *
   * @param[in] freeMesh free mesh
   * @param[in] fixedMesh fixed mesh
   */
  virtual void setProblemData(
      const SurfaceMesh3D< LO, SC > * freeMesh,
      const SurfaceMesh3D< LO, SC > * fixedMesh
      );

  void setCostType(
      costType cost
      ) {
    this->myCost = cost;
  }

protected:

private:

  //! Default constructor
  BernoulliSubproblem( );

  void elem2NodalAreaWeighted(
      const Vector< LO, SC > & elem,
      Vector< LO, SC > & nodal
      ) const;

  void updateCost( );

  void getShapeGradientL2Norm(
      const Vector< LO, SCVT > & perturbationX1,
      const Vector< LO, SCVT > & perturbationX2,
      const Vector< LO, SCVT > & perturbationX3,
      SCVT & dx
      ) const;

  void getShapeGradientVNorm(
      const Vector< LO, SCVT > & perturbationX1,
      const Vector< LO, SCVT > & perturbationX2,
      const Vector< LO, SCVT > & perturbationX3,
      SCVT & dx
      ) const;

  SC dirichletDataFree;

  Vector< LO, SC > * dirichletDataFixed;

  Vector< LO, SC > * primalDirichletData;

  Vector< LO, SC > * primalNeumannData;

  Vector< LO, SC > * dualDirichletData;

  Vector< LO, SC > * dualNeumannData;

  FullMatrix< LO, SC > * V;

  FullMatrix< LO, SC > * K;
  
  IdentityOperator< LO, SC > * id;

  SurfaceMesh3D< LO, SC > * mesh;

  LO nElemsFree;

  LO nElemsFixed;

  LO nNodesFree;

  LO nNodesFixed;

  SCVT Q;

  SCVT cost;

  basisType dirBasis;

  basisType neuBasis;

  costType myCost;

  inline double timeDiff( timeval start, timeval end ) {
    return (( end.tv_sec - start.tv_sec ) * 1000000u +
        end.tv_usec - start.tv_usec ) / 1.e6;
  }

};

}

// include .cpp file to overcome linking problems due to templates
#include "BernoulliSubproblem.cpp"

#endif	/* BERNOULLISUBPROBLEM_H */

