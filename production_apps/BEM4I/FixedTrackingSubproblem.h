/*!
 * @file    FixedTrackingSubproblem.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    August 12, 2014
 * @brief   Header file for class FixedTrackingSubproblem
 * 
 */

#ifndef FIXEDTRACKINGSUBPROBLEM_H
#define	FIXEDTRACKINGSUBPROBLEM_H

#include "SurfaceMesh3D.h"
#include "Vector.h"
#include "OptimizationSubproblem.h"

namespace bem4i {

template< class LO, class SC >
class FixedTrackingSubproblem : public OptimizationSubproblem< LO, SC > {

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  /*!
   * Constructs FixedTrackingSubproblem from the fixed and free meshes,
   * given Dirichlet data, and the fitting constant Q.
   * 
   * @param[in] dirichletDataFree Dirichlet data for the free component
   * @param[in] dirichletDataFixed Dirichlet data for the fixed component
   * @param[in] Q fitting constant
   */
  FixedTrackingSubproblem(
      SC dirichletDataFree,
      SC dirichletDataFixed,
      const Vector< LO, SC > & targetNeumannDataFixed
      );

  //! Destructor
  virtual ~FixedTrackingSubproblem( );

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

protected:

private:

  //! Default constructor
  FixedTrackingSubproblem( );

  void elem2NodalAreaWeighted(
      const Vector< LO, SC > & elem,
      Vector< LO, SC > & nodal
      ) const;

  SC dirichletDataFree;

  SC dirichletDataFixed;
  
  Vector< LO, SC > * targetNeumannDataFixed;

  Vector< LO, SC > * primalDirichletData;

  Vector< LO, SC > * primalNeumannData;

  Vector< LO, SC > * dualDirichletData;

  Vector< LO, SC > * dualNeumannData;

  SurfaceMesh3D< LO, SC > * mesh;

  LO nElemsFree;

  LO nElemsFixed;

  LO nNodesFree;

  LO nNodesFixed;

  SCVT cost;

};

}

// include .cpp file to overcome linking problems due to templates
#include "FixedTrackingSubproblem.cpp"

#endif	/* FIXEDTRACKINGSUBPROBLEM_H */

