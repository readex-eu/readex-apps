/*!
 * @file    OptimizationSubproblem.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    May 23, 2014
 * @brief   Header file for class OptimizationSubproblem
 * 
 */

#ifndef OPTIMIZATIONSUBPROBLEM_H
#define	OPTIMIZATIONSUBPROBLEM_H

#include "Vector.h"
#include "Problem.h"
#include "SurfaceMesh3D.h"

namespace bem4i {

template< class LO, class SC >
class OptimizationSubproblem : public Problem< LO, SC > {

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  virtual ~OptimizationSubproblem( ) {
  };

  /*!
   * Solves the boundary value subproblem and fills in data necessary to 
   * compute the shape gradient
   */
  virtual bool solve( ) = 0;

  /*!
   * returns the value of the cost functional (need to solve first)
   */
  virtual void getCost(
      SCVT & cost
      ) const = 0;

  /*!
   * Returns the shape gradient corresponding to the perturbation given
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
      ) const = 0;

  /*!
   * Returns the shape gradient corresponding to the perturbation given
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
      ) const = 0;

  /*!
   * Returns the shape gradient kernel in each node of the free mesh
   *
   * @param[in,out] grad shape gradient kernel
   */
  virtual void getShapeGradient(
      Vector< LO, SC > & grad
      ) const = 0;

  /*!
   * Updates current mesh and mesh data
   *
   * @param[in] freeMesh free mesh
   * @param[in] fixedMesh fixed mesh
   */
  virtual void setProblemData(
      const SurfaceMesh3D< LO, SC > * freeMesh,
      const SurfaceMesh3D< LO, SC > * fixedMesh
      ) = 0;


  //! Prints basic info about this.
  virtual void printInfo( ) const = 0;

  //! Prints available data to Paraview vtu file.
  virtual void printVtu(
      const string & filename
      ) const = 0;

};

}

#endif	/* OPTIMIZATIONSUBPROBLEM_H */
