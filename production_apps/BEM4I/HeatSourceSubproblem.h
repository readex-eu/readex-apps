/*!
 * @file    HeatSourceSubproblem.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    September 1, 2014
 * @brief   Header file for class HeatSourceSubproblem
 * 
 */

#ifndef HEATSOURCESUBPROBLEM_H
#define	HEATSOURCESUBPROBLEM_H

#include "SurfaceMesh3D.h"
#include "Vector.h"
#include "OptimizationSubproblem.h"
#include "BEBilinearFormLaplace1Layer.h"
#include "BEBilinearFormLaplace2Layer.h"
#include "RepresentationFormulaLaplace.h"

namespace bem4i {

template< class LO, class SC >
class HeatSourceSubproblem : public OptimizationSubproblem< LO, SC > {

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  /*!
   * Constructs HeatSourceSubproblem from the fixed and free meshes,
   * given Dirichlet data, and the fitting constant Q.
   * 
   * @param[in] targetNeumannDataDetector Neumann data on the sensor
   */
  HeatSourceSubproblem(
      const Vector< LO, SC > & targetNeumannDataSensor
      );

  //! Destructor
  virtual ~HeatSourceSubproblem( );

  //! Prints basic info about this.
  virtual void printInfo( ) const;

  //! Prints available data to Paraview vtu file.
  virtual void printVtu(
      const string & filename
      ) const;

  void setCostMultiplier(
      SCVT costMultiplier
      ) {
    this->costMultiplier = costMultiplier;
  }

  void setRegularizationParameter(
      SCVT regularizationPar
      ) {
    this->regularizationPar = regularizationPar;
  }

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
      const SurfaceMesh3D< LO, SC > * meshSource,
      const SurfaceMesh3D< LO, SC > * meshSensor
      );

protected:

private:

  //! Default constructor
  HeatSourceSubproblem( );

  void elem2NodalAreaWeighted(
      const Vector< LO, SC > & elem,
      Vector< LO, SC > & nodal
      ) const;

  void setUpPrimalRHS(
      Vector< LO, SC > & rhs
      ) const;

  SC evalNewtonKernel(
      const SCVT * x,
      const SCVT * y,
      const SCVT * n )
  const;

  void getQuadratureNodes(
      const SCVT * x1,
      const SCVT * x2,
      const SCVT * x3,
      int quadratureOrder,
      SCVT * nodes )
  const;


  Vector< LO, SC > * targetNeumannDataSensor;

  Vector< LO, SC > * primalDirichletData;

  Vector< LO, SC > * primalNeumannData;

  Vector< LO, SC > * adjointDirichletData;

  Vector< LO, SC > * adjointNeumannData;

  Vector< LO, SC > * adjointOnSource;

  SurfaceMesh3D< LO, SC > * meshSensor;

  SurfaceMesh3D< LO, SC > * meshSource;

  LO nElemsFree;

  LO nElemsFixed;

  LO nNodesFree;

  LO nNodesFixed;

  SCVT cost;

  SCVT costMultiplier;

  SCVT regularizationPar;

};

}

// include .cpp file to overcome linking problems due to templates
#include "HeatSourceSubproblem.cpp"

#endif	/* HEATSOURCESUBPROBLEM_H */

