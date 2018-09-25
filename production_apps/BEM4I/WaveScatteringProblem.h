/*!
 * @file    WaveScatteringProblem.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    June 10, 2014
 * @brief   Header file for class WaveScatteringProblem
 * 
 */

#ifndef WAVESCATTERINGPROBLEM_H
#define	WAVESCATTERINGPROBLEM_H

#include "Problem.h"
#include "Macros.h"
#include "BEIntegratorWave.h"
#include "BEBilinearFormWave1Layer.h"
#include "BEBilinearFormWaveHypersingular.h"
#include <algorithm>


namespace bem4i {

enum WaveSettings {

  PROBLEM_TYPE, //!< type of problem (Dirichlet, Neumann, Mixed) (enum)
  INPUT_MESH_FILE, //!< string with path to the mesh file (string)
  INPUT_MESH, //!< input mesh (pointer to mesh)
  N_REFINE_INPUT, //!< number of refines of the input mesh (int)
  SCALE_INPUT, //!< factor to scale input mesh (SCVT)
  END_TIME, //!< end time (SCVT)
  INCIDENT_WAVE, //!< incident wave (functio pointer)
  INCIDENT_WAVE_DU_DN, //!< du/dn of inc. wave on boundary (for neu., f. point.)
  N_TIME_STEPS, //!< number of time-steps (int)
  LEGENDRE_ORDER, //!< order of Legendre polynomial (int)
  SPACE_QUAD_ORDER, //!< order of spatial quadrature (pointer to array of 4 int)
  TEMP_QUAD_ORDER, //!< order of temporal quadrature (int)
  SOLVER_TYPE, //!< type of solver (enum)
  ITERATIVE_SOLVER, //!< type of iterative solver (GMRES, DGMRES)
  SOLVER_PRECISION, //!< precision of iterative solver (SCVT)
  SOLVER_MAX_IT, //!< maximum number of iteration of iterative solver (int)
  SOLVER_RESTARTS, //!< restarts of GMRES iterative solver (int)
  DGMRES_MAX_DIM, //!< maximum dimension of an invariant subspace (int)
  DGMRES_MAX_DIM_IT, //!< max. n. of vectors to add to subspace per iter. (int)
  EVAL_MESH_FILE, //!< mesh file for evaluating result (string)
  EVAL_MESH, //!< evaluation mesh (pointer to mesh)
  N_REFINE_EVAL, //!< number of refines of the evaluation mesh (LO)
  SCALE_EVAL, //!< factor to scale the evaluation mesh (SCVT)
  OUTPUT_FILE //!< name base of the output vtu file (string)
};

enum ProblemType {

  DIRICHLET, //!< Dirichlet problem
  NEUMANN, //!< Neumann problem
  MIXED //!< Dirichlet and Neumann problem
};

enum SolverType {
  DIRECT, //!< direct solver (based on Stewart)
  ITERATIVE //!< GMRES
};

enum IterativeSolverType {
  GMRES, //!< GMRES
  DGMRES //!< GMRES with deflations
};

/*! 
 * Class representing problem of time dependent scattering
 * 
 */
template<class LO, class SC>
class WaveScatteringProblem : public Problem<LO, SC> {

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:
  //! default constructor
  WaveScatteringProblem( );

  //! copy constructor
  WaveScatteringProblem( const WaveScatteringProblem& orig );

  //! destructor
  virtual ~WaveScatteringProblem( );

  virtual bool solve( );

  /*!
   * Sets parameters of a solver
   */
  virtual void set(
      int name,
      SCVT value
      );

  /*!
   * Sets parameters of a solver
   */
  virtual void set(
      int name,
      const string & value
      );

  /*!
   * Sets parameters of a solver
   */
  virtual void set(
      int name,
      void * value
      );

  /*!
   * Sets parameters of a solver
   */
  virtual void set(
      int name,
      int value
      );

  /*!
   * Sets parameter of a solver
   */
  virtual void set(
      int name,
      SC( *incWave )( SCVT t, SCVT* x, SCVT *n )
      );


  /*!
   * Sets parameters of a solver
   */
  virtual void set(
      int name,
      SC( *incWave )( SCVT t, SCVT* x )
      );

  /*!
   * Sets an incident wave
   * 
   * @param[in]   *incWave pointer to the function representing incident wave
   *               (time, point, normal to the surface in point)
   */
  void setIncidentWaveDuDn( SC( *incWave )( SCVT t, SCVT* x, SCVT *n ) ) {
    this->incWaveDuDn = incWave;
    parameters[INCIDENT_WAVE_DU_DN] = true;
  }

  /*!
   * Sets an incident wave
   * 
   * @param[in]   *infWave pointer to the function representing incident wave
   *               (time, point, normal to the surface in point)
   */
  void setIncidentWave( SC( *incWave )( SCVT t, SCVT* x ) ) {
    this->incWave = incWave;
    parameters[INCIDENT_WAVE] = true;
  }

protected:

  //! sets default values
  void setDefaultParameters( );

  //! checks if all necessary parameters are set (in this case, return true)
  bool checkParameters( ) const;

  //! creates system matrix and RHS for the Dirichlet problem
  void prepareDirichletSystem(
      MPIBlockMatrix<LO, SC> &A,
      Vector<LO, SC> &rhs
      );

  //! creates system matrix and RHS for the Neumann problem
  void prepareNeumannSystem(
      MPIBlockMatrix<LO, SC> &A,
      Vector<LO, SC> &rhs
      );

  //! solves the problem using iterative solver
  bool solveIteratively(
      MPIBlockMatrix<LO, SC> &matrix,
      Vector<LO, SC> &rhs,
      Vector<LO, SC> &solution
      );

  //! solves the problem using direct solver
  bool solveDirectly(
      MPIBlockMatrix<LO, SC> &matrix,
      Vector<LO, SC> &rhs,
      Vector<LO, SC> &solution
      );

  //! saves to vtu files
  void saveVtu( Vector<LO, SC> &solution );


  //! type of problem (Dirichlet, Neumann, Mixed)
  ProblemType problemType;

  //! type of solver (direct, iterative)
  SolverType solverType;
  
  //! type of iterative solver (GMRES, DGMRES)
  IterativeSolverType iterativeSolver;

  //! number of refines of input surface mesh
  int nInputRefines;

  //! number of refines of evaluation mesh
  int nEvalRefines;

  //! number of time-steps
  int nTimeSteps;

  //! order of Legendre polynomial
  int legendreOrder;

  //! order of temporal quadrature
  int tempQuadOrder;

  //! order of spatial quadrature
  int *spatialQuad;

  //! path to the surface mesh
  std::string inputMeshFile;

  //! path to the evaluation mesh
  std::string evalMeshFile;

  //! path to the evaluation mesh
  std::string outputFile;

  //! function pointer to incident wave function (time, space, normal)
  SC( *incWaveDuDn )( SCVT, SCVT*, SCVT* );

  //! function pointer to du/dn of incident wave function (time, space)
  SC( *incWave )( SCVT, SCVT* );

  //! end time
  SCVT endTime;

  //! iterative solver precision
  SCVT solverPrecision;

  //! scale the evaluation mesh by given factor
  SCVT evalScaleFactor;
  
  //! scale the input mesh by given factor
  SCVT meshScaleFactor;

  //! maximum number of iterations
  int maxIt;

  //! restarts of GMRES
  int gmresRestarts;
  
  //! max. dimension of invariant subspace for DGMRES
  int dgmresMaxDim;

  //! max. eigs to add to subspace per iteration (DGMRES)
  int dgmresMaxDimIt;
  
  //! surface mesh pointer
  SurfaceMesh3D<LO, SC>* mesh;

  //! evaluation mesh pointer
  SurfaceMesh3D<LO, SC>* evalMesh;

  //! array of booleans indicating set parameters
  bool* parameters;
};

} // end of namespace bem4i

// include .cpp file to overcome linking problems due to templates
#include "WaveScatteringProblem.cpp"

#endif	/* WAVESCATTERINGPROBLEM_H */

