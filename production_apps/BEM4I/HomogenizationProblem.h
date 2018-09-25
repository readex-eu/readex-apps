/*!
 * @file    HomogenizationProblem.h 
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    December 2, 2014
 * @brief   Header file for class HomogenizationProblem
 * 
 */

#ifndef HOMOGENIZATIONPROBLEM_H
#define	HOMOGENIZATIONPROBLEM_H

#include "Problem.h"
#include "BEBilinearFormLaplace1Layer.h"
#include "BEBilinearFormLaplace2Layer.h"
#include "BEBilinearFormLaplaceHypersingular.h"

namespace bem4i {

/*! 
 * Class representing direct homoganization problem
 * 
 */
template<class LO, class SC>
class HomogenizationProblem : public Problem<LO, SC> {

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  enum HomogenizationSettings {

    MESH_INCLUSION, //!< mesh of the incluson (pointer to SurfaceMesh3D)
    N_CELL_REFINE, //!< number of refines for the basic cell mesh  
    SYSTEM_TYPE, //!< system matrix type

    SETTINGS_SIZE //!< number of settings, always the last in this enum!!!
  };

  typedef typename HomogenizationProblem< LO, SC >::HomogenizationSettings
  Settings;

  enum SystemType {

    FULL_STEKLOV_POINCARE, //!< Steklov-Poincare formulation with full matrices

    SYSTEM_TYPE_SIZE //!< number of settings, always the last in this enum!!!     
  };

  HomogenizationProblem( );

  //! destructor
  virtual ~HomogenizationProblem( );

  virtual bool solve( );

  //  /*!
  //   * Sets parameters of a solver
  //   */
  //  template< class argType >
  //  void set(
  //      Settings name,
  //      argType & value
  //      );

  void getHomoMatrix(
      FullMatrix<LO, SC> & homoMatrix
      ) const {
    if ( this->homoMatrix ) {
      this->homoMatrix->copy( homoMatrix );
    }
  }

  void setMeshInclusion(
      SurfaceMesh3D< LO, SC > & meshInclusion
      ) {
    this->meshInclusion = &meshInclusion;
  }

  int getNSegsPerEdge( ) const {
    return nSegsPerEdge;
  }

  void setNSegsPerEdge(
      unsigned nSegsPerEdge
      ) {
    this->nSegsPerEdge = nSegsPerEdge;
  }

  //! prints available mesh data
  void printVtu(
      const std::string & fileName
      );

protected:

  //! surface mesh pointer
  SurfaceMesh3D<LO, SC> * meshCell;

  //! evaluation mesh pointer
  SurfaceMesh3D<LO, SC> * meshInclusion;

  //! concatenated mesh
  SurfaceMesh3D<LO, SC> * mesh;

  //! matrix ensuring periodicity
  SparseMatrix<LO, SC> * P;

  //! homogenized coefficients
  FullMatrix<LO, SC> * homoMatrix;

  //! auxiliary function on G1
  FullMatrix<LO, SC> * chi1;

  //! auxiliary function on G2
  FullMatrix<LO, SC> * chi2;

  //! Type of the system to be solved
  SystemType systemType;

  //! number of segments per cube edge
  int nSegsPerEdge;

  //! coefficient for inclusion
  SCVT a1;

  //! coefficient for matrix
  SCVT a2;

  //! order of quadrature
  int quadOrder;

  //! type of quadrature
  quadratureType quadType;

  //! array of booleans indicating set parameters
  bool parameters[ SETTINGS_SIZE ];

  HomogenizationProblem(
      const HomogenizationProblem & orig
      );

  //! sets default values
  void setDefaultParameters( );

  //! checks if all necessary parameters are set (in this case, return true)
  bool checkParameters( ) const;

  bool solveFullSteklovPoincare( );

  void assembleHomoMatrix( );

  void prepareMesh( );

  void assembleIdentity(
      SparseMatrix<LO, SC> & M,
      Vector<LO, SC> & a
      );

  void assembleIdentityG1G1(
      SparseMatrix<LO, SC> & MG1G1
      );

  void copyIntoSubmatrix(
      const FullMatrix< LO, SC > & M,
      LO rowsStart,
      LO rowsEnd,
      LO colsStart,
      LO colsEnd,
      FullMatrix< LO, SC > & MG1G1
      );

  void copyFromSubmatrix(
      FullMatrix< LO, SC > & M,
      LO rowsStart,
      LO colsStart,
      const FullMatrix< LO, SC > & MG1G1
      );

  void addLeftUpper(
      const FullMatrix< LO, SC > & M,
      LO nRows,
      LO nCols,
      FullMatrix< LO, SC > & MG1G1
      );

  void prepareFullSteklovPoincareMatrix(
      FullMatrix<LO, SC> & A
      );

  void prepareFullSteklovPoincareRHS(
      FullMatrix<LO, SC> & rhs
      );

};

} // end of namespace bem4i

// include .cpp file to overcome linking problems due to templates
#include "HomogenizationProblem.cpp"

#endif	/* HOMOGENIZATIONPROBLEM_H */

