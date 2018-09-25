/*!
 * @file    MultiresolutionOptimizer.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    May 19, 2014
 * @brief   Header file for class MultiresolutionOptimizer
 * 
 */

#ifndef MULTIRESOLUTIONOPTIMIZER_H
#define	MULTIRESOLUTIONOPTIMIZER_H

#include <string>
#include <vector>

#include "SurfaceMesh3D.h"
#include "OpenMeshWrapper.h"
#include "OptimizationSubproblem.h"
#include "BernoulliSubproblem.h"
#include "IOHelper.h"

#include "nlopt.hpp"

#ifdef IPOPT
#include "IpTNLP.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpIpoptApplication.hpp"
#endif

namespace bem4i {

enum MOType {

  freeForm,
  fixedRef
};

/*! 
 * Class for multiresolution shape optimization
 * 
 */
template<class LO, class SC>
#ifdef IPOPT
class MultiresolutionOptimizer : public Ipopt::TNLP {

#else
class MultiresolutionOptimizer {

#endif

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! default constructor
  MultiresolutionOptimizer( MOType type = freeForm );

  //  /*!
  //   * Constructs MultiresolutionOptimizer from SurfaceMesh3D (meshes are copied!)
  //   * 
  //   * @param[in] fixedPart fixed mesh
  //   * @param[in] freePart (coarse) free mesh 
  //   * @param[in] analysisLevel subdivision level for analysis 
  //   * @param[in] maxOptimLevel maximal optimization level 
  //   */
  //  MultiresolutionOptimizer(
  //      const SurfaceMesh3D< LO, SC > * fixedPart,
  //      const SurfaceMesh3D< LO, SC > * freePart,
  //      int analysisLevel,
  //      int maxOptimLevel,
  //      OptimizationSubproblem< LO, SC > & bvp,
  //      const string & outputPath,
  //      const string & tagFile = "",
  //      SCVT nloptFRel = 1e-2,
  //      int jumpTol = 1
  //      );

  //! destructor
  ~MultiresolutionOptimizer( );

  void setFixedPart(
      const SurfaceMesh3D< LO, SC > & fixedPart
      ) {
    if ( this->fixedPart ) delete this->fixedPart;
    this->fixedPart = new SurfaceMesh3D< LO, SC >( fixedPart );
  }

  void setFreePart(
      const SurfaceMesh3D< LO, SC > & freePart,
      int analysisLevel
      ) {
    this->setFreePart( freePart, analysisLevel, "", "" );
  }

  void setFreePart(
      const SurfaceMesh3D< LO, SC > & freePart,
      int analysisLevel,
      const std::string & tagFile,
      const std::string & bndFile
      ) {

    if ( this->freePartOptim ) delete this->freePartOptim;
    this->freePartOptim = new OpenMeshWrapper< LO, SC >( freePart );
    this->analysisLevel = analysisLevel;
    this->setFreePartAnalAndInitMatrices( );

    if ( !tagFile.empty( ) ) {
      this->freePartOptim->readTags( tagFile );
    }

    this->resetBounds( );

    if ( !bndFile.empty( ) ) {
      this->readRelativeBounds( bndFile );
    }

    this->alphaOpt.clear( );
    this->alphaOpt.resize( this->freePartOptim->getNNodes( ), 0.0 );

    SCVT node[ 3 ];
    this->xOpt.clear( );
    for ( LO i = 0; i < this->freePartOptim->getNNodes( ); ++i ) {
      this->freePartOptim->getNode( i, node );
      this->xOpt.push_back( node[ 0 ] );
      this->xOpt.push_back( node[ 1 ] );
      this->xOpt.push_back( node[ 2 ] );
    }

    // save current reference domain
    if ( this->freePartOptimRef ) delete this->freePartOptimRef;
    this->freePartOptimRef =
        new OpenMeshWrapper< LO, SC >( *this->freePartOptim );
  }

  void setMaxOptimLevel(
      int level
      ) {
    this->maxOptimLevel = level;
  }

  void setBVP(
      OptimizationSubproblem< LO, SC > & bvp
      ) {
    this->bvp = &bvp;
  }

  void setNloptFRel(
      SCVT val
      ) {
    this->nloptFRel = val;
  }

  void setNloptJumpTol(
      int tol
      ) {
    this->nloptJumpTol = tol;
  }

  void setOutputPath(
      const std::string & outputPath
      ) {
    const std::string createDir = "mkdir -p " + outputPath;
    if ( system( createDir.c_str( ) ) ) exit( -1 );
    this->outputPath = outputPath + "/";
  }

  /*!
   * Propagates data from coarse to fine level using subdivision matrices
   * 
   * @param[in] coarse data defined on the coarse level
   * @param[in,out] fine data defined on the fine level (to be computed)
   */
  void coarseToFinePropagate(
      const Vector< LO, SC > & coarse,
      Vector< LO, SC > & fine
      ) const;

  /*!
   * Propagates data from fine to coarse level using least square
   * fitting Eigen (QR decomposition) with the subdivision matrix
   * 
   * @param[in] fine data defined on the fine level
   * @param[in,out] coarse data defined on the coarse level (to be computed)
   */
  void fineToCoarsePropagate(
      const Vector< LO, SC > & fine,
      Vector< LO, SC > & coarse
      );

  //! returns number of nodes of the free part on current optimization level

  inline LO getFreePartOptimNNodes( ) const {
    return this->freePartOptim->getNNodes( );
  }

  //! returns number of nodes of the free part on analysis level

  inline LO getFreePartAnalNNodes( ) const {
    return this->freePartAnal->getNNodes( );
  }

  //! returns number of nodes of the fixed part

  inline LO getFixedPartNNodes( ) const {
    return this->fixedPart->getNNodes( );
  }

  //! sets normalizeCost flag

  void setNormalizeCost(
      bool normalizeCost
      ) {
    this->normalizeCost = normalizeCost;
  }

  //! sets normalizeGrad flag

  void setNormalizeGrad(
      bool normalizeGrad
      ) {
    this->normalizeGrad = normalizeGrad;
  }

  //! sets normalizeOnEachLevel flag

  void setNormalizeOnEachLevel(
      bool normalizeOnEachLevel
      ) {
    this->normalizeOnEachLevel = normalizeOnEachLevel;
  }

  //! prints basic info about the object
  void printInfo( ) const;

  /*!
   * Prints fixed mesh to the xml paraview file format
   * 
   * @param[in] filename string with the target file name
   */
  void printFixedPartVtu(
      const string & filename
      ) const;

  /*!
   * Prints free mesh on analysis level to the xml paraview file format
   * 
   * @param[in] filename string with the target file name
   */
  void printFreePartAnalVtu(
      const string & filename
      ) const;

  /*!
   * Prints free mesh on analysis level to the xml paraview file format
   * 
   * @param[in] meshFile string with the target file name
   * @param[in] nodeNames vector of names for nodal data
   * @param[in] nodalData vector of data defined nodewise
   * @param[in] elemNames vector of names for element data
   * @param[in] elemData vector of data defined elementwise
   */
  void printFreePartAnalVtu(
      const string & meshFile,
      const std::vector< string > * nodeNames,
      const std::vector< Vector< LO, SC >* > * nodalData,
      const std::vector< string > * elemNames,
      const std::vector< Vector< LO, SC >* > * elemData
      ) const;

  /*!
   * Prints free mesh on analysis level to the xml paraview file format
   * 
   * @param[in] meshFile string with the target file name
   * @param[in] nodeNames vector of names for nodal data
   * @param[in] nodalData vector of data defined nodewise
   * @param[in] elemNames vector of names for element data
   * @param[in] elemData vector of data defined elementwise
   */
  void printFreePartOptimVtu(
      const string & meshFile,
      const std::vector< string > * nodeNames,
      const std::vector< Vector< LO, SC >* > * nodalData,
      const std::vector< string > * elemNames,
      const std::vector< Vector< LO, SC >* > * elemData
      ) const;

  /*!
   * Prints mesh to the xml paraview file format
   * 
   * @param[in] filename string with the target file name
   */
  void printFreePartOptimVtu(
      const string & filename
      ) const;

  /*!
   * Performs optimization on current level
   */
  void optimize( );

  /*!
   * Cost function for NLopt
   * 
   * @param[in] x current point
   * @param[in,out] grad current gradient
   * @param[in,out] data user data
   */
  static double nloptCostFunctionFreeForm(
      const std::vector< double > & x,
      std::vector< double > & grad,
      void * data
      );

  /*!
   * Cost function for NLopt
   * 
   * @param[in] alpha current coefficients
   * @param[in,out] grad current gradient
   * @param[in,out] data user data
   */
  static double nloptCostFunctionFixedReference(
      const std::vector< double > & alpha,
      std::vector< double > & grad,
      void * data
      );

protected:

private:

  /*!
   * Sets relative bounds on optimization parameters
   * 
   * @param[in] boundFile file describing the lower and upper bounds
   */
  bool readRelativeBounds(
      const std::string & boundsFile
      );

  /*!
   * Sets relative bounds on optimization parameters
   * 
   * @param[in] bound relative bound for all parameters
   */
  void setRelativeBounds(
      SCVT bound
      );


  //! resets bounds to -inf, +inf
  void resetBounds( );

  /*!
   * Sets absolute bounds ([x - lb, x + ub])
   * 
   * @param[in,out] lb lower bound
   * @param[in,out] ub upper bound
   */
  void setUpAbsoluteBounds(
      std::vector< double > & lb,
      std::vector< double > & ub
      ) const;

  //! sets optimization bounds for next optimization level
  // needs to be called before the subdivision of freePartOptim!!!
  void setBoundsForNextLevel( );

  /*!
   * Computes gradient for free movement of free nodes
   * 
   * @param[in,out] grad output gradient
   */
  void getGradientFree(
      Vector< LO, SC > & grad
      ) const;

  /*!
   * Computes gradient for free normal movement of free nodes
   * 
   * @param[in,out] grad output gradient
   */
  void getGradientNormal(
      Vector< LO, SC > & grad
      ) const;

  /*!
   * Computes gradient for normal movement of free nodes
   * 
   * @param[in,out] grad output gradient
   */
  void getGradientNormalFixedReference(
      Vector< LO, SC > & grad
      ) const;

  /*!
   * Computes gradient from the fine kernel
   * 
   * @param[in,out] grad output gradient
   */
  void getGradientKernel(
      Vector< LO, SC > & grad
      );

  //  /*!
  //   * Get lower and upper bounds on coarse mesh coordinates
  //   * 
  //   * @param[in,out] lb lower bounds
  //   * @param[in,out] ub upper bounds
  //   */
  //  void getBoundsOnCoarseMesh( );

  /*!
   * Performs free form optimization on current level
   * every iteration changes reference domain to the latest configuration
   */
  void optimizeFreeForm( );

  /*!
   * Performs free form optimization on current level
   * every iteration is obtained by a perturbation of original domain
   */
  void optimizeFixedReference( );

  /*!
   * Increases the subdivision level if the optimized mesh
   */
  void increaseOptimizationLevel( );

  void increaseOptimizationLevelFreeForm( );

  void increaseOptimizationLevelFixedRef( );

  void updateFreePartAnal( );

  void setFreePartAnalAndInitMatrices( );

  void computeSubdivMatrix( );

  void convergenceMonitor( );

  void updateFreePartOptimFromReference(
      const std::vector< double > & alpha
      );

  int currOptimLevel;

  int analysisLevel;

  int maxOptimLevel;

  std::string outputPath;

  SurfaceMesh3D< LO, SC > * fixedPart;

  SurfaceMesh3D< LO, SC > * freePartAnal;

  OpenMeshWrapper< LO, SC > * freePartOptim;

  OpenMeshWrapper< LO, SC > * freePartOptimRef;

  std::vector< SparseMatrix< LO, SCVT > * > subdivMatrices;

  SparseMatrix< LO, SCVT > * subdivMatrix;

  OptimizationSubproblem< LO, SC > * bvp;

  std::vector< double > relativeLowerBounds;

  std::vector< double > relativeUpperBounds;

  MOType motype;

  //---optimization-helpers---------------------------------------------------//

  nlopt::opt * optimizer;

  SCVT nloptFRel;

  SCVT cost;

  SCVT costOpt;

  SCVT costFirst;

  SCVT gradNormFirst;

  bool normalizeCost;

  bool normalizeGrad;

  bool normalizeOnEachLevel;

  std::vector< double > xOpt;

  std::vector< double > alphaOpt;

  int iter;

  int iterLevel;

  int jump;

  int nloptJumpTol;


  // IPOPT
#ifdef IPOPT

public:
  //! Method to return some info about the nlp
  virtual bool get_nlp_info(
      Ipopt::Index & n,
      Ipopt::Index & m,
      Ipopt::Index & nnz_jac_g,
      Ipopt::Index & nnz_h_lag,
      IndexStyleEnum & index_style
      );

  //! Method to return the bounds for my problem
  virtual bool get_bounds_info(
      Ipopt::Index n,
      Ipopt::Number * x_l,
      Ipopt::Number * x_u,
      Ipopt::Index m,
      Ipopt::Number * g_l,
      Ipopt::Number * g_u
      );

  //! Method to return the starting point for the algorithm
  virtual bool get_starting_point(
      Ipopt::Index n,
      bool init_x,
      Ipopt::Number * x,
      bool init_z,
      Ipopt::Number * z_L,
      Ipopt::Number * z_U,
      Ipopt::Index m,
      bool init_lambda,
      Ipopt::Number * lambda
      );

  //! Method to return the objective value
  virtual bool eval_f(
      Ipopt::Index n,
      const Ipopt::Number * x,
      bool new_x,
      Ipopt::Number & obj_value
      );

  //! Method to return the gradient of the objective
  virtual bool eval_grad_f(
      Ipopt::Index n,
      const Ipopt::Number * x,
      bool new_x,
      Ipopt::Number * grad_f
      );

  //! Method to return the constraint residuals
  virtual bool eval_g(
      Ipopt::Index n,
      const Ipopt::Number * x,
      bool new_x,
      Ipopt::Index m,
      Ipopt::Number * g
      );

  /*! 
   * Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(
      Ipopt::Index n,
      const Ipopt::Number * x,
      bool new_x,
      Ipopt::Index m,
      Ipopt::Index nele_jac,
      Ipopt::Index * iRow,
      Ipopt::Index * jCol,
      Ipopt::Number * values
      );

  /*!
   *  This method is called when the algorithm is complete so the TNLP can 
   * store/write the solution
   */
  virtual void finalize_solution(
      Ipopt::SolverReturn status,
      Ipopt::Index n,
      const Ipopt::Number * x,
      const Ipopt::Number * z_L,
      const Ipopt::Number * z_U,
      Ipopt::Index m,
      const Ipopt::Number * g,
      const Ipopt::Number * lambda,
      Ipopt::Number obj_value,
      const Ipopt::IpoptData * ip_data,
      Ipopt::IpoptCalculatedQuantities * ip_cq
      );

  virtual bool intermediate_callback(
      Ipopt::AlgorithmMode mode,
      Ipopt::Index iter,
      Ipopt::Number obj_value,
      Ipopt::Number inf_pr,
      Ipopt::Number inf_du,
      Ipopt::Number mu,
      Ipopt::Number d_norm,
      Ipopt::Number regularization_size,
      Ipopt::Number alpha_du,
      Ipopt::Number alpha_pr,
      Ipopt::Index ls_trials,
      const Ipopt::IpoptData * ip_data,
      Ipopt::IpoptCalculatedQuantities * ip_cq
      );

  void setIpoptTol(
      double ipopt_tol
      ) {
    this->ipopt_tol = ipopt_tol;
  }

  void setIpoptAcceptableIter(
      int ipopt_acceptable_iter
      ) {
    this->ipopt_acceptable_iter = ipopt_acceptable_iter;
  }

  void setIpoptAcceptableObjChangeTol(
      double ipopt_acceptable_obj_change_tol
      ) {
    this->ipopt_acceptable_obj_change_tol = ipopt_acceptable_obj_change_tol;
  }

  void setIpoptAcceptableTol(
      double ipopt_acceptable_tol
      ) {
    this->ipopt_acceptable_tol = ipopt_acceptable_tol;
  }

  void setIpoptMaxIter(
      int ipopt_max_iter
      ) {
    this->ipopt_max_iter = ipopt_max_iter;
  }

private:

  // fixme: dirty, dirty hack
  MultiresolutionOptimizer< LO, SC > * getShallowCopy( );

  void shallowCopy(
      MultiresolutionOptimizer< LO, SC > & copy
      ) const;

  /*!
   * Performs free form optimization on current level
   * every iteration is obtained by a perturbation of original domain
   */
  void optimizeFixedReferenceIpopt( );

  bool iAmShallowCopy;

  double ipopt_tol;

  double ipopt_acceptable_tol;

  int ipopt_max_iter;

  int ipopt_acceptable_iter;

  double ipopt_acceptable_obj_change_tol;

#endif

};

}

// include .cpp file to overcome linking problems due to templates
#include "MultiresolutionOptimizer.cpp"

#endif	/* MULTIRESOLUTIONOPTIMIZER_H */
