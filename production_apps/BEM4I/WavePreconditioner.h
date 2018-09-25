/*!
 * @file    WavePreconditioner.h
 * @author  Michal Merta 
 * @date    March 12, 2014
 * @brief   Header file for abstract class for preconditioner for wave equation
 * 
 */

#ifndef WAVEPRECONDITIONER_H
#define	WAVEPRECONDITIONER_H

#include "LeftPreconditioner.h"
#include "WavePreconditionerTriangular.h"
#include "BlockMatrix.h"
#include "MPIBlockMatrix.h"

namespace bem4i {

template<class LO, class SC>
class WavePreconditioner : public LeftPreconditioner<LO, SC> {

public:

  //! default constructor
  WavePreconditioner( );

  /*!
   * @brief Constructor taking system matrix for wave equation as argument
   * 
   * Constructor taking hypersingular/single layer matrix for wave equation as argument
   * @param A         system matrix
   * @param maxLevel  maximum level of recursion of preconditioner
   * @param maxIters  maximum number of iterations on each level (from top)
   */
  WavePreconditioner(
      BlockMatrix<LO, SC> *A,
      int maxLevel = -1, 
      int *maxIters = nullptr
      );

  /*!
   * @brief Constructor taking system matrix for wave equation as argument
   * 
   * Constructor taking hypersingular/single layer matrix for wave equation as argument
   * @param A     system matrix
   * @param maxLevel  maximum level of recursion of preconditioner
   * @param maxIters  maximum number of iterations on each level (from top)
   */
  WavePreconditioner(
      MPIBlockMatrix<LO, SC> *A,
      int maxLevel = -1,
      int *maxIters = nullptr
      );

  //! copy constructor
  WavePreconditioner(
      const WavePreconditioner& orig
      );

  //! destructor
  virtual ~WavePreconditioner( );

  /*!
   * @brief Performs a matrix-vector multiplication
   * 
   * Computes y = beta*y + alpha*this*x
   * @param A 
   * @param x
   * @param y
   * @param alpha
   * @param beta
   */
  virtual void apply(
      Vector<LO, SC> const &x,
      Vector<LO, SC> &y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

  //! sets current level of preconditioner application recursion

  inline void setCurrentLevel(
      int level
      ) {
    this->level = level;
  }

  //! gets current level of preconditioner application recursion

  inline int getCurrentLevel( ) {
    return this->level;
  }

private:

  //! system matrix for wave equation
  BlockMatrix<LO, SC> * sysMatrix;

  //! lower block triangular part of a system matrix
  BlockMatrix<LO, SC> * preconditioner;

  //! system matrix for wave equation
  MPIBlockMatrix<LO, SC> * sysMatrixMPI;

  //! lower block triangular part of a system matrix
  MPIBlockMatrix<LO, SC> * preconditionerMPI;

  //! applies preconditioner based on block matrix
  void applyPreconditioner(
      BlockMatrix<LO, SC> &matrix,
      Vector<LO, SC> const &x,
      Vector<LO, SC> &y );

  //! applies preconditioner based on block matrix
  void applyPreconditioner(
      MPIBlockMatrix<LO, SC> &matrix,
      Vector<LO, SC> const &x,
      Vector<LO, SC> &y );

//  void applyGSPreconditioner(
//      BlockMatrix<LO, SC> &matrix,
//      Vector<LO, SC> const &x,
//      Vector<LO, SC> &y );
//
//  void applyLInverse(
//      BlockMatrix<LO, SC> &matrix,
//      Vector<LO, SC> const &x,
//      Vector<LO, SC> &y );

  //! current level of recursion
  int level;

  //! maximum level of recursion
  int maxLevel;
  
  //! maximum number of iterations on each level of recursion (start from top)
  int *maxIters;

  //! default max iterations on each level
  const int defaultMaxIters = 20;
  
  //! whether deallocate maxIters
  bool deleteMaxIters;


};

}

// include .cpp file to overcome linking problems due to templates
#include "WavePreconditioner.cpp"

#endif	/* WAVEPRECONDITIONER_H */

