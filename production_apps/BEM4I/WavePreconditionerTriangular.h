/*!
 * @file    WavePreconditionerTriangular.h
 * @author  Michal Merta 
 * @date    March 12, 2014
 * @brief   Header file for abstract class for preconditioner for wave equation
 * 
 */

#ifndef WAVEPRECONDITIONERTRIANGULAR_H
#define	WAVEPRECONDITIONERTRIANGULAR_H

#include "LeftPreconditioner.h"
#include "BlockMatrix.h"

namespace bem4i {

template<class LO, class SC>
class WavePreconditionerTriangular : public LeftPreconditioner<LO, SC> {

public:

  //! default constructor
  WavePreconditionerTriangular( );

  /*!
   * @brief Constructor taking system matrix for wave equation as argument
   * 
   * Constructor taking hypersingular/single layer matrix for wave equation as argument
   * @param A     system matrix
   */
  WavePreconditionerTriangular( BlockMatrix<LO, SC> *A, int maxLevel = -1 );

  //! copy constructor
  WavePreconditionerTriangular( const WavePreconditionerTriangular& orig );

  //! destructor
  virtual ~WavePreconditionerTriangular( );

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
  virtual void apply( Vector<LO, SC> const &x, Vector<LO, SC> &y, bool transA = false,
      SC alpha = 1.0, SC beta = 0.0 );
  
  //! sets current level of preconditioner application recursion
  inline void setCurrentLevel(int level) {
    this->level = level;
  }
  
  //! gets current level of preconditioner application recursion
  inline int getCurrentLevel() {
    return this->level;
  }

private:

  //! system matrix for wave equation
  BlockMatrix<LO, SC> *sysMatrix;

  //! lower block triangular part of a system matrix
  BlockMatrix<LO, SC> *preconditioner;
  

};

}

// include .cpp file to overcome linking problems due to templates
#include "WavePreconditionerTriangular.cpp"

#endif	/* WAVEPRECONDITIONERTRIANGULAR_H */

