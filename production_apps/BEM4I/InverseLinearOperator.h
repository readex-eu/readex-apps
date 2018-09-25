/*!
 * @file    InverseLinearOperator.h
 * @author  Jan Zapletal
 * @date    April 26, 2016
 * @brief   Header file for the class InverseLinearOperator
 * 
 */

#ifndef INVERSELINEAROPERATOR_H
#define	INVERSELINEAROPERATOR_H

#include "LinearOperator.h"
#include "Vector.h"
#include <vector>

namespace bem4i {

template< class LO, class SC >
class InverseLinearOperator;

template< class LO, class SC >
struct SolverParameters {
  typedef typename GetType<LO, SC>::SCVT SCVT;

  enum IterativeSolver {
    CG,
    GMRES,
    FGMRES,
    DGMRES
  };

  IterativeSolver solver;
  SCVT precision;
  LO maxIter;
  std::string msg;

  SolverParameters( ) {
    solver = CG;
    precision = 1e-8;
    maxIter = 1000;
    msg = "";
  }
  
};

template<class LO, class SC>
class InverseLinearOperator : public LinearOperator<LO, SC> {
public:

  typedef typename GetType<LO, SC>::SCVT SCVT;

  InverseLinearOperator(
      LinearOperator< LO, SC > * op,
      SolverParameters< LO, SC > * params,
      LinearOperator< LO, SC > * precond = nullptr
      ) {

    this->op = op;
    this->params = params;
    this->precond = precond;
    
    this->dimDomain = op->getDimRange( );
    this->dimRange = op->getDimDomain( );
  };

  //! destructor

  virtual ~InverseLinearOperator( ) {
  };

  /*!
   * @brief Applies operator on a vector
   * 
   * Computes y = beta*y + alpha*this*x
   * @param A 
   * @param x
   * @param y
   * @param alpha
   * @param beta
   */
  virtual void apply(
      Vector< LO, SC > const & x,
      Vector< LO, SC > & y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

protected:

  LinearOperator< LO, SC > * op;
  LinearOperator< LO, SC > * precond;
  SolverParameters< LO, SC > * params;

private:

};

}

#include "InverseLinearOperator.cpp"

#endif	/* INVERSELINEAROPERATOR_H */
