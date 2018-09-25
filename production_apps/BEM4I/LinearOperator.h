/*!
 * @file    LinearOperator.h
 * @author  Michal Merta 
 * @date    July 8, 2013
 * @brief   Header file for abstract class LinearOperator
 * 
 */

#ifndef LINEAROPERATOR_H
#define	LINEAROPERATOR_H

#include "Vector.h"
#include "FullMatrix.h"
#include "LeftPreconditioner.h"
#include "IterativeSolver.h"

namespace bem4i {

template <class LO, class SC>
class FullMatrix;

template <class LO, class SC>
class Vector;

/*! 
 * Abstract class representing arbitrary linear operator
 * 
 */
template<class LO, class SC>
class LinearOperator {
public:

  typedef typename GetType<LO, SC>::SCVT SCVT;

  //! default constructor
  // todo: implement this constructor for every subclass, 
  // make the default private

  LinearOperator( ) {
    this->dimRange = 0;
    this->dimDomain = 0;
  }

  LinearOperator(
      LO dimRange,
      LO dimDomain
      ) {
  };

  //! destructor

  virtual ~LinearOperator( ) {
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
      ) = 0;

  inline LO getDimRange( ) {
    return this->dimRange;
  }

  inline LO getDimDomain( ) {
    return this->dimDomain;
  }

  /*!
   * @brief CG solver for symmetric, positive definite matrices
   * 
   * @param rhs   right-hand side vector
   * @param sol   solution vector
   * @param prec  solver relative precision
   */
  bool CGSolve(
      const Vector<LO, SC> & rhs,
      Vector<LO, SC> & sol,
      SCVT prec = 1e-6,
      LO maxIt = 1000,
      LinearOperator<LO, SC> * M = nullptr,
      const std::string & msg = ""
      ) {
    return IterativeSolver<LO, SC>::CGSolve( *this, rhs, sol, prec, maxIt, M, msg );
  }

  /*!
   * @brief GMRES solver for general nonsymmetric systems
   * 
   * @param rhs[in]       right-hand side vector
   * @param sol[in,out]   solution vector
   * @param prec[in]      solver relative precision
   * @param restarts[in]  number of iterations per before restart
   * @param maxIt[in]     maximum number of iterations
   */
  bool GMRESSolve(
      Vector<LO, SC> const &rhs,
      Vector<LO, SC> &sol,
      SCVT precision,
      LO maxIt,
      LO restarts = 1000,
      LinearOperator<LO, SC> *M = nullptr,
      const std::string & msg = ""
      ) {
    return IterativeSolver<LO, SC>::GMRESSolve( *this, rhs, sol, precision,
        maxIt, restarts, M, msg );
  }

  /*!
   * @brief FGMRES solver for general nonsymmetric systems
   * 
   * @param rhs[in]       right-hand side vector
   * @param sol[in,out]   solution vector
   * @param prec[in]      solver relative precision
   * @param restarts[in]  number of iterations per before restart
   * @param maxIt[in]     maximum number of iterations
   */
  bool FGMRESSolve(
      Vector<LO, SC> const &rhs,
      Vector<LO, SC> &sol,
      SCVT precision,
      LO maxIt,
      LO restarts = 1000,
      LinearOperator<LO, SC> *M = nullptr,
      const std::string & msg = ""
      ) {
    return IterativeSolver<LO, SC>::FGMRESSolve( *this, rhs, sol, precision,
        maxIt, restarts, M, msg );
  }

  /*!
   * @brief DGMRES with defaltions solver for general nonsymmetric systems
   * 
   * @param rhs[in]       right-hand side vector
   * @param sol[in,out]   solution vector
   * @param prec[in]      solver relative precision
   * @param restarts[in]  number of iterations per before restart
   * @param maxInvSpaceDim[in]     maximum number of vector in the basis of the 
   *                               invariant subspace
   * @param maxEigsPerRestart[in]  maximum number of vector to add in one 
   *                               restart           
   * @param maxIt[in]     maximum number of iterations
   */
  bool DGMRESSolve(
      Vector<LO, SC> const &rhs,
      Vector<LO, SC> &sol,
      SCVT precision,
      LO maxIt,
      LO restarts = 10,
      LO maxInvSpaceDim = 1,
      LO maxEigsPerRestart = 1,
      LinearOperator<LO, SC> *M = nullptr,
      const std::string & msg = ""
      ) {

    return IterativeSolver<LO, SC>::DGMRESSolve( *this, rhs, sol, precision,
        maxIt, restarts, maxInvSpaceDim, maxEigsPerRestart, M, msg );
  }


protected:

  LO dimRange;
  LO dimDomain;

private:

};

template < class LO, class SC >
class Zero : public LinearOperator< LO, SC > {
public:

  Zero(
      LO dimRange,
      LO dimDomain
      ) {

    this->dimRange = dimRange;
    this->dimDomain = dimDomain;
  }

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
      Vector<LO, SC> const & x,
      Vector<LO, SC> & y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      ) {

    y.scale( beta );
  }

};

}

#endif	/* LINEAROPERATOR_H */
