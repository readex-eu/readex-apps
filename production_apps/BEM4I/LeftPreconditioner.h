/*!
 * @file    LeftPreconditioner.h
 * @author  Michal Merta 
 * @date    March 12, 2014
 * @brief   Header file for abstract class LeftPreconditioner
 * 
 */

#ifndef LEFTPRECONDITIONER_H
#define	LEFTPRECONDITIONER_H

#include "Preconditioner.h"


namespace bem4i {

/*! 
 * Abstract class representing a left preconditioner of a system of lin. eq.
 * 
 * the class applies the inverse of a matrix spectrally equivallent to a given system matrix
 */
template<class LO, class SC>
class LeftPreconditioner : public Preconditioner<LO, SC> {

public:

  typedef typename GetType<LO, SC>::SCVT SCVT;

private:

};

/*! 
 * Class applies identity operator as a "preconditioner"
 * 
 */
template<class LO, class SC>
class LeftIdentityPreconditioner : public LeftPreconditioner<LO, SC> {

public:

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
    y.add( x, alpha );
  }
};

}

#endif	/* LEFTPRECONDITIONER_H */
