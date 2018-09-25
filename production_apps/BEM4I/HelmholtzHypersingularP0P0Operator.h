/*!
 * @file    HelmholtzHypersingularP0P0Operator.h
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    July 23, 2015
 * @brief   Header file for class HelmholtzHypersingularP0P0Operator
 */

#ifndef HELMHOLTZHYPERSINGULARP0P0OPERATOR_H
#define	HELMHOLTZHYPERSINGULARP0P0OPERATOR_H

#include <vector>

#include "LinearOperator.h"
#include "Vector.h"
#include "BESpace.h"
#include "SparseMatrix.h"
#include "IterativeSolver.h"

namespace bem4i {

template<class LO, class SC>
class HelmholtzHypersingularP0P0Operator :
public LinearOperator<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  /*!
   * @param H1     
   * @param H2     
   */
  HelmholtzHypersingularP0P0Operator(
      LinearOperator< LO, SC > * H1,
      LinearOperator< LO, SC > * H2,
      SurfaceMesh3D<LO, SC> & mesh
      );

  /*!
   * @brief Applies operator on a vector ( y = beta * y + alpha * this * x )
   *
   * @param x      operator argument
   * @param y      result (user pre-allocated)
   * @param alpha  y = beta * y + alpha * this * x
   * @param beta   y = beta * y + alpha * this * x
   */
  void apply(
      Vector<LO, SC> const & x,
      Vector<LO, SC> & y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

 

  //! Destructor
  virtual ~HelmholtzHypersingularP0P0Operator( );
  
private:

  
  LinearOperator< LO, SC > * H1;

  LinearOperator< LO, SC > * H2;
  
  Vector<LO, SC> *diagElems1;
  
  Vector<LO, SC> *diagElems2;

 
  //! Default constructor (empty)
  HelmholtzHypersingularP0P0Operator( );

  //! Copy constructor (empty)
  HelmholtzHypersingularP0P0Operator(
      const HelmholtzHypersingularP0P0Operator & orig
      );

};

}
// include .cpp file to overcome linking problems due to templates
#include "HelmholtzHypersingularP0P0Operator.cpp"

#endif	/* HELMHOLTZHYPERSINGULARP0P0OPERATOR_H */
