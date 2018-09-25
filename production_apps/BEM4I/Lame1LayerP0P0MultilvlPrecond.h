/*!
 * @file    Lame1LayerP0P0MultilvlPrecond.h
 * @author  Lukas Maly 
 * @date    July 20, 2016
 * @brief   Header file for class Lame1LayerP0P0MultilvlPrecond
 * 
 */


#ifndef LAME1LAYERP0P0MULTILVLPRECOND_H
#define LAME1LAYERP0P0MULTILVLPRECOND_H

#include <vector>

#include "Laplace1LayerP0P0MultilvlPrecond.h"
#include "LinearOperator.h"
#include "BlockLinearOperator.h"
#include "Vector.h"
#include "Tree.h"
#include "BESpace.h"

namespace bem4i {

/*!
 * Class representing an Artificial Multilevel Preconditioner
 */
template<class LO, class SC>
class Lame1LayerP0P0MultilvlPrecond : public LinearOperator<LO, SC> {
  
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! constructor taking a mesh as argument
 Lame1LayerP0P0MultilvlPrecond(
      SurfaceMesh3D<LO, SC> * mesh,
      LO maxElems = 10
      );

  /*
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
 
  virtual  ~Lame1LayerP0P0MultilvlPrecond( );

private:
  
 Lame1LayerP0P0MultilvlPrecond( ) {
    this->dimRange = 0;
    this->dimDomain = 0;
  }

  void init( );

  // preconditioner for V
  Laplace1LayerP0P0MultilvlPrecond<LO, SC> * Vprecond;

  // block linerar operator 
  BlockLinearOperator<LO, SC> * LameBlockPrecond;

  // zero operator 
  Zero<LO, SC> * O;

};

} //end of namespace
#include "Lame1LayerP0P0MultilvlPrecond.cpp"

#endif /* Lame1LAYERP0P0MULTILVLPRECOND_H */
