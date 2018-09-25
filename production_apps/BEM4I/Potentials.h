/*!
 * @file    Potentials.h
 * @author  Michal Merta 
 * @date    May 20, 2014
 * @brief   Header file for pure virtual class Potentials
 * 
 */

#ifndef POTENTIALS_H
#define	POTENTIALS_H

#include "Vector.h"
//#include "BESpace.h"
#include "Macros.h"
#include "Mesh.h"

namespace bem4i {
  
  // forward declaration
  template<class LO3, class SC3>
  class BESpace;

/*! 
 * Abstract class for potential operators
 * 
 */
template<class LO, class SC>
class Potentials {
  


  typedef typename GetType<LO, SC>::SCVT SCVT;

public:


  /*!
   * function evaluates single layer potential in given array of points
   *  
   * @param[in]       x
   * @param[in,out]   result
   */
  virtual void singleLayerPotential(
      const SCVT * x, 
      LO n,
      Vector<LO, SC> & values,
      SCVT t = 0.0
      ) const = 0;

  /*!
   * function evaluates double layer potential in given array of points
   *  
   * @param[in]       x
   * @param[in,out]   result
   */
  virtual void doubleLayerPotential(
      const SCVT * x,
      LO n,
      Vector<LO, SC> & values,
      SCVT t = 0.0
      ) const = 0;

  /*!
   * function evaluates single layer potential in nodes of given mesh
   *  
   * @param[in]       mesh input mesh
   * @param[in]       interior interior/exterior flag
   * @param[in,out]   values
   */
  virtual void singleLayerPotential(
      Mesh<LO, SC> & mesh,
      Vector<LO, SC> & values,
      SCVT t = 0.0
      ) const = 0;

  /*!
   * function evaluates double layer potential in nodes of given mesh
   *  
   * @param[in]       mesh input mesh
   * @param[in]       interior interior/exterior flag
   * @param[in,out]   values
   */
  virtual void doubleLayerPotential(
      Mesh<LO, SC> & mesh,
      Vector<LO, SC> & values,
      SCVT t = 0.0
      ) const = 0;


protected:

  //! boundary element space
  BESpace<LO, SC> * space;

  //! density function
  Vector<LO, SC> * density;


};

} // end of namespace bem4i

#endif	/* POTENTIALS_H */

