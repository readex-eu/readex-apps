/*!
 * @file    PotentialsLaplace.h
 * @author  Jan Zapletal
 * @date    February 17, 2015
 * @brief   Header file for class PotentialsLaplace
 * 
 */

#ifndef POTENTIALSLAPLACE_H
#define	POTENTIALSLAPLACE_H

#include "Vector.h"
//#include "BESpace.h"
#include "Macros.h"
#include "Mesh.h"
#include "Potentials.h"
#include "BEIntegratorLaplace.h"

namespace bem4i {

  // forward declaration
  template<class LO3, class SC3>
  class BESpace;

/*! 
 * Class representing single and double layer potentials for Laplace equation
 * 
 */
template<class LO, class SC>
class PotentialsLaplace : public Potentials<LO, SC> {

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  /*!
   * constructor taking Dirichlet and Neumann data as arguments
   * 
   * @param[in]     space
   * @param[in]     dirichlet - vector of dirichlet data
   * @param[in]     neumann - vector of neumann data
   * 
   */
  PotentialsLaplace(
      BESpace<LO, SC> * space,
      Vector<LO, SCVT> * density
      );

  //! destructor
  virtual ~PotentialsLaplace( );

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
      ) const {
  };

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
      ) const;

  /*!
   * function evaluates single layer potential in nodes of given mesh
   *  
   * @param[in]       mesh input mesh
   * @param[in,out]   values
   */
  virtual void singleLayerPotential(
      Mesh<LO, SC> & mesh,
      Vector<LO, SC> & values,
      SCVT t = 0.0
      ) const {
  };

  /*!
   * function evaluates double layer potential in nodes of given mesh
   *  
   * @param[in]       mesh input mesh
   * @param[in,out]   values
   */
  virtual void doubleLayerPotential(
      Mesh<LO, SC> & mesh,
      Vector<LO, SC> & values,
      SCVT t = 0.0
      ) const;

private:

  PotentialsLaplace( ) {
  };

  //! copy constructor
  PotentialsLaplace(
      const PotentialsLaplace & orig
      );

};

} // end of namespace bem4i

// include .cpp file to overcome linking problems due to templates
#include "PotentialsLaplace.cpp"

#endif	/* POTENTIALSLAPLACE_H */

