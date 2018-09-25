/*!
 * @file    PotentialsWave.h
 * @author  Michal Merta 
 * @date    May 20, 2014
 * @brief   Header file for class PotentialsWave
 * 
 */

#ifndef POTENTIALSWAVE_H
#define	POTENTIALSWAVE_H

#include "Vector.h"
#include "BESpace.h"
#include "Macros.h"
#include "Mesh.h"
#include "Potentials.h"

namespace bem4i {

/*! 
 * Class representing single and double layer potentials for wave equation
 * 
 */
template<class LO, class SC>
class PotentialsWave : public Potentials<LO, SC> {

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! copy constructor

  PotentialsWave( ) {

  };

  //! copy constructor
  PotentialsWave(
      const PotentialsWave& orig
      );

  /*!
   * constructor taking Dirichlet and Neumann data as arguments
   * 
   * @param[in]     space
   * @param[in]     dirichlet - vector of dirichlet data
   * @param[in]     neumann - vector of neumann data
   * 
   */
  PotentialsWave(
      BESpace<LO, SC> * space,
      Vector<LO, SCVT> * density
      );

  //! destructor
  virtual ~PotentialsWave( );

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
   * @param[in]       interior interior/exterior flag
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
   * @param[in]       interior interior/exterior flag
   * @param[in,out]   values
   */
  virtual void doubleLayerPotential(
      Mesh<LO, SC> & mesh,
      Vector<LO, SC> & values,
      SCVT t = 0.0
      ) const;


};

} // end of namespace bem4i

// include .cpp file to overcome linking problems due to templates
#include "PotentialsWave.cpp"

#endif	/* POTENTIALSWAVE_H */

