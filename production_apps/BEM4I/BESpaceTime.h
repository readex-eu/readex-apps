/*!
 * @file    BESpaceTime.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    December 12, 2013
 * @brief   Header file for class BESpaceTime
 * 
 */

#ifndef BESPACETIME_H
#define	BESPACETIME_H

#include "BESpace.h"

namespace bem4i {

template<class LO, class SC>
class BESpaceTime : public BESpace<LO, SC> {

  // to get inner type of complex numbers (for Helmholtz)
  typedef typename GetType<LO, SC>::SCVT SCVT;


public:
  //! default constructor
  BESpaceTime( );

  //! copy constructor
  BESpaceTime( const BESpaceTime& orig );

  /*!
   * constructor taking the final time and the number of timesteps
   * 
   * @param[in]   mesh
   * @param[in]   innerAnsatzType   ansatz type of the inner integral (constant/linear)
   * @param[in]   outerAnsatzType   ansatz type of the outer integral (constant/linear)
   */
  BESpaceTime(
      SurfaceMesh3D<LO, SC>* mesh,
      basisType ansatzFunctionType,
      basisType testFunctionType,
      int legendreOrder,
      SCVT finalTime,
      LO nTimeSteps
      ) ;

  //! destructor
  virtual ~BESpaceTime( );
  
  //! returns number of time-steps
  inline LO getNTimeSteps() const {
    return nTimeSteps;
  }
  
  //! returns final time
  inline SC getFinalTime() const {
    return finalTime;
  }
  
  inline SCVT getDt() const {
    return dt;
  }
  
  inline int getLegendreOrder() const {
    return legendreOrder;
  }

private:
  
  protected:
    
    //! T = <0, finalTime>
    SCVT finalTime;
    
    //! number of time-steps
    LO nTimeSteps;
    
    SCVT dt;
    
    //! order of Legendre polynomials used for the construction of temporal basis functions
    int legendreOrder;

};

}

// include .cpp file to overcome linking problems due to templates
#include "BESpaceTime.cpp"

#endif	/* BESPACETIME_H */

