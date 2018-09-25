/*!
 * @file    BESpaceTime.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    December 12, 2013
 */

#ifdef BESPACETIME_H

namespace bem4i {

template<class LO, class SC>
BESpaceTime<LO, SC>::BESpaceTime( ) {
}

template<class LO, class SC>
BESpaceTime<LO, SC>::BESpaceTime( const BESpaceTime& orig ) {
}

template<class LO, class SC>
BESpaceTime<LO, SC>::~BESpaceTime( ) {
}

template<class LO, class SC>
BESpaceTime<LO, SC>::BESpaceTime(
    SurfaceMesh3D<LO, SC>* mesh,
    basisType ansatzFunctionType,
    basisType testFunctionType,
    int legendreOrder,
    SCVT finalTime,
    LO nTimeSteps
    ) : BESpace<LO, SC>(mesh, ansatzFunctionType, testFunctionType) {
  this->legendreOrder = legendreOrder;
  this->finalTime = finalTime;
  this->nTimeSteps = nTimeSteps;
  this->dt = finalTime / (nTimeSteps - 1);
}

}

#endif
