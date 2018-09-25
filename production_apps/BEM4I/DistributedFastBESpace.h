/*!
 * @file    DistributedFastBESpace.h
 * @author  Michal Merta 
 * @date    February 17, 2016
 * @brief   Header file for class DistributedFastBESpace
 * 
 */

#ifndef DISTRIBUTEDFASTBESPACE_H
#define	DISTRIBUTEDFASTBESPACE_H

#include <algorithm>

#include "BESpace.h"


namespace bem4i {


/*!
 * BE space for fast boundary element methods with distributed matrices
 * 
 * the class is responsible mainly for decomposition of input mesh by METIS
 */
template<class LO, class SC>
class DistributedFastBESpace : public FastBESpace<LO, SC> {

typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  /*!
   * constructor taking a mesh and ansatzs types
   * constructor then creates a cluster tree based on given mesh, 
   * sets (possibly default) fast BEM parameters, creates a block cluster tree, and
   * creates (non)admissible pairs
   * 
   * @param[in]   mesh
   * @param[in]   innerAnsatzType   ansatz type of the inner integral (constant/linear)
   * @param[in]   outerAnsatzType   ansatz type of the outer integral (constant/linear)
   * @param[in]   eta               admissibility constant
   */
  DistributedFastBESpace(
      SurfaceMesh3D<LO, SC>* mesh,
      LO nSubmeshes,
      basisType ansatzFunctionType,
      basisType testFunctionType,
      SC eta = 1.0,
      bool groupAdmClusters = true,
      LO maxElems = 0
      );


  //! destructor
  virtual ~DistributedFastBESpace( );

  //! setter for eta

  inline void setEta( SC eta ) {
    this->eta = eta;
  }

  //! getter for eta
  inline SC getEta( ) const {
    return eta;
  }

  //! setter for nMax
  inline void setNMax( int nMax ) {
    this->nMax = nMax;
  }

  //! getter for nMax
  inline int getNMax( ) const {
    return nMax;
  }

  //! setter for quadrature orders of fmm integrals
  inline void setQuadOrder( int quadOrder ) {
    this->quadOrder = quadOrder;
  }

  //! getter for quadrature order
  inline int getQuadOrder( ) const {
    return quadOrder;
  }

  //! setter for starting level for passing tree
  inline void setStartLevel( int startLevel ) {
    this->startLevel = startLevel;
  }

  //! getter for starting level for passing tree
  inline int getStartLevel( ) const {
    return startLevel;
  }

  //! setter for fast bem parameters
  inline void setFastBEMParams(
      SCVT eta = 1.0,
      int nMax = 6,
      int quadOrder = 7,
      int startLevel = 0
      ) {
    this->eta = eta;
    this->nMax = nMax;
    this->quadOrder = quadOrder;
    this->startLevel = startLevel;
  }

  //! setter for aca epsilon
  inline void setEpsilonACA(
      SCVT epsilon ) {
    this->epsilonACA = epsilon;
  }

  //! setter for aca epsilon
  inline SCVT getEpsilonACA( ) {
    return this->epsilonACA;
  }

  //! setter for ACA zero scale factor
  inline void setScaleACA(
      SCVT scale ) {
    this->scaleACA = scale;
  }

  //! setter for aca epsilon
  inline SCVT getScaleACA( ) {
    return this->scaleACA;
  }
  
  //! setter for max elements per cluster
  inline void setMaxElemsPerCluster( LO maxElemsPerCluster ) {
    this->maxElemsPerCluster = maxElemsPerCluster;
  }

  
protected:

  //! addmissibility constant
  SCVT eta;

  //! default maximum number of elements per cluster
  LO defaultMaxElemsPerCluster = 15;

  //! max index of fmm expansion
  int nMax;

  //! starting level for passing tree
  int startLevel;

  //! quadrature orders of integrals in fmm
  int quadOrder;

  //! epsilon for ACA
  SCVT epsilonACA;
  
  //! scaling factor for approximation of zero in ACA
  SCVT scaleACA;

  //! group admissible clusters together
  bool groupAdmClusters;
  
  //! maximum number of elements per cluster
  LO maxElemsPerCluster;

  //! number of submeshes to split an input mesh into
  LO nSubmeshes;
  
  //! array of pointers to submeshes
  SurfaceMesh3D<LO, SC> ** submeshes;
  
  //! array of FastBESpaces associated with individual submeshes
  FastBESpace<LO, SC> ** bespaces;
  
  //! global indices corresponding to local indic
  std::vector<LO> ** outerDOFs;

  std::vector<LO> ** innerDOFs;
  
  //! array of trees (one per submesh)
  Tree<BECluster<LO, SC>*, LO> ** trees; 

  //! calls metis to decompose an input mesh into submeshes
  void decomposeMesh();
  

private:
    //! default constructor
  DistributedFastBESpace( );

  //! copy constructor
  DistributedFastBESpace( const DistributedFastBESpace& orig );

};
}

// include .cpp file to overcome linking problems due to templates
#include "DistributedFastBESpace.cpp"

#endif	/* DISTRIBUTEDFASTBESPACE_H */

