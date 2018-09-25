/*!
 * @file    FastBESpace.h
 * @author  Michal Merta 
 * @date    August 14, 2013
 * @brief   Header file for class FastBESpace
 * 
 */

#ifndef FASTBESPACE_H
#define	FASTBESPACE_H

#include <algorithm>

#include "BESpace.h"
#include "Tree.h"


namespace bem4i {

//! struct representing pairs of boundary element clusters

template<class LO, class SC>
struct BEBlockCluster {

  BEBlockCluster( ) {
    leftCluster = nullptr;
    rightCluster = nullptr;
    outerDOFs = nullptr;
    innerDOFs = nullptr;
  }

  ~BEBlockCluster( ) {
    if ( outerDOFs ) {
      delete outerDOFs;
    }
    if ( innerDOFs ) {
      delete innerDOFs;
    }
  }

  //! left (row) cluster
  BECluster<LO, SC>* leftCluster;

  //! right (column) cluster
  BECluster<LO, SC>* rightCluster;

  //! degrees of freedom associated with left cluster
  std::vector<LO>* outerDOFs;

  //! degrees of freedom associated with right cluster
  std::vector<LO>* innerDOFs;

  //! whether the clusters in pair are mutually admissible
  bool admissible;

  //! whether the cluster is a leaf
  bool leaf;
};

/*!
 * BE space for fast boundary element methods
 * 
 * the class is responsible for finding (non)admissible pair of BE clusters
 */
template<class LO, class SC>
class FastBESpace : public BESpace<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! default constructor
  FastBESpace( );

  //! copy constructor
  FastBESpace(
      const FastBESpace & orig
      );

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
  FastBESpace(
      SurfaceMesh3D<LO, SC>* mesh,
      basisType ansatzFunctionType,
      basisType testFunctionType,
      SC eta = 1.0,
      int nMax = 6,
      int quadOrder = 7,
      int startLevel = 0,
  //    bool groupAdmClusters = true,
      LO maxElems = 0
      );

  /*!
   * constructor taking a mesh and ansatzs types and left and right cluster tree
   * sets (possibly default) fast BEM parameters, creates a block cluster tree, and
   * creates (non)admissible pairs
   * @param[in]   mesh
   * @param[in]   innerAnsatzType   ansatz type of the inner integral (constant/linear)
   * @param[in]   outerAnsatzType   ansatz type of the outer integral (constant/linear)
   * @param[in]   leftClusterTree   cluster tree associated with row part of the mesh
   * @param[in]   rightClusterTree  cluster tree associated with column spart of the mesh
   * @param[in]   eta               admissibility constant
   */
  FastBESpace(
      SurfaceMesh3D<LO, SC>* mesh,
      basisType ansatzFunctionType,
      basisType testFunctionType,
      Tree<BECluster<LO, SC>*, LO >* leftClusterTree,
      Tree<BECluster<LO, SC>*, LO >* rightClusterTree,
      SC eta = 1.0,
      int nMax = 6,
      int quadOrder = 7,
      int startLevel = 0,
//      bool groupAdmClusters = true,
      LO maxElems = 0
      );

  /*!
   * constructor taking a mesh and ansatzs types and left and right cluster tree
   * sets (possibly default) fast BEM parameters, creates a block cluster tree, and
   * creates (non)admissible pairs
   * @param[in]   leftMesh
   * @param[in]   rightMesh
   * @param[in]   innerAnsatzType   ansatz type of the inner integral (constant/linear)
   * @param[in]   outerAnsatzType   ansatz type of the outer integral (constant/linear)
   * @param[in]   leftClusterTree   cluster tree associated with row part of the mesh
   * @param[in]   rightClusterTree  cluster tree associated with column spart of the mesh
   * @param[in]   eta               admissibility constant
   */
  FastBESpace(
      SurfaceMesh3D<LO, SC>* leftMesh,
      SurfaceMesh3D<LO, SC>* rightMesh,
      basisType ansatzFunctionType,
      basisType testFunctionType,
      Tree<BECluster<LO, SC>*, LO >* leftClusterTree,
      Tree<BECluster<LO, SC>*, LO >* rightClusterTree,
      SC eta = 1.0,
  //    bool groupAdmClusters = true,
      LO maxElems = 0
      );

  /*!
   * constructor taking a mesh and ansatzs types, left=right cluster tree,
   * sets (possibly default) fast BEM parameters, creates a block cluster tree, and
   * creates (non)admissible pairs 
   * @param[in]   mesh
   * @param[in]   innerAnsatzType   ansatz type of the inner integral (constant/linear)
   * @param[in]   outerAnsatzType   ansatz type of the outer integral (constant/linear)
   * @param[in]   clusterTree       cluster tree associated with row and colum part of the mesh
   * @param[in]   eta               admissibility constant
   */
  FastBESpace(
      SurfaceMesh3D<LO, SC>* mesh,
      basisType ansatzFunctionType,
      basisType testFunctionType,
      Tree<BECluster<LO, SC>*, LO >* clusterTree,
      SCVT eta = 1.0,
      int nMax = 6,
      int quadOrder = 7,
      int startLevel = 0,
//      bool groupAdmClusters = true,
      LO maxElems = 0
      );

  //! destructor
  virtual ~FastBESpace( );

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

  //! sets the maximum number of elements allowed in block cluster

  inline void setMaxAggElemsPerCluster( LO maxElems ) {
    this->maxAggregatedElemsPerCluster = maxElems;
  }

  //! returns the maximum number of elements allowed in block cluster

  inline LO getMaxAggElemsPerCluster( ) {
    return this->maxAggregatedElemsPerCluster;
  }

  //! returns a tree associated with rows of a matrix

  inline Tree<BECluster<LO, SC>*, LO >* getLeftClusterTree( ) {
    return leftClusterTree;
  }

  //! returns a tree associated with columns of a matrix

  inline Tree<BECluster<LO, SC>*, LO >* getRightClusterTree( ) {
    return rightClusterTree;
  }

  //! returns pointer to the list of admissible blocks

  inline std::vector<BEBlockCluster<LO, SC>*>* getAdmissibleLeaves( ) {
    return &( this->admissibleLeaves );
  }

  //! returns pointer to the list of non-admissible blocks

  inline std::vector<BEBlockCluster<LO, SC>*>* getNonadmissibleLeaves( ) {
    return &( this->nonadmissibleLeaves );
  }

  //! method sets up block cluster tree
  void createBlockClusterTree( );

  /*
  //! collects a vector containing distribution of elements in clusters
  void elems2Clusters(
      std::vector<LO> &clusters
      ) const;
   */

protected:


  //! cluster tree associated with rows of matrix
  Tree<BECluster<LO, SC>*, LO >* leftClusterTree;

  //! cluster tree associated with columns of matrix
  Tree<BECluster<LO, SC>*, LO >* rightClusterTree;

  //! tree of boudnary element cluster pairs
  Tree<BEBlockCluster<LO, SC>*, LO >* blockClusterTree;

  //! method does recursively the actual block cluster constuction
  void doCreateBlockClusterTree( TreeMember<BECluster<LO, SC>*>* leftPointer, TreeMember<BECluster<LO, SC>*>* rightPointer );

  //! the method collects all leaves of block cluster tree and sorts them to admissible and nonadmissible
  void collectLeaves( );

  //! addmissibility constant
  SCVT eta;

  //! default maximum number of elements per cluster
  int defaultMaxElemsPerCluster = 15;

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
//  bool groupAdmClusters;

  //! maximum number of elements in aggregated clusters
  LO maxAggregatedElemsPerCluster;


  //! list of nonadmissible leaves
  std::vector<BEBlockCluster<LO, SC>*> nonadmissibleLeaves;

  //! list of admissible leaves
  std::vector<BEBlockCluster<LO, SC>*> admissibleLeaves;

  void doElems2Clusters(
      std::vector<LO> &clusters,
      LO & currentCluster
      ) const;

  void decomposeMesh( );

private:

};
}

// include .cpp file to overcome linking problems due to templates
#include "FastBESpace.cpp"

#endif	/* FASTBESPACE_H */

