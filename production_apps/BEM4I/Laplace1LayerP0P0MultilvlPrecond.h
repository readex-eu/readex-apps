/*!
 * @file    Laplace1LayerP0P0MultilvlPrecond.h
 * @author  Lukas Maly 
 * @date    July 20, 2016
 * @brief   Header file for class Laplace1LayerP0P0MultilvlPrecond
 * 
 */


#ifndef LAPLACE1LAYERP0P0MULTILVLPRECOND_H
#define LAPLACE1LAYERP0P0MULTILVLPRECOND_H

#include <vector>

#include "LinearOperator.h"
#include "Vector.h"
#include "Tree.h"
#include "BESpace.h"

namespace bem4i {

/*!
 * Class representing an Artificial Multilevel Preconditioner
 */
template<class LO, class SC>
class Laplace1LayerP0P0MultilvlPrecond : public LinearOperator<LO, SC> {
  
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! constructor taking a mesh as argument
 Laplace1LayerP0P0MultilvlPrecond(
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

  void applyProjection(
      Tree<BECluster<LO, SC>*, LO > & clusterTree,
      Vector< LO, SC > const & x,
      Vector< LO, SC > & y,
      bool transA = false
      );
 
  virtual  ~Laplace1LayerP0P0MultilvlPrecond( ) {
    for(LO i = 0; i < this->projectionMappingElem2Cluster.size(); ++i )
      delete this->projectionMappingElem2Cluster[i];
  };

private:
  
 Laplace1LayerP0P0MultilvlPrecond( ) {
    this->dimRange = 0;
    this->dimDomain = 0;
  }

  void init( );

  // Tree of BEClusters .. multilevel mesh clusters
  Tree<BECluster<LO, SC> *, LO > tree;

  // maximal number of elems within the tree-leaves
  LO maxElemsPerLeave; 

  // Mesh
  SurfaceMesh3D<LO, SC> * mesh; 

  // Starting level: diam of cluster < 1
  LO startLevel;

  // Ending level = minimal depth
  LO endLevel;

  // projectionMappingElem2Cluster from leave to root
  std::vector< Vector< LO, SC > * > projectionMappingElem2Cluster;  

  // Max diameters of cluster on each level
  std::vector< SC > clusterLevelDiam;


};

} //end of namespace
#include "Laplace1LayerP0P0MultilvlPrecond.cpp"

#endif /* LAPLACE1LAYERP0P0MULTILVLPRECOND_H */
