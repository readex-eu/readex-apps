/*!
 * @file    Laplace1LayerP0P0MultilvlPrecond.cpp
 * @author  Lukas Maly
 * @date    June 20, 2016
 */

//#include "Laplace1LayerP0P0MultilvlPrecond.h"

#ifdef LAPLACE1LAYERP0P0MULTILVLPRECOND_H

namespace bem4i {

template<class LO, class SC>
Laplace1LayerP0P0MultilvlPrecond<LO, SC>::
Laplace1LayerP0P0MultilvlPrecond(
    SurfaceMesh3D<LO, SC> * mesh,
    LO maxElems
    ) {
  this->mesh = mesh;
  this->maxElemsPerLeave = maxElems;

  this->dimDomain = mesh->getNElements( );
  this->dimRange = mesh->getNElements( );

  this->init( );

}

template<class LO, class SC>
void Laplace1LayerP0P0MultilvlPrecond<LO, SC>::init( ) 
{
  //! Do nested dissection/clustering
  mesh->nestedDissection( this->tree, this->maxElemsPerLeave );
  
  //BECluster<LO, SC> *becluster = tree.get();
  
  /*
  //! Print leave nodes 
  std::vector< LO > cluster( mesh->getNElements(), 0 ); 
  tree.printLeaves( cluster );

  Vector< LO, SC > leaves( cluster.size( ) );
  SC * leavesData = leaves.getData( );
  for( LO i = 0; i < cluster.size(); ++i )
    leaves.set( i , SC(cluster[i]) );  

  std::vector< std::string > veName{ "leaves" };
  std::vector< Vector< LO, SC > * > veData{ &leaves };
  
  mesh->printParaviewVtu( "output/clustering.vtu", nullptr, nullptr, &veName,
    &veData, nullptr, nullptr, nullptr, nullptr );
  */

  this->endLevel = tree.getMinDepth( tree.getRootNode() );
  LO maxDepth = tree.getMaxDepth( tree.getRootNode() );
  
  // TODO: predelat zpatky na pole poli, vektor vektoru
  this->projectionMappingElem2Cluster.resize( maxDepth+1 );
  for(LO i = 0; i < this->projectionMappingElem2Cluster.size(); ++i )
    this->projectionMappingElem2Cluster[i] = 
      new Vector<LO,SC>( mesh->getNElements( ) );

  std::vector< std::vector<SC> > hj( maxDepth+1 );
  mesh->getNestedProjectionMassRatio( tree, 
    this->projectionMappingElem2Cluster );
  mesh->getNestedProjectionDiam( tree, hj ); 

  //std::cout << "clusterDiameters: " << endl;
  //for( LO i = 0; i < hj.size(); i++ ) {
  //  for( LO j = 0; j < hj[i].size(); j++ ) {
  //    std::cout << hj[i][j] << " ";
  //  }
  //  std::cout << std::endl;
  //}
  
  this->startLevel = -1;
  //! use maximum on each level of clDiams as h_j within the preconditioner
  this->clusterLevelDiam.resize(maxDepth+1);
  std::cout << "clusterLevelDiam.size = " << maxDepth+1 << std::endl;
  for( LO i = 0; i < clusterLevelDiam.size(); i++)
    this->clusterLevelDiam[i] =  0.5 * (*std::max_element( hj[i].begin() , 
        hj[i].end() ) + *std::min_element( hj[i].begin() , hj[i].end() ));
  this->clusterLevelDiam.push_back( mesh->getDiscretizationParameter( ) );

  /*
  //! diameter on next level set to be a half on the previous one
  for (LO i = 1; i < clusterLevelDiam.size(); i++)
  {
    clusterLevelDiam[i] = clusterLevelDiam[i-1]/2.0;
  } 
  */

  for( LO i = 0; i < clusterLevelDiam.size()-1; i++) {
    if( this->clusterLevelDiam[i] < 1 ) { 
      this->startLevel = i;
      break;
    }
  }
  
  if( this->startLevel > 0 || this->startLevel == -1 ) {
    std::cout << 
    "Diameter of coarse levels is over 1, the preconditioner might be unefficient!" <<
    std::endl;

    for (LO i = 0; i < this->endLevel; i++)
    {
      std::cout << "level " << i << ": " << clusterLevelDiam[i] << std::endl;
    }
 
    // zacina pokazde odshora
    this->startLevel = 0;
  }

  std::cout << "starlevel = " << this->startLevel << std::endl;
  std::cout << "endlevel = " << this->endLevel << std::endl;

  if( this->startLevel == -1 ) {
    std::cout << 
    "Mesh needs to be scaled, diameter on all levels is greater than 1!" <<
    std::endl;
    return;
  }

  if( this->startLevel > this-> endLevel ) {
    std::cout << "Your mesh is not suitable for this preconditioner, ";
    std::cout << "try to use finer clustering ... " << std::endl;
    return;
  }
}

//! implements aplication of y = C^{-1}*x = M^{-1} * A^{-1/2} * M^{-1} * x
template<class LO, class SC> 
void Laplace1LayerP0P0MultilvlPrecond< LO, SC >::apply(
      Vector< LO, SC > const & x,
      Vector< LO, SC > & y,
      bool transA,
      SC alpha,
      SC beta
      ) {

  Vector< LO, SC > yCopy;
  Vector< LO, SC > z( x );
  SC * zData = z.getData( );
  SC * yData = y.getData( );

  if ( std::abs( beta ) > 0.0 ) {
    y.copy( yCopy );
  }
  y.setAll( 0.0 );

  tree.setPointerToRoot( );

  // M^{-1} * x * alpha
  for ( LO i = 0; i < this->mesh->getNElements( ); ++i ) {
    zData[ i ] *= alpha / this->mesh->getElemArea( i );
  }
  
  // A^{-1/2} * x :
  applyProjection( tree, z, y, transA );

  // M^{-1} * x
  for ( LO i = 0; i < this->mesh->getNElements( ); ++i ) {
    yData[ i ] /= this->mesh->getElemArea( i );
  }

  if ( std::abs( beta ) > (SCVT) 0.0 ) {
    y.add( yCopy, beta );
  }
}

//! [Aw]_k = < sum_j hj*(Q_j - Q_j-1)*sum_i*w_i*phi_i , phi_k >
template<class LO, class SC> 
void Laplace1LayerP0P0MultilvlPrecond< LO, SC >::applyProjection(
      Tree<BECluster<LO, SC>*, LO > & clusterTree,
      Vector< LO, SC > const & x,
      Vector< LO, SC > & y,
      bool transA
      ) {

    BECluster<LO, SC> *currCluster;
    LO currLevel = clusterTree.getLevel( );
    if ( currLevel > this->endLevel ) return;
    //std::cout << "level: " << currLevel << std::endl;
    LO idx;
    SC diff_h;
    SC res = 0.0;

    // -Q_0 * h_0 * w_i * phi_i
    if( startLevel > 0 && currLevel == (startLevel-1) ) {
      currCluster = clusterTree.get( );
      //std::cout << " ----- extra level bellow start level ----- " << std::endl;
      //diff_h = (-1.0) * (1.0/(this->clusterLevelDiam[currLevel + 1]));
      diff_h = (-1.0) * (1.0/std::sqrt(this->clusterLevelDiam[currLevel + 1])); 
      //diff_h = (-1.0) * this->clusterLevelDiam[currLevel + 1];
    
      for (LO i = 0; i < currCluster->nelems; ++i) {
        idx = currCluster->elems->at( i );
        //std::cout << idx << ",";
        res += projectionMappingElem2Cluster[ currLevel ]->get( idx ) * 
          x.get( idx );
      }
      //std::cout << std::endl;
      //std::cout << "res = " << res << std::endl;
      for (LO i = 0; i < currCluster->nelems; ++i) {
        idx = currCluster->elems->at( i );
        y.add( idx, res * diff_h * mesh->getElemArea( idx ) );
      }
    }

    if( currLevel >= startLevel && currLevel <= endLevel ) {
      currCluster = clusterTree.get( );
      if( currLevel == endLevel ) {
        //diff_h = 1.0/(this->clusterLevelDiam[currLevel]) 
        //      - 1.0/(this->clusterLevelDiam.back( ));
        diff_h = 1.0/std::sqrt(this->clusterLevelDiam[currLevel]) 
              - 1.0/std::sqrt(this->clusterLevelDiam.back( ));
        //diff_h = this->clusterLevelDiam[currLevel] 
        //      - this->clusterLevelDiam.back( );
      } else {
        //diff_h = 1.0/(this->clusterLevelDiam[currLevel]) 
        //      - 1.0/(this->clusterLevelDiam[currLevel + 1]);
        diff_h = 1.0/std::sqrt(this->clusterLevelDiam[currLevel]) 
              - 1.0/std::sqrt(this->clusterLevelDiam[currLevel + 1]);
        //diff_h = this->clusterLevelDiam[currLevel] 
        //      - this->clusterLevelDiam[currLevel + 1];
      } 
    
      for (LO i = 0; i < currCluster->nelems; ++i) {
        idx = currCluster->elems->at( i );
        //std::cout << idx << ", ";
        res += projectionMappingElem2Cluster[ currLevel ]->get( idx ) * 
          x.get( idx );
      }
      //std::cout << std::endl;
      //std::cout << "res = " << res << std::endl;
      for (LO i = 0; i < currCluster->nelems; ++i) {
        idx = currCluster->elems->at( i );
        y.add( idx, res * diff_h * mesh->getElemArea( idx ) );
      }     
    }

    // Q_J * h_J * w_i * phi_i ... virtual element level
    if( currLevel == endLevel ) {
      currCluster = clusterTree.get( );
      //diff_h = 1.0/(this->clusterLevelDiam.back( )); 
      diff_h = 1.0/std::sqrt(this->clusterLevelDiam.back( )); 
      //diff_h = this->clusterLevelDiam.back( );
    
      for (LO i = 0; i < currCluster->nelems; ++i) {
        idx = currCluster->elems->at( i );
        y.add( idx, x.get( idx ) * diff_h * mesh->getElemArea( idx ) );
      }
    }

    LO nSons = clusterTree.getNSons( );
    for (LO i = 0; i < nSons; ++i) {
      clusterTree.setPointerToSon( i );
      applyProjection( clusterTree, x, y, transA );
      clusterTree.setPointerToParent( );
    }
  }

/*
#ifdef BLAS_INT
template class Laplace1LayerP0P0MultilvlPrecond<int, double>;
template class Laplace1LayerP0P0MultilvlPrecond<int, float>;
#endif

#ifdef BLAS_LONG
template class Laplace1LayerP0P0MultilvlPrecond<long, double>;
template class Laplace1LayerP0P0MultilvlPrecond<long, float>;
#endif
*/

} // end of namespace

#endif
