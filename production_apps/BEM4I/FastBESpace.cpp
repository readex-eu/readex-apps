/*!
 * @file    FastBESpace.cpp
 * @author  Michal Merta 
 * @date    August 14, 2013
 * 
 */

#ifdef FASTBESPACE_H

namespace bem4i {

template<class LO, class SC>
FastBESpace<LO, SC>::FastBESpace( ) {
}

template<class LO, class SC>
FastBESpace<LO, SC>::FastBESpace(
    const FastBESpace & orig
    ) {
}

template<class LO, class SC>
FastBESpace<LO, SC>::~FastBESpace( ) {
  delete blockClusterTree;
}

template<class LO, class SC>
FastBESpace<LO, SC>::FastBESpace(
    SurfaceMesh3D<LO, SC>* mesh,
    basisType ansatzFunctionType,
    basisType testFunctionType,
    SC eta,
    int nMax,
    int quadOrder,
    int startLevel,
    //    bool groupAdmClusters,
    LO maxElems
    ) {
  this->mesh = this->leftMesh = this->rightMesh = mesh;
  this->ansatzFunctionType = ansatzFunctionType;
  this->testFunctionType = testFunctionType;

  this->nInnerElems = mesh->getNElements( );
  this->nOuterElems = mesh->getNElements( );

  this->innerElems.resize( this->nInnerElems );
  this->outerElems.resize( this->nOuterElems );

  for ( int i = 0; i < this->nInnerElems; i++ ) {
    this->innerElems[i] = i;
    this->outerElems[i] = i;
  }

  this->totalDOFsPerInnerElem = basisType2functions[ansatzFunctionType];

  this->totalDOFsPerOuterElem = basisType2functions[testFunctionType];

  setFastBEMParams( eta, nMax, quadOrder, startLevel );

  // set default values for ACA
  this->epsilonACA = 1e-4;
  this->scaleACA = 1.0;

  //  this->groupAdmClusters = groupAdmClusters;
  if ( maxElems == 0 ) {
    this->maxAggregatedElemsPerCluster = this->nInnerElems;
  } else {
    this->maxAggregatedElemsPerCluster = maxElems;
  }

  // create a cluster tree and block cluster tree for fast bem
  this->leftClusterTree = new Tree<BECluster<LO, SC>*, LO >( );
  this->rightClusterTree = this->leftClusterTree;
  this->mesh->nestedDissection( *this->leftClusterTree, defaultMaxElemsPerCluster );

  // create a mapping between elements and appropriate degrees of freedom in global matrix
  //this->defineMappingToDOFs( );

  this->createBlockClusterTree( );
}

template<class LO, class SC>
FastBESpace<LO, SC>::FastBESpace(
    SurfaceMesh3D<LO, SC>* mesh,
    basisType ansatzFunctionType,
    basisType testFunctionType,
    Tree<BECluster<LO, SC>*, LO >* clusterTree,
    SCVT eta,
    int nMax,
    int quadOrder,
    int startLevel,
    //    bool groupAdmClusters,
    LO maxElems
    ) {
  this->mesh = this->leftMesh = this->rightMesh = mesh;
  this->ansatzFunctionType = ansatzFunctionType;
  this->testFunctionType = testFunctionType;

  this->leftClusterTree = clusterTree;
  this->rightClusterTree = clusterTree;

  this->nInnerElems = mesh->getNElements( );
  this->nOuterElems = mesh->getNElements( );

  this->innerElems.resize( this->nInnerElems );
  this->outerElems.resize( this->nOuterElems );

  for ( int i = 0; i < this->nInnerElems; i++ ) {
    this->innerElems[i] = i;
    this->outerElems[i] = i;
  }

  this->totalDOFsPerInnerElem = basisType2functions[ansatzFunctionType];

  this->totalDOFsPerOuterElem = basisType2functions[testFunctionType];

  setFastBEMParams( eta, nMax, quadOrder, startLevel );

  // set default values for ACA
  this->epsilonACA = 1e-4;
  this->scaleACA = 1.0;

  //  this->groupAdmClusters = groupAdmClusters;
  if ( maxElems == 0 ) {
    this->maxAggregatedElemsPerCluster = this->nInnerElems;
  } else {
    this->maxAggregatedElemsPerCluster = maxElems;
  }

  // create a mapping between elements and appropriate degrees of freedom in global matrix
  //this->defineMappingToDOFs( );

  this->createBlockClusterTree( );


}

template<class LO, class SC>
FastBESpace<LO, SC>::FastBESpace(
    SurfaceMesh3D<LO, SC>* mesh,
    basisType ansatzFunctionType,
    basisType testFunctionType,
    Tree<BECluster<LO, SC>*, LO >* leftClusterTree,
    Tree<BECluster<LO, SC>*, LO >* rightClusterTree,
    SC eta,
    int nMax,
    int quadOrder,
    int startLevel,
    //    bool groupAdmClusters,
    LO maxElems
    ) {
  this->mesh = this->leftMesh = this->rightMesh = mesh;
  this->ansatzFunctionType = ansatzFunctionType;
  this->testFunctionType = testFunctionType;

  this->leftClusterTree = leftClusterTree;
  this->rightClusterTree = rightClusterTree;

  this->nInnerElems = mesh->getNElements( );
  this->nOuterElems = mesh->getNElements( );

  this->innerElems.resize( this->nInnerElems );
  this->outerElems.resize( this->nOuterElems );

  for ( int i = 0; i < this->nInnerElems; i++ ) {
    this->innerElems[i] = i;
    this->outerElems[i] = i;
  }

  this->totalDOFsPerInnerElem = basisType2functions[ansatzFunctionType];

  this->totalDOFsPerOuterElem = basisType2functions[testFunctionType];

  setFastBEMParams( eta, nMax, quadOrder, startLevel );

  //  this->groupAdmClusters = groupAdmClusters;
  if ( maxElems == 0 ) {
    this->maxAggregatedElemsPerCluster = this->nInnerElems;
  } else {
    this->maxAggregatedElemsPerCluster = maxElems;
  }

  // create a mapping between elements and appropriate degrees of freedom in global matrix
  //this->defineMappingToDOFs( );

  this->createBlockClusterTree( );
}

template<class LO, class SC>
FastBESpace<LO, SC>::FastBESpace(
    SurfaceMesh3D<LO, SC>* leftMesh,
    SurfaceMesh3D<LO, SC>* rightMesh,
    basisType ansatzFunctionType,
    basisType testFunctionType,
    Tree<BECluster<LO, SC>*, LO >* leftClusterTree,
    Tree<BECluster<LO, SC>*, LO >* rightClusterTree,
    SC eta,
    //    bool groupAdmClusters,
    LO maxElems
    ) {
  this->leftMesh = leftMesh;
  this->rightMesh = rightMesh;

  this->ansatzFunctionType = ansatzFunctionType;
  this->testFunctionType = testFunctionType;

  this->leftClusterTree = leftClusterTree;
  this->rightClusterTree = rightClusterTree;

  this->nInnerElems = rightMesh->getNElements( );
  this->nOuterElems = leftMesh->getNElements( );

  this->innerElems.resize( this->nInnerElems );
  this->outerElems.resize( this->nOuterElems );

  for ( int i = 0; i < this->nInnerElems; i++ ) {
    this->innerElems[i] = i;
  }

  for ( int i = 0; i < this->nOuterElems; i++ ) {
    this->outerElems[i] = i;
  }

  this->totalDOFsPerInnerElem = basisType2functions[ansatzFunctionType];

  this->totalDOFsPerOuterElem = basisType2functions[testFunctionType];

  this->eta = eta;

  //  this->groupAdmClusters = groupAdmClusters;
  if ( maxElems == 0 ) {
    this->maxAggregatedElemsPerCluster =
        ( this->nInnerElems + this->nOuterElems ) / 2;
  } else {
    this->maxAggregatedElemsPerCluster = maxElems;
  }


  // create a mapping between elements and appropriate degrees of freedom in global matrix
  //this->defineMappingToDOFs( );

  this->createBlockClusterTree( );
}

template<class LO, class SC>
void FastBESpace<LO, SC>::createBlockClusterTree( ) {

  // create tree of cluster pairs
  blockClusterTree = new Tree< BEBlockCluster< LO, SC > *, LO >( );
  BEBlockCluster< LO, SC > * root = new BEBlockCluster< LO, SC >( );
  blockClusterTree->createRoot( root );
  blockClusterTree->setPointerToRoot( );

  // start with roots of both trees
  TreeMember< BECluster< LO, SC > * > * leftPointer =
      leftClusterTree->getRootNode( );
  TreeMember< BECluster< LO, SC > * > * rightPointer =
      rightClusterTree->getRootNode( );

  BECluster<LO, SC> * leftCluster = leftPointer->data;
  BECluster<LO, SC> * rightCluster = rightPointer->data;
  root->leftCluster = leftCluster;
  root->rightCluster = rightCluster;

  doCreateBlockClusterTree( leftPointer, rightPointer );

  blockClusterTree->setPointerToRoot( );
  collectLeaves( );
}

template<class LO, class SC>
void FastBESpace<LO, SC>::doCreateBlockClusterTree(
    TreeMember< BECluster< LO, SC > * > * leftPointer,
    TreeMember< BECluster< LO, SC > * > * rightPointer
    ) {

  BEBlockCluster< LO, SC > * newBlockCluster;
  int idx = 0;

  BEBlockCluster< LO, SC > * currentCluster = blockClusterTree->get( );
  BECluster< LO, SC > * leftCluster = currentCluster->leftCluster;
  BECluster< LO, SC > * rightCluster = currentCluster->rightCluster;

  // compute distance between centroids
  SCVT distance = 0.0;
  distance = std::sqrt(
      ( leftCluster->centroid[ 0 ] - rightCluster->centroid[ 0 ] ) *
      ( leftCluster->centroid[ 0 ] - rightCluster->centroid[ 0 ] ) +
      ( leftCluster->centroid[ 1 ] - rightCluster->centroid[ 1 ] ) *
      ( leftCluster->centroid[ 1 ] - rightCluster->centroid[ 1 ] ) +
      ( leftCluster->centroid[ 2 ] - rightCluster->centroid[ 2 ] ) *
      ( leftCluster->centroid[ 2 ] - rightCluster->centroid[ 2 ] ) );

  if ( 2.0 * std::min( leftCluster->radius, rightCluster->radius ) <=
      eta * ( distance - leftCluster->radius - rightCluster->radius )
      && ( ( leftCluster->nelems + rightCluster->nelems ) / 2 ) <
      maxAggregatedElemsPerCluster ) {
    //&& (groupAdmClusters || leftClusterTree->getNSons( leftPointer ) == 0 || rightClusterTree->getNSons( rightPointer ) == 0) ) {
    //      std::cout << leftCluster->centroid[0] << ", " << leftCluster->centroid[1] << ", " << leftCluster->centroid[2] << std::endl;
    //std::cout << rightCluster->centroid[0] << ", " << rightCluster->centroid[1] << ", " << rightCluster->centroid[2] << std::endl;
    //std::cout << "dist: " << distance << ", left r.: " << leftCluster->radius << ", right r.: " << rightCluster->radius << std::endl;
    //std::cout<<2 * std::min(leftCluster->radius, rightCluster->radius) << ", " <<eta * (distance - leftCluster->radius - rightCluster->radius) << std::endl << std::endl;;
    currentCluster->admissible = true;
    currentCluster->leaf = true;

    // add right cluster to the list of admissible clusters of the left cluster
    if ( leftCluster->admClusters == nullptr ) {
      leftCluster->admClusters = new std::vector< BECluster< LO, SC > * >( );
    }
    leftCluster->admClusters->push_back( rightCluster );
  } else {

    currentCluster->admissible = false;
    int nLeftSons = leftClusterTree->getNSons( leftPointer );
    int nRightSons = rightClusterTree->getNSons( rightPointer );
    int nNonadmissibleLeaves = 0;
    int nAdmissibleLeaves = 0;

    for ( int i = 0; i < nLeftSons; i++ ) {
      leftPointer = leftClusterTree->getSonNode( leftPointer, i );
      for ( int j = 0; j < nRightSons; j++, idx++ ) {
        rightPointer = rightClusterTree->getSonNode( rightPointer, j );
        newBlockCluster = new BEBlockCluster< LO, SC >( );
        newBlockCluster->leftCluster = leftClusterTree->get( leftPointer );
        newBlockCluster->rightCluster = rightClusterTree->get( rightPointer );
        blockClusterTree->addSon( newBlockCluster );
        blockClusterTree->setPointerToSon( idx );
        doCreateBlockClusterTree( leftPointer, rightPointer );

        if ( !newBlockCluster->admissible && newBlockCluster->leaf ) {
          nNonadmissibleLeaves++;
        }

        if ( newBlockCluster->admissible && newBlockCluster->leaf ) {
          nAdmissibleLeaves++;
        }

        blockClusterTree->setPointerToParent( );
        rightPointer = rightClusterTree->getParentNode( rightPointer );
      }
      leftPointer = leftClusterTree->getParentNode( leftPointer );
    }

    if ( idx == 0 ) {

      currentCluster->leaf = true;

    } else if ( nNonadmissibleLeaves == nLeftSons * nRightSons
        //        && groupAdmClusters
        && ( ( ( leftCluster->nelems + rightCluster->nelems ) / 2 ) <
        maxAggregatedElemsPerCluster )
        ) {

      currentCluster->leaf = true;

    } else if ( nAdmissibleLeaves == nLeftSons * nRightSons
        //        && groupAdmClusters
        && ( ( ( leftCluster->nelems + rightCluster->nelems ) / 2 ) <
        maxAggregatedElemsPerCluster ) ) {

      currentCluster->leaf = true;
      currentCluster->admissible = true;

    } else {

      currentCluster->leaf = false;

    }
  }
}

//  template<class LO, class SC>
//  void FastBESpace<LO, SC>::detectFarfieldOnAllLevels(TreeMember<BECluster<LO, SC>*>* leftPointer, TreeMember<BECluster<LO, SC>*>* rightPointer) {
//    
//    BECluster<LO, SC>* leftCluster = leftClusterTree->get(leftPointer);
//    BECluster<LO, SC>* rightCluster = rightClusterTree->get(leftPointer);
//    
//    double distance = std::sqrt((leftCluster->centroid[0] - rightCluster->centroid[0]) *
//            (leftCluster->centroid[0] - rightCluster->centroid[0]) +
//            (leftCluster->centroid[1] - rightCluster->centroid[1]) *
//            (leftCluster->centroid[1] - rightCluster->centroid[1]) +
//            (leftCluster->centroid[2] - rightCluster->centroid[2]) *
//            (leftCluster->centroid[2] - rightCluster->centroid[2]));
//    
//    if (2 * std::min(leftCluster->radius, rightCluster->radius) <=
//            eta * (distance - leftCluster->radius - rightCluster->radius)) {
//      if (leftCluster->admClusters == nullptr) {
//        leftCluster->admClusters = new std::vector<BECluster<LO, SC>* >();
//      }
//      
//      leftCluster->admClusters->push_back(rightCluster);
//    } else {
//      for (int i = 0; i < leftClusterTree->getNSons(leftPointer); i++) {
//        leftPointer = leftClusterTree->getSonNode(leftPointer, i);
//        for (int j = 0; j < rightClusterTree->getNSons(rightPointer); j++) {
//          rightPointer = rightClusterTree->getSonNode(rightPointer, j);
//          detectFarfieldOnAllLevels(leftPointer, rightPointer);
//          rightPointer = rightClusterTree->getParentNode(rightPointer);
//        }
//        leftPointer = leftClusterTree->getParentNode(leftPointer);
//      }
//    }
//  }

template<class LO, class SC>
void FastBESpace<LO, SC>::collectLeaves( ) {
  int nSons;
  BEBlockCluster<LO, SC>* blockCluster;

  blockCluster = blockClusterTree->get( );

  // set degrees of freedom for left and right cluster
  if ( blockCluster->leaf ) {

    blockCluster->innerDOFs = new std::vector<LO>;
    blockCluster->outerDOFs = new std::vector<LO>;

    // for right cluster
    if ( this->ansatzFunctionType == p0 ) {
      std::copy( blockCluster->rightCluster->elems->begin( ),
          blockCluster->rightCluster->elems->end( ),
          std::back_inserter( *blockCluster->innerDOFs ) );
    } else if ( this->ansatzFunctionType == p1 ) {
      std::copy( blockCluster->rightCluster->nodes->begin( ),
          blockCluster->rightCluster->nodes->end( ),
          std::back_inserter( *blockCluster->innerDOFs ) );
    } else if ( this->ansatzFunctionType == p1dis ) {
      blockCluster->innerDOFs->clear( );
      blockCluster->innerDOFs->reserve(
          3 * blockCluster->rightCluster->nelems );
      for ( LO i = 0; i < blockCluster->rightCluster->nelems; ++i ) {
        blockCluster->innerDOFs->push_back(
            3 * blockCluster->rightCluster->elems->at( i ) );
        blockCluster->innerDOFs->push_back(
            3 * blockCluster->rightCluster->elems->at( i ) + 1 );
        blockCluster->innerDOFs->push_back(
            3 * blockCluster->rightCluster->elems->at( i ) + 2 );
      }
    }

    // left cluster
    if ( this->testFunctionType == p0 ) {
      std::copy( blockCluster->leftCluster->elems->begin( ),
          blockCluster->leftCluster->elems->end( ),
          std::back_inserter( *blockCluster->outerDOFs ) );
    } else if ( this->testFunctionType == p1 ) {
      std::copy( blockCluster->leftCluster->nodes->begin( ),
          blockCluster->leftCluster->nodes->end( ),
          std::back_inserter( *blockCluster->outerDOFs ) );
    } else if ( this->testFunctionType == p1dis ) {
      blockCluster->outerDOFs->clear( );
      blockCluster->outerDOFs->reserve(
          3 * blockCluster->leftCluster->nelems );
      for ( LO i = 0; i < blockCluster->leftCluster->nelems; ++i ) {
        blockCluster->outerDOFs->push_back(
            3 * blockCluster->leftCluster->elems->at( i ) );
        blockCluster->outerDOFs->push_back(
            3 * blockCluster->leftCluster->elems->at( i ) + 1 );
        blockCluster->outerDOFs->push_back(
            3 * blockCluster->leftCluster->elems->at( i ) + 2 );
      }
    }
    /*
    std::cout << blockCluster->leftCluster->elems->size( ) << " "
        << blockCluster->rightCluster->elems->size( ) << " "
        << blockCluster->admissible << std::endl;
    */
  }

  if ( !blockCluster->admissible && blockCluster->leaf ) {
    nonadmissibleLeaves.push_back( blockCluster );
  } else if ( blockCluster->admissible && blockCluster->leaf ) {
    admissibleLeaves.push_back( blockCluster );
  } else {
    nSons = blockClusterTree->getNSons( );
    for ( int i = 0; i < nSons; i++ ) {
      blockClusterTree->setPointerToSon( i );
      collectLeaves( );
      blockClusterTree->setPointerToParent( );
    }
  }
}

/*
template<class LO, class SC>
void FastBESpace<LO, SC>::elems2Clusters(
    std::vector<LO> &clusters
    ) const {
  clusters.clear( );
  clusters.resize( this->mesh->getNElements( ), -1 );
  LO currentCluster = 0;
  doElems2Clusters( clusters, currentCluster );
}

template<class LO, class SC>
void FastBESpace<LO, SC>::doElems2Clusters(
    std::vector<LO> &clusters,
    LO & currentCluster
    ) const {

  BECluster<LO, SC>* cluster;
  std::vector<LO>* clusterElems;

  cluster = leftClusterTree->get( );
  int nSons = leftClusterTree->getNSons( );

  if ( nSons == 0 ) {
    //std::cout << blockCluster->leftCluster->nelems << std::endl;
    //std::cout << blockCluster->rightCluster->nelems << std::endl;
    clusterElems = cluster->elems;
    if ( clusters[( *clusterElems )[0]] == -1 ) {
      for ( auto it = clusterElems->begin( );
          it != clusterElems->end( ); ++it ) {
        clusters[ *it ] = currentCluster;
      }
      ++currentCluster;
    }
  } else {
    nSons = leftClusterTree->getNSons( );
    for ( int i = 0; i < nSons; i++ ) {
      leftClusterTree->setPointerToSon( i );
      doElems2Clusters( clusters, currentCluster );
      leftClusterTree->setPointerToParent( );
    }
  }

  //  int nSons;
  //  BEBlockCluster<LO, SC>* blockCluster;
  //  std::vector<LO>* clusterElems;
  //
  //  blockCluster = blockClusterTree->get( );
  //
  //  if ( blockCluster->leaf ) {
  //    //std::cout << blockCluster->leftCluster->nelems << std::endl;
  //    //std::cout << blockCluster->rightCluster->nelems << std::endl;
  //    clusterElems = blockCluster->leftCluster->elems;
  //    if ( clusters[( *clusterElems )[0]] == -1 ) {
  //      for ( auto it = clusterElems->begin( );
  //          it != clusterElems->end( ); ++it ) {
  //        clusters[ *it ] = currentCluster;
  //      }
  //      ++currentCluster;
  //    }
  //    
  //    clusterElems = blockCluster->rightCluster->elems;
  //    if ( clusters[( *clusterElems )[0]] == -1 ) {
  //      for ( auto it = clusterElems->begin( );
  //          it != clusterElems->end( ); ++it ) {
  //        clusters[ *it ] = currentCluster;
  //      }
  //      ++currentCluster;
  //    }
  //  } else {
  //    nSons = blockClusterTree->getNSons( );
  //    for ( int i = 0; i < nSons; i++ ) {
  //      blockClusterTree->setPointerToSon( i );
  //      doElems2Clusters( clusters, currentCluster );
  //      blockClusterTree->setPointerToParent( );
  //    }
  //  }
}
 */

}

#endif

