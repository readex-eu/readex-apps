/*!
 * @file    DistributedFastBESpace.cpp
 * @author  Michal Merta 
 * @date    February 17, 2016
 * 
 */

#ifdef DISTRIBUTEDFASTBESPACE_H

namespace bem4i {

template<class LO, class SC>
DistributedFastBESpace<LO, SC>::DistributedFastBESpace(
    ) {
}

template<class LO, class SC>
DistributedFastBESpace<LO, SC>::DistributedFastBESpace(
    const DistributedFastBESpace& orig
    ) {
}

template<class LO, class SC>
DistributedFastBESpace<LO, SC>::~DistributedFastBESpace( ) {
  for ( LO i = 0; i < nSubmeshes; i++ ) {
    if ( submeshes[i] != nullptr ) {
      delete submeshes[i];
    }

    if ( this->innerDOFs != nullptr && this->innerDOFs[i] != nullptr ) {
      delete this->innerDOFs[i];
    }
    if ( this->outerDOFs != nullptr && this->outerDOFs[i] != nullptr ) {
      delete this->outerDOFs[i];
    }
    if (this->outerDOFs != nullptr) {
      delete this->outerDOFs;
    }
    if (this->innerDOFs != nullptr) {
      delete this->innerDOFs;
    }

  }
  delete [] submeshes;

}

template<class LO, class SC>
DistributedFastBESpace<LO, SC>::DistributedFastBESpace(
    SurfaceMesh3D<LO, SC>* mesh,
    LO nSubmeshes,
    basisType ansatzFunctionType,
    basisType testFunctionType,
    SC eta,
    bool groupAdmClusters,
    LO maxElems
    ) {
  this->mesh = mesh;
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

  this->eta = eta;
  this->maxElemsPerCluster = defaultMaxElemsPerCluster;

  // set default values for ACA
  this->epsilonACA = 1e-4;
  this->scaleACA = 1.0;

  this->groupAdmClusters = groupAdmClusters;
  this->maxAggregatedElemsPerCluster = maxElems;

  this->nSubmeshes = nSubmeshes;

  // create array of pointers to new submeshes (create HERE by new and 
  // deallocated in destructor)
  submeshes = new SurfaceMesh3D<LO, SC>*[nSubmeshes];
  bespaces = new FastBESpace<LO, SC>*[nSubmeshes * nSubmeshes];
  trees = new Tree<BECluster<LO, SC>*, LO >*[nSubmeshes];
  this->outerDOFs = new std::vector<LO>*[nSubmeshes * nSubmeshes];
  this->innerDOFs = new std::vector<LO>*[nSubmeshes * nSubmeshes];

  for ( LO i = 0; i < nSubmeshes; ++i ) {
    submeshes[i] = nullptr;
    trees[i] = nullptr;
  }
  for ( LO i = 0; i < nSubmeshes * nSubmeshes; ++i ) {
    bespaces[i] = nullptr;
    this->outerDOFs[ i ] = nullptr;
    this->innerDOFs[ i ] = nullptr;
  }

  this->decomposeMesh();
  
}

template<class LO, class SC>
void DistributedFastBESpace<LO, SC>::decomposeMesh( ) {
  // call METIS and decompose mesh. store result in array submeshes.

  // ... fill here 

  // create trees and FastBESpaces 
  for ( LO i = 0; i < nSubmeshes; i++ ) {
    if ( trees[i] ) delete trees[i];
    trees[i] = new Tree<BECluster<LO, SC>*, LO>( );
    submeshes[i]->nestedDissection( *trees[i], maxElemsPerCluster );
  }
  
  for ( LO i = 0; i < nSubmeshes; i++ ) {
    for ( LO j = 0; j < nSubmeshes; j++ ) {
      if ( bespaces[j * nSubmeshes + i] ) delete bespaces[j * nSubmeshes + i];
      bespaces[i * nSubmeshes + j] = new FastBESpace<LO, SC>( submeshes[j], 
          submeshes[i], this->ansatzFunctionType, this->testFunctionType, 
          trees[j], trees[i], this->eta, this->groupAdmClusters );
    }
  }
}

}
#endif