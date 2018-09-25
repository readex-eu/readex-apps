/*!
 * @file    BESpace.cpp
 * @author  Michal Merta 
 * @date    July 9, 2013
 * 
 */



#ifdef BESPACE_H

namespace bem4i {

template<class LO, class SC>
BESpace<LO, SC>::BESpace( ) {
}

template<class LO, class SC>
BESpace<LO, SC>::BESpace(
    const BESpace & orig
    ) {
}

template<class LO, class SC>
BESpace<LO, SC>::BESpace(
    SurfaceMesh3D<LO, SC> * mesh,
    basisType ansatzFunctionType,
    basisType testFunctionType
    ) {

  this->rightMesh = this->leftMesh = this->mesh = mesh;
  this->ansatzFunctionType = ansatzFunctionType;
  this->testFunctionType = testFunctionType;

  nInnerElems = mesh->getNElements( );
  nOuterElems = mesh->getNElements( );

  innerElems.resize( nInnerElems );
  outerElems.resize( nOuterElems );

  for ( int i = 0; i < nInnerElems; i++ ) {
    innerElems[i] = i;
    outerElems[i] = i;
  }

  totalDOFsPerInnerElem = basisType2functions[ansatzFunctionType];

  totalDOFsPerOuterElem = basisType2functions[testFunctionType];

  // create a mapping between elements and appropriate degrees of freedom in global matrix
  //defineMappingToDOFs( );
}

// for now just assume we work with whole mesh

/*
template<class LO, class SC>
BESpace<LO, SC>::BESpace(SurfaceMesh3D<LO>& mesh, vector<LO> innerElems, vector<LO> outerElems,
          basisType innerAnsatz, basisType outerAnsatz) {
  this->mesh = mesh;
  this->innerAnsatz = innerAnsatz;
  this->outerAnsatz = outerAnsatz;
  
  nInnerElems = mesh->getNElements();
  nOuterElems = mesh->getNElements();
  
  innerElems.resize(nInnerElems);
  outerElems.resize(nOuterElems);
  
  this->innerElems = innerElems;
  this->outerElems = outerElems;
  
  totalDOFsPerInnerElem = basisType2functions[ansatzFunctionType];
  
  totalDOFsPerOuterElem = basisType2functions[testFunctionType];
  
  // create a mapping between elements and appropriate degrees of freedom in global matrix
  defineMappingToDOFs();
}
 */


template<class LO, class SC>
BESpace<LO, SC>::~BESpace( ) {
}

/*
template<class LO, class SC>
void BESpace<LO, SC>::defineMappingToDOFs( ) {

  innerIdx2DOFs.resize( basisType2functions[ansatzFunctionType] * nInnerElems );
  outerIdx2DOFs.resize( basisType2functions[testFunctionType] * nOuterElems );

  int idx;

  // for now, just set this only for scalar pw constant and linear functions
  // TODO later we need to extend this for general scalar/vector functions (but in some subclass)
  if ( ansatzFunctionType == p0 ) {
    memcpy( &innerIdx2DOFs[0], &innerElems[0], nInnerElems * sizeof (LO ) );
  } else if ( ansatzFunctionType == p1 ) {
    for ( int i = 0; i < nInnerElems; i++ ) {
      idx = innerElems[i];
      rightMesh->getElement( idx, 
          &innerIdx2DOFs[ i * rightMesh->getNNodesPerElem( )] );
    }
  } else if ( ansatzFunctionType == p1dis ) {
    for ( LO i = 0; i < nInnerElems; i++ ) {
      idx = innerElems[i];
      innerIdx2DOFs[3 * i] = 3 * idx;
      innerIdx2DOFs[3 * i + 1] = 3 * idx + 1;
      innerIdx2DOFs[3 * i + 2] = 3 * idx + 2;
    }
  }

  if ( testFunctionType == p0 ) {
    memcpy( &outerIdx2DOFs[0], &outerElems[0], nOuterElems * sizeof (LO ) );
  } else if ( testFunctionType == p1 ) {
    for ( int i = 0; i < nOuterElems; i++ ) {
      idx = outerElems[i];
      leftMesh->getElement( idx, 
          &outerIdx2DOFs[i * leftMesh->getNNodesPerElem( )] );
    }
  } else if ( testFunctionType == p1dis ) {
    for ( LO i = 0; i < nOuterElems; i++ ) {
      idx = outerElems[i];
      outerIdx2DOFs[3 * i] = 3 * idx;
      outerIdx2DOFs[3 * i + 1] = 3 * idx + 1;
      outerIdx2DOFs[3 * i + 2] = 3 * idx + 2;
    }
  }
}
 */

template<class LO, class SC>
void BESpace<LO, SC>::getOuterElemDOFs( LO elemIdx, LO* DOFs ) const {

  if ( testFunctionType == p0 ) {
    DOFs[0] = elemIdx;
  } else if ( testFunctionType == p1 ) {
    LO idx = outerElems[elemIdx];
    leftMesh->getElement( idx, DOFs );
  } else if ( testFunctionType == p1dis ) {
    LO idx = outerElems[elemIdx];
    DOFs[0] = 3 * idx;
    DOFs[1] = 3 * idx + 1;
    DOFs[2] = 3 * idx + 2;
  }

  //memcpy( DOFs, &outerIdx2DOFs[elemIdx * totalDOFsPerOuterElem], totalDOFsPerOuterElem * sizeof (LO ) );
}

template<class LO, class SC>
void BESpace<LO, SC>::getInnerElemDOFs( LO elemIdx, LO* DOFs ) const {

  if ( ansatzFunctionType == p0 ) {
    DOFs[0] = elemIdx;
  } else if ( ansatzFunctionType == p1 ) {
    LO idx = innerElems[elemIdx];
    rightMesh->getElement( idx, DOFs );
  } else if ( ansatzFunctionType == p1dis ) {
    LO idx = innerElems[elemIdx];
    DOFs[0] = 3 * idx;
    DOFs[1] = 3 * idx + 1;
    DOFs[2] = 3 * idx + 2;
  }

  //memcpy( DOFs, &innerIdx2DOFs[elemIdx * totalDOFsPerInnerElem], totalDOFsPerInnerElem * sizeof (LO ) );
}

template<class LO, class SC>
void BESpace<LO, SC>::getClusterOuterElemDOFs( BECluster<LO, SC> const *cluster, LO elemIdx, LO *DOFs ) const {
  // get local degrees of freedom associated with given cluster element
  LO idx;
  switch ( testFunctionType ) {
    case p0:
      DOFs[0] = elemIdx;
      break;
    case p1:
      //memcpy( DOFs, &outerIdx2DOFs[( *cluster->elems )[elemIdx] * totalDOFsPerOuterElem], totalDOFsPerOuterElem * sizeof (LO ) );
      idx = outerElems[( *cluster->elems )[elemIdx]];
      leftMesh->getElement( idx, DOFs );
      DOFs[0] = std::distance(
          cluster->nodes->begin( ), cluster->nodes->find( DOFs[0] ) );
      DOFs[1] = std::distance(
          cluster->nodes->begin( ), cluster->nodes->find( DOFs[1] ) );
      DOFs[2] = std::distance(
          cluster->nodes->begin( ), cluster->nodes->find( DOFs[2] ) );
      break;
    case p1dis:
      DOFs[0] = 3 * elemIdx;
      DOFs[1] = 3 * elemIdx + 1;
      DOFs[2] = 3 * elemIdx + 2;
      break;
  }
}

template<class LO, class SC>
void BESpace<LO, SC>::getClusterInnerElemDOFs( BECluster<LO, SC> const *cluster, LO elemIdx, LO *DOFs ) const {
  // get local degrees of freedom associated with given cluster element
  LO idx;
  switch ( ansatzFunctionType ) {
    case p0:
      DOFs[0] = elemIdx;
      break;
    case p1:
      //memcpy( DOFs, &innerIdx2DOFs[( *cluster->elems )[elemIdx] * totalDOFsPerInnerElem], totalDOFsPerInnerElem * sizeof (LO ) );
      idx = innerElems[( *cluster->elems )[elemIdx]];
      rightMesh->getElement( idx, DOFs );

      DOFs[0] = std::distance(
          cluster->nodes->begin( ), cluster->nodes->find( DOFs[0] ) );
      DOFs[1] = std::distance(
          cluster->nodes->begin( ), cluster->nodes->find( DOFs[1] ) );
      DOFs[2] = std::distance(
          cluster->nodes->begin( ), cluster->nodes->find( DOFs[2] ) );
      break;
    case p1dis:
      DOFs[0] = 3 * elemIdx;
      DOFs[1] = 3 * elemIdx + 1;
      DOFs[2] = 3 * elemIdx + 2;
      break;
  }
}


}

#endif
