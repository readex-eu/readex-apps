/*!
 * @file    FMMMatrix.cpp
 * @author  Michal Merta 
 * @date    August 19, 2013
 * 
 */

#ifdef FMMMATRIX_H

namespace bem4i {

template<class LO, class SC>
FMMMatrix<LO, SC>::FMMMatrix( ) {
  this->nonadmBlocksSize = 0;
}

template<class LO, class SC>
FMMMatrix<LO, SC>::FMMMatrix( const FMMMatrix& orig ) {
}

template<class LO, class SC>
FMMMatrix<LO, SC>::FMMMatrix( LO nRows, LO nCols, FMMKernel<LO, SC>* kernel ) {
  this->kernel = kernel;
  this->nRows = nRows;
  this->nCols = nCols;
  this->nonadmBlocksSize = 0;
}

template<class LO, class SC>
FMMMatrix<LO, SC>::~FMMMatrix( ) {
  for ( int i = 0; i < nonAdmissibleBlocks.size( ); i++ ) {
    if ( nonAdmissibleBlocks[i] ) {
      delete nonAdmissibleBlocks[i];
    }
  }
}

template<class LO, class SC>
void FMMMatrix<LO, SC>::apply(
    Vector<LO, SC> const &x,
    Vector<LO, SC> &y,
    bool transA,
    SC alpha,
    SC beta
    ) {
  //TODO: multiply by alpha!
  y.scale( beta );

  // apply nonadmissible blocks
  for ( int i = 0; i < nonAdmissibleBlocks.size( ); i++ ) {
    nonAdmissibleBlocks[i]->apply( x, *( nonadmissibleLeaves[i]->innerDOFs ),
        y, *( nonadmissibleLeaves[i]->outerDOFs ), transA, alpha, 1.0 );
  }

  this->kernel->setVector( &x );
  this->kernel->upward( );
  this->kernel->downward( );
  this->kernel->getLeftClusterTree( )->setPointerToRoot( );
  applyFarfieldBlocks( y );
  //this->kernel->getLeftClusterTree()->setPointerToRoot();
}

template<class LO, class SC>
void FMMMatrix<LO, SC>::applyFarfieldBlocks( Vector<LO, SC> &result ) {
  SC res = 0.0;

  Tree<BECluster<LO, SC>*, LO > *clusterTree = this->kernel->getLeftClusterTree( );

  for ( int i = 0; i < clusterTree->getNSons( ); i++ ) {
    clusterTree->setPointerToSon( i );
    applyFarfieldBlocks( result );
    clusterTree->setPointerToParent( );
  }

  if ( clusterTree->getNSons( ) == 0 ) {
    for ( int i = 0; i < clusterTree->get( )->nelems; i++ ) {
      res = this->kernel->computeApproximation( clusterTree->get( ), i );
      //std::cout << res << std::endl;
      result.add( ( *clusterTree->get( )->elems )[i], res );
    }
  }
}

}
#endif
