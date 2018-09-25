/*!
 * @file    MPIACAMatrix.cpp
 * @author  Michal Merta 
 * @date    March 13, 2015
 * 
 */

#ifdef MPIACAMATRIX_H

namespace bem4i {

template<class LO, class SC>
MPIACAMatrix<LO, SC>::MPIACAMatrix( ) {
  this->nonadmBlocksSize = 0;
  this->admBlocksSize = 0;
  this->communicator = MPI_COMM_WORLD;
}

template<class LO, class SC>
MPIACAMatrix<LO, SC>::MPIACAMatrix(
    const MPIACAMatrix& orig
    ) {
}

template<class LO, class SC>
MPIACAMatrix<LO, SC>::MPIACAMatrix(
    LO nRows,
    LO nCols,
    MPI_Comm communicator
    ) {

  this->nRows = nRows;
  this->nCols = nCols;
  this->nonadmBlocksSize = 0;
  this->admBlocksSize = 0;
  this->communicator = communicator;
}

template<class LO, class SC>
MPIACAMatrix<LO, SC>::~MPIACAMatrix( ) {

}

template<class LO, class SC>
void MPIACAMatrix<LO, SC>::apply(
    Vector<LO, SC> const &x,
    Vector<LO, SC> &y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  MPI_Barrier( communicator );

  //y.scale( beta );
  std::vector<LO>* innerDOFs;
  std::vector<LO>* outerDOFs;

  LO yLength = transA ? this->nCols : this->nRows;

  Vector<LO, SC> *myY = new Vector<LO, SC>( yLength );
  Vector<LO, SC> *globalY = new Vector<LO, SC>( yLength );

  myY->setAll( 0.0 );

  if ( !transA ) {
    // apply nonadmissible blocks
    for ( int i = 0; i < this->nonAdmissibleBlocks.size( ); i++ ) {
      if ( this->nonAdmissibleBlocks[i] ) {
        this->nonAdmissibleBlocks[i]->apply( x, *( this->nonadmissibleInnerDOFs[i] ), *myY,
            *( this->nonadmissibleOuterDOFs[i] ), transA, alpha, 1.0 );
      }
    }

    // apply admissible blocks
    for ( int i = 0; i < this->admissibleBlocks.size( ); i++ ) {
      if ( !this->admissibleBlocks[i].first ) continue;
      innerDOFs = this->admissibleInnerDOFs[i];
      outerDOFs = this->admissibleOuterDOFs[i];
      Vector<LO, SC> localX( innerDOFs->size( ) );
      Vector<LO, SC> localY( outerDOFs->size( ), true );
      for ( LO j = 0; j < innerDOFs->size( ); j++ ) {
        localX.set( j, x.get( ( *innerDOFs )[j] ) );
      }

      if ( this->admissibleBlocks[i].second ) {
        Vector<LO, SC> tmp( this->admissibleBlocks[i].second->getNCols( ), true );
        this->admissibleBlocks[i].second->apply( localX, tmp, true );
        this->admissibleBlocks[i].first->apply( tmp, localY, false, alpha );
      } else {
        this->admissibleBlocks[i].first->apply( localX, localY, false, alpha );
      }
      for ( LO j = 0; j < outerDOFs->size( ); j++ ) {
        myY->add( ( *outerDOFs )[j], localY.get( j ) );
      }
    }
  } else {
    // apply nonadmissible blocks
    for ( int i = 0; i < this->nonAdmissibleBlocks.size( ); i++ ) {
      if ( this->nonAdmissibleBlocks[i] != nullptr ) {
        this->nonAdmissibleBlocks[i]->apply( x, *( this->nonadmissibleOuterDOFs[i] ), *myY,
            *( this->nonadmissibleInnerDOFs[i] ), transA, alpha, 1.0 );
      }
    }

    // apply admissible blocks
    for ( int i = 0; i < this->admissibleBlocks.size( ); i++ ) {
      if ( this->admissibleBlocks[i].first == nullptr ) continue;
      innerDOFs = this->admissibleOuterDOFs[i];
      outerDOFs = this->admissibleInnerDOFs[i];
      Vector<LO, SC> localX( innerDOFs->size( ) );
      Vector<LO, SC> localY( outerDOFs->size( ), true );
      for ( LO j = 0; j < innerDOFs->size( ); j++ ) {
        localX.set( j, x.get( ( *innerDOFs )[j] ) );
      }

      if ( this->admissibleBlocks[i].second ) {
        Vector<LO, SC> tmp( this->admissibleBlocks[i].first->getNCols( ), true );
        this->admissibleBlocks[i].first->apply( localX, tmp, true );
        this->admissibleBlocks[i].second->apply( tmp, localY, false, alpha );
      } else {
        this->admissibleBlocks[i].first->apply( localX, localY, true, alpha );
      }
      for ( LO j = 0; j < outerDOFs->size( ); j++ ) {
        myY->add( ( *outerDOFs )[j], localY.get( j ) );
      }
    }
  }

  MPI_Barrier( communicator );
  MPI_Allreduce( myY->getData( ), globalY->getData( ), yLength,
      GetType<LO, SC>::MPI_SC( ), MPI_SUM, communicator );

  y.scale( beta );
  y.add( *globalY );

  delete myY;
  delete globalY;

}





}
#endif

