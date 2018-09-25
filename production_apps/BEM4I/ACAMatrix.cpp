/*!
 * @file    ACAMatrix.cpp
 * @author  Michal Merta
 * @date    August 19, 2013
 *
 */

#ifdef ACAMATRIX_H

namespace bem4i {

template<class LO, class SC>
ACAMatrix<LO, SC>::ACAMatrix( ) {
  this->nonadmBlocksSize = 0;
  this->admBlocksSize = 0;
  this->maxBlockSize = 0;
  this->p12p1dis = false;
  this->p1dis2p1 = false;
  this->p12p1disMat = nullptr;
  this->p1dis2p1Mat = nullptr;
}

template<class LO, class SC>
ACAMatrix<LO, SC>::ACAMatrix( const ACAMatrix& orig ) {
}

template<class LO, class SC>
ACAMatrix<LO, SC>::ACAMatrix( LO nRows, LO nCols ) {

  this->nRows = nRows;
  this->nCols = nCols;
  this->nonadmBlocksSize = 0;
  this->admBlocksSize = 0;
  this->maxBlockSize = 0;
  this->p12p1dis = false;
  this->p1dis2p1 = false;
  this->p12p1disMat = nullptr;
  this->p1dis2p1Mat = nullptr;
}

template<class LO, class SC>
ACAMatrix<LO, SC>::~ACAMatrix( ) {
  for ( int i = 0; i < admissibleBlocks.size( ); i++ ) {
    if ( admissibleBlocks[i].first ) {
      delete admissibleBlocks[i].first;
    }
    if ( admissibleBlocks[i].second ) {
      delete admissibleBlocks[i].second;
    }
  }
  for ( int i = 0; i < nonAdmissibleBlocks.size( ); i++ ) {
    if ( nonAdmissibleBlocks[i] ) {
      delete nonAdmissibleBlocks[i];
    }
  }
  if ( this->p12p1disMat ) {
    delete this->p12p1disMat;
  }

  if ( this->p1dis2p1Mat ) {
    delete this->p1dis2p1Mat;
  }



  for ( auto it = this->admissibleInnerDOFs.begin( );
      it != this->admissibleInnerDOFs.end( ); ++it ) {
    delete *it;
  }
  this->admissibleInnerDOFs.clear( );
  for ( auto it = this->admissibleOuterDOFs.begin( );
      it != this->admissibleOuterDOFs.end( ); ++it ) {
    delete *it;
  }
  this->admissibleOuterDOFs.clear( );
  for ( auto it = this->nonadmissibleInnerDOFs.begin( );
      it != this->nonadmissibleInnerDOFs.end( ); ++it ) {
    delete *it;
  }
  this->nonadmissibleInnerDOFs.clear( );
  for ( auto it = this->nonadmissibleOuterDOFs.begin( );
      it != this->nonadmissibleOuterDOFs.end( ); ++it ) {
    delete *it;
  }
  this->nonadmissibleOuterDOFs.clear( );
}
/*
template<class LO, class SC>
void ACAMatrix<LO, SC>::apply(
    Vector<LO, SC> const & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  y.scale( beta );
  Vector<LO, SC> * auxY = &y;
  Vector<LO, SC> * auxX = nullptr;

  if ( ( !transA && this->p12p1dis ) || ( transA && this->p1dis2p1 ) ) {
    auxX = new Vector<LO, SC>( this->p12p1disMat->getNRows( ) );
    this->p12p1disMat->apply( x, *auxX );
  } else {
    // TODO performance issue for copying
    auxX = new Vector<LO, SC>( x );
  }

  if ( ( !transA && this->p1dis2p1 ) || ( transA && this->p12p1dis ) ) {
    auxY = new Vector<LO, SC>( this->p12p1disMat->getNRows( ) );
  } else {
    //auxY = new Vector<LO, SC>( y.getLength( ) );
  }

  std::vector<LO>* innerDOFs;
  std::vector<LO>* outerDOFs;


  if ( !transA ) {
    // apply nonadmissible blocks
    for ( int i = 0; i < nonAdmissibleBlocks.size( ); i++ ) {
      if ( nonAdmissibleBlocks[i] != nullptr ) {
        nonAdmissibleBlocks[i]->apply(
 *auxX, *nonadmissibleInnerDOFs[i], *auxY,
 *( nonadmissibleOuterDOFs[i] ), transA, alpha, 1.0 );
      }
    }

    // apply admissible blocks
    for ( int i = 0; i < admissibleBlocks.size( ); i++ ) {
      if ( admissibleBlocks[i].first == nullptr ) continue;
      innerDOFs = admissibleInnerDOFs[i];
      outerDOFs = admissibleOuterDOFs[i];
      Vector<LO, SC> localX( innerDOFs->size( ) );
      Vector<LO, SC> localY( outerDOFs->size( ), true );
      for ( LO j = 0; j < innerDOFs->size( ); j++ ) {
        localX.set( j, auxX->get( ( *innerDOFs )[j] ) );
      }

      if ( admissibleBlocks[i].second ) {
        Vector<LO, SC> tmp( admissibleBlocks[i].second->getNCols( ), true );
        admissibleBlocks[i].second->apply( localX, tmp, true );
        admissibleBlocks[i].first->apply( tmp, localY, false, alpha );
      } else {
        admissibleBlocks[i].first->apply( localX, localY, false, alpha );
      }
      for ( LO j = 0; j < outerDOFs->size( ); j++ ) {
        auxY->add( ( *outerDOFs )[j], localY.get( j ) );
      }
    }
  } else {
    // apply nonadmissible blocks
    for ( int i = 0; i < nonAdmissibleBlocks.size( ); i++ ) {
      if ( nonAdmissibleBlocks[i] != nullptr ) {
        nonAdmissibleBlocks[i]->apply(
 *auxX, *( nonadmissibleOuterDOFs[i] ), *auxY,
 *( nonadmissibleInnerDOFs[i] ), transA, alpha, 1.0 );
      }
    }

    // apply admissible blocks
    for ( int i = 0; i < admissibleBlocks.size( ); i++ ) {
      if ( admissibleBlocks[i].first == nullptr ) continue;
      innerDOFs = admissibleOuterDOFs[i];
      outerDOFs = admissibleInnerDOFs[i];
      Vector<LO, SC> localX( innerDOFs->size( ) );
      Vector<LO, SC> localY( outerDOFs->size( ), true );
      for ( LO j = 0; j < innerDOFs->size( ); j++ ) {
        localX.set( j, auxX->get( ( *innerDOFs )[j] ) );
      }

      if ( admissibleBlocks[i].second ) {
        Vector<LO, SC> tmp( admissibleBlocks[i].first->getNCols( ), true );
        admissibleBlocks[i].first->apply( localX, tmp, true );
        admissibleBlocks[i].second->apply( tmp, localY, false, alpha );
      } else {
        admissibleBlocks[i].first->apply( localX, localY, true, alpha );
      }
      for ( LO j = 0; j < outerDOFs->size( ); j++ ) {
        auxY->add( ( *outerDOFs )[j], localY.get( j ) );
      }
    }
  }

  if ( ( !transA && this->p1dis2p1 ) || ( transA && this->p12p1dis ) ) {
    this->p12p1disMat->apply( *auxY, y, true );
    delete auxY;
  }
  delete auxX;
}
 */
///*

template<class LO, class SC>
void ACAMatrix<LO, SC>::apply(
    Vector<LO, SC> const & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  y.scale( beta );
  Vector<LO, SC> * auxY = &y;
  Vector<LO, SC> * auxX = nullptr;

  if ( !transA && this->p12p1dis ) {
    auxX = new Vector<LO, SC>( this->p12p1disMat->getNRows( ) );
    this->p12p1disMat->apply( x, *auxX );
  } else if ( transA && this->p1dis2p1 ) {
    auxX = new Vector<LO, SC>( this->p1dis2p1Mat->getNCols( ) );
    this->p1dis2p1Mat->apply( x, *auxX, true );
  } else {
    // TODO performance issue for copying
    auxX = new Vector<LO, SC>( x );
  }
  if ( !transA && this->p1dis2p1 ) {
    auxY = new Vector<LO, SC>( this->p1dis2p1Mat->getNCols( ) );
  } else if ( transA && this->p12p1dis ) {
    auxY = new Vector<LO, SC>( this->p12p1disMat->getNRows( ) );
  }

#pragma omp parallel
  {
    std::vector< LO > * innerDOFs;
    std::vector< LO > * outerDOFs;
    SC * localXData = new SC[ this->maxBlockSize ];
    SC * localYData = new SC[ this->maxBlockSize ];
    SC * tmpData = new SC[ this->maxBlockSize ];
    Vector< LO, SC > localX;
    Vector< LO, SC > localY( this->maxBlockSize, localYData );
    Vector< LO, SC > tmp( this->maxBlockSize, tmpData );
    localY.setAll( 0.0 ); // has to be done for the result of multiplication!
    tmp.setAll( 0.0 ); // has to be done for the result of multiplication!

    if ( !transA ) {
      // apply nonadmissible blocks
#pragma omp for schedule( dynamic, 16 )
      for ( int i = 0; i < this->nonAdmissibleBlocks.size( ); ++i ) {
        if ( !this->nonAdmissibleBlocks[ i ] ) continue;

        innerDOFs = this->nonadmissibleInnerDOFs[ i ];
        outerDOFs = this->nonadmissibleOuterDOFs[ i ];
        localX.setData( innerDOFs->size( ), localXData );
        localY.setData( outerDOFs->size( ), localYData );

        for ( LO j = 0; j < innerDOFs->size( ); ++j ) {
          localX.set( j, auxX->get( ( *innerDOFs )[ j ] ) );
        }

        this->nonAdmissibleBlocks[ i ]->apply( localX, localY, false, alpha );

        for ( LO j = 0; j < outerDOFs->size( ); ++j ) {
          auxY->addAtomic( ( *outerDOFs )[ j ], localY.get( j ) );
        }
      }

      // apply admissible blocks
#pragma omp for schedule( dynamic, 16 )
      for ( int i = 0; i < this->admissibleBlocks.size( ); ++i ) {
        if ( !this->admissibleBlocks[ i ].first ) continue;

        innerDOFs = this->admissibleInnerDOFs[ i ];
        outerDOFs = this->admissibleOuterDOFs[ i ];
        localX.setData( innerDOFs->size( ), localXData );
        localY.setData( outerDOFs->size( ), localYData );

        for ( LO j = 0; j < innerDOFs->size( ); ++j ) {
          localX.set( j, auxX->get( ( *innerDOFs )[ j ] ) );
        }

        if ( this->admissibleBlocks[ i ].second ) {
          tmp.setData( this->admissibleBlocks[ i ].second->getNCols( ),
              tmpData );
          this->admissibleBlocks[ i ].second->apply( localX, tmp, true );
          this->admissibleBlocks[ i ].first->apply( tmp, localY, false, alpha );
        } else {
          this->admissibleBlocks[ i ].first->apply( localX, localY, false,
              alpha );
        }

        for ( LO j = 0; j < outerDOFs->size( ); ++j ) {
          auxY->addAtomic( ( *outerDOFs )[ j ], localY.get( j ) );
        }
      }
    } else {

      // apply nonadmissible blocks
#pragma omp for schedule( dynamic, 16 )
      for ( int i = 0; i < this->nonAdmissibleBlocks.size( ); ++i ) {
        if ( !this->nonAdmissibleBlocks[ i ] ) continue;

        innerDOFs = this->nonadmissibleOuterDOFs[ i ];
        outerDOFs = this->nonadmissibleInnerDOFs[ i ];
        localX.setData( innerDOFs->size( ), localXData );
        localY.setData( outerDOFs->size( ), localYData );

        for ( LO j = 0; j < innerDOFs->size( ); ++j ) {
          localX.set( j, auxX->get( ( *innerDOFs )[ j ] ) );
        }

        this->nonAdmissibleBlocks[ i ]->apply( localX, localY, true, alpha );

        for ( LO j = 0; j < outerDOFs->size( ); ++j ) {
          auxY->addAtomic( ( *outerDOFs )[ j ], localY.get( j ) );
        }
      }

      // apply admissible blocks
#pragma omp for schedule( dynamic, 16 )
      for ( int i = 0; i < this->admissibleBlocks.size( ); ++i ) {
        if ( !this->admissibleBlocks[ i ].first ) continue;

        innerDOFs = this->admissibleOuterDOFs[ i ];
        outerDOFs = this->admissibleInnerDOFs[ i ];
        localX.setData( innerDOFs->size( ), localXData );
        localY.setData( outerDOFs->size( ), localYData );

        for ( LO j = 0; j < innerDOFs->size( ); ++j ) {
          localX.set( j, auxX->get( ( *innerDOFs )[ j ] ) );
        }

        if ( this->admissibleBlocks[ i ].second ) {
          tmp.setData( this->admissibleBlocks[ i ].first->getNCols( ),
              tmpData );
          this->admissibleBlocks[ i ].first->apply( localX, tmp, true );
          this->admissibleBlocks[ i ].second->apply( tmp, localY, false,
              alpha );
        } else {
          this->admissibleBlocks[ i ].first->apply( localX, localY, true,
              alpha );
        }

        for ( LO j = 0; j < outerDOFs->size( ); ++j ) {
          auxY->addAtomic( ( *outerDOFs )[ j ], localY.get( j ) );
        }
      }
    }

    delete [] localXData;
    delete [] localYData;
    delete [] tmpData;
  }

  if ( !transA && this->p1dis2p1 ) {
    this->p1dis2p1Mat->apply( *auxY, y, false, 1.0, 1.0 );
    delete auxY;
  } else if ( transA && this->p12p1dis ) {
    this->p12p1disMat->apply( *auxY, y, true, 1.0, 1.0 );
    delete auxY;
  }

  delete auxX;
}
//*/

template<class LO, class SC>
void ACAMatrix<LO, SC>::setP12p1disMatFromTriplets(
    LO nRows,
    LO nCols,
    std::vector<LO> & vecI,
    std::vector<LO> & vecJ,
    std::vector<SC> & vecV
    ) {
  if ( p12p1disMat != nullptr ) {
    delete p12p1disMat;
  }
  this->p12p1disMat = new SparseMatrix<LO, SC>( );
  this->p12p1disMat->setFromTriplets( nRows, nCols, vecI, vecJ, vecV );
}

template<class LO, class SC>
void ACAMatrix<LO, SC>::setP1dis2p1MatFromTriplets(
    LO nRows,
    LO nCols,
    std::vector<LO> & vecI,
    std::vector<LO> & vecJ,
    std::vector<SC> & vecV
    ) {
  if ( p1dis2p1Mat != nullptr ) {
    delete p1dis2p1Mat;
  }
  this->p1dis2p1Mat = new SparseMatrix<LO, SC>( );
  this->p1dis2p1Mat->setFromTriplets( nRows, nCols, vecI, vecJ, vecV );
}




}
#endif
