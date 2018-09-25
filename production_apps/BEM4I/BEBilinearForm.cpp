/*!
 * @file    BEBilinearForm.cpp
 * @author  Michal Merta 
 * @date    July 12, 2013
 * 
 */



#ifdef BEBILINEARFORM_H

namespace bem4i {

template<class LO, class SC>
BEBilinearForm<LO, SC>::BEBilinearForm( ) {
  
  this->counterRow = 0;
  this->counterCol = 0;
}

template<class LO, class SC>
BEBilinearForm<LO, SC>::BEBilinearForm(
    const BEBilinearForm& orig
    ) {
}

template<class LO, class SC>
BEBilinearForm<LO, SC>::~BEBilinearForm( ) {
}

//template<class LO, class SC>
//void BEBilinearForm<LO, SC>::assembleACABlock(
//    BEBlockCluster<LO, SC>* block,
//    FullMatrix<LO, SC> & U,
//    FullMatrix<LO, SC> & V,
//    void * voidIntegrator
//    ) const {
//
//  FastBESpace<LO, SC>* fastSpace =
//      static_cast<FastBESpace<LO, SC>*> ( this->space );
//
//  SCVT zeroEps = 1e-14 * fastSpace->getScaleACA( );
//
//  BECluster<LO, SC> * leftCluster = block->leftCluster;
//  BECluster<LO, SC> * rightCluster = block->rightCluster;
//
//  LO nRows = fastSpace->getClusterOuterDOFs( leftCluster );
//  LO nCols = fastSpace->getClusterInnerDOFs( rightCluster );
//  std::vector<LO> idxI;
//  idxI.reserve( nRows );
//
//  std::vector<LO> idxJ;
//  idxJ.reserve( nCols );
//  std::vector< Vector< LO, SC > * > u;
//  u.reserve( nRows );
//  std::vector< Vector< LO, SC > * > v;
//  v.reserve( nCols );
//
//  std::vector<bool> idxICompl( nRows, true );
//  std::vector<bool> idxJCompl( nCols, true );
//
//  Vector<LO, SCVT> c( nRows, true );
//  Vector<LO, SCVT> r( nCols, true );
//
//  LO iters = 0;
//  LO rank = 0;
//
//  bool typeIsRow = true;
//  bool searchForRow = true;
//
//  LO pivotRow = 0;
//  LO pivotCol = 0;
//  SC gamma;
//  SCVT squareNormS = 0.0;
//  SCVT normUV = 0.0;
//  SC maxRowRes = 0.0;
//  SC maxColRes = 0.0;
//
//  // preallocate memory for rows and columns
//  SCVT percentPrealloc = 0.1;
//  LO preallocSize = (LO) ( percentPrealloc * (SCVT) std::min( nRows, nCols ) );
//  SC * rowBuffer = new SC[ preallocSize * nCols ];
//  SC * colBuffer = new SC[ preallocSize * nRows ];
//  SC * currRowData = nullptr;
//  SC * currColData = nullptr;
//  Vector< LO, SC > * currentRow = nullptr;
//  Vector< LO, SC > * currentColumn = nullptr;
//  bool deleteData = false;
//
//  // stop if all rows and cols used
//  while ( idxI.size( ) != nRows && idxJ.size( ) != nCols ) {
//
//    // find pivot row
//    if ( searchForRow ) {
//      pivotRow = std::find( idxICompl.begin( ), idxICompl.end( ), true )
//          - idxICompl.begin( );
//      typeIsRow = true;
//    }
//    // no possible row found (should not happen, loop should end sooner)
//    if ( pivotRow == nRows ) {
//      break;
//    }
//    if ( typeIsRow ) {
//      // if buffer is not full, use memory from it
//      if ( u.size( ) < preallocSize ) {
//        currRowData = rowBuffer + nCols * u.size( );
//        deleteData = false;
//      } else {
//        currRowData = new SC[ nCols ];
//        deleteData = true;
//      }
//      //currentRow->setData( nCols, currRowData, false );
//      currentRow = new Vector< LO, SC >( nCols, currRowData, false );
//      currentRow->setDeleteData( deleteData );
//
//      this->assembleRow( *block, ( *block->outerDOFs )[pivotRow],
//          *block->innerDOFs, *currentRow, voidIntegrator );
//
//      idxI.push_back( pivotRow );
//      idxICompl[ pivotRow ] = false;
//
//      if ( currentRow->norm2( ) < zeroEps ) {
//        delete currentRow;
//        currentRow = nullptr;
//        searchForRow = true;
//        ++iters;
//        continue;
//      }
//
//      for ( LO i = 0; i < nCols; ++i ) {
//        r.add( i, std::abs( currentRow->get( i ) ) );
//      }
//
//      for ( LO i = 0; i < u.size( ); ++i ) {
//        currentRow->add( *v[ i ], -u[ i ]->get( pivotRow ) );
//      }
//      pivotCol = currentRow->findAbsMax( ); // <-- vyjmout sloupec, ktery uz byl pouzit
//      maxRowRes = currentRow->get( pivotCol );
//
//      if ( std::abs( maxRowRes ) > zeroEps ) {
//        //if buffer is not full, use memory from it
//        if ( u.size( ) < preallocSize ) {
//          currColData = colBuffer + nRows * u.size( );
//          deleteData = false;
//        } else {
//          currColData = new SC[ nRows ];
//          deleteData = true;
//        }
//        //currentColumn->setData( nRows, currColData, false );
//        currentColumn = new Vector< LO, SC >( nRows, currColData, false );
//        currentColumn->setDeleteData( deleteData );
//        gamma = ( (SCVT) 1.0 ) / maxRowRes;
//
//        this->assembleColumn( *block, ( *block->innerDOFs )[pivotCol],
//            *block->outerDOFs, *currentColumn, voidIntegrator );
//
//        idxJ.push_back( pivotCol );
//        idxJCompl[ pivotCol ] = false;
//
//        for ( LO i = 0; i < nRows; ++i ) {
//          c.add( i, std::abs( currentColumn->get( i ) ) );
//        }
//
//        for ( LO i = 0; i < v.size( ); ++i ) {
//          currentColumn->add( *u[ i ], -v[ i ]->get( pivotCol ) );
//        }
//
//        pivotRow = currentColumn->findAbsMax( );
//        u.push_back( currentColumn );
//        currentRow->scale( gamma );
//        v.push_back( currentRow );
//        ++rank;
//      } else {
//        delete currentRow;
//        currentRow = nullptr;
//        ++iters;
//        searchForRow = true;
//        continue;
//      }
//    } else { // looking for column
//
//      // if buffer is not full, use memory from it
//      if ( u.size( ) < preallocSize ) {
//        currColData = colBuffer + nRows * u.size( );
//        deleteData = false;
//      } else {
//        currColData = new SC[ nRows ];
//        deleteData = true;
//      }
//      //currentColumn->setData( nRows, currColData, false );
//      currentColumn = new Vector< LO, SC >( nRows, currColData, false );
//      currentColumn->setDeleteData( deleteData );
//
//      this->assembleColumn( *block, ( *block->innerDOFs )[pivotCol],
//          *block->outerDOFs, *currentColumn, voidIntegrator );
//
//      idxJ.push_back( pivotCol );
//      idxJCompl[ pivotCol ] = false;
//
//      if ( currentColumn->norm2( ) < zeroEps ) {
//        delete currentColumn;
//        currentColumn = nullptr;
//        searchForRow = true;
//        ++iters;
//        continue;
//      }
//
//      for ( LO i = 0; i < nRows; ++i ) {
//        c.add( i, std::abs( currentColumn->get( i ) ) );
//      }
//
//      for ( LO i = 0; i < v.size( ); ++i ) {
//        currentColumn->add( *u[ i ], -v[ i ]->get( pivotCol ) );
//      }
//
//      pivotRow = currentColumn->findAbsMax( );
//      maxColRes = currentColumn->get( pivotRow );
//
//      if ( std::abs( maxColRes ) > zeroEps ) {
//
//        if ( u.size( ) < preallocSize ) {
//          currRowData = rowBuffer + nCols * u.size( );
//          deleteData = false;
//        } else {
//          currRowData = new SC[ nCols ];
//          deleteData = true;
//        }
//        //currentRow->setData( nCols, currRowData, false );
//        currentRow = new Vector< LO, SC >( nCols, currRowData, false );
//        currentRow->setDeleteData( deleteData );
//
//        gamma = ( (SCVT) 1.0 ) / maxColRes;
//        this->assembleRow( *block, ( *block->outerDOFs )[ pivotRow ],
//            *block->innerDOFs, *currentRow, voidIntegrator );
//
//        idxI.push_back( pivotRow );
//        idxICompl[ pivotRow ] = false;
//
//        for ( LO i = 0; i < nCols; ++i ) {
//          r.add( i, std::abs( currentRow->get( i ) ) );
//        }
//
//        for ( LO i = 0; i < u.size( ); ++i ) {
//          currentRow->add( *v[ i ], -u[ i ]->get( pivotRow ) );
//        }
//
//        pivotCol = currentRow->findAbsMax( );
//        currentColumn->scale( gamma );
//        u.push_back( currentColumn );
//        v.push_back( currentRow );
//        ++rank;
//      } else {
//        delete currentColumn;
//        currentColumn = nullptr;
//        searchForRow = true;
//        ++iters;
//        continue;
//      }
//    } // end column
//
//    // check approximation
//    //if ( maxRowRes > zeroEps && maxColRes > zeroEps ) {
//    if ( currentColumn && currentRow ) {
//      normUV = currentColumn->norm2( ) * currentRow->norm2( );
//      squareNormS += normUV * normUV;
//      for ( LO i = 0; i < u.size( ) - 1; ++i ) {
//        squareNormS += 2.0 * std::real( u.back( )->dot( *u[ i ] ) *
//            v.back( )->dot( *v[ i ] ) );
//      }
//    }
//    // stopping criterion
//    // this stopping criterion modified to sqrt( eps )!
//    if ( currentColumn && currentRow &&
//        normUV > fastSpace->getEpsilonACA( ) * std::sqrt( squareNormS ) ) {
//
//      searchForRow = false;
//      ++iters;
//      continue;
//    } else {
//      searchForRow = true;
//
//      for ( LO i = 0; i < idxICompl.size( ); ++i ) {
//        if ( idxICompl[ i ] && std::abs( c.get( i ) ) < zeroEps ) {
//          pivotRow = i;
//          typeIsRow = true;
//          searchForRow = false;
//          break;
//        }
//      }
//      if ( !searchForRow ) {
//        ++iters;
//        continue;
//      }
//
//      for ( LO i = 0; i < idxJCompl.size( ); ++i ) {
//        if ( idxJCompl[ i ] && std::abs( r.get( i ) ) < zeroEps ) {
//          pivotCol = i;
//          typeIsRow = false;
//          searchForRow = false;
//          break;
//        }
//      }
//      if ( !searchForRow ) {
//        ++iters;
//        continue;
//      }
//
//      break;
//    }
//  } // end while
//
//  // setting up matrices U (columns) and V (rows)
//  U.resize( nRows, u.size( ) );
//  V.resize( nCols, v.size( ) );
//
//  for ( LO i = 0; i < u.size( ); ++i ) {
//    memcpy( U.getData( ) + i * nRows, u[ i ]->getData( ),
//        nRows * sizeof ( SC ) );
//    memcpy( V.getData( ) + i * nCols, v[ i ]->getData( ),
//        nCols * sizeof ( SC ) );
//    delete u[ i ];
//    delete v[ i ];
//  }
//
//  delete [] rowBuffer;
//  delete [] colBuffer;
//} // end of ACA()

template<class LO, class SC>
void BEBilinearForm<LO, SC>::assembleACABlock(
    BEBlockCluster<LO, SC>* block,
    FullMatrix<LO, SC> & U,
    FullMatrix<LO, SC> & V,
    void * voidIntegrator
    ) const {

  FastBESpace<LO, SC>* fastSpace =
      static_cast<FastBESpace<LO, SC>*> ( this->space );

  SCVT zero = (SCVT) 1e-14;
  SCVT zeroEps = zero * fastSpace->getScaleACA( );

  BECluster<LO, SC> * leftCluster = block->leftCluster;
  BECluster<LO, SC> * rightCluster = block->rightCluster;

  LO nRows = fastSpace->getClusterOuterDOFs( leftCluster );
  LO nCols = fastSpace->getClusterInnerDOFs( rightCluster );
  std::vector<LO> idxI;
  idxI.reserve( nRows );

  std::vector<LO> idxJ;
  idxJ.reserve( nCols );
  std::vector< Vector< LO, SC > * > u;
  u.reserve( nRows );
  std::vector< Vector< LO, SC > * > v;
  v.reserve( nCols );

  std::vector<bool> idxICompl( nRows, true );
  std::vector<bool> idxJCompl( nCols, true );

  Vector<LO, SCVT> c( nRows, true );
  Vector<LO, SCVT> r( nCols, true );

  LO iters = 0;
  LO rank = 0;

  bool typeIsRow = true;
  bool searchForRow = true;

  LO pivotRow = 0;
  LO pivotCol = 0;
  SC gamma;
  SCVT squareNormS = 0.0;
  SCVT normUV = 0.0;
  SC maxRowRes = 0.0;
  SC maxColRes = 0.0;

  // preallocate memory for rows and columns
  SCVT percentPrealloc = 0.1;
  LO preallocSize = (LO) ( percentPrealloc * (SCVT) std::min( nRows, nCols ) );
  SC * rowBuffer = new SC[ preallocSize * nCols ];
  SC * colBuffer = new SC[ preallocSize * nRows ];
  SC * rowCopy = new SC[ nCols ];
  SC * colCopy = new SC[ nRows ];
  SC * currRowData = nullptr;
  SC * currColData = nullptr;
  Vector< LO, SC > * currentRow = nullptr;
  Vector< LO, SC > * currentColumn = nullptr;
  bool deleteData = false;

  // stop if all rows and cols used
  while ( idxI.size( ) != nRows && idxJ.size( ) != nCols ) {

    // find pivot row
    if ( searchForRow ) {
      pivotRow = std::find( idxICompl.begin( ), idxICompl.end( ), true )
          - idxICompl.begin( );
      typeIsRow = true;
    }
    // no possible row found (should not happen, loop should end sooner)
    if ( pivotRow == nRows ) {
      break;
    }
    if ( typeIsRow ) {
      // if buffer is not full, use memory from it
      if ( u.size( ) < preallocSize ) {
        currRowData = rowBuffer + nCols * u.size( );
        deleteData = false;
      } else {
        currRowData = new SC[ nCols ];
        deleteData = true;
      }

      currentRow = new Vector< LO, SC >( nCols, currRowData, false );
      currentRow->setDeleteData( deleteData );

      this->assembleRow( *block, ( *block->outerDOFs )[pivotRow],
          *block->innerDOFs, *currentRow, voidIntegrator );

      idxI.push_back( pivotRow );
      idxICompl[ pivotRow ] = false;

      memcpy( rowCopy, currentRow->getData( ), nCols * sizeof ( SC ) );

      for ( LO i = 0; i < u.size( ); ++i ) {
        currentRow->add( *v[ i ], -u[ i ]->get( pivotRow ) );
      }
      pivotCol = currentRow->findAbsMax( ); // <-- vyjmout sloupec, ktery uz byl pouzit
      maxRowRes = currentRow->get( pivotCol );

      if ( std::abs( maxRowRes ) > zeroEps ) {

        for ( LO i = 0; i < nCols; ++i ) {
          r.add( i, std::abs( rowCopy[ i ] ) );
        }

        //if buffer is not full, use memory from it
        if ( u.size( ) < preallocSize ) {
          currColData = colBuffer + nRows * u.size( );
          deleteData = false;
        } else {
          currColData = new SC[ nRows ];
          deleteData = true;
        }

        currentColumn = new Vector< LO, SC >( nRows, currColData, false );
        currentColumn->setDeleteData( deleteData );
        gamma = ( (SCVT) 1.0 ) / maxRowRes;

        this->assembleColumn( *block, ( *block->innerDOFs )[pivotCol],
            *block->outerDOFs, *currentColumn, voidIntegrator );

        idxJ.push_back( pivotCol );
        idxJCompl[ pivotCol ] = false;

        for ( LO i = 0; i < nRows; ++i ) {
          c.add( i, std::abs( currentColumn->get( i ) ) );
        }

        for ( LO i = 0; i < v.size( ); ++i ) {
          currentColumn->add( *u[ i ], -v[ i ]->get( pivotCol ) );
        }

        pivotRow = currentColumn->findAbsMax( );
        u.push_back( currentColumn );
        currentRow->scale( gamma );
        v.push_back( currentRow );
        ++rank;
      } else {
        delete currentRow;
        currentRow = nullptr;
        ++iters;
        searchForRow = true;
        continue;
      }
    } else { // looking for column

      // if buffer is not full, use memory from it
      if ( u.size( ) < preallocSize ) {
        currColData = colBuffer + nRows * u.size( );
        deleteData = false;
      } else {
        currColData = new SC[ nRows ];
        deleteData = true;
      }

      currentColumn = new Vector< LO, SC >( nRows, currColData, false );
      currentColumn->setDeleteData( deleteData );

      this->assembleColumn( *block, ( *block->innerDOFs )[pivotCol],
          *block->outerDOFs, *currentColumn, voidIntegrator );

      idxJ.push_back( pivotCol );
      idxJCompl[ pivotCol ] = false;

      memcpy( colCopy, currentColumn->getData( ), nRows * sizeof ( SC ) );

      for ( LO i = 0; i < v.size( ); ++i ) {
        currentColumn->add( *u[ i ], -v[ i ]->get( pivotCol ) );
      }

      pivotRow = currentColumn->findAbsMax( );
      maxColRes = currentColumn->get( pivotRow );

      if ( std::abs( maxColRes ) > zeroEps ) {

        for ( LO i = 0; i < nRows; ++i ) {
          c.add( i, std::abs( colCopy[ i ] ) );
        }

        if ( u.size( ) < preallocSize ) {
          currRowData = rowBuffer + nCols * u.size( );
          deleteData = false;
        } else {
          currRowData = new SC[ nCols ];
          deleteData = true;
        }

        currentRow = new Vector< LO, SC >( nCols, currRowData, false );
        currentRow->setDeleteData( deleteData );

        gamma = ( (SCVT) 1.0 ) / maxColRes;
        this->assembleRow( *block, ( *block->outerDOFs )[ pivotRow ],
            *block->innerDOFs, *currentRow, voidIntegrator );

        idxI.push_back( pivotRow );
        idxICompl[ pivotRow ] = false;

        for ( LO i = 0; i < nCols; ++i ) {
          r.add( i, std::abs( currentRow->get( i ) ) );
        }

        for ( LO i = 0; i < u.size( ); ++i ) {
          currentRow->add( *v[ i ], -u[ i ]->get( pivotRow ) );
        }

        pivotCol = currentRow->findAbsMax( );
        currentColumn->scale( gamma );
        u.push_back( currentColumn );
        v.push_back( currentRow );
        ++rank;
      } else {
        delete currentColumn;
        currentColumn = nullptr;
        searchForRow = true;
        ++iters;
        continue;
      }
    } // end column

    // check approximation
    if ( currentColumn && currentRow ) {
      normUV = currentColumn->norm2( ) * currentRow->norm2( );
      squareNormS += normUV * normUV;
      for ( LO i = 0; i < u.size( ) - 1; ++i ) {
        squareNormS += 2.0 * std::real( u.back( )->dot( *u[ i ] ) *
            v.back( )->dot( *v[ i ] ) );
      }

      //zeroEps = zero * normUV / std::sqrt( nRows * nCols );
    }
    // stopping criterion
    if ( currentColumn && currentRow &&
        normUV > fastSpace->getEpsilonACA( ) * std::sqrt( squareNormS ) ) {

      searchForRow = false;
      ++iters;
      continue;
    } else {
      searchForRow = true;

      for ( LO i = 0; i < idxICompl.size( ); ++i ) {
        if ( idxICompl[ i ] && std::abs( c.get( i ) ) < zeroEps ) {
          pivotRow = i;
          typeIsRow = true;
          searchForRow = false;
          break;
        }
      }
      if ( !searchForRow ) {
        ++iters;
        continue;
      }

      for ( LO i = 0; i < idxJCompl.size( ); ++i ) {
        if ( idxJCompl[ i ] && std::abs( r.get( i ) ) < zeroEps ) {
          pivotCol = i;
          typeIsRow = false;
          searchForRow = false;
          break;
        }
      }
      if ( !searchForRow ) {
        ++iters;
        continue;
      }

      break;
    }
  } // end while

  // setting up matrices U (columns) and V (rows)
  U.resize( nRows, u.size( ) );
  V.resize( nCols, v.size( ) );

  for ( LO i = 0; i < u.size( ); ++i ) {
    memcpy( U.getData( ) + i * nRows, u[ i ]->getData( ),
        nRows * sizeof ( SC ) );
    memcpy( V.getData( ) + i * nCols, v[ i ]->getData( ),
        nCols * sizeof ( SC ) );
    delete u[ i ];
    delete v[ i ];
  }

  delete [] rowBuffer;
  delete [] colBuffer;
  delete [] rowCopy;
  delete [] colCopy;
} // end of ACA()

template<class LO, class SC>
void BEBilinearForm<LO, SC>::leavesOuterP1dis2p1( ) const {

  FastBESpace<LO, SC>* fastSpace =
      static_cast<FastBESpace<LO, SC> *> ( this->space );
  std::vector< BEBlockCluster<LO, SC> * > * leaves =
      fastSpace->getNonadmissibleLeaves( );

  for ( LO i = 0; i < leaves->size( ); ++i ) {
    leaves->at( i )->outerDOFs->clear( );
    leaves->at( i )->outerDOFs->reserve( leaves->at( i )->leftCluster->nnodes );
    std::copy( leaves->at( i )->leftCluster->nodes->begin( ),
        leaves->at( i )->leftCluster->nodes->end( ),
        std::back_inserter( *( leaves->at( i )->outerDOFs ) ) );
  }

  leaves = fastSpace->getAdmissibleLeaves( );

  for ( LO i = 0; i < leaves->size( ); ++i ) {
    leaves->at( i )->outerDOFs->clear( );
    leaves->at( i )->outerDOFs->reserve( leaves->at( i )->leftCluster->nnodes );
    std::copy( leaves->at( i )->leftCluster->nodes->begin( ),
        leaves->at( i )->leftCluster->nodes->end( ),
        std::back_inserter( *( leaves->at( i )->outerDOFs ) ) );
  }
}

template<class LO, class SC>
void BEBilinearForm<LO, SC>::leavesInnerP1dis2p1( ) const {

  FastBESpace<LO, SC>* fastSpace =
      static_cast<FastBESpace<LO, SC> *> ( this->space );
  std::vector< BEBlockCluster<LO, SC> * > * leaves =
      fastSpace->getNonadmissibleLeaves( );

  for ( LO i = 0; i < leaves->size( ); ++i ) {
    leaves->at( i )->innerDOFs->clear( );
    leaves->at( i )->innerDOFs->reserve(
        leaves->at( i )->rightCluster->nnodes );
    std::copy( leaves->at( i )->rightCluster->nodes->begin( ),
        leaves->at( i )->rightCluster->nodes->end( ),
        std::back_inserter( *( leaves->at( i )->innerDOFs ) ) );
  }

  leaves = fastSpace->getAdmissibleLeaves( );

  for ( LO i = 0; i < leaves->size( ); ++i ) {
    leaves->at( i )->innerDOFs->clear( );
    leaves->at( i )->innerDOFs->reserve(
        leaves->at( i )->rightCluster->nnodes );
    std::copy( leaves->at( i )->rightCluster->nodes->begin( ),
        leaves->at( i )->rightCluster->nodes->end( ),
        std::back_inserter( *( leaves->at( i )->innerDOFs ) ) );
  }
}

template<class LO, class SC>
void BEBilinearForm<LO, SC>::leavesOuterP12p1dis( ) const {

  FastBESpace<LO, SC>* fastSpace =
      static_cast<FastBESpace<LO, SC> *> ( this->space );
  std::vector< BEBlockCluster<LO, SC> * > * leaves =
      fastSpace->getNonadmissibleLeaves( );

  for ( LO i = 0; i < leaves->size( ); ++i ) {
    leaves->at( i )->outerDOFs->clear( );
    leaves->at( i )->outerDOFs->reserve(
        3 * leaves->at( i )->leftCluster->nelems );
    for ( LO j = 0; j < leaves->at( i )->leftCluster->nelems; ++j ) {
      leaves->at( i )->outerDOFs->push_back(
          3 * leaves->at( i )->leftCluster->elems->at( j ) );
      leaves->at( i )->outerDOFs->push_back(
          3 * leaves->at( i )->leftCluster->elems->at( j ) + 1 );
      leaves->at( i )->outerDOFs->push_back(
          3 * leaves->at( i )->leftCluster->elems->at( j ) + 2 );
    }
  }

  leaves = fastSpace->getAdmissibleLeaves( );

  for ( LO i = 0; i < leaves->size( ); ++i ) {
    leaves->at( i )->outerDOFs->clear( );
    leaves->at( i )->outerDOFs->reserve(
        3 * leaves->at( i )->leftCluster->nelems );
    for ( LO j = 0; j < leaves->at( i )->leftCluster->nelems; ++j ) {
      leaves->at( i )->outerDOFs->push_back(
          3 * leaves->at( i )->leftCluster->elems->at( j ) );
      leaves->at( i )->outerDOFs->push_back(
          3 * leaves->at( i )->leftCluster->elems->at( j ) + 1 );
      leaves->at( i )->outerDOFs->push_back(
          3 * leaves->at( i )->leftCluster->elems->at( j ) + 2 );
    }
  }
}

template<class LO, class SC>
void BEBilinearForm<LO, SC>::leavesInnerP12p1dis( ) const {

  FastBESpace<LO, SC>* fastSpace =
      static_cast<FastBESpace<LO, SC> *> ( this->space );
  std::vector< BEBlockCluster<LO, SC> * > * leaves =
      fastSpace->getNonadmissibleLeaves( );

  for ( LO i = 0; i < leaves->size( ); ++i ) {
    leaves->at( i )->innerDOFs->clear( );
    leaves->at( i )->innerDOFs->reserve(
        3 * leaves->at( i )->rightCluster->nelems );
    for ( LO j = 0; j < leaves->at( i )->rightCluster->nelems; ++j ) {
      leaves->at( i )->innerDOFs->push_back(
          3 * leaves->at( i )->rightCluster->elems->at( j ) );
      leaves->at( i )->innerDOFs->push_back(
          3 * leaves->at( i )->rightCluster->elems->at( j ) + 1 );
      leaves->at( i )->innerDOFs->push_back(
          3 * leaves->at( i )->rightCluster->elems->at( j ) + 2 );
    }
  }

  leaves = fastSpace->getAdmissibleLeaves( );

  for ( LO i = 0; i < leaves->size( ); ++i ) {
    leaves->at( i )->innerDOFs->clear( );
    leaves->at( i )->innerDOFs->reserve(
        3 * leaves->at( i )->rightCluster->nelems );
    for ( LO j = 0; j < leaves->at( i )->rightCluster->nelems; ++j ) {
      leaves->at( i )->innerDOFs->push_back(
          3 * leaves->at( i )->rightCluster->elems->at( j ) );
      leaves->at( i )->innerDOFs->push_back(
          3 * leaves->at( i )->rightCluster->elems->at( j ) + 1 );
      leaves->at( i )->innerDOFs->push_back(
          3 * leaves->at( i )->rightCluster->elems->at( j ) + 2 );
    }
  }
}

template<class LO, class SC>
void BEBilinearForm<LO, SC>::assemble(
    ACAMatrix<LO, SC>& matrix
    ) const {
  
  this->counterRow = 0;
  this->counterCol = 0;

  // assemble nonadmissible blocks
  FastBESpace<LO, SC>* fastSpace =
      static_cast<FastBESpace<LO, SC>*> ( this->space );

  matrix.setNRows( this->space->getOuterDOFs( ) );
  matrix.setNCols( this->space->getInnerDOFs( ) );

  bool resetInnerP1dis2p1 = false;
  bool resetOuterP1dis2p1 = false;
  if ( fastSpace->getAnsatzFunctionType( ) == p1 ) {
    fastSpace->setAnsatzFunctionType( p1dis );
    this->leavesInnerP12p1dis( );
    resetInnerP1dis2p1 = true;
  }
  if ( fastSpace->getTestFunctionType( ) == p1 ) {
    fastSpace->setTestFunctionType( p1dis );
    this->leavesOuterP12p1dis( );
    resetOuterP1dis2p1 = true;
  }

  // iterate through nonadmissible leaves and assemble appropriate part of 
  // matrix
  std::vector<BEBlockCluster<LO, SC>* >* nonAdmLeaves =
      fastSpace->getNonadmissibleLeaves( );

  matrix.setNonadmissibleDOFs( *nonAdmLeaves );

  int nNonAdm = nonAdmLeaves->size( );
  matrix.resizeNonAdmBlocks( nNonAdm );

  std::vector<BEBlockCluster<LO, SC>* >* admLeaves =
      fastSpace->getAdmissibleLeaves( );

  matrix.setAdmissibleDOFs( *admLeaves );

  int nAdm = admLeaves->size( );
  matrix.resizeAdmBlocks( nAdm );

#ifdef VERBOSE
  ProgressMonitor::init( "Assembling nonadmissible blocks", nNonAdm );
#endif

#pragma omp parallel 
  {

    BEBlockCluster<LO, SC>* block;
    BECluster<LO, SC>* leftCluster;
    BECluster<LO, SC>* rightCluster;

    void * voidIntegrator = this->createIntegrator( );

#pragma omp for schedule(dynamic)
    for ( int i = 0; i < nNonAdm; i++ ) {
      block = ( *nonAdmLeaves )[i];
      leftCluster = block->leftCluster;
      rightCluster = block->rightCluster;
      FullMatrix<LO, SC> *fullBlock = new FullMatrix<LO, SC>( );
      assemble( leftCluster, rightCluster, *fullBlock, voidIntegrator );

      matrix.addNonadmissibleBlock( fullBlock, i );
#ifdef VERBOSE
      ProgressMonitor::step( );
#endif

    }

    this->destroyIntegrator( voidIntegrator );
  }

#ifdef VERBOSE
  ProgressMonitor::init( "Assembling admissible blocks", nAdm );
#endif

#pragma omp parallel 
  {
    BEBlockCluster<LO, SC>* block;
    void * voidIntegrator = this->createIntegrator( );

    // assemble admissible blocks
#pragma omp for schedule(dynamic)
    for ( LO i = 0; i < nAdm; i++ ) {
      block = ( *admLeaves )[i];

      FullMatrix<LO, SC> *acaU = new FullMatrix<LO, SC>( );
      FullMatrix<LO, SC> *acaV = new FullMatrix<LO, SC>( );
      this->assembleACABlock( block, *acaU, *acaV, voidIntegrator );

      if ( acaU->getNCols( ) > 0 ) {
        if ( acaU->getNRows( ) * acaV->getNRows( ) <
            acaU->getNRows( ) * acaU->getNCols( ) +
            acaV->getNRows( ) * acaV->getNCols( ) ) {

          // if the size of UV' < (size of U) + (size of V) multiply matrices

          FullMatrix<LO, SC> *UV = new FullMatrix<LO, SC>(
              acaU->getNRows( ), acaV->getNRows( ) );
          UV->multiply( *acaU, *acaV, false, true );

          matrix.addAdmissibleBlock( std::pair<
              FullMatrix<LO, SC>*, FullMatrix<LO, SC>* >( UV, nullptr ), i );

          delete acaU;
          delete acaV;
        } else {
          matrix.addAdmissibleBlock( std::pair<
              FullMatrix<LO, SC>*, FullMatrix<LO, SC>* >( acaU, acaV ), i );
        }
      } else {
        matrix.addAdmissibleBlock( std::pair<
            FullMatrix<LO, SC>*, FullMatrix<LO, SC>* >( nullptr, nullptr ), i );
      }

#ifdef VERBOSE
      ProgressMonitor::step( );
#endif
    }

    this->destroyIntegrator( voidIntegrator );
  }

  if ( resetInnerP1dis2p1 ) {
    this->assembleP12p1disMat( matrix );
    matrix.setP12p1dis( true );
    fastSpace->setAnsatzFunctionType( p1 );
    this->leavesInnerP1dis2p1( );
  }
  if ( resetOuterP1dis2p1 ) {
    this->assembleP1dis2p1Mat( matrix );
    matrix.setP1dis2p1( true );
    fastSpace->setTestFunctionType( p1 );
    this->leavesOuterP1dis2p1( );
  }

}

template<class LO, class SC>
void BEBilinearForm<LO, SC>::assembleP12p1disMat(
    ACAMatrix<LO, SC> & matrix
    ) const {
  std::vector<LO> veci;
  std::vector<LO> vecj;
  std::vector<SC> vecv;

  SurfaceMesh3D<LO, SC> * mesh =
      static_cast<FastBESpace<LO, SC>*> ( this->space )->getRightMesh( );

  LO nElems = mesh->getNElements( );

  veci.reserve( 3 * nElems );
  vecj.reserve( 3 * nElems );
  vecv.reserve( 3 * nElems );

  for ( LO i = 0; i < nElems; ++i ) {
    veci.push_back( 3 * i );
    veci.push_back( 3 * i + 1 );
    veci.push_back( 3 * i + 2 );
    LO elemental[3];
    mesh->getElement( i, elemental );
    vecj.push_back( elemental[0] );
    vecj.push_back( elemental[1] );
    vecj.push_back( elemental[2] );
    vecv.push_back( 1 );
    vecv.push_back( 1 );
    vecv.push_back( 1 );
  }

  matrix.setP12p1disMatFromTriplets( 3 * nElems, mesh->getNNodes( ), veci, vecj,
      vecv );
}

template<class LO, class SC>
void BEBilinearForm<LO, SC>::assembleP1dis2p1Mat(
    ACAMatrix<LO, SC> & matrix
    ) const {
  std::vector<LO> veci;
  std::vector<LO> vecj;
  std::vector<SC> vecv;

  SurfaceMesh3D<LO, SC> * mesh =
      static_cast<FastBESpace<LO, SC>*> ( this->space )->getLeftMesh( );

  LO nElems = mesh->getNElements( );

  veci.reserve( 3 * nElems );
  vecj.reserve( 3 * nElems );
  vecv.reserve( 3 * nElems );

  for ( LO i = 0; i < nElems; ++i ) {
    veci.push_back( 3 * i );
    veci.push_back( 3 * i + 1 );
    veci.push_back( 3 * i + 2 );
    LO elemental[3];
    mesh->getElement( i, elemental );
    vecj.push_back( elemental[0] );
    vecj.push_back( elemental[1] );
    vecj.push_back( elemental[2] );
    vecv.push_back( 1 );
    vecv.push_back( 1 );
    vecv.push_back( 1 );
  }

  matrix.setP1dis2p1MatFromTriplets( mesh->getNNodes( ), 3 * nElems, vecj, veci,
      vecv );
}

template<class LO, class SC>
void BEBilinearForm<LO, SC>::assemble(
    MPIACAMatrix<LO, SC>& matrix
    ) const {

  // assemble nonadmissible blocks

  matrix.setNRows( this->space->getOuterDOFs( ) );
  matrix.setNCols( this->space->getInnerDOFs( ) );

  FastBESpace<LO, SC>* fastSpace =
      static_cast<FastBESpace<LO, SC>*> ( this->space );

  // iterate through nonadmissible leaves and assemble appropriate part of 
  // matrix
  std::vector<BEBlockCluster<LO, SC>* >* nonAdmLeaves =
      fastSpace->getNonadmissibleLeaves( );

  matrix.setNonadmissibleDOFs( *nonAdmLeaves );

  int nNonAdm = nonAdmLeaves->size( );
  matrix.resizeNonAdmBlocks( nNonAdm );

  std::vector<BEBlockCluster<LO, SC>* >* admLeaves =
      fastSpace->getAdmissibleLeaves( );

  matrix.setAdmissibleDOFs( *admLeaves );

  int nAdm = admLeaves->size( );
  matrix.resizeAdmBlocks( nAdm );

  std::cout << "Number of nonadmissible blocks: " << nNonAdm << std::endl;
  std::cout << "Number of admissible blocks: " <<
      fastSpace->getAdmissibleLeaves( )->size( ) << std::endl;

  int size, rank;
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // prepare the MPI loop over nonadmissible blocks
  int myNonAdmSize = nNonAdm / size;
  if ( rank < nNonAdm % size ) myNonAdmSize++;
  int *myNonAdmLeaves = new int[myNonAdmSize];
  for ( int i = 0; i < myNonAdmSize; i++ ) {
    myNonAdmLeaves[i] = i * size + rank;
  }

#ifdef VERBOSE
  ProgressMonitor::init( "Assembling nonadmissible blocks", nNonAdm );
#endif

#pragma omp parallel 
  {

    BEBlockCluster<LO, SC>* block;
    BECluster<LO, SC>* leftCluster;
    BECluster<LO, SC>* rightCluster;
    void * voidIntegrator = this->createIntegrator( );

#pragma omp for schedule(dynamic)
    for ( int i = 0; i < myNonAdmSize; i++ ) {
      block = ( *nonAdmLeaves )[myNonAdmLeaves[i]];
      leftCluster = block->leftCluster;
      rightCluster = block->rightCluster;
      FullMatrix<LO, SC> *fullBlock = new FullMatrix<LO, SC>( );
      assemble( leftCluster, rightCluster, *fullBlock, voidIntegrator );

      matrix.addNonadmissibleBlock( fullBlock, myNonAdmLeaves[i] );
#ifdef VERBOSE
      ProgressMonitor::step( );
#endif

      //#pragma omp atomic
      //      currIt++;
    }

    this->destroyIntegrator( voidIntegrator );
    ;
  }

  std::vector<LO> sortedIdx( admLeaves->size( ) );
  for ( LO i = 0; i < admLeaves->size( ); i++ ) {
    sortedIdx[i] = i;
  }
  std::sort( sortedIdx.begin( ), sortedIdx.end( ),
      [admLeaves]( size_t i1, size_t i2 ) {
        return (*admLeaves )[i1]->leftCluster->nelems *
            ( *admLeaves )[i1]->rightCluster->nelems >
            ( *admLeaves )[i2]->leftCluster->nelems *
            ( *admLeaves )[i2]->rightCluster->nelems;
      }
  );

  // prepare the MPI loop over admissible blocks
  int myAdmSize = nAdm / size;
  if ( rank < nAdm % size ) myAdmSize++;
  int *myAdmLeaves = new int[myAdmSize];
  for ( int i = 0; i < myAdmSize; i++ ) {
    myAdmLeaves[i] = sortedIdx[i] * size + rank;
  }


#ifdef VERBOSE
  ProgressMonitor::init( "Assembling admissible blocks", nAdm );
#endif

#pragma omp parallel 
  {
    BEBlockCluster<LO, SC>* block;
    void * voidIntegrator = this->createIntegrator( );

    // assemble admissible blocks
#pragma omp for schedule(dynamic)
    for ( LO i = 0; i < myAdmSize; i++ ) {
      block = ( *admLeaves )[myAdmLeaves[i]];

      FullMatrix<LO, SC> *acaU = new FullMatrix<LO, SC>( );
      FullMatrix<LO, SC> *acaV = new FullMatrix<LO, SC>( );
      this->assembleACABlock( block, *acaU, *acaV, voidIntegrator );


      //      if ( acaU->getNRows( ) == 1 && acaU->getNCols( ) == 1 ) {
      //        
      //        // if ACA could not approximate block, assemble it as a full 
      //        
      //        FullMatrix<LO, SC> *fullBlock = new FullMatrix<LO, SC>( );
      //        assemble( block->leftCluster, block->rightCluster, *fullBlock );
      //        matrix.addAdmissibleBlock(
      //            std::pair<FullMatrix<LO, SC>*, FullMatrix<LO, SC>* >( fullBlock,
      //            nullptr ),
      //            myAdmLeaves[i] );
      //        delete acaU;
      //        delete acaV;
      //      } else 
      if ( acaU->getNRows( ) * acaV->getNRows( ) <
          acaU->getNRows( ) * acaU->getNCols( ) +
          acaV->getNRows( ) * acaV->getNCols( ) ) {
        //if (true) {
        //std::cout << acaU->getNRows( ) << " " << acaV->getNRows( ) << std::endl;
        // if the size of UV' < (size of U) + (size of V) multiply matrices

        FullMatrix<LO, SC> *UV = new FullMatrix<LO, SC>(
            acaU->getNRows( ), acaV->getNRows( ) );
        UV->multiply( *acaU, *acaV, false, true );

        matrix.addAdmissibleBlock(
            std::pair<FullMatrix<LO, SC>*, FullMatrix<LO, SC>* >( UV, nullptr ),
            myAdmLeaves[i] );

        delete acaU;
        delete acaV;
      } else {
        matrix.addAdmissibleBlock(
            std::pair<FullMatrix<LO, SC>*, FullMatrix<LO, SC>* >( acaU, acaV ),
            myAdmLeaves[i] );
      }

#ifdef VERBOSE
      ProgressMonitor::step( );
#endif
    }
  }

  delete [] myAdmLeaves;
  delete [] myNonAdmLeaves;

  // send nonadmissible blocks to master (process with rank 0)
  //  MPI_Barrier( MPI_COMM_WORLD );
  //  MPI_Status status;
  //  if ( rank == 0 ) {
  //    for ( LO i = 0; i < nNonAdm - myNonAdmSize; i++ ) {
  //      int index = 0;
  //      int nRows = 0;
  //      int nCols = 0;
  //      MPI_Recv( &index, 1, MPI_LO, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
  //          &status );
  //      MPI_Recv( &nRows, 1, MPI_LO, status.MPI_SOURCE, MPI_ANY_TAG,
  //          MPI_COMM_WORLD, &status );
  //      MPI_Recv( &nCols, 1, MPI_LO, status.MPI_SOURCE, MPI_ANY_TAG,
  //          MPI_COMM_WORLD, &status );
  //      SC *data = new SC[nRows * nCols];
  //      MPI_Recv( data, nRows*nCols, MPI_SC, status.MPI_SOURCE, MPI_ANY_TAG,
  //          MPI_COMM_WORLD, &status );
  //      FullMatrix<LO, SC> *fullBlock = new FullMatrix<LO, SC>( nRows, nCols,
  //          data );
  //      matrix.addNonadmissibleBlock( fullBlock, index );
  //    }
  //  } else {
  //    for ( LO i = 0; i < myNonAdmSize; i++ ) {
  //      int index = myNonAdmLeaves[i];
  //      FullMatrix<LO, SC>* currentMatrix = matrix.getNonAdmissibleBlock( index );
  //      int nRows = currentMatrix->getNRows( );
  //      int nCols = currentMatrix->getNCols( );
  //      MPI_Send( &index, 1, MPI_LO, 0, 0, MPI_COMM_WORLD );
  //      MPI_Send( &nRows, 1, MPI_LO, 0, 0, MPI_COMM_WORLD );
  //      MPI_Send( &nCols, 1, MPI_LO, 0, 0, MPI_COMM_WORLD );
  //      MPI_Send( currentMatrix->getData( ), nRows*nCols, MPI_SC, 0, 0,
  //          MPI_COMM_WORLD );
  //    }
  //  }
  //  
  //  MPI_Barrier(MPI_COMM_WORLD);
  //  // send admissible blocks to master
  //  if ( rank == 0 ) {
  //    for ( LO i = 0; i < nAdm - myAdmSize; i++ ) {
  //      int index = 0;
  //      int nRowsU = 0;
  //      int nColsU = 0;
  //      int nRowsV = 0;
  //      int nColsV = 0;
  //      // first receive U
  //      MPI_Recv( &index, 1, MPI_LO, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
  //          &status );
  //      MPI_Recv( &nRowsU, 1, MPI_LO, status.MPI_SOURCE, MPI_ANY_TAG,
  //          MPI_COMM_WORLD, &status );
  //      MPI_Recv( &nColsU, 1, MPI_LO, status.MPI_SOURCE, MPI_ANY_TAG,
  //          MPI_COMM_WORLD, &status );
  //      SC *dataU = new SC[nRowsU * nColsU];
  //      MPI_Recv( dataU, nRowsU*nColsU, MPI_SC, status.MPI_SOURCE, MPI_ANY_TAG,
  //          MPI_COMM_WORLD, &status );
  //      FullMatrix<LO, SC> *U = new FullMatrix<LO, SC>( nRowsU, nColsU,
  //          dataU );
  //      // then receive V (if not null)
  //      FullMatrix<LO, SC> *V = nullptr;
  //      MPI_Recv( &nRowsV, 1, MPI_LO, status.MPI_SOURCE, MPI_ANY_TAG,
  //          MPI_COMM_WORLD, &status );
  //      if ( nRowsV != 0 ) {
  //        MPI_Recv( &nColsV, 1, MPI_LO, status.MPI_SOURCE, MPI_ANY_TAG,
  //            MPI_COMM_WORLD, &status );
  //        SC *dataV = new SC[nRowsV * nColsV];
  //        MPI_Recv( dataV, nRowsV*nColsV, MPI_SC, status.MPI_SOURCE, MPI_ANY_TAG,
  //            MPI_COMM_WORLD, &status );
  //        V = new FullMatrix<LO, SC>( nRowsV, nColsV, dataV );
  //      }
  //
  //      matrix.addAdmissibleBlock(
  //          std::pair<FullMatrix<LO, SC>*, FullMatrix<LO, SC>* >( U, V ), index );
  //    }
  //  } else {
  //    for ( LO i = 0; i < myAdmSize; i++ ) {
  //      int index = myAdmLeaves[i];
  //      std::pair<FullMatrix<LO, SC>*, FullMatrix<LO, SC>* > currentPair =
  //          matrix.getAdmissibleBlock( index );
  //      int nRowsU = currentPair.first->getNRows( );
  //      int nColsU = currentPair.first->getNCols( );
  //      MPI_Send( &index, 1, MPI_LO, 0, 0, MPI_COMM_WORLD );
  //      MPI_Send( &nRowsU, 1, MPI_LO, 0, 0, MPI_COMM_WORLD );
  //      MPI_Send( &nColsU, 1, MPI_LO, 0, 0, MPI_COMM_WORLD );
  //      MPI_Send( currentPair.first->getData( ), nRowsU*nColsU, MPI_SC, 0, 0,
  //          MPI_COMM_WORLD );
  //      if ( currentPair.second == nullptr ) {
  //        int zero = 0;
  //        MPI_Send( &zero, 1, MPI_LO, 0, 0, MPI_COMM_WORLD );
  //      } else {
  //        int nRowsV = currentPair.second->getNRows( );
  //        int nColsV = currentPair.second->getNCols( );
  //        MPI_Send( &nRowsV, 1, MPI_LO, 0, 0, MPI_COMM_WORLD );
  //        MPI_Send( &nColsV, 1, MPI_LO, 0, 0, MPI_COMM_WORLD );
  //        MPI_Send( currentPair.second->getData( ), nRowsV*nColsV, MPI_SC, 0, 0,
  //            MPI_COMM_WORLD );
  //      }
  //    }
  //  }


}





}


#endif
