/*!
 * @file    BEBilinearFormHelmholtzHypersingular.cpp
 * @author  Jan Zapletal
 * @date    September 3, 2013
 * 
 */

#ifdef BEBILINEARFORMHELMHOLTZHYPERSINGULAR_H

namespace bem4i {

template<class LO, class SC>
BEBilinearFormHelmholtzHypersingular<LO, SC>::
BEBilinearFormHelmholtzHypersingular( ) {
  initCurl( );
}

template<class LO, class SC>
BEBilinearFormHelmholtzHypersingular<LO, SC>::
BEBilinearFormHelmholtzHypersingular(
    const BEBilinearFormHelmholtzHypersingular & orig
    ) {
}

template<class LO, class SC>
BEBilinearFormHelmholtzHypersingular<LO, SC>::
BEBilinearFormHelmholtzHypersingular(
    BESpace<LO, SC>* space,
    int* quadratureOrder,
    SC kappa,
    quadratureType quadrature,
    int* quadratureOrderDisjointElems
    ) {

  this->quadrature = quadrature;
  this->space = space;
  this->kappa = kappa;
  if ( quadratureOrder ) {
    this->quadratureOrder = quadratureOrder;
  } else {
    switch ( quadrature ) {
      case SauterSchwab:
        this->quadratureOrder = defaultQuadraturesSauterSchwab;
        break;
      case Steinbach:
        this->quadratureOrder = defaultQuadraturesSteinbach;
    }
  }
  initCurl( );

  this->quadratureOrderDisjointElems = quadratureOrderDisjointElems;

}

template<class LO, class SC>
BEBilinearFormHelmholtzHypersingular<LO, SC>::
~BEBilinearFormHelmholtzHypersingular( ) {
}

template <class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assemble(
    BECluster<LO, SC> const *leftCluster,
    BECluster<LO, SC> const *rightCluster,
    FullMatrix<LO, SC>& matrix,
    void * voidIntegrator
    ) const {

  if ( this->space->getAnsatzFunctionType( ) == p1dis &&
      this->space->getTestFunctionType( ) == p1dis ) {

    this->assembleP1DisP1Dis( leftCluster, rightCluster, matrix,
        voidIntegrator );
    return;
  }

  vector<LO>* innerElems = rightCluster->elems;
  vector<LO>* outerElems = leftCluster->elems;

  // allocate local element matrices
  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  matrix.resize( this->space->getClusterOuterDOFs( leftCluster ),
      this->space->getClusterInnerDOFs( rightCluster ) );

  LO iMax = outerElems->size( );
  LO jMax = innerElems->size( );

  BEIntegratorHelmholtz<LO, SC> * integrator =
      static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );

  // vector of indices where to put values in the global matrix
  vector<LO> rowIndices( nLocalRows );
  vector<LO> colIndices( nLocalCols );
  FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );

  for ( LO j = 0; j < jMax; j++ ) {
    for ( LO i = 0; i < iMax; i++ ) {
      integrator->getElemMatrixHypersingular( ( *outerElems )[i],
          ( *innerElems )[j], elemMatrix );
      this->space->getClusterOuterElemDOFs( leftCluster, i, &rowIndices[0] );
      this->space->getClusterInnerElemDOFs( rightCluster, j, &colIndices[0] );

      matrix.addToPositions( rowIndices, colIndices, elemMatrix );

    }
  }

}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleP1DisP1Dis(
    const BECluster< LO, SC > * leftCluster,
    const BECluster< LO, SC > * rightCluster,
    FullMatrix< LO, SC > & matrix,
    void * voidIntegrator
    ) const {

  BEIntegratorHelmholtz<LO, SC> * integrator =
      static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );

  const vector<LO> * innerElems = rightCluster->elems;
  const vector<LO> * outerElems = leftCluster->elems;
  LO iMax = outerElems->size( );
  LO jMax = innerElems->size( );

  matrix.resize( 3 * iMax, 3 * jMax );
  SC * matrixData = matrix.getData( );

  FullMatrix<LO, SC> elemMatrix( 3, 3 );
  SC * elemMatrixData = elemMatrix.getData( );

  for ( LO i = 0; i < iMax; ++i ) {
    for ( LO j = 0; j < jMax; ++j ) {

      integrator->getElemMatrixHypersingular( outerElems->at( i ),
          innerElems->at( j ), elemMatrix );
      for ( int oRot = 0; oRot < 3; ++oRot ) {
        for ( int iRot = 0; iRot < 3; ++iRot ) {
          matrixData[ ( 3 * j + iRot ) * 3 * iMax + 3 * i + oRot ] =
              elemMatrixData[ iRot * 3 + oRot ];
        }
      }
    }
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleRow(
    const BEBlockCluster< LO, SC > & block,
    LO idx,
    const std::vector<LO> & cols,
    Vector<LO, SC> & row,
    void * voidIntegrator
    ) const {

#pragma omp atomic update
  ++this->counterRow;

  if ( this->space->getTestFunctionType( ) == p1dis &&
      this->space->getAnsatzFunctionType( ) == p1dis ) {
    this->assembleRowP1DisP1Dis( block, idx, cols, row, voidIntegrator );
    return;
  }

  BEIntegratorHelmholtz<LO, SC> * integrator =
      static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );

  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );
  // vector of indices where to put values in the global matrix
  vector<LO> rowIndices( nLocalRows );
  vector<LO> colIndices( nLocalCols );

  LO nCols = cols.size( );
  //row.resize( nCols, true );
  row.setAll( 0.0 );
  vector< LO > testSupport;
  vector< LO > ansatzSupport;
  this->space->getOuterSupport( idx, testSupport,
      *( block.leftCluster->elems ) );

  for ( LO i = 0; i < testSupport.size( ); i++ ) {
    this->space->getOuterElemDOFs( testSupport[i], &rowIndices[0] );
    LO localRowIdx = 0;
    for ( LO auxI = 0; auxI < nLocalRows; auxI++ ) {
      if ( rowIndices[auxI] == idx ) {
        localRowIdx = auxI;
      }
    }
    for ( LO j = 0; j < nCols; j++ ) {
      this->space->getInnerSupport( cols[ j ], ansatzSupport,
          *( block.rightCluster->elems ) );

      for ( LO k = 0; k < ansatzSupport.size( ); k++ ) {
        integrator->getElemMatrixHypersingular( testSupport[i],
            ansatzSupport[k], elemMatrix );

        this->space->getInnerElemDOFs( ansatzSupport[k], &colIndices[0] );
        LO localColIdx = 0;
        for ( LO auxJ = 0; auxJ < nLocalCols; auxJ++ ) {
          if ( colIndices[auxJ] == cols[ j ] ) {
            localColIdx = auxJ;
          }
        }
        row.add( j, elemMatrix.get( localRowIdx, localColIdx ) );
      }
    }
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleRowP1DisP1Dis(
    const BEBlockCluster< LO, SC > & block,
    LO idx,
    const std::vector<LO> & cols,
    Vector<LO, SC> & row,
    void * voidIntegrator
    ) const {

  BEIntegratorHelmholtz<LO, SC> * integrator =
      static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );

  FullMatrix<LO, SC> elemMatrix( 3, 3 );

  for ( LO i = 0; i < block.rightCluster->nelems; ++i ) {
    integrator->getElemMatrixHypersingular( idx / 3,
        block.rightCluster->elems->at( i ), elemMatrix );
    row.set( 3 * i, elemMatrix.get( idx % 3, 0 ) );
    row.set( 3 * i + 1, elemMatrix.get( idx % 3, 1 ) );
    row.set( 3 * i + 2, elemMatrix.get( idx % 3, 2 ) );
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleColumn(
    const BEBlockCluster< LO, SC > & block,
    LO idx,
    const std::vector< LO > & rows,
    Vector<LO, SC>& col,
    void * voidIntegrator
    ) const {

#pragma omp atomic update
  ++this->counterCol;

  if ( this->space->getTestFunctionType( ) == p1dis &&
      this->space->getAnsatzFunctionType( ) == p1dis ) {
    this->assembleColumnP1DisP1Dis( block, idx, rows, col, voidIntegrator );
    return;
  }

  BEIntegratorHelmholtz<LO, SC> * integrator =
      static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );

  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );
  // vector of indices where to put values in the global matrix
  vector<LO> rowIndices( nLocalRows );
  vector<LO> colIndices( nLocalCols );

  LO nRows = rows.size( );
  //col.resize( nRows, true );
  col.setAll( 0.0 );
  vector< LO > testSupport;
  vector< LO > ansatzSupport;
  this->space->getInnerSupport( idx, ansatzSupport,
      *( block.rightCluster->elems ) );

  for ( LO i = 0; i < ansatzSupport.size( ); i++ ) {
    this->space->getInnerElemDOFs( ansatzSupport[i], &colIndices[0] );
    LO localColIdx = 0;
    for ( LO auxI = 0; auxI < nLocalCols; auxI++ ) {
      if ( colIndices[auxI] == idx ) {
        localColIdx = auxI;
      }
    }
    for ( LO j = 0; j < nRows; j++ ) {
      this->space->getOuterSupport( rows[ j ], testSupport,
          *( block.leftCluster->elems ) );

      for ( LO k = 0; k < testSupport.size( ); k++ ) {
        integrator->getElemMatrixHypersingular( testSupport[k],
            ansatzSupport[i], elemMatrix );

        this->space->getOuterElemDOFs( testSupport[k], &rowIndices[0] );
        LO localRowIdx = 0;
        for ( LO auxJ = 0; auxJ < nLocalRows; auxJ++ ) {
          if ( rowIndices[auxJ] == rows[ j ] ) {
            localRowIdx = auxJ;
          }
        }
        col.add( j, elemMatrix.get( localRowIdx, localColIdx ) );
      }
    }
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleColumnP1DisP1Dis(
    const BEBlockCluster< LO, SC > & block,
    LO idx,
    const std::vector<LO> & rows,
    Vector<LO, SC> & col,
    void * voidIntegrator
    ) const {

  BEIntegratorHelmholtz<LO, SC> * integrator =
      static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );

  FullMatrix<LO, SC> elemMatrix( 3, 3 );

  for ( LO i = 0; i < block.leftCluster->nelems; ++i ) {
    integrator->getElemMatrixHypersingular(
        block.leftCluster->elems->at( i ), idx / 3, elemMatrix );
    col.set( 3 * i, elemMatrix.get( 0, idx % 3 ) );
    col.set( 3 * i + 1, elemMatrix.get( 1, idx % 3 ) );
    col.set( 3 * i + 2, elemMatrix.get( 2, idx % 3 ) );
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assemble(
    FullMatrix<LO, SC> & matrix
    ) const {

  if ( this->space->getAnsatzFunctionType( ) == p1 &&
      this->space->getTestFunctionType( ) == p1 ) {

    this->assembleP1P1( matrix );
    return;
  }

  vector<LO> innerElems = this->space->getInnerElems( );
  vector<LO> outerElems = this->space->getOuterElems( );

  // allocate local element matrices
  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  matrix.resize( this->space->getOuterDOFs( ),
      this->space->getInnerDOFs( ) );

  LO iMax = outerElems.size( );
  LO jMax = innerElems.size( );

#ifdef VERBOSE
  ProgressMonitor::init( "Assembling Helmholtz D", iMax * jMax );
#endif

#pragma omp parallel shared(matrix)
  {
    BEIntegratorHelmholtz<LO, SC> integrator( this->space,
        this->quadratureOrder, this->kappa, this->quadrature,
        this->quadratureOrderDisjointElems );
    // vector of indices where to put values in the global matrix
    FullMatrix<LO, SC>* matrixBuffer[BUFFER_SIZE];
    vector<LO>* rowIndicesBuffer[BUFFER_SIZE];
    vector<LO>* colIndicesBuffer[BUFFER_SIZE];
    for ( int i = 0; i < BUFFER_SIZE; i++ ) {
      matrixBuffer[i] = new FullMatrix<LO, SC>( nLocalRows, nLocalCols );
      rowIndicesBuffer[i] = new vector<LO>( nLocalRows );
      colIndicesBuffer[i] = new vector<LO>( nLocalCols );
    }
    int counter = 0;

#pragma omp for collapse(2)
    for ( LO i = 0; i < iMax; i++ ) {
      for ( LO j = 0; j < jMax; j++ ) {
        integrator.getElemMatrixHypersingular( outerElems[i], innerElems[j],
            *( matrixBuffer[counter] ) );

        this->space->getOuterElemDOFs( outerElems[i],
            &( *( rowIndicesBuffer[counter] ) )[0] );
        this->space->getInnerElemDOFs( innerElems[j],
            &( *( colIndicesBuffer[counter] ) )[0] );
        counter++;
        if ( counter == BUFFER_SIZE ) {
          counter = 0;

          for ( int k = 0; k < BUFFER_SIZE; k++ ) {
            matrix.addToPositionsAtomic( *( rowIndicesBuffer[k] ),
                *( colIndicesBuffer[k] ), *( matrixBuffer[k] ) );
          }

#ifdef VERBOSE
          ProgressMonitor::step( BUFFER_SIZE );
#endif
        }
      }
    }

    for ( int i = 0; i < counter; i++ ) {
      matrix.addToPositionsAtomic( *( rowIndicesBuffer[i] ),
          *( colIndicesBuffer[i] ), *( matrixBuffer[i] ) );
    }

#ifdef VERBOSE
    ProgressMonitor::step( counter );
#endif
    for ( int i = 0; i < BUFFER_SIZE; i++ ) {
      delete matrixBuffer[i];
      delete rowIndicesBuffer[i];
      delete colIndicesBuffer[i];
    }
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleP1P1(
    FullMatrix<LO, SC> & matrix
    ) const {

  LO nNodes = this->space->getMesh( )->getNNodes( );
  LO nElems = this->space->getMesh( )->getNElements( );

  matrix.resize( nNodes, nNodes, false );
  SC * matrixData = matrix.getData( );

#pragma omp parallel for schedule( static, 32 )
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nNodes; ++i ) {
      matrixData[ j * nNodes + i ] = 0.0;
    }
  }

#pragma omp parallel
  {
    BEIntegratorHelmholtz<LO, SC> integrator( this->space,
        this->quadratureOrder, this->kappa, this->quadrature,
        this->quadratureOrderDisjointElems );

    FullMatrix< LO, SC > elemMatrix( 3, 3 );
    SC * elemMatrixData = elemMatrix.getData( );
    LO iElem[ 3 ], oElem[ 3 ];

#pragma omp for collapse( 2 ) schedule( dynamic, 32 )
    for ( LO i = 0; i < nElems; ++i ) { // outer (row)
      for ( LO j = 0; j < nElems; ++j ) { // inner (col)

        integrator.getElemMatrixHypersingular( i, j, elemMatrix );

        this->space->getMesh( )->getElement( i, oElem );
        this->space->getMesh( )->getElement( j, iElem );

        for ( int oRot = 0; oRot < 3; ++oRot ) {
          for ( int iRot = 0; iRot < 3; ++iRot ) {
#pragma omp atomic update
            reinterpret_cast<SCVT *> ( matrixData )
                [ 2 * ( iElem[ iRot ] * nNodes + oElem[ oRot ] ) ] +=
                std::real( elemMatrixData[ iRot * 3 + oRot ] );
#pragma omp atomic update
            reinterpret_cast<SCVT *> ( matrixData )
                [ 2 * ( iElem[ iRot ] * nNodes + oElem[ oRot ] ) + 1 ] +=
                std::imag( elemMatrixData[ iRot * 3 + oRot ] );
          }
        }
      }
    }
  } // end omp parallel
}
/*
template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assemble( FullMatrix<LO, SC>& matrix, FullMatrix<LO, SC>& V ) {

  vector<LO> innerElems = this->space->getInnerElems( );
  vector<LO> outerElems = this->space->getOuterElems( );

  // allocate local element matrices
  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  matrix.resize( this->space->getOuterDOFs( ),
      this->space->getInnerDOFs( ) );

  LO iMax = outerElems.size( );
  LO jMax = innerElems.size( );

#pragma omp parallel shared(matrix)
  {
    BEIntegratorHelmholtz<LO, SC> integrator( this->space, this->quadratureOrder,
        this->kappa, this->quadrature, this->quadratureOrderDisjointElems );
    // vector of indices where to put values in the global matrix
    vector<LO> rowIndices( nLocalRows );
    vector<LO> colIndices( nLocalCols );
    FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );

#pragma omp for collapse(2)
    for ( LO i = 0; i < iMax; i++ ) {
      for ( LO j = 0; j < jMax; j++ ) {
        integrator.getElemMatrixHypersingular( outerElems[i], innerElems[j], V, elemMatrix );
        this->space->getOuterElemDOFs( outerElems[i], &rowIndices[0] );
        this->space->getInnerElemDOFs( innerElems[j], &colIndices[0] );

        matrix.addToPositionsAtomic( rowIndices, colIndices, elemMatrix );

      }
    }
  }
}
*/

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleP0P0(
    FullMatrix<LO, SC>& matrix
    ) const {

  LO nElems = this->space->getMesh( )->getNElements( );

  this->assembleH2P0P0( matrix );

  Vector< LO, SC > ones( nElems );
  ones.setAll( 1.0 );
  Vector< LO, SC > diag( nElems );
  matrix.apply( ones, diag );
  matrix.scale( -1.0 );
  for ( LO i = 0; i < nElems; ++i ) {
    matrix.add( i, i, diag.get( i ) );
  }

  FullMatrix< LO, SC > H1;
  this->assembleH1P0P0( H1 );

  H1.apply( ones, diag );
  for ( LO i = 0; i < nElems; ++i ) {
    matrix.add( i, i, -diag.get( i ) );
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleH1P0P0(
    FullMatrix<LO, SC>& matrix
    ) const {

  vector<LO> innerElems = this->space->getInnerElems( );
  vector<LO> outerElems = this->space->getOuterElems( );

  // allocate local element matrices
  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  matrix.resize( this->space->getOuterDOFs( ),
      this->space->getInnerDOFs( ) );

  LO iMax = outerElems.size( );
  LO jMax = innerElems.size( );

#ifdef VERBOSE
  ProgressMonitor::init( "Assembling Helmholtz H1", iMax * jMax );
#endif

#pragma omp parallel shared(matrix)
  {
    BEIntegratorHelmholtz<LO, SC> integrator( this->space,
        this->quadratureOrder, this->kappa, this->quadrature,
        this->quadratureOrderDisjointElems );
    // vector of indices where to put values in the global matrix
    FullMatrix<LO, SC>* matrixBuffer[BUFFER_SIZE];
    vector<LO>* rowIndicesBuffer[BUFFER_SIZE];
    vector<LO>* colIndicesBuffer[BUFFER_SIZE];
    for ( int i = 0; i < BUFFER_SIZE; i++ ) {
      matrixBuffer[i] = new FullMatrix<LO, SC>( nLocalRows, nLocalCols );
      rowIndicesBuffer[i] = new vector<LO>( nLocalRows );
      colIndicesBuffer[i] = new vector<LO>( nLocalCols );
    }
    int counter = 0;

#pragma omp for collapse(2)
    for ( LO i = 0; i < iMax; i++ ) {
      //#pragma omp for
      for ( LO j = 0; j < jMax; j++ ) {
        integrator.computeElemMatrixH1SauterSchwabP0P0( outerElems[i],
            innerElems[j], *( matrixBuffer[counter] ) );
        this->space->getOuterElemDOFs( outerElems[i],
            &( *( rowIndicesBuffer[counter] ) )[0] );
        this->space->getInnerElemDOFs( innerElems[j],
            &( *( colIndicesBuffer[counter] ) )[0] );
        counter++;

        if ( counter == BUFFER_SIZE ) {
          counter = 0;

          for ( int k = 0; k < BUFFER_SIZE; k++ ) {
            matrix.addToPositionsAtomic( *( rowIndicesBuffer[k] ),
                *( colIndicesBuffer[k] ), *( matrixBuffer[k] ) );
          }

#ifdef VERBOSE
          ProgressMonitor::step( BUFFER_SIZE );
#endif
        }
      }
    }

    for ( int i = 0; i < counter; i++ ) {
      matrix.addToPositionsAtomic( *( rowIndicesBuffer[i] ),
          *( colIndicesBuffer[i] ), *( matrixBuffer[i] ) );
    }

#ifdef VERBOSE
    ProgressMonitor::step( counter );
#endif
    for ( int i = 0; i < BUFFER_SIZE; i++ ) {
      delete matrixBuffer[i];
      delete rowIndicesBuffer[i];
      delete colIndicesBuffer[i];
    }
  }

}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleH2P0P0(
    FullMatrix<LO, SC>& matrix
    ) const {

  vector<LO> innerElems = this->space->getInnerElems( );
  vector<LO> outerElems = this->space->getOuterElems( );

  // allocate local element matrices
  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  matrix.resize( this->space->getOuterDOFs( ),
      this->space->getInnerDOFs( ) );

  LO iMax = outerElems.size( );
  LO jMax = innerElems.size( );

#ifdef VERBOSE
  ProgressMonitor::init( "Assembling Helmholtz H2", iMax * jMax );
#endif

#pragma omp parallel shared(matrix)
  {
    BEIntegratorHelmholtz<LO, SC> integrator( this->space,
        this->quadratureOrder, this->kappa, this->quadrature,
        this->quadratureOrderDisjointElems );
    // vector of indices where to put values in the global matrix
    FullMatrix<LO, SC>* matrixBuffer[BUFFER_SIZE];
    vector<LO>* rowIndicesBuffer[BUFFER_SIZE];
    vector<LO>* colIndicesBuffer[BUFFER_SIZE];
    for ( int i = 0; i < BUFFER_SIZE; i++ ) {
      matrixBuffer[i] = new FullMatrix<LO, SC>( nLocalRows, nLocalCols );
      rowIndicesBuffer[i] = new vector<LO>( nLocalRows );
      colIndicesBuffer[i] = new vector<LO>( nLocalCols );
    }
    int counter = 0;

#pragma omp for collapse(2)
    for ( LO i = 0; i < iMax; i++ ) {
      //#pragma omp for
      for ( LO j = 0; j < jMax; j++ ) {
        integrator.computeElemMatrixH2SauterSchwabP0P0( outerElems[i], innerElems[j], *( matrixBuffer[counter] ) );
        this->space->getOuterElemDOFs( outerElems[i], &( *( rowIndicesBuffer[counter] ) )[0] );
        this->space->getInnerElemDOFs( innerElems[j], &( *( colIndicesBuffer[counter] ) )[0] );
        counter++;

        if ( counter == BUFFER_SIZE ) {
          counter = 0;

          for ( int k = 0; k < BUFFER_SIZE; k++ ) {
            matrix.addToPositionsAtomic( *( rowIndicesBuffer[k] ),
                *( colIndicesBuffer[k] ), *( matrixBuffer[k] ) );
          }

#ifdef VERBOSE
          ProgressMonitor::step( BUFFER_SIZE );
#endif
        }
      }
    }

    for ( int i = 0; i < counter; i++ ) {
      matrix.addToPositionsAtomic( *( rowIndicesBuffer[i] ),
          *( colIndicesBuffer[i] ), *( matrixBuffer[i] ) );
    }

#ifdef VERBOSE
    ProgressMonitor::step( counter );
#endif
    for ( int i = 0; i < BUFFER_SIZE; i++ ) {
      delete matrixBuffer[i];
      delete rowIndicesBuffer[i];
      delete colIndicesBuffer[i];
    }
  }

}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleH1P0P0(
    ACAMatrix<LO, SC>& matrix
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

#ifdef VERBOSE
  ProgressMonitor::init( "Assembling nonadmissible blocks", nNonAdm );
#endif

#pragma omp parallel 
  {

    BEBlockCluster<LO, SC>* block;
    BECluster<LO, SC>* leftCluster;
    BECluster<LO, SC>* rightCluster;

#pragma omp for schedule(dynamic)
    for ( int i = 0; i < nNonAdm; i++ ) {
      block = ( *nonAdmLeaves )[i];
      leftCluster = block->leftCluster;
      rightCluster = block->rightCluster;
      FullMatrix<LO, SC> *fullBlock = new FullMatrix<LO, SC>( );
      assembleH1P0P0( leftCluster, rightCluster, *fullBlock );

      matrix.addNonadmissibleBlock( fullBlock, i );
#ifdef VERBOSE
      ProgressMonitor::step( );
#endif

      //#pragma omp atomic
      //      currIt++;
    }
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
      this->assembleACABlockH1P0P0( block, *acaU, *acaV, voidIntegrator );


      //      if ( acaU->getNRows( ) == 1 && acaU->getNCols( ) == 1 ) {
      //        
      //        // if ACA could not approximate block, assemble it as a full 
      //        
      //        FullMatrix<LO, SC> *fullBlock = new FullMatrix<LO, SC>( );
      //        assemble( block->leftCluster, block->rightCluster, *fullBlock );
      //        matrix.addAdmissibleBlock(
      //            std::pair<FullMatrix<LO, SC>*, FullMatrix<LO, SC>* >( fullBlock,
      //            nullptr ),
      //            i );
      //        delete acaU;
      //        delete acaV;        
      //      } else 
      if ( acaU->getNRows( ) * acaV->getNRows( ) <
          acaU->getNRows( ) * acaU->getNCols( ) +
          acaV->getNRows( ) * acaV->getNCols( ) ) {

        // if the size of UV' < (size of U) + (size of V) multiply matrices

        FullMatrix<LO, SC> *UV = new FullMatrix<LO, SC>(
            acaU->getNRows( ), acaV->getNRows( ) );
        UV->multiply( *acaU, *acaV, false, true );

        matrix.addAdmissibleBlock(
            std::pair<FullMatrix<LO, SC>*, FullMatrix<LO, SC>* >( UV, nullptr ),
            i );

        delete acaU;
        delete acaV;
      } else {
        matrix.addAdmissibleBlock(
            std::pair<FullMatrix<LO, SC>*, FullMatrix<LO, SC>* >( acaU, acaV ),
            i );
      }

#ifdef VERBOSE
      ProgressMonitor::step( );
#endif
    }

    this->destroyIntegrator( voidIntegrator );
    ;
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleH2P0P0(
    ACAMatrix<LO, SC>& matrix
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

#ifdef VERBOSE
  ProgressMonitor::init( "Assembling nonadmissible blocks", nNonAdm );
#endif

#pragma omp parallel 
  {

    BEBlockCluster<LO, SC>* block;
    BECluster<LO, SC>* leftCluster;
    BECluster<LO, SC>* rightCluster;

#pragma omp for schedule(dynamic)
    for ( int i = 0; i < nNonAdm; i++ ) {
      block = ( *nonAdmLeaves )[i];
      leftCluster = block->leftCluster;
      rightCluster = block->rightCluster;
      FullMatrix<LO, SC> *fullBlock = new FullMatrix<LO, SC>( );
      assembleH2P0P0( leftCluster, rightCluster, *fullBlock );

      matrix.addNonadmissibleBlock( fullBlock, i );
#ifdef VERBOSE
      ProgressMonitor::step( );
#endif

      //#pragma omp atomic
      //      currIt++;
    }
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
      this->assembleACABlockH2P0P0( block, *acaU, *acaV, voidIntegrator );


      //      if ( acaU->getNRows( ) == 1 && acaU->getNCols( ) == 1 ) {
      //        
      //        // if ACA could not approximate block, assemble it as a full 
      //        
      //        FullMatrix<LO, SC> *fullBlock = new FullMatrix<LO, SC>( );
      //        assemble( block->leftCluster, block->rightCluster, *fullBlock );
      //        matrix.addAdmissibleBlock(
      //            std::pair<FullMatrix<LO, SC>*, FullMatrix<LO, SC>* >( fullBlock,
      //            nullptr ),
      //            i );
      //        delete acaU;
      //        delete acaV;        
      //      } else 
      if ( acaU->getNRows( ) * acaV->getNRows( ) <
          acaU->getNRows( ) * acaU->getNCols( ) +
          acaV->getNRows( ) * acaV->getNCols( ) ) {

        // if the size of UV' < (size of U) + (size of V) multiply matrices

        FullMatrix<LO, SC> *UV = new FullMatrix<LO, SC>(
            acaU->getNRows( ), acaV->getNRows( ) );
        UV->multiply( *acaU, *acaV, false, true );

        matrix.addAdmissibleBlock(
            std::pair<FullMatrix<LO, SC>*, FullMatrix<LO, SC>* >( UV, nullptr ),
            i );

        delete acaU;
        delete acaV;
      } else {
        matrix.addAdmissibleBlock(
            std::pair<FullMatrix<LO, SC>*, FullMatrix<LO, SC>* >( acaU, acaV ),
            i );
      }

#ifdef VERBOSE
      ProgressMonitor::step( );
#endif
    }

    this->destroyIntegrator( voidIntegrator );
    ;
  }
}

template <class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleH1P0P0(
    BECluster<LO, SC> const *leftCluster,
    BECluster<LO, SC> const *rightCluster,
    FullMatrix<LO, SC>& matrix
    ) const {

  vector<LO>* innerElems = rightCluster->elems;
  vector<LO>* outerElems = leftCluster->elems;

  // allocate local element matrices
  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  matrix.resize( this->space->getClusterOuterDOFs( leftCluster ),
      this->space->getClusterInnerDOFs( rightCluster ) );

  LO iMax = outerElems->size( );
  LO jMax = innerElems->size( );

  BEIntegratorHelmholtz<LO, SC> integrator( this->space, this->quadratureOrder,
      kappa, this->quadrature, this->quadratureOrderDisjointElems );
  // vector of indices where to put values in the global matrix
  vector<LO> rowIndices( nLocalRows );
  vector<LO> colIndices( nLocalCols );
  FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );

  for ( LO j = 0; j < jMax; j++ ) {
    for ( LO i = 0; i < iMax; i++ ) {
      integrator.computeElemMatrixH1SauterSchwabP0P0( ( *outerElems )[i],
          ( *innerElems )[j], elemMatrix );
      this->space->getClusterOuterElemDOFs( leftCluster, i, &rowIndices[0] );
      this->space->getClusterInnerElemDOFs( rightCluster, j, &colIndices[0] );

      matrix.addToPositions( rowIndices, colIndices, elemMatrix );

    }
  }

}

template <class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleH2P0P0(
    BECluster<LO, SC> const *leftCluster,
    BECluster<LO, SC> const *rightCluster,
    FullMatrix<LO, SC>& matrix
    ) const {

  vector<LO>* innerElems = rightCluster->elems;
  vector<LO>* outerElems = leftCluster->elems;

  // allocate local element matrices
  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  matrix.resize( this->space->getClusterOuterDOFs( leftCluster ),
      this->space->getClusterInnerDOFs( rightCluster ) );

  LO iMax = outerElems->size( );
  LO jMax = innerElems->size( );

  BEIntegratorHelmholtz<LO, SC> integrator( this->space, this->quadratureOrder,
      kappa, this->quadrature, this->quadratureOrderDisjointElems );
  // vector of indices where to put values in the global matrix
  vector<LO> rowIndices( nLocalRows );
  vector<LO> colIndices( nLocalCols );
  FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );

  for ( LO j = 0; j < jMax; j++ ) {
    for ( LO i = 0; i < iMax; i++ ) {
      integrator.computeElemMatrixH2SauterSchwabP0P0( ( *outerElems )[i],
          ( *innerElems )[j], elemMatrix );
      this->space->getClusterOuterElemDOFs( leftCluster, i, &rowIndices[0] );
      this->space->getClusterInnerElemDOFs( rightCluster, j, &colIndices[0] );

      matrix.addToPositions( rowIndices, colIndices, elemMatrix );

    }
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleRowH1P0P0(
    const BEBlockCluster< LO, SC > & block,
    LO idx,
    const std::vector<LO> & cols,
    Vector<LO, SC> & row,
    void * voidIntegrator
    ) const {

  BEIntegratorHelmholtz<LO, SC> * integrator =
      static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );

  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );
  // vector of indices where to put values in the global matrix
  vector<LO> rowIndices( nLocalRows );
  vector<LO> colIndices( nLocalCols );

  LO nCols = cols.size( );
  // row.resize( nCols, true );
  row.setAll( 0.0 );
  vector< LO > testSupport;
  vector< LO > ansatzSupport;
  this->space->getOuterSupport( idx, testSupport,
      *( block.leftCluster->elems ) );

  for ( LO i = 0; i < testSupport.size( ); i++ ) {
    this->space->getOuterElemDOFs( testSupport[i], &rowIndices[0] );
    LO localRowIdx = 0;
    for ( LO auxI = 0; auxI < nLocalRows; auxI++ ) {
      if ( rowIndices[auxI] == idx ) {
        localRowIdx = auxI;
      }
    }
    for ( LO j = 0; j < nCols; j++ ) {
      this->space->getInnerSupport( cols[ j ], ansatzSupport,
          *( block.rightCluster->elems ) );

      for ( LO k = 0; k < ansatzSupport.size( ); k++ ) {
        integrator->computeElemMatrixH1SauterSchwabP0P0( testSupport[i],
            ansatzSupport[k], elemMatrix );

        this->space->getInnerElemDOFs( ansatzSupport[k], &colIndices[0] );
        LO localColIdx = 0;
        for ( LO auxJ = 0; auxJ < nLocalCols; auxJ++ ) {
          if ( colIndices[auxJ] == cols[ j ] ) {
            localColIdx = auxJ;
          }
        }
        row.add( j, elemMatrix.get( localRowIdx, localColIdx ) );
      }
    }
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleColumnH1P0P0(
    const BEBlockCluster< LO, SC > & block,
    LO idx,
    const std::vector< LO > & rows,
    Vector<LO, SC>& col,
    void * voidIntegrator
    ) const {

  BEIntegratorHelmholtz<LO, SC> * integrator =
      static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );

  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );
  // vector of indices where to put values in the global matrix
  vector<LO> rowIndices( nLocalRows );
  vector<LO> colIndices( nLocalCols );

  LO nRows = rows.size( );
  //col.resize( nRows, true );
  col.setAll( 0.0 );
  vector< LO > testSupport;
  vector< LO > ansatzSupport;
  this->space->getInnerSupport( idx, ansatzSupport,
      *( block.rightCluster->elems ) );

  for ( LO i = 0; i < ansatzSupport.size( ); i++ ) {
    this->space->getInnerElemDOFs( ansatzSupport[i], &colIndices[0] );
    LO localColIdx = 0;
    for ( LO auxI = 0; auxI < nLocalCols; auxI++ ) {
      if ( colIndices[auxI] == idx ) {
        localColIdx = auxI;
      }
    }
    for ( LO j = 0; j < nRows; j++ ) {
      this->space->getOuterSupport( rows[ j ], testSupport,
          *( block.leftCluster->elems ) );

      for ( LO k = 0; k < testSupport.size( ); k++ ) {
        integrator->computeElemMatrixH1SauterSchwabP0P0( testSupport[k],
            ansatzSupport[i], elemMatrix );

        this->space->getOuterElemDOFs( testSupport[k], &rowIndices[0] );
        LO localRowIdx = 0;
        for ( LO auxJ = 0; auxJ < nLocalRows; auxJ++ ) {
          if ( rowIndices[auxJ] == rows[ j ] ) {
            localRowIdx = auxJ;
          }
        }
        col.add( j, elemMatrix.get( localRowIdx, localColIdx ) );
      }
    }
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleRowH2P0P0(
    const BEBlockCluster< LO, SC > & block,
    LO idx,
    const std::vector<LO> & cols,
    Vector<LO, SC> & row,
    void * voidIntegrator
    ) const {

  BEIntegratorHelmholtz<LO, SC> * integrator =
      static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );

  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );
  // vector of indices where to put values in the global matrix
  vector<LO> rowIndices( nLocalRows );
  vector<LO> colIndices( nLocalCols );

  LO nCols = cols.size( );
  // row.resize( nCols, true );
  row.setAll( 0.0 );
  vector< LO > testSupport;
  vector< LO > ansatzSupport;
  this->space->getOuterSupport( idx, testSupport,
      *( block.leftCluster->elems ) );

  for ( LO i = 0; i < testSupport.size( ); i++ ) {
    this->space->getOuterElemDOFs( testSupport[i], &rowIndices[0] );
    LO localRowIdx = 0;
    for ( LO auxI = 0; auxI < nLocalRows; auxI++ ) {
      if ( rowIndices[auxI] == idx ) {
        localRowIdx = auxI;
      }
    }
    for ( LO j = 0; j < nCols; j++ ) {
      this->space->getInnerSupport( cols[ j ], ansatzSupport,
          *( block.rightCluster->elems ) );

      for ( LO k = 0; k < ansatzSupport.size( ); k++ ) {
        integrator->computeElemMatrixH2SauterSchwabP0P0( testSupport[i],
            ansatzSupport[k], elemMatrix );

        this->space->getInnerElemDOFs( ansatzSupport[k], &colIndices[0] );
        LO localColIdx = 0;
        for ( LO auxJ = 0; auxJ < nLocalCols; auxJ++ ) {
          if ( colIndices[auxJ] == cols[ j ] ) {
            localColIdx = auxJ;
          }
        }
        row.add( j, elemMatrix.get( localRowIdx, localColIdx ) );
      }
    }
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleColumnH2P0P0(
    const BEBlockCluster< LO, SC > & block,
    LO idx,
    const std::vector< LO > & rows,
    Vector<LO, SC>& col,
    void * voidIntegrator
    ) const {

  BEIntegratorHelmholtz<LO, SC> * integrator =
      static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );

  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );
  // vector of indices where to put values in the global matrix
  vector<LO> rowIndices( nLocalRows );
  vector<LO> colIndices( nLocalCols );

  LO nRows = rows.size( );
  //col.resize( nRows, true );
  col.setAll( 0.0 );
  vector< LO > testSupport;
  vector< LO > ansatzSupport;
  this->space->getInnerSupport( idx, ansatzSupport,
      *( block.rightCluster->elems ) );

  for ( LO i = 0; i < ansatzSupport.size( ); i++ ) {
    this->space->getInnerElemDOFs( ansatzSupport[i], &colIndices[0] );
    LO localColIdx = 0;
    for ( LO auxI = 0; auxI < nLocalCols; auxI++ ) {
      if ( colIndices[auxI] == idx ) {
        localColIdx = auxI;
      }
    }
    for ( LO j = 0; j < nRows; j++ ) {
      this->space->getOuterSupport( rows[ j ], testSupport,
          *( block.leftCluster->elems ) );

      for ( LO k = 0; k < testSupport.size( ); k++ ) {
        integrator->computeElemMatrixH2SauterSchwabP0P0( testSupport[k],
            ansatzSupport[i], elemMatrix );

        this->space->getOuterElemDOFs( testSupport[k], &rowIndices[0] );
        LO localRowIdx = 0;
        for ( LO auxJ = 0; auxJ < nLocalRows; auxJ++ ) {
          if ( rowIndices[auxJ] == rows[ j ] ) {
            localRowIdx = auxJ;
          }
        }
        col.add( j, elemMatrix.get( localRowIdx, localColIdx ) );
      }
    }
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleACABlockH1P0P0(
    BEBlockCluster<LO, SC>* block,
    FullMatrix<LO, SC> & U,
    FullMatrix<LO, SC> & V,
    void * voidIntegrator
    ) const {

  FastBESpace<LO, SC>* fastSpace =
      static_cast<FastBESpace<LO, SC>*> ( this->space );

  SCVT zeroEps = 1e-14 * fastSpace->getScaleACA( );

  BECluster<LO, SC>* leftCluster = block->leftCluster;
  BECluster<LO, SC>* rightCluster = block->rightCluster;

  LO nRows = fastSpace->getClusterOuterDOFs( leftCluster );
  LO nCols = fastSpace->getClusterInnerDOFs( rightCluster );
  //std::cout << nRows << " " << nCols<< std::endl;
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

  Vector<LO, SC> * ru = nullptr;
  Vector<LO, SC> * rv = nullptr;

  LO iters = 0;
  LO maxIters = 30;
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

  // stop if all rows and cols used
  while ( idxI.size( ) != nRows && idxJ.size( ) != nCols && iters < maxIters ) {

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
      Vector<LO, SC> * currentRow = new Vector<LO, SC>( nCols );
      this->assembleRowH1P0P0( *block, ( *block->outerDOFs )[pivotRow],
          *block->innerDOFs, *currentRow, voidIntegrator );

      idxI.push_back( pivotRow );
      idxICompl[ pivotRow ] = false;

      if ( std::abs( currentRow->norm2( ) ) < zeroEps ) {
        delete currentRow;
        searchForRow = true;
        ++iters;
        continue;
      }

      for ( LO i = 0; i < nCols; ++i ) {
        r.add( i, std::abs( currentRow->get( i ) ) );
      }

      rv = currentRow; //new Vector< LO, SC >( *currentRow );
      for ( LO i = 0; i < u.size( ); ++i ) {
        rv->add( *v[ i ], -u[ i ]->get( pivotRow ) );
      }
      pivotCol = rv->findAbsMax( ); // <-- vyjmout sloupec, ktery uz byl pouzit
      maxRowRes = rv->get( pivotCol );

      if ( std::abs( maxRowRes ) > zeroEps ) {
        Vector< LO, SC > * currentColumn =
            new Vector<LO, SC>( nRows );
        gamma = ( (SCVT) 1.0 ) / maxRowRes;
        this->assembleColumnH1P0P0( *block, ( *block->innerDOFs )[pivotCol],
            *block->outerDOFs, *currentColumn, voidIntegrator );

        idxJ.push_back( pivotCol );
        idxJCompl[ pivotCol ] = false;

        for ( LO i = 0; i < nRows; ++i ) {
          c.add( i, std::abs( currentColumn->get( i ) ) );
        }

        ru = currentColumn; //new Vector< LO, SC >( currentColumn );
        for ( LO i = 0; i < v.size( ); ++i ) {
          ru->add( *u[ i ], -v[ i ]->get( pivotCol ) );
        }

        pivotRow = ru->findAbsMax( );
        u.push_back( ru );
        rv->scale( gamma );
        v.push_back( rv );
        ++rank;
      } else {
        delete rv;
        rv = nullptr;
        ++iters;
        searchForRow = true;
        continue;
      }
    } else { // looking for column
      Vector< LO, SC > * currentColumn =
          new Vector<LO, SC>( nRows );

      this->assembleColumnH1P0P0( *block, ( *block->innerDOFs )[pivotCol],
          *block->outerDOFs, *currentColumn, voidIntegrator );
      idxJ.push_back( pivotCol );
      idxJCompl[ pivotCol ] = false;

      if ( currentColumn->norm2( ) < zeroEps ) {

        delete currentColumn;
        searchForRow = true;
        ++iters;
        continue;
      }

      for ( LO i = 0; i < nRows; ++i ) {
        c.add( i, std::abs( currentColumn->get( i ) ) );
      }

      ru = currentColumn; //new Vector< LO, SC >( currentColumn );
      for ( LO i = 0; i < v.size( ); ++i ) {
        ru->add( *u[ i ], -v[ i ]->get( pivotCol ) );
      }

      pivotRow = ru->findAbsMax( );
      maxColRes = ru->get( pivotRow );

      if ( std::abs( maxColRes ) > zeroEps ) {

        Vector< LO, SC > * currentRow =
            new Vector<LO, SC>( nCols );

        gamma = ( (SCVT) 1.0 ) / maxColRes;
        this->assembleRowH1P0P0( *block, ( *block->outerDOFs )[pivotRow], // originally pivotCol
            *block->innerDOFs, *currentRow, voidIntegrator );
        idxI.push_back( pivotRow );
        idxICompl[ pivotRow ] = false;

        for ( LO i = 0; i < nCols; ++i ) {
          r.add( i, std::abs( currentRow->get( i ) ) );
        }

        rv = currentRow; // new Vector< LO, SC >( currentRow );
        for ( LO i = 0; i < u.size( ); ++i ) {
          rv->add( *v[ i ], -u[ i ]->get( pivotRow ) );
        }

        pivotCol = rv->findAbsMax( );
        ru->scale( gamma );
        u.push_back( ru );
        v.push_back( rv );
        ++rank;
      } else {
        //        if ( u.size( ) >= preallocSize ) {
        //          delete [] currColData;
        //        }
        delete ru;
        ru = nullptr;
        searchForRow = true;
        ++iters;
        continue;
      }
    } // end column

    // check approximation
    //if ( maxRowRes > zeroEps && maxColRes > zeroEps ) {
    if ( ru && rv ) {
      normUV = ru->norm2( ) * rv->norm2( );
      squareNormS += normUV * normUV;
      for ( LO i = 0; i < u.size( ) - 1; ++i ) {
        squareNormS += 2.0 * std::real( u.back( )->dot( *u[ i ] ) *
            v.back( )->dot( *v[ i ] ) );
      }
    }
    // stopping criterion
    // this stopping criterion modified to sqrt( eps )!
    if ( ru && rv &&
        normUV > fastSpace->getEpsilonACA( ) * std::sqrt( squareNormS ) ) {
      //normUV > std::sqrt( fastSpace->getEpsilonACA( ) ) *
      //std::sqrt( squareNormS ) ) {
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

} // end of ACA()

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleACABlockH2P0P0(
    BEBlockCluster<LO, SC>* block,
    FullMatrix<LO, SC> & U,
    FullMatrix<LO, SC> & V,
    void * voidIntegrator
    ) const {

  FastBESpace<LO, SC>* fastSpace =
      static_cast<FastBESpace<LO, SC>*> ( this->space );

  SCVT zeroEps = 1e-14 * fastSpace->getScaleACA( );

  BECluster<LO, SC>* leftCluster = block->leftCluster;
  BECluster<LO, SC>* rightCluster = block->rightCluster;

  LO nRows = fastSpace->getClusterOuterDOFs( leftCluster );
  LO nCols = fastSpace->getClusterInnerDOFs( rightCluster );
  //std::cout << nRows << " " << nCols<< std::endl;
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

  Vector<LO, SC> * ru = nullptr;
  Vector<LO, SC> * rv = nullptr;

  LO iters = 0;
  LO maxIters = 30;
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

  // stop if all rows and cols used
  while ( idxI.size( ) != nRows && idxJ.size( ) != nCols && iters < maxIters ) {

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
      Vector<LO, SC> * currentRow = new Vector<LO, SC>( nCols );
      this->assembleRowH2P0P0( *block, ( *block->outerDOFs )[pivotRow],
          *block->innerDOFs, *currentRow, voidIntegrator );

      idxI.push_back( pivotRow );
      idxICompl[ pivotRow ] = false;

      if ( std::abs( currentRow->norm2( ) ) < zeroEps ) {
        delete currentRow;
        searchForRow = true;
        ++iters;
        continue;
      }

      for ( LO i = 0; i < nCols; ++i ) {
        r.add( i, std::abs( currentRow->get( i ) ) );
      }

      rv = currentRow; //new Vector< LO, SC >( *currentRow );
      for ( LO i = 0; i < u.size( ); ++i ) {
        rv->add( *v[ i ], -u[ i ]->get( pivotRow ) );
      }
      pivotCol = rv->findAbsMax( ); // <-- vyjmout sloupec, ktery uz byl pouzit
      maxRowRes = rv->get( pivotCol );

      if ( std::abs( maxRowRes ) > zeroEps ) {
        Vector< LO, SC > * currentColumn =
            new Vector<LO, SC>( nRows );
        gamma = ( (SCVT) 1.0 ) / maxRowRes;
        this->assembleColumnH2P0P0( *block, ( *block->innerDOFs )[pivotCol],
            *block->outerDOFs, *currentColumn, voidIntegrator );

        idxJ.push_back( pivotCol );
        idxJCompl[ pivotCol ] = false;

        for ( LO i = 0; i < nRows; ++i ) {
          c.add( i, std::abs( currentColumn->get( i ) ) );
        }

        ru = currentColumn; //new Vector< LO, SC >( currentColumn );
        for ( LO i = 0; i < v.size( ); ++i ) {
          ru->add( *u[ i ], -v[ i ]->get( pivotCol ) );
        }

        pivotRow = ru->findAbsMax( );
        u.push_back( ru );
        rv->scale( gamma );
        v.push_back( rv );
        ++rank;
      } else {
        delete rv;
        rv = nullptr;
        ++iters;
        searchForRow = true;
        continue;
      }
    } else { // looking for column
      Vector< LO, SC > * currentColumn =
          new Vector<LO, SC>( nRows );

      this->assembleColumnH2P0P0( *block, ( *block->innerDOFs )[pivotCol],
          *block->outerDOFs, *currentColumn, voidIntegrator );
      idxJ.push_back( pivotCol );
      idxJCompl[ pivotCol ] = false;

      if ( currentColumn->norm2( ) < zeroEps ) {

        delete currentColumn;
        searchForRow = true;
        ++iters;
        continue;
      }

      for ( LO i = 0; i < nRows; ++i ) {
        c.add( i, std::abs( currentColumn->get( i ) ) );
      }

      ru = currentColumn; //new Vector< LO, SC >( currentColumn );
      for ( LO i = 0; i < v.size( ); ++i ) {
        ru->add( *u[ i ], -v[ i ]->get( pivotCol ) );
      }

      pivotRow = ru->findAbsMax( );
      maxColRes = ru->get( pivotRow );

      if ( std::abs( maxColRes ) > zeroEps ) {

        Vector< LO, SC > * currentRow =
            new Vector<LO, SC>( nCols );

        gamma = ( (SCVT) 1.0 ) / maxColRes;
        this->assembleRowH2P0P0( *block, ( *block->outerDOFs )[pivotRow], // originally pivotCol
            *block->innerDOFs, *currentRow, voidIntegrator );
        idxI.push_back( pivotRow );
        idxICompl[ pivotRow ] = false;

        for ( LO i = 0; i < nCols; ++i ) {
          r.add( i, std::abs( currentRow->get( i ) ) );
        }

        rv = currentRow; // new Vector< LO, SC >( currentRow );
        for ( LO i = 0; i < u.size( ); ++i ) {
          rv->add( *v[ i ], -u[ i ]->get( pivotRow ) );
        }

        pivotCol = rv->findAbsMax( );
        ru->scale( gamma );
        u.push_back( ru );
        v.push_back( rv );
        ++rank;
      } else {
        //        if ( u.size( ) >= preallocSize ) {
        //          delete [] currColData;
        //        }
        delete ru;
        ru = nullptr;
        searchForRow = true;
        ++iters;
        continue;
      }
    } // end column

    // check approximation
    //if ( maxRowRes > zeroEps && maxColRes > zeroEps ) {
    if ( ru && rv ) {
      normUV = ru->norm2( ) * rv->norm2( );
      squareNormS += normUV * normUV;
      for ( LO i = 0; i < u.size( ) - 1; ++i ) {
        squareNormS += 2.0 * std::real( u.back( )->dot( *u[ i ] ) *
            v.back( )->dot( *v[ i ] ) );
      }
    }
    // stopping criterion
    // this stopping criterion modified to sqrt( eps )!
    if ( ru && rv &&
        normUV > fastSpace->getEpsilonACA( ) * std::sqrt( squareNormS ) ) {
      //normUV > std::sqrt( fastSpace->getEpsilonACA( ) ) *
      //std::sqrt( squareNormS ) ) {
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

} // end of ACA()

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleH1P0P0(
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

#pragma omp for schedule(dynamic)
    for ( int i = 0; i < myNonAdmSize; i++ ) {
      block = ( *nonAdmLeaves )[myNonAdmLeaves[i]];
      leftCluster = block->leftCluster;
      rightCluster = block->rightCluster;
      FullMatrix<LO, SC> *fullBlock = new FullMatrix<LO, SC>( );
      assembleH1P0P0( leftCluster, rightCluster, *fullBlock );

      matrix.addNonadmissibleBlock( fullBlock, myNonAdmLeaves[i] );
#ifdef VERBOSE
      ProgressMonitor::step( );
#endif

      //#pragma omp atomic
      //      currIt++;
    }
  }

  // prepare the MPI loop over admissible blocks
  int myAdmSize = nAdm / size;
  if ( rank < nAdm % size ) myAdmSize++;
  int *myAdmLeaves = new int[myAdmSize];
  for ( int i = 0; i < myAdmSize; i++ ) {
    myAdmLeaves[i] = i * size + rank;
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
      this->assembleACABlockH1P0P0( block, *acaU, *acaV, voidIntegrator );

      if ( acaU->getNRows( ) * acaV->getNRows( ) <
          acaU->getNRows( ) * acaU->getNCols( ) +
          acaV->getNRows( ) * acaV->getNCols( ) ) {

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

    this->destroyIntegrator( voidIntegrator );
    ;
  }
  delete [] myAdmLeaves;
  delete [] myNonAdmLeaves;
}

template<class LO, class SC>
void BEBilinearFormHelmholtzHypersingular<LO, SC>::assembleH2P0P0(
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

#pragma omp for schedule(dynamic)
    for ( int i = 0; i < myNonAdmSize; i++ ) {
      block = ( *nonAdmLeaves )[myNonAdmLeaves[i]];
      leftCluster = block->leftCluster;
      rightCluster = block->rightCluster;
      FullMatrix<LO, SC> *fullBlock = new FullMatrix<LO, SC>( );
      assembleH2P0P0( leftCluster, rightCluster, *fullBlock );

      matrix.addNonadmissibleBlock( fullBlock, myNonAdmLeaves[i] );
#ifdef VERBOSE
      ProgressMonitor::step( );
#endif

      //#pragma omp atomic
      //      currIt++;
    }
  }

  // prepare the MPI loop over admissible blocks
  int myAdmSize = nAdm / size;
  if ( rank < nAdm % size ) myAdmSize++;
  int *myAdmLeaves = new int[myAdmSize];
  for ( int i = 0; i < myAdmSize; i++ ) {
    myAdmLeaves[i] = i * size + rank;
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
      this->assembleACABlockH2P0P0( block, *acaU, *acaV, voidIntegrator );

      if ( acaU->getNRows( ) * acaV->getNRows( ) <
          acaU->getNRows( ) * acaU->getNCols( ) +
          acaV->getNRows( ) * acaV->getNCols( ) ) {

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

    this->destroyIntegrator( voidIntegrator );
    ;
  }
  delete [] myAdmLeaves;
  delete [] myNonAdmLeaves;
}

}

#endif
