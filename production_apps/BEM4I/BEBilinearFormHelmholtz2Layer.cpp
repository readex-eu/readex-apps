/*!;
 * @file    BEBilinearFormHelmholtz2Layer.cpp
 * @author  Jan Zapletal
 * @date    August 14, 2013
 * 
 */

#ifdef BEBILINEARFORMHELMHOLTZ2LAYER_H

namespace bem4i {

template<class LO, class SC>
BEBilinearFormHelmholtz2Layer<LO, SC>::BEBilinearFormHelmholtz2Layer( ) {
}

template<class LO, class SC>
BEBilinearFormHelmholtz2Layer<LO, SC>::BEBilinearFormHelmholtz2Layer( const BEBilinearFormHelmholtz2Layer& orig ) {
  this->space = orig.space;
  this->quadratureOrder = orig.quadratureOrder;
}

template<class LO, class SC>
BEBilinearFormHelmholtz2Layer<LO, SC>::BEBilinearFormHelmholtz2Layer(
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

  this->quadratureOrderDisjointElems = quadratureOrderDisjointElems;
}

template<class LO, class SC>
BEBilinearFormHelmholtz2Layer<LO, SC>::~BEBilinearFormHelmholtz2Layer( ) {
}

template<class LO, class SC>
void BEBilinearFormHelmholtz2Layer<LO, SC>::assemble(
    BECluster<LO, SC> const *leftCluster,
    BECluster<LO, SC> const *rightCluster,
    FullMatrix<LO, SC>& matrix,
    void * voidIntegrator
    ) const {

  if ( this->space->getAnsatzFunctionType( ) == p1dis &&
      this->space->getTestFunctionType( ) == p0 ) {

    this->assembleP0P1Dis( leftCluster, rightCluster, matrix, voidIntegrator );
    return;
  } else if ( this->space->getAnsatzFunctionType( ) == p1dis &&
      this->space->getTestFunctionType( ) == p1dis ) {

    this->assembleP1DisP1Dis( leftCluster, rightCluster, matrix,
        voidIntegrator );
    return;
  } else if ( this->space->getAnsatzFunctionType( ) == p0 &&
      this->space->getTestFunctionType( ) == p0 ) {

    this->assembleP0P0( leftCluster, rightCluster, matrix, voidIntegrator );
    return;
  }

  BEIntegratorHelmholtz<LO, SC> * integrator =
      static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );

  const vector<LO>* innerElems = rightCluster->elems;
  const vector<LO>* outerElems = leftCluster->elems;

  // allocate local element matrices
  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  matrix.resize( this->space->getClusterOuterDOFs( leftCluster ),
      this->space->getClusterInnerDOFs( rightCluster ) );

  LO iMax = outerElems->size( );
  LO jMax = innerElems->size( );

  // vector of indices where to put values in the global matrix
  vector<LO> rowIndices( nLocalRows );
  vector<LO> colIndices( nLocalCols );
  FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );

  for ( LO j = 0; j < jMax; j++ ) {
    for ( LO i = 0; i < iMax; i++ ) {
      integrator->getElemMatrix2Layer( ( *outerElems )[i], ( *innerElems )[j],
          elemMatrix );
      this->space->getClusterOuterElemDOFs( leftCluster, i, &rowIndices[0] );
      this->space->getClusterInnerElemDOFs( rightCluster, j, &colIndices[0] );
      matrix.addToPositions( rowIndices, colIndices, elemMatrix );
    }
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtz2Layer<LO, SC>::assembleP0P0(
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

  matrix.resize( iMax, jMax );

  FullMatrix<LO, SC> elemMatrix( 1, 1 );

  for ( LO j = 0; j < jMax; ++j ) {
    for ( LO i = 0; i < iMax; ++i ) {
      integrator->getElemMatrix2Layer( outerElems->at( i ), innerElems->at( j ),
          elemMatrix );
      matrix.set( i, j, elemMatrix.get( 0, 0 ) );
    }
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtz2Layer<LO, SC>::assembleP0P1Dis(
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

  matrix.resize( iMax, 3 * jMax );

  FullMatrix<LO, SC> elemMatrix( 1, 3 );

  for ( LO j = 0; j < jMax; ++j ) {
    for ( LO i = 0; i < iMax; ++i ) {

      integrator->getElemMatrix2Layer( outerElems->at( i ), innerElems->at( j ),
          elemMatrix );
      matrix.set( i, 3 * j, elemMatrix.get( 0, 0 ) );
      matrix.set( i, 3 * j + 1, elemMatrix.get( 0, 1 ) );
      matrix.set( i, 3 * j + 2, elemMatrix.get( 0, 2 ) );
    }
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtz2Layer<LO, SC>::assembleP1DisP1Dis(
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

  for ( LO j = 0; j < jMax; ++j ) {
    for ( LO i = 0; i < iMax; ++i ) {

      integrator->getElemMatrix2Layer( outerElems->at( i ), innerElems->at( j ),
          elemMatrix );
      for ( int iRot = 0; iRot < 3; ++iRot ) {
        for ( int oRot = 0; oRot < 3; ++oRot ) {
          matrixData[ ( 3 * j + iRot ) * 3 * iMax + 3 * i + oRot ] =
              elemMatrixData[ iRot * 3 + oRot ];
        }
      }
    }
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtz2Layer<LO, SC>::assembleRow(
    const BEBlockCluster< LO, SC > & block,
    LO idx,
    const std::vector<LO> & cols,
    Vector<LO, SC> & row,
    void * voidIntegrator
    ) const {

#pragma omp atomic update
  ++this->counterRow;

  if ( this->space->getTestFunctionType( ) == p0 &&
      this->space->getAnsatzFunctionType( ) == p1dis ) {
    this->assembleRowP0P1Dis( block, idx, cols, row, voidIntegrator );
    return;
  } else if ( this->space->getTestFunctionType( ) == p0 &&
      this->space->getAnsatzFunctionType( ) == p0 ) {
    this->assembleRowP0P0( block, idx, cols, row, voidIntegrator );
    return;
  } else if ( this->space->getTestFunctionType( ) == p1dis &&
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
        integrator->getElemMatrix2Layer( testSupport[i], ansatzSupport[k],
            elemMatrix );

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
void BEBilinearFormHelmholtz2Layer<LO, SC>::assembleRowP0P0(
    const BEBlockCluster< LO, SC > & block,
    LO idx,
    const std::vector<LO> & cols,
    Vector<LO, SC> & row,
    void * voidIntegrator
    ) const {

  BEIntegratorHelmholtz<LO, SC> * integrator =
      static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );

  FullMatrix<LO, SC> elemMatrix( 1, 1 );

  for ( LO i = 0; i < block.rightCluster->nelems; ++i ) {
    integrator->getElemMatrix2Layer( idx, block.rightCluster->elems->at( i ),
        elemMatrix );
    row.set( i, elemMatrix.get( 0, 0 ) );
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtz2Layer<LO, SC>::assembleRowP0P1Dis(
    const BEBlockCluster< LO, SC > & block,
    LO idx,
    const std::vector<LO> & cols,
    Vector<LO, SC> & row,
    void * voidIntegrator
    ) const {

  BEIntegratorHelmholtz<LO, SC> * integrator =
      static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );

  FullMatrix<LO, SC> elemMatrix( 1, 3 );

  for ( LO i = 0; i < block.rightCluster->nelems; ++i ) {
    integrator->getElemMatrix2Layer( idx, block.rightCluster->elems->at( i ),
        elemMatrix );
    row.set( 3 * i, elemMatrix.get( 0, 0 ) );
    row.set( 3 * i + 1, elemMatrix.get( 0, 1 ) );
    row.set( 3 * i + 2, elemMatrix.get( 0, 2 ) );
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtz2Layer<LO, SC>::assembleRowP1DisP1Dis(
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
    integrator->getElemMatrix2Layer( idx / 3,
        block.rightCluster->elems->at( i ), elemMatrix );
    row.set( 3 * i, elemMatrix.get( idx % 3, 0 ) );
    row.set( 3 * i + 1, elemMatrix.get( idx % 3, 1 ) );
    row.set( 3 * i + 2, elemMatrix.get( idx % 3, 2 ) );
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtz2Layer<LO, SC>::assembleColumn(
    const BEBlockCluster< LO, SC > & block,
    LO idx,
    const std::vector< LO > & rows,
    Vector<LO, SC>& col,
    void * voidIntegrator
    ) const {

#pragma omp atomic update
  ++this->counterCol;

  if ( this->space->getTestFunctionType( ) == p0 &&
      this->space->getAnsatzFunctionType( ) == p1dis ) {
    this->assembleColumnP0P1Dis( block, idx, rows, col, voidIntegrator );
    return;
  } else if ( this->space->getTestFunctionType( ) == p0 &&
      this->space->getAnsatzFunctionType( ) == p0 ) {
    this->assembleColumnP0P0( block, idx, rows, col, voidIntegrator );
    return;
  } else if ( this->space->getTestFunctionType( ) == p1dis &&
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
  col.resize( nRows, true );
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
        integrator->getElemMatrix2Layer( testSupport[k], ansatzSupport[i],
            elemMatrix );

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
void BEBilinearFormHelmholtz2Layer<LO, SC>::assembleColumnP0P0(
    const BEBlockCluster< LO, SC > & block,
    LO idx,
    const std::vector<LO> & rows,
    Vector<LO, SC> & col,
    void * voidIntegrator
    ) const {

  BEIntegratorHelmholtz<LO, SC> * integrator =
      static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );

  FullMatrix<LO, SC> elemMatrix( 1, 1 );

  for ( LO i = 0; i < block.leftCluster->nelems; ++i ) {
    integrator->getElemMatrix2Layer( block.leftCluster->elems->at( i ), idx,
        elemMatrix );
    col.set( i, elemMatrix.get( 0, 0 ) );
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtz2Layer<LO, SC>::assembleColumnP0P1Dis(
    const BEBlockCluster< LO, SC > & block,
    LO idx,
    const std::vector<LO> & rows,
    Vector<LO, SC> & col,
    void * voidIntegrator
    ) const {

  BEIntegratorHelmholtz<LO, SC> * integrator =
      static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );

  FullMatrix<LO, SC> elemMatrix( 1, 3 );

  for ( LO i = 0; i < block.leftCluster->nelems; ++i ) {
    integrator->getElemMatrix2Layer( block.leftCluster->elems->at( i ), idx / 3,
        elemMatrix );
    col.set( i, elemMatrix.get( 0, idx % 3 ) );
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtz2Layer<LO, SC>::assembleColumnP1DisP1Dis(
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
    integrator->getElemMatrix2Layer(
        block.leftCluster->elems->at( i ), idx / 3, elemMatrix );
    col.set( 3 * i, elemMatrix.get( 0, idx % 3 ) );
    col.set( 3 * i + 1, elemMatrix.get( 1, idx % 3 ) );
    col.set( 3 * i + 2, elemMatrix.get( 2, idx % 3 ) );
  }
}

template<class LO, class SC>
void BEBilinearFormHelmholtz2Layer<LO, SC>::assemble(
    FullMatrix<LO, SC> & matrix
    ) const {

  //READEX_REGION_DEFINE(assemble_helmholtz_k);
  //READEX_REGION_START(assemble_helmholtz_k,"assemble_helmholtz_k",SCOREP_USER_REGION_TYPE_COMMON);
  
  if ( this->space->getAnsatzFunctionType( ) == p1 &&
      this->space->getTestFunctionType( ) == p0 ) {

    this->assembleP0P1( matrix );
    return;
  } else if ( this->space->getAnsatzFunctionType( ) == p1 &&
      this->space->getTestFunctionType( ) == p1 ) {

    this->assembleP1P1( matrix );
    return;
  } else if ( this->space->getAnsatzFunctionType( ) == p0 &&
      this->space->getTestFunctionType( ) == p0 ) {

    this->assembleP0P0( matrix );
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
  ProgressMonitor::init( "Assembling Helmholtz K", iMax * jMax );
#endif

#pragma omp parallel shared(matrix)
  {
    BEIntegratorHelmholtz<LO, SC> integrator( this->space, this->quadratureOrder,
        this->kappa, this->quadrature, this->quadratureOrderDisjointElems );

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
        integrator.getElemMatrix2Layer( outerElems[i], innerElems[j], *( matrixBuffer[counter] ) );
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
  
  //READEX_REGION_STOP(assemble_helmholtz_k); 
}

template<class LO, class SC>
void BEBilinearFormHelmholtz2Layer<LO, SC>::assembleP0P1(
    FullMatrix<LO, SC> & matrix
    ) const {

  //READEX_REGION_DEFINE(assemble_helmholtz_k01);
  //READEX_REGION_START(assemble_helmholtz_k01,"assemble_helmholtz_k01",SCOREP_USER_REGION_TYPE_COMMON);

  LO nNodes = this->space->getMesh( )->getNNodes( );
  LO nElems = this->space->getMesh( )->getNElements( );

  matrix.resize( nElems, nNodes, false );
  SC * matrixData = matrix.getData( );

#pragma omp parallel for schedule( static, 32 )
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nElems; ++i ) {
      matrixData[ j * nElems + i ] = 0.0;
    }
  }

#pragma omp parallel
  {
    BEIntegratorHelmholtz<LO, SC> integrator( this->space,
        this->quadratureOrder, this->kappa,
        this->quadrature, this->quadratureOrderDisjointElems );

    FullMatrix< LO, SC > elemMatrix( 1, 3 );
    LO elem[ 3 ];

    SC * elemMatrixData = elemMatrix.getData( );

#pragma omp for schedule( dynamic, 8 ) // no collapse, no atomic update
    for ( LO i = 0; i < nElems; ++i ) {
      for ( LO j = 0; j < nElems; ++j ) {
        integrator.getElemMatrix2Layer( i, j, elemMatrix );
        this->space->getMesh( )->getElement( j, elem );

        matrixData[ elem[ 0 ] * nElems + i ] += elemMatrixData[ 0 ];
        matrixData[ elem[ 1 ] * nElems + i ] += elemMatrixData[ 1 ];
        matrixData[ elem[ 2 ] * nElems + i ] += elemMatrixData[ 2 ];
      }
    }
  } // end omp parallel
  
  //READEX_REGION_STOP(assemble_helmholtz_k01); 
}

template<class LO, class SC>
void BEBilinearFormHelmholtz2Layer<LO, SC>::assembleP1P1(
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
        this->quadratureOrder, this->kappa,
        this->quadrature, this->quadratureOrderDisjointElems );

    FullMatrix< LO, SC > elemMatrix( 3, 3 );
    LO iElem[ 3 ], oElem[ 3 ];

    SC * elemMatrixData = elemMatrix.getData( );

#pragma omp for collapse( 2 ) schedule( dynamic, 32 )
    for ( LO i = 0; i < nElems; ++i ) { // i outer row
      for ( LO j = 0; j < nElems; ++j ) { // j inner col
        integrator.getElemMatrix2Layer( i, j, elemMatrix );
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

template<class LO, class SC>
void BEBilinearFormHelmholtz2Layer<LO, SC>::assembleP0P0(
    FullMatrix<LO, SC> & matrix
    ) const {

  LO nElems = this->space->getMesh( )->getNElements( );

  matrix.resize( nElems, nElems, false );
  SC * matrixData = matrix.getData( );

#pragma omp parallel
  {
    BEIntegratorHelmholtz<LO, SC> integrator( this->space,
        this->quadratureOrder, this->kappa, this->quadrature,
        this->quadratureOrderDisjointElems );

    FullMatrix< LO, SC > elemMatrix( 1, 1 );

#pragma omp for schedule( dynamic, 8 )
    for ( LO j = 0; j < nElems; ++j ) {
      for ( LO i = 0; i < nElems; ++i ) {
        integrator.getElemMatrix2Layer( i, j, elemMatrix );
        matrixData[ j * nElems + i ] = elemMatrix.get( 0, 0 );
      }
    }
  } // end omp parallel
}

}

#endif
