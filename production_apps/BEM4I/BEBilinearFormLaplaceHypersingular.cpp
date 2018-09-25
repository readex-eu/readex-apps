/*!
 * @file    BEBilinearFormLaplaceHypersingular.cpp
 * @author  Jan Zapletal
 * @date    October 17, 2013
 * 
 */

#ifdef BEBILINEARFORMLAPLACEHYPERSINGULAR_H

namespace bem4i {

template<class LO, class SC>
BEBilinearFormLaplaceHypersingular<LO, SC>::
BEBilinearFormLaplaceHypersingular( ) {
  this->initCurl( );
}

template<class LO, class SC>
BEBilinearFormLaplaceHypersingular<LO, SC>::BEBilinearFormLaplaceHypersingular(
    const BEBilinearFormLaplaceHypersingular & orig
    ) {
  this->space = orig.space;
  memcpy( &( this->x1[0] ), &( orig.x1[0] ), 3 * sizeof (SC ) );
  memcpy( &( this->x2[0] ), &( orig.x2[0] ), 3 * sizeof (SC ) );
  memcpy( &( this->x3[0] ), &( orig.x3[0] ), 3 * sizeof (SC ) );
}

template<class LO, class SC>
BEBilinearFormLaplaceHypersingular<LO, SC>::BEBilinearFormLaplaceHypersingular(
    BESpace<LO, SC>* space,
    int* quadratureOrder,
    quadratureType quadrature,
    int* quadratureOrderDisjointElems
    ) {
  this->quadrature = quadrature;
  this->space = space;
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


  initCurl( );
}

template<class LO, class SC>
BEBilinearFormLaplaceHypersingular<LO, SC>::
~BEBilinearFormLaplaceHypersingular( ) {

}

template<class LO, class SC>
void BEBilinearFormLaplaceHypersingular<LO, SC>::assemble(
    LaplaceHypersingularOperator<LO, SC> & op
    ) const {

#ifdef VERBOSE
  ProgressMonitor::init( "Assembling Laplace D" );
#endif

  op.assemble( );

#ifdef VERBOSE
  ProgressMonitor::step( );
#endif

}

template<class LO, class SC>
void BEBilinearFormLaplaceHypersingular<LO, SC>::assemble(
    FullMatrix<LO, SC>& matrix
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
  ProgressMonitor::init( "Assembling Laplace D", iMax * jMax );
#endif

#pragma omp parallel shared(matrix)
  {
    BEIntegratorLaplace<LO, SC> integrator( this->space, this->quadratureOrder,
        this->quadrature, this->quadratureOrderDisjointElems );

    // bufferd vector of indices where to put values in the global matrix
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
void BEBilinearFormLaplaceHypersingular<LO, SC>::assembleP1P1(
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
    BEIntegratorLaplace<LO, SC> integrator( this->space,
        this->quadratureOrder, this->quadrature,
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
            matrixData[ iElem[ iRot ] * nNodes + oElem[ oRot ] ] +=
                elemMatrixData[ iRot * 3 + oRot ];
          }
        }
      }
    }
  } // end omp parallel
}
/*
template<class LO, class SC>
void BEBilinearFormLaplaceHypersingular<LO, SC>::assemble(
    FullMatrix<LO, SC> & matrix,
    FullMatrix<LO, SC> & V
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
  ProgressMonitor::init( "Assembling Laplace D", iMax * jMax );
#endif

#pragma omp parallel shared(matrix)
  {
    BEIntegratorLaplace<LO, SC> integrator( this->space, this->quadratureOrder,
        this->quadrature );

    // bufferd vector of indices where to put values in the global matrix
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
    for ( LO j = 0; j < jMax; j++ ) {
      for ( LO i = 0; i < iMax; i++ ) {
        integrator.getElemMatrixHypersingular( outerElems[i], innerElems[j],
            V, *( matrixBuffer[counter] ) );
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
*/
template<class LO, class SC>
void BEBilinearFormLaplaceHypersingular<LO, SC>::assemble(
    const BECluster< LO, SC > * leftCluster,
    const BECluster< LO, SC > * rightCluster,
    FullMatrix< LO, SC > & matrix,
    void * voidIntegrator
    ) const {

  if ( this->space->getAnsatzFunctionType( ) == p1dis &&
      this->space->getTestFunctionType( ) == p1dis ) {

    this->assembleP1DisP1Dis( leftCluster, rightCluster, matrix,
        voidIntegrator );
    return;
  }

  vector<LO> * innerElems = rightCluster->elems;
  vector<LO> * outerElems = leftCluster->elems;

  // allocate local element matrices
  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  matrix.resize( this->space->getClusterOuterDOFs( leftCluster ),
      this->space->getClusterInnerDOFs( rightCluster ) );

  LO iMax = outerElems->size( );
  LO jMax = innerElems->size( );

  BEIntegratorLaplace<LO, SC> * integrator =
      static_cast<BEIntegratorLaplace< LO, SC > *> ( voidIntegrator );

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
void BEBilinearFormLaplaceHypersingular<LO, SC>::assembleP1DisP1Dis(
    const BECluster<LO, SC> * leftCluster,
    const BECluster<LO, SC> * rightCluster,
    FullMatrix<LO, SC> & matrix,
    void * voidIntegrator
    ) const {

  BEIntegratorLaplace<LO, SC> * integrator =
      static_cast<BEIntegratorLaplace< LO, SC > *> ( voidIntegrator );

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
void BEBilinearFormLaplaceHypersingular<LO, SC>::assembleRow(
    const BEBlockCluster<LO, SC>& block,
    LO idx,
    const std::vector<LO>& cols,
    Vector<LO, SC>& row,
    void* voidIntegrator
    ) const {

#pragma omp atomic update
  ++this->counterRow;

  if ( this->space->getTestFunctionType( ) == p1dis &&
      this->space->getAnsatzFunctionType( ) == p1dis ) {
    this->assembleRowP1DisP1Dis( block, idx, cols, row, voidIntegrator );
    return;
  }

  BEIntegratorLaplace<LO, SC> * integrator =
      static_cast<BEIntegratorLaplace< LO, SC > *> ( voidIntegrator );

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
void BEBilinearFormLaplaceHypersingular<LO, SC>::assembleRowP1DisP1Dis(
    const BEBlockCluster<LO, SC> & block,
    LO idx,
    const std::vector<LO> & cols,
    Vector<LO, SC> & row,
    void * voidIntegrator
    ) const {

  BEIntegratorLaplace<LO, SC> * integrator =
      static_cast<BEIntegratorLaplace< LO, SC > *> ( voidIntegrator );

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
void BEBilinearFormLaplaceHypersingular<LO, SC>::assembleColumn(
    const BEBlockCluster<LO, SC> & block,
    LO idx,
    const std::vector<LO> & rows,
    Vector<LO, SC> & col,
    void * voidIntegrator
    ) const {

#pragma omp atomic update
  ++this->counterCol;

  if ( this->space->getTestFunctionType( ) == p1dis &&
      this->space->getAnsatzFunctionType( ) == p1dis ) {
    this->assembleColumnP1DisP1Dis( block, idx, rows, col, voidIntegrator );
    return;
  }

  BEIntegratorLaplace<LO, SC> * integrator =
      static_cast<BEIntegratorLaplace< LO, SC > *> ( voidIntegrator );

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
void BEBilinearFormLaplaceHypersingular<LO, SC>::assembleColumnP1DisP1Dis(
    const BEBlockCluster<LO, SC>& block,
    LO idx,
    const std::vector<LO>& rows,
    Vector<LO, SC>& col,
    void * voidIntegrator
    ) const {

  BEIntegratorLaplace<LO, SC> * integrator =
      static_cast<BEIntegratorLaplace< LO, SC > *> ( voidIntegrator );

  FullMatrix<LO, SC> elemMatrix( 3, 3 );

  for ( LO i = 0; i < block.leftCluster->nelems; ++i ) {
    integrator->getElemMatrixHypersingular(
        block.leftCluster->elems->at( i ), idx / 3, elemMatrix );
    col.set( 3 * i, elemMatrix.get( 0, idx % 3 ) );
    col.set( 3 * i + 1, elemMatrix.get( 1, idx % 3 ) );
    col.set( 3 * i + 2, elemMatrix.get( 2, idx % 3 ) );
  }
}

}

#endif
