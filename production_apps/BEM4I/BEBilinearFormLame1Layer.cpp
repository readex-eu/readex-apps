/*!
 * @file    BEBilinearFormLame1Layer.cpp
 * @author  Michal Merta
 * @author  Jan Zapletal
 * @date    April 29, 2015
 *
 */

#ifdef BEBILINEARFORMLAME1LAYER_H

namespace bem4i {

template<class LO, class SC>
BEBilinearFormLame1Layer<LO, SC>::BEBilinearFormLame1Layer( ) {
}

template<class LO, class SC>
BEBilinearFormLame1Layer<LO, SC>::BEBilinearFormLame1Layer(
    const BEBilinearFormLame1Layer& orig
    ) {
}

template<class LO, class SC>
BEBilinearFormLame1Layer<LO, SC>::BEBilinearFormLame1Layer(
    BESpace<LO, SC>* space,
    int* quadratureOrder,
    quadratureType quadrature,
    int* quadratureOrderDisjointElems,
    bool employSymmetricity
    ) {
  this->space = space;
  this->quadrature = quadrature;

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

  // default values = steel
  this->nu = 0.29;
  this->E = 2.0e11;

  this->quadratureOrderDisjointElems = quadratureOrderDisjointElems;
  this->employSymmetricity = employSymmetricity;
}

template<class LO, class SC>
BEBilinearFormLame1Layer<LO, SC>::~BEBilinearFormLame1Layer( ) {
}

/*
template<class LO, class SC>
void BEBilinearFormLame1Layer<LO, SC>::assemble(
    FullMatrix<LO, SC>& matrix,
    FullMatrix<LO, SC>& V
    ) const {

  const vector<LO> innerElems = this->space->getInnerElems( );
  const vector<LO> outerElems = this->space->getOuterElems( );

  // allocate local element matrices
  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  LO dim = 3;

  LO nRowsScalar = this->space->getOuterDOFs( );
  LO nColsScalar = this->space->getInnerDOFs( );

  matrix.resize( dim * nRowsScalar, dim * nColsScalar );
  V.resize( nRowsScalar, nColsScalar );

  LO iMax = outerElems.size( );
  LO jMax = innerElems.size( );

  bool symmetry = this->space->getTestFunctionType( ) == p0 &&
      this->space->getAnsatzFunctionType( ) == p0 &&
      this->employSymmetricity;

  LO sizeVector = 0;
  if ( symmetry ) {
    sizeVector = 18;
  } else {
    sizeVector = 9;
  }
  // vector of global row indices in matrix
  std::vector<LO> glRInd( sizeVector );
  // vector of global column indices in matrix
  std::vector<LO> glCInd( sizeVector );
  // vector of global values in matrix
  std::vector<SC> glVal( sizeVector );

#ifdef VERBOSE
  ProgressMonitor::init( "Assembling Laplace V", iMax * jMax );
#endif

#pragma omp parallel shared(matrix)
  {
    BEIntegratorLame<LO, SC> integrator( this->space, this->quadratureOrder,
        this->quadrature, this->quadratureOrderDisjointElems );
    // vector of indices where to put values in the global matrix
    FullMatrix<LO, SC>* matrixBuffer[BUFFER_SIZE];
    vector<LO>* rowIndicesBuffer[BUFFER_SIZE];
    vector<LO>* colIndicesBuffer[BUFFER_SIZE];
    for ( int i = 0; i < BUFFER_SIZE; i++ ) {
      matrixBuffer[i] = new FullMatrix<LO, SC>( nLocalRows, 7 * nLocalCols );
      rowIndicesBuffer[i] = new vector<LO>( nLocalRows );
      colIndicesBuffer[i] = new vector<LO>( nLocalCols );
    }
    int counter = 0;

    if ( symmetry ) {
#pragma omp for schedule(dynamic)
      for ( LO i = 0; i < iMax; i++ ) {
        for ( LO j = 0; j <= i; j++ ) {

          integrator.getElemMatrix1Layer( outerElems[i], innerElems[j],
 *( matrixBuffer[counter] ) );

          this->space->getOuterElemDOFs( outerElems[i],
              &( *( rowIndicesBuffer[counter] ) )[0] );
          this->space->getInnerElemDOFs( innerElems[j],
              &( *( colIndicesBuffer[counter] ) )[0] );

          counter++;

          if ( counter == BUFFER_SIZE ) {
            counter = 0;
#pragma omp critical
            {
              for ( int k = 0; k < BUFFER_SIZE; k++ ) {
                this->getGlobalIndices( rowIndicesBuffer[k]->at( 0 ),
                    colIndicesBuffer[k]->at( 0 ), *( matrixBuffer[k] ),
                    glRInd, glCInd, glVal );

                matrix.addToPositions( glRInd, glCInd, glVal );
                V.add( rowIndicesBuffer[k]->at( 0 ),
                    colIndicesBuffer[k]->at( 0 ),
                    matrixBuffer[k]->get( 0, 0 ) );
                V.add( colIndicesBuffer[k]->at( 0 ),
                    rowIndicesBuffer[k]->at( 0 ),
                    matrixBuffer[k]->get( 0, 0 ) );
              }
            }
#ifdef VERBOSE
            ProgressMonitor::step( BUFFER_SIZE );
#endif
          }
        }
      }
#pragma omp critical
      {
        for ( int i = 0; i < counter; i++ ) {
          //( matrixBuffer[i] )->print();
          this->getGlobalIndices( rowIndicesBuffer[i]->at( 0 ),
              colIndicesBuffer[i]->at( 0 ), *( matrixBuffer[i] ),
              glRInd, glCInd, glVal );

          matrix.addToPositions( glRInd, glCInd, glVal );
          V.add( rowIndicesBuffer[i]->at( 0 ), colIndicesBuffer[i]->at( 0 ),
              matrixBuffer[i]->get( 0, 0 ) );
          V.add( colIndicesBuffer[i]->at( 0 ), rowIndicesBuffer[i]->at( 0 ),
              matrixBuffer[i]->get( 0, 0 ) );
        }
      }
    } else {
#pragma omp for collapse(2)
      for ( LO i = 0; i < iMax; i++ ) {
        for ( LO j = 0; j < jMax; j++ ) {

          integrator.getElemMatrix1Layer( outerElems[i], innerElems[j],
 *( matrixBuffer[counter] ) );

          this->space->getOuterElemDOFs( outerElems[i],
              &( *( rowIndicesBuffer[counter] ) )[0] );
          this->space->getInnerElemDOFs( innerElems[j],
              &( *( colIndicesBuffer[counter] ) )[0] );

          counter++;

          if ( counter == BUFFER_SIZE ) {
            counter = 0;
#pragma omp critical
            {
              for ( int k = 0; k < BUFFER_SIZE; k++ ) {

                this->getGlobalIndices( rowIndicesBuffer[k]->at( 0 ),
                    colIndicesBuffer[k]->at( 0 ), *( matrixBuffer[k] ),
                    glRInd, glCInd, glVal );

                matrix.addToPositions( glRInd, glCInd, glVal );
                V.add( rowIndicesBuffer[k]->at( 0 ), colIndicesBuffer[k]->at( 0 ),
                    matrixBuffer[k]->get( 0, 0 ) );
              }
            }
          }
#ifdef VERBOSE
          ProgressMonitor::step( BUFFER_SIZE );
#endif
        }
      }
#pragma omp critical
      {
        for ( int i = 0; i < counter; i++ ) {
          this->getGlobalIndices( rowIndicesBuffer[i]->at( 0 ),
              colIndicesBuffer[i]->at( 0 ), *( matrixBuffer[i] ),
              glRInd, glCInd, glVal );

          matrix.addToPositions( glRInd, glCInd, glVal );
          V.add( rowIndicesBuffer[i]->at( 0 ), colIndicesBuffer[i]->at( 0 ),
              matrixBuffer[i]->get( 0, 0 ) );
        }
      }
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

  if ( symmetry ) {
    for ( LO i = 0; i < nRowsScalar; i++ ) {
      V.set( i, i, V.get( i, i ) / 2.0 );
    }
    for ( LO i = 0; i < matrix.getNRows( ); i++ ) {
      matrix.set( i, i, matrix.get( i, i ) / 2.0 );

    }
    for ( LO i = 0; i < matrix.getNRows( ) - nRowsScalar; i++ ) {
      matrix.set( i + nRowsScalar, i, matrix.get( i + nRowsScalar, i ) / 2.0 );
      matrix.set( i, i + nRowsScalar, matrix.get( i, i + nRowsScalar ) / 2.0 );
    }
    for ( LO i = 0; i < matrix.getNRows( ) - 2 * nRowsScalar; i++ ) {
      matrix.set( i + 2 * nRowsScalar, i, matrix.get( i + 2 * nRowsScalar, i )
          / 2.0 );
      matrix.set( i, i + 2 * nRowsScalar, matrix.get( i, i + 2 * nRowsScalar )
          / 2.0 );
    }
  }

  // add the Laplace matrix V to the diagonal of the result
  SC val = 0.0;
#pragma omp parallel for shared(matrix) private(val)
  for ( LO j = 0; j < nColsScalar; ++j ) {
    for ( LO i = 0; i < nRowsScalar; ++i ) {
      val = V.get( i, j ) * ( 3.0 - 4.0 * this->nu );
      matrix.add( i, j, val );
      matrix.add( i + nRowsScalar, j + nColsScalar, val );
      matrix.add( i + 2 * nRowsScalar, j + 2 * nColsScalar, val );
    }
  }

  matrix.scale( ( 1.0 + this->nu ) / ( 2.0 * E * ( 1.0 - this->nu ) ) );

  //matrix.print( );
}
 */
///*

template<class LO, class SC>
void BEBilinearFormLame1Layer<LO, SC>::assemble(
    FullMatrix<LO, SC> & matrix,
    FullMatrix<LO, SC> & V
    ) const {

  if ( this->space->getTestFunctionType( ) == p0 &&
      this->space->getAnsatzFunctionType( ) == p0 ) {

    this->assembleP0P0( matrix, V );
    return;
  } else if ( this->space->getTestFunctionType( ) == p1 &&
      this->space->getAnsatzFunctionType( ) == p1 ) {

    this->assembleP1P1( matrix, V );
    return;
  }

}

template<class LO, class SC>
void BEBilinearFormLame1Layer<LO, SC>::assembleP0P0(
    FullMatrix<LO, SC> & matrix,
    FullMatrix<LO, SC> & V
    ) const {

  LO nElems = this->space->getMesh( )->getNElements( );
  LO nElems3 = 3 * nElems;

  matrix.resize( nElems3, nElems3, false );
  V.resize( nElems, nElems, false );

  SC * matrixData = matrix.getData( );
  SC * VData = V.getData( );

#pragma omp parallel
  {
    BEIntegratorLame<LO, SC> integrator( this->space, this->quadratureOrder,
        this->quadrature, this->quadratureOrderDisjointElems );

    FullMatrix< LO, SC > elemMatrix( 1, 7 );
    SC * elemMatrixData = elemMatrix.getData( );

    if ( this->employSymmetricity ) {
#pragma omp for schedule( dynamic, 4 )
      for ( LO i = 0; i < nElems; ++i ) {
        for ( LO j = 0; j <= i; ++j ) {

          integrator.getElemMatrix1Layer( i, j, elemMatrix );

          matrixData[ j * nElems3 + i ] =
              elemMatrixData[ 1 ];
          matrixData[ j * nElems3 + i + nElems ] =
              elemMatrixData[ 4 ];
          matrixData[ j * nElems3 + i + 2 * nElems ] =
              elemMatrixData[ 5 ];
          matrixData[ ( j + nElems ) * nElems3 + i ] =
              elemMatrixData[ 4 ];
          matrixData[ ( j + nElems ) * nElems3 + i + nElems ] =
              elemMatrixData[ 2 ];
          matrixData[ ( j + nElems ) * nElems3 + i + 2 * nElems ] =
              elemMatrixData[ 6 ];
          matrixData[ ( j + 2 * nElems ) * nElems3 + i ] =
              elemMatrixData[ 5 ];
          matrixData[ ( j + 2 * nElems ) * nElems3 + i + nElems ] =
              elemMatrixData[ 6 ];
          matrixData[ ( j + 2 * nElems ) * nElems3 + i + 2 * nElems ] =
              elemMatrixData[ 3 ];

          if ( i != j ) {
            matrixData[ i * nElems3 + j ] =
                elemMatrixData[ 1 ];
            matrixData[ ( i + nElems ) * nElems3 + j ] =
                elemMatrixData[ 4 ];
            matrixData[ ( i + 2 * nElems ) * nElems3 + j ] =
                elemMatrixData[ 5 ];
            matrixData[ i * nElems3 + j + nElems ] =
                elemMatrixData[ 4 ];
            matrixData[ ( i + nElems ) * nElems3 + j + nElems ] =
                elemMatrixData[ 2 ];
            matrixData[ ( i + 2 * nElems ) * nElems3 + j + nElems ] =
                elemMatrixData[ 6 ];
            matrixData[ i * nElems3 + j + 2 * nElems ] =
                elemMatrixData[ 5 ];
            matrixData[ ( i + nElems ) * nElems3 + j + 2 * nElems ] =
                elemMatrixData[ 6 ];
            matrixData[ ( i + 2 * nElems ) * nElems3 + j + 2 * nElems ] =
                elemMatrixData[ 3 ];
          }

          VData[ j * nElems + i ] = elemMatrixData[ 0 ];
          if ( i != j ) {
            VData[ i * nElems + j ] = elemMatrixData[ 0 ];
          }
        }
      }
    } else {
#pragma omp for schedule( dynamic, 8 )
      for ( LO j = 0; j < nElems; ++j ) {
        for ( LO i = 0; i < nElems; ++i ) {

          integrator.getElemMatrix1Layer( i, j, elemMatrix );

          matrixData[ j * nElems3 + i ] =
              elemMatrixData[ 1 ];
          matrixData[ j * nElems3 + i + nElems ] =
              elemMatrixData[ 4 ];
          matrixData[ j * nElems3 + i + 2 * nElems ] =
              elemMatrixData[ 5 ];
          matrixData[ ( j + nElems ) * nElems3 + i ] =
              elemMatrixData[ 4 ];
          matrixData[ ( j + nElems ) * nElems3 + i + nElems ] =
              elemMatrixData[ 2 ];
          matrixData[ ( j + nElems ) * nElems3 + i + 2 * nElems ] =
              elemMatrixData[ 6 ];
          matrixData[ ( j + 2 * nElems ) * nElems3 + i ] =
              elemMatrixData[ 5 ];
          matrixData[ ( j + 2 * nElems ) * nElems3 + i + nElems ] =
              elemMatrixData[ 6 ];
          matrixData[ ( j + 2 * nElems ) * nElems3 + i + 2 * nElems ] =
              elemMatrixData[ 3 ];

          VData[ j * nElems + i ] = elemMatrixData[ 0 ];
        }
      }
    }
  }

  // add the Laplace matrix V to the diagonal of the result
  SC val = 0.0;
#pragma omp parallel for schedule( dynamic, 32 ) private( val )
  for ( LO j = 0; j < nElems; ++j ) {
    for ( LO i = 0; i < nElems; ++i ) {
      val = VData[ j * nElems + i ] * ( 3.0 - 4.0 * this->nu );
      matrixData[ j * nElems3 + i ] += val;
      matrixData[ ( j + nElems ) * nElems3 + i + nElems ] += val;
      matrixData[ ( j + 2 * nElems ) * nElems3 + i + 2 * nElems ] += val;
    }
  }

  SCVT mult = ( (SCVT) 1.0 + this->nu ) /
      ( (SCVT) 2.0 * E * ( (SCVT) 1.0 - this->nu ) );

#pragma omp parallel for simd
  for ( LO i = 0; i < nElems3 * nElems3; ++i ) {
    matrixData[ i ] *= mult;
  }
}

template<class LO, class SC>
void BEBilinearFormLame1Layer<LO, SC>::assembleP1P1(
    FullMatrix<LO, SC> & matrix,
    FullMatrix<LO, SC> & V
    ) const {

  LO nNodes = this->space->getMesh( )->getNNodes( );
  LO nElems = this->space->getMesh( )->getNElements( );
  LO nNodes3 = 3 * nNodes;

  matrix.resize( nNodes3, nNodes3, false );
  V.resize( nNodes, nNodes, false );

  SC * matrixData = matrix.getData( );
  SC * VData = V.getData( );

#pragma omp parallel for schedule( static, 32 )
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nNodes; ++i ) {
      VData[ j * nNodes + i ] = 0.0;
    }
  }

#pragma omp parallel for schedule( static, 32 )
  for ( LO j = 0; j < nNodes3; ++j ) {
    for ( LO i = 0; i < nNodes3; ++i ) {
      matrixData[ j * nNodes3 + i ] = 0.0;
    }
  }

#pragma omp parallel
  {
    BEIntegratorLame<LO, SC> integrator( this->space, this->quadratureOrder,
        this->quadrature, this->quadratureOrderDisjointElems );

    FullMatrix< LO, SC > elemMatrix( 3, 21 );
    SC * elemMatrixData = elemMatrix.getData( );
    LO iElem[ 3 ], oElem[ 3 ];

#pragma omp for schedule( dynamic, 8 )
    for ( LO j = 0; j < nElems; ++j ) {
      for ( LO i = 0; i < nElems; ++i ) {

        integrator.getElemMatrix1Layer( i, j, elemMatrix );

        this->space->getMesh( )->getElement( i, oElem );
        this->space->getMesh( )->getElement( j, iElem );

        for ( int oRot = 0; oRot < 3; ++oRot ) {
          for ( int iRot = 0; iRot < 3; ++iRot ) {
            //Vlap
#pragma omp atomic update
            VData[ iElem[ iRot ] * nNodes + oElem[ oRot ] ] +=
                elemMatrixData[ iRot * 3 + oRot ];
            //V11
#pragma omp atomic update
            matrixData[ iElem[ iRot ] * nNodes3 + oElem[ oRot ] ] +=
                elemMatrixData[ iRot * 3 + oRot + 9 ];
            //V22
#pragma omp atomic update
            matrixData[ ( iElem[ iRot ] + nNodes ) * nNodes3 +
                oElem[ oRot ] + nNodes ] +=
                elemMatrixData[ iRot * 3 + oRot + 18 ];
            //V33
#pragma omp atomic update
            matrixData[ ( iElem[ iRot ] + 2 * nNodes ) * nNodes3 +
                oElem[ oRot ] + 2 * nNodes ] +=
                elemMatrixData[ iRot * 3 + oRot + 27 ];
            //V12
#pragma omp atomic update
            matrixData[ ( iElem[ iRot ] + nNodes ) * nNodes3 +
                oElem[ oRot ] ] +=
                elemMatrixData[ iRot * 3 + oRot + 36 ];
#pragma omp atomic update
            matrixData[ iElem[ iRot ] * nNodes3 +
                oElem[ oRot ] + nNodes ] +=
                elemMatrixData[ iRot * 3 + oRot + 36 ];
            //V13
#pragma omp atomic update
            matrixData[ ( iElem[ iRot ] + 2 * nNodes ) * nNodes3 +
                oElem[ oRot ] ] +=
                elemMatrixData[ iRot * 3 + oRot + 45 ];
#pragma omp atomic update
            matrixData[ iElem[ iRot ] * nNodes3 +
                oElem[ oRot ] + 2 * nNodes ] +=
                elemMatrixData[ iRot * 3 + oRot + 45 ];
            //V23
#pragma omp atomic update
            matrixData[ ( iElem[ iRot ] + 2 * nNodes ) * nNodes3 +
                oElem[ oRot ] + nNodes ] +=
                elemMatrixData[ iRot * 3 + oRot + 54 ];
#pragma omp atomic update
            matrixData[ ( iElem[ iRot ] + nNodes ) * nNodes3 +
                oElem[ oRot ] + 2 * nNodes ] +=
                elemMatrixData[ iRot * 3 + oRot + 54 ];
          }
        }
      }
    }
  }

  // add the Laplace matrix V to the diagonal of the result
  SC val = 0.0;
#pragma omp parallel for schedule( dynamic, 32 ) private( val )
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nNodes; ++i ) {
      val = VData[ j * nNodes + i ] * ( 3.0 - 4.0 * this->nu );
      matrixData[ j * nNodes3 + i ] += val;
      matrixData[ ( j + nNodes ) * nNodes3 + i + nNodes ] += val;
      matrixData[ ( j + 2 * nNodes ) * nNodes3 + i + 2 * nNodes ] += val;
    }
  }

  SCVT mult = ( (SCVT) 1.0 + this->nu ) /
      ( (SCVT) 2.0 * E * ( (SCVT) 1.0 - this->nu ) );

#pragma omp parallel for simd
  for ( LO i = 0; i < nNodes3 * nNodes3; ++i ) {
    matrixData[ i ] *= mult;
  }
}

template<class LO, class SC>
void BEBilinearFormLame1Layer<LO, SC>::getGlobalIndices(
    LO rowInd,
    LO colInd,
    const FullMatrix<LO, SC> & values,
    std::vector<LO> & glRows,
    std::vector<LO> & glCols,
    std::vector<SC> & glValues
    ) const {

  LO nRowsScalar = this->space->getOuterDOFs( );
  LO nColsScalar = this->space->getInnerDOFs( );

  if ( this->space->getTestFunctionType( ) == p0 &&
      this->space->getAnsatzFunctionType( ) == p0 &&
      this->employSymmetricity ) {

    // V11 upper + lower
    glRows[0] = rowInd;
    glCols[0] = colInd;
    glRows[1] = colInd;
    glCols[1] = rowInd;
    // V21 upper + lower
    glRows[2] = rowInd + nRowsScalar;
    glCols[2] = colInd;
    glRows[3] = colInd + nRowsScalar;
    glCols[3] = rowInd;
    // V31 upper + lower
    glRows[4] = rowInd + 2 * nRowsScalar;
    glCols[4] = colInd;
    glRows[5] = colInd + 2 * nRowsScalar;
    glCols[5] = rowInd;
    // V12 upper + lower
    glRows[6] = rowInd;
    glCols[6] = colInd + nColsScalar;
    glRows[7] = colInd;
    glCols[7] = rowInd + nColsScalar;
    // V22 upper + lower
    glRows[8] = rowInd + nRowsScalar;
    glCols[8] = colInd + nColsScalar;
    glRows[9] = colInd + nRowsScalar;
    glCols[9] = rowInd + nColsScalar;
    // V32 upper + lower
    glRows[10] = rowInd + 2 * nRowsScalar;
    glCols[10] = colInd + nColsScalar;
    glRows[11] = colInd + 2 * nRowsScalar;
    glCols[11] = rowInd + nColsScalar;
    // V13 upper + lower
    glRows[12] = rowInd;
    glCols[12] = colInd + 2 * nColsScalar;
    glRows[13] = colInd;
    glCols[13] = rowInd + 2 * nColsScalar;
    // V32 upper + lower
    glRows[14] = rowInd + nRowsScalar;
    glCols[14] = colInd + 2 * nColsScalar;
    glRows[15] = colInd + nRowsScalar;
    glCols[15] = rowInd + 2 * nColsScalar;
    // V33 upper + lower
    glRows[16] = rowInd + 2 * nRowsScalar;
    glCols[16] = colInd + 2 * nColsScalar;
    glRows[17] = colInd + 2 * nRowsScalar;
    glCols[17] = rowInd + 2 * nColsScalar;

    // values are stored in matrix in order Vlap, V11, V22, V33, V12, V13, V23
    // V11
    glValues[0] = values.get( 0, 1 );
    glValues[1] = values.get( 0, 1 );
    // V21
    glValues[2] = values.get( 0, 4 );
    glValues[3] = values.get( 0, 4 );
    // V31
    glValues[4] = values.get( 0, 5 );
    glValues[5] = values.get( 0, 5 );
    // V12
    glValues[6] = values.get( 0, 4 );
    glValues[7] = values.get( 0, 4 );
    // V22
    glValues[8] = values.get( 0, 2 );
    glValues[9] = values.get( 0, 2 );
    // V32
    glValues[10] = values.get( 0, 6 );
    glValues[11] = values.get( 0, 6 );
    // V13
    glValues[12] = values.get( 0, 5 );
    glValues[13] = values.get( 0, 5 );
    // V23
    glValues[14] = values.get( 0, 6 );
    glValues[15] = values.get( 0, 6 );
    // V33
    glValues[16] = values.get( 0, 3 );
    glValues[17] = values.get( 0, 3 );
  } else {
    // V11
    glRows[0] = rowInd;
    glCols[0] = colInd;
    // V21
    glRows[1] = rowInd + nRowsScalar;
    glCols[1] = colInd;
    // V31
    glRows[2] = rowInd + 2 * nRowsScalar;
    glCols[2] = colInd;
    // V12
    glRows[3] = rowInd;
    glCols[3] = colInd + nColsScalar;
    // V22
    glRows[4] = rowInd + nRowsScalar;
    glCols[4] = colInd + nColsScalar;
    // V32
    glRows[5] = rowInd + 2 * nRowsScalar;
    glCols[5] = colInd + nColsScalar;
    // V13
    glRows[6] = rowInd;
    glCols[6] = colInd + 2 * nColsScalar;
    // V32
    glRows[7] = rowInd + nRowsScalar;
    glCols[7] = colInd + 2 * nColsScalar;
    // V33
    glRows[8] = rowInd + 2 * nRowsScalar;
    glCols[8] = colInd + 2 * nColsScalar;

    // values are stored in matrix in order Vlap, V11, V22, V33, V12, V13, V23
    // V11
    glValues[0] = values.get( 0, 1 );
    // V21
    glValues[1] = values.get( 0, 4 );
    // V31
    glValues[2] = values.get( 0, 5 );
    // V12
    glValues[3] = values.get( 0, 4 );
    // V22
    glValues[4] = values.get( 0, 2 );
    // V32
    glValues[5] = values.get( 0, 6 );
    // V13
    glValues[6] = values.get( 0, 5 );
    // V23
    glValues[7] = values.get( 0, 6 );
    // V33
    glValues[8] = values.get( 0, 3 );
  }
}

template<class LO, class SC>
void BEBilinearFormLame1Layer<LO, SC>::getGlobalIndices(
    LO rowInd,
    LO colInd,
    const FullMatrix<LO, SC> & values,
    LO * glRows,
    LO * glCols,
    SC * glValues
    ) const {

  LO nElems = this->space->getMesh( )->getNElements( );

  // V11
  glRows[0] = rowInd;
  glCols[0] = colInd;
  // V21
  glRows[1] = rowInd + nElems;
  glCols[1] = colInd;
  // V31
  glRows[2] = rowInd + 2 * nElems;
  glCols[2] = colInd;
  // V12
  glRows[3] = rowInd;
  glCols[3] = colInd + nElems;
  // V22
  glRows[4] = rowInd + nElems;
  glCols[4] = colInd + nElems;
  // V32
  glRows[5] = rowInd + 2 * nElems;
  glCols[5] = colInd + nElems;
  // V13
  glRows[6] = rowInd;
  glCols[6] = colInd + 2 * nElems;
  // V32
  glRows[7] = rowInd + nElems;
  glCols[7] = colInd + 2 * nElems;
  // V33
  glRows[8] = rowInd + 2 * nElems;
  glCols[8] = colInd + 2 * nElems;

  // values are stored in matrix in order Vlap, V11, V22, V33, V12, V13, V23
  // V11
  glValues[0] = values.get( 0, 1 );
  // V21
  glValues[1] = values.get( 0, 4 );
  // V31
  glValues[2] = values.get( 0, 5 );
  // V12
  glValues[3] = values.get( 0, 4 );
  // V22
  glValues[4] = values.get( 0, 2 );
  // V32
  glValues[5] = values.get( 0, 6 );
  // V13
  glValues[6] = values.get( 0, 5 );
  // V23
  glValues[7] = values.get( 0, 6 );
  // V33
  glValues[8] = values.get( 0, 3 );
}

template<class LO, class SC>
void BEBilinearFormLame1Layer<LO, SC>::assembleAllMIC(
    FullMatrix<LO, SC>& VLame,
    FullMatrix<LO, SC>& VLap,
    FullMatrix<LO, SC>& KLaplace,
    BESpace< LO, SC > & bespaceK
    ) const {

#if N_MIC > 0

  // preallocate data for numerical integration
  int qOrderOuter = this->quadratureOrderDisjointElems[ 0 ];
  int qOrderInner = this->quadratureOrderDisjointElems[ 1 ];
  int outerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 0 ] ];
  int innerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 1 ] ];

  SCVT * vOutW = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * vInW = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * outerX1ref = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * outerX2ref = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * innerX1ref = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * innerX2ref = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * phi1Values = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * phi2Values = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * phi3Values = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );

  int counter = 0;
  double * qPointsIn = quadPoints[ this->quadratureOrderDisjointElems[ 1 ] ];
  for ( int i = 0; i < outerPoints; ++i ) {
    for ( int j = 0; j < innerPoints; ++j ) {
      vOutW[ counter ] =
          (SCVT) quadWeights[ this->quadratureOrderDisjointElems[ 0 ] ][ i ];
      vInW[ counter ] =
          (SCVT) quadWeights[ this->quadratureOrderDisjointElems[ 1 ] ][ j ];
      outerX1ref[ counter ] = (SCVT) quadPoints_X1[ qOrderOuter ][ i ];
      outerX2ref[ counter ] = (SCVT) quadPoints_X2[ qOrderOuter ][ i ];
      innerX1ref[ counter ] = (SCVT) quadPoints_X1[ qOrderInner ][ j ];
      innerX2ref[ counter ] = (SCVT) quadPoints_X2[ qOrderInner ][ j ];
      phi1Values[ counter ] = 1.0 - qPointsIn[ 2 * j ] - qPointsIn[ 2 * j + 1 ];
      phi2Values[ counter ] = qPointsIn[ 2 * j ];
      phi3Values[ counter ] = qPointsIn[ 2 * j + 1 ];
      ++counter;
    }
  }

  // prepare data to be sent to mic
  LO nNodes = this->space->getMesh( )->getNNodes( );
  LO nElems = this->space->getMesh( )->getNElements( );
  SCVT * nodes = this->space->getMesh( )->getNodes( )->data( );
  LO * elems = this->space->getMesh( )->getElements( )->data( );
  SCVT * areas = this->space->getMesh( )->getAreas( )->data( );
  SCVT * normals = this->space->getMesh( )->getNormals( )->data( );

  LO dim = 3;

  LO nRowsScalar = this->space->getOuterDOFs( );
  LO nColsScalar = this->space->getInnerDOFs( );

  VLame.resize( dim * nRowsScalar, dim * nColsScalar, false );
  KLaplace.resize( nElems, nNodes, false );
  VLap.resize( nElems, nElems, false );

  //FullMatrix<LO, SC> VLap( nElems, nElems, false );
  //FullMatrix<LO, SC> V11( nElems, nElems, false );
  //FullMatrix<LO, SC> V22( nElems, nElems, false );
  //FullMatrix<LO, SC> V33( nElems, nElems, false );
  //FullMatrix<LO, SC> V12( nElems, nElems, false );
  //FullMatrix<LO, SC> V13( nElems, nElems, false );
  //FullMatrix<LO, SC> V23( nElems, nElems, false );
  //FullMatrix<LO, SC> KLap( nElems, nNodes, false );

  SC * VLameData = VLame.getData( );
  SC * VLapData = VLap.getData( );
  //SC * V11Data = V11.getData( );
  //SC * V22Data = V22.getData( );
  //SC * V33Data = V33.getData( );
  //SC * V12Data = V12.getData( );
  //SC * V13Data = V13.getData( );
  //SC * V23Data = V23.getData( );
  SC * KData = KLaplace.getData( );

#pragma omp parallel for schedule( static, 32 )
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nElems; ++i ) {
      KData[ j * nElems + i ] = 0.0;
    }
  }

  // divide the global matrix among available MICs
  double CPUrelative2MIC = 1.3;
  double totalPower = N_MIC + CPUrelative2MIC;
  LO nMICRows = nElems / totalPower;
  //std::cout << nElems << " " << nMICRows << std::endl;
  LO nRowsPerMIC[ N_MIC ];
  long dataLength[ N_MIC ];
  long dataInBytes[ N_MIC ];

  for ( LO i = 0; i < N_MIC; ++i ) {
    nRowsPerMIC[ i ] = nMICRows;
    dataLength[ i ] = nMICRows * nElems;
    dataInBytes[ i ] = dataLength[ i ] * sizeof ( SC );
  }
  //nRowsPerMIC[ N_MIC - 1 ] += nElems % N_MIC;
  //dataLength[ N_MIC - 1] += ( nElems % N_MIC ) * nElems;
  //dataInBytes[ N_MIC - 1] = dataLength[ N_MIC - 1 ] * sizeof ( SC );
  LO nSubmatrices =
      std::ceil( ( (SCVT) dataInBytes[ N_MIC - 1 ] ) / MIC_CHUNK );

  LO nRowsPerSubmatrix[ N_MIC ];
  LO remainingRowsPerMIC[ N_MIC ];
  for ( LO i = 0; i < N_MIC; ++i ) {
    nRowsPerSubmatrix[ i ] = nRowsPerMIC[ i ] / nSubmatrices;
    remainingRowsPerMIC[ i ] = nRowsPerMIC[ i ] % nSubmatrices;
  }

  SC * matrixBuffersVLap[ 2 ][ N_MIC ];
  SC * matrixBuffersV11[ 2 ][ N_MIC ];
  SC * matrixBuffersV22[ 2 ][ N_MIC ];
  SC * matrixBuffersV33[ 2 ][ N_MIC ];
  SC * matrixBuffersV12[ 2 ][ N_MIC ];
  SC * matrixBuffersV13[ 2 ][ N_MIC ];
  SC * matrixBuffersV23[ 2 ][ N_MIC ];
  SC * matrixBuffersK[ 2 ][ N_MIC ];
  long bufferLengthsV[ N_MIC ];
  long bufferLengthsK[ N_MIC ];
  for ( LO i = 0; i < N_MIC; ++i ) {
    bufferLengthsV[ i ] = ( nRowsPerSubmatrix[ i ] + remainingRowsPerMIC[ i ] )
        * nElems;
    bufferLengthsK[ i ] = ( nRowsPerSubmatrix[ i ] + remainingRowsPerMIC[ i ] )
        * nNodes;
    matrixBuffersVLap[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersVLap[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV11[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV11[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV22[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV22[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV33[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV33[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV12[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV12[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV13[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV13[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV23[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV23[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersK[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsK[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersK[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsK[ i ] * sizeof ( SC ), DATA_ALIGN );
  }

  LO numCPUThreads = 0;

#pragma omp parallel
  {
#pragma omp single
    numCPUThreads = omp_get_num_threads( ) - N_MIC;
    if ( numCPUThreads < 0 ) numCPUThreads = 1;
  }

  std::vector< std::vector< LO > * > rowIdxV;
  rowIdxV.resize( numCPUThreads );
  std::vector<std::vector< LO > * > colIdxV;
  colIdxV.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesVLap;
  valuesVLap.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesV11;
  valuesV11.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesV22;
  valuesV22.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesV33;
  valuesV33.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesV12;
  valuesV12.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesV13;
  valuesV13.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesV23;
  valuesV23.resize( numCPUThreads );

  std::vector< std::vector< LO > * > rowIdxK;
  rowIdxK.resize( numCPUThreads );
  std::vector< std::vector< LO > * > colIdxK;
  colIdxK.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesK;
  valuesK.resize( numCPUThreads );

  for ( int i = 0; i < numCPUThreads; ++i ) {
    rowIdxV[ i ] = new std::vector< LO >;
    colIdxV[ i ] = new std::vector< LO >;
    valuesVLap[ i ] = new std::vector< SC >;
    valuesV11[ i ] = new std::vector< SC >;
    valuesV22[ i ] = new std::vector< SC >;
    valuesV33[ i ] = new std::vector< SC >;
    valuesV12[ i ] = new std::vector< SC >;
    valuesV13[ i ] = new std::vector< SC >;
    valuesV23[ i ] = new std::vector< SC >;
    rowIdxV[ i ]->reserve( 200 * nElems / numCPUThreads );
    colIdxV[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesVLap[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesV11[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesV22[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesV33[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesV12[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesV13[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesV23[ i ]->reserve( 200 * nElems / numCPUThreads );

    rowIdxK[ i ] = new std::vector< LO >;
    colIdxK[ i ] = new std::vector< LO >;
    valuesK[ i ] = new std::vector< SC >;
    rowIdxK[ i ]->reserve( 200 * nElems / numCPUThreads );
    colIdxK[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesK[ i ]->reserve( 200 * nElems / numCPUThreads );
  }





#pragma omp parallel
  {
    if ( omp_get_thread_num( ) < N_MIC ) {
 int device = omp_get_thread_num( );
    SC * myDataVLap1 = matrixBuffersVLap[ 0 ][ device ];
    SC * myDataVLap2 = matrixBuffersVLap[ 1 ][ device ];
    SC * myDataV11_1 = matrixBuffersV11[ 0 ][ device ];
    SC * myDataV11_2 = matrixBuffersV11[ 1 ][ device ];
    SC * myDataV22_1 = matrixBuffersV22[ 0 ][ device ];
    SC * myDataV22_2 = matrixBuffersV22[ 1 ][ device ];
    SC * myDataV33_1 = matrixBuffersV33[ 0 ][ device ];
    SC * myDataV33_2 = matrixBuffersV33[ 1 ][ device ];
    SC * myDataV12_1 = matrixBuffersV12[ 0 ][ device ];
    SC * myDataV12_2 = matrixBuffersV12[ 1 ][ device ];
    SC * myDataV13_1 = matrixBuffersV13[ 0 ][ device ];
    SC * myDataV13_2 = matrixBuffersV13[ 1 ][ device ];
    SC * myDataV23_1 = matrixBuffersV23[ 0 ][ device ];
    SC * myDataV23_2 = matrixBuffersV23[ 1 ][ device ];
    SC * myDataK1 = matrixBuffersK[ 0 ][ device ];
    SC * myDataK2 = matrixBuffersK[ 1 ][ device ];

#pragma offload_transfer target( mic : device ) \
    nocopy( myDataVLap1 : length( bufferLengthsV[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataVLap2 : length( bufferLengthsV[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV11_1 : length( bufferLengthsV[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV11_2 : length( bufferLengthsV[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV22_1 : length( bufferLengthsV[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV22_2 : length( bufferLengthsV[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV33_1 : length( bufferLengthsV[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV33_2 : length( bufferLengthsV[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV12_1 : length( bufferLengthsV[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV12_2 : length( bufferLengthsV[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV13_1 : length( bufferLengthsV[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV13_2 : length( bufferLengthsV[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV23_1 : length( bufferLengthsV[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV23_2 : length( bufferLengthsV[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataK1 : length( bufferLengthsK[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataK2 : length( bufferLengthsK[ device ] ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    in( nodes : length( 3 * nNodes ) alloc_if( 1 ) free_if( 0 ) ) \
    in( elems : length( 3 * nElems ) alloc_if( 1 ) free_if( 0 ) ) \
    in( areas : length( nElems ) alloc_if( 1 ) free_if( 0 ) ) \
    in( normals : length( 3 * nElems ) alloc_if( 1 ) free_if( 0 ) ) \
    in( phi1Values : length( outerPoints * innerPoints ) \
    alloc_if(1) free_if(0)) \
    in( phi2Values : length( outerPoints * innerPoints ) \
    alloc_if(1) free_if(0)) \
    in( phi3Values : length( outerPoints * innerPoints ) \
    alloc_if(1) free_if(0)) \
    in( vOutW : length( outerPoints * innerPoints ) \
    alloc_if( 1 ) free_if( 0 ) ) \
    in( vInW : length( outerPoints * innerPoints ) \
    alloc_if( 1 ) free_if( 0 ) ) \
    in( nRowsPerSubmatrix : alloc_if( 1 ) free_if( 0 ) ) \
    in( remainingRowsPerMIC : alloc_if( 1 ) free_if( 0 ) ) \
    in( outerX1ref : length( outerPoints * innerPoints ) \
    alloc_if( 1 ) free_if( 0 ) ) \
    in( outerX2ref : length( outerPoints * innerPoints ) \
    alloc_if( 1 ) free_if( 0 ) ) \
    in( innerX1ref : length( outerPoints * innerPoints ) \
    alloc_if( 1 ) free_if( 0 ) ) \
    in( innerX2ref : length( outerPoints * innerPoints ) \
    alloc_if( 1 ) free_if( 0 ) ) \
    in( this : alloc_if( 1 ) free_if( 0 ) ) 

      // double buffering (one buffer for computation, one buffer sent to cpu)
      for ( LO s = 0; s <= nSubmatrices; ++s ) {
//std::cout << s  << std::endl;
        if ( s != nSubmatrices ) {
          SC * computationDataVLap = matrixBuffersVLap[ s % 2 ][ device ];
          SC * computationDataV11 = matrixBuffersV11[ s % 2 ][ device ];
          SC * computationDataV22 = matrixBuffersV22[ s % 2 ][ device ];
          SC * computationDataV33 = matrixBuffersV33[ s % 2 ][ device ];
          SC * computationDataV12 = matrixBuffersV12[ s % 2 ][ device ];
          SC * computationDataV13 = matrixBuffersV13[ s % 2 ][ device ];
          SC * computationDataV23 = matrixBuffersV23[ s % 2 ][ device ];
          SC * computationDataK = matrixBuffersK[ s % 2 ][ device ];

      //signal( matrixBuffersVLap[ s % 2 ][ device ] ) 
#pragma offload target( mic : device ) \
      in( computationDataVLap : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( computationDataV11 : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( computationDataV22 : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( computationDataV33 : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( computationDataV12 : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( computationDataV13 : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( computationDataV23 : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( computationDataK : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( nodes : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( elems : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( areas : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( normals : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( phi1Values : length( 0 ) alloc_if(0) free_if(0)) \
      in( phi2Values : length( 0 ) alloc_if(0) free_if(0)) \
      in( phi3Values : length( 0 ) alloc_if(0) free_if(0)) \
      in( vOutW : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( vInW : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( nRowsPerSubmatrix :  alloc_if( 0 ) free_if( 0 ) ) \
      in( remainingRowsPerMIC : alloc_if( 0 ) free_if( 0 ) ) \
      in( outerX1ref : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( outerX2ref : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( innerX1ref : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( innerX2ref : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( s ) \
      signal( computationDataVLap ) \
      in( this : length( 0 ) alloc_if( 0 ) free_if( 0 ) )
          {
#pragma omp parallel for schedule( static, 32 ) num_threads( N_MIC_THREADS )
#pragma vector aligned(computationDataK)
            for ( LO i = 0; i < bufferLengthsK[ device ]; ++i ) {
              computationDataK[ i ] = 0.0;
            }

#pragma omp parallel num_threads( N_MIC_THREADS )
            {

              SCVT * outerX1 = (SCVT *) _mm_malloc(
                  outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
              SCVT * outerX2 = (SCVT *) _mm_malloc(
                  outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
              SCVT * outerX3 = (SCVT *) _mm_malloc(
                  outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
              SCVT * innerX1 = (SCVT *) _mm_malloc(
                  outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
              SCVT * innerX2 = (SCVT *) _mm_malloc(
                  outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
              SCVT * innerX3 = (SCVT *) _mm_malloc(
                  outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );

              SC elemMatrixVLap = 0.0;
              SC elemMatrixV11 = 0.0;
              SC elemMatrixV22 = 0.0;
              SC elemMatrixV33 = 0.0;
              SC elemMatrixV12 = 0.0;
              SC elemMatrixV13 = 0.0;
              SC elemMatrixV23 = 0.0;
              SC elemMatrixK[3] = { 0.0, 0.0, 0.0 };
              LO startRow = device * nMICRows + s * nRowsPerSubmatrix[ device ];
              LO endRow = startRow + nRowsPerSubmatrix[ device ];
              endRow +=
                  ( s == nSubmatrices - 1 ) ? remainingRowsPerMIC[ device ] : 0;

//#pragma omp for collapse( 2 ) schedule( dynamic, 32 )

//#pragma novector
#pragma omp for schedule( dynamic, 16 )
                for ( LO j = 0; j < nElems; ++j ) {
                    
              for ( LO i = startRow; i < endRow; ++i ) {
                  if ( !BEIntegrator<LO, SC, BEIntegratorLaplace<LO, SC> >::
                      areElementsDisjointOnMIC( i, j, nodes, elems ) ) {
                    continue;
                  }

                  BEIntegratorLame<LO, SC >::
                      computeElemMatrixAllMIC(
                      nodes, elems, areas, normals, i, j, qOrderOuter,
                      qOrderInner, outerX1ref, outerX2ref, innerX1ref, innerX2ref,
                      outerX1, outerX2, outerX3, innerX1, innerX2, innerX3,
                      phi1Values, phi2Values, phi3Values, vOutW, vInW,
                      &elemMatrixV11, &elemMatrixV22,
                      &elemMatrixV33, &elemMatrixV12, &elemMatrixV13,
                      &elemMatrixV23, &elemMatrixVLap, elemMatrixK );

                  computationDataVLap[ j * ( endRow - startRow ) + i - startRow ] =
                      elemMatrixVLap;
                  computationDataV11[ j * ( endRow - startRow ) + i - startRow ] =
                      elemMatrixV11;
                  computationDataV22[ j * ( endRow - startRow ) + i - startRow ] =
                      elemMatrixV22;
                  computationDataV33[ j * ( endRow - startRow ) + i - startRow ] =
                      elemMatrixV33;
                  computationDataV12[ j * ( endRow - startRow ) + i - startRow ] =
                      elemMatrixV12;
                  computationDataV13[ j * ( endRow - startRow ) + i - startRow ] =
                      elemMatrixV13;
                  computationDataV23[ j * ( endRow - startRow ) + i - startRow ] =
                      elemMatrixV23;
#pragma omp atomic update
                  computationDataK[ elems[ 3 * j ] *
                      ( endRow - startRow ) + i - startRow ] += elemMatrixK[ 0 ];
#pragma omp atomic update
                  computationDataK[ elems[ 3 * j + 1 ] *
                      ( endRow - startRow ) + i - startRow ] += elemMatrixK[ 1 ];
#pragma omp atomic update
                  computationDataK[ elems[ 3 * j + 2 ] *
                      ( endRow - startRow ) + i - startRow ] += elemMatrixK[ 2 ];
              
                }
              }
              _mm_free( outerX1 );
              _mm_free( outerX2 );
              _mm_free( outerX3 );
              _mm_free( innerX1 );
              _mm_free( innerX2 );
              _mm_free( innerX3 );
            }
          }
        }
        // send data to CPU
        if ( s > 0 ) {
          SC * outputBufferVLap = matrixBuffersVLap[ ( s - 1 ) % 2 ][ device ];
          SC * outputBufferV11 = matrixBuffersV11[ ( s - 1 ) % 2 ][ device ];
          SC * outputBufferV22 = matrixBuffersV22[ ( s - 1 ) % 2 ][ device ];
          SC * outputBufferV33 = matrixBuffersV33[ ( s - 1 ) % 2 ][ device ];
          SC * outputBufferV12 = matrixBuffersV12[ ( s - 1 ) % 2 ][ device ];
          SC * outputBufferV13 = matrixBuffersV13[ ( s - 1 ) % 2 ][ device ];
          SC * outputBufferV23 = matrixBuffersV23[ ( s - 1 ) % 2 ][ device ];
          SC * outputBufferK = matrixBuffersK[ ( s - 1 ) % 2 ][ device ];
          
#pragma offload_transfer target( mic : device ) \
        wait(outputBufferVLap) \
        out( outputBufferVLap : length( bufferLengthsV[ device ] ) free_if(0)) \
        out( outputBufferV11 : length( bufferLengthsV[ device ] ) free_if(0)) \
        out( outputBufferV22 : length( bufferLengthsV[ device ] ) free_if(0)) \
        out( outputBufferV33 : length( bufferLengthsV[ device ] ) free_if(0)) \
        out( outputBufferV12 : length( bufferLengthsV[ device ] ) free_if(0)) \
        out( outputBufferV13 : length( bufferLengthsV[ device ] ) free_if(0)) \
        out( outputBufferV23 : length( bufferLengthsV[ device ] ) free_if(0)) \
        out( outputBufferK : length( bufferLengthsK[ device ] ) free_if(0))


          for ( LO j = 0; j < nElems; ++j ) {
            LO rem = ( s == nSubmatrices - 1 ) ? remainingRowsPerMIC[ device ] : 0;
            memcpy( VLapData + j * nElems + device * nMICRows +
                ( s - 1 ) * nRowsPerSubmatrix[ device ],
                matrixBuffersVLap[ ( s - 1 ) % 2 ][ device ]
                + j * ( nRowsPerSubmatrix[ device ] + rem ),
                ( nRowsPerSubmatrix[ device ] + rem ) * sizeof ( SC ) );
            memcpy( VLameData + j * nElems * dim + device * nMICRows +
                ( s - 1 ) * nRowsPerSubmatrix[ device ],
                matrixBuffersV11[ ( s - 1 ) % 2 ][ device ]
                + j * ( nRowsPerSubmatrix[ device ] + rem ),
                ( nRowsPerSubmatrix[ device ] + rem ) * sizeof ( SC ) );
            memcpy( VLameData + ( 3 * nElems * nElems + nElems ) +
                j * nElems * dim + device * nMICRows +
                ( s - 1 ) * nRowsPerSubmatrix[ device ],
                matrixBuffersV22[ ( s - 1 ) % 2 ][ device ]
                + j * ( nRowsPerSubmatrix[ device ] + rem ),
                ( nRowsPerSubmatrix[ device ] + rem ) * sizeof ( SC ) );
            memcpy( VLameData + ( 6 * nElems * nElems + 2 * nElems ) +
                j * nElems * dim + device * nMICRows +
                ( s - 1 ) * nRowsPerSubmatrix[ device ],
                matrixBuffersV33[ ( s - 1 ) % 2 ][ device ]
                + j * ( nRowsPerSubmatrix[ device ] + rem ),
                ( nRowsPerSubmatrix[ device ] + rem ) * sizeof ( SC ) );
            memcpy( VLameData + ( 3 * nElems * nElems ) +
                j * nElems * dim + device * nMICRows +
                ( s - 1 ) * nRowsPerSubmatrix[ device ],
                matrixBuffersV12[ ( s - 1 ) % 2 ][ device ]
                + j * ( nRowsPerSubmatrix[ device ] + rem ),
                ( nRowsPerSubmatrix[ device ] + rem ) * sizeof ( SC ) );
            memcpy( VLameData + ( 1 * nElems ) +
                j * nElems * dim + device * nMICRows +
                ( s - 1 ) * nRowsPerSubmatrix[ device ],
                matrixBuffersV12[ ( s - 1 ) % 2 ][ device ]
                + j * ( nRowsPerSubmatrix[ device ] + rem ),
                ( nRowsPerSubmatrix[ device ] + rem ) * sizeof ( SC ) );
            memcpy( VLameData + ( 6 * nElems * nElems ) +
                j * nElems * dim + device * nMICRows +
                ( s - 1 ) * nRowsPerSubmatrix[ device ],
                matrixBuffersV13[ ( s - 1 ) % 2 ][ device ]
                + j * ( nRowsPerSubmatrix[ device ] + rem ),
                ( nRowsPerSubmatrix[ device ] + rem ) * sizeof ( SC ) );
            memcpy( VLameData + ( 2 * nElems ) +
                j * nElems * dim + device * nMICRows +
                ( s - 1 ) * nRowsPerSubmatrix[ device ],
                matrixBuffersV13[ ( s - 1 ) % 2 ][ device ]
                + j * ( nRowsPerSubmatrix[ device ] + rem ),
                ( nRowsPerSubmatrix[ device ] + rem ) * sizeof ( SC ) );
            memcpy( VLameData + ( 6 * nElems * nElems + nElems ) +
                j * nElems * dim + device * nMICRows +
                ( s - 1 ) * nRowsPerSubmatrix[ device ],
                matrixBuffersV23[ ( s - 1 ) % 2 ][ device ]
                + j * ( nRowsPerSubmatrix[ device ] + rem ),
                ( nRowsPerSubmatrix[ device ] + rem ) * sizeof ( SC ) );
            memcpy( VLameData + ( 3 * nElems * nElems + 2 * nElems ) +
                j * nElems * dim + device * nMICRows +
                ( s - 1 ) * nRowsPerSubmatrix[ device ],
                matrixBuffersV23[ ( s - 1 ) % 2 ][ device ]
                + j * ( nRowsPerSubmatrix[ device ] + rem ),
                ( nRowsPerSubmatrix[ device ] + rem ) * sizeof ( SC ) );
          }
          for ( LO k = 0; k < nNodes; ++k ) {
            LO rem = ( s == nSubmatrices  - 1) ? remainingRowsPerMIC[ device ] : 0;
            memcpy( KData + k * nElems + device * nMICRows +
                ( s - 1 ) * nRowsPerSubmatrix[ device ],
                matrixBuffersK[ ( s - 1 ) % 2 ][ device ] +
                k * ( nRowsPerSubmatrix[ device ] + rem ),
                ( nRowsPerSubmatrix[ device ] + rem ) * sizeof ( SC ) );
          }
       }
      }

      //SC * myDataVLap1 = matrixBuffersVLap[ 0 ][ device ];
      //SC * myDataVLap2 = matrixBuffersVLap[ 1 ][ device ];
      //SC * myDataV11_1 = matrixBuffersV11[ 0 ][ device ];
      //SC * myDataV11_2 = matrixBuffersV11[ 1 ][ device ];
      //SC * myDataV22_1 = matrixBuffersV22[ 0 ][ device ];
      //SC * myDataV22_2 = matrixBuffersV22[ 1 ][ device ];
      //SC * myDataV33_1 = matrixBuffersV33[ 0 ][ device ];
      //SC * myDataV33_2 = matrixBuffersV33[ 1 ][ device ];
      //SC * myDataV12_1 = matrixBuffersV12[ 0 ][ device ];
      //SC * myDataV12_2 = matrixBuffersV12[ 1 ][ device ];
      //SC * myDataV13_1 = matrixBuffersV13[ 0 ][ device ];
      //SC * myDataV13_2 = matrixBuffersV13[ 1 ][ device ];
      //SC * myDataV23_1 = matrixBuffersV23[ 0 ][ device ];
      //SC * myDataV23_2 = matrixBuffersV23[ 1 ][ device ];
      //SC * myDataK1 = matrixBuffersK[ 0 ][ device ];
      //SC * myDataK2 = matrixBuffersK[ 1 ][ device ];
#pragma offload_transfer target( mic : device ) \
      nocopy( myDataVLap1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataVLap2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV11_1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV11_2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV22_1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV22_2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV33_1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV33_2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV12_1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV12_2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV13_1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV13_2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV23_1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV23_2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataK1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataK2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( nodes : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( elems : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( areas : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( vOutW : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( vInW : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( outerX1ref : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( outerX2ref : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( innerX1ref : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( innerX2ref : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( phi1Values : length( 0 ) alloc_if(0) free_if(1)) \
      nocopy( phi2Values : length( 0 ) alloc_if(0) free_if(1)) \
      nocopy( phi3Values : length( 0 ) alloc_if(0) free_if(1))

      std::cout << "MIC FINISHED" << std::endl;
    } else {
#pragma omp single nowait
      {
          std::cout << "CPU STARTED"<< std::endl;
        for ( LO i = 0; i < nElems; ++i ) {
#pragma omp task
          {
            int myThreadNum = omp_get_thread_num( );
            BEIntegratorLame<LO, SC> integratorV( this->space,
                this->quadratureOrder, this->quadrature,
                this->quadratureOrderDisjointElems );
            BEIntegratorLaplace<LO, SC> integratorK( &bespaceK,
                this->quadratureOrder, this->quadrature,
                this->quadratureOrderDisjointElems );
            FullMatrix< LO, SC > elemMatrixV( 1, 7, true );
            FullMatrix< LO, SC > elemMatrixK( 1, 3, true );

            for ( LO j = 0; j < nElems; ++j ) {
              if ( BEIntegrator<LO, SC, BEIntegratorLaplace<LO, SC> >::
                  areElementsDisjointOnMIC( i, j, nodes, elems ) ) {
                continue;
              }
              integratorV.getElemMatrix1Layer( i, j, elemMatrixV );
              integratorK.getElemMatrix2Layer( i, j, elemMatrixK );

              rowIdxV[ myThreadNum - N_MIC ]->push_back( i );
              colIdxV[ myThreadNum - N_MIC ]->push_back( j );
              valuesVLap[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixV.get( 0, 0 ) );
              valuesV11[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixV.get( 0, 1 ) );
              valuesV22[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixV.get( 0, 2 ) );
              valuesV33[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixV.get( 0, 3 ) );
              valuesV12[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixV.get( 0, 4 ) );
              valuesV13[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixV.get( 0, 5 ) );
              valuesV23[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixV.get( 0, 6 ) );

              rowIdxK[ myThreadNum - N_MIC ]->push_back( i );
              rowIdxK[ myThreadNum - N_MIC ]->push_back( i );
              rowIdxK[ myThreadNum - N_MIC ]->push_back( i );
              colIdxK[ myThreadNum - N_MIC ]->push_back( elems[ 3 * j ] );
              colIdxK[ myThreadNum - N_MIC ]->push_back( elems[ 3 * j + 1 ] );
              colIdxK[ myThreadNum - N_MIC ]->push_back( elems[ 3 * j + 2 ] );
              valuesK[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixK.get( 0, 0 ) );
              valuesK[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixK.get( 0, 1 ) );
              valuesK[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixK.get( 0, 2 ) );
            } // end for j elems
          } // end omp task
        } // end i elems
        LO startRow = N_MIC * nMICRows;
            for ( LO j = 0; j < nElems; ++j ) {
#pragma omp task
          {
            SCVT * outerX1 = (SCVT *) _mm_malloc(
                outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
            SCVT * outerX2 = (SCVT *) _mm_malloc(
                outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
            SCVT * outerX3 = (SCVT *) _mm_malloc(
                outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
            SCVT * innerX1 = (SCVT *) _mm_malloc(
                outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
            SCVT * innerX2 = (SCVT *) _mm_malloc(
                outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
            SCVT * innerX3 = (SCVT *) _mm_malloc(
                outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );

            SC elemMatrixVLap = 0.0;
            SC elemMatrixV11 = 0.0;
            SC elemMatrixV22 = 0.0;
            SC elemMatrixV33 = 0.0;
            SC elemMatrixV12 = 0.0;
            SC elemMatrixV13 = 0.0;
            SC elemMatrixV23 = 0.0;
            SC elemMatrixK[3] = { 0.0, 0.0, 0.0 };

        for ( LO i = startRow; i < nElems; ++i ) {

              if ( !BEIntegrator<LO, SC, BEIntegratorLaplace<LO, SC> >::
                  areElementsDisjointOnMIC( i, j, nodes, elems ) ) {
                continue;
              }

              BEIntegratorLame<LO, SC >::
                  computeElemMatrixAllMIC(
                  nodes, elems, areas, normals, i, j, qOrderOuter,
                  qOrderInner, outerX1ref, outerX2ref, innerX1ref, innerX2ref,
                  outerX1, outerX2, outerX3, innerX1, innerX2, innerX3,
                  phi1Values, phi2Values, phi3Values, vOutW, vInW,
                  &elemMatrixV11, &elemMatrixV22,
                  &elemMatrixV33, &elemMatrixV12, &elemMatrixV13,
                  &elemMatrixV23, &elemMatrixVLap, elemMatrixK );

              VLap.set( i, j, elemMatrixVLap );

#pragma omp atomic update
              KData[ elems[ 3 * j ] *
                  nElems + i ] += elemMatrixK[0];
#pragma omp atomic update
              KData[ elems[ 3 * j + 1 ] *
                  nElems + i ] += elemMatrixK[1];
#pragma omp atomic update
              KData[ elems[ 3 * j + 2 ] *
                  nElems + i ] += elemMatrixK[2];

              FullMatrix<LO, SC> tmpMatrix( 1, 7 );
              // vector of global row indices in matrix
              std::vector<LO> glRInd( 9 );
              // vector of global column indices in matrix
              std::vector<LO> glCInd( 9 );
              // vector of global values in matrix
              std::vector<SC> glVal( 9 );

              tmpMatrix.set( 0, 0, elemMatrixVLap );
              tmpMatrix.set( 0, 1, elemMatrixV11 );
              tmpMatrix.set( 0, 2, elemMatrixV22 );
              tmpMatrix.set( 0, 3, elemMatrixV33 );
              tmpMatrix.set( 0, 4, elemMatrixV12 );
              tmpMatrix.set( 0, 5, elemMatrixV13 );
              tmpMatrix.set( 0, 6, elemMatrixV23 );
              this->getGlobalIndices( i,
                  j, tmpMatrix,
                  glRInd, glCInd, glVal );
              for ( LO i = 0; i < 9; i++ ) {
                VLame.set( glRInd[i], glCInd[i], glVal[i] );
              }

            }

            _mm_free( outerX1 );
            _mm_free( outerX2 );
            _mm_free( outerX3 );
            _mm_free( innerX1 );
            _mm_free( innerX2 );
            _mm_free( innerX3 );
          } // end omp task
        } // end i elems
        std::cout << "CPU  FINISHED" << std::endl;
      } // end omp single nowait
    } // end else
#pragma omp barrier
  } // end omp parallel
#pragma omp parallel num_threads( numCPUThreads )
  {
    int myThreadNum = omp_get_thread_num( );
    FullMatrix<LO, SC> tmpMatrix( 1, 7 );
    // vector of global row indices in matrix
    std::vector<LO> glRInd( 9 );
    // vector of global column indices in matrix
    std::vector<LO> glCInd( 9 );
    // vector of global values in matrix
    std::vector<SC> glVal( 9 );
    for ( LO j = 0; j < rowIdxV[ myThreadNum]->size( ); ++j ) {
      tmpMatrix.set( 0, 0, valuesVLap[ myThreadNum ]->at( j ) );
      tmpMatrix.set( 0, 1, valuesV11[ myThreadNum ]->at( j ) );
      tmpMatrix.set( 0, 2, valuesV22[ myThreadNum ]->at( j ) );
      tmpMatrix.set( 0, 3, valuesV33[ myThreadNum ]->at( j ) );
      tmpMatrix.set( 0, 4, valuesV12[ myThreadNum ]->at( j ) );
      tmpMatrix.set( 0, 5, valuesV13[ myThreadNum ]->at( j ) );
      tmpMatrix.set( 0, 6, valuesV23[ myThreadNum ]->at( j ) );
      this->getGlobalIndices( rowIdxV[ myThreadNum ]->at( j ),
          colIdxV[ myThreadNum ]->at( j ), tmpMatrix,
          glRInd, glCInd, glVal );
      for ( LO i = 0; i < 9; i++ ) {
        VLame.set( glRInd[i], glCInd[i], glVal[i] );
      }
      //VLame.addToPositions( glRInd, glCInd, glVal );
      VLap.set( rowIdxV[ myThreadNum ]->at( j ), colIdxV[ myThreadNum ]->at( j ),
          valuesVLap[ myThreadNum ]->at( j ) );
    }
    for ( LO j = 0; j < rowIdxK[myThreadNum]->size( ); ++j ) {
#pragma omp atomic update
      KData[ colIdxK[ myThreadNum ]->at( j ) * nElems +
          rowIdxK[ myThreadNum ]->at( j ) ] += valuesK[ myThreadNum ]->at( j );
    }
  }

  // add the Laplace matrix V to the diagonal of the result
  SC val = 0.0;
  LO nRows = 3 * nRowsScalar;
  LO i = 0;
#pragma omp parallel for private(val)
  for ( LO j = 0; j < nColsScalar; ++j ) {

#pragma omp simd \
linear( i : 1 ) \
private( val ) \
simdlen( DATA_WIDTH )
    for ( i = 0; i < nRowsScalar; ++i ) {
      val = VLapData[ j * nRowsScalar + i ] * ( 3.0 - 4.0 * this->nu ); 
      VLameData[ j * nRows + i ] += val;
      VLameData[ (j+nColsScalar) * nRows + i + nRowsScalar ] += val;
      VLameData[ (j+2 * nColsScalar) * nRows + i + 2 * nRowsScalar ] += val;
    }
  }
  VLame.scale( ( 1.0 + this->nu ) / ( 2.0 * E * ( 1.0 - this->nu ) ) );
  /*SCVT mult = ( (SCVT) 1.0 + this->nu ) /
      ( (SCVT) 2.0 * E * ( (SCVT) 1.0 - this->nu ) );
      
#pragma omp parallel for simd 
  for ( LO i = 0; i < nRows * nRows; ++i ) {
    VLameData[ i ] *= mult;
  }*/
  for ( LO i = 0; i < N_MIC; ++i ) {
    _mm_free( matrixBuffersVLap[ 0 ][ i ] );
    _mm_free( matrixBuffersVLap[ 1 ][ i ] );
    _mm_free( matrixBuffersV11[ 0 ][ i ] );
    _mm_free( matrixBuffersV11[ 1 ][ i ] );
    _mm_free( matrixBuffersV22[ 0 ][ i ] );
    _mm_free( matrixBuffersV22[ 1 ][ i ] );
    _mm_free( matrixBuffersV33[ 0 ][ i ] );
    _mm_free( matrixBuffersV33[ 1 ][ i ] );
    _mm_free( matrixBuffersV12[ 0 ][ i ] );
    _mm_free( matrixBuffersV12[ 1 ][ i ] );
    _mm_free( matrixBuffersV13[ 0 ][ i ] );
    _mm_free( matrixBuffersV13[ 1 ][ i ] );
    _mm_free( matrixBuffersV23[ 0 ][ i ] );
    _mm_free( matrixBuffersV23[ 1 ][ i ] );
    _mm_free( matrixBuffersK[ 0 ][ i ] );
    _mm_free( matrixBuffersK[ 1 ][ i ] );
  }

  _mm_free( vOutW );
  _mm_free( vInW );
  _mm_free( outerX1ref );
  _mm_free( outerX2ref );
  _mm_free( innerX1ref );
  _mm_free( innerX2ref );
  _mm_free( phi1Values );
  _mm_free( phi2Values );
  _mm_free( phi3Values );

  for ( int i = 0; i < numCPUThreads; ++i ) {
    delete rowIdxV[ i ];
    delete colIdxV[ i ];
    delete valuesVLap[ i ];
    delete valuesV11[ i ];
    delete valuesV22[ i ];
    delete valuesV33[ i ];
    delete valuesV12[ i ];
    delete valuesV13[ i ];
    delete valuesV23[ i ];
    delete rowIdxK[ i ];
    delete colIdxK[ i ];
    delete valuesK[ i ];
  }

#endif

}

template<class LO, class SC>
void BEBilinearFormLame1Layer<LO, SC>::assembleAllMIC(
    BlockMatrix<LO, SC>& VLame,
    FullMatrix<LO, SC>& VLap,
    FullMatrix<LO, SC>& KLaplace,
    BESpace< LO, SC > & bespaceK
    ) const {

#if N_MIC > 0
  // preallocate data for numerical integration
  int qOrderOuter = this->quadratureOrderDisjointElems[ 0 ];
  int qOrderInner = this->quadratureOrderDisjointElems[ 1 ];
  int outerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 0 ] ];
  int innerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 1 ] ];

  SCVT * vOutW = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * vInW = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * outerX1ref = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * outerX2ref = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * innerX1ref = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * innerX2ref = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * phi1Values = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * phi2Values = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
  SCVT * phi3Values = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );

  int counter = 0;
  double * qPointsIn = quadPoints[ this->quadratureOrderDisjointElems[ 1 ] ];
  for ( int i = 0; i < outerPoints; ++i ) {
    for ( int j = 0; j < innerPoints; ++j ) {
      vOutW[ counter ] =
          (SCVT) quadWeights[ this->quadratureOrderDisjointElems[ 0 ] ][ i ];
      vInW[ counter ] =
          (SCVT) quadWeights[ this->quadratureOrderDisjointElems[ 1 ] ][ j ];
      outerX1ref[ counter ] = (SCVT) quadPoints_X1[ qOrderOuter ][ i ];
      outerX2ref[ counter ] = (SCVT) quadPoints_X2[ qOrderOuter ][ i ];
      innerX1ref[ counter ] = (SCVT) quadPoints_X1[ qOrderInner ][ j ];
      innerX2ref[ counter ] = (SCVT) quadPoints_X2[ qOrderInner ][ j ];
      phi1Values[ counter ] = 1.0 - qPointsIn[ 2 * j ] - qPointsIn[ 2 * j + 1 ];
      phi2Values[ counter ] = qPointsIn[ 2 * j ];
      phi3Values[ counter ] = qPointsIn[ 2 * j + 1 ];
      ++counter;
    }
  }
  // prepare data to be sent to mic
  LO nNodes = this->space->getMesh( )->getNNodes( );
  LO nElems = this->space->getMesh( )->getNElements( );
  SCVT * nodes = this->space->getMesh( )->getNodes( )->data( );
  LO * elems = this->space->getMesh( )->getElements( )->data( );
  SCVT * areas = this->space->getMesh( )->getAreas( )->data( );
  SCVT * normals = this->space->getMesh( )->getNormals( )->data( );

  LO dim = 3;

  LO nRowsScalar = this->space->getOuterDOFs( );
  LO nColsScalar = this->space->getInnerDOFs( );

  LO rows[3];
  rows[0] = nElems; 
  rows[1] = nElems;
  rows[2] = nElems;

  VLame.resize( dim, dim, rows, rows );
  FullMatrix<LO, SC> * VLame11 = new FullMatrix<LO, SC>( nElems, nElems, false );
  FullMatrix<LO, SC> * VLame22 = new FullMatrix<LO, SC>( nElems, nElems, false );
  FullMatrix<LO, SC> * VLame33 = new FullMatrix<LO, SC>( nElems, nElems, false );
  FullMatrix<LO, SC> * VLame12 = new FullMatrix<LO, SC>( nElems, nElems, false );
  FullMatrix<LO, SC> * VLame13 = new FullMatrix<LO, SC>( nElems, nElems, false );
  FullMatrix<LO, SC> * VLame23 = new FullMatrix<LO, SC>( nElems, nElems, false );

  KLaplace.resize( nElems, nNodes, false );
  VLap.resize( nElems, nElems, false );

  SC * VLapData = VLap.getData( );
  SC * VLame11Data = VLame11->getData();
  SC * VLame22Data = VLame22->getData();
  SC * VLame33Data = VLame33->getData();
  SC * VLame12Data = VLame12->getData();
  SC * VLame13Data = VLame13->getData();
  SC * VLame23Data = VLame23->getData();
  SC * KData = KLaplace.getData( );

#pragma omp parallel for schedule( static, 32 )
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nElems; ++i ) {
      KData[ j * nElems + i ] = 0.0;
    }
  }

  // divide the global matrix among available MICs
  double CPUrelative2MIC = 0.7;
  double totalPower = N_MIC + CPUrelative2MIC;
  LO nMICRows = nElems / totalPower;
  LO nRowsPerMIC[ N_MIC ];
  long dataLength[ N_MIC ];
  long dataInBytes[ N_MIC ];
  for ( LO i = 0; i < N_MIC; ++i ) {
    nRowsPerMIC[ i ] = nMICRows;
    dataLength[ i ] = nMICRows * nElems;
    dataInBytes[ i ] = dataLength[ i ] * sizeof ( SC );
  }
  //nRowsPerMIC[ N_MIC - 1 ] += nElems % N_MIC;
  //dataLength[ N_MIC - 1] += ( nElems % N_MIC ) * nElems;
  //dataInBytes[ N_MIC - 1] = dataLength[ N_MIC - 1 ] * sizeof ( SC );
  LO nSubmatrices =
      std::ceil( ( (SCVT) dataInBytes[ N_MIC - 1 ] ) / MIC_CHUNK );

  LO nRowsPerSubmatrix[ N_MIC ];
  LO remainingRowsPerMIC[ N_MIC ];
  for ( LO i = 0; i < N_MIC; ++i ) {
    nRowsPerSubmatrix[ i ] = nRowsPerMIC[ i ] / nSubmatrices;
    remainingRowsPerMIC[ i ] = nRowsPerMIC[ i ] % nSubmatrices;
  }

  // total lenght of i,j buffers for close elements
  LO closeBuffLen = nRowsPerMIC[ 0 ] * 24;
  std::vector<LO> closeI;
  std::vector<LO> closeJ;

  SC * matrixBuffersVLap[ 2 ][ N_MIC ];
  SC * matrixBuffersV11[ 2 ][ N_MIC ];
  SC * matrixBuffersV22[ 2 ][ N_MIC ];
  SC * matrixBuffersV33[ 2 ][ N_MIC ];
  SC * matrixBuffersV12[ 2 ][ N_MIC ];
  SC * matrixBuffersV13[ 2 ][ N_MIC ];
  SC * matrixBuffersV23[ 2 ][ N_MIC ];
  SC * matrixBuffersK[ 2 ][ N_MIC ];
  LO * closeElemsBufferI[ 2 ][ N_MIC ];
  LO * closeElemsBufferJ[ 2 ][ N_MIC ];
  long bufferLengthsV[ N_MIC ];
  long bufferLengthsK[ N_MIC ];
  for ( LO i = 0; i < N_MIC; ++i ) {
    bufferLengthsV[ i ] = ( nRowsPerSubmatrix[ i ] + remainingRowsPerMIC[ i ] )
        * nElems;
    bufferLengthsK[ i ] = ( nRowsPerSubmatrix[ i ] + remainingRowsPerMIC[ i ] )
        * nNodes;
    matrixBuffersVLap[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersVLap[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV11[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV11[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV22[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV22[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV33[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV33[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV12[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV12[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV13[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV13[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV23[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersV23[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsV[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersK[ 0 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsK[ i ] * sizeof ( SC ), DATA_ALIGN );
    matrixBuffersK[ 1 ][ i ] =
        (SC *) _mm_malloc( bufferLengthsK[ i ] * sizeof ( SC ), DATA_ALIGN );
    closeElemsBufferI[ 0 ][ i ] = 
        (LO *) _mm_malloc( closeBuffLen * sizeof( LO), DATA_ALIGN );
    closeElemsBufferI[ 1 ][ i ] = 
        (LO *) _mm_malloc( closeBuffLen * sizeof( LO), DATA_ALIGN );
    closeElemsBufferJ[ 0 ][ i ] = 
        (LO *) _mm_malloc( closeBuffLen * sizeof( LO), DATA_ALIGN );
    closeElemsBufferJ[ 1 ][ i ] = 
        (LO *) _mm_malloc( closeBuffLen * sizeof( LO), DATA_ALIGN );
  }

  
  LO numCPUThreads = 0;
  LO numClosePanels[] = {0, 0}; 

#pragma omp parallel
  {
#pragma omp single
    numCPUThreads = omp_get_num_threads( ) - N_MIC;
    if ( numCPUThreads < 0 ) numCPUThreads = 1;
  }

  std::vector< std::vector< LO > * > rowIdxV;
  rowIdxV.resize( numCPUThreads );
  std::vector<std::vector< LO > * > colIdxV;
  colIdxV.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesVLap;
  valuesVLap.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesV11;
  valuesV11.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesV22;
  valuesV22.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesV33;
  valuesV33.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesV12;
  valuesV12.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesV13;
  valuesV13.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesV23;
  valuesV23.resize( numCPUThreads );

  std::vector< std::vector< LO > * > rowIdxK;
  rowIdxK.resize( numCPUThreads );
  std::vector< std::vector< LO > * > colIdxK;
  colIdxK.resize( numCPUThreads );
  std::vector< std::vector< SC > * > valuesK;
  valuesK.resize( numCPUThreads );

  for ( int i = 0; i < numCPUThreads; ++i ) {
    rowIdxV[ i ] = new std::vector< LO >;
    colIdxV[ i ] = new std::vector< LO >;
    valuesVLap[ i ] = new std::vector< SC >;
    valuesV11[ i ] = new std::vector< SC >;
    valuesV22[ i ] = new std::vector< SC >;
    valuesV33[ i ] = new std::vector< SC >;
    valuesV12[ i ] = new std::vector< SC >;
    valuesV13[ i ] = new std::vector< SC >;
    valuesV23[ i ] = new std::vector< SC >;
    rowIdxV[ i ]->reserve( 200 * nElems / numCPUThreads );
    colIdxV[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesVLap[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesV11[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesV22[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesV33[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesV12[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesV13[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesV23[ i ]->reserve( 200 * nElems / numCPUThreads );

    rowIdxK[ i ] = new std::vector< LO >;
    colIdxK[ i ] = new std::vector< LO >;
    valuesK[ i ] = new std::vector< SC >;
    rowIdxK[ i ]->reserve( 200 * nElems / numCPUThreads );
    colIdxK[ i ]->reserve( 200 * nElems / numCPUThreads );
    valuesK[ i ]->reserve( 200 * nElems / numCPUThreads );
  }





#pragma omp parallel
  {
    if ( omp_get_thread_num( ) < N_MIC ) {
 int device = omp_get_thread_num( );
    SC * myDataVLap1 = matrixBuffersVLap[ 0 ][ device ];
    SC * myDataVLap2 = matrixBuffersVLap[ 1 ][ device ];
    SC * myDataV11_1 = matrixBuffersV11[ 0 ][ device ];
    SC * myDataV11_2 = matrixBuffersV11[ 1 ][ device ];
    SC * myDataV22_1 = matrixBuffersV22[ 0 ][ device ];
    SC * myDataV22_2 = matrixBuffersV22[ 1 ][ device ];
    SC * myDataV33_1 = matrixBuffersV33[ 0 ][ device ];
    SC * myDataV33_2 = matrixBuffersV33[ 1 ][ device ];
    SC * myDataV12_1 = matrixBuffersV12[ 0 ][ device ];
    SC * myDataV12_2 = matrixBuffersV12[ 1 ][ device ];
    SC * myDataV13_1 = matrixBuffersV13[ 0 ][ device ];
    SC * myDataV13_2 = matrixBuffersV13[ 1 ][ device ];
    SC * myDataV23_1 = matrixBuffersV23[ 0 ][ device ];
    SC * myDataV23_2 = matrixBuffersV23[ 1 ][ device ];
    SC * myDataK1 = matrixBuffersK[ 0 ][ device ];
    SC * myDataK2 = matrixBuffersK[ 1 ][ device ];
    LO * myCloseElemsI_1 = closeElemsBufferI[ 0 ][ device ];
    LO * myCloseElemsI_2 = closeElemsBufferI[ 1 ][ device ];
    LO * myCloseElemsJ_1 = closeElemsBufferJ[ 0 ][ device ];
    LO * myCloseElemsJ_2 = closeElemsBufferJ[ 1 ][ device ];
    LO lengthV = bufferLengthsV[ device ];
    LO lengthK = bufferLengthsK[ device ];

#pragma offload_transfer target( mic : device ) \
    nocopy( myDataVLap1 : length( lengthV ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataVLap2 : length( lengthV ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV11_1 : length( lengthV ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV11_2 : length( lengthV ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV22_1 : length( lengthV ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV22_2 : length( lengthV ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV33_1 : length( lengthV ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV33_2 : length( lengthV ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV12_1 : length( lengthV ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV12_2 : length( lengthV ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV13_1 : length( lengthV ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV13_2 : length( lengthV ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV23_1 : length( lengthV ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataV23_2 : length( lengthV ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataK1 : length( lengthK ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myDataK2 : length( lengthK ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myCloseElemsI_1 : length( closeBuffLen ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myCloseElemsI_2 : length( closeBuffLen ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myCloseElemsJ_1 : length( closeBuffLen ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    nocopy( myCloseElemsJ_2 : length( closeBuffLen ) alloc_if( 1 ) \
    free_if( 0 ) ) \
    in( nodes : length( 3 * nNodes ) alloc_if( 1 ) free_if( 0 ) ) \
    in( elems : length( 3 * nElems ) alloc_if( 1 ) free_if( 0 ) ) \
    in( areas : length( nElems ) alloc_if( 1 ) free_if( 0 ) ) \
    in( normals : length( 3 * nElems ) alloc_if( 1 ) free_if( 0 ) ) \
    in( phi1Values : length( outerPoints * innerPoints ) \
    alloc_if(1) free_if(0)) \
    in( phi2Values : length( outerPoints * innerPoints ) \
    alloc_if(1) free_if(0)) \
    in( phi3Values : length( outerPoints * innerPoints ) \
    alloc_if(1) free_if(0)) \
    in( vOutW : length( outerPoints * innerPoints ) \
    alloc_if( 1 ) free_if( 0 ) ) \
    in( vInW : length( outerPoints * innerPoints ) \
    alloc_if( 1 ) free_if( 0 ) ) \
    in( nRowsPerSubmatrix : alloc_if( 1 ) free_if( 0 ) ) \
    in( remainingRowsPerMIC : alloc_if( 1 ) free_if( 0 ) ) \
    in( outerX1ref : length( outerPoints * innerPoints ) \
    alloc_if( 1 ) free_if( 0 ) ) \
    in( outerX2ref : length( outerPoints * innerPoints ) \
    alloc_if( 1 ) free_if( 0 ) ) \
    in( innerX1ref : length( outerPoints * innerPoints ) \
    alloc_if( 1 ) free_if( 0 ) ) \
    in( innerX2ref : length( outerPoints * innerPoints ) \
    alloc_if( 1 ) free_if( 0 ) ) \
    in( this : alloc_if( 1 ) free_if( 0 ) ) 
      // double buffering (one buffer for computation, one buffer sent to cpu)
      for ( LO s = 0; s <= nSubmatrices; ++s ) {
        if ( s != nSubmatrices ) {
          SC * computationDataVLap = matrixBuffersVLap[ s % 2 ][ device ];
          SC * computationDataV11 = matrixBuffersV11[ s % 2 ][ device ];
          SC * computationDataV22 = matrixBuffersV22[ s % 2 ][ device ];
          SC * computationDataV33 = matrixBuffersV33[ s % 2 ][ device ];
          SC * computationDataV12 = matrixBuffersV12[ s % 2 ][ device ];
          SC * computationDataV13 = matrixBuffersV13[ s % 2 ][ device ];
          SC * computationDataV23 = matrixBuffersV23[ s % 2 ][ device ];
          SC * computationDataK = matrixBuffersK[ s % 2 ][ device ];
          LO * bufferI = closeElemsBufferI[ s % 2 ][ device ];
          LO * bufferJ = closeElemsBufferJ[ s % 2 ][ device ];

      //signal( matrixBuffersVLap[ s % 2 ][ device ] ) 
#pragma offload target( mic : device ) \
      in( computationDataVLap : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( computationDataV11 : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( computationDataV22 : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( computationDataV33 : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( computationDataV12 : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( computationDataV13 : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( computationDataV23 : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( computationDataK : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( bufferI : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( bufferJ : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( nodes : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( elems : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( areas : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( normals : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( phi1Values : length( 0 ) alloc_if(0) free_if(0)) \
      in( phi2Values : length( 0 ) alloc_if(0) free_if(0)) \
      in( phi3Values : length( 0 ) alloc_if(0) free_if(0)) \
      in( vOutW : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( vInW : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( nRowsPerSubmatrix :  alloc_if( 0 ) free_if( 0 ) ) \
      in( remainingRowsPerMIC : alloc_if( 0 ) free_if( 0 ) ) \
      in( outerX1ref : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( outerX2ref : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( innerX1ref : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( innerX2ref : length( 0 ) alloc_if( 0 ) free_if( 0 ) ) \
      in( s ) \
      signal( computationDataVLap ) \
      in( this : length( 0 ) alloc_if( 0 ) free_if( 0 ) )
          {
              LO ** closePanelsI = new LO*[N_MIC_THREADS];
              LO ** closePanelsJ = new LO*[N_MIC_THREADS];
              LO * counterClose = new LO[N_MIC_THREADS];

              LO allocClose = 24*nRowsPerSubmatrix[device];
              #pragma omp parallel for schedule(dynamic) num_threads(N_MIC_THREADS)
              for (LO i = 0; i < N_MIC_THREADS; ++i) {
                counterClose[i] = 0;
                closePanelsI[i] = (LO*) _mm_malloc(allocClose * 
                    sizeof(LO), DATA_ALIGN);
                 closePanelsJ[i] = (LO*) _mm_malloc(allocClose * 
                    sizeof(LO), DATA_ALIGN);
              }
#pragma omp parallel for schedule( static, 32 ) num_threads( N_MIC_THREADS )
#pragma vector aligned(computationDataK)
            for ( LO i = 0; i < bufferLengthsK[ device ]; ++i ) {
              computationDataK[ i ] = 0.0;
            }

#pragma omp parallel num_threads( N_MIC_THREADS )
            {

              int myThreadNum = omp_get_thread_num();  
              SCVT * outerX1 = (SCVT *) _mm_malloc(
                  outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
              SCVT * outerX2 = (SCVT *) _mm_malloc(
                  outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
              SCVT * outerX3 = (SCVT *) _mm_malloc(
                  outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
              SCVT * innerX1 = (SCVT *) _mm_malloc(
                  outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
              SCVT * innerX2 = (SCVT *) _mm_malloc(
                  outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
              SCVT * innerX3 = (SCVT *) _mm_malloc(
                  outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );

              SC elemMatrixVLap = 0.0;
              SC elemMatrixV11 = 0.0;
              SC elemMatrixV22 = 0.0;
              SC elemMatrixV33 = 0.0;
              SC elemMatrixV12 = 0.0;
              SC elemMatrixV13 = 0.0;
              SC elemMatrixV23 = 0.0;
              SC elemMatrixK[3] = { 0.0, 0.0, 0.0 };
              LO startRow = device * nMICRows + s * nRowsPerSubmatrix[ device ];
              LO endRow = startRow + nRowsPerSubmatrix[ device ];
              endRow +=
                  ( s == nSubmatrices - 1 ) ? remainingRowsPerMIC[ device ] : 0;

//#pragma omp for collapse( 2 ) schedule( dynamic, 32 )

//#pragma novector
#pragma omp for schedule( dynamic, 16 )
                for ( LO j = 0; j < nElems; ++j ) {
                    
              for ( LO i = startRow; i < endRow; ++i ) {
                  if ( !BEIntegrator<LO, SC, BEIntegratorLaplace<LO, SC> >::
                      areElementsDisjointOnMIC( i, j, nodes, elems ) ) {
                      counterClose[myThreadNum]++;
                      if (counterClose[myThreadNum] < allocClose) {
                        closePanelsI[myThreadNum][counterClose[myThreadNum] - 1] = i;
                        closePanelsJ[myThreadNum][counterClose[myThreadNum] - 1] = j;
                      } else {
                        std::cout << "JE TO V PRDELI, KRETENE!" << std::endl;
                        exit(-1);
                      }
                    continue;
                  }

                  BEIntegratorLame<LO, SC >::
                      computeElemMatrixAllMIC(
                      nodes, elems, areas, normals, i, j, qOrderOuter,
                      qOrderInner, outerX1ref, outerX2ref, innerX1ref, innerX2ref,
                      outerX1, outerX2, outerX3, innerX1, innerX2, innerX3,
                      phi1Values, phi2Values, phi3Values, vOutW, vInW,
                      &elemMatrixV11, &elemMatrixV22,
                      &elemMatrixV33, &elemMatrixV12, &elemMatrixV13,
                      &elemMatrixV23, &elemMatrixVLap, elemMatrixK );

                  computationDataVLap[ j * ( endRow - startRow ) + i - startRow ] =
                      elemMatrixVLap;
                  computationDataV11[ j * ( endRow - startRow ) + i - startRow ] =
                      elemMatrixV11;
                  computationDataV22[ j * ( endRow - startRow ) + i - startRow ] =
                      elemMatrixV22;
                  computationDataV33[ j * ( endRow - startRow ) + i - startRow ] =
                      elemMatrixV33;
                  computationDataV12[ j * ( endRow - startRow ) + i - startRow ] =
                      elemMatrixV12;
                  computationDataV13[ j * ( endRow - startRow ) + i - startRow ] =
                      elemMatrixV13;
                  computationDataV23[ j * ( endRow - startRow ) + i - startRow ] =
                      elemMatrixV23;
#pragma omp atomic update
                  computationDataK[ elems[ 3 * j ] *
                      ( endRow - startRow ) + i - startRow ] += elemMatrixK[ 0 ];
#pragma omp atomic update
                  computationDataK[ elems[ 3 * j + 1 ] *
                      ( endRow - startRow ) + i - startRow ] += elemMatrixK[ 1 ];
#pragma omp atomic update
                  computationDataK[ elems[ 3 * j + 2 ] *
                      ( endRow - startRow ) + i - startRow ] += elemMatrixK[ 2 ];
              
                }
              }
              _mm_free( outerX1 );
              _mm_free( outerX2 );
              _mm_free( outerX3 );
              _mm_free( innerX1 );
              _mm_free( innerX2 );
              _mm_free( innerX3 );
            }
            LO totalClose = 0;
            for (LO i = 0 ; i < N_MIC_THREADS; ++i) {
                totalClose += counterClose[i];
            }
            numClosePanels[s] = totalClose;
            if ( totalClose < closeBuffLen ) {
                LO copiedIdx = 0;
              for (LO i = 0 ; i < N_MIC_THREADS; ++i) {
                memcpy( bufferI + copiedIdx, closePanelsI[i], counterClose[i] * sizeof(LO) );
                memcpy( bufferJ + copiedIdx, closePanelsJ[i], counterClose[i] * sizeof(LO) );
                copiedIdx += counterClose[i];
              }
            } else {
                std::cout << "JE TO V JESTE VETSI PRDELI, KRETENE!" << std::endl;
                exit(-1);
            }

              for (LO i = 0; i < N_MIC_THREADS; ++i) {
                _mm_free(closePanelsI[i]);
                _mm_free(closePanelsJ[i]);
              }
              delete [] closePanelsI;
              delete [] closePanelsJ;
              delete [] counterClose;
          }
        }
        // send data to CPU
        if ( s > 0 ) {
          SC * outputBufferVLap = matrixBuffersVLap[ ( s - 1 ) % 2 ][ device ];
          SC * outputBufferV11 = matrixBuffersV11[ ( s - 1 ) % 2 ][ device ];
          SC * outputBufferV22 = matrixBuffersV22[ ( s - 1 ) % 2 ][ device ];
          SC * outputBufferV33 = matrixBuffersV33[ ( s - 1 ) % 2 ][ device ];
          SC * outputBufferV12 = matrixBuffersV12[ ( s - 1 ) % 2 ][ device ];
          SC * outputBufferV13 = matrixBuffersV13[ ( s - 1 ) % 2 ][ device ];
          SC * outputBufferV23 = matrixBuffersV23[ ( s - 1 ) % 2 ][ device ];
          SC * outputBufferK = matrixBuffersK[ ( s - 1 ) % 2 ][ device ];
          LO * outputBufferCloseI = closeElemsBufferI[ ( s - 1 ) % 2 ][ device ];
          LO * outputBufferCloseJ = closeElemsBufferJ[ ( s - 1 ) % 2 ][ device ];
          
          SC * VLapHost = VLapData + device * nMICRows*nElems + (s-1) * nRowsPerSubmatrix[ device ] * nElems;
          SC * VLame11Host = VLame11Data + device * nMICRows*nElems + (s-1) * nRowsPerSubmatrix[ device ] * nElems;
          SC * VLame22Host = VLame22Data + device * nMICRows*nElems + (s-1) * nRowsPerSubmatrix[ device ] * nElems;
          SC * VLame33Host = VLame33Data + device * nMICRows*nElems + (s-1) * nRowsPerSubmatrix[ device ] * nElems;
          SC * VLame12Host = VLame12Data + device * nMICRows*nElems + (s-1) * nRowsPerSubmatrix[ device ] * nElems;
          SC * VLame13Host = VLame13Data + device * nMICRows*nElems + (s-1) * nRowsPerSubmatrix[ device ] * nElems;
          SC * VLame23Host = VLame23Data + device * nMICRows*nElems + (s-1) * nRowsPerSubmatrix[ device ] * nElems;
          LO rem = ( s == nSubmatrices - 1) ? remainingRowsPerMIC[ device ] : 0;
          LO bufferLength = (nRowsPerSubmatrix[ device ] + rem) * nElems; 
          
#pragma offload_transfer target( mic : device ) \
        wait(outputBufferVLap) \
        out( outputBufferVLap : length( bufferLength ) alloc_if(0) free_if(0) into(VLapHost) ) \
        out( outputBufferV11 : length( bufferLength ) alloc_if(0) free_if(0) into(VLame11Host)) \
        out( outputBufferV22 : length( bufferLength ) alloc_if(0) free_if(0) into(VLame22Host)) \
        out( outputBufferV33 : length( bufferLength ) alloc_if(0) free_if(0) into(VLame33Host)) \
        out( outputBufferV12 : length( bufferLength ) alloc_if(0) free_if(0) into(VLame12Host)) \
        out( outputBufferV13 : length( bufferLength ) alloc_if(0) free_if(0) into(VLame13Host)) \
        out( outputBufferV23 : length( bufferLength ) alloc_if(0) free_if(0) into(VLame23Host)) \
        out( outputBufferK : length( bufferLengthsK[ device ] ) alloc_if(0) free_if(0)) \
        out( outputBufferCloseI : length( numClosePanels[(s-1) % 2 ]) alloc_if(0) free_if(0) ) \
        out( outputBufferCloseJ : length( numClosePanels[(s-1) % 2 ]) alloc_if(0) free_if(0) )

         for ( LO k = 0; k < nNodes; ++k ) {
            LO rem = ( s == nSubmatrices - 1) ? remainingRowsPerMIC[ device ] : 0;
            memcpy( KData + k * nElems + device * nMICRows +
                ( s - 1 ) * nRowsPerSubmatrix[ device ],
                matrixBuffersK[ ( s - 1 ) % 2 ][ device ] +
                k * ( nRowsPerSubmatrix[ device ] + rem ),
                ( nRowsPerSubmatrix[ device ] + rem ) * sizeof ( SC ) );
          }
        std::vector<LO> currCloseI( outputBufferCloseI, outputBufferCloseI + numClosePanels[(s-1)%2] );
        std::vector<LO> currCloseJ( outputBufferCloseJ, outputBufferCloseJ + numClosePanels[(s-1)%2] );
        for (LO k = 0 ; k < numClosePanels[(s-1)%2]; ++k) {
            closeI.push_back(outputBufferCloseI[k]);
            closeJ.push_back(outputBufferCloseJ[k]);
        }
       }
      }

      #pragma offload_transfer target( mic : device ) \
      nocopy( myDataVLap1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataVLap2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV11_1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV11_2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV22_1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV22_2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV33_1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV33_2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV12_1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV12_2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV13_1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV13_2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV23_1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataV23_2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataK1 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( myDataK2 : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( nodes : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( elems : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( areas : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( vOutW : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( vInW : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( outerX1ref : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( outerX2ref : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( innerX1ref : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( innerX2ref : alloc_if( 0 ) free_if( 1 ) ) \
      nocopy( phi1Values : length( 0 ) alloc_if(0) free_if(1)) \
      nocopy( phi2Values : length( 0 ) alloc_if(0) free_if(1)) \
      nocopy( phi3Values : length( 0 ) alloc_if(0) free_if(1))

      std::cout << "MIC FINISHED" << std::endl;
    } else {
#pragma omp single nowait
      {
          std::cout << "CPU STARTED"<< std::endl;

       /* for ( LO i = 0; i < nElems; ++i ) {
#pragma omp task
          {
            int myThreadNum = omp_get_thread_num( );
            BEIntegratorLame<LO, SC> integratorV( this->space,
                this->quadratureOrder, this->quadrature,
                this->quadratureOrderDisjointElems );
            BEIntegratorLaplace<LO, SC> integratorK( &bespaceK,
                this->quadratureOrder, this->quadrature,
                this->quadratureOrderDisjointElems );
            FullMatrix< LO, SC > elemMatrixV( 1, 7, true );
            FullMatrix< LO, SC > elemMatrixK( 1, 3, true );

            for ( LO j = 0; j < nElems; ++j ) {
              if ( BEIntegrator<LO, SC, BEIntegratorLaplace<LO, SC> >::
                  areElementsDisjointOnMIC( j, i, nodes, elems ) ) {
                continue;
              }
              integratorV.getElemMatrix1Layer( i, j, elemMatrixV );
              integratorK.getElemMatrix2Layer( i, j, elemMatrixK );

              rowIdxV[ myThreadNum - N_MIC ]->push_back( i );
              colIdxV[ myThreadNum - N_MIC ]->push_back( j );
              valuesVLap[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixV.get( 0, 0 ) );
              valuesV11[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixV.get( 0, 1 ) );
              valuesV22[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixV.get( 0, 2 ) );
              valuesV33[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixV.get( 0, 3 ) );
              valuesV12[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixV.get( 0, 4 ) );
              valuesV13[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixV.get( 0, 5 ) );
              valuesV23[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixV.get( 0, 6 ) );

              rowIdxK[ myThreadNum - N_MIC ]->push_back( i );
              rowIdxK[ myThreadNum - N_MIC ]->push_back( i );
              rowIdxK[ myThreadNum - N_MIC ]->push_back( i );
              colIdxK[ myThreadNum - N_MIC ]->push_back( elems[ 3 * j ] );
              colIdxK[ myThreadNum - N_MIC ]->push_back( elems[ 3 * j + 1 ] );
              colIdxK[ myThreadNum - N_MIC ]->push_back( elems[ 3 * j + 2 ] );
              valuesK[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixK.get( 0, 0 ) );
              valuesK[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixK.get( 0, 1 ) );
              valuesK[ myThreadNum - N_MIC ]
                  ->push_back( elemMatrixK.get( 0, 2 ) );
            } // end for j elems
          } // end omp task
        }*/ // end i elems
        LO startRow = N_MIC * nMICRows;
std::cout << "start" <<std::endl;
        for ( LO i = startRow; i < nElems; ++i ) {
#pragma omp task
          {

            SCVT * outerX1 = (SCVT *) _mm_malloc(
                outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
            SCVT * outerX2 = (SCVT *) _mm_malloc(
                outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
            SCVT * outerX3 = (SCVT *) _mm_malloc(
                outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
            SCVT * innerX1 = (SCVT *) _mm_malloc(
                outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
            SCVT * innerX2 = (SCVT *) _mm_malloc(
                outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
            SCVT * innerX3 = (SCVT *) _mm_malloc(
                outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );

            SC elemMatrixVLap = 0.0;
            SC elemMatrixV11 = 0.0;
            SC elemMatrixV22 = 0.0;
            SC elemMatrixV33 = 0.0;
            SC elemMatrixV12 = 0.0;
            SC elemMatrixV13 = 0.0;
            SC elemMatrixV23 = 0.0;
            SC elemMatrixK[3] = { 0.0, 0.0, 0.0 };



        for ( LO j = 0; j < nElems; ++j ) {
              if ( !BEIntegrator<LO, SC, BEIntegratorLaplace<LO, SC> >::
                  areElementsDisjointOnMIC( i, j, nodes, elems ) ) {
                continue;
              }

              BEIntegratorLame<LO, SC >::
                  computeElemMatrixAllMIC(
                  nodes, elems, areas, normals, i, j, qOrderOuter,
                  qOrderInner, outerX1ref, outerX2ref, innerX1ref, innerX2ref,
                  outerX1, outerX2, outerX3, innerX1, innerX2, innerX3,
                  phi1Values, phi2Values, phi3Values, vOutW, vInW,
                  &elemMatrixV11, &elemMatrixV22,
                  &elemMatrixV33, &elemMatrixV12, &elemMatrixV13,
                  &elemMatrixV23, &elemMatrixVLap, elemMatrixK );

              VLap.set( j, i, elemMatrixVLap );

//#pragma omp atomic update
              KData[ elems[ 3 * j ] *
                  nElems + i ] += elemMatrixK[0];
//#pragma omp atomic update
             KData[ elems[ 3 * j + 1 ] *
                  nElems + i ] += elemMatrixK[1];
//#pragma omp atomic update
              KData[ elems[ 3 * j + 2 ] *
                  nElems + i ] += elemMatrixK[2];

              VLame11->set(j,i, elemMatrixV11);
              VLame22->set(j,i, elemMatrixV22);
              VLame33->set(j,i, elemMatrixV33);
              VLame12->set(j,i, elemMatrixV12);
              VLame13->set(j,i, elemMatrixV13);
              VLame23->set(j,i, elemMatrixV23);
            }

            _mm_free( outerX1 );
            _mm_free( outerX2 );
            _mm_free( outerX3 );
            _mm_free( innerX1 );
            _mm_free( innerX2 );
            _mm_free( innerX3 );
          } // end omp task
        } // end i elems
        std::cout << "CPU  FINISHED" << std::endl;
      } // end omp single nowait
    } // end else
#pragma omp barrier
  } // end omp parallel
#pragma omp parallel num_threads( numCPUThreads )
  {
    int myThreadNum = omp_get_thread_num( );
    for ( LO j = 0; j < rowIdxV[ myThreadNum]->size( ); ++j ) {
      VLap.set( rowIdxV[ myThreadNum ]->at( j ), colIdxV[ myThreadNum ]->at( j ),
          valuesVLap[ myThreadNum ]->at( j ) );
      VLame11->set( rowIdxV[ myThreadNum ]->at( j ), colIdxV[ myThreadNum ]->at( j ),
          valuesV11[ myThreadNum ]->at( j ) );
      VLame22->set( rowIdxV[ myThreadNum ]->at( j ), colIdxV[ myThreadNum ]->at( j ),
          valuesV22[ myThreadNum ]->at( j ) );
      VLame33->set( rowIdxV[ myThreadNum ]->at( j ), colIdxV[ myThreadNum ]->at( j ),
          valuesV33[ myThreadNum ]->at( j ) );
      VLame12->set( rowIdxV[ myThreadNum ]->at( j ), colIdxV[ myThreadNum ]->at( j ),
          valuesV12[ myThreadNum ]->at( j ) );
      VLame13->set( rowIdxV[ myThreadNum ]->at( j ), colIdxV[ myThreadNum ]->at( j ),
          valuesV13[ myThreadNum ]->at( j ) );
      VLame23->set( rowIdxV[ myThreadNum ]->at( j ), colIdxV[ myThreadNum ]->at( j ),
          valuesV23[ myThreadNum ]->at( j ) );
    }
    for ( LO j = 0; j < rowIdxK[myThreadNum]->size( ); ++j ) {
#pragma omp atomic update
      KData[ colIdxK[ myThreadNum ]->at( j ) * nElems +
          rowIdxK[ myThreadNum ]->at( j ) ] += valuesK[ myThreadNum ]->at( j );
    }
  }

  // add the Laplace matrix V to the diagonal of the result
  SC val = 0.0;
  LO i = 0;
#pragma omp parallel for private(val)
  for ( LO j = 0; j < nColsScalar; ++j ) {

#pragma omp simd \
linear( i : 1 ) \
private( val ) \
simdlen( DATA_WIDTH )
    for ( i = 0; i < nRowsScalar; ++i ) {
      val = VLapData[ j * nRowsScalar + i ] * ( 3.0 - 4.0 * this->nu ); 
      VLame11Data[ j * nRowsScalar + i ] += val;
      VLame22Data[ j * nRowsScalar + i ] += val;
      VLame33Data[ j * nRowsScalar + i ] += val;
    }
  }
  
  VLame11->scale( ( 1.0 + this->nu ) / ( 2.0 * E * ( 1.0 - this->nu ) ) );
  VLame22->scale( ( 1.0 + this->nu ) / ( 2.0 * E * ( 1.0 - this->nu ) ) );
  VLame33->scale( ( 1.0 + this->nu ) / ( 2.0 * E * ( 1.0 - this->nu ) ) );
  VLame12->scale( ( 1.0 + this->nu ) / ( 2.0 * E * ( 1.0 - this->nu ) ) );
  VLame13->scale( ( 1.0 + this->nu ) / ( 2.0 * E * ( 1.0 - this->nu ) ) );
  VLame23->scale( ( 1.0 + this->nu ) / ( 2.0 * E * ( 1.0 - this->nu ) ) );
  VLame.setBlock( 0, 0, VLame11, true );
  VLame.setBlock( 1, 1, VLame22, true );
  VLame.setBlock( 2, 2, VLame33, true );
  VLame.setBlock( 1, 0, VLame12, true );
  VLame.setBlock( 2, 0, VLame13, true );
  VLame.setBlock( 2, 1, VLame23, true );
  VLame.setBlock( 0, 1, VLame13, true );
  VLame.setBlock( 0, 2, VLame13, true );
  VLame.setBlock( 1, 2, VLame23, true );

  /*SCVT mult = ( (SCVT) 1.0 + this->nu ) /
      ( (SCVT) 2.0 * E * ( (SCVT) 1.0 - this->nu ) );
      
#pragma omp parallel for simd 
  for ( LO i = 0; i < nRows * nRows; ++i ) {
    VLameData[ i ] *= mult;
  }*/
  for ( LO i = 0; i < N_MIC; ++i ) {
    _mm_free( matrixBuffersVLap[ 0 ][ i ] );
    _mm_free( matrixBuffersVLap[ 1 ][ i ] );
    _mm_free( matrixBuffersV11[ 0 ][ i ] );
    _mm_free( matrixBuffersV11[ 1 ][ i ] );
    _mm_free( matrixBuffersV22[ 0 ][ i ] );
    _mm_free( matrixBuffersV22[ 1 ][ i ] );
    _mm_free( matrixBuffersV33[ 0 ][ i ] );
    _mm_free( matrixBuffersV33[ 1 ][ i ] );
    _mm_free( matrixBuffersV12[ 0 ][ i ] );
    _mm_free( matrixBuffersV12[ 1 ][ i ] );
    _mm_free( matrixBuffersV13[ 0 ][ i ] );
    _mm_free( matrixBuffersV13[ 1 ][ i ] );
    _mm_free( matrixBuffersV23[ 0 ][ i ] );
    _mm_free( matrixBuffersV23[ 1 ][ i ] );
    _mm_free( matrixBuffersK[ 0 ][ i ] );
    _mm_free( matrixBuffersK[ 1 ][ i ] );
  }

  _mm_free( vOutW );
  _mm_free( vInW );
  _mm_free( outerX1ref );
  _mm_free( outerX2ref );
  _mm_free( innerX1ref );
  _mm_free( innerX2ref );
  _mm_free( phi1Values );
  _mm_free( phi2Values );
  _mm_free( phi3Values );

  for ( int i = 0; i < numCPUThreads; ++i ) {
    delete rowIdxV[ i ];
    delete colIdxV[ i ];
    delete valuesVLap[ i ];
    delete valuesV11[ i ];
    delete valuesV22[ i ];
    delete valuesV33[ i ];
    delete valuesV12[ i ];
    delete valuesV13[ i ];
    delete valuesV23[ i ];
    delete rowIdxK[ i ];
    delete colIdxK[ i ];
    delete valuesK[ i ];
  }
#endif

}

}

#endif
