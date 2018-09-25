/*!
 * @file    BEBilinearFormWave1Layer.cpp
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    December 10, 2013
 * 
 */


#ifdef BEBILINEARFORMWAVE1LAYER_H

namespace bem4i {

template<class LO, class SC>
BEBilinearFormWave1Layer<LO, SC>::BEBilinearFormWave1Layer(
    ) {
}

template<class LO, class SC>
BEBilinearFormWave1Layer<LO, SC>::BEBilinearFormWave1Layer(
    const BEBilinearFormWave1Layer& orig
    ) {
}

template<class LO, class SC>
BEBilinearFormWave1Layer<LO, SC>::BEBilinearFormWave1Layer(
    BESpace<LO, SC>* space,
    int* quadratureOrder,
    int timeQuadratureOrder,
    int* quadratureOrderDisjointElems
    ) {
  this->space = space;

  if ( quadratureOrder ) {
    this->quadratureOrder = quadratureOrder;
  } else {
    this->quadratureOrder = defaultQuadraturesSauterSchwab;
  }

  this->timeQuadOrder = timeQuadratureOrder;

  this->quadratureOrderDisjointElems = quadratureOrderDisjointElems;
}

template<class LO, class SC>
BEBilinearFormWave1Layer<LO, SC>::~BEBilinearFormWave1Layer( ) {
}

template<class LO, class SC>
void BEBilinearFormWave1Layer<LO, SC>::assemble(
    FullMatrix<LO, SC>& matrix
    ) const {

  BESpaceTime<LO, SC> *spaceTime =
      static_cast<BESpaceTime<LO, SC>*> ( this->space );
  BEIntegratorWave<LO, SC> integrator(
      this->space, this->quadratureOrder, this->timeQuadOrder );
  const vector<LO> innerElems = this->space->getInnerElems( );
  const vector<LO> outerElems = this->space->getOuterElems( );

  // allocate local element matrices
  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  // max order of Legendre polynomials in temporal basis
  int legOrder = spaceTime->getLegendreOrder( );

  int nTimeSteps = spaceTime->getNTimeSteps( );

  FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );

  matrix.resize( ( nTimeSteps ) * ( legOrder + 1 ) *
      this->space->getOuterDOFs( ),
      ( nTimeSteps ) * ( legOrder + 1 ) * this->space->getInnerDOFs( ) );

  // vector of indices where to put values in the global matrix
  vector<LO> rowIndices( nLocalRows );
  vector<LO> colIndices( nLocalCols );

  // outer two loops over time blocks
  for ( int i = 0; i < nTimeSteps; i++ ) {
    for ( int j = 0; j < nTimeSteps; j++ ) {
      integrator.setCurrentFunctions( i, j );

      // two loops over Legendre polynomial order
      for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
        for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {
          integrator.setCurrentLegendreOrder( legBasis, legTest );

          if ( i <= j + 1 ) {
            // inner two loops over space
            for ( auto itOut = outerElems.begin( ); itOut != outerElems.end( );
                itOut++ ) {
              for ( auto itIn = innerElems.begin( ); itIn != innerElems.end( );
                  itIn++ ) {
                integrator.getElemMatrix1Layer( *itOut, *itIn, elemMatrix );
                this->space->getOuterElemDOFs( *itOut, &rowIndices[0] );
                this->space->getInnerElemDOFs( *itIn, &colIndices[0] );
                for ( int k = 0; k < nLocalRows; k++ ) {
                  rowIndices[k] += ( j * ( legOrder + 1 ) + legTest ) *
                      this->space->getOuterDOFs( );
                }
                for ( int k = 0; k < nLocalCols; k++ ) {
                  colIndices[k] += ( i * ( legOrder + 1 ) + legBasis ) *
                      this->space->getInnerDOFs( );
                }
                matrix.addToPositions( rowIndices, colIndices, elemMatrix );
              }
            }
          }
        }
      }
    }
  }
}

template<class LO, class SC>
void BEBilinearFormWave1Layer<LO, SC>::assemble(
    SparseMatrix<LO, SC>& matrix
    ) const {

  //    BESpaceTime<LO, SC> *spaceTime = static_cast<BESpaceTime<LO, SC>*> ( this->space );
  //    BEIntegratorWave<LO, SC> integrator( this->space, this->quadratureOrder, this->timeQuadOrder );
  //    const vector<LO> innerElems = this->space->getInnerElems( );
  //    const vector<LO> outerElems = this->space->getOuterElems( );
  //  
  //    // allocate local element matrices
  //    LO nLocalRows = this->space->getDOFsPerOuterElem( );
  //    LO nLocalCols = this->space->getDOFsPerInnerElem( );
  //  
  //    int nTimeSteps = spaceTime->getNTimeSteps( );
  //  
  //    FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );
  //    
  //    matrix.resize( nTimeSteps * this->space->getOuterDOFs( ),
  //        nTimeSteps * this->space->getInnerDOFs( ) );
  //    matrix.getEigenSparseMatrix()->uncompress();
  //    
  //    Eigen::VectorXi resSizes( nTimeSteps * this->space->getInnerDOFs( ) );
  //    for ( int i = 0; i < nTimeSteps; i++ ) {
  //      for ( int j = 0; j < this->space->getInnerDOFs( ); j++ ) {
  //        resSizes( i * this->space->getInnerDOFs( ) + j ) = ( nTimeSteps -i) * this->space->getInnerDOFs( );
  //  
  //     }
  //    }
  //  
  //    matrix.getEigenSparseMatrix( )->reserve( resSizes );
  //    // vector of indices where to put values in the global matrix
  //    vector<LO> rowIndices( nLocalRows );
  //    vector<LO> colIndices( nLocalCols );
  //  
  //  
  //    // outer two loops over time blocks
  //    for ( int i = 0; i < nTimeSteps; i++ ) { // basis
  //      for ( int j = 0; j < nTimeSteps; j++ ) { // test
  //        std::cout << " basis: " << i << ", test: " << j << std::endl;
  //        if ( i <= j + 1 ) {
  //          integrator.setCurrentFunctions( i, j );
  //  
  //          // inner two loops over space
  //          for ( auto itOut = outerElems.begin( ); itOut != outerElems.end( ); itOut++ ) {
  //            for ( auto itIn = innerElems.begin( ); itIn != innerElems.end( ); itIn++ ) {
  //              integrator.getElemMatrix1Layer( *itOut, *itIn, elemMatrix );
  //              this->space->getOuterElemDOFs( *itOut, &rowIndices[0] );
  //              this->space->getInnerElemDOFs( *itIn, &colIndices[0] );
  //              for ( int k = 0; k < nLocalRows; k++ ) {
  //                rowIndices[k] += j * this->space->getOuterDOFs( );
  //              }
  //              for ( int k = 0; k < nLocalCols; k++ ) {
  //                colIndices[k] += i * this->space->getInnerDOFs( );
  //              }
  //              matrix.addToPositions( rowIndices, colIndices, elemMatrix );
  //            }
  //          }
  //        }
  //      }
  //    }
  //    matrix.makeCompressed( );


  BESpaceTime<LO, SC> *spaceTime =
      static_cast<BESpaceTime<LO, SC>*> ( this->space );
  BEIntegratorWave<LO, SC> integrator( this->space, this->quadratureOrder,
      this->timeQuadOrder );
  const vector<LO> innerElems = this->space->getInnerElems( );
  const vector<LO> outerElems = this->space->getOuterElems( );
  int offset = 0;

  // max order of Legendre polynomials in temporal basis
  int legOrder = spaceTime->getLegendreOrder( );

  // allocate local element matrices
  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  int nTimeSteps = spaceTime->getNTimeSteps( );

  FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );

  matrix.resize( ( nTimeSteps ) * ( legOrder + 1 ) *
      this->space->getOuterDOFs( ),
      ( nTimeSteps ) * ( legOrder + 1 ) * this->space->getInnerDOFs( ) );

  // vector of indices where to put values in the global matrix
  vector<LO> rowIndices( nLocalRows );
  vector<LO> colIndices( nLocalCols );

  // create Eigen triplets and insert values into matrix (duplicated indices are summed!)
  std::vector<Eigen::Triplet<SC, LO> > tripletList;
  tripletList.reserve( ( nTimeSteps + legOrder + 1 ) *
      this->space->getInnerDOFs( ) * this->space->getInnerDOFs( ) );

  // now we have to divide computation of the matrix into several parts - 
  // computation of the first column, second column except its first and last block is copied to the rest of the matrix, 
  // the last column and the last row

  // first column
  for ( int j = 0; j < nTimeSteps; j++ ) { // test
    //std::cout << " basis: 0" << ", test: " << j << std::endl;
    integrator.setCurrentFunctions( 0, j );

    // two loops over Legendre polynomials in basis
    for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
      for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {
        integrator.setCurrentLegendreOrder( legBasis, legTest );

        // inner two loops over space
        for ( auto itOut = outerElems.begin( ); itOut != outerElems.end( );
            itOut++ ) {
          for ( auto itIn = innerElems.begin( ); itIn != innerElems.end( );
              itIn++ ) {
            integrator.getElemMatrix1Layer( *itOut, *itIn, elemMatrix );
            this->space->getOuterElemDOFs( *itOut, &rowIndices[0] );
            this->space->getInnerElemDOFs( *itIn, &colIndices[0] );
            for ( int k = 0; k < nLocalRows; k++ ) {
              rowIndices[k] += ( j * ( legOrder + 1 ) + legTest ) *
                  this->space->getOuterDOFs( );
            }
            for ( int k = 0; k < nLocalCols; k++ ) {
              colIndices[k] += legBasis * this->space->getInnerDOFs( );
            }
            //matrix.addToPositions( rowIndices, colIndices, elemMatrix );
            for ( int k = 0; k < nLocalRows; k++ ) {
              for ( int l = 0; l < nLocalCols; l++ ) {
                if ( fabs( elemMatrix.get( k, l ) ) > EPS ) {
                  tripletList.push_back( Eigen::Triplet<SC, LO>(
                      rowIndices[k], colIndices[l], elemMatrix.get( k, l ) ) );
                }
              }
            }
          }
        }
      }
    }
  }

  // second column copied to the rest of the matrix (except the last column and row)
  for ( int j = 0; j < nTimeSteps - 1; j++ ) { // test
    //std::cout << " basis: 1" << ", test: " << j << std::endl;
    integrator.setCurrentFunctions( 1, j );

    if ( j == 0 ) {
      // do not copy block on position (0, 1)
      offset = nTimeSteps - 2;
    } else {
      offset = 0;
    }

    // two loops over Legendre polynomials in basis
    for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
      for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {
        integrator.setCurrentLegendreOrder( legBasis, legTest );

        // inner two loops over space
        for ( auto itOut = outerElems.begin( ); itOut != outerElems.end( );
            itOut++ ) {
          for ( auto itIn = innerElems.begin( ); itIn != innerElems.end( );
              itIn++ ) {
            integrator.getElemMatrix1Layer( *itOut, *itIn, elemMatrix );
            this->space->getOuterElemDOFs( *itOut, &rowIndices[0] );
            this->space->getInnerElemDOFs( *itIn, &colIndices[0] );
            for ( int k = 0; k < nLocalRows; k++ ) {
              rowIndices[k] += ( j * ( legOrder + 1 ) + legTest ) *
                  this->space->getOuterDOFs( );
            }
            for ( int k = 0; k < nLocalCols; k++ ) {
              colIndices[k] += ( legOrder + 1 + legBasis ) *
                  this->space->getInnerDOFs( );
            }
            //matrix.addToPositions( rowIndices, colIndices, elemMatrix );
            for ( int t = 0; t + j + offset < nTimeSteps - 1; t++ ) {
              for ( int k = 0; k < nLocalRows; k++ ) {
                for ( int l = 0; l < nLocalCols; l++ ) {
                  if ( fabs( elemMatrix.get( k, l ) ) > EPS ) {
                    tripletList.push_back( Eigen::Triplet<SC, LO>(
                        rowIndices[k] + t * ( legOrder + 1 ) *
                        this->space->getOuterDOFs( ),
                        colIndices[l] + t * ( legOrder + 1 ) *
                        this->space->getInnerDOFs( ),
                        elemMatrix.get( k, l ) ) );
                  }
                }
              }
            }
          }
        }
      }
    }
  }



  // upper diagonal blocks - starting in the third column
  //std::cout << " basis: 1" << ", test: 2" << std::endl;
  integrator.setCurrentFunctions( 2, 1 );

  // two loops over Legendre polynomials in basis
  for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
    for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {
      integrator.setCurrentLegendreOrder( legBasis, legTest );

      // inner two loops over space
      for ( auto itOut = outerElems.begin( ); itOut != outerElems.end( );
          itOut++ ) {
        for ( auto itIn = innerElems.begin( ); itIn != innerElems.end( );
            itIn++ ) {
          integrator.getElemMatrix1Layer( *itOut, *itIn, elemMatrix );
          this->space->getOuterElemDOFs( *itOut, &rowIndices[0] );
          this->space->getInnerElemDOFs( *itIn, &colIndices[0] );
          for ( int k = 0; k < nLocalRows; k++ ) {
            rowIndices[k] += ( legOrder + 1 + legTest ) *
                this->space->getOuterDOFs( );
          }
          for ( int k = 0; k < nLocalCols; k++ ) {
            colIndices[k] += ( 2 * ( legOrder + 1 ) + legBasis ) *
                this->space->getInnerDOFs( );
          }
          //matrix.addToPositions( rowIndices, colIndices, elemMatrix );
          for ( int t = 0; t + 2 < nTimeSteps - 1; t++ ) {
            for ( int k = 0; k < nLocalRows; k++ ) {
              for ( int l = 0; l < nLocalCols; l++ ) {
                if ( fabs( elemMatrix.get( k, l ) ) > EPS ) {
                  tripletList.push_back( Eigen::Triplet<SC, LO>( rowIndices[k]
                      + t * ( legOrder + 1 ) * this->space->getOuterDOFs( ),
                      colIndices[l] + t * ( legOrder + 1 ) *
                      this->space->getInnerDOFs( ), elemMatrix.get( k, l ) ) );
                }
              }
            }
          }
        }
      }
    }
  }


  // last column
  for ( int j = nTimeSteps - 2; j < nTimeSteps; j++ ) { // test
    //std::cout << " basis: " << nTimeSteps - 1 << ", test: " << j << std::endl;
    integrator.setCurrentFunctions( nTimeSteps - 1, j );

    // two loops over Legendre polynomials in basis
    for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
      for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {
        integrator.setCurrentLegendreOrder( legBasis, legTest );

        // inner two loops over space
        for ( auto itOut = outerElems.begin( ); itOut != outerElems.end( );
            itOut++ ) {
          for ( auto itIn = innerElems.begin( ); itIn != innerElems.end( );
              itIn++ ) {
            integrator.getElemMatrix1Layer( *itOut, *itIn, elemMatrix );
            this->space->getOuterElemDOFs( *itOut, &rowIndices[0] );
            this->space->getInnerElemDOFs( *itIn, &colIndices[0] );
            for ( int k = 0; k < nLocalRows; k++ ) {
              rowIndices[k] += ( j * ( legOrder + 1 ) + legTest ) *
                  this->space->getOuterDOFs( );
            }
            for ( int k = 0; k < nLocalCols; k++ ) {
              colIndices[k] += ( ( nTimeSteps - 1 )*( legOrder + 1 ) +
                  legBasis ) * this->space->getInnerDOFs( );
            }
            //matrix.addToPositions( rowIndices, colIndices, elemMatrix );
            for ( int k = 0; k < nLocalRows; k++ ) {
              for ( int l = 0; l < nLocalCols; l++ ) {
                if ( fabs( elemMatrix.get( k, l ) ) > EPS ) {
                  tripletList.push_back( Eigen::Triplet<SC, LO>( rowIndices[k],
                      colIndices[l], elemMatrix.get( k, l ) ) );
                }
              }
            }
          }
        }
      }
    }
  }

  // last row
  for ( int i = 1; i < nTimeSteps - 1; i++ ) { // basis
    //std::cout << " basis: " << i << ", test: " << nTimeSteps - 1 << std::endl;
    integrator.setCurrentFunctions( i, nTimeSteps - 1 );

    // two loops over Legendre polynomials in basis
    for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
      for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {
        integrator.setCurrentLegendreOrder( legBasis, legTest );

        // inner two loops over space
        for ( auto itOut = outerElems.begin( ); itOut != outerElems.end( );
            itOut++ ) {
          for ( auto itIn = innerElems.begin( ); itIn != innerElems.end( );
              itIn++ ) {
            integrator.getElemMatrix1Layer( *itOut, *itIn, elemMatrix );
            this->space->getOuterElemDOFs( *itOut, &rowIndices[0] );
            this->space->getInnerElemDOFs( *itIn, &colIndices[0] );
            for ( int k = 0; k < nLocalRows; k++ ) {
              rowIndices[k] += ( ( nTimeSteps - 1 )*( legOrder + 1 ) + legTest )
                  * this->space->getOuterDOFs( );
            }
            for ( int k = 0; k < nLocalCols; k++ ) {
              colIndices[k] += ( i * ( legOrder + 1 ) + legBasis ) *
                  this->space->getInnerDOFs( );
            }
            //matrix.addToPositions( rowIndices, colIndices, elemMatrix );
            for ( int k = 0; k < nLocalRows; k++ ) {
              for ( int l = 0; l < nLocalCols; l++ ) {
                if ( fabs( elemMatrix.get( k, l ) ) > EPS ) {
                  tripletList.push_back( Eigen::Triplet<SC, LO>( rowIndices[k],
                      colIndices[l], elemMatrix.get( k, l ) ) );
                }
              }
            }

          }
        }
      }
    }
  }

  matrix.getEigenSparseMatrix( )->setFromTriplets( tripletList.begin( ),
      tripletList.end( ) );
  matrix.makeCompressed( );



}

template<class LO, class SC>
void BEBilinearFormWave1Layer<LO, SC>::assemble(
    BlockMatrix<LO, SC>& matrix
    ) const {

  // NOT WORKING !
  std::cout << "WARNING: The method assemble(BlockMatrix) is not yet functional" << std::endl;

  BESpaceTime<LO, SC> *spaceTime =
      static_cast<BESpaceTime<LO, SC>*> ( this->space );
  BEIntegratorWave<LO, SC> integrator( this->space, this->quadratureOrder,
      this->timeQuadOrder );
  const vector<LO> innerElems = this->space->getInnerElems( );
  const vector<LO> outerElems = this->space->getOuterElems( );

  SparseMatrix<LO, SC> *block;

  // max order of Legendre polynomials in temporal basis
  int legOrder = spaceTime->getLegendreOrder( );

  // allocate local element matrices
  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  int nTimeSteps = spaceTime->getNTimeSteps( );

  FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );

  LO *numbersOfRows = new LO[nTimeSteps];
  LO *numbersOfCols = new LO[nTimeSteps];

  for ( LO i = 0; i < nTimeSteps; i++ ) {
    numbersOfRows[i] = ( legOrder + 1 ) * this->space->getOuterDOFs( );
    numbersOfCols[i] = ( legOrder + 1 ) * this->space->getOuterDOFs( );
  }

  //LO blockRows = this->space->getOuterDOFs( ) * ( legOrder + 1 );
  //LO blockCols = this->space->getInnerDOFs( ) * ( legOrder + 1 );

  matrix.resize( nTimeSteps, nTimeSteps, numbersOfRows, numbersOfCols );

  // vector of indices where to put values in the global matrix
  vector<LO> rowIndices( nLocalRows );
  vector<LO> colIndices( nLocalCols );

  delete [] numbersOfCols;
  delete [] numbersOfRows;

  // at first, assemble the first column
  for ( int i = 0; i < nTimeSteps; i++ ) {

    block = new SparseMatrix<LO, SC>( ( legOrder + 1 ) *
        this->space->getOuterDOFs( ),
        ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
    assembleBlock( i, 0, *block );
    block->makeCompressed( );
    matrix.setBlock( i, 0, block, true );
  }

  // block on position (0,1) is single
  block = new SparseMatrix<LO, SC>( ( legOrder + 1 ) *
      this->space->getOuterDOFs( ),
      ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
  assembleBlock( 0, 1, *block );
  block->makeCompressed( );
  matrix.setBlock( 0, 1, block, true );

  // blocks in the second column 
  for ( int i = 1; i < nTimeSteps - 1; i++ ) {
    block = new SparseMatrix<LO, SC>( ( legOrder + 1 ) *
        this->space->getOuterDOFs( ),
        ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
    assembleBlock( i, 1, *block );
    block->makeCompressed( );
    matrix.setBlock( i, 1, block, true );
  }

  // the last row
  for ( int i = 1; i < nTimeSteps - 1; i++ ) {
    block = new SparseMatrix<LO, SC>( ( legOrder + 1 ) *
        this->space->getOuterDOFs( ),
        ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
    assembleBlock( nTimeSteps - 1, i, *block );
    block->makeCompressed( );
    matrix.setBlock( nTimeSteps - 1, i, block, true );
  }

  // the upper diagonal
  block = new SparseMatrix<LO, SC>( ( legOrder + 1 ) *
      this->space->getOuterDOFs( ),
      ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
  assembleBlock( 1, 2, *block );
  matrix.setBlock( 1, 2, block, true );
  for ( int j = 2; j < nTimeSteps - 1; j++ ) {
    matrix.setBlock( j, j + 1, block, false );
  }

  // copy pointers to the second column blocks on the appropriate positions
  for ( int i = 1; i < nTimeSteps - 2; i++ ) {
    for ( int j = 1; j < nTimeSteps - i - 1; j++ ) {
      matrix.setBlock( i + j, 1 + j, matrix.getBlock( i, 1 ), false );
    }
  }
  // the last column
  block = new SparseMatrix<LO, SC>( ( legOrder + 1 ) *
      this->space->getOuterDOFs( ),
      ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
  assembleBlock( nTimeSteps - 2, nTimeSteps - 1, *block );
  block->makeCompressed( );
  matrix.setBlock( nTimeSteps - 2, nTimeSteps - 1, block, true );

  block = new SparseMatrix<LO, SC>( ( legOrder + 1 ) *
      this->space->getOuterDOFs( ),
      ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
  assembleBlock( nTimeSteps - 1, nTimeSteps - 1, *block );
  block->makeCompressed( );
  matrix.setBlock( nTimeSteps - 1, nTimeSteps - 1, block, true );
  //  
  //  // we need to assemble only the first two and the last column
  //  // outer two loops over time blocks
  //  for ( int i = 0; i < 2; i++ ) { // basis
  //    for ( int j = 0; j < nTimeSteps; j++ ) { // test
  //      SparseMatrix<LO, SC> *block = new SparseMatrix<LO, SC>( blockRows,
  //          blockCols );
  //      integrator.setCurrentFunctions( i, j );
  //
  //      // two loops over Legendre polynomials in basis
  //      for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
  //        for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {
  //          integrator.setCurrentLegendreOrder( legBasis, legTest );
  //
  //          // inner two loops over space
  //          for ( auto itOut = outerElems.begin( ); itOut != outerElems.end( ); itOut++ ) {
  //            for ( auto itIn = innerElems.begin( ); itIn != innerElems.end( ); itIn++ ) {
  //              integrator.getElemMatrix1Layer( *itOut, *itIn, elemMatrix );
  //              this->space->getOuterElemDOFs( *itOut, &rowIndices[0] );
  //              this->space->getInnerElemDOFs( *itIn, &colIndices[0] );
  //              for ( int k = 0; k < nLocalRows; k++ ) {
  //                rowIndices[k] += legTest * this->space->getOuterDOFs( );
  //              }
  //              for ( int k = 0; k < nLocalCols; k++ ) {
  //                colIndices[k] += legBasis * this->space->getInnerDOFs( );
  //              }
  //              block->addToPositions( rowIndices, colIndices, elemMatrix );
  //            }
  //          }
  //        }
  //      }
  //      matrix.setBlock( j, i, block, true );
  //      block->makeCompressed( );
  //    }
  //  }
  //
  //
  //
  //
  //  // copy pointers to the second column blocks on the appropriate positions
  //  for ( int i = 1; i < nTimeSteps - 2; i++ ) {
  //    for ( int j = 1; j < nTimeSteps - i - 1; j++ ) {
  //      matrix.setBlock( i + j, 1 + j, matrix.getBlock( i, 1 ), false );
  //    }
  //  }
  //
  //
  //  // the last row
  //  for ( int i = 2; i < nTimeSteps - 1; i++ ) { // basis
  //
  //    SparseMatrix<LO, SC> *block = new SparseMatrix<LO, SC>( blockRows,
  //        blockCols );
  //
  //    integrator.setCurrentFunctions( i, nTimeSteps - 1 );
  //
  //    // two loops over Legendre polynomials in basis
  //    for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
  //      for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {
  //        integrator.setCurrentLegendreOrder( legBasis, legTest );
  //
  //        // inner two loops over space
  //        for ( auto itOut = outerElems.begin( ); itOut != outerElems.end( ); itOut++ ) {
  //          for ( auto itIn = innerElems.begin( ); itIn != innerElems.end( ); itIn++ ) {
  //            integrator.getElemMatrix1Layer( *itOut, *itIn, elemMatrix );
  //            this->space->getOuterElemDOFs( *itOut, &rowIndices[0] );
  //            this->space->getInnerElemDOFs( *itIn, &colIndices[0] );
  //            for ( int k = 0; k < nLocalRows; k++ ) {
  //              rowIndices[k] += legTest * this->space->getOuterDOFs( );
  //            }
  //            for ( int k = 0; k < nLocalCols; k++ ) {
  //              colIndices[k] += legBasis * this->space->getInnerDOFs( );
  //            }
  //            block->addToPositions( rowIndices, colIndices, elemMatrix );
  //          }
  //        }
  //      }
  //    }
  //    matrix.setBlock( nTimeSteps - 1, i, block, true );
  //    block->makeCompressed( );
  //
  //  }
  //
  //  // the upper diagonal
  //  SparseMatrix<LO, SC> *block = new SparseMatrix<LO, SC>( blockRows,
  //      blockCols );
  //  integrator.setCurrentFunctions( 2, 1 );
  //
  //  // two loops over Legendre polynomials in basis
  //  for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
  //    for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {
  //      integrator.setCurrentLegendreOrder( legBasis, legTest );
  //
  //      for ( auto itOut = outerElems.begin( ); itOut != outerElems.end( ); itOut++ ) {
  //        for ( auto itIn = innerElems.begin( ); itIn != innerElems.end( ); itIn++ ) {
  //          integrator.getElemMatrix1Layer( *itOut, *itIn, elemMatrix );
  //          this->space->getOuterElemDOFs( *itOut, &rowIndices[0] );
  //          this->space->getInnerElemDOFs( *itIn, &colIndices[0] );
  //          for ( int k = 0; k < nLocalRows; k++ ) {
  //            rowIndices[k] += legTest * this->space->getOuterDOFs( );
  //          }
  //          for ( int k = 0; k < nLocalCols; k++ ) {
  //            colIndices[k] += legBasis * this->space->getInnerDOFs( );
  //          }
  //          block->addToPositions( rowIndices, colIndices, elemMatrix );
  //        }
  //      }
  //    }
  //  }
  //  matrix.setBlock( 1, 2, block, true );
  //  block->makeCompressed( );
  //
  //  // copy pointer to upper diagonal
  //  for ( int j = 2; j < nTimeSteps - 1; j++ ) {
  //    matrix.setBlock( j, j + 1, block, false );
  //  }
  //
  //  // the last column
  //  for ( int j = nTimeSteps - 2; j < nTimeSteps; j++ ) { // test
  //    SparseMatrix<LO, SC> *block = new SparseMatrix<LO, SC>( blockRows,
  //        blockCols );
  //
  //    integrator.setCurrentFunctions( nTimeSteps - 1, j );
  //
  //    // two loops over Legendre polynomials in basis
  //    for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
  //      for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {
  //        integrator.setCurrentLegendreOrder( legBasis, legTest );
  //
  //        // inner two loops over space
  //        for ( auto itOut = outerElems.begin( ); itOut != outerElems.end( ); itOut++ ) {
  //          for ( auto itIn = innerElems.begin( ); itIn != innerElems.end( ); itIn++ ) {
  //            integrator.getElemMatrix1Layer( *itOut, *itIn, elemMatrix );
  //            this->space->getOuterElemDOFs( *itOut, &rowIndices[0] );
  //            this->space->getInnerElemDOFs( *itIn, &colIndices[0] );
  //            for ( int k = 0; k < nLocalRows; k++ ) {
  //              rowIndices[k] += legTest * this->space->getOuterDOFs( );
  //            }
  //            for ( int k = 0; k < nLocalCols; k++ ) {
  //              colIndices[k] += legBasis * this->space->getInnerDOFs( );
  //            }
  //            block->addToPositions( rowIndices, colIndices, elemMatrix );
  //          }
  //        }
  //      }
  //    }
  //    block->makeCompressed( );
  //    matrix.setBlock( j, nTimeSteps - 1, block, true );
  //  }


}

template<class LO, class SC>
void BEBilinearFormWave1Layer<LO, SC>::assembleBlock(
    LO row,
    LO column,
    SparseMatrix<LO, SC> &localBlock
    ) const {

  // cast to bespace to bespacetime
  BESpaceTime<LO, SC> *spaceTime = static_cast<BESpaceTime<LO, SC>*> ( this->space );
  BEIntegratorWave<LO, SC> integrator( this->space, this->quadratureOrder, this->timeQuadOrder );
  const vector<LO> innerElems = this->space->getInnerElems( );
  const vector<LO> outerElems = this->space->getOuterElems( );

  LO iMax = outerElems.size( );
  LO jMax = innerElems.size( );

  integrator.setCurrentFunctions( column, row );

  // max order of Legendre polynomials in temporal basis
  int legOrder = spaceTime->getLegendreOrder( );

  // allocate local element matrices
  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  // preallocate the matrix
  localBlock.resize( ( legOrder + 1 ) * this->space->getOuterDOFs( ),
      ( legOrder + 1 ) * this->space->getInnerDOFs( ) );


  // two loops over Legendre polynomials in basis
  for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
    for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {
      integrator.setCurrentLegendreOrder( legBasis, legTest );

#pragma omp parallel shared(localBlock)
      {
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

#pragma omp for schedule(dynamic)
        for ( LO i = 0; i < iMax; i++ ) {
          for ( LO j = 0; j < jMax; j++ ) {
            integrator.getElemMatrix1Layer( outerElems[i], innerElems[j],
                *( matrixBuffer[counter] ) );
            this->space->getOuterElemDOFs( outerElems[i],
                &( *( rowIndicesBuffer[counter] ) )[0] );
            this->space->getInnerElemDOFs( innerElems[j],
                &( *( colIndicesBuffer[counter] ) )[0] );
            // correct for the current Legendre polynomial
            for ( int k = 0; k < nLocalRows; k++ ) {
              ( *( rowIndicesBuffer[counter] ) )[k] +=
                  legTest * this->space->getOuterDOFs( );
            }
            for ( int k = 0; k < nLocalCols; k++ ) {
              ( *( colIndicesBuffer[counter] ) )[k] +=
                  legBasis * this->space->getInnerDOFs( );
            }
            counter++;
            if ( counter == BUFFER_SIZE ) {
              counter = 0;
#pragma omp critical
              {
                for ( int m = 0; m < BUFFER_SIZE; m++ ) {
                  for ( int k = 0; k < nLocalRows; k++ ) {
                    for ( int l = 0; l < nLocalCols; l++ ) {
                      if ( fabs( ( *( matrixBuffer[m] ) ).get( k, l ) ) > EPS ) {
                        localBlock.getEigenSparseMatrix( )->coeffRef(
                            ( *( rowIndicesBuffer[m] ) )[k],
                            ( *( colIndicesBuffer[m] ) )[l] ) +=
                            ( *( matrixBuffer[m] ) ).get( k, l );
                        //localTripletList.push_back( Eigen::Triplet<SC, LO>( ( *( rowIndicesBuffer[m] ) )[k], ( *( colIndicesBuffer[m] ) )[l], ( *( matrixBuffer[m] ) ).get( k, l ) ) );
                      }
                    }
                  }
                }
              }
            }
          }
        }
#pragma omp critical
        {
          for ( int m = 0; m < counter; m++ ) {
            for ( int k = 0; k < nLocalRows; k++ ) {
              for ( int l = 0; l < nLocalCols; l++ ) {
                if ( fabs( ( *( matrixBuffer[m] ) ).get( k, l ) ) > EPS ) {
                  localBlock.getEigenSparseMatrix( )->coeffRef(
                      ( *( rowIndicesBuffer[m] ) )[k],
                      ( *( colIndicesBuffer[m] ) )[l] ) +=
                      ( *( matrixBuffer[m] ) ).get( k, l );
                  //localTripletList.push_back( Eigen::Triplet<SC, LO>( ( *( rowIndicesBuffer[m] ) )[k], ( *( colIndicesBuffer[m] ) )[l], ( *( matrixBuffer[m] ) ).get( k, l ) ) );
                }
              }
            }
          }
        }
        for ( int i = 0; i < BUFFER_SIZE; i++ ) {
          delete matrixBuffer[i];
          delete rowIndicesBuffer[i];
          delete colIndicesBuffer[i];
        }
      }
    }
  }



}

}
#endif
