/*!
 * @file    BEBilinearFormWaveHypersingular.cpp
 * @author  Michal Merta
 * @date    January 1, 2014
 * 
 */

#ifdef BEBILINEARFORMWAVEHYPERSINGULAR_H

// use ordinary basis functions described in Sauter, Veit: Retarded time-domain...
#ifndef EXPERIMENTAL_WAVE

namespace bem4i {

template<class LO, class SC>
BEBilinearFormWaveHypersingular<LO, SC>::BEBilinearFormWaveHypersingular( ) {
}

template<class LO, class SC>
BEBilinearFormWaveHypersingular<LO, SC>::BEBilinearFormWaveHypersingular(
    const BEBilinearFormWaveHypersingular& orig
    ) {
}

template<class LO, class SC>
BEBilinearFormWaveHypersingular<LO, SC>::~BEBilinearFormWaveHypersingular( ) {
}

template<class LO, class SC>
BEBilinearFormWaveHypersingular<LO, SC>::BEBilinearFormWaveHypersingular(
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

  this->quadratureOrderDisjointElems = quadratureOrderDisjointElems;

  this->timeQuadOrder = timeQuadratureOrder;

  this->nPos = 0;
  this->nPre = 0;
}

template<class LO, class SC>
void BEBilinearFormWaveHypersingular<LO, SC>::assemble(
    SparseMatrix<LO, SC>& matrix
    ) const {

  BESpaceTime<LO, SC> *spaceTime = static_cast<BESpaceTime<LO, SC>*> (
      this->space );
  // max order of Legendre polynomials in temporal basis
  int legOrder = spaceTime->getLegendreOrder( );

  int nTimeSteps = spaceTime->getNTimeSteps( );
  LO nnz;
  LO length = 0;

  matrix.resize( nTimeSteps * this->space->getOuterDOFs( ) * ( legOrder + 1 ),
      nTimeSteps * this->space->getInnerDOFs( ) * ( legOrder + 1 ) );

  // create Eigen triplets and insert values into matrix
  // (duplicated indices are summed!)
  std::vector<Eigen::Triplet<SC, LO> > tripletList;
  tripletList.reserve( nTimeSteps * this->space->getInnerDOFs( ) *
      this->space->getInnerDOFs( ) * ( legOrder + 1 ) );

  // at first, assemble the first column
  //assembleFirstColumn(tripletList);
  for ( int i = 0; i < nTimeSteps; i++ ) {
    assembleBlock( i, 0, tripletList, nnz );
  }


  // block on position (0,1) is single
  assembleBlock( 0, 1, tripletList, nnz );

  // blocks in the second column (except the first and the last) are copied to
  // the rest of the matrix
  for ( int i = 1; i < nTimeSteps - 1; i++ ) {
    assembleBlock( i, 1, tripletList, nnz );
    length = tripletList.size( );
    for ( int j = 1; j < nTimeSteps - i - 1; j++ ) {
      for ( int k = 0; k < nnz; k++ ) {
        Eigen::Triplet<SC, LO> newTriplet(
            tripletList[length - nnz + k].row( ) + j *
            this->space->getOuterDOFs( ) * ( legOrder + 1 ),
            tripletList[length - nnz + k].col( ) + j *
            this->space->getInnerDOFs( ) * ( legOrder + 1 ),
            tripletList[length - nnz + k].value( ) );
        tripletList.push_back( newTriplet );
      }
    }
  }

  // upper diagonal blocks
  assembleBlock( 1, 2, tripletList, nnz );
  length = tripletList.size( );
  for ( int j = 1; j < nTimeSteps - 3; j++ ) {
    for ( int k = 0; k < nnz; k++ ) {
      Eigen::Triplet<SC, LO> newTriplet(
          tripletList[length - nnz + k].row( ) + j *
          this->space->getOuterDOFs( ) * ( legOrder + 1 ),
          tripletList[length - nnz + k].col( ) + j *
          this->space->getInnerDOFs( ) * ( legOrder + 1 ),
          tripletList[length - nnz + k].value( ) );
      tripletList.push_back( newTriplet );
    }
  }

  // the last column
  assembleBlock( nTimeSteps - 2, nTimeSteps - 1, tripletList, nnz );
  assembleBlock( nTimeSteps - 1, nTimeSteps - 1, tripletList, nnz );

  // the last row
  for ( int i = 1; i < nTimeSteps - 1; i++ ) {
    assembleBlock( nTimeSteps - 1, i, tripletList, nnz );
  }
  std::cout << " assembled local triplets" << std::endl;
  matrix.getEigenSparseMatrix( )->setFromTriplets( tripletList.begin( ),
      tripletList.end( ) );
  matrix.makeCompressed( );
}

template<class LO, class SC>
void BEBilinearFormWaveHypersingular<LO, SC>::assemble(
    BlockMatrix<LO, SC>& matrix
    ) const {

  BESpaceTime<LO, SC> *spaceTime =
      static_cast<BESpaceTime<LO, SC>*> ( this->space );

  // max order of Legendre polynomials in temporal basis
  int legOrder = spaceTime->getLegendreOrder( );

  SparseMatrix<LO, SC> *block;

  int nTimeSteps = spaceTime->getNTimeSteps( );

  LO *numbersOfRows = new LO[nTimeSteps];
  LO *numbersOfCols = new LO[nTimeSteps];
  for ( LO i = 0; i < nTimeSteps; i++ ) {
    numbersOfRows[i] = ( legOrder + 1 ) * this->space->getOuterDOFs( );
    numbersOfCols[i] = ( legOrder + 1 ) * this->space->getInnerDOFs( );
  }

  matrix.resize( nTimeSteps, nTimeSteps, numbersOfRows, numbersOfCols );

  delete [] numbersOfCols;
  delete [] numbersOfRows;

  // at first, assemble the first column
  for ( int i = 0; i < nTimeSteps; i++ ) {

    block = new SparseMatrix<LO, SC>(
        ( legOrder + 1 ) * this->space->getOuterDOFs( ),
        ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
    assembleBlock( i, 0, *block );
    block->makeCompressed( );
    matrix.setBlock( i, 0, block, true );
  }

  // block on position (0,1) is single
  block = new SparseMatrix<LO, SC>(
      ( legOrder + 1 ) * this->space->getOuterDOFs( ),
      ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
  assembleBlock( 0, 1, *block );
  block->makeCompressed( );
  matrix.setBlock( 0, 1, block, true );

  // blocks in the second column 
  for ( int i = 1; i < nTimeSteps - 1; i++ ) {
    block = new SparseMatrix<LO, SC>(
        ( legOrder + 1 ) * this->space->getOuterDOFs( ),
        ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
    assembleBlock( i, 1, *block );
    block->makeCompressed( );
    matrix.setBlock( i, 1, block, true );
  }



  // the last row
  for ( int i = 1; i < nTimeSteps - 1; i++ ) {
    block = new SparseMatrix<LO, SC>(
        ( legOrder + 1 ) * this->space->getOuterDOFs( ),
        ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
    assembleBlock( nTimeSteps - 1, i, *block );
    block->makeCompressed( );
    matrix.setBlock( nTimeSteps - 1, i, block, true );
  }

  // the upper diagonal
  block = new SparseMatrix<LO, SC>(
      ( legOrder + 1 ) * this->space->getOuterDOFs( ),
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
  block = new SparseMatrix<LO, SC>(
      ( legOrder + 1 ) * this->space->getOuterDOFs( ),
      ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
  assembleBlock( nTimeSteps - 2, nTimeSteps - 1, *block );
  block->makeCompressed( );
  matrix.setBlock( nTimeSteps - 2, nTimeSteps - 1, block, true );

  block = new SparseMatrix<LO, SC>(
      ( legOrder + 1 ) * this->space->getOuterDOFs( ),
      ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
  assembleBlock( nTimeSteps - 1, nTimeSteps - 1, *block );
  block->makeCompressed( );
  matrix.setBlock( nTimeSteps - 1, nTimeSteps - 1, block, true );


}

template<class LO, class SC>
void BEBilinearFormWaveHypersingular<LO, SC>::assemble(
    MPIBlockMatrix<LO, SC>& matrix
    ) const {
  // timeval start, stop;
  BESpaceTime<LO, SC> *spaceTime = static_cast<BESpaceTime<LO, SC>*>
      ( this->space );

  SparseMatrix<LO, SC> *block = nullptr;

  // max order of Legendre polynomials in temporal basis
  int legOrder = spaceTime->getLegendreOrder( );

  int nTimeSteps = spaceTime->getNTimeSteps( );

  LO *numbersOfRows = new LO[nTimeSteps];
  LO *numbersOfCols = new LO[nTimeSteps];
  for ( LO i = 0; i < nTimeSteps; i++ ) {
    numbersOfRows[i] = ( legOrder + 1 ) * this->space->getOuterDOFs( );
    numbersOfCols[i] = ( legOrder + 1 ) * this->space->getInnerDOFs( );
  }

  // prepare the map of blocks per rank
  int *ranks = new int[nTimeSteps * nTimeSteps];
  //memset( ranks, 0, nTimeSteps * nTimeSteps * sizeof (int ) );
  for ( int i = 0; i < nTimeSteps * nTimeSteps; i++ ) {
    ranks[i] = -1;
  }
  int m, n, nProc;
  MPI_Comm_size( MPI_COMM_WORLD, &nProc );
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // compute number of nonzero rows in the first column
  SC diam = 2.0 * this->space->getMesh( )->getRadius( );
  int nZBlocks = ceil( diam / spaceTime->getDt( ) + 1.0 ) + 1;

  if ( nZBlocks > nTimeSteps ) {
    nZBlocks = nTimeSteps;
  }
  //compute the total number of non-zero blocks (excpt first col. and last row)
  int totNZBlocks = ( (double) ( ( nTimeSteps - 1 ) * ( nTimeSteps - 1 ) ) ) -
      ( (double) ( nTimeSteps - 1 )*( nTimeSteps - 2 ) ) / 2.0 -
      ( (double) std::max( 0, ( nTimeSteps - nZBlocks - 1 ) ) *
      ( nTimeSteps - nZBlocks - 2 ) ) / 2.0;

  int blocksPerProc = ceil( (double) totNZBlocks / (double) nProc );

  int owner = 0, counter = 0, totCounter = 0;
  bool breakLoop = false;
  int *nBlocksPerRank = new int[nProc];
  memset( nBlocksPerRank, 0, nProc * sizeof (int ) );

  for ( int i = 0; i < nTimeSteps - 1; i++ ) {
    for ( int j = 0; j < nTimeSteps - 1 - i; j++ ) {
      m = i + j;
      n = 1 + j;

      ranks[n * nTimeSteps + m] = owner;

      if ( ( n == 1 ) || ( m == 1 && n == 2 ) || ( ( m == nTimeSteps - 2 ) &&
          ( n == nTimeSteps - 1 ) ) || ( counter == 0 ) ) {
        // count the number of assembled blocks on each process
        nBlocksPerRank[owner]++;
      }

      counter++;
      totCounter++;

      if ( counter == blocksPerProc ) {
        counter = 0;
        owner++;
      }
      if ( totCounter == totNZBlocks ) {
        breakLoop = true;
        break;
      }
    }
    if ( breakLoop ) {
      break;
    }
  }
  // the remaining blocks are put on the processes with the lowest amount of blocks
  for ( int i = 0; i < nZBlocks; i++ ) {
    owner = std::distance( nBlocksPerRank, std::min_element(
        nBlocksPerRank, nBlocksPerRank + nProc ) );
    ranks[i] = owner;
    nBlocksPerRank[owner]++;
  }

  for ( int i = nTimeSteps - nZBlocks; i < nTimeSteps; i++ ) {
    owner = std::distance( nBlocksPerRank, std::min_element(
        nBlocksPerRank, nBlocksPerRank + nProc ) );
    ranks[i * nTimeSteps + nTimeSteps - 1] = owner;
    nBlocksPerRank[owner]++;
  }


  if ( rank == 0 )
    for ( int i = 0; i < nTimeSteps; i++ ) {
      for ( int j = 0; j < nTimeSteps; j++ ) {
        // blocks distribution
        std::cout << ranks[j * nTimeSteps + i] << " ";
      }
      std::cout << std::endl;
    }


  for ( int i = 0; i < nProc; i++ ) {
    // numbers of blocks
    std::cout << nBlocksPerRank[i] << " ";
  }

  matrix.resize( nTimeSteps, nTimeSteps, numbersOfRows, numbersOfCols, ranks );
  //
  delete [] numbersOfCols;
  delete [] numbersOfRows;

  // if the matrix belongs to the current process, assemble it, otherwise set blocks to null
  double t1, t2;
  t1 = MPI_Wtime( );
  // at first, assemble the first column
  for ( int i = 0; i < nTimeSteps; i++ ) {
    if ( matrix.getOwner( i, 0 ) != -1 ) {
      block = new SparseMatrix<LO, SC>(
          ( legOrder + 1 ) * this->space->getOuterDOFs( ),
          ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
      std::cout << " before " << std::endl;
      assembleMPIBlock( i, 0, *block, matrix.getOwner( i, 0 ) );
      std::cout << i << std::endl;
      // throw away blocks on nonowner processors
      if ( rank != matrix.getOwner( i, 0 ) ) {
        delete block;
      }
    }
    if ( rank == matrix.getOwner( i, 0 ) ) {
      block->makeCompressed( );
      matrix.setBlock( i, 0, block, true );
    } else {
      matrix.setBlock( i, 0, nullptr, false );
    }
  }
  t2 = MPI_Wtime( );
  printf( "First column assembly time: %1.2f\n", t2 - t1 );
  fflush( stdout );

  // block on position (0,1) is single
  t1 = MPI_Wtime( );
  block = new SparseMatrix<LO, SC>(
      ( legOrder + 1 ) * this->space->getOuterDOFs( ),
      ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
  assembleMPIBlock( 0, 1, *block, matrix.getOwner( 0, 1 ) );
  if ( rank == matrix.getOwner( 0, 1 ) ) {
    block->makeCompressed( );
    matrix.setBlock( 0, 1, block, true );
  } else {
    matrix.setBlock( 0, 1, nullptr, false );
    delete block;
  }

  // blocks in the second column 
  for ( int i = 1; i < nTimeSteps - 1; i++ ) {
    if ( matrix.getOwner( i, 1 ) != -1 ) {
      block = new SparseMatrix<LO, SC>(
          ( legOrder + 1 ) * this->space->getOuterDOFs( ),
          ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
      assembleMPIBlock( i, 1, *block, matrix.getOwner( i, 1 ) );
      // throw away blocks on nonowner processors
      if ( rank != matrix.getOwner( i, 1 ) ) {
        delete block;
      }
    }
    if ( rank == matrix.getOwner( i, 1 ) ) {
      block->makeCompressed( );
      matrix.setBlock( i, 1, block, true );
    } else {
      matrix.setBlock( i, 1, nullptr, false );
    }
  }
  //
  t2 = MPI_Wtime( );
  printf( "Second column assembly time: %1.2f\n", t2 - t1 );
  fflush( stdout );
  t1 = MPI_Wtime( );
  // the last row
  for ( int i = 1; i < nTimeSteps - 1; i++ ) {
    if ( matrix.getOwner( nTimeSteps - 1, i ) != -1 ) {
      block = new SparseMatrix<LO, SC>(
          ( legOrder + 1 ) * this->space->getOuterDOFs( ),
          ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
      assembleMPIBlock( nTimeSteps - 1, i, *block,
          matrix.getOwner( nTimeSteps - 1, i ) );
      // throw away blocks on nonowner processors
      if ( rank != matrix.getOwner( nTimeSteps - 1, i ) ) {
        delete block;
      }
    }
    if ( rank == matrix.getOwner( nTimeSteps - 1, i ) ) {
      block->makeCompressed( );
      matrix.setBlock( nTimeSteps - 1, i, block, true );
    } else {
      matrix.setBlock( nTimeSteps - 1, i, nullptr, false );
    }
  }
  t2 = MPI_Wtime( );
  printf( "Last row assembly time: %1.2f\n", t2 - t1 );
  fflush( stdout );

  t1 = MPI_Wtime( );
  // the upper diagonal
  block = new SparseMatrix<LO, SC>(
      ( legOrder + 1 ) * this->space->getOuterDOFs( ),
      ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
  assembleMPIBlock( 1, 2, *block, matrix.getOwner( 1, 2 ) );
  if ( matrix.getOwner( 1, 2 ) == rank ) {
    block->makeCompressed( );
    matrix.setBlock( 1, 2, block, true );
  } else {
    delete block;
    matrix.setBlock( 1, 2, nullptr, false );
  }
  block = ( SparseMatrix<LO, SC>* ) matrix.getBlock( 1, 2 );
  int sentFrom = -1;
  int sentTo = -1;
  for ( int j = 2; j < nTimeSteps - 1; j++ ) {
    if ( ( rank == matrix.getOwner( 1, 2 )
        || rank == matrix.getOwner( j, j + 1 ) )
        && ( matrix.getOwner( 1, 2 ) != matrix.getOwner( j, j + 1 ) )
        && ( sentTo != matrix.getOwner( j, j + 1 )
        || sentFrom != matrix.getOwner( 1, 2 ) ) ) {
      sendBlock( matrix, 1, 2, j, j + 1, matrix.getOwner( 1, 2 ),
          matrix.getOwner( j, j + 1 ) );

      sentTo = matrix.getOwner( j, j + 1 );
      sentFrom = matrix.getOwner( 1, 2 );
      if ( rank == matrix.getOwner( j, j + 1 ) ) {
        block = ( SparseMatrix<LO, SC>* ) matrix.getBlock( j, j + 1 );
      }
    } else {
      if ( rank == matrix.getOwner( j, j + 1 ) ) {
        matrix.setBlock( j, j + 1, block, false );
        //std::cout << "setting block " << j << " " << 1+j << " to " << "1 2 " << " on rank " << rank <<  std::endl;
      }
    }
  }
  t2 = MPI_Wtime( );
  printf( "Upper diagonal assembly time: %1.2f\n", t2 - t1 );
  fflush( stdout );
  MPI_Barrier( MPI_COMM_WORLD );

  sentFrom = -1;
  sentTo = -1;
  //
  //
  // copy pointers to the second column blocks on the appropriate positions
  t1 = MPI_Wtime( );
  for ( int i = 1; i < nTimeSteps - 2; i++ ) {
    block = ( SparseMatrix<LO, SC>* ) matrix.getBlock( i, 1 );
    for ( int j = 1; j < nTimeSteps - i - 1; j++ ) {

      if ( ( rank == matrix.getOwner( i, 1 )
          || rank == matrix.getOwner( i + j, 1 + j ) )
          && matrix.getOwner( i, 1 ) != matrix.getOwner( i + j, 1 + j )
          && ( sentTo != matrix.getOwner( i + j, j + 1 )
          || sentFrom != matrix.getOwner( i, 1 ) ) ) {
        //std::cout << "source: "<<matrix.getOwner( i, 1 ) << std::endl;
        //std::cout << "desst: "<<matrix.getOwner(i+j, 1+j ) << std::endl;
        sendBlock( matrix, i, 1, i + j, 1 + j, matrix.getOwner( i, 1 ),
            matrix.getOwner( i + j, 1 + j ) );
        sentTo = matrix.getOwner( i + j, 1 + j );
        sentFrom = matrix.getOwner( i, 1 );
        if ( rank == matrix.getOwner( i + j, 1 + j ) ) {
          block = ( SparseMatrix<LO, SC>* ) matrix.getBlock( i + j, j + 1 );
        }
      } else {
        if ( rank == matrix.getOwner( i + j, 1 + j ) ) {
          matrix.setBlock( i + j, 1 + j, block, false );
          //std::cout << "setting block " << i+j << " " << 1+j << " to " << i << " " <<"1" <<" on rank " << rank <<   std::endl;
        }
      }
    }
  }
  t2 = MPI_Wtime( );
  printf( "Second column copy time: %1.2f\n", t2 - t1 );
  fflush( stdout );
  MPI_Barrier( MPI_COMM_WORLD );


  // the last column
  t1 = MPI_Wtime( );
  block = new SparseMatrix<LO, SC>(
      ( legOrder + 1 ) * this->space->getOuterDOFs( ),
      ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
  assembleMPIBlock( nTimeSteps - 2, nTimeSteps - 1, *block,
      matrix.getOwner( nTimeSteps - 2, nTimeSteps - 1 ) );
  if ( matrix.getOwner( nTimeSteps - 2, nTimeSteps - 1 ) == rank ) {
    block->makeCompressed( );
    matrix.setBlock( nTimeSteps - 2, nTimeSteps - 1, block, true );
  } else {
    matrix.setBlock( nTimeSteps - 2, nTimeSteps - 1, nullptr, false );
    delete block;
  }

  block = new SparseMatrix<LO, SC>(
      ( legOrder + 1 ) * this->space->getOuterDOFs( ),
      ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
  assembleMPIBlock( nTimeSteps - 1, nTimeSteps - 1, *block,
      matrix.getOwner( nTimeSteps - 1, nTimeSteps - 1 ) );
  if ( matrix.getOwner( nTimeSteps - 1, nTimeSteps - 1 ) == rank ) {
    block->makeCompressed( );
    matrix.setBlock( nTimeSteps - 1, nTimeSteps - 1, block, true );
  } else {
    matrix.setBlock( nTimeSteps - 1, nTimeSteps - 1, nullptr, false );
    delete block;
  }
  t2 = MPI_Wtime( );
  printf( "Last column assembly time: %1.2f\n", t2 - t1 );
  fflush( stdout );

  MPI_Barrier( MPI_COMM_WORLD );
  /*  for (int r = 0; r < nProc; r++) {
    if (rank == r) {
    for (int i = 0; i < 20; i++) {
      for (int j = 0; j < 20; j++) {
        if (rank == matrix.getOwner(i, j)) {
        if (matrix.getBlock(i, j ) != nullptr) {
        std::cout << i << " " <<j <<std::endl;
          Vector<LO, SC> tst(matrix.getBlock(i,j)->getNRows());
          Vector<LO, SC> res(matrix.getBlock(i,j)->getNRows());
          tst.setAll(1.0);
            matrix.getBlock(i,j)->apply(tst, res);
            res.print();
            }
      }
    }
    }
  }  
  MPI_Barrier(MPI_COMM_WORLD);
    }
   */

  delete [] ranks;
  delete [] nBlocksPerRank;
}

template<class LO, class SC>
void BEBilinearFormWaveHypersingular<LO, SC>::sendBlock(
    MPIBlockMatrix<LO, SC>& matrix,
    int sourceRow,
    int sourceColumn,
    int destRow,
    int destColumn,
    int source,
    int destination ) const {

  int rank, nProc;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &nProc );

  MPI_Datatype MPI_LO = GetType<LO, SC>::MPI_LO( );
  MPI_Datatype MPI_SC = GetType<LO, SC>::MPI_SC( );

  double t1, t2;

  if ( rank == source ) {

    SparseMatrix<LO, SC> *localBlock =
        ( SparseMatrix<LO, SC>* ) matrix.getBlock( sourceRow, sourceColumn );
    // serialize my block and prepare it for sending
    if ( localBlock->getEigenSparseMatrix( ) != nullptr
        && localBlock->getEigenSparseMatrix( )->nonZeros( ) > 0 ) {

      int nnz = localBlock->getEigenSparseMatrix( )->nonZeros( );
      LO nRows = localBlock->getNRows( );
      LO nCols = localBlock->getNCols( );
      int outSize = localBlock->getEigenSparseMatrix( )->outerSize( );
      int inSize = localBlock->getEigenSparseMatrix( )->innerSize( );

      //      for ( int k = 0; k < localBlock->getEigenSparseMatrix( )->outerSize( );
      //          ++k ) {
      //        for ( typename Eigen::SparseMatrix<SC>::InnerIterator it(
      //            *localBlock->getEigenSparseMatrix( ), k ); it; ++it ) {
      //          //std::cout << counter << std::endl;
      //          vals[counter] = it.value( );
      //          rows[counter] = it.row( );
      //          cols[counter] = it.col( );
      //          counter++;
      //        }
      //      }

      // send block to destination process
      MPI_Send( &nnz, 1, MPI_INT, destination, 0, MPI_COMM_WORLD );

      t1 = MPI_Wtime( );
      MPI_Send( &nRows, 1, MPI_LO, destination, 0, MPI_COMM_WORLD );
      MPI_Send( &nCols, 1, MPI_LO, destination, 0, MPI_COMM_WORLD );
      MPI_Send( &outSize, 1, MPI_INT, destination, 0, MPI_COMM_WORLD );
      MPI_Send( &inSize, 1, MPI_INT, destination, 0, MPI_COMM_WORLD );
      MPI_Send( localBlock->getEigenSparseMatrix( )->valuePtr( ),
          nnz, MPI_SC, destination, 0, MPI_COMM_WORLD );
      MPI_Send( localBlock->getEigenSparseMatrix( )->outerIndexPtr( ),
          outSize, MPI_LO, destination, 0, MPI_COMM_WORLD );
      MPI_Send( localBlock->getEigenSparseMatrix( )->innerIndexPtr( ),
          nnz, MPI_LO, destination, 0, MPI_COMM_WORLD );
      t2 = MPI_Wtime( );
      printf( "Send time: %1.2f\n", t2 - t1 );
      fflush( stdout );
      std::cout << "Sending block " << sourceRow << ", " << sourceColumn << " from rank " << rank << " to rank " << destination << " on position " << destRow << ", " << destColumn << std::endl;

    } else {
      int zero = 0;
      MPI_Send( &zero, 1, MPI_INT, destination, 0, MPI_COMM_WORLD );
    }
  } else if ( rank == destination ) {

    SparseMatrix<LO, SC> *localBlock =
        ( SparseMatrix<LO, SC>* ) matrix.getBlock( destRow, destColumn );

    LO nnz, nRows, nCols, outSize, inSize;
    MPI_Status status;

    MPI_Recv( &nnz, 1, MPI_INT, source, 0, MPI_COMM_WORLD, &status );

    if ( nnz > 0 ) {
      if ( localBlock != nullptr ) {
        delete localBlock;
      }

      MPI_Recv( &nRows, 1, MPI_LO, source, 0, MPI_COMM_WORLD, &status );
      MPI_Recv( &nCols, 1, MPI_LO, source, 0, MPI_COMM_WORLD, &status );

      MPI_Recv( &outSize, 1, MPI_INT, source, 0, MPI_COMM_WORLD, &status );
      MPI_Recv( &inSize, 1, MPI_INT, source, 0, MPI_COMM_WORLD, &status );

      //      LO *rows = new LO[nnz];
      //      LO *cols = new LO[nnz];
      //      SCVT *vals = new SCVT[nnz];

      localBlock = new SparseMatrix<LO, SC>( nRows, nCols );
      matrix.setBlock( destRow, destColumn, localBlock, true );

      localBlock->getEigenSparseMatrix( )->resize( nRows, nCols );
      localBlock->getEigenSparseMatrix( )->makeCompressed( );
      localBlock->getEigenSparseMatrix( )->resizeNonZeros( nnz );

      MPI_Recv( localBlock->getEigenSparseMatrix( )->valuePtr( ),
          nnz, MPI_SC, source, 0, MPI_COMM_WORLD, &status );
      MPI_Recv( localBlock->getEigenSparseMatrix( )->outerIndexPtr( ),
          outSize, MPI_LO, source, 0, MPI_COMM_WORLD, &status );
      MPI_Recv( localBlock->getEigenSparseMatrix( )->innerIndexPtr( ),
          nnz, MPI_LO, source, 0, MPI_COMM_WORLD, &status );
      localBlock->getEigenSparseMatrix( )->finalize( );
    } else {
      matrix.setBlock( destRow, destColumn, nullptr, false );
    }
  }
}


//template<class LO, class SC>
//void BEBilinearFormWaveHypersingular<LO, SC>::assemble( MPIBlockMatrix<LO, SC>& matrix ) const {
//
//  BESpaceTime<LO, SC> *spaceTime = static_cast<BESpaceTime<LO, SC>*> ( this->space );
//  BEIntegratorWave<LO, SC> integrator( this->space, this->quadratureOrder, this->timeQuadOrder );
//  const vector<LO> innerElems = this->space->getInnerElems( );
//  const vector<LO> outerElems = this->space->getOuterElems( );
//  SparseMatrix<LO, SC> *block;
//
//  int nTimeSteps = spaceTime->getNTimeSteps( );
//  
//  LO *nnz = new LO[nTimeSteps];
//  estimateNnz(nnz);
//  
// 
//  LO *numbersOfRows = new LO[nTimeSteps];
//  LO *numbersOfCols = new LO[nTimeSteps];
//  for ( LO i = 0; i < nTimeSteps; i++ ) {
//    numbersOfRows[i] = this->space->getOuterDOFs( );
//    numbersOfCols[i] = this->space->getInnerDOFs( );
//  }
//
//  // prepare the map of blocks per rank
//  int *ranks = new int[nTimeSteps * nTimeSteps];
//  memset( ranks, 0, nTimeSteps * nTimeSteps * sizeof (int ) );
//  int m, n, nProc;
//  MPI_Comm_size( MPI_COMM_WORLD, &nProc );
//  int rank;
//  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
//
//  // compute number of nonzero rows in the first column
//  SC diam = 2.0 * this->space->getMesh( )->getRadius( );
//  int nZBlocks = ceil( diam / spaceTime->getDt( ) + 1.0 ) + 1;
//
//  if ( nZBlocks > nTimeSteps ) {
//    nZBlocks = nTimeSteps;
//  }
//  // compute the total number of non-zero blocks
//  int totNZBlocks = ( ( (double) ( nTimeSteps * nTimeSteps ) ) / 2.0 - ( (double) nTimeSteps ) / 2.0 )
//      - ( (double) ( nTimeSteps - nZBlocks - 1 ) * ( nTimeSteps - nZBlocks - 2 ) ) / 2.0;
//
//  int blocksPerProc = ceil( (double) totNZBlocks / (double) nProc );
//
//  int owner = 0, counter = 0, totCounter = 0;
//  bool breakLoop = false;
//  int *nBlocksPerRank = new int[nProc];
//  memset( nBlocksPerRank, 0, nProc * sizeof (int ) );
//  int *nnzPerRank = new int[nProc];
//  memset( nnzPerRank, 0, nProc * sizeof (int ) );  
//  
//  LO estTotNnz = 0; 
//  for (int i = 0 ; i< nTimeSteps - 1; i++) {
//    estTotNnz += nnz[i]*(nTimeSteps-i - 1);
//  }
//  int nnzPerCore = estTotNnz / nProc;
//  std::cout << nnzPerCore << std::endl;
//  
//  
//  for ( int i = 0; i < nTimeSteps - 1; i++ ) {
//    for ( int j = 0; j < nTimeSteps - 1 - i; j++ ) {
//      m = i + j;
//      n = 1 + j;
//
//      ranks[n * nTimeSteps + m] = owner;
//      
//      nnzPerRank[owner]+=nnz[i];
//      std::cout << nnzPerRank[owner ] << std::endl;
//
//      if ( ( n == 1 ) || ( m == 1 && n == 2 ) || ( ( m == nTimeSteps - 2 ) && ( n == nTimeSteps - 1 ) ) || ( counter == 0 ) ) {
//        // count the number of assembled blocks on each process
//        nBlocksPerRank[owner]++;
//      }
//
//      counter++;
//      totCounter++;
//
//      if ( nnzPerRank[owner] > nnzPerCore) {
//        owner++;
//      }
//      if ( totCounter == totNZBlocks ) {
//        breakLoop = true;
//        break;
//      }
//    }
//    if ( breakLoop ) {
//      break;
//    }
//  }
//  // the remaining blocks are put on the processes with the lowest amount of blocks
//  for ( int i = 0; i < nZBlocks; i++ ) {
//    owner = std::distance( nnzPerRank, std::min_element( nnzPerRank, nnzPerRank + nProc ) );
//    ranks[i] = owner;
//    nnzPerRank[owner]+=nnz[i];
//  }
//
//  for ( int i = nTimeSteps - nZBlocks; i < nTimeSteps; i++ ) {
//    owner = std::distance( nnzPerRank, std::min_element( nnzPerRank, nnzPerRank + nProc ) );
//    ranks[i * nTimeSteps + nTimeSteps - 1] = owner;
//    nnzPerRank[owner]+=nnz[i];
//  }
//
//
//  if ( rank == 0 )
//    for ( int i = 0; i < nTimeSteps; i++ ) {
//      for ( int j = 0; j < nTimeSteps; j++ ) {
//        // blocks distribution
//        std::cout << ranks[j * nTimeSteps + i] << " ";
//      }
//      std::cout << std::endl;
//    }
//
//
//  for ( int i = 0; i < nProc; i++ ) {
//    // numbers of blocks
//    std::cout << nBlocksPerRank[i] << " ";
//  }
//
//
//
//  matrix.resize( nTimeSteps, nTimeSteps, numbersOfRows, numbersOfCols, ranks );
//  //
//  delete [] numbersOfCols;
//  delete [] numbersOfRows;
//
//  // if the matrix belongs to the current process, assemble it, otherwise set blocks to null
//
//  // at first, assemble the first column
//  for ( int i = 0; i < nTimeSteps; i++ ) {
//    if ( matrix.getOwner( i, 0 ) == rank ) {
//      block = new SparseMatrix<LO, SC>( this->space->getOuterDOFs( ),
//          this->space->getInnerDOFs( ) );
//      assembleBlock( i, 0, *block );
//      block->makeCompressed( );
//      matrix.setBlock( i, 0, block, true );
//    } else {
//      matrix.setBlock( i, 0, nullptr, false );
//    }
//  }
//
//  // block on position (0,1) is single
//  if ( matrix.getOwner( 0, 1 ) == rank ) {
//    block = new SparseMatrix<LO, SC>( this->space->getOuterDOFs( ),
//        this->space->getInnerDOFs( ) );
//    assembleBlock( 0, 1, *block );
//    block->makeCompressed( );
//    matrix.setBlock( 0, 1, block, true );
//  } else {
//    matrix.setBlock( 0, 1, nullptr, false );
//  }
//
//  // blocks in the second column 
//  for ( int i = 1; i < nTimeSteps - 1; i++ ) {
//    if ( matrix.getOwner( i, 1 ) == rank ) {
//      block = new SparseMatrix<LO, SC>( this->space->getOuterDOFs( ),
//          this->space->getInnerDOFs( ) );
//      assembleBlock( i, 1, *block );
//      block->makeCompressed( );
//      matrix.setBlock( i, 1, block, true );
//    } else {
//      matrix.setBlock( i, 1, nullptr, false );
//    }
//  }
//  //
//  // the last column
//  if ( matrix.getOwner( nTimeSteps - 2, nTimeSteps - 1 ) == rank ) {
//    block = new SparseMatrix<LO, SC>( this->space->getOuterDOFs( ),
//        this->space->getInnerDOFs( ) );
//    assembleBlock( nTimeSteps - 2, nTimeSteps - 1, *block );
//    block->makeCompressed( );
//    matrix.setBlock( nTimeSteps - 2, nTimeSteps - 1, block, true );
//  } else {
//    matrix.setBlock( nTimeSteps - 2, nTimeSteps - 1, nullptr, false );
//  }
//
//  if ( matrix.getOwner( nTimeSteps - 1, nTimeSteps - 1 ) == rank ) {
//    block = new SparseMatrix<LO, SC>( this->space->getOuterDOFs( ),
//        this->space->getInnerDOFs( ) );
//    assembleBlock( nTimeSteps - 1, nTimeSteps - 1, *block );
//    block->makeCompressed( );
//    matrix.setBlock( nTimeSteps - 1, nTimeSteps - 1, block, true );
//  } else {
//    matrix.setBlock( nTimeSteps - 1, nTimeSteps - 1, nullptr, false );
//  }
//
//
//  // the last row
//  for ( int i = 1; i < nTimeSteps - 1; i++ ) {
//    if ( matrix.getOwner( nTimeSteps - 1, i ) == rank ) {
//      block = new SparseMatrix<LO, SC>( this->space->getOuterDOFs( ),
//          this->space->getInnerDOFs( ) );
//      assembleBlock( nTimeSteps - 1, i, *block );
//      block->makeCompressed( );
//      matrix.setBlock( nTimeSteps - 1, i, block, true );
//    } else {
//      matrix.setBlock( nTimeSteps - 1, i, nullptr, false );
//    }
//  }
//
//  // the upper diagonal
//  if ( matrix.getOwner( 1, 2 ) == rank ) {
//    block = new SparseMatrix<LO, SC>( this->space->getOuterDOFs( ),
//        this->space->getInnerDOFs( ) );
//    assembleBlock( 1, 2, *block );
//    block->makeCompressed( );
//    matrix.setBlock( 1, 2, block, true );
//  } else {
//    matrix.setBlock( 1, 2, nullptr, false );
//  }
//  block = ( SparseMatrix<LO, SC>* ) matrix.getBlock( 1, 2 );
//  for ( int j = 2; j < nTimeSteps - 1; j++ ) {
//    if ( matrix.getOwner( j, j + 1 ) == rank ) {
//      if ( block == nullptr ) {
//        block = new SparseMatrix<LO, SC>( this->space->getOuterDOFs( ), this->space->getInnerDOFs( ) );
//        assembleBlock( j, j + 1, *block );
//        block->makeCompressed( );
//        matrix.setBlock( j, j + 1, block, true );
//      } else {
//        matrix.setBlock( j, j + 1, block, false );
//      }
//    }
//  }
//
//
//  // copy pointers to the second column blocks on the appropriate positions
//  for ( int i = 1; i < nTimeSteps - 2; i++ ) {
//    block = ( SparseMatrix<LO, SC>* ) matrix.getBlock( i, 1 );
//    for ( int j = 1; j < nTimeSteps - i - 1; j++ ) {
//      if ( matrix.getOwner( i + j, 1 + j ) == rank ) {
//        if ( block == nullptr ) {
//          block = new SparseMatrix<LO, SC>( this->space->getOuterDOFs( ), this->space->getInnerDOFs( ) );
//          assembleBlock( i + j, 1 + j, *block );
//          block->makeCompressed( );
//          matrix.setBlock( i + j, 1 + j, block, true );
//        } else {
//          matrix.setBlock( i + j, 1 + j, block, false );
//        }
//      }
//    }
//  }
//  delete [] ranks;
//  delete [] nBlocksPerRank;
//}

template<class LO, class SC>
void BEBilinearFormWaveHypersingular<LO, SC>::estimateNnz(
    LO * nnz
    ) const {

  BESpaceTime<LO, SC> *spaceTime =
      static_cast<BESpaceTime<LO, SC>*> ( this->space );
  SC node1[3];
  SC node2[3];
  SC dist;
  for ( int i = 0; i < spaceTime->getNTimeSteps( ); i++ ) {
    nnz[i] = 0;
  }
  // watch out! this is O(n^2)!
  for ( int i = 0; i < this->space->getMesh( )->getNNodes( ); i++ ) {
    for ( int j = 0; j < this->space->getMesh( )->getNNodes( ); j++ ) {
      this->space->getMesh( )->getNode( i, node1 );
      this->space->getMesh( )->getNode( j, node2 );
      dist = std::sqrt( ( node1[0] - node2[0] )*( node1[0] - node2[0] )
          +( node1[1] - node2[1] )*( node1[1] - node2[1] ) +
          ( node1[2] - node2[2] )*( node1[2] - node2[2] ) );
      int idx = floor( dist / spaceTime->getDt( ) );
      nnz[idx]++;
      nnz[idx + 1]++;
      if ( idx != 0 ) {
        nnz[idx - 1]++;
      }
    }
  }
}

template<class LO, class SC>
void BEBilinearFormWaveHypersingular<LO, SC>::assemble(
    FullMatrix<LO, SC>& matrix
    ) const {

  BESpaceTime<LO, SC> *spaceTime =
      static_cast<BESpaceTime<LO, SC>*> ( this->space );
  BEIntegratorWave<LO, SC> integrator( this->space, this->quadratureOrder,
      this->timeQuadOrder, this->quadratureOrderDisjointElems );
  const vector<LO> innerElems = this->space->getInnerElems( );
  const vector<LO> outerElems = this->space->getOuterElems( );

  matrix.setAll( 0.0 );

  // allocate local element matrices
  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  int nTimeSteps = spaceTime->getNTimeSteps( );

  FullMatrix<LO, SC> elemMatrix( nLocalRows, nLocalCols );

  matrix.resize( nTimeSteps * this->space->getOuterDOFs( ),
      nTimeSteps * this->space->getInnerDOFs( ) );

  // vector of indices where to put values in the global matrix
  vector<LO> rowIndices( nLocalRows );
  vector<LO> colIndices( nLocalCols );

  // outer two loops over time blocks
  for ( int i = 0; i < nTimeSteps; i++ ) {
    for ( int j = 0; j < nTimeSteps; j++ ) {
      if ( i <= j + 1 ) {
        integrator.setCurrentFunctions( i, j );

        // inner two loops over space
        for ( auto itOut = outerElems.begin( ); itOut != outerElems.end( );
            itOut++ ) {
          for ( auto itIn = innerElems.begin( ); itIn != innerElems.end( );
              itIn++ ) {
            integrator.getElemMatrixHypersingular( *itOut, *itIn, elemMatrix );
            this->space->getOuterElemDOFs( *itOut, &rowIndices[0] );
            this->space->getInnerElemDOFs( *itIn, &colIndices[0] );
            for ( int k = 0; k < nLocalRows; k++ ) {
              rowIndices[k] += j * this->space->getOuterDOFs( );
            }
            for ( int k = 0; k < nLocalCols; k++ ) {
              colIndices[k] += i * this->space->getInnerDOFs( );
            }
            matrix.addToPositions( rowIndices, colIndices, elemMatrix );
          }
        }
      }
    }
  }
}

template<class LO, class SC>
void BEBilinearFormWaveHypersingular<LO, SC>::assembleBlock(
    LO row,
    LO column,
    std::vector<Eigen::Triplet<SC, LO> > &tripletList,
    LO &nnz
    ) const {

  // set number of nonzeros to zero
  nnz = 0;

  // set local triplet list
  std::vector<Eigen::Triplet<SC, LO> > localTripletList;

  BESpaceTime<LO, SC> *spaceTime = static_cast<BESpaceTime<LO, SC>*> (
      this->space );
  int nTimeSteps = spaceTime->getNTimeSteps( );
  int legOrder = spaceTime->getLegendreOrder( );

#ifdef EXPERIMENTAL_WAVE
  BEIntegratorWave<LO, SC> integrator( this->space, this->quadratureOrder, this->timeQuadOrder, nPre, nPos );
#else
  BEIntegratorWave<LO, SC> integrator( this->space, this->quadratureOrder,
      this->timeQuadOrder, this->quadratureOrderDisjointElems );
#endif
  const vector<LO> innerElems = this->space->getInnerElems( );
  const vector<LO> outerElems = this->space->getOuterElems( );

  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  LO iMax = outerElems.size( );
  LO jMax = innerElems.size( );

  integrator.setCurrentFunctions( column - nPre, row - nPre );

  // two loops over Legendre polynomials in basis
  for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
    for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {
      integrator.setCurrentLegendreOrder( legBasis, legTest );
#pragma omp parallel shared(localTripletList)
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

        // inner two loops over space
#pragma omp for collapse(2)
        for ( LO i = 0; i < iMax; i++ ) {
          for ( LO j = 0; j < jMax; j++ ) {
            integrator.getElemMatrixHypersingular( outerElems[i], innerElems[j],
                *( matrixBuffer[counter] ) );
            this->space->getOuterElemDOFs( outerElems[i],
                &( *( rowIndicesBuffer[counter] ) )[0] );
            this->space->getInnerElemDOFs( innerElems[j],
                &( *( colIndicesBuffer[counter] ) )[0] );
            for ( int k = 0; k < nLocalRows; k++ ) {
#ifndef EXPERIMENTAL_WAVE
              ( *( rowIndicesBuffer[counter] ) )[k] +=
                  row * this->space->getOuterDOFs( ) * ( legOrder + 1 ) +
                  this->space->getOuterDOFs( ) * legTest;
#else
              ( *( rowIndicesBuffer[counter] ) )[k] += ( row ) * this->space->getOuterDOFs( );
#endif
            }
            for ( int k = 0; k < nLocalCols; k++ ) {
#ifndef EXPERIMENTAL_WAVE
              ( *( colIndicesBuffer[counter] ) )[k] += column *
                  this->space->getInnerDOFs( )* ( legOrder + 1 ) +
                  this->space->getInnerDOFs( ) * legBasis;
#else
              ( *( colIndicesBuffer[counter] ) )[k] += ( column ) * this->space->getInnerDOFs( );
#endif
            }
            counter++;
            //matrix.addToPositions( rowIndices, colIndices, elemMatrix );
            if ( counter == BUFFER_SIZE ) {
              counter = 0;
#pragma omp critical
              {
                for ( int m = 0; m < BUFFER_SIZE; m++ ) {
                  for ( int k = 0; k < nLocalRows; k++ ) {
                    for ( int l = 0; l < nLocalCols; l++ ) {
                      if ( fabs( ( *( matrixBuffer[m] ) ).get( k, l ) ) > EPS ) {
                        localTripletList.push_back( Eigen::Triplet<SC, LO>(
                            ( *( rowIndicesBuffer[m] ) )[k], ( *(
                            colIndicesBuffer[m] ) )[l],
                            ( *( matrixBuffer[m] ) ).get( k, l ) ) );
                        //nnz++;
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
                  localTripletList.push_back( Eigen::Triplet<SC, LO>( ( *(
                      rowIndicesBuffer[m] ) )[k],
                      ( *( colIndicesBuffer[m] ) )[l],
                      ( *( matrixBuffer[m] ) ).get( k, l ) ) );
                  //nnz++;
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

  Eigen::SparseMatrix<SC, Eigen::ColMajor, LO> localMatrix( ( nTimeSteps + nPos
      + nPre ) * this->space->getOuterDOFs( ) * ( legOrder + 1 ),
      ( nTimeSteps + nPos + nPre ) * this->space->getInnerDOFs( ) *
      ( legOrder + 1 ) );
  localMatrix.setFromTriplets( localTripletList.begin( ),
      localTripletList.end( ) );
  nnz += localMatrix.nonZeros( );
  for ( LO k = 0; k < localMatrix.outerSize( ); ++k ) {
    typename Eigen::SparseMatrix<SC, Eigen::ColMajor, LO>::InnerIterator it(
        localMatrix, k );
    for (; it; ++it ) {
      tripletList.push_back( Eigen::Triplet<SC, LO>( it.row( ), it.col( ),
          it.value( ) ) );
    }
  }
}

template<class LO, class SC>
void BEBilinearFormWaveHypersingular<LO, SC>::assembleBlock(
    LO row,
    LO column,
    SparseMatrix<LO, SC> &localBlock
    ) const {

  // cast to bespace to bespacetime
  BESpaceTime<LO, SC> *spaceTime =
      static_cast<BESpaceTime<LO, SC>*> ( this->space );
  SurfaceMesh3D<LO, SC> *mesh = this->space->getMesh( );


  // max order of Legendre polynomials in temporal basis
  int legOrder = spaceTime->getLegendreOrder( );

  // set local triplet list
  std::vector<Eigen::Triplet<SC, LO> > localTripletList;

#ifdef EXPERIMENTAL_WAVE
  BEIntegratorWave<LO, SC> integrator( this->space, this->quadratureOrder, this->timeQuadOrder, nPre, nPos );
#else
  BEIntegratorWave<LO, SC> integrator( this->space, this->quadratureOrder,
      this->timeQuadOrder, this->quadratureOrderDisjointElems );
#endif
  const vector<LO> innerElems = this->space->getInnerElems( );
  const vector<LO> outerElems = this->space->getOuterElems( );

  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  LO iMax = outerElems.size( );
  //LO jMax = innerElems.size( );

  integrator.setCurrentFunctions( column, row );

  // precompute nonzero pattern

  std::vector<Eigen::Triplet<SC, LO> > nnzTripletList;
  nnzTripletList.reserve( ( legOrder + 1 )*10 * mesh->getNNodes( ) );

  //std::cout << ( legOrder + 1 )*50 * mesh->getNNodes( ) << std::endl;
  SCVT testStart, testEnd;
  SCVT basisStart, basisEnd;
  SCVT dt = spaceTime->getDt( );
  LO N = spaceTime->getNTimeSteps( );
  if ( row == 0 ) {
    testStart = 0.0;
    testEnd = dt;
  } else if ( row == N - 1 ) {
    testStart = ( N - 2 ) * dt;
    testEnd = ( N - 1 ) * dt;
  } else {
    testStart = ( row - 1 ) * dt;
    testEnd = ( row + 1 ) * dt;
  }
  if ( column == 0 ) {
    basisStart = 0.0;
    basisEnd = dt;
  } else if ( column == N - 1 ) {
    basisStart = ( N - 2 ) * dt;
    basisEnd = ( N - 1 ) * dt;
  } else {
    basisStart = ( column - 1 ) * dt;
    basisEnd = ( column + 1 ) * dt;
  }

  std::vector<std::set<LO> > lightCone( mesh->getNElements( ) );

#pragma omp parallel shared(lightCone, nnzTripletList)
  {
    std::vector<LO> elems1;
    std::vector<LO> elems2;
    std::vector<std::set<LO> > lightConeLocal( mesh->getNElements( ) );
    std::vector<Eigen::Triplet<SC, LO> > nnzTripletListLocal;
    nnzTripletListLocal.reserve( ( legOrder + 1 )*10 * mesh->getNNodes( ) );
    SCVT centroid1[3];
    SCVT centroid2[3];
    SCVT radius1 = 0.0;
    SCVT radius2 = 0.0;
    SCVT x1[3], x2[3], x3[3];
    SCVT dist, distCentr;
    SCVT minDist = 0.0;
    SCVT maxDist = 0.0;
#pragma omp for
    for ( LO i = 0; i < mesh->getNNodes( ); ++i ) {
      mesh->getElements( i, elems1 );
      mesh->getNode( i, centroid1 );
      radius1 = 0.0;
      for ( auto it = elems1.begin( ); it != elems1.end( ); it++ ) {
        mesh->getNodes( *it, x1, x2, x3 );
        dist = DIST3( x1, centroid1 );
        if ( dist > radius1 ) {
          radius1 = dist;
        }
        dist = DIST3( x2, centroid1 );
        if ( dist > radius1 ) {
          radius1 = dist;
        }
        dist = DIST3( x3, centroid1 );
        if ( dist > radius1 ) {
          radius1 = dist;
        }
      }
      for ( LO j = 0; j < mesh->getNNodes( ); ++j ) {
        mesh->getElements( j, elems2 );
        mesh->getNode( j, centroid2 );
        radius2 = 0.0;
        for ( auto it = elems2.begin( ); it != elems2.end( ); it++ ) {
          mesh->getNodes( *it, x1, x2, x3 );
          dist = DIST3( x1, centroid2 );
          if ( dist > radius2 ) {
            radius2 = dist;
          }
          dist = DIST3( x2, centroid2 );
          if ( dist > radius2 ) {
            radius2 = dist;
          }
          dist = DIST3( x3, centroid2 );
          if ( dist > radius2 ) {
            radius2 = dist;
          }
        }
        distCentr = DIST3( centroid1, centroid2 );
        minDist = distCentr - radius1 - radius2;
        maxDist = distCentr + radius1 + radius2;
        if ( ( maxDist > ( testStart - basisEnd ) ) &&
            ( minDist < ( testEnd - basisStart ) ) ) {
          for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
            for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {
              nnzTripletListLocal.push_back( Eigen::Triplet<SC, LO>(
                  legBasis * this->space->getInnerDOFs( ) + i,
                  legTest * this->space->getOuterDOFs( ) + j, 0.0 ) );
            }
          }
          for ( auto it = elems2.begin( ); it != elems2.end( ); it++ ) {
            for ( auto it2 = elems1.begin( ); it2 != elems1.end( ); it2++ ) {
              lightConeLocal[*it2].insert( *it );
            }
          }
        }
      }
    }
#pragma omp critical
    {
      for ( LO i = 0; i < mesh->getNElements( ); i++ ) {
        //        for ( auto it2 = lightConeLocal[i].begin( ); it2 != lightConeLocal[i].end( ); it2++ ) {
        //          lightCone[i].insert( *it2 );
        //        }
        lightCone[i].insert( lightConeLocal[i].begin( ),
            lightConeLocal[i].end( ) );
      }
      for ( auto it2 = nnzTripletListLocal.begin( ); it2 !=
          nnzTripletListLocal.end( ); it2++ ) {
        nnzTripletList.push_back( *it2 );
      }
      //nnzTripletList.insert( nnzTripletList.begin(), nnzTripletListLocal.begin( ), nnzTripletListLocal.end( ) );
    }
  }

  //std::cout << row << " " << column << ": " << nnzTripletList.size( ) << std::endl;


  // preallocate the matrix
  localBlock.resize( ( legOrder + 1 ) * this->space->getOuterDOFs( ),
      ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
  localBlock.getEigenSparseMatrix( )->setFromTriplets( nnzTripletList.begin( ),
      nnzTripletList.end( ) );


  // two loops over Legendre polynomials in basis
  for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
    for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {

      if ( ( legTest > legBasis ) && ( row != 0 )&& ( column != 0 ) &&
          ( row != N - 1 )&&( column != N - 1 ) ) {
        continue;
      }

      integrator.setCurrentLegendreOrder( legBasis, legTest );

      SC multiplier = 0.0;
      if ( ( ( legBasis + legTest ) % 2 ) == 0 ) {
        multiplier = 1.0;
      } else {
        multiplier = -1.0;
      }

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

        // inner two loops over space
#pragma omp for schedule(dynamic)//collapse(2)
        for ( LO i = 0; i < iMax; i++ ) {
          //for ( LO j = 0; j < lightCone[i].size(); j++ ) {
          for ( auto it = lightCone[i].begin( ); it != lightCone[i].end( );
              it++ ) {
            //            std::cout << i << ": " << std::endl;
            //            std::cout << *it << " " << innerElems[*it] << std::endl;
            integrator.getElemMatrixHypersingular( outerElems[i],
                innerElems[*it], *( matrixBuffer[counter] ) );
            this->space->getOuterElemDOFs( outerElems[i],
                &( *( rowIndicesBuffer[counter] ) )[0] );
            this->space->getInnerElemDOFs( innerElems[*it],
                &( *( colIndicesBuffer[counter] ) )[0] );

            // correct for the current Legendre polynomial
            for ( int k = 0; k < nLocalRows; k++ ) {
              ( *( rowIndicesBuffer[counter] ) )[k] += legTest *
                  this->space->getOuterDOFs( );
            }
            for ( int k = 0; k < nLocalCols; k++ ) {
              ( *( colIndicesBuffer[counter] ) )[k] += legBasis *
                  this->space->getInnerDOFs( );
            }

            counter++;
            //matrix.addToPositions( rowIndices, colIndices, elemMatrix );
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
                        if ( legBasis != legTest && ( row != 0 )&&
                            ( column != 0 ) && ( row != N - 1 )&&
                            ( column != N - 1 ) ) {

                          LO lrow = ( *( rowIndicesBuffer[m] ) )[k] -
                              legTest * this->space->getOuterDOFs( ) +
                              legBasis * this->space->getOuterDOFs( );

                          LO lcol = ( *( colIndicesBuffer[m] ) )[l] -
                              legBasis * this->space->getInnerDOFs( ) +
                              legTest * this->space->getInnerDOFs( );


                          localBlock.getEigenSparseMatrix( )->coeffRef(
                              lrow, lcol ) +=
                              multiplier * ( *( matrixBuffer[m] ) ).get( k, l );
                        }
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

                  if ( legBasis != legTest && ( row != 0 )&& ( column != 0 ) &&
                      ( row != N - 1 )&&( column != N - 1 ) ) {
                    LO lrow = ( *( rowIndicesBuffer[m] ) )[k] -
                        legTest * this->space->getOuterDOFs( ) +
                        legBasis * this->space->getOuterDOFs( );

                    LO lcol = ( *( colIndicesBuffer[m] ) )[l] -
                        legBasis * this->space->getInnerDOFs( ) +
                        legTest * this->space->getInnerDOFs( );

                    localBlock.getEigenSparseMatrix( )->coeffRef(
                        lrow, lcol ) +=
                        multiplier * ( *( matrixBuffer[m] ) ).get( k, l );

                  }
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

  // assemble the block from the triplets
  //localBlock.resize( ( legOrder + 1 ) * this->space->getOuterDOFs( ), ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
  //localBlock.getEigenSparseMatrix( )->setFromTriplets( localTripletList.begin( ), localTripletList.end( ) );

}

template<class LO, class SC>
void BEBilinearFormWaveHypersingular<LO, SC>::assembleMPIBlock(
    LO row,
    LO column,
    SparseMatrix<LO, SC> &localBlock,
    int owner
    ) const {

  int rank, nProc;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &nProc );

  MPI_Datatype MPI_LO = GetType<LO, SC>::MPI_LO( );
  MPI_Datatype MPI_SC = GetType<LO, SC>::MPI_SC( );

  // cast the bespace to bespacetime
  BESpaceTime<LO, SC> *spaceTime =
      static_cast<BESpaceTime<LO, SC>*> ( this->space );
  SurfaceMesh3D<LO, SC> *mesh = this->space->getMesh( );


  // max order of Legendre polynomials in temporal basis
  int legOrder = spaceTime->getLegendreOrder( );

  const vector<LO> innerElems = this->space->getInnerElems( );
  const vector<LO> outerElems = this->space->getOuterElems( );

  LO nLocalRows = this->space->getDOFsPerOuterElem( );
  LO nLocalCols = this->space->getDOFsPerInnerElem( );

  //LO iMax = outerElems.size( );
  //LO jMax = innerElems.size( );



  // precompute nonzero pattern

  std::vector<Eigen::Triplet<SC, LO> > nnzTripletList;
  //nnzTripletList.reserve( ( legOrder + 1 )*10 * mesh->getNNodes( ) );

  //std::cout << ( legOrder + 1 )*50 * mesh->getNNodes( ) << std::endl;
  SCVT testStart, testEnd;
  SCVT basisStart, basisEnd;
  SCVT dt = spaceTime->getDt( );
  LO N = spaceTime->getNTimeSteps( );
  if ( row == 0 ) {
    testStart = 0.0;
    testEnd = dt;
  } else if ( row == N - 1 ) {
    testStart = ( N - 2 ) * dt;
    testEnd = ( N - 1 ) * dt;
  } else {
    testStart = ( row - 1 ) * dt;
    testEnd = ( row + 1 ) * dt;
  }
  if ( column == 0 ) {
    basisStart = 0.0;
    basisEnd = dt;
  } else if ( column == N - 1 ) {
    basisStart = ( N - 2 ) * dt;
    basisEnd = ( N - 1 ) * dt;
  } else {
    basisStart = ( column - 1 ) * dt;
    basisEnd = ( column + 1 ) * dt;
  }

  // divide nodes among processors
  LO *rankStarts = new LO[nProc + 1];
  LO nodesPerProc = mesh->getNNodes( ) / nProc;
  LO remainder = mesh->getNNodes( ) % nProc;
  rankStarts[0] = 0;
  for ( int i = 1; i < nProc + 1; i++ ) {
    rankStarts[i] = rankStarts[i - 1] + nodesPerProc;
    if ( remainder > 0 ) {
      rankStarts[i]++;
      remainder--;
    }
  }

  std::vector< std::vector<LO> > processLightCone( mesh->getNElements( ) );
  std::vector< std::set<LO> > lightCone( mesh->getNElements( ) );

  // local (i,j,0.0) on each MPI process to be allreduced
  std::vector< LO > indILocal;
  std::vector< LO > indJLocal;
  indILocal.reserve( ( legOrder + 1 )*10 * mesh->getNNodes( ) );
  indJLocal.reserve( ( legOrder + 1 )*10 * mesh->getNNodes( ) );

#pragma omp parallel shared(processLightCone, nnzTripletList)
  {
    std::vector<LO> elems1;
    std::vector<LO> elems2;
    std::vector<std::set<LO> > lightConeLocal( mesh->getNElements( ) );
    std::vector<Eigen::Triplet<SC, LO> > nnzTripletListLocal;

    SCVT centroid1[3];
    SCVT centroid2[3];
    SCVT radius1 = 0.0;
    SCVT radius2 = 0.0;
    SCVT x1[3], x2[3], x3[3];
    SCVT dist, distCentr;
    SCVT minDist = 0.0;
    SCVT maxDist = 0.0;
#pragma omp for
    for ( LO i = rankStarts[rank]; i < rankStarts[rank + 1]; i++ ) {

      mesh->getElements( i, elems1 );
      mesh->getNode( i, centroid1 );
      radius1 = 0.0;
      for ( auto it = elems1.begin( ); it != elems1.end( ); it++ ) {
        mesh->getNodes( *it, x1, x2, x3 );
        dist = DIST3( x1, centroid1 );
        if ( dist > radius1 ) {
          radius1 = dist;
        }
        dist = DIST3( x2, centroid1 );
        if ( dist > radius1 ) {
          radius1 = dist;
        }
        dist = DIST3( x3, centroid1 );
        if ( dist > radius1 ) {
          radius1 = dist;
        }
      }
      for ( LO j = 0; j < mesh->getNNodes( ); ++j ) {
        mesh->getElements( j, elems2 );
        mesh->getNode( j, centroid2 );
        radius2 = 0.0;
        for ( auto it = elems2.begin( ); it != elems2.end( ); it++ ) {
          mesh->getNodes( *it, x1, x2, x3 );
          dist = DIST3( x1, centroid2 );
          if ( dist > radius2 ) {
            radius2 = dist;
          }
          dist = DIST3( x2, centroid2 );
          if ( dist > radius2 ) {
            radius2 = dist;
          }
          dist = DIST3( x3, centroid2 );
          if ( dist > radius2 ) {
            radius2 = dist;
          }
        }
        distCentr = DIST3( centroid1, centroid2 );
        minDist = distCentr - radius1 - radius2;
        maxDist = distCentr + radius1 + radius2;
        if ( ( maxDist > ( testStart - basisEnd ) ) &&
            ( minDist < ( testEnd - basisStart ) ) ) {
          for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
            for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {
              nnzTripletListLocal.push_back( Eigen::Triplet<SC, LO>(
                  legBasis * this->space->getInnerDOFs( ) + i,
                  legTest * this->space->getOuterDOFs( ) + j, 0.0 ) );
            }
          }
          for ( auto it = elems2.begin( ); it != elems2.end( ); it++ ) {
            for ( auto it2 = elems1.begin( ); it2 != elems1.end( ); it2++ ) {
              lightConeLocal[*it2].insert( *it );
            }
          }
        }
      }
    }
#pragma omp critical
    {
      for ( LO i = 0; i < mesh->getNElements( ); i++ ) {
        //        for ( auto it2 = lightConeLocal[i].begin( ); it2 != lightConeLocal[i].end( ); it2++ ) {
        //          lightCone[i].insert( *it2 );
        //        }
        for ( auto it = lightConeLocal[i].begin( ); it !=
            lightConeLocal[i].end( ); it++ ) {
          typename std::vector<LO>::iterator it2
              = std::lower_bound( processLightCone[i].begin( ),
              processLightCone[i].end( ), *it );
          if ( it2 == processLightCone[i].end( ) || *it < *it2 ) {
            processLightCone[i].insert( it2, *it );
          }
        }
      }
      for ( auto it2 = nnzTripletListLocal.begin( ); it2 !=
          nnzTripletListLocal.end( ); it2++ ) {
        indILocal.push_back( ( *it2 ).row( ) );
        indJLocal.push_back( ( *it2 ).col( ) );
        //nnzTripletList.push_back( *it2 );

      }
      //nnzTripletList.insert( nnzTripletList.begin(), nnzTripletListLocal.begin( ), nnzTripletListLocal.end( ) );
    }
  }

  LO *lightConeSize = new LO[mesh->getNElements( )];
  LO *tmpLightConeSize = new LO[mesh->getNElements( )];
  int totalSize = 0;
  for ( LO i = 0; i < mesh->getNElements( ); i++ ) {
    lightConeSize[i] = processLightCone[i].size( );
    totalSize += lightConeSize[i];
  }

  LO *lsToSend = new LO[totalSize];
  LO offset = 0;
  for ( LO i = 0; i < mesh->getNElements( ); i++ ) {
    if ( lightConeSize[i] > 0 ) {
      memcpy( lsToSend + offset, &( ( processLightCone[i] )[0] ),
          lightConeSize[i] * sizeof (LO ) );
      offset += lightConeSize[i];
    }
  }

  for ( int i = 0; i < nProc; i++ ) {
    int nElems = 0;
    if ( rank == i ) {
      MPI_Bcast( &totalSize, 1, MPI_INT, i, MPI_COMM_WORLD );
      MPI_Bcast( lsToSend, totalSize, MPI_LO, i, MPI_COMM_WORLD );
      MPI_Bcast( lightConeSize, (int) mesh->getNElements( ), MPI_LO, i,
          MPI_COMM_WORLD );

      offset = 0;

      for ( LO j = 0; j < mesh->getNElements( ); j++ ) {
        for ( LO k = 0; k < lightConeSize[j]; k++ ) {
          lightCone[j].insert( lsToSend[offset + k] );
        }
        offset += lightConeSize[j];
      }
    } else {
      MPI_Bcast( &nElems, 1, MPI_INT, i, MPI_COMM_WORLD );
      LO *recvBuff = new LO[nElems];
      MPI_Bcast( recvBuff, nElems, MPI_LO, i, MPI_COMM_WORLD );
      MPI_Bcast( tmpLightConeSize, (int) mesh->getNElements( ), MPI_LO, i,
          MPI_COMM_WORLD );

      offset = 0;

      for ( LO j = 0; j < mesh->getNElements( ); j++ ) {
        for ( LO k = 0; k < tmpLightConeSize[j]; k++ ) {
          lightCone[j].insert( recvBuff[offset + k] );
        }
        offset += tmpLightConeSize[j];
      }
      delete [] recvBuff;
    }
  }

  delete [] lightConeSize;
  // gather data to all processes
  int nLocalTriplets = indILocal.size( );
  int *sizes = new int[nProc];
  int *displ = new int[nProc];
  LO nTriplets = 0;
  MPI_Allgather( &nLocalTriplets, 1, MPI_INT, sizes, 1, MPI_INT,
      MPI_COMM_WORLD );
  for ( LO i = 0; i < nProc; i++ ) {
    nTriplets += sizes[i];
  }

  displ[0] = 0;
  for ( LO i = 1; i < nProc; i++ ) {
    displ[i] = displ[i - 1] + sizes[i - 1];
  }
  LO *indIGlobal = new LO[nTriplets];
  LO *indJGlobal = new LO[nTriplets];
  MPI_Allgatherv( &indILocal[0], sizes[rank], MPI_LO, indIGlobal, sizes,
      displ, MPI_LO, MPI_COMM_WORLD );
  MPI_Allgatherv( &indJLocal[0], sizes[rank], MPI_LO, indJGlobal, sizes,
      displ, MPI_LO, MPI_COMM_WORLD );

  nnzTripletList.reserve( nTriplets );
  for ( LO i = 0; i < nTriplets; i++ ) {
    nnzTripletList.push_back( Eigen::Triplet<SC, LO>(
        indIGlobal[i], indJGlobal[i], 0.0 ) );
  }

  delete [] sizes;
  delete [] displ;
  delete [] indIGlobal;
  delete [] indJGlobal;

  std::vector<std::pair<LO, LO> > elemPairs;
  for ( LO i = 0; i < mesh->getNElements( ); i++ ) {
    for ( auto it = lightCone[i].begin( ); it != lightCone[i].end( ); it++ ) {
      elemPairs.push_back( std::pair<LO, LO>( i, *it ) );
    }
  }


  LO pairsPerProc = elemPairs.size( ) / nProc;
  remainder = elemPairs.size( ) % nProc;
  rankStarts[0] = 0;
  for ( int i = 1; i < nProc + 1; i++ ) {
    rankStarts[i] = rankStarts[i - 1] + pairsPerProc;
    if ( remainder > 0 ) {
      rankStarts[i]++;
      remainder--;
    }
  }

  // preallocate the matrix
  localBlock.resize( ( legOrder + 1 ) * this->space->getOuterDOFs( ),
      ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
  localBlock.getEigenSparseMatrix( )->setFromTriplets( nnzTripletList.begin( ),
      nnzTripletList.end( ) );


  // two loops over Legendre polynomials in basis
  //  for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
  //    for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {
  for ( int legBasis = 0; legBasis < legOrder + 1; legBasis++ ) {
    for ( int legTest = 0; legTest < legOrder + 1; legTest++ ) {

      if ( ( legTest > legBasis ) && ( row != 0 )&& ( column != 0 ) &&
          ( row != N - 1 )&&( column != N - 1 ) ) {
        continue;
      }



      SC multiplier = 0.0;
      if ( ( legBasis + legTest ) % 2 == 0 ) {
        multiplier = 1.0;
      } else {
        multiplier = -1.0;
      }
#pragma omp parallel shared(localBlock)
      {

        BEIntegratorWave<LO, SC> integrator( this->space, this->quadratureOrder,
            this->timeQuadOrder, this->quadratureOrderDisjointElems );
        integrator.setCurrentFunctions( column, row );
        integrator.setCurrentLegendreOrder( legBasis, legTest );

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

        // inner two loops over space
#pragma omp for schedule(dynamic)//collapse(2)
        for ( LO pair = rankStarts[rank]; pair < rankStarts[rank + 1]; pair++ ) {
          LO i = elemPairs[pair].first;
          LO j = elemPairs[pair].second;

          //for ( LO j = 0; j < lightCone[i].size(); j++ ) {
          //for ( auto it = lightCone[i].begin( ); it != lightCone[i].end( ); it++ ) {
          //std::cout << i << ": "<< std::endl;
          //std::cout << *it << " " <<innerElems[*it]  << std::endl;
          integrator.getElemMatrixHypersingular( outerElems[j], innerElems[i],
              *( matrixBuffer[counter] ) );
          this->space->getOuterElemDOFs( outerElems[j],
              &( *( rowIndicesBuffer[counter] ) )[0] );
          this->space->getInnerElemDOFs( innerElems[i],
              &( *( colIndicesBuffer[counter] ) )[0] );

          // correct for the current Legendre polynomial
          for ( int k = 0; k < nLocalRows; k++ ) {
            ( *( rowIndicesBuffer[counter] ) )[k] += legTest *
                this->space->getOuterDOFs( );
          }
          for ( int k = 0; k < nLocalCols; k++ ) {
            ( *( colIndicesBuffer[counter] ) )[k] += legBasis *
                this->space->getInnerDOFs( );
          }


          counter++;
          //matrix.addToPositions( rowIndices, colIndices, elemMatrix );
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
                      if ( legBasis != legTest && ( row != 0 )&&
                          ( column != 0 ) && ( row != N - 1 )&&
                          ( column != N - 1 ) ) {

                        LO lrow = ( *( rowIndicesBuffer[m] ) )[k] -
                            legTest * this->space->getOuterDOFs( ) +
                            legBasis * this->space->getOuterDOFs( );

                        LO lcol = ( *( colIndicesBuffer[m] ) )[l] -
                            legBasis * this->space->getInnerDOFs( ) +
                            legTest * this->space->getInnerDOFs( );


                        localBlock.getEigenSparseMatrix( )->coeffRef(
                            lrow, lcol ) +=
                            multiplier * ( *( matrixBuffer[m] ) ).get( k, l );
                      }
                      //localTripletList.push_back( Eigen::Triplet<SC, LO>( ( *( rowIndicesBuffer[m] ) )[k], ( *( colIndicesBuffer[m] ) )[l], ( *( matrixBuffer[m] ) ).get( k, l ) ) );
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
                  if ( legBasis != legTest && ( row != 0 )&&
                      ( column != 0 ) && ( row != N - 1 )&&
                      ( column != N - 1 ) ) {

                    LO lrow = ( *( rowIndicesBuffer[m] ) )[k] -
                        legTest * this->space->getOuterDOFs( ) +
                        legBasis * this->space->getOuterDOFs( );

                    LO lcol = ( *( colIndicesBuffer[m] ) )[l] -
                        legBasis * this->space->getInnerDOFs( ) +
                        legTest * this->space->getInnerDOFs( );


                    localBlock.getEigenSparseMatrix( )->coeffRef(
                        lrow, lcol ) +=
                        multiplier * ( *( matrixBuffer[m] ) ).get( k, l );
                  }
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

  if ( rank != owner ) {
    //localBlock.getEigenSparseMatrix( )->prune( 0.0 );
    // now send every part to the block's owner 
    LO nnz = localBlock.getEigenSparseMatrix( )->nonZeros( );
    //    LO *rows = new LO[nnz];
    //    LO *cols = new LO[nnz];
    //    SCVT *vals = new SCVT[nnz];
    //    LO counter = 0;
    //std::cout << nnz << std::endl;

    //    for ( int k = 0; k < localBlock.getEigenSparseMatrix( )->outerSize( ); ++k ) {
    //      for ( typename Eigen::SparseMatrix<SC>::InnerIterator it(
    //          *localBlock.getEigenSparseMatrix( ), k ); it; ++it ) {
    //        //std::cout << counter << std::endl;
    //        vals[counter] = it.value( );
    //        rows[counter] = it.row( );
    //        cols[counter] = it.col( );
    //        counter++;
    //      }
    //    }
    //    MPI_Send( &counter, 1, MPI_INT, owner, 0, MPI_COMM_WORLD );
    //    MPI_Send( vals, counter, MPI_DOUBLE, owner, 0, MPI_COMM_WORLD );
    //    MPI_Send( rows, counter, MPI_INT, owner, 0, MPI_COMM_WORLD );
    //    MPI_Send( cols, counter, MPI_INT, owner, 0, MPI_COMM_WORLD );

    //    delete [] rows;
    //    delete [] cols;
    //    delete [] vals;

    MPI_Reduce( localBlock.getEigenSparseMatrix( )->valuePtr( ), nullptr,
        (int) nnz, MPI_SC, MPI_SUM, owner, MPI_COMM_WORLD );

  } else {
    LO nnz = localBlock.getEigenSparseMatrix( )->nonZeros( );

    MPI_Reduce( MPI_IN_PLACE, localBlock.getEigenSparseMatrix( )->valuePtr( ),
        (int) nnz, MPI_SC, MPI_SUM, owner, MPI_COMM_WORLD );

    //    MPI_Status status;
    //    for ( int i = 0; i < nProc; i++ ) {
    //      if ( i != owner ) {
    //        MPI_Recv( &nnz, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
    //        LO *rows = new LO[nnz];
    //        LO *cols = new LO[nnz];
    //        SCVT *vals = new SCVT[nnz];
    //        MPI_Recv( vals, nnz, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
    //        MPI_Recv( rows, nnz, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
    //        MPI_Recv( cols, nnz, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
    //
    //        for ( LO j = 0; j < nnz; j++ ) {
    //          localBlock.getEigenSparseMatrix( )->coeffRef(
    //              rows[j], cols[j] ) += vals[j];
    //        }
    //        delete [] rows;
    //        delete [] cols;
    //        delete [] vals;
    //      }
    //    }
    localBlock.getEigenSparseMatrix( )->prune( 0.0 );
  }
  //   assemble the block from the triplets
  //  localBlock.resize( ( legOrder + 1 ) * this->space->getOuterDOFs( ), ( legOrder + 1 ) * this->space->getInnerDOFs( ) );
  //  localBlock.getEigenSparseMatrix( )->setFromTriplets( localTripletList.begin( ), localTripletList.end( ) );


  delete [] rankStarts;

}


}



#else 
// use experimental basis functions which are not a partition of unity

namespace bem4i {

template<class LO, class SC>
void BEBilinearFormWaveHypersingular<LO, SC>::assemble( SparseMatrix<LO, SC>& matrix ) const {


  BESpaceTime<LO, SC> *spaceTime = static_cast<BESpaceTime<LO, SC>*> ( this->space );

  int nTimeSteps = spaceTime->getNTimeSteps( );
  LO nnz;
  LO length = 0;

  matrix.resize( ( nTimeSteps + nPos + nPre ) * this->space->getOuterDOFs( ),
      ( nTimeSteps + nPos + nPre ) * this->space->getInnerDOFs( ) );

  // create Eigen triplets and insert values into matrix (duplicated indices are summed!)
  std::vector<Eigen::Triplet<SC, LO> > tripletList;
  tripletList.reserve( ( nTimeSteps + nPos + nPre ) * this->space->getInnerDOFs( ) * this->space->getInnerDOFs( ) );


  // at first, assemble the first column and copy it to the rest of the matrix
  for ( int i = 0; i < nTimeSteps + nPos + nPre; i++ ) {
    assembleBlock( i, 0, tripletList, nnz );
    length = tripletList.size( );


    for ( int j = 1; j < nTimeSteps + nPos + nPre - i; j++ ) {

      for ( int k = 0; k < nnz; k++ ) {
        Eigen::Triplet<SC, LO> newTriplet(
            tripletList[length - nnz + k].row( ) + j * this->space->getOuterDOFs( ),
            tripletList[length - nnz + k].col( ) + j * this->space->getInnerDOFs( ),
            tripletList[length - nnz + k].value( ) );
        tripletList.push_back( newTriplet );

      }
    }
  }

  // block on position (0,1) is copied on the upper diagonal
  assembleBlock( 0, 1, tripletList, nnz );
  length = tripletList.size( );
  for ( int j = 1; j < nTimeSteps + nPos + nPre - 1; j++ ) {
    for ( int k = 0; k < nnz; k++ ) {
      Eigen::Triplet<SC, LO> newTriplet(
          tripletList[length - nnz + k].row( ) + j * this->space->getOuterDOFs( ),
          tripletList[length - nnz + k].col( ) + j * this->space->getInnerDOFs( ),
          tripletList[length - nnz + k].value( ) );
      tripletList.push_back( newTriplet );

    }
  }
  matrix.getEigenSparseMatrix( )->setFromTriplets( tripletList.begin( ), tripletList.end( ) );
  matrix.makeCompressed( );


  //  Eigen::MatrixXd m;
  //  m = (*(matrix.getEigenSparseMatrix()));
  //  std::cout << m << std::endl;


  //
  //  BESpaceTime<LO, SC> *spaceTime = static_cast<BESpaceTime<LO, SC>*> ( this->space );
  //
  //  int nTimeSteps = spaceTime->getNTimeSteps( );
  //  LO nnz;
  //  LO length = 0;
  //
  //  matrix.resize( ( nTimeSteps + nPos ) * this->space->getOuterDOFs( ),
  //      ( nTimeSteps + nPos ) * this->space->getInnerDOFs( ) );
  //
  //  // create Eigen triplets and insert values into matrix (duplicated indices are summed!)
  //  std::vector<Eigen::Triplet<SC, LO> > tripletList;
  //  tripletList.reserve( ( nTimeSteps + nPos ) * this->space->getInnerDOFs( ) * this->space->getInnerDOFs( ) );
  //
  //  // at first, assemble the first column and copy it to the rest of the matrix
  //  for ( int i = 1; i < nTimeSteps + nPos; i++ ) {
  //    assembleBlock( i, 1, tripletList, nnz );
  //    length = tripletList.size( );
  //    for ( int j = 1; j < nTimeSteps + nPos - i; j++ ) {
  //
  //      for ( int k = 0; k < nnz; k++ ) {
  //        Eigen::Triplet<SC, LO> newTriplet(
  //            tripletList[length - nnz + k].row( ) + j * this->space->getOuterDOFs( ),
  //            tripletList[length - nnz + k].col( ) + j * this->space->getInnerDOFs( ),
  //            tripletList[length - nnz + k].value( ) );
  //        tripletList.push_back( newTriplet );
  //      }
  //    }
  //  }
  //
  //  // block on position (0,1) is copied on the upper diagonal
  //  assembleBlock( 1, 2, tripletList, nnz );
  //  length = tripletList.size( );
  //  for ( int j = 1; j < nTimeSteps + nPos - 2; j++ ) {
  //    for ( int k = 0; k < nnz; k++ ) {
  //      Eigen::Triplet<SC, LO> newTriplet(
  //          tripletList[length - nnz + k].row( ) + j * this->space->getOuterDOFs( ),
  //          tripletList[length - nnz + k].col( ) + j * this->space->getInnerDOFs( ),
  //          tripletList[length - nnz + k].value( ) );
  //      tripletList.push_back( newTriplet );
  //    }
  //  }
  //
  //  matrix.getEigenSparseMatrix( )->setFromTriplets( tripletList.begin( ), tripletList.end( ) );
  //  matrix.makeCompressed( );
}

}
#endif


#endif
