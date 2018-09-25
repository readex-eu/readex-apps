/*!
 * @file    MPIBlockMatrix.cpp
 * @author  Michal Merta 
 * @date    November 25, 2013
 * 
 */

#ifdef MPIBLOCKMATRIX_H

namespace bem4i {

template<class LO, class SC>
MPIBlockMatrix<LO, SC>::MPIBlockMatrix( ) {
  this->blocks = nullptr;
  this->numbersOfCols = nullptr;
  this->numbersOfRows = nullptr;
  this->delBlocks = nullptr;
  this->ranks = nullptr;
  this->communicator = MPI_COMM_WORLD;
}

template<class LO, class SC>
MPIBlockMatrix<LO, SC>::MPIBlockMatrix(
    int nBlockRows,
    int nBlockCols,
    LO * numbersOfRows,
    LO * numbersOfCols,
    int * ranks,
    MPI_Comm communicator
    ) {

  // allocate memory for pointers
  this->blocks = new Matrix<LO, SC>*[nBlockRows * nBlockCols];
  this->numbersOfRows = new LO[nBlockRows];
  this->numbersOfCols = new LO[nBlockCols];
  this->delBlocks = new bool[nBlockRows * nBlockCols];
  this->ranks = new int[nBlockRows * nBlockCols];
  this->nBlockRows = nBlockRows;
  this->nBlockCols = nBlockCols;
  this->communicator = communicator;
  memset( this->delBlocks, 0, nBlockRows * nBlockCols * sizeof ( bool ) );

  // copy numbers of rows and columns of each block
  memcpy( this->numbersOfRows, numbersOfRows, nBlockRows * sizeof ( LO ) );
  memcpy( this->numbersOfCols, numbersOfCols, nBlockCols * sizeof ( LO ) );

  memcpy( this->ranks, ranks, nBlockRows * nBlockCols * sizeof ( int ) );

  for ( int i = 0; i < this->nBlockRows * this->nBlockCols; i++ ) {
    this->blocks[i] = nullptr;
  }

  this->nRows = 0;
  this->nCols = 0;
  for ( int i = 0; i < this->nBlockRows; i++ ) {
    this->nRows += this->numbersOfRows[i];
  }

  for ( int i = 0; i < this->nBlockCols; i++ ) {
    this->nCols += this->numbersOfCols[i];
  }

}

template<class LO, class SC>
MPIBlockMatrix<LO, SC>::~MPIBlockMatrix( ) {

  if ( ranks ) {
    delete [] ranks;
  }

  this->ranks = nullptr;
}

template<class LO, class SC>
void MPIBlockMatrix<LO, SC>::apply(
    Vector<LO, SC> const & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  MPI_Barrier( communicator );

  // iterate through matrix blocks and multiply
  Vector<LO, SC> *yi;
  Vector<LO, SC> *xj;
  LO yLength = transA ? this->nCols : this->nRows;

  Vector<LO, SC> *localY = new Vector<LO, SC>( yLength );
  Vector<LO, SC> *globalY = new Vector<LO, SC>( yLength );

  localY->setAll( 0.0 );
  LO yPos = 0;
  LO xPos = 0;

  if ( !transA ) {
    for ( int i = 0; i < this->nBlockRows; i++ ) {
      yi = new Vector<LO, SC>( this->numbersOfRows[i],
          localY->getData( ) + yPos, false );
      yPos += this->numbersOfRows[i];
      xPos = 0;
      for ( int j = 0; j < this->nBlockCols; j++ ) {
        if ( this->getBlock( i, j ) && amIOwner( i, j ) ) {
          xj = new Vector<LO, SC>( this->numbersOfCols[j], x.getData( ) + xPos,
              false );
          this->getBlock( i, j )->apply( *xj, *yi, transA, alpha, 1.0 );
          delete xj;
        }
        xPos += this->numbersOfCols[j];

      }
      delete yi;
    }
  } else {
    for ( int i = 0; i < this->nBlockCols; i++ ) {
      yi = new Vector<LO, SC>( this->numbersOfCols[i],
          localY->getData( ) + yPos, false );
      yPos += this->numbersOfCols[i];
      xPos = 0;
      for ( int j = 0; j < this->nBlockRows; j++ ) {
        if ( this->getBlock( j, i ) && amIOwner( j, i ) ) {
          xj = new Vector<LO, SC>( this->numbersOfRows[j], x.getData( ) + xPos,
              false );
          this->getBlock( j, i )->apply( *xj, *yi, transA, alpha, 1.0 );
          delete xj;
        }
        xPos += this->numbersOfRows[j];

      }
      delete yi;
    }
  }
  MPI_Barrier( communicator );
  MPI_Allreduce( localY->getData( ), globalY->getData( ), yLength, 
      GetType<LO, SC>::MPI_SC(), MPI_SUM, communicator );
  
  y.scale( beta );
  y.add( *globalY );

  delete localY;
  delete globalY;

}

template<class LO, class SC>
void MPIBlockMatrix<LO, SC>::resize(
    int nBlockRows,
    int nBlockCols,
    LO * numbersOfRows,
    LO * numbersOfCols,
    int * ranks ) {

  if ( this->ranks ) {
    delete [] ranks;
  }

  if ( this->numbersOfRows ) {
    delete [] numbersOfRows;
  }
  if ( this->numbersOfCols ) {
    delete [] numbersOfCols;
  }

  if ( this->delBlocks ) {
    for ( int i = 0; i< this->getNBlocks( ); i++ ) {
      if ( this->delBlocks[i] == true ) {
        delete this->blocks[i];
      }
    }
    delete [] this->delBlocks;
    delete [] this->blocks;
  }

  // allocate memory for pointers
  this->blocks = new Matrix<LO, SC>*[nBlockRows * nBlockCols];
  this->numbersOfRows = new LO[nBlockRows];
  this->numbersOfCols = new LO[nBlockCols];
  this->ranks = new int[nBlockRows * nBlockCols ];
  this->delBlocks = new bool[nBlockRows * nBlockCols];
  this->nBlockRows = nBlockRows;
  this->nBlockCols = nBlockCols;
  memset( this->delBlocks, 0, nBlockRows * nBlockCols * sizeof ( bool ) );

  // copy numbers of rows and columns of each block
  memcpy( this->numbersOfRows, numbersOfRows, nBlockRows * sizeof ( LO ) );
  memcpy( this->numbersOfCols, numbersOfCols, nBlockCols * sizeof ( LO ) );
  memcpy( this->ranks, ranks, nBlockRows * nBlockCols * sizeof ( int ) );

  this->nRows = 0;
  this->nCols = 0;
  for ( int i = 0; i < this->nBlockRows; i++ ) {
    this->nRows += this->numbersOfRows[i];
  }

  for ( int i = 0; i < this->nBlockCols; i++ ) {
    this->nCols += this->numbersOfCols[i];
  }

  for ( int i = 0; i < nBlockCols * nBlockRows; i++ ) {
    this->blocks[i] = nullptr;
  }

}

template<class LO, class SC>
void MPIBlockMatrix<LO, SC>::resize(
    int nBlockRows,
    int nBlockCols,
    LO * numbersOfRows,
    LO * numbersOfCols ) {

  int * ranks = new int[nBlockRows * nBlockCols];
  memset( ranks, 0, nBlockRows * nBlockCols * sizeof ( int ) );
  this->resize( nBlockRows, nBlockCols, numbersOfRows, numbersOfCols, ranks );
  delete [] ranks;
}

}

#endif
