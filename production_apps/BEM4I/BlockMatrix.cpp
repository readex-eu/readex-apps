/*!
 * @file    BlockMatrix.cpp
 * @author  Michal Merta 
 * @date    November 25, 2013
 * 
 */

#ifdef BLOCKMATRIX_H

namespace bem4i {

template<class LO, class SC>
BlockMatrix<LO, SC>::BlockMatrix( ) {
  this->blocks = nullptr;
  this->numbersOfCols = nullptr;
  this->numbersOfRows = nullptr;
  this->delBlocks = nullptr;
}

template<class LO, class SC>
BlockMatrix<LO, SC>::BlockMatrix(
    const BlockMatrix& orig
    ) {

  this->nBlockRows = orig.getNBlockRows( );
  this->nBlockCols = orig.getNBlockCols( );

  this->numbersOfRows = new LO[ nBlockRows ];
  this->numbersOfCols = new LO[ nBlockCols ];
  this->delBlocks = new bool[ nBlockRows * nBlockCols ];

  memcpy( this->delBlocks, orig.delBlocks,
      nBlockRows * nBlockCols * sizeof ( bool ) );
  memcpy( this->numbersOfRows, orig.numbersOfRows, nBlockRows * sizeof ( LO ) );
  memcpy( this->numbersOfCols, orig.numbersOfCols, nBlockCols * sizeof ( LO ) );

  // perform shallow copy of data
  this->blocks = new Matrix<LO, SC>*[nBlockRows * nBlockCols];
  for ( int i = 0; i < nBlockRows * nBlockCols; i++ ) {
    this->blocks[i] = orig.blocks[i];
    this->delBlocks[i] = false;
  }

}

template<class LO, class SC>
BlockMatrix<LO, SC>::BlockMatrix(
    int nBlockRows,
    int nBlockCols,
    LO * numbersOfRows,
    LO * numbersOfCols
    ) {
  // allocate memory for pointers
  this->blocks = new Matrix<LO, SC>*[nBlockRows * nBlockCols];
  // set all pointers to nullptr
  for ( int i = 0; i < nBlockCols * nBlockRows; i++ ) {
    this->blocks[i] = nullptr;
  }
  this->numbersOfRows = new LO[nBlockRows];
  this->numbersOfCols = new LO[nBlockCols];
  this->delBlocks = new bool[nBlockRows * nBlockCols];
  this->nBlockRows = nBlockRows;
  this->nBlockCols = nBlockCols;

  memset( this->delBlocks, 0, nBlockRows * nBlockCols * sizeof ( bool ) );

  // copy numbers of rows and columns of each block
  memcpy( this->numbersOfRows, numbersOfRows, nBlockRows * sizeof ( LO ) );
  memcpy( this->numbersOfCols, numbersOfCols, nBlockCols * sizeof ( LO ) );

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
BlockMatrix<LO, SC>::~BlockMatrix( ) {

  if ( numbersOfRows ) {
    delete [] numbersOfRows;
  }
  if ( numbersOfCols ) {
    delete [] numbersOfCols;
  }

  if ( delBlocks ) {
    for ( int i = 0; i< this->getNBlocks( ); i++ ) {
      if ( delBlocks[i] == true ) {
        delete blocks[i];
      }
    }
  }

  if ( delBlocks ) {
    delete [] delBlocks;
  }

  if ( blocks ) {
    delete [] blocks;
  }

  this->blocks = nullptr;
  this->numbersOfCols = nullptr;
  this->numbersOfRows = nullptr;
  this->delBlocks = nullptr;

}

template<class LO, class SC>
void BlockMatrix<LO, SC>::apply(
    Vector<LO, SC> const & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  // iterate through matrix blocks and multiply
  Vector<LO, SC> *yi;
  Vector<LO, SC> *xj;
  LO yPos = 0;
  LO xPos = 0;

  if ( !transA ) {
    for ( int i = 0; i < nBlockRows; i++ ) {
      yi = new Vector<LO, SC>( numbersOfRows[i], y.getData( ) + yPos, false );
      yi->scale( beta );
      yPos += numbersOfRows[i];
      xPos = 0;
      for ( int j = 0; j < nBlockCols; j++ ) {
        xj = new Vector<LO, SC>( numbersOfCols[j], x.getData( ) + xPos, false );
        if ( this->getBlock( i, j ) ) {
          this->getBlock( i, j )->apply( *xj, *yi, transA, alpha, 1.0 );
        }
        xPos += numbersOfCols[j];
        delete xj;
      }
      delete yi;
    }
  } else {
    for ( int i = 0; i < nBlockCols; i++ ) {
      yi = new Vector<LO, SC>( numbersOfCols[i], y.getData( ) + yPos, false );
      yi->scale( beta );
      yPos += numbersOfCols[i];
      xPos = 0;
      for ( int j = 0; j < nBlockRows; j++ ) {
        xj = new Vector<LO, SC>( numbersOfRows[j], x.getData( ) + xPos, false );
        if ( this->getBlock( j, i ) ) {
          this->getBlock( j, i )->apply( *xj, *yi, transA, alpha, 1.0 );
        }
        xPos += numbersOfRows[j];
        delete xj;
      }
      delete yi;
    }
  }

}

template<class LO, class SC>
void BlockMatrix<LO, SC>::resize(
    int nBlockRows,
    int nBlockCols,
    LO * numbersOfRows,
    LO * numbersOfCols
    ) {

  if ( this->numbersOfRows ) {
    delete [] numbersOfRows;
  }
  if ( this->numbersOfCols ) {
    delete [] numbersOfCols;
  }

  if ( this->delBlocks ) {
    for ( int i = 0; i< this->getNBlocks( ); i++ ) {
      if ( delBlocks[i] == true ) {
        delete blocks[i];
      }
    }
    delete [] delBlocks;
    delete [] blocks;
  }

  // allocate memory for pointers
  this->blocks = new Matrix<LO, SC>*[nBlockRows * nBlockCols];
  this->numbersOfRows = new LO[nBlockRows];
  this->numbersOfCols = new LO[nBlockCols];
  this->delBlocks = new bool[nBlockRows * nBlockCols];
  this->nBlockRows = nBlockRows;
  this->nBlockCols = nBlockCols;
  memset( this->delBlocks, 0, nBlockRows * nBlockCols * sizeof ( bool ) );

  // copy numbers of rows and columns of each block
  memcpy( this->numbersOfRows, numbersOfRows, nBlockRows * sizeof ( LO ) );
  memcpy( this->numbersOfCols, numbersOfCols, nBlockCols * sizeof ( LO ) );

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

}
#endif
