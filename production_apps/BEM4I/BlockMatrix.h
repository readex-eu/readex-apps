/*!
 * @file    BlockMatrix.h
 * @author  Michal Merta 
 * @date    November 25, 2013
 * @brief   Header file for the BlockMatrix class
 * 
 */

#ifndef BLOCKMATRIX_H
#define	BLOCKMATRIX_H

// for some reason this has to be included bcs of SparseMatrix compilation
#include "FullMatrix.h"  
#include "Matrix.h"
#include "Vector.h"

namespace bem4i {

/*! 
 * Class representing a block matrix
 * 
 * 
 */
template<class LO, class SC>
class BlockMatrix : public Matrix<LO, SC> {
public:

  //! default constructor
  BlockMatrix( );

  //! copy constructor
  BlockMatrix(
      const BlockMatrix& orig
      );

  /*! Consturctor taking number of blocks and their sizes as arguments
   * 
   * @param[in]   nBlockRows  number of block row
   * @param[in]   nBlockCols  number of block columns
   * @param[in]   numbersOfRows array of legth nBlockRows containing numbers of rows of each block row
   * @param[in]   numbersOfCols array of legth nBlockCols containing numbers of columns of each block column
   */
  BlockMatrix(
      int nBlockRows,
      int nBlockCols,
      LO *numbersOfRows,
      LO *numbersOfCols
      );

  //! destructor
  virtual ~BlockMatrix( );

  //! returns total number of blocks

  inline int getNBlocks( ) {
    return nBlockRows * nBlockCols;
  }

  //! returns number of row blocks

  inline int getNBlockRows( ) const {
    return nBlockRows;
  }

  //! returns number of column blocks

  inline int getNBlockCols( ) const {
    return nBlockCols;
  }

  //! returns block on position (i, j)

  inline Matrix<LO, SC>* getBlock(
      LO i,
      LO j
      ) const {
    return this->blocks[j * nBlockRows + i];
  }

  /*! Inserts a block to a given location of the matrix
   * 
   * @param[in]   i  row
   * @param[in]   j  column
   * @param[in]   block block of the matrix to insert (can be nullptr - considered as a zero matrix)
   * @param[in]   del bool flag indicating whether the matrix is responsible for the deletion of a given block
   */
  inline void setBlock(
      LO i,
      LO j,
      Matrix<LO, SC>* block,
      bool del = false
      ) {
    blocks[j * nBlockRows + i] = block;
    delBlocks[j * nBlockRows + i] = del;
  }


  /*!
   * @brief Performs a matrix-vector multiplication
   * 
   * Computes y = beta*y + alpha*this*x
   * @param A 
   * @param x
   * @param y
   * @param alpha
   * @param beta
   */
  virtual void apply(
      Vector<LO, SC> const & x,
      Vector<LO, SC> & y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

  /*! Method resizes a matrix according to the given parameters
   * 
   * @param[in]   nBlockRows  number of block row
   * @param[in]   nBlockCols  number of block columns
   * @param[in]   numbersOfRows array of legth nBlockRows containing numbers of rows of each block row
   * @param[in]   numbersOfCols array of legth nBlockCols containing numbers of columns of each block column
   */
  virtual void resize(
      int nBlockRows,
      int nBlockCols,
      LO *numbersOfRows,
      LO *numbersOfCols
      );

  LO getNRowsOfBlock( LO i ) {
    return numbersOfRows[i];
  }

  LO getNColsOfBlock( LO i ) {
    return numbersOfCols[i];
  }

  void print( std::ostream &stream = std::cout ) const {
    for ( LO i = 0; i < nBlockCols; i++ ) {
      for ( LO j = 0; j < nBlockRows; j++ ) {
        stream << std::endl << "Block (" << j << ", " << i << "): "
            << std::endl;
        if ( getBlock( j, i ) != nullptr ) {
          stream << "Local rows: " << getBlock( j, i )->getNRows( )
              << std::endl;
          stream << "Local columns: " << getBlock( j, i )->getNCols( )
              << std::endl;
          getBlock( j, i )->print( stream );
        } else {
          stream << "Local rows: 0" << std::endl;
          stream << "Local columns: 0" << std::endl;
        }
      }
    }
    stream << std::endl;
  }

protected:

  //! number of block rows
  int nBlockRows;

  //! number of block columns
  int nBlockCols;

  //! 1D array of pointers to matrix blocks (columnwise)
  Matrix<LO, SC>** blocks;

  //! number of rows in each block
  LO * numbersOfRows;

  //! number of cols in each block
  LO * numbersOfCols;

  //! array of booleans indicating the blocks which matrix has to delete
  bool * delBlocks;

  LO getMaxRows( ) {
    LO maxRows = 0;
    if ( numbersOfRows != nullptr ) {
      for ( LO i = 0; i < nBlockRows; ++i ) {
        if ( numbersOfRows[i] > maxRows ) {
          maxRows = numbersOfRows[i];
        }
      }
    }
    return maxRows;
  }

  LO getMaxCols( ) {
    LO maxCols = 0;
    if ( numbersOfCols != nullptr ) {
      for ( LO i = 0; i < nBlockCols; ++i ) {
        if ( numbersOfCols[i] > maxCols ) {
          maxCols = numbersOfCols[i];
        }
      }
    }
    return maxCols;
  }

};

} // end of namespace bem4i

// include .cpp file to overcome linking problems due to templates
#include "BlockMatrix.cpp"

#endif	/* BLOCKMATRIX_H */

