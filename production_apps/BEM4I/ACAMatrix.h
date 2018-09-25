/*!
 * @file    ACAKernel.h
 * @author  Dalibor Lukas
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    June 16, 2014
 * @brief   Header file for class ACAMatrix
 * 
 */

#ifndef ACAMATRIX_H
#define	ACAMATRIX_H

#include "FullMatrix.h"
#include "FastBESpace.h"

namespace bem4i {

template<class LO, class SC>
class ACAMatrix : public Matrix<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! default constructor
  ACAMatrix( );

  /*!
   * Constructor allocating a full matrix
   * 
   * @param[in]   nRows number of rows
   * @param[in]   nCols number of columns
   */
  ACAMatrix( LO nRows, LO nCols );



  //! destructor
  virtual ~ACAMatrix( );

  //! adds full nonadmissible block to the list of nonadmissible blocks
/*
  inline void addNonadmissibleBlock(
      FullMatrix<LO, SC>* block
      ) {
    nonAdmissibleBlocks.push_back( block );
#pragma omp atomic update
    nonadmBlocksSize += block->getNCols( ) * block->getNRows( );
  }
*/
  //! adds full nonadmissible block associated with idx-th nonadmissible leaf to the list of nonadmissible blocks

  inline void addNonadmissibleBlock(
      FullMatrix<LO, SC>* block,
      LO idx
      ) {
    nonAdmissibleBlocks.at( idx ) = block;
#pragma omp atomic update
    nonadmBlocksSize += block->getNCols( ) * block->getNRows( );
  }

  //! returns pointer to a nonadmissible matrix with index idx

  inline FullMatrix<LO, SC>* getNonAdmissibleBlock(
      LO idx
      ) {
    return nonAdmissibleBlocks[idx];
  }

  //! adds list of nonadmissible cluster pairs

  inline void setNonadmissibleDOFs(
      std::vector< BEBlockCluster< LO, SC > * > & leaves
      ) {

    this->nonadmissibleInnerDOFs.reserve( leaves.size( ) );
    this->nonadmissibleOuterDOFs.reserve( leaves.size( ) );

    for ( LO i = 0; i < leaves.size( ); ++i ) {
      this->nonadmissibleInnerDOFs.push_back(
          new std::vector< LO >( *( leaves[ i ]->innerDOFs ) ) );
      if ( leaves.at( i )->innerDOFs->size( ) > this->maxBlockSize ) {
        this->maxBlockSize = leaves.at( i )->innerDOFs->size( );
      }
      this->nonadmissibleOuterDOFs.push_back(
          new std::vector< LO >( *( leaves[ i ]->outerDOFs ) ) );
      if ( leaves.at( i )->outerDOFs->size( ) > this->maxBlockSize ) {
        this->maxBlockSize = leaves.at( i )->outerDOFs->size( );
      }
    }
  }

  //! adds list of admissible cluster pairs

  inline void setAdmissibleDOFs(
      std::vector< BEBlockCluster< LO, SC > * > & leaves
      ) {

    this->admissibleInnerDOFs.reserve( leaves.size( ) );
    this->admissibleOuterDOFs.reserve( leaves.size( ) );

    for ( LO i = 0; i < leaves.size( ); ++i ) {
      this->admissibleInnerDOFs.push_back(
          new std::vector< LO >( *( leaves[ i ]->innerDOFs ) ) );
      if ( leaves.at( i )->innerDOFs->size( ) > this->maxBlockSize ) {
        this->maxBlockSize = leaves.at( i )->innerDOFs->size( );
      }
      this->admissibleOuterDOFs.push_back(
          new std::vector< LO >( *( leaves[ i ]->outerDOFs ) ) );
      if ( leaves.at( i )->outerDOFs->size( ) > this->maxBlockSize ) {
        this->maxBlockSize = leaves.at( i )->outerDOFs->size( );
      }
    }
  }

  //! adds a pair of U, V matrices from ACA 

  inline void addAdmissibleBlock(
      std::pair<FullMatrix<LO, SC>*, FullMatrix<LO, SC>* > blockUV,
      LO idx
      ) {
    admissibleBlocks[ idx ] = blockUV;

    if ( blockUV.first ) {
#pragma omp atomic update
      admBlocksSize += blockUV.first->getNRows( ) * blockUV.first->getNCols( );
    }
    if ( blockUV.second ) {
#pragma omp atomic update
      admBlocksSize += blockUV.second->getNRows( )
          * blockUV.second->getNCols( );
    }
  }

  //! returns pointer to a pair of matrices U,V with index idx

  inline std::pair<FullMatrix<LO, SC>*, FullMatrix<LO, SC>* >
  getAdmissibleBlock(
      LO idx
      ) {
    return admissibleBlocks[idx];
  }


  //! applies matrix to a vector
  virtual void apply(
      Vector<LO, SC> const &x,
      Vector<LO, SC> &y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

  inline SCVT getCompressionRatio( ) {

    LO nCols = this->getNCols( );
    LO nRows = this->getNRows( );

    if ( nCols * nRows == 0 ) return 1.0;

    if ( this->p12p1disMat ) nCols = this->p12p1disMat->getNRows( );
    if ( this->p1dis2p1Mat ) nRows = this->p1dis2p1Mat->getNCols( );

    return (SCVT) ( ( nonadmBlocksSize + admBlocksSize ) /
        ( (SCVT) nCols * nRows ) );
  }

  inline void resizeNonAdmBlocks(
      LO size
      ) {
    this->nonAdmissibleBlocks.resize( size );
  }

  inline void resizeAdmBlocks(
      LO size
      ) {
    this->admissibleBlocks.resize( size );
  }

  void print(
      std::ostream &stream = std::cout
      ) const {
    std::cout << "ACA Matrix\n";
    std::cout << "Number of rows: " << this->nRows << std::endl;
    std::cout << "Number of cols: " << this->nCols << std::endl;
  };

  void setP12p1dis(
      bool p12p1dis
      ) {
    this->p12p1dis = p12p1dis;
  }

  void setP1dis2p1(
      bool p1dis2p1
      ) {
    this->p1dis2p1 = p1dis2p1;
  }

  bool getP12p1dis(
      ) {
    return this->p12p1dis;
  }

  bool getP1dis2p1(
      ) {
    return this->p1dis2p1;
  }

  void setP12p1disMatFromTriplets(
      LO nRows,
      LO nCols,
      std::vector<LO> & vecI,
      std::vector<LO> & vecJ,
      std::vector<SC> & vecV
      );

  void setP1dis2p1MatFromTriplets(
      LO nRows,
      LO nCols,
      std::vector<LO> & vecI,
      std::vector<LO> & vecJ,
      std::vector<SC> & vecV
      );

protected:

  //! full size of approximated blocks
  LO nonadmBlocksSize;

  //! size of approximated blocks
  LO admBlocksSize;

  //! vector of nonadmissible leaves
  //std::vector<BEBlockCluster<LO, SC>*> nonadmissibleLeaves;

  //! vector of admissible leaves
  //std::vector<BEBlockCluster<LO, SC>*> admissibleLeaves;

  //! vector of nonadmissible matrix blocks
  std::vector<FullMatrix<LO, SC>* > nonAdmissibleBlocks;

  std::vector<std::vector<LO>*> admissibleInnerDOFs;

  std::vector<std::vector<LO>*> admissibleOuterDOFs;

  std::vector<std::vector<LO>*> nonadmissibleInnerDOFs;

  std::vector<std::vector<LO>*> nonadmissibleOuterDOFs;

  LO maxBlockSize;

  //! vector of admissible matrix blocks
  std::vector<std::pair< FullMatrix<LO, SC>*, FullMatrix<LO, SC>*> >
  admissibleBlocks;

  //! whether to use transformation matrix from p1dis to p1 elements
  bool p12p1dis;

  bool p1dis2p1;

  SparseMatrix<LO, SC> * p12p1disMat;

  SparseMatrix<LO, SC> * p1dis2p1Mat;

  // bool deleteInnerDOFsLists;

  // bool deleteOuterDOFsLists;

private:
  //! copy constructor
  ACAMatrix( const ACAMatrix& orig );

};

}

// include .cpp file to overcome linking problems due to templates
#include "ACAMatrix.cpp"

#endif	/* ACAMATRIX_H */

