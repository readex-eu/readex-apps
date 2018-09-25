/*!
 * @file    FMMKernel.h
 * @author  Michal Merta 
 * @date    August 21, 2013
 * @brief   Header file for class FMMMatrix
 * 
 */
#ifndef FMMMATRIX_H
#define	FMMMATRIX_H

#include "FullMatrix.h"
#include "FMMKernel.h"

namespace bem4i {

template<class LO, class SC>
class FMMMatrix : public Matrix<LO, SC> {

public:

  //! default constructor
  FMMMatrix( );

  /*!
   * Constructor allocating a full matrix
   * 
   * @param[in]   nRows number of rows
   * @param[in]   nCols number of columns
   * @param[in]   kernel fmm kernel to use
   */
  FMMMatrix( LO nRows, LO nCols, FMMKernel<LO, SC>* kernel );

  //! copy constructor
  FMMMatrix( const FMMMatrix& orig );

  //! destructor
  virtual ~FMMMatrix( );

  //! adds full nonadmissible block to the list of nonadmissible blocks

  inline void addNonadmissibleBlock( FullMatrix<LO, SC>* block ) {
    nonAdmissibleBlocks.push_back( block );
    nonadmBlocksSize += block->getNCols( ) * block->getNRows( );
  }

  //! adds full nonadmissible block associated with idx-th nonadmissible leaf to the list of nonadmissible blocks

  inline void addNonadmissibleBlock( FullMatrix<LO, SC>* block, LO idx ) {
    nonAdmissibleBlocks.at( idx ) = block;
    nonadmBlocksSize += block->getNCols( ) * block->getNRows( );
  }

  //! adds list of nonadmissible cluster pairs

  inline void setNonadmissibleLeaves( std::vector<BEBlockCluster<LO, SC>*> leaves ) {
    nonadmissibleLeaves = leaves;
  }

  //! adds list of nonadmissible cluster pairs

  inline void addAdmissibleLeaves( BEBlockCluster<LO, SC>* leaf ) {
    admissibleLeaves.push_back( leaf );
  }

  //! sets fast multipole kernel

  inline void setKernel( FMMKernel<LO, SC>* kernel ) {
    this->kernel = kernel;
  }

  //! applies matrix to a vector
  virtual void apply(
      Vector<LO, SC> const & x,
      Vector<LO, SC> & y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

  inline SC getCompressionRatio( ) {
    return (SC) ( nonadmBlocksSize / ( ( SC ) this->getNCols( ) * this->getNRows( ) ) );
  }

  inline void resizeNonAdmBlocks( LO size ) {
    this->nonAdmissibleBlocks.resize( size );
  }

  void print( std::ostream &stream = std::cout ) const {
    std::cout << "Fast Multipole Method Matrix\n";
    std::cout << "Number of rows: " << this->nRows << std::endl;
    std::cout << "Number of cols: " << this->nCols << std::endl;
  };

private:

  //! kernel associated with mesh and problem type
  FMMKernel<LO, SC>* kernel;

  //! full size of approximated blocks
  LO nonadmBlocksSize;

  //! vector of nonadmissible leaves
  std::vector<BEBlockCluster<LO, SC>*> nonadmissibleLeaves;

  //! vector of admissible leaves
  std::vector<BEBlockCluster<LO, SC>*> admissibleLeaves;

  //! vector of nonadmissible matrix blocks
  std::vector<FullMatrix<LO, SC>* > nonAdmissibleBlocks;

  //! method applies admissible blocks (approximated by FMM)
  void applyFarfieldBlocks(
      Vector<LO, SC> & result
      );

};

}

// include .cpp file to overcome linking problems due to templates
#include "FMMMatrix.cpp"

#endif	/* FMMMATRIX_H */

