/*!
 * @file    BlockLinearOperator.h
 * @author  Jan Zapletal
 * @date    June 22, 2016
 * @brief   Header file for the class BlockLinearOperator
 * 
 */

#ifndef BLOCKLINEAROPERATOR_H
#define	BLOCKLINEAROPERATOR_H

#include "LinearOperator.h"
#include "Vector.h"
#include <vector>

namespace bem4i {

template<class LO, class SC>
class BlockLinearOperator : public LinearOperator<LO, SC> {
public:

  typedef typename GetType<LO, SC>::SCVT SCVT;

  BlockLinearOperator(
      LO blockDimDomain = 0,
      LO blockDimRange = 0
      ) {

    this->blockDimDomain = blockDimDomain;
    this->blockDimRange = blockDimRange;

    LO size = blockDimDomain * blockDimRange;
    this->block.resize( size, nullptr );
    this->trans.resize( size, false );
    this->alpha.resize( size, 1.0 );
  };

  inline void setBlock(
      LO i,
      LO j,
      LinearOperator< LO, SC > * op,
      bool trans = false,
      SC alpha = 1.0
      ) {

    LO ind = i + j * this->blockDimRange;
    this->block[ ind ] = op;
    this->trans[ ind ] = trans;
    this->alpha[ ind ] = alpha;
  }

  inline LinearOperator< LO, SC > * getBlock(
      LO i,
      LO j
      ) const {

    return this->block[ i + j * this->blockDimRange ];
  }

  inline SC getAlpha(
      LO i,
      LO j
      ) const {

    return this->alpha[ i + j * this->blockDimRange ];
  }

  inline bool getTrans(
      LO i,
      LO j
      ) const {

    return this->trans[ i + j * this->blockDimRange ];
  }

  bool isValid( ) {

    LO myDimRange, myDimDomain, test;
    // fixme dirty
    this->dimRange = 0;
    this->dimDomain = 0;

    for ( LO i = 0; i < this->blockDimRange; ++i ) {

      if ( !this->getBlock( i, 0 ) ) return false;

      myDimRange = ( !this->getTrans( i, 0 ) ) ?
          this->getBlock( i, 0 )->getDimRange( ) :
          this->getBlock( i, 0 )->getDimDomain( );
      
      this->dimRange += myDimRange;

      for ( LO j = 1; j < this->blockDimDomain; ++j ) {

        if ( !this->getBlock( i, j ) ) return false;

        test = ( !this->getTrans( i, j ) ) ?
            this->getBlock( i, j )->getDimRange( ) :
            this->getBlock( i, j )->getDimDomain( );

        if ( test != myDimRange ) return false;

      }
    }

    for ( LO j = 0; j < this->blockDimDomain; ++j ) {

      myDimDomain = ( !this->getTrans( 0, j ) ) ?
          this->getBlock( 0, j )->getDimDomain( ) :
          this->getBlock( 0, j )->getDimRange( );
      
      this->dimDomain += myDimDomain;

      for ( LO i = 1; i < this->blockDimRange; ++i ) {

        test = ( !this->getTrans( i, j ) ) ?
            this->getBlock( i, j )->getDimDomain( ) :
            this->getBlock( i, j )->getDimRange( );

        if ( test != myDimDomain ) return false;

      }
    }

    return true;
  }

  //! destructor

  virtual ~BlockLinearOperator( ) {
  };

  /*!
   * @brief Applies operator on a vector
   * 
   * Computes y = beta*y + alpha*this*x
   * @param A 
   * @param x
   * @param y
   * @param alpha
   * @param beta
   */
  virtual void apply(
      Vector< LO, SC > const & x,
      Vector< LO, SC > & y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );

protected:

  BlockLinearOperator( ) {
  };

  LO blockDimDomain;
  LO blockDimRange;
  std::vector< LinearOperator<LO, SC> * > block;
  std::vector< bool > trans;
  std::vector< SC > alpha;

private:

};

}

#include "BlockLinearOperator.cpp"

#endif	/* BLOCKLINEAROPERATOR_H */
