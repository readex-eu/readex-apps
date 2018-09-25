/*!
 * @file    CompoundLinearOperator.h
 * @author  Jan Zapletal
 * @date    April 26, 2016
 * @brief   Header file for the class CompoundLinearOperator
 * 
 */

#ifndef COMPOUNDLINEAROPERATOR_H
#define	COMPOUNDLINEAROPERATOR_H

#include "LinearOperator.h"
#include "Vector.h"
#include <vector>

namespace bem4i {

template<class LO, class SC>
class CompoundLinearOperator : public LinearOperator<LO, SC> {
public:

  typedef typename GetType<LO, SC>::SCVT SCVT;

  CompoundLinearOperator(
      LO reserveSize = 0
      ) {

    this->compound.reserve( reserveSize );
    this->trans.reserve( reserveSize );
    this->alpha.reserve( reserveSize );
    this->maxSize = 0;
  };

  bool isValid( ) const {

    LO size = this->compound.size( );
    LO dimRange, dimDomainNext;

    for ( LO i = 0; i < size - 1; ++i ) {
      if ( !this->compound[ i ] || !this->compound[ i + 1 ] ) return false;

      dimRange = ( !this->trans[ i ] ) ? this->compound[ i ]->getDimRange( ) :
          this->compound[ i ]->getDimDomain( );

      dimDomainNext = ( !this->trans[ i + 1 ] ) ?
          this->compound[ i + 1 ]->getDimDomain( ) :
          this->compound[ i + 1 ]->getDimRange( );
      
      if( dimRange != dimDomainNext ) return false;
    }

    return true;
  }

  //! destructor

  virtual ~CompoundLinearOperator( ) {
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

  void push_back(
      LinearOperator< LO, SC> * op,
      bool trans = false,
      SC alpha = 1.0
      ) {

    if ( op ) {
      this->compound.push_back( op );
      this->trans.push_back( trans );
      this->alpha.push_back( alpha );

      if ( this->compound.size( ) == 1 ) {
        if ( !trans ) {
          this->dimDomain = op->getDimDomain( );
        } else {
          this->dimDomain = op->getDimRange( );
        }
      }

      if ( !trans ) {
        this->dimRange = op->getDimRange( );
      } else {
        this->dimRange = op->getDimDomain( );
      }

      if ( op->getDimRange( ) > this->maxSize ) {
        this->maxSize = op->getDimRange( );
      }

      if ( op->getDimDomain( ) > this->maxSize ) {
        this->maxSize = op->getDimDomain( );
      }
    }
  };

protected:

  std::vector< LinearOperator<LO, SC> * > compound;
  std::vector< bool > trans;
  std::vector< SC > alpha;

  LO maxSize;

private:

};

}

#include "CompoundLinearOperator.cpp"

#endif	/* COMPOUNDLINEAROPERATOR_H */
