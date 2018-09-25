/*!
 * @file    SumLinearOperator.h
 * @author  Jan Zapletal
 * @date    June 23, 2016
 * @brief   Header file for the class SumLinearOperator
 * 
 */

#ifndef SUMLINEAROPERATOR_H
#define	SUMLINEAROPERATOR_H

#include "LinearOperator.h"
#include "Vector.h"
#include <vector>

namespace bem4i {

template<class LO, class SC>
class SumLinearOperator : public LinearOperator<LO, SC> {
public:

  typedef typename GetType<LO, SC>::SCVT SCVT;

  SumLinearOperator(
      LO reserveSize = 0
      ) {
    
    this->sum.reserve( reserveSize );
    this->trans.reserve( reserveSize );
    this->alpha.reserve( reserveSize );
  };

  bool isValid( ) const {

    LO size = this->sum.size( );

    if ( size == 0 ) return true;
    if ( !this->sum[ 0 ] ) return false;

    LO dimDomain, dimRange, test;
    if ( !this->trans[ 0 ] ) {
      dimDomain = this->sum[ 0 ]->getDimDomain( );
      dimRange = this->sum[ 0 ]->getDimRange( );
    } else {
      dimDomain = this->sum[ 0 ]->getDimRange( );
      dimRange = this->sum[ 0 ]->getDimDomain( );
    }

    for ( LO i = 1; i < size; ++i ) {
      if ( !this->sum[ i ] ) return false;

      test = ( !this->trans[ i ] ) ? this->sum[ i ]->getDimDomain( ) :
          this->sum[ i ]->getDimRange( );

      if ( test != dimDomain ) return false;

      test = ( !this->trans[ i ] ) ? this->sum[ i ]->getDimRange( ) :
          this->sum[ i ]->getDimDomain( );

      if ( test != dimRange ) return false;

    }

    return true;
  }

  //! destructor

  virtual ~SumLinearOperator( ) {
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
      this->sum.push_back( op );
      this->trans.push_back( trans );
      this->alpha.push_back( alpha );

      if ( this->sum.size( ) == 1 ) {
        if ( !trans ) {
          this->dimDomain = op->getDimDomain( );
          this->dimRange = op->getDimRange( );
        } else {
          this->dimDomain = op->getDimRange( );
          this->dimRange = op->getDimDomain( );
        }
      }
    }
  };

protected:

  std::vector< LinearOperator<LO, SC> * > sum;
  std::vector< bool > trans;
  std::vector< SC > alpha;

private:

};

}

#include "SumLinearOperator.cpp"

#endif	/* SUMLINEAROPERATOR_H */
