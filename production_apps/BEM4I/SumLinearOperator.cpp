/*!
 * @file    SumLinearOperator.cpp
 * @author  Jan Zapletal
 * @date    June 23, 2016
 * 
 */


#ifdef SUMLINEAROPERATOR_H

namespace bem4i {

template< class LO, class SC >
void SumLinearOperator< LO, SC >::apply(
    Vector< LO, SC > const & x,
    Vector< LO, SC > & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  LO size = this->sum.size( );

  Vector< LO, SC > yCopy( y.getLength( ), false );
  if ( beta != 0.0 ) {
    y.copy( yCopy );
  }
  y.setAll( 0.0 );

  LO i;
  for ( i = 0; i < size; ++i ) {
    this->sum[ i ]->apply( x, y, transA != this->trans[ i ],
        alpha * this->alpha[ i ], 1.0 );
  }

  if ( beta != 0.0 ) {
    y.add( yCopy, beta );
  }

}

}
#endif
