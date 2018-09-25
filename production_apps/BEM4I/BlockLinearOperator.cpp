/*!
 * @file    BlockLinearOperator.cpp
 * @author  Jan Zapletal
 * @date    June 22, 2016
 * 
 */


#ifdef BLOCKLINEAROPERATOR_H

namespace bem4i {

template< class LO, class SC >
void BlockLinearOperator< LO, SC >::apply(
    Vector< LO, SC > const & x,
    Vector< LO, SC > & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  Vector<LO, SC> yi;
  Vector<LO, SC> xj;
  LO yPos = 0;
  LO xPos = 0;
  SC * yData = y.getData( );
  SC * xData = x.getData( );
  LO dimRange, dimDomain;
  LinearOperator< LO, SC > * block;

  if ( !transA ) {

    for ( int i = 0; i < this->blockDimRange; ++i ) {
      dimRange = ( !this->getTrans( i, 0 ) ) ?
          this->getBlock( i, 0 )->getDimRange( ) :
          this->getBlock( i, 0 )->getDimDomain( );
      yi.setData( dimRange, yData + yPos, false );
      yi.scale( beta );
      yPos += dimRange;
      xPos = 0;
      for ( int j = 0; j < this->blockDimDomain; ++j ) {
        block = this->getBlock( i, j );
        dimDomain = ( !this->getTrans( i, j ) ) ?
            block->getDimDomain( ) : block->getDimRange( );
        xj.setData( dimDomain, xData + xPos, false );

        block->apply( xj, yi, this->getTrans( i, j ),
            alpha * this->getAlpha( i, j ), 1.0 );
        xPos += dimDomain;
      }
    }
  } else {

    for ( int i = 0; i < this->blockDimDomain; ++i ) {
      dimDomain = ( !this->getTrans( 0, i ) ) ?
          this->getBlock( 0, i )->getDimDomain( ) :
          this->getBlock( 0, i )->getDimRange( );
      yi.setData( dimDomain, yData + yPos, false );
      yi.scale( beta );
      yPos += dimDomain;
      xPos = 0;
      for ( int j = 0; j < this->blockDimRange; ++j ) {
        block = this->getBlock( j, i );
        dimRange = ( !this->getTrans( j, i ) ) ?
            block->getDimRange( ) : block->getDimDomain( );
        xj.setData( dimRange, xData + xPos, false );

        block->apply( xj, yi, !this->getTrans( j, i ),
            alpha * this->getAlpha( j, i ), 1.0 );
        xPos += dimRange;
      }
    }
  }

}

}
#endif
