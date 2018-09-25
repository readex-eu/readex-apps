/*!
 * @file    CompoundLinearOperator.cpp
 * @author  Jan Zapletal
 * @date    April 26, 2016
 * 
 */


#ifdef COMPOUNDLINEAROPERATOR_H

namespace bem4i {

template< class LO, class SC >
void CompoundLinearOperator< LO, SC >::apply(
    Vector< LO, SC > const & x,
    Vector< LO, SC > & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  LO size = this->compound.size( );

  if ( size == 0 ) {
    return;
  } else if ( size == 1 ) {
    this->compound[ 0 ]->apply( x, y, transA != this->trans[ 0 ],
        alpha * this->alpha[ 0 ], beta );
    return;
  }

  std::vector< SC * > auxData;
  auxData.reserve( 2 );
  auxData.push_back( new SC[ this->maxSize ] );
  auxData.push_back( new SC[ this->maxSize ] );
  Vector< LO, SC > src( this->maxSize, auxData[ 0 ] );
  Vector< LO, SC > tgt( this->maxSize, auxData[ 1 ] );
  src.setAll( 0.0 );
  tgt.setAll( 0.0 );
  LO i;
  LO tgtSize;
  LO srcSize;

  if ( !transA ) {

    tgtSize = ( !this->trans[ 0 ] ) ? this->compound[ 0 ]->getDimRange( ) :
        this->compound[ 0 ]->getDimDomain( );
    tgt.setData( tgtSize, auxData[ 0 ], false );
    this->compound[ 0 ]->apply( x, tgt, this->trans[ 0 ],
        alpha * this->alpha[ 0 ], 0.0 );

    for ( i = 1; i < size - 1; ++i ) {
      srcSize = ( !this->trans[ i ] ) ? this->compound[ i ]->getDimDomain( ) :
          this->compound[ i ]->getDimRange( );
      src.setData( srcSize, auxData[ ( i + 1 ) % 2 ], false );
      tgtSize = ( !this->trans[ i ] ) ? this->compound[ i ]->getDimRange( ) :
          this->compound[ i ]->getDimDomain( );
      tgt.setData( tgtSize, auxData[ i % 2 ], false );
      this->compound[ i ]->apply( src, tgt, this->trans[ i ], this->alpha[ i ],
          0.0 );
    }

    srcSize = ( !this->trans[ i ] ) ? this->compound[ i ]->getDimDomain( ) :
        this->compound[ i ]->getDimRange( );
    src.setData( srcSize, auxData[ ( i + 1 ) % 2 ], false );
    this->compound[ i ]->apply( src, y, this->trans[ i ], this->alpha[ i ],
        beta );
  } else {

    tgtSize = ( !this->trans[ size - 1 ] ) ?
        this->compound[ size - 1 ]->getDimDomain( ) :
        this->compound[ size - 1 ]->getDimRange( );
    tgt.setData( tgtSize, auxData[ 0 ], false );
    this->compound[ size - 1 ]->apply( x, tgt, !this->trans[ size - 1 ],
        alpha * this->alpha[ size - 1 ], 0.0 );
    LO ind;

    for ( i = 1; i < size - 1; ++i ) {
      ind = size - i - 1;
      srcSize = ( !this->trans[ ind ] ) ?
          this->compound[ ind ]->getDimRange( ) :
          this->compound[ ind ]->getDimDomain( );
      src.setData( srcSize, auxData[ ( i + 1 ) % 2 ], false );
      tgtSize = ( !this->trans[ ind ] ) ?
          this->compound[ ind ]->getDimDomain( ) :
          this->compound[ ind ]->getDimRange( );
      tgt.setData( tgtSize, auxData[ i % 2 ], false );
      this->compound[ ind ]->apply( src, tgt, !this->trans[ ind ],
          this->alpha[ ind ], 0.0 );
    }

    srcSize = ( !this->trans[ 0 ] ) ?
        this->compound[ 0 ]->getDimRange( ) :
        this->compound[ 0 ]->getDimDomain( );
    src.setData( srcSize, auxData[ ( i + 1 ) % 2 ], false );
    this->compound[ 0 ]->apply( src, y, !this->trans[ 0 ], this->alpha[ 0 ],
        beta );
  }

  delete [] auxData[ 0 ];
  delete [] auxData[ 1 ];
}

}
#endif
