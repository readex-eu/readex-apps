/*!
 * @file    MathFun.cpp
 * @author  Jan Zapletal
 * @date    February 20, 2015
 * 
 */

#ifdef MATHFUN_H

namespace bem4i {

template<class LO, class SC>
void MathFun<LO, SC>::legendre(
    SCVT x,
    LO maxL,
    std::vector< SCVT > & res
    ) {

  if ( std::abs( x ) > 1.0 || maxL < 0 ) return;

  res.clear( );
  res.resize( (LO) ( 0.5 * ( maxL + 1 ) * ( maxL + 2 ) ) );

  SCVT sqroot = std::sqrt( 1.0 - x * x );
  SCVT p00 = 1;
  SCVT p10 = x;
  SCVT p11 = -sqroot;

  res[ 0 ] = p00;
  if ( maxL < 1 ) return;
  res[ 1 ] = p10;
  res[ 2 ] = p11;

  SCVT p_l_m;
  SCVT p_lm_m;
  for ( LO l = 2; l <= maxL; ++l ) {
    // P_l^0
    p_l_m = res[ lmtoiP( l - 1, 0 ) ];
    p_lm_m = res[ lmtoiP( l - 2, 0 ) ];
    res[ lmtoiP( l, 0 ) ] = legendreIncrementL( l - 1, 0, x, p_l_m, p_lm_m );

    // P_l^m
    //if ( std::abs( x ) > 1.0 - EPS ) {
    if ( sqroot < EPS ) {
      for ( LO m = 1; m <= l; ++m ) {
        res[ lmtoiP( l, m ) ] = 0.0;
      }
    } else {
      for ( LO m = 1; m <= l; ++m ) {

        p_l_m = res[ lmtoiP( l, m - 1 ) ];
        p_lm_m = res[ lmtoiP( l - 1, m - 1 ) ];
        res[ lmtoiP( l, m ) ] =
            legendreIncrementM( l, m - 1, x, p_l_m, p_lm_m, sqroot );
      }
    }
  }

}

template<class LO, class SC>
typename bem4i::MathFun<LO, SC>::SCVT MathFun<LO, SC>::legendreIncrementL(
    LO l,
    LO m,
    SCVT x,
    SCVT p_l_m,
    SCVT p_lm_m
    ) {

  SCVT p_lp_m =
      ( ( 2.0 * l + 1 ) * x * p_l_m - ( l + m ) * p_lm_m ) / ( l - m + 1.0 );

  return p_lp_m;
}

template<class LO, class SC>
typename bem4i::MathFun<LO, SC>::SCVT MathFun<LO, SC>::legendreIncrementM(
    LO l,
    LO m,
    SCVT x,
    SCVT p_l_m,
    SCVT p_lm_m,
    SCVT sqroot
    ) {

  SCVT p_l_mp =
      ( ( l - m ) * x * p_l_m - ( l + m ) * p_lm_m ) / sqroot;

  return p_l_mp;
}

template<class LO, class SC>
void MathFun<LO, SC>::sphericalHarmonicsRealBasis(
    const SCVT * x,
    LO maxL,
    std::vector< SCVT > & res
    ) {

  if ( maxL < 0 ) return;

  SCVT r = NORM3( x );
  if ( std::abs( r - 1.0 ) >= EPS ) {
    std::cout << "Point not lying on the unit sphere!" << std::endl;
    return;
  }

  res.clear( );
  res.resize( ( maxL + 1 ) * ( maxL + 1 ) );

  // convert to spherical coordinates
  //SCVT theta = std::acos( x[ 2 ] );
  SCVT phi = std::atan2( x[ 1 ], x[ 0 ] );

  std::vector< SCVT > leg( 0.5 * ( maxL + 1 ) * ( maxL + 2 ) );
  legendre( x[ 2 ], maxL, leg );

  // precompute factorials 0! ... 2L!
  LO * fact = new LO[ 2 * maxL + 1 ];
  fact[ 0 ] = 1;
  for ( LO l = 1; l <= 2 * maxL; ++l ) {
    fact[ l ] = fact[ l - 1 ] * l;
  }

  SCVT factor;
  //SCVT sqrt2 = std::sqrt( 2.0 );
  SCVT sqrt2 = 1.0; // not weighing
  SCVT pifact = 1.0 / ( 4.0 * M_PI );
  for ( LO l = 0; l <= maxL; ++l ) {
    // Y_l^0
    res[ lmtoiY( l, 0 ) ] =
        std::sqrt( pifact * ( 2.0 * l + 1.0 ) ) * leg[ lmtoiP( l, 0 ) ];

    for ( LO m = 1; m <= l; ++m ) {
      factor = sqrt2 * std::sqrt( pifact * ( 2.0 * l + 1.0 ) *
          ( (SCVT) fact[ l - m ] / fact[ l + m ] ) ) * leg[ lmtoiP( l, m ) ];
      res[ lmtoiY( l, m ) ] = factor * std::cos( m * phi );
      res[ lmtoiY( l, -m ) ] = factor * std::sin( m * phi );
    }
  }

  delete [] fact;

}

template<class LO, class SC>
void MathFun<LO, SC>::sphericalHarmonicsRealBasis(
    const SCVT * x,
    LO nPoints,
    LO maxL,
    std::vector< Vector< LO, SCVT > * > & res
    ) {

  if ( res.size( ) != ( maxL + 1 ) * ( maxL + 1 ) ) return;
  for ( auto it = res.begin( ); it != res.end( ); ++it ) {
    if ( *it == nullptr ) return;
    ( *it )->resize( nPoints, false );
  }

  std::vector< std::vector< SCVT > > aux( nPoints );
  for ( LO i = 0; i < nPoints; ++i ) {
    sphericalHarmonicsRealBasis( x + 3 * i, maxL, aux[ i ] );
  }

  for ( LO i = 0; i < nPoints; ++i ) {
    for ( LO l = 0; l <= maxL; ++l ) {
      for ( LO m = -l; m <= l; ++m ) {
        res[ lmtoiY( l, m ) ]->set( i, aux[ i ][ lmtoiY( l, m ) ] );
      }
    }
  }

}

}

#endif
