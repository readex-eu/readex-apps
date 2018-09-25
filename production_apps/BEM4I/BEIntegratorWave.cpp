/*!
 * @file    BEIntegratorWave.cpp
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    December 12, 2013
 * 
 */

#ifdef BEINTEGRATORWAVE_H

#ifndef EXPERIMENTAL_WAVE

// use ordinary basis functions described in Sauter, Veit: Retarded time-domain...

namespace bem4i {

template<class LO, class SC>
BEIntegratorWave<LO, SC>::BEIntegratorWave( ) {

}

template<class LO, class SC>
BEIntegratorWave<LO, SC>::BEIntegratorWave( const BEIntegratorWave& orig ) {
}

template<class LO, class SC>
BEIntegratorWave<LO, SC>::BEIntegratorWave(
    BESpace<LO, SC>* space,
    int* quadratureOrder,
    int timeQuadOrder,
    int* quadratureOrderDisjointElems,
    int nChebIntervals,
    int nChebPoints
    ) {

  this->space = space;
  this->quadratureOrder = quadratureOrder;
  this->timeQuadOrder = timeQuadOrder;
  this->quadrature = SauterSchwab;
  this->quadratureOrderDisjointElems = quadratureOrderDisjointElems;
  
  this->currLegOrderBasis = 0;
  this->currLegOrderTest = 0;

  N = static_cast<BESpaceTime<LO, SC>*> ( this->space )->getNTimeSteps( );
  dt = static_cast<BESpaceTime<LO, SC>*> ( this->space )->getDt( );

  // initialize arrays for Chebyshev interpolation of temporal basis functions
  this->nChebIntervals = nChebIntervals;
  this->nChebPoints = nChebPoints;
  chebIntervalStarts = new SC[nChebIntervals + 1];
  psiChebCoeffs = new SCVT*[nChebIntervals];
  psiTildeChebCoeffs = new SCVT*[nChebIntervals];
  psi1LayerChebCoeffs = new SCVT*[nChebIntervals];
  for ( int i = 0; i < nChebIntervals; i++ ) {
    psiChebCoeffs[i] = new SCVT[nChebPoints];
    psiTildeChebCoeffs[i] = new SCVT[nChebPoints];
    psi1LayerChebCoeffs[i] = new SCVT[nChebPoints];
  }

  if ( quadratureOrderDisjointElems ) {
    this->initDisjointQuadratureData( quadratureOrderDisjointElems );
  }
  
  this->initSauterSchwabQuadratureData( quadratureOrder );

}

template<class LO, class SC>
BEIntegratorWave<LO, SC>::~BEIntegratorWave( ) {
  delete [] chebIntervalStarts;
  for ( int i = 0; i < nChebIntervals; i++ ) {
    delete [] psiChebCoeffs[i];
    delete [] psiTildeChebCoeffs[i];
    delete [] psi1LayerChebCoeffs[i];
  }
  delete [] psiChebCoeffs;
  delete [] psiTildeChebCoeffs;
  delete [] psi1LayerChebCoeffs;



}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::computeElemMatrix1Layer(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const {

  if ( this->space->getTestFunctionType( ) == p0 &&
      this->space->getAnsatzFunctionType( ) == p0 ) {
    switch ( this->quadrature ) {
      case SauterSchwab:
        this->computeElemMatrix1LayerSauterSchwabP0P0( outerElem, innerElem,
            matrix );
        break;
      case Steinbach:
        std::cout << "Not implemented!" << std::endl;
        break;
    }
  } else if ( this->space->getTestFunctionType( ) == p1 &&
      this->space->getAnsatzFunctionType( ) == p1 ) {
    //computeElemMatrix1LayerP1P1( outerElem, innerElem, matrix );
    std::cout << "Not implemented!" << std::endl;
  } else {
    std::cout << "Not implemented!" << std::endl;
  }
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::computeElemMatrix2Layer(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const {

  if ( this->space->getTestFunctionType( ) == p0 &&
      this->space->getAnsatzFunctionType( ) == p1 ) {
    if ( this->quadratureOrderDisjointElems &&
        this->areElementsDisjoint( outerElem, innerElem ) ) {
      this->computeElemMatrix1LayerDisjointP0P0( outerElem, innerElem, matrix );
    }
    switch ( this->quadrature ) {
      case SauterSchwab:
        this->computeElemMatrix2LayerSauterSchwabP0P1( outerElem, innerElem,
            matrix );
        break;
      case Steinbach:
        std::cout << "Not implemented!" << std::endl;
        break;
    }
  } else if ( this->space->getTestFunctionType( ) == p0 &&
      this->space->getAnsatzFunctionType( ) == p0 ) {
    //computeElemMatrix2LayerP0P0( outerElem, innerElem, matrix );
    std::cout << "Not implemented!" << std::endl;
  } else if ( this->space->getTestFunctionType( ) == p1 &&
      this->space->getAnsatzFunctionType( ) == p1 ) {
    //computeElemMatrix2LayerP1P1( outerElem, innerElem, matrix );
    std::cout << "Not implemented!" << std::endl;
  } else {
    std::cout << "Not implemented!" << std::endl;
  }
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::computeElemMatrixHypersingular(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const {

  if ( this->space->getTestFunctionType( ) == p1 &&
      this->space->getAnsatzFunctionType( ) == p1 ) {
    if ( this->quadratureOrderDisjointElems &&
        this->areElementsDisjoint( outerElem, innerElem ) ) {
      this->computeElemMatrixHypersingularDisjointP1P1( outerElem, innerElem,
          matrix );
      //std::cout << " bb " << std::endl;
      return;
    } else {
      computeElemMatrixHypersingularSauterSchwabP1P1( outerElem, innerElem,
          matrix );
    }
  } else {
    std::cout << "Not implemented!" << std::endl;
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::evalErf(
    SCVT t,
    SCVT a,
    SCVT b
    ) const {
  SCVT arg = 2.0 * ( t - a ) / ( b - a ) - 1.0;
  if ( arg <= -1.0 ) {
    return 0.0;
  } else if ( arg >= 1.0 ) {
    return 1.0;
  } else {
    return 0.5 * ( erf( 2.0 * atanh( arg ) ) + 1.0 );
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::
evalErfDot(
    SCVT t,
    SCVT a,
    SCVT b
    ) const {

  SCVT arg = ( a + b - 2.0 * t ) / ( a - b );
  if ( arg <= -1.0 ) {
    return 0.0;
  } else if ( arg >= 1.0 ) {
    return 0.0;
  } else {
    SCVT ath = atanh( arg );
    return ( ( ( a - b ) * std::exp( -4.0 * ath * ath ) ) /
        ( std::sqrt( M_PI )*( a - t )*( b - t ) ) );
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::
evalErfDDot(
    SCVT t,
    SCVT a,
    SCVT b
    ) const {
  SCVT arg = ( a + b - 2.0 * t ) / ( a - b );
  if ( arg <= -1.0 ) {
    return 0.0;
  } else if ( arg >= 1.0 ) {
    return 0.0;
  } else {
    SCVT ath = atanh( arg );
    return ( ( ( a - b ) * std::exp( -4.0 * ath * ath ) *
        ( -4.0 * ( a - b ) * ath + a + b - 2.0 * t ) ) /
        ( std::sqrt( M_PI )*( a - t )*( a - t )*( b - t )*( b - t ) ) );
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::
evalLegendrePolynomial(
    SCVT t,
    int order
    ) const {
  if ( t<-1.0 || t > 1.0 ) {
    return 0.0;
  }
  switch ( order ) {
    case 0:
      return 1.0;
    case 1:
      return t;
    case 2:
      return ( 1.5 * t * t - 0.5 );
    case 3:
      return (0.5 * ( 5.0 * t * t * t - 3.0 * t ) );
    case 4:
      return (0.125 * ( 35.0 * t * t * t * t - 30.0 * t * t + 3.0 ) );
  }
  return 0.0;
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::
evalLegendrePolynomialDot(
    SCVT t,
    int order
    ) const {
  if ( t<-1.0 || t > 1.0 ) {
    return 0.0;
  }
  switch ( order ) {
    case 0:
      return 0.0;
    case 1:
      return 1.0;
    case 2:
      return ( 3.0 * t );
    case 3:
      return ( 7.5 * t * t - 1.5 );
    case 4:
      return (17.5 * t * t * t - 7.5 * t );
  }
  return 0.0;
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::
evalLegendrePolynomialDDot(
    SCVT t,
    int order
    ) const {
  if ( t<-1.0 || t > 1.0 ) {
    return 0.0;
  }
  switch ( order ) {
    case 0:
      return 0.0;
    case 1:
      return 0.0;
    case 2:
      return ( 3.0 );
    case 3:
      return ( 15.0 * t );
    case 4:
      return ( 52.5 * t * t - 7.5 );
  }
  return 0.0;
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::
evalPUM(
    SCVT t,
    LO i
    ) const {

  SCVT tIMinus1 = dt * ( i - 1 );
  SCVT tI = dt * i;
  SCVT tIPlus1 = dt * ( i + 1 );

  if ( i == 0 ) {
    return 1.0 - evalErf( t, tI, tIPlus1 );
  } else if ( i == N - 1 ) {
    return evalErf( t, tIMinus1, tI );
  } else {
    if ( t <= tI ) {
      return evalErf( t, tIMinus1, tI );
    } else {
      return 1.0 - evalErf( t, tI, tIPlus1 );
    }
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::
evalPUMDot(
    SCVT t,
    LO i
    ) const {

  SCVT tIMinus1 = dt * ( i - 1 );
  SCVT tI = dt * i;
  SCVT tIPlus1 = dt * ( i + 1 );

  if ( i == 0 ) {
    return (-evalErfDot( t, tI, tIPlus1 ) );
  } else if ( i == N - 1 ) {
    return evalErfDot( t, tIMinus1, tI );
  } else {
    if ( fabs( t - tI ) < EPS ) {
      return 0.0;
    } else if ( t < tI ) {
      return evalErfDot( t, tIMinus1, tI );
    } else {
      return (-evalErfDot( t, tI, tIPlus1 ) );
    }
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::
evalPUMDDot(
    SCVT t,
    LO i
    ) const {

  SCVT tIMinus1 = dt * ( i - 1 );
  SCVT tI = dt * i;
  SCVT tIPlus1 = dt * ( i + 1 );

  if ( i == 0 ) {
    return ( -evalErfDDot( t, tI, tIPlus1 ) );
  } else if ( i == N - 1 ) {
    return evalErfDDot( t, tIMinus1, tI );
  } else {
    if ( fabs( t - tI ) < EPS ) {
      return 0.0;
    } else if ( t < tI ) {
      return evalErfDDot( t, tIMinus1, tI );
    } else {
      return (-evalErfDDot( t, tI, tIPlus1 ) );
    }
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::
evalB(
    SCVT t,
    LO i,
    LO order
    ) const {

  SCVT tIMinus1 = dt * ( i - 1 );
  SCVT tI = dt * i;
  SCVT tIPlus1 = dt * ( i + 1 );

  if ( i == 0 ) {
    //order = ( order < 2 ) ? 2 : order;
    return evalPUM( t, i ) * t * t *
        evalLegendrePolynomial( 2.0 * t / dt - 1.0, order ) *8.0 / ( dt * dt );
  } else if ( i == N - 1 ) {
    return evalPUM( t, i ) *
        evalLegendrePolynomial( 2.0 * ( t - tIMinus1 ) / ( dt ) - 1.0, order );
  } else {
    return evalPUM( t, i ) *
        evalLegendrePolynomial( ( t - tIMinus1 ) / ( dt ) - 1.0, order );
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::
evalBDot(
    SCVT t,
    LO i,
    LO order
    ) const {

  SCVT tIMinus1 = dt * ( i - 1 );
  SCVT tI = dt * i;
  SCVT tIPlus1 = dt * ( i + 1 );

  if ( i == 0 ) {
    //return ( t * t * evalPUMDot( t, i ) + 2.0 * t * evalPUM( t, i ) ) *
    //8.0 / ( dt * dt );
    //order = ( order < 2 ) ? 2 : order;
    return ( t * ( 2.0 * evalPUM( t, i ) *
        ( ( t * evalLegendrePolynomialDot( 2.0 * t / tIPlus1 - 1.0, order ) ) /
        ( tIPlus1 ) + evalLegendrePolynomial( 2.0 * t / tIPlus1 - 1.0, order ) )
        + t * evalLegendrePolynomial( 2.0 * t / tIPlus1 - 1.0, order ) *
        evalPUMDot( t, i ) ) ) *8.0 / ( dt * dt );
  } else if ( i == N - 1 ) {
    //return evalPUMDot( t, i );
    return (2.0 * evalPUM( t, i ) *
        evalLegendrePolynomialDot( ( tI + tIMinus1 - 2.0 * t ) /
        ( tIMinus1 - tI ), order ) / ( tI - tIMinus1 ) +
        evalPUMDot( t, i ) * evalLegendrePolynomial( ( tI + tIMinus1 - 2.0 * t )
        / ( tIMinus1 - tI ), order ) );
  } else {
    //return evalPUMDot( t, i );
    return (2.0 * evalPUM( t, i ) *
        evalLegendrePolynomialDot( ( tIPlus1 + tIMinus1 - 2.0 * t ) /
        ( tIMinus1 - tIPlus1 ), order ) / ( tIPlus1 - tIMinus1 ) +
        evalPUMDot( t, i ) *
        evalLegendrePolynomial( ( tIPlus1 + tIMinus1 - 2.0 * t ) /
        ( tIMinus1 - tIPlus1 ), order ) );
  }

}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::
evalBDDot(
    SCVT t,
    LO i,
    LO order
    ) const {

  SCVT tIMinus1 = dt * ( i - 1 );
  SCVT tI = dt * i;
  SCVT tIPlus1 = dt * ( i + 1 );

  if ( i == 0 ) {
    //return (t * t * evalPUMDDot( t, i ) + 
    //4.0 * t * evalPUMDot( t, i ) + 2.0 * evalPUM( t, i ) ) *8 / ( dt * dt );
    //order = ( order < 2 ) ? 2 : order;
    return ( ( 4.0 * t * t * evalPUM( t, i ) *
        evalLegendrePolynomialDDot( 2 * t / tIPlus1 - 1.0, order ) ) /
        ( tIPlus1 * tIPlus1 ) +
        ( 4.0 * t * evalLegendrePolynomialDot( 2.0 * t / tIPlus1 - 1.0, order )
        * ( 2.0 * evalPUM( t, i ) + t * evalPUMDot( t, i ) ) ) / ( tIPlus1 ) +
        evalLegendrePolynomial( 2.0 * t / tIPlus1 - 1.0, order )*
        ( 2.0 * evalPUM( t, i ) + t * ( 4.0 * evalPUMDot( t, i ) + t *
        evalPUMDDot( t, i ) ) ) )*8 / ( dt * dt );
  } else if ( i == N - 1 ) {
    //return evalPUMDDot( t, i );
    return ( ( 4.0 * evalPUM( t, i ) *
        evalLegendrePolynomialDDot( ( tI + tIMinus1 - 2.0 * t ) /
        ( tIMinus1 - tI ), order ) ) / ( ( tI - tIMinus1 )*( tI - tIMinus1 ) ) +
        ( 4.0 * evalPUMDot( t, i ) *
        evalLegendrePolynomialDot( ( tI + tIMinus1 - 2.0 * t ) /
        ( tIMinus1 - tI ), order ) ) / ( tI - tIMinus1 ) + evalPUMDDot( t, i ) *
        evalLegendrePolynomial( ( tI + tIMinus1 - 2.0 * t ) /
        ( tIMinus1 - tI ), order ) );
  } else {
    //return evalPUMDDot( t, i );
    return ( ( 4.0 * evalPUM( t, i ) *
        evalLegendrePolynomialDDot( ( tIPlus1 + tIMinus1 - 2.0 * t ) /
        ( tIMinus1 - tIPlus1 ), order ) ) / ( ( tIPlus1 - tIMinus1 )*
        ( tIPlus1 - tIMinus1 ) ) + ( 4.0 * evalPUMDot( t, i ) *
        evalLegendrePolynomialDot( ( tIPlus1 + tIMinus1 - 2.0 * t ) /
        ( tIMinus1 - tIPlus1 ), order ) ) / ( tIPlus1 - tIMinus1 ) +
        evalPUMDDot( t, i ) * evalLegendrePolynomial(
        ( tIPlus1 + tIMinus1 - 2.0 * t ) / ( tIMinus1 - tIPlus1 ), order ) );
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::
integrate1Layer(
    SCVT start,
    SCVT end,
    SCVT r
    ) const {
  SCVT ret = 0.0;
  LO qSize = lineQuadSizes[timeQuadOrder];
  SCVT *qPoints = lineQuadPoints[timeQuadOrder];
  SCVT *qWeights = lineQuadWeights[timeQuadOrder];

  SCVT length = end - start;
  SCVT t = 0.0;

  for ( int i = 0; i < qSize; i++ ) {
    t = start + length * qPoints[i];
    ret += qWeights[i] * evalBDot( t - r, currBasisF, currLegOrderBasis ) *
        evalB( t, currTestF, currLegOrderTest );
  }
  ret *= length;
  return ret;
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::
integrateHypersingularPart2(
    SCVT start,
    SCVT end,
    SCVT r
    ) const {
  SCVT ret = 0.0;
  LO qSize = lineQuadSizes[timeQuadOrder];
  SCVT *qPoints = lineQuadPoints[timeQuadOrder];
  SCVT *qWeights = lineQuadWeights[timeQuadOrder];

  SCVT length = end - start;
  SCVT t = 0.0;

  for ( int i = 0; i < qSize; i++ ) {
    t = start + length * qPoints[i];
    ret += qWeights[i] * evalBDDot( t - r, currBasisF, currLegOrderBasis ) *
        evalBDot( t, currTestF, currLegOrderTest );
  }
  ret *= length;
  return ret;
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::
integrateHypersingularPart1(
    SCVT start,
    SCVT end,
    SCVT r
    ) const {
  // evaluates psi tilde
  SCVT ret = 0.0;
  LO qSize = lineQuadSizes[timeQuadOrder];
  SCVT *qPoints = lineQuadPoints[timeQuadOrder];
  SCVT *qWeights = lineQuadWeights[timeQuadOrder];

  SCVT length = end - start;
  SCVT t = 0.0;

  for ( int i = 0; i < qSize; i++ ) {
    t = start + length * qPoints[i];
    ret += qWeights[i] * evalB( t - r, currBasisF, currLegOrderBasis ) *
        evalBDot( t, currTestF, currLegOrderTest );
  }
  ret *= length;
  return ret;
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::getDirichletRHS(
    Vector<LO, SC> &rhs
    ) const {

  BESpaceTime<LO, SC> *spaceTime =
      static_cast<BESpaceTime<LO, SC>*> ( this->space );

  // max order of Legendre polynomials in temporal basis
  int legOrder = spaceTime->getLegendreOrder( );

  LO nElems = this->space->getMesh( )->getNElements( );
  rhs.resize( nElems * N * ( legOrder + 1 ) );

  SCVT ret = 0.0;
  LO qSize = lineQuadSizes[timeQuadOrder];
  SCVT *qPoints = lineQuadPoints[timeQuadOrder];
  SCVT *qWeights = lineQuadWeights[timeQuadOrder];

  SCVT t = 0.0;
  SC g;
  SCVT length, start, end;

  for ( int i = 0; i < N; i++ ) {

    for ( int currentLeg = 0; currentLeg < legOrder + 1; currentLeg++ ) {

      if ( i == 0 ) {
        start = 0.0;
        end = dt;
      } else if ( i == N - 1 ) {
        start = ( N - 2 ) * dt;
        end = ( N - 1 ) * dt;
      } else {
        start = ( i - 1 ) * dt;
        end = ( i + 1 ) * dt;
      }
      length = end - start;
      ret = 0.0;
      for ( int q = 0; q < qSize; q++ ) {
        t = start + length * qPoints[q];
        g = ( 4.0 * t * t * t - 2.0 * t * t * t * t ) * std::exp( -2.0 * t );
        ret += qWeights[q] * evalB( t, i, currentLeg ) * g;
      }
      ret *= length;
      for ( LO j = 0; j < nElems; j++ ) {
        rhs.set( ( i * ( legOrder + 1 ) + currentLeg ) * nElems + j,
            this->space->getMesh( )->getElemArea( j ) * ret );
      }
    }
  }
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::computeChebyshevForHypersingular(
    ) {

  // at first, compute supports of psidot and psiddot depending on current
  // test and basis functions
  SCVT minTest, maxTest, minBasis, maxBasis;
  SCVT minSupp, maxSupp;
  if ( currTestF == 0 ) {
    minTest = 0.0;
    maxTest = dt;
  } else if ( currTestF == N - 1 ) {
    minTest = ( N - 2 ) * dt;
    maxTest = ( N - 1 ) * dt;
  } else {
    minTest = ( currTestF - 1 ) * dt;
    maxTest = ( currTestF + 1 ) * dt;
  }

  if ( currBasisF == 0 ) {
    minBasis = 0.0;
    maxBasis = dt;
  } else if ( currBasisF == N - 1 ) {
    minBasis = ( N - 2 ) * dt;
    maxBasis = ( N - 1 ) * dt;
  } else {
    minBasis = ( currBasisF - 1 ) * dt;
    maxBasis = ( currBasisF + 1 ) * dt;
  }

  minSupp = std::max( minTest - maxBasis, 0.0 );
  maxSupp = maxTest - minBasis;

  // roots of Chebyshev polynomial
  SC *zeros = new SC[nChebPoints];
  SC *psi = new SC[nChebPoints];
  SC *psiTilde = new SC[nChebPoints];
  SC t;
  chebDelta = ( maxSupp - minSupp ) / ( (SC) nChebIntervals );
  for ( int i = 0; i <= nChebIntervals; i++ ) {
    chebIntervalStarts[i] = minSupp + i * chebDelta;
  }
  for ( int i = 0; i < nChebIntervals; i++ ) {
    // compute roots of Chebyshev pol
    chebyshevZeros( chebIntervalStarts[i], chebIntervalStarts[i + 1], zeros );

    // compute values of psi and psitilde in zeros of Cheb. interpolant
    for ( int j = 0; j < nChebPoints; j++ ) {
      t = zeros[j];
      SCVT intStart = std::min( std::max( minTest, minBasis + t ), ( N - 1 ) *
          dt );
      SCVT intEnd = std::min( maxBasis + t, maxTest );

      if ( intEnd <= intStart ) {
        psi[j] = 0.0;
        psiTilde[j] = 0.0;
      } else {
        psi[j] = integrateHypersingularPart2( intStart, intEnd, t );
        psiTilde[j] = integrateHypersingularPart1( intStart, intEnd, t );
      }

    }
    chebyshevCoefficients( psi, psiChebCoeffs[i] );
    chebyshevCoefficients( psiTilde, psiTildeChebCoeffs[i] );
  }

  delete [] zeros;
  delete [] psi;
  delete [] psiTilde;
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::computeChebyshevFor1Layer(
    ) {

  // at first, compute supports of psidot and psiddot depending on
  // current test and basis functions
  SCVT minTest, maxTest, minBasis, maxBasis;
  SCVT minSupp, maxSupp;
  if ( currTestF == 0 ) {
    minTest = 0.0;
    maxTest = dt;
  } else if ( currTestF == N - 1 ) {
    minTest = ( N - 2 ) * dt;
    maxTest = ( N - 1 ) * dt;
  } else {
    minTest = ( currTestF - 1 ) * dt;
    maxTest = ( currTestF + 1 ) * dt;
  }

  if ( currBasisF == 0 ) {
    minBasis = 0.0;
    maxBasis = dt;
  } else if ( currBasisF == N - 1 ) {
    minBasis = ( N - 2 ) * dt;
    maxBasis = ( N - 1 ) * dt;
  } else {
    minBasis = ( currBasisF - 1 ) * dt;
    maxBasis = ( currBasisF + 1 ) * dt;
  }

  minSupp = std::max( minTest - maxBasis, 0.0 );
  maxSupp = maxTest - minBasis;

  // roots of Chebyshev polynomial
  SC *zeros = new SC[nChebPoints];
  SC *psi = new SC[nChebPoints];
  SC t;
  for ( int i = 0; i <= nChebIntervals; i++ ) {
    chebIntervalStarts[i] = minSupp + i * ( maxSupp - minSupp ) /
        ( (SC) nChebIntervals );
  }
  for ( int i = 0; i < nChebIntervals; i++ ) {
    // compute roots of Chebyshev pol
    chebyshevZeros( chebIntervalStarts[i], chebIntervalStarts[i + 1], zeros );

    // compute values of psi and psitilde in zeros of Cheb. interpolant
    for ( int j = 0; j < nChebPoints; j++ ) {
      t = zeros[j];
      SCVT intStart = std::min( std::max( minTest, minBasis + t ), ( N - 1 ) *
          dt );
      SCVT intEnd = std::min( maxBasis + t, maxTest );

      if ( intEnd <= intStart ) {
        psi[j] = 0.0;
      } else {
        psi[j] = integrate1Layer( intStart, intEnd, t );
      }

    }
    chebyshevCoefficients( psi, psi1LayerChebCoeffs[i] );
  }

  delete [] zeros;
  delete [] psi;
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::chebyshevZeros(
    SC a,
    SC b,
    SC *zeros
    ) const {
  SC angle;
  int i;
  for ( i = 0; i < nChebPoints; i++ ) {
    angle = (SC) ( 2 * i + 1 ) * M_PI / (SC) ( 2 * nChebPoints );
    zeros[i] = std::cos( angle );
    zeros[i] = 0.5 * ( a + b ) + zeros[i]*0.5 * ( b - a );
  }
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::chebyshevCoefficients(
    SC *fx,
    SC *coeffs
    ) const {
  SC angle;
  memset( coeffs, 0, nChebPoints * sizeof (SC ) );
  for ( int i = 0; i < nChebPoints; i++ ) {
    for ( int j = 0; j < nChebPoints; j++ ) {
      angle = (SC) ( i * ( 2 * j + 1 ) ) * M_PI / (SC) ( 2 * nChebPoints );
      coeffs[i] += fx[j] * std::cos( angle );
    }
    coeffs[i] *= 2.0 / (SC) ( nChebPoints );
  }
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::getNeumannRHS(
    Vector<LO, SC> &rhs
    ) const {

  BESpaceTime<LO, SC> *spaceTime =
      static_cast<BESpaceTime<LO, SC>*> ( this->space );

  // max order of Legendre polynomials in temporal basis
  int legOrder = spaceTime->getLegendreOrder( );

  LO nElems = this->space->getMesh( )->getNElements( );
  LO nNodes = this->space->getMesh( )->getNNodes( );
  rhs.resize( nNodes * N * ( legOrder + 1 ) );
  rhs.setAll( 0.0 );

  SCVT ret = 0.0;
  LO qSize = lineQuadSizes[timeQuadOrder];
  SCVT *qPoints = lineQuadPoints[timeQuadOrder];
  SCVT *qWeights = lineQuadWeights[timeQuadOrder];

  // variables for Gaussian quadrature 
  int quadOrder = 8;
  int quad_size = quadSizes[quadOrder];
  double *quad_points = quadPoints[quadOrder];
  double *quad_weight = quadWeights[quadOrder];

  SCVT t = 0.0;
  SC gVal;
  SCVT length, start, end;
  //SCVT x1[3], x2[3], x3[3];
  //SCVT* qNodes = new SCVT[ quadSizes[ this->quadratureOrder[ 0 ] ] * 3 ];
  SCVT phi1, phi2, phi3, intPx, intPy;
  vector<LO> colIndices( 3 );
  SCVT x1[3], x2[3], x3[3], n[3];
  SCVT * qNodes = new SCVT[ 3 * quad_size ];
  SCVT * xLocal;

  for ( int i = 0; i < N; i++ ) {

    for ( int currentLeg = 0; currentLeg < legOrder + 1; currentLeg++ ) {

      if ( i == 0 ) {
        start = 0.0;
        end = dt;
      } else if ( i == N - 1 ) {
        start = ( N - 2 ) * dt;
        end = ( N - 1 ) * dt;
      } else {
        start = ( i - 1 ) * dt;
        end = ( i + 1 ) * dt;
      }
      length = end - start;
      ret = 0.0;


      for ( int j = 0; j < nElems; j++ ) {

        this->space->getMesh( )->getNodes( j, x1, x2, x3 );
        this->getQuadratureNodes( x1, x2, x3, quadOrder, qNodes );
        this->space->getOuterElemDOFs( j, &colIndices[0] );
        this->space->getMesh( )->getNormal( j, n );

        SCVT cent[3];
        this->space->getMesh( )->getCentroid( j, cent );


        for ( int k = 0; k < quad_size; k++ ) {
          intPx = quad_points[ 2 * k ];
          intPy = quad_points[ 2 * k + 1 ];
          phi1 = 1.0 - intPx - intPy;
          phi2 = intPx;
          phi3 = intPy;
          xLocal = qNodes + 3 * k;
          //phi1 = 1.0 - intPx;
          //phi2 = intPx - intPy;
          //phi3 = intPy;
          ret = 0.0;

          for ( int q = 0; q < qSize; q++ ) {

            t = start + length * qPoints[q];
            gVal = uIncNeu( t, xLocal, n, 0 );
            //g = ( t * t * t * t ) * std::exp( -2.0 * t );

            //g = std::sin( 4.0 * t * t * t * t ) * t * std::exp( -t );
            //g = std::sin( 3.0 * t ) * t * t * std::exp( -t );
            //g = ( std::sin( 2.0 * t ) )*( std::sin( 2.0 * t ) ) * t * 
            //  std::exp( -t ) * 0.5 * sqrt( 3.0 / M_PI ) * cent[2];

            if ( t < 0.0 ) {
              gVal = 0.0;
            }

            //}
            SCVT scaleFact = 1.0;
            if ( legOrder == 0 ) {
              scaleFact = ( std::sqrt( M_PI ) );
            }
            ret += qWeights[q] * gVal * evalBDot( t, i, currentLeg ) *
                scaleFact;

          }

          ret *= length;
          //ret = 1.0;

          rhs.add( ( i * ( legOrder + 1 ) + currentLeg ) * nNodes +
              colIndices[0], quad_weight[k] * phi1 * ret *
              this->space->getMesh( )->getElemArea( j ) );
          rhs.add( ( i * ( legOrder + 1 ) + currentLeg ) * nNodes +
              colIndices[1], quad_weight[k] * phi2 * ret *
              this->space->getMesh( )->getElemArea( j ) );
          rhs.add( ( i * ( legOrder + 1 ) + currentLeg ) * nNodes +
              colIndices[2], quad_weight[k] * phi3 * ret *
              this->space->getMesh( )->getElemArea( j ) );
        }
      }
    }
  }
  //rhs.scale( 1.0 / ( sqrt( M_PI )*2 ) );
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::getNeumannRHS(
    Vector<LO, SC> &rhs,
    SC( *incWave ) ( SCVT, SCVT*, SCVT* )
    ) const {

  BESpaceTime<LO, SC> *spaceTime =
      static_cast<BESpaceTime<LO, SC>*> ( this->space );

  // max order of Legendre polynomials in temporal basis
  int legOrder = spaceTime->getLegendreOrder( );

  LO nElems = this->space->getMesh( )->getNElements( );
  LO nNodes = this->space->getMesh( )->getNNodes( );
  rhs.resize( nNodes * N * ( legOrder + 1 ) );
  rhs.setAll( 0.0 );

  SCVT ret = 0.0;
  LO qSize = lineQuadSizes[timeQuadOrder];
  SCVT *qPoints = lineQuadPoints[timeQuadOrder];
  SCVT *qWeights = lineQuadWeights[timeQuadOrder];

  // variables for Gaussian quadrature 
  int quadOrder = 8;
  int quad_size = quadSizes[quadOrder];
  double *quad_points = quadPoints[quadOrder];
  double *quad_weight = quadWeights[quadOrder];

  SCVT t = 0.0;
  SC gVal;
  SCVT length, start, end;
  //SCVT x1[3], x2[3], x3[3];
  //SCVT* qNodes = new SCVT[ quadSizes[ this->quadratureOrder[ 0 ] ] * 3 ];
  SCVT phi1, phi2, phi3, intPx, intPy;
  vector<LO> colIndices( 3 );
  SCVT x1[3], x2[3], x3[3], n[3];
  SCVT * qNodes = new SCVT[ 3 * quad_size ];
  SCVT * xLocal;

  for ( int i = 0; i < N; i++ ) {

    for ( int currentLeg = 0; currentLeg < legOrder + 1; currentLeg++ ) {

      if ( i == 0 ) {
        start = 0.0;
        end = dt;
      } else if ( i == N - 1 ) {
        start = ( N - 2 ) * dt;
        end = ( N - 1 ) * dt;
      } else {
        start = ( i - 1 ) * dt;
        end = ( i + 1 ) * dt;
      }
      length = end - start;
      ret = 0.0;


      for ( int j = 0; j < nElems; j++ ) {

        this->space->getMesh( )->getNodes( j, x1, x2, x3 );
        this->getQuadratureNodes( x1, x2, x3, quadOrder, qNodes );
        this->space->getOuterElemDOFs( j, &colIndices[0] );
        this->space->getMesh( )->getNormal( j, n );

        SCVT cent[3];
        this->space->getMesh( )->getCentroid( j, cent );


        for ( int k = 0; k < quad_size; k++ ) {
          intPx = quad_points[ 2 * k ];
          intPy = quad_points[ 2 * k + 1 ];
          phi1 = 1.0 - intPx - intPy;
          phi2 = intPx;
          phi3 = intPy;
          xLocal = qNodes + 3 * k;
          //phi1 = 1.0 - intPx;
          //phi2 = intPx - intPy;
          //phi3 = intPy;
          ret = 0.0;

          for ( int q = 0; q < qSize; q++ ) {

            t = start + length * qPoints[q];
            gVal = incWave( t, xLocal, n );
            //g = ( t * t * t * t ) * std::exp( -2.0 * t );

            //g = std::sin( 4.0 * t * t * t * t ) * t * std::exp( -t );
            //g = std::sin( 3.0 * t ) * t * t * std::exp( -t );
            //g = ( std::sin( 2.0 * t ) )*( std::sin( 2.0 * t ) ) * t * 
            //  std::exp( -t ) * 0.5 * sqrt( 3.0 / M_PI ) * cent[2];

            if ( t < 0.0 ) {
              gVal = 0.0;
            }

            //}
            SCVT scaleFact = 1.0;
            if ( legOrder == 0 ) {
              scaleFact = ( std::sqrt( M_PI ) );
            }
            ret += qWeights[q] * gVal * evalBDot( t, i, currentLeg ) * scaleFact;

          }

          ret *= length;
          //ret = 1.0;

          rhs.add( ( i * ( legOrder + 1 ) + currentLeg ) * nNodes +
              colIndices[0], quad_weight[k] * phi1 * ret *
              this->space->getMesh( )->getElemArea( j ) );
          rhs.add( ( i * ( legOrder + 1 ) + currentLeg ) * nNodes +
              colIndices[1], quad_weight[k] * phi2 * ret *
              this->space->getMesh( )->getElemArea( j ) );
          rhs.add( ( i * ( legOrder + 1 ) + currentLeg ) * nNodes +
              colIndices[2], quad_weight[k] * phi3 * ret *
              this->space->getMesh( )->getElemArea( j ) );
        }
      }
    }
  }
  //rhs.scale( 1.0 / ( sqrt( M_PI )*2 ) );
}

template<class LO, class SC>
SC BEIntegratorWave<LO, SC>::uIncNeu(
    SCVT t,
    SCVT* x,
    SCVT *n,
    int type
    ) const {

  switch ( type ) {
    case 0:
      return ( 1.0 / ( std::sqrt( M_PI )*2 ) )*
          ( std::sin( 3.0 * t ) * t * t * std::exp( -t ) );
    case 1:
    {
      SCVT A = 1.0;
      SCVT k[3] = { M_PI, 0.0, 0.0 };
      SCVT omega = std::sqrt( k[0] * k[0] + k[1] * k[1] + k[2] * k[2] );
      SCVT mf = 6.0 * M_PI;
      SCVT mt = 8.0 * M_PI;
      SCVT phi0 = 3.5 * M_PI;
      SCVT kx = DOT3( k, x );
      SCVT kn = DOT3( k, n );
      if ( ( omega * t - mf >= kx ) && ( kx >= omega * t - mt ) ) {
        return A * std::sin( kx + phi0 - omega * t ) * kn;
      } else {
        return 0.0;
      }
    }
    case 2:
      return ( std::sin( 2.0 * t ) * std::sin( 2.0 * t ) * t *
          std::exp( -t ) )*0.5 * std::sqrt( 3 / M_PI ) * x[2];
    case 3:
      return ( std::sin( 2.0 * M_PI * t ) * t * t * t *
          std::exp( -2.0 * t ) )*0.5 * std::sqrt( 3 / M_PI ) * x[2];
    case 4:
    {
      SCVT s[3] = { 0.0, 3.0, 3.0 };
      SCVT R = 0.05;
      SCVT r = DIST3( s, x );


      if ( r - t + 12 * R >= 0.0 && r - t - R <= 0.0 ) {
        SCVT dfdx1 = ( ( M_PI * std::sin( ( M_PI * ( 3 * R - t + r ) ) / R ) *
            ( 2 * s[0] - 2 * x[0] ) ) / ( 8 * R * r ) -
            ( M_PI * std::sin( ( M_PI * ( 3 * R - t + r ) ) /
            ( 2 * R ) )*( 2 * s[0] - 2 * x[0] ) ) / ( 4 * R * r ) ) / r +
            ( ( 2 * s[0] - 2 * x[0] )*( std::cos( ( M_PI * ( 3 * R - t + r ) )
            / R ) / 4 - std::cos( ( M_PI * ( 3 * R - t + r ) ) / ( 2 * R ) ) +
            3 / 4 ) ) / ( 2 * r * std::sqrt( r ) );
        SCVT dfdx2 = ( ( M_PI * std::sin( ( M_PI * ( 3 * R - t + r ) ) / R ) *
            ( 2 * s[1] - 2 * x[1] ) ) / ( 8 * R * r ) -
            ( M_PI * std::sin( ( M_PI * ( 3 * R - t + r ) ) /
            ( 2 * R ) )*( 2 * s[1] - 2 * x[1] ) ) / ( 4 * R * r ) ) / r +
            ( ( 2 * s[1] - 2 * x[1] )*( std::cos( ( M_PI * ( 3 * R - t + r ) )
            / R ) / 4 - std::cos( ( M_PI * ( 3 * R - t + r ) ) / ( 2 * R ) ) +
            3 / 4 ) ) / ( 2 * r * std::sqrt( r ) );
        SCVT dfdx3 = ( ( M_PI * std::sin( ( M_PI * ( 3 * R - t + r ) ) / R ) *
            ( 2 * s[2] - 2 * x[2] ) ) / ( 8 * R * r ) -
            ( M_PI * std::sin( ( M_PI * ( 3 * R - t + r ) ) / ( 2 * R ) )*
            ( 2 * s[2] - 2 * x[2] ) ) / ( 4 * R * r ) ) / r +
            ( ( 2 * s[2] - 2 * x[2] )*( std::cos( ( M_PI * ( 3 * R - t + r ) )
            / R ) / 4 - std::cos( ( M_PI * ( 3 * R - t + r ) ) / ( 2 * R ) ) +
            3 / 4 ) ) / ( 2 * r * std::sqrt( r ) );

        SCVT dfdn = dfdx1 * n[0] + dfdx2 * n[1] + dfdx3 * n[2];
        return dfdn;
      } else {
        return 0;
      }
    }
    default:
      return 0.0;
  }
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::computeElemMatrixHypersingularDisjointP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const {

  SCVT x1[3], x2[3], x3[3], y1[3], y2[3], y3[3];
  SCVT intPxInner, intPyInner, intPxOuter, intPyOuter;
  SCVT phi1x, phi2x, phi3x, phi1y, phi2y, phi3y;

  int qOrderOut = this->quadratureOrderDisjointElems[0];
  int qOrderIn = this->quadratureOrderDisjointElems[1];

  // getting outer quadrature points
  SCVT* outerQuadratureNodes = new SCVT[ quadSizes[ qOrderOut ] * 3 ];
  SCVT* innerQuadratureNodes = new SCVT[ quadSizes[ qOrderIn ] * 3 ];
  this->space->getMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->space->getMesh( )->getNodes( innerElem, y1, y2, y3 );

  // q. points in current elements
  this->getQuadratureNodes( x1, x2, x3, qOrderOut,
      outerQuadratureNodes );
  this->getQuadratureNodes( y1, y2, y3, qOrderIn,
      innerQuadratureNodes );

  int outerPoints = quadSizes[ qOrderOut ];
  int innerPoints = quadSizes[ qOrderIn ];

  // q. points on the reference triangle
  double *qPointsOut = quadPoints[qOrderOut];
  double *qPointsIn = quadPoints[qOrderIn];

  SCVT* outerCurl = this->space->getMesh( )->getCurls( )->getData( ) +
      outerElem * 9;
  SCVT* innerCurl = this->space->getMesh( )->getCurls( )->getData( ) +
      innerElem * 9;
  FullMatrix< LO, SC > curlMatrix( 3, 3 );

  for ( int outerRot = 0; outerRot < 3; outerRot++ ) {
    for ( int innerRot = 0; innerRot < 3; innerRot++ ) {
      curlMatrix.set( outerRot, innerRot, DOT3( ( outerCurl + 3 * outerRot ),
          ( innerCurl + 3 * innerRot ) ) );
    }
  }

  SCVT nInner[3], nOuter[3];
  this->getSpace( )->getMesh( )->getNormal( innerElem, nInner );
  this->getSpace( )->getMesh( )->getNormal( outerElem, nOuter );
  SCVT normalDot = DOT3( nOuter, nInner );
  SCVT innerArea = this->getSpace( )->getMesh( )->getElemArea( innerElem );
  SCVT outerArea = this->getSpace( )->getMesh( )->getElemArea( outerElem );

  matrix.setAll( 0.0 );

  SC kernel = 0.0;
  SC kernelDDot = 0.0;
  SC entry = 0.0;

  for ( int i = 0; i < outerPoints; i++ ) {
    intPxOuter = qPointsOut[ 2 * i ];
    intPyOuter = qPointsOut[ 2 * i + 1 ];
    phi1x = 1.0 - intPxOuter - intPyOuter;
    phi2x = intPxOuter;
    phi3x = intPyOuter;

    for ( int j = 0; j < innerPoints; j++ ) {
      intPxInner = qPointsIn[ 2 * j ];
      intPyInner = qPointsIn[ 2 * j + 1 ];

      // evaluate basis functions
      phi1y = 1.0 - intPxInner - intPyInner;
      phi2y = intPxInner;
      phi3y = intPyInner;


      kernel = evalHypersingularPart1(
          ( outerQuadratureNodes + 3 * i ),
          ( innerQuadratureNodes + 3 * j ) ) *
          quadWeights[ qOrderOut ][ i ] *
          quadWeights[ qOrderIn ][ j ];
      kernelDDot = evalHypersingularPart2(
          ( outerQuadratureNodes + 3 * i ),
          ( innerQuadratureNodes + 3 * j ) ) *
          quadWeights[ qOrderOut ][ i ] *
          quadWeights[ qOrderIn ][ j ];
      entry += kernel;
      matrix.add( 0, 0, kernelDDot * ( normalDot * phi1y * phi1x ) );
      matrix.add( 0, 1, kernelDDot * ( normalDot * phi2y * phi1x ) );
      matrix.add( 0, 2, kernelDDot * ( normalDot * phi3y * phi1x ) );
      matrix.add( 1, 0, kernelDDot * ( normalDot * phi1y * phi2x ) );
      matrix.add( 1, 1, kernelDDot * ( normalDot * phi2y * phi2x ) );
      matrix.add( 1, 2, kernelDDot * ( normalDot * phi3y * phi2x ) );
      matrix.add( 2, 0, kernelDDot * ( normalDot * phi1y * phi3x ) );
      matrix.add( 2, 1, kernelDDot * ( normalDot * phi2y * phi3x ) );
      matrix.add( 2, 2, kernelDDot * ( normalDot * phi3y * phi3x ) );

    }
  }
  matrix.add( curlMatrix, entry );
  matrix.scale( innerArea * outerArea );
  //matrix.print();

  delete [] innerQuadratureNodes;
  delete [] outerQuadratureNodes;

}

#pragma omp declare simd simdlen( DATA_WIDTH )

template<class LO, class SC>
SC BEIntegratorWave<LO, SC>::chebyshevInterpolant(
    const SC a,
    const SC b,
    const SC * const coeffs,
    const SC t
    ) const {

  SC di, dip1, dip2, y;

  dip1 = 0.0;
  di = 0.0;
  y = ( 2.0 * t - a - b ) / ( b - a );
  for ( int i = nChebPoints - 1; 1 <= i; i-- ) {
    dip2 = dip1;
    dip1 = di;
    di = 2.0 * y * dip1 - dip2 + coeffs[i];
  }
  return (y * di - dip1 + 0.5 * coeffs[0] );

}

//template<class LO, class SC>
//SC BEIntegratorWave<LO, SC>::chebyshevInterpolant(
//    const SC a,
//    const SC b,
//    const SC* const coeffs,
//    const SC t
//    ) const {
//
//  SC di, dip1, dip2, y;
//
//  dip1 = 0.0;
//  di = 0.0;
//  y = ( 2.0 * t - a - b ) / ( b - a );
//  for ( int i = nChebPoints - 1; 1 <= i; i-- ) {
//    dip2 = dip1;
//    dip1 = di;
//    di = 2.0 * y * dip1 - dip2 + coeffs[i];
//  }
//  return (y * di - dip1 + 0.5 * coeffs[0] );
//
//}

//template<class LO, class SC>
//void BEIntegratorWave<LO, SC>::computeElemMatrixHypersingularSauterSchwabP1P1(
//    LO outerElem,
//    LO innerElem,
//    FullMatrix<LO, SC>& matrix
//    ) const {
//
//  SCVT* outerCurl = this->space->getMesh( )->getCurls( )->getData( ) +
//      outerElem * 9;
//  SCVT* innerCurl = this->space->getMesh( )->getCurls( )->getData( ) +
//      innerElem * 9;
//  FullMatrix< LO, SC > curlMatrix( 3, 3 );
//
//  for ( int outerRot = 0; outerRot < 3; outerRot++ ) {
//    for ( int innerRot = 0; innerRot < 3; innerRot++ ) {
//      curlMatrix.set( outerRot, innerRot, DOT3( ( outerCurl + 3 * outerRot ),
//          ( innerCurl + 3 * innerRot ) ) );
//    }
//  }
//
//  matrix.setAll( 0.0 );
//
//  SCVT nInner[3], nOuter[3];
//  this->getSpace( )->getMesh( )->getNormal( innerElem, nInner );
//  this->getSpace( )->getMesh( )->getNormal( outerElem, nOuter );
//  SCVT normalDot = DOT3( nOuter, nInner );
//
//  int outerRot;
//  int innerRot;
//  int type;
//
//  this->getSauterSchwabQuadratureType( innerElem, outerElem, type, innerRot,
//      outerRot ); //, innerSwap, outerSwap );
//
//  int qOrder0 = this->quadratureOrder[ 0 ];
//  int qOrder1 = this->quadratureOrder[ 1 ];
//  int qOrder2 = this->quadratureOrder[ 2 ];
//  int qOrder3 = this->quadratureOrder[ 3 ];
//
//  int qSize0 = lineQuadSizes[ qOrder0 ];
//  int qSize1 = lineQuadSizes[ qOrder1 ];
//  int qSize2 = lineQuadSizes[ qOrder2 ];
//  int qSize3 = lineQuadSizes[ qOrder3 ];
//
//  SCVT xRef[2], yRef[2];
//  SCVT x[3], y[3];
//  SCVT w[4];
//  SCVT x1[3], x2[3], x3[3];
//  SCVT y1[3], y2[3], y3[3];
//  SCVT weights[4];
//  SCVT jacobian;
//  SCVT phi1x, phi2x, phi3x, phi1y, phi2y, phi3y;
//  SC kernel = 0.0;
//  SC kernelDDot = 0.0;
//  SC entry = 0.0;
//
//  this->getSpace( )->getMesh( )->getNodes( outerElem, x1, x2, x3 );
//  this->getSpace( )->getMesh( )->getNodes( innerElem, y1, y2, y3 );
//
//  SCVT innerArea = this->getSpace( )->getMesh( )->getElemArea( innerElem );
//  SCVT outerArea = this->getSpace( )->getMesh( )->getElemArea( outerElem );
//
//
//  for ( int dksi = 0; dksi < qSize0; dksi++ ) {
//    w[0] = lineQuadPoints[qOrder0][dksi];
//    weights[0] = lineQuadWeights[qOrder0][dksi];
//    for ( int deta3 = 0; deta3 < qSize1; deta3++ ) {
//      w[3] = lineQuadPoints[qOrder1][deta3];
//      weights[1] = lineQuadWeights[qOrder1][deta3];
//      for ( int deta2 = 0; deta2 < qSize2; deta2++ ) {
//        w[2] = lineQuadPoints[qOrder2][deta2];
//        weights[2] = lineQuadWeights[qOrder2][deta2];
//        for ( int deta1 = 0; deta1 < qSize3; deta1++ ) {
//          w[1] = lineQuadPoints[qOrder3][deta1];
//          weights[3] = lineQuadWeights[qOrder3][deta1];
//
//          switch ( type ) {
//              // identical panels
//            case 0:
//              for ( int simplex = 0; simplex < 6; simplex++ ) {
//                this->cube2triIdentical( w, simplex, xRef, yRef, jacobian );
//                this->tri2element( x1, x2, x3, xRef, x );
//                this->tri2element( y1, y2, y3, yRef, y );
//                evalHypersingular( x[0], x[1], x[2], y[0], y[1], y[2], 
//                    kernel, kernelDDot);
//                kernel *=  weights[0] * weights[1] * weights[2] * weights[3] *
//                    jacobian;
//                kernelDDot *=  weights[0] * weights[1] * weights[2] * weights[3] *
//                    jacobian;
////                kernel = evalHypersingularPart1( x, y ) *
////                    weights[0] * weights[1] * weights[2] * weights[3] *
////                    jacobian;
////                kernelDDot = evalHypersingularPart2( x, y ) *
////                    weights[0] * weights[1] * weights[2] * weights[3] *
////                    jacobian;
//
//                entry += kernel;
//                phi1x = 1.0 - xRef[0];
//                phi2x = xRef[0] - xRef[1];
//                phi3x = xRef[1];
//                phi1y = 1.0 - yRef[0];
//                phi2y = yRef[0] - yRef[1];
//                phi3y = yRef[1];
//
//                matrix.add( 0, 0, kernelDDot * ( normalDot * phi1y * phi1x ) );
//                matrix.add( 0, 1, kernelDDot * ( normalDot * phi2y * phi1x ) );
//                matrix.add( 0, 2, kernelDDot * ( normalDot * phi3y * phi1x ) );
//                matrix.add( 1, 0, kernelDDot * ( normalDot * phi1y * phi2x ) );
//                matrix.add( 1, 1, kernelDDot * ( normalDot * phi2y * phi2x ) );
//                matrix.add( 1, 2, kernelDDot * ( normalDot * phi3y * phi2x ) );
//                matrix.add( 2, 0, kernelDDot * ( normalDot * phi1y * phi3x ) );
//                matrix.add( 2, 1, kernelDDot * ( normalDot * phi2y * phi3x ) );
//                matrix.add( 2, 2, kernelDDot * ( normalDot * phi3y * phi3x ) );
//              }
//              break;
//
//              // common edge
//            case 1:
//              for ( int simplex = 0; simplex < 5; simplex++ ) {
//                this->cube2triCommonEdge( w, simplex, xRef, yRef, jacobian );
//                // <("): outer element swaps first two nodes, so that the edges agree
//                this->tri2element( x1, x2, x3, xRef, outerRot, x, true );
//                this->tri2element( y1, y2, y3, yRef, innerRot, y );
//                evalHypersingular( x[0], x[1], x[2], y[0], y[1], y[2], 
//                    kernel, kernelDDot);
//                kernel *=  weights[0] * weights[1] * weights[2] * weights[3] *
//                    jacobian;
//                kernelDDot *=  weights[0] * weights[1] * weights[2] * weights[3] *
//                    jacobian;
////                kernel = jacobian * evalHypersingularPart1( x, y ) *
////                    weights[0] * weights[1] * weights[2] * weights[3];
////                kernelDDot = evalHypersingularPart2( x, y ) *
////                    weights[0] * weights[1] * weights[2] * weights[3] *
////                    jacobian;
//
//                entry += kernel;
//                phi1x = 1.0 - xRef[0];
//                phi2x = xRef[0] - xRef[1];
//                phi3x = xRef[1];
//                phi1y = 1.0 - yRef[0];
//                phi2y = yRef[0] - yRef[1];
//                phi3y = yRef[1];
//                // outer indices swap due to comment <(")
//                matrix.add( this->mod3( 1 + outerRot ), this->mod3( innerRot ),
//                    kernelDDot * ( normalDot * phi1y * phi1x ) );
//                matrix.add( this->mod3( 1 + outerRot ),
//                    this->mod3( 1 + innerRot ), kernelDDot *
//                    ( normalDot * phi2y * phi1x ) );
//                matrix.add( this->mod3( 1 + outerRot ), this->mod3( 2 +
//                    innerRot ), kernelDDot * ( normalDot * phi3y * phi1x ) );
//                matrix.add( this->mod3( outerRot ), this->mod3( innerRot ),
//                    kernelDDot * ( normalDot * phi1y * phi2x ) );
//                matrix.add( this->mod3( outerRot ), this->mod3( 1 + innerRot ),
//                    kernelDDot * ( normalDot * phi2y * phi2x ) );
//                matrix.add( this->mod3( outerRot ), this->mod3( 2 + innerRot ),
//                    kernelDDot * ( normalDot * phi3y * phi2x ) );
//                matrix.add( this->mod3( 2 + outerRot ), this->mod3( innerRot ),
//                    kernelDDot * ( normalDot * phi1y * phi3x ) );
//                matrix.add( this->mod3( 2 + outerRot ), this->mod3( 1 +
//                    innerRot ), kernelDDot * ( normalDot * phi2y * phi3x ) );
//                matrix.add( this->mod3( 2 + outerRot ), this->mod3( 2 +
//                    innerRot ), kernelDDot * ( normalDot * phi3y * phi3x ) );
//              }
//              break;
//
//              // common vertex
//            case 2:
//              for ( int simplex = 0; simplex < 2; simplex++ ) {
//                this->cube2triCommonVertex( w, simplex, xRef, yRef, jacobian );
//                this->tri2element( x1, x2, x3, xRef, outerRot, x );
//                this->tri2element( y1, y2, y3, yRef, innerRot, y );
//                evalHypersingular( x[0], x[1], x[2], y[0], y[1], y[2], 
//                    kernel, kernelDDot);
//                kernel *=  weights[0] * weights[1] * weights[2] * weights[3] *
//                    jacobian;
//                kernelDDot *=  weights[0] * weights[1] * weights[2] * weights[3] *
//                    jacobian;
////                kernel = evalHypersingularPart1( x, y ) *
////                    weights[0] * weights[1] * weights[2] * weights[3] *
////                    jacobian;
////                kernelDDot = evalHypersingularPart2( x, y ) *
////                    weights[0] * weights[1] * weights[2] * weights[3] *
////                    jacobian;
//                entry += kernel;
//                phi1x = 1.0 - xRef[0];
//                phi2x = xRef[0] - xRef[1];
//                phi3x = xRef[1];
//                phi1y = 1.0 - yRef[0];
//                phi2y = yRef[0] - yRef[1];
//                phi3y = yRef[1];
//                matrix.add( this->mod3( outerRot ), this->mod3( innerRot ),
//                    kernelDDot * ( normalDot * phi1y * phi1x ) );
//                matrix.add( this->mod3( outerRot ), this->mod3( 1 + innerRot ),
//                    kernelDDot * ( normalDot * phi2y * phi1x ) );
//                matrix.add( this->mod3( outerRot ), this->mod3( 2 + innerRot ),
//                    kernelDDot * ( normalDot * phi3y * phi1x ) );
//                matrix.add( this->mod3( 1 + outerRot ), this->mod3( innerRot ),
//                    kernelDDot * ( normalDot * phi1y * phi2x ) );
//                matrix.add( this->mod3( 1 + outerRot ), this->mod3( 1 +
//                    innerRot ), kernelDDot * ( normalDot * phi2y * phi2x ) );
//                matrix.add( this->mod3( 1 + outerRot ), this->mod3( 2 +
//                    innerRot ), kernelDDot * ( normalDot * phi3y * phi2x ) );
//                matrix.add( this->mod3( 2 + outerRot ), this->mod3( innerRot ),
//                    kernelDDot * ( normalDot * phi1y * phi3x ) );
//                matrix.add( this->mod3( 2 + outerRot ), this->mod3( 1 +
//                    innerRot ), kernelDDot * ( normalDot * phi2y * phi3x ) );
//                matrix.add( this->mod3( 2 + outerRot ), this->mod3( 2 +
//                    innerRot ), kernelDDot * ( normalDot * phi3y * phi3x ) );
//              }
//              break;
//
//              // disjoint triangles  
//            case 3:
//              this->cube2triDisjoint( w, xRef, yRef, jacobian );
//              this->tri2element( x1, x2, x3, xRef, x );
//              this->tri2element( y1, y2, y3, yRef, y );
//              evalHypersingular( x[0], x[1], x[2], y[0], y[1], y[2], 
//                    kernel, kernelDDot);
//                kernel *=  weights[0] * weights[1] * weights[2] * weights[3] *
//                    jacobian;
//                kernelDDot *=  weights[0] * weights[1] * weights[2] * weights[3] *
//                    jacobian;
////              kernel = evalHypersingularPart1( x, y ) *
////                  weights[0] * weights[1] * weights[2] * weights[3] *
////                  jacobian;
////              kernelDDot = evalHypersingularPart2( x, y ) *
////                  weights[0] * weights[1] * weights[2] * weights[3] *
////                  jacobian;
//              entry += kernel;
//              phi1x = 1.0 - xRef[0];
//              phi2x = xRef[0] - xRef[1];
//              phi3x = xRef[1];
//              phi1y = 1.0 - yRef[0];
//              phi2y = yRef[0] - yRef[1];
//              phi3y = yRef[1];
//              matrix.add( 0, 0, kernelDDot * ( normalDot * phi1y * phi1x ) );
//              matrix.add( 0, 1, kernelDDot * ( normalDot * phi2y * phi1x ) );
//              matrix.add( 0, 2, kernelDDot * ( normalDot * phi3y * phi1x ) );
//              matrix.add( 1, 0, kernelDDot * ( normalDot * phi1y * phi2x ) );
//              matrix.add( 1, 1, kernelDDot * ( normalDot * phi2y * phi2x ) );
//              matrix.add( 1, 2, kernelDDot * ( normalDot * phi3y * phi2x ) );
//              matrix.add( 2, 0, kernelDDot * ( normalDot * phi1y * phi3x ) );
//              matrix.add( 2, 1, kernelDDot * ( normalDot * phi2y * phi3x ) );
//              matrix.add( 2, 2, kernelDDot * ( normalDot * phi3y * phi3x ) );
//              break;
//          }
//        }
//      }
//    }
//  }
//
//  matrix.add( curlMatrix, entry );
//  matrix.scale( 4.0 * innerArea * outerArea );
//
//  //std::cout << type << std::endl;
//  matrix.print();
//}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::computeElemMatrixHypersingularSauterSchwabP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const {

  SCVT* outerCurl = this->space->getLeftMesh( )->getCurls( )->getData( ) +
      outerElem * 9;
  SCVT* innerCurl = this->space->getRightMesh( )->getCurls( )->getData( ) +
      innerElem * 9;
  FullMatrix< LO, SC > curlMatrix( 3, 3 );

  for ( int outerRot = 0; outerRot < 3; outerRot++ ) {
    for ( int innerRot = 0; innerRot < 3; innerRot++ ) {
      curlMatrix.set( outerRot, innerRot, DOT3( ( outerCurl + 3 * outerRot ),
          ( innerCurl + 3 * innerRot ) ) );
    }
  }

  matrix.setAll( 0.0 );

  SCVT nInner[3], nOuter[3];
  this->getSpace( )->getRightMesh( )->getNormal( innerElem, nInner );
  this->getSpace( )->getLeftMesh( )->getNormal( outerElem, nOuter );
  SCVT normalDot = DOT3( nOuter, nInner );

  int outerRot;
  int innerRot;
  int type;

  this->getSauterSchwabQuadratureType( innerElem, outerElem, type, innerRot,
      outerRot ); //, innerSwap, outerSwap );

  int qSize0 = lineQuadSizes[ this->quadratureOrder[ 0 ] ];
  int qSize1 = lineQuadSizes[ this->quadratureOrder[ 1 ] ];
  int qSize2 = lineQuadSizes[ this->quadratureOrder[ 2 ] ];
  int qSize3 = lineQuadSizes[ this->quadratureOrder[ 3 ] ];
  int totalSize = qSize0 * qSize1 * qSize2 * qSize3;

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->getSpace( )->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->getSpace( )->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );

  int nSimplex;
  SCVT ** jacobianP;
  SCVT ** x1refP;
  SCVT ** x2refP;
  SCVT ** y1refP;
  SCVT ** y2refP;

  this->getReferenceSauterSchwabData( type, nSimplex, jacobianP, x1refP, x2refP,
      y1refP, y2refP );

  SCVT * jacobian;
  SCVT * x1ref;
  SCVT * x2ref;
  SCVT * y1ref;
  SCVT * y2ref;

  SC entry11 = 0.0;
  SC entry21 = 0.0;
  SC entry31 = 0.0;
  SC entry12 = 0.0;
  SC entry22 = 0.0;
  SC entry32 = 0.0;
  SC entry13 = 0.0;
  SC entry23 = 0.0;
  SC entry33 = 0.0;

  SC kernel = 0.0;
  SC kernelDDot = 0.0;
  SC entry = 0.0;

  SCVT mult, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y;
  int i;

  SCVT * x1ss = this->x1ss;
  SCVT * x2ss = this->x2ss;
  SCVT * x3ss = this->x3ss;
  SCVT * y1ss = this->y1ss;
  SCVT * y2ss = this->y2ss;
  SCVT * y3ss = this->y3ss;

  const int align = DATA_ALIGN;

  for ( int simplex = 0; simplex < nSimplex; ++simplex ) {
    this->updateSauterSchwabQuadratureNodes( type, simplex, x1, x2, x3,
        outerRot, y1, y2, y3, innerRot );
    jacobian = jacobianP[ simplex ];
    x1ref = x1refP[ simplex ];
    x2ref = x2refP[ simplex ];
    y1ref = y1refP[ simplex ];
    y2ref = y2refP[ simplex ];
/*
#pragma omp simd linear( i : 1 ) reduction( + : entry11, entry21, entry31,\
entry12, entry22, entry32, entry13, entry23, entry33, entry ) \
aligned( x1ss, x2ss, x3ss, y1ss, y2ss, y3ss, jacobian : align ) \
aligned( x1ref, x2ref, y1ref, y2ref : align )
*/
    for ( i = 0; i < totalSize; ++i ) {
      this->evalHypersingular( x1ss[ i ], x2ss[ i ],
          x3ss[ i ], y1ss[ i ], y2ss[ i ], y3ss[ i ], kernel, kernelDDot );
      kernel *= jacobian[ i ];
      kernelDDot *= jacobian[ i ];
      
      entry += kernel;

      phi1x = (SCVT) 1.0 - x1ref[ i ];
      phi2x = x1ref[ i ] - x2ref[ i ];
      phi3x = x2ref[ i ];
      phi1y = (SCVT) 1.0 - y1ref[ i ];
      phi2y = y1ref[ i ] - y2ref[ i ];
      phi3y = y2ref[ i ];

      entry11 += kernelDDot * normalDot * phi1x * phi1y;
      entry21 += kernelDDot * normalDot * phi2x * phi1y;
      entry31 += kernelDDot * normalDot * phi3x * phi1y;
      entry12 += kernelDDot * normalDot * phi1x * phi2y;
      entry22 += kernelDDot * normalDot * phi2x * phi2y;
      entry32 += kernelDDot * normalDot * phi3x * phi2y;
      entry13 += kernelDDot * normalDot * phi1x * phi3y;
      entry23 += kernelDDot * normalDot * phi2x * phi3y;
      entry33 += kernelDDot * normalDot * phi3x * phi3y;
    }
  }

  // first two y indices swapped for common edge!
  if ( type == 1 ) {
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ),
        entry11 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ),
        entry21 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ),
        entry31 );
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ),
        entry12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ),
        entry22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ),
        entry32 );
  } else {
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ),
        entry11 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ),
        entry21 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ),
        entry31 );
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ),
        entry12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ),
        entry22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ),
        entry32 );
  }
  matrix.set( this->mod3( outerRot ), this->mod3( 2 + innerRot ),
      entry13 );
  matrix.set( this->mod3( 1 + outerRot ), this->mod3( 2 + innerRot ),
      entry23 );
  matrix.set( this->mod3( 2 + outerRot ), this->mod3( 2 + innerRot ),
      entry33 );

  matrix.add( curlMatrix, entry );

  SCVT areaMult = (SCVT) 4.0 *
      this->getSpace( )->getRightMesh( )->getElemArea( innerElem ) *
      this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );
  matrix.scale( areaMult );

    

  //std::cout << type << std::endl;
  //matrix.print();
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::doubleLayerPotentialP1(
    const SCVT* xCoord,
    LO numPoints,
    SCVT t,
    const Vector<LO, SC> & density,
    Vector<LO, SC> & values
    ) const {

  BESpaceTime<LO, SC> *spaceTime =
      static_cast<BESpaceTime<LO, SC>*> ( this->space );

  LO nElems = this->space->getMesh( )->getNElements( );
  LO nNodes = this->space->getMesh( )->getNNodes( );

  //  std::cout << nElems << std::endl;
  //  std::cout << nNodes << std::endl;

  //  SC * x;
  SC * y;
  SC alpha1, alpha2;
  SCVT x1[3], x2[3], x3[3], n[ 3 ];
  SCVT normxy;
  LO interval;
  SCVT dot;
  SCVT retardedT;

  LO ind[3];
  SC val;
  SCVT intPx, intPy;

  // variables for Gaussian quadrature 
  int qOrder = 5;
  int quad_size = quadSizes[qOrder];
  double *quad_points = quadPoints[qOrder];
  double *quad_weight = quadWeights[qOrder];
  SCVT *qNodes = new SCVT[3 * quad_size];

  SCVT * linValues = new SCVT[3 * quad_size];
  SCVT phi;
  for ( int i = 0; i < quad_size; i++ ) {
    intPx = quad_points[ 2 * i ];
    intPy = quad_points[ 2 * i + 1 ];
    linValues[i] = 1.0 - intPx - intPy;
    linValues[i + quad_size] = intPx;
    linValues[i + 2 * quad_size] = intPy;
  }

  int legOrder = spaceTime->getLegendreOrder( );
  //  std::cout << legOrder << std::endl;


  for ( LO i = 0; i < numPoints; i++ ) {
    val = 0.0;
    const SCVT * x = xCoord + 3 * i;
    //std::cout << i << std::endl;
    for ( LO j = 0; j < nElems; j++ ) {
      this->space->getMesh( )->getNodes( j, x1, x2, x3 );
      this->space->getMesh( )->getElement( j, ind );
      this->getQuadratureNodes( x1, x2, x3, qOrder, qNodes );
      this->space->getMesh( )->getNormal( j, n );

      for ( int spatialDof = 0; spatialDof < 3; ++spatialDof ) {
        for ( int temporalDof = 0; temporalDof < legOrder + 1; ++temporalDof ) {

          for ( int k = 0; k < quad_size; k++ ) {
            phi = linValues[ k + spatialDof * quad_size ];
            y = qNodes + 3 * k;
            normxy = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) +
                ( x[1] - y[1] )*( x[1] - y[1] ) +
                ( x[2] - y[2] )*( x[2] - y[2] ) );
            dot = n[0]*( x[0] - y[0] ) + n[1]*( x[1] - y[1] ) +
                n[2]*( x[2] - y[2] );
            retardedT = t - normxy;
            if ( retardedT <= 0 ) {
              val = 0.0;
              break;
            }
            interval = std::floor( retardedT / spaceTime->getDt( ) );

            alpha1 = density.get( interval * ( legOrder + 1 ) * nNodes +
                temporalDof * nNodes + ind[spatialDof] );
            alpha2 = density.get( ( interval + 1 ) * ( legOrder + 1 ) * nNodes +
                temporalDof * nNodes + ind[spatialDof] );

            //            std::cout << interval << std::endl;
            //            std::cout << normxy << std::endl;
            //            std::cout << dot << std::endl;
            //            std::cout << retardedT << std::endl;
            //            std::cou  t << alpha1 << std::endl;
            //            std::cout << alpha2 << std::endl;
            //            std::cout << spaceTime->getMesh( )->getElemArea( j ) << std::endl;
            //            std::cout << std::endl;

            val += spaceTime->getMesh( )->getElemArea( j ) *
                quad_weight[k]*( dot / normxy ) *
                ( phi * ( alpha1 * evalB( retardedT, interval, temporalDof ) +
                alpha2 * evalB( retardedT, interval + 1, temporalDof ) )
                / ( normxy * normxy ) +
                phi * ( alpha1 * evalBDot( retardedT, interval, temporalDof ) +
                alpha2 * evalBDot( retardedT, interval + 1, temporalDof ) )
                / normxy );
          }
        }
      }
    }

    values.set( i, -PI_FACT * val );
    //std::cout << val << std::endl;
  }
  delete [] qNodes;
}

}

#else

// use experimental basis functions which are not partition of unity


namespace bem4i {

template<class LO, class SC>
BEIntegratorWave<LO, SC>::BEIntegratorWave( ) {
}

template<class LO, class SC>
BEIntegratorWave<LO, SC>::BEIntegratorWave( const BEIntegratorWave& orig ) {
}

template<class LO, class SC>
BEIntegratorWave<LO, SC>::BEIntegratorWave(
    BESpace<LO, SC>* space,
    int* quadratureOrder,
    int timeQuadOrder,
    int nPre,
    int nPos,
    int nChebIntervals,
    int nChebPoints
    ) {

  this->space = space;
  this->quadratureOrder = quadratureOrder;
  this->timeQuadOrder = timeQuadOrder;
  this->quadrature = SauterSchwab;

  N = static_cast<BESpaceTime<LO, SC>*> ( this->space )->getNTimeSteps( );
  dt = static_cast<BESpaceTime<LO, SC>*> ( this->space )->getDt( );

  // initialize arrays for Chebyshev interpolation of temporal basis functions
  this->nChebIntervals = nChebIntervals;
  this->nChebPoints = nChebPoints;
  chebIntervalStarts = new SC[nChebIntervals + 1];
  psiChebCoeffs = new SCVT*[nChebIntervals];
  psiTildeChebCoeffs = new SCVT*[nChebIntervals];
  for ( int i = 0; i < nChebIntervals; i++ ) {
    psiChebCoeffs[i] = new SCVT[nChebPoints];
    psiTildeChebCoeffs[i] = new SCVT[nChebPoints];
  }

  this->nPre = nPre;
  this->nPos = nPos;

}

template<class LO, class SC>
BEIntegratorWave<LO, SC>::~BEIntegratorWave( ) {
  delete [] chebIntervalStarts;
  for ( int i = 0; i < nChebIntervals; i++ ) {
    delete [] psiChebCoeffs[i];
    delete [] psiTildeChebCoeffs[i];
  }
  delete [] psiChebCoeffs;
  delete [] psiTildeChebCoeffs;
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::computeElemMatrix1Layer(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const {

  if ( this->space->getTestFunctionType( ) == p0 && this->space->getAnsatzFunctionType( ) == p0 ) {
    switch ( this->quadrature ) {
      case SauterSchwab:
        this->computeElemMatrix1LayerSauterSchwabP0P0( outerElem, innerElem, matrix );
        break;
      case Steinbach:
        std::cout << "Not implemented!" << std::endl;
        break;
    }
  } else if ( this->space->getTestFunctionType( ) == p1 && this->space->getAnsatzFunctionType( ) == p1 ) {
    //computeElemMatrix1LayerP1P1( outerElem, innerElem, matrix );
    std::cout << "Not implemented!" << std::endl;
  } else {
    std::cout << "Not implemented!" << std::endl;
  }
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::computeElemMatrix2Layer(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const {

  if ( this->space->getTestFunctionType( ) == p0 && this->space->getAnsatzFunctionType( ) == p1 ) {
    switch ( this->quadrature ) {
      case SauterSchwab:
        this->computeElemMatrix2LayerSauterSchwabP0P1( outerElem, innerElem, matrix );
        break;
      case Steinbach:
        std::cout << "Not implemented!" << std::endl;
        break;
    }
  } else if ( this->space->getTestFunctionType( ) == p0 && this->space->getAnsatzFunctionType( ) == p0 ) {
    //computeElemMatrix2LayerP0P0( outerElem, innerElem, matrix );
    std::cout << "Not implemented!" << std::endl;
  } else if ( this->space->getTestFunctionType( ) == p1 && this->space->getAnsatzFunctionType( ) == p1 ) {
    //computeElemMatrix2LayerP1P1( outerElem, innerElem, matrix );
    std::cout << "Not implemented!" << std::endl;
  } else {
    std::cout << "Not implemented!" << std::endl;
  }
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::computeElemMatrixHypersingular(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const {

  if ( this->space->getTestFunctionType( ) == p1 && this->space->getAnsatzFunctionType( ) == p1 ) {
    computeElemMatrixHypersingularSauterSchwabP1P1( outerElem, innerElem, matrix );
  } else {
    std::cout << "Not implemented!" << std::endl;
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::evalErf( SCVT t, SCVT a, SCVT b ) const {
  SCVT arg = 2.0 * ( t - a ) / ( b - a ) - 1.0;
  if ( arg <= -1.0 ) {
    return 0.0;
  } else if ( arg >= 1.0 ) {
    return 1.0;
  } else {
    return 0.5 * ( erf( 2.0 * atanh( arg ) ) + 1.0 );
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::evalErfDot( SCVT t, SCVT a, SCVT b ) const {

  SCVT arg = ( a + b - 2.0 * t ) / ( a - b );
  if ( arg <= -1.0 ) {
    return 0.0;
  } else if ( arg >= 1.0 ) {
    return 0.0;
  } else {
    SCVT ath = atanh( arg );
    return ( ( ( a - b ) * std::exp( -4.0 * ath * ath ) ) / ( sqrt( M_PI )*( a - t )*( b - t ) ) );
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::evalErfDDot( SCVT t, SCVT a, SCVT b ) const {
  SCVT arg = ( a + b - 2.0 * t ) / ( a - b );
  if ( arg <= -1.0 ) {
    return 0.0;
  } else if ( arg >= 1.0 ) {
    return 0.0;
  } else {
    SCVT ath = atanh( arg );
    return ( ( ( a - b ) * std::exp( -4.0 * ath * ath ) * ( -4.0 * ( a - b ) * ath + a + b - 2.0 * t ) ) / ( sqrt( M_PI )*( a - t )*( a - t )*( b - t )*( b - t ) ) );
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::evalPUM( SCVT t, LO i ) const {

  SCVT tIMinus1 = dt * ( i - 1 );
  SCVT tI = dt * i;
  SCVT tIPlus1 = dt * ( i + 1 );

  if ( t <= tI ) {
    return evalErf( t, tIMinus1, tI );
  } else {
    return 1.0 - evalErf( t, tI, tIPlus1 );
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::evalPUMDot( SCVT t, LO i ) const {

  SCVT tIMinus1 = dt * ( i - 1 );
  SCVT tI = dt * i;
  SCVT tIPlus1 = dt * ( i + 1 );

  if ( fabs( t - tI ) < EPS ) {
    return 0.0;
  } else if ( t < tI ) {
    return evalErfDot( t, tIMinus1, tI );
  } else {
    return (-evalErfDot( t, tI, tIPlus1 ) );
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::evalPUMDDot( SCVT t, LO i ) const {

  SCVT tIMinus1 = dt * ( i - 1 );
  SCVT tI = dt * i;
  SCVT tIPlus1 = dt * ( i + 1 );

  if ( fabs( t - tI ) < EPS ) {
    return 0.0;
  } else if ( t < tI ) {
    return evalErfDDot( t, tIMinus1, tI );
  } else {
    return (-evalErfDDot( t, tI, tIPlus1 ) );
  }
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::evalB( SCVT t, LO i, LO order ) const {

  SCVT tIMinus1 = dt * ( i - 1 );
  SCVT tI = dt * i;
  SCVT tIPlus1 = dt * ( i + 1 );

  return evalPUM( t, i ) * evalLegendrePolynomial( 2.0 * ( t - tIMinus1 ) / ( tIPlus1 - tIMinus1 ) - 1.0, order );
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::evalBDot( SCVT t, LO i, LO order ) const {

  // let's not make it too complicated for now - suppose Legendres are of order 0
  return evalPUMDot( t, i );
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::evalBDDot( SCVT t, LO i, LO order ) const {

  // let's not make it too complicated for now - suppose Legendres are of order 0
  return evalPUMDDot( t, i );
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::integrate1Layer(
    SCVT start,
    SCVT end,
    SCVT r
    ) const {
  SCVT ret = 0.0;
  LO qSize = lineQuadSizes[timeQuadOrder];
  SCVT *qPoints = lineQuadPoints[timeQuadOrder];
  SCVT *qWeights = lineQuadWeights[timeQuadOrder];

  SCVT length = end - start;
  SCVT t = 0.0;

  for ( int i = 0; i < qSize; i++ ) {
    t = start + length * qPoints[i];
    ret += qWeights[i] * evalBDot( t - r, currBasisF ) *
        evalB( t, currTestF );
  }
  ret *= length;
  return ret;
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::integrateHypersingularPart2(
    SCVT start,
    SCVT end,
    SCVT r
    ) const {
  SCVT ret = 0.0;
  LO qSize = lineQuadSizes[timeQuadOrder];
  SCVT *qPoints = lineQuadPoints[timeQuadOrder];
  SCVT *qWeights = lineQuadWeights[timeQuadOrder];

  SCVT length = end - start;
  SCVT t = 0.0;

  for ( int i = 0; i < qSize; i++ ) {
    t = start + length * qPoints[i];
    ret += qWeights[i] * evalBDDot( t - r, currBasisF ) *
        evalBDot( t, currTestF );
  }
  ret *= length;
  return ret;
}

template<class LO, class SC>
typename bem4i::BEIntegratorWave<LO, SC>::SCVT BEIntegratorWave<LO, SC>::integrateHypersingularPart1(
    SCVT start,
    SCVT end,
    SCVT r
    ) const {
  // evaluates psi tilde
  SCVT ret = 0.0;
  LO qSize = lineQuadSizes[timeQuadOrder];
  SCVT *qPoints = lineQuadPoints[timeQuadOrder];
  SCVT *qWeights = lineQuadWeights[timeQuadOrder];

  SCVT length = end - start;
  SCVT t = 0.0;

  for ( int i = 0; i < qSize; i++ ) {
    t = start + length * qPoints[i];
    ret += qWeights[i] * evalB( t - r, currBasisF ) *
        evalBDot( t, currTestF );
  }
  ret *= length;
  return ret;
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::getDirichletRHS( Vector<LO, SC> &rhs ) const {
  LO nElems = this->space->getMesh( )->getNElements( );
  rhs.resize( nElems * N );

  SCVT ret = 0.0;
  LO qSize = lineQuadSizes[timeQuadOrder];
  SCVT *qPoints = lineQuadPoints[timeQuadOrder];
  SCVT *qWeights = lineQuadWeights[timeQuadOrder];

  SCVT t = 0.0;
  SC g;
  SCVT length, start, end;

  for ( int i = 0; i < N; i++ ) {

    if ( i == 0 ) {
      start = 0.0;
      end = dt;
    } else if ( i == N - 1 ) {
      start = ( N - 2 ) * dt;
      end = ( N - 1 ) * dt;
    } else {
      start = ( i - 1 ) * dt;
      end = ( i + 1 ) * dt;
    }
    length = end - start;
    ret = 0.0;
    for ( int q = 0; q < qSize; q++ ) {
      t = start + length * qPoints[q];
      g = ( 4.0 * t * t * t - 2.0 * t * t * t * t ) * std::exp( -2.0 * t );
      ret += qWeights[q] * evalB( t, i, 0 ) * g;
    }
    ret *= length;
    for ( LO j = 0; j < nElems; j++ ) {
      rhs.set( i * nElems + j, this->space->getMesh( )->getElemArea( j ) * ret );
    }
  }
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::getNeumannRHS( Vector<LO, SC> &rhs ) const {

  LO nElems = this->space->getMesh( )->getNElements( );
  LO nNodes = this->space->getMesh( )->getNNodes( );
  //  rhs.resize( nNodes * ( N + nPos  ) );
  rhs.resize( nNodes * ( N + nPos + nPre ) );
  rhs.setAll( 0.0 );

  SCVT ret = 0.0;
  LO qSize = lineQuadSizes[timeQuadOrder];
  SCVT *qPoints = lineQuadPoints[timeQuadOrder];
  SCVT *qWeights = lineQuadWeights[timeQuadOrder];

  // variables for Gaussian quadrature 
  int quad_size = quadSizes[7];
  double *quad_points = quadPoints[7];
  double *quad_weight = quadWeights[7];

  SCVT t = 0.0;
  SC g;
  SCVT length, start, end;
  //SCVT x1[3], x2[3], x3[3];
  //SCVT* qNodes = new SCVT[ quadSizes[ this->quadratureOrder[ 0 ] ] * 3 ];
  SCVT phi1, phi2, phi3, intPx, intPy;
  vector<LO> colIndices( 3 );

  for ( int i = 0; i < N + nPos + nPre; i++ ) {


    start = ( -nPre + i - 1 ) * dt;
    end = ( -nPre + i + 1 ) * dt;
    length = end - start;
    ret = 0.0;


    for ( int j = 0; j < nElems; j++ ) {

      //this->space->getMesh( )->getNodes( j, x1, x2, x3 );
      //this->getQuadratureNodes( x1, x2, x3, this->quadratureOrder[ 0 ], qNodes );
      this->space->getOuterElemDOFs( j, &colIndices[0] );


      for ( int k = 0; k < quad_size; k++ ) {
        intPx = quad_points[ 2 * k ];
        intPy = quad_points[ 2 * k + 1 ];
        phi1 = 1.0 - intPx - intPy;
        phi2 = intPx;
        phi3 = intPy;
        //phi1 = 1.0 - intPx;
        //phi2 = intPx - intPy;
        //phi3 = intPy;
        ret = 0.0;

        for ( int q = 0; q < qSize; q++ ) {

          t = start + length * qPoints[q];

          //g = ( t * t * t * t ) * std::exp( -2.0 * t );

          //g = std::sin( 4.0 * t * t * t * t ) * t * std::exp( -t );
          g = std::sin( 3.0 * t ) * t * t * std::exp( -t );
          if ( t < dt ) {
            g = 0.0;
          }
          ret += qWeights[q] * g * evalBDot( t, -nPre + i, 0 ) / ( sqrt( M_PI )*2.0 );
        }
        ret *= length;
        //ret = 1.0;

        rhs.add( i * nNodes + colIndices[0], quad_weight[k] * phi1 * ret * this->space->getMesh( )->getElemArea( j ) );
        rhs.add( i * nNodes + colIndices[1], quad_weight[k] * phi2 * ret * this->space->getMesh( )->getElemArea( j ) );
        rhs.add( i * nNodes + colIndices[2], quad_weight[k] * phi3 * ret * this->space->getMesh( )->getElemArea( j ) );
      }
    }
  }
  //rhs.scale(1.0/(sqrt(M_PI)*2));

}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::computeElemMatrixHypersingularSauterSchwabP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const {

  SCVT* outerCurl = this->space->getMesh( )->getCurls( )->getData( ) + outerElem * 9;
  SCVT* innerCurl = this->space->getMesh( )->getCurls( )->getData( ) + innerElem * 9;
  FullMatrix< LO, SC > curlMatrix( 3, 3 );

  for ( int outerRot = 0; outerRot < 3; outerRot++ ) {
    for ( int innerRot = 0; innerRot < 3; innerRot++ ) {
      curlMatrix.set( outerRot, innerRot, DOT3( ( outerCurl + 3 * outerRot ), ( innerCurl + 3 * innerRot ) ) );
    }
  }

  matrix.setAll( 0.0 );

  SCVT nInner[3], nOuter[3];
  this->getSpace( )->getMesh( )->getNormal( innerElem, nInner );
  this->getSpace( )->getMesh( )->getNormal( outerElem, nOuter );
  SCVT normalDot = DOT3( nOuter, nInner );

  int outerRot;
  int innerRot;
  int type;

  this->getSauterSchwabQuadratureType( innerElem, outerElem, type, innerRot, outerRot ); //, innerSwap, outerSwap );

  int qOrder0 = this->quadratureOrder[ 0 ];
  int qOrder1 = this->quadratureOrder[ 1 ];
  int qOrder2 = this->quadratureOrder[ 2 ];
  int qOrder3 = this->quadratureOrder[ 3 ];

  int qSize0 = lineQuadSizes[ qOrder0 ];
  int qSize1 = lineQuadSizes[ qOrder1 ];
  int qSize2 = lineQuadSizes[ qOrder2 ];
  int qSize3 = lineQuadSizes[ qOrder3 ];

  SCVT xRef[2], yRef[2];
  SCVT x[3], y[3];
  SCVT w[4];
  SCVT x1[3], x2[3], x3[3];
  SCVT y1[3], y2[3], y3[3];
  SCVT weights[4];
  SCVT jacobian;
  SCVT phi1x, phi2x, phi3x, phi1y, phi2y, phi3y;
  SC kernel = 0.0;
  SC kernelDDot = 0.0;
  SC entry = 0.0;

  this->getSpace( )->getMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->getSpace( )->getMesh( )->getNodes( innerElem, y1, y2, y3 );

  SCVT innerArea = this->getSpace( )->getMesh( )->getElemArea( innerElem );
  SCVT outerArea = this->getSpace( )->getMesh( )->getElemArea( outerElem );


  for ( int dksi = 0; dksi < qSize0; dksi++ ) {
    w[0] = lineQuadPoints[qOrder0][dksi];
    weights[0] = lineQuadWeights[qOrder0][dksi];
    for ( int deta3 = 0; deta3 < qSize1; deta3++ ) {
      w[3] = lineQuadPoints[qOrder1][deta3];
      weights[1] = lineQuadWeights[qOrder1][deta3];
      for ( int deta2 = 0; deta2 < qSize2; deta2++ ) {
        w[2] = lineQuadPoints[qOrder2][deta2];
        weights[2] = lineQuadWeights[qOrder2][deta2];
        for ( int deta1 = 0; deta1 < qSize3; deta1++ ) {
          w[1] = lineQuadPoints[qOrder3][deta1];
          weights[3] = lineQuadWeights[qOrder3][deta1];

          switch ( type ) {
              // identical panels
            case 0:
              for ( int simplex = 0; simplex < 6; simplex++ ) {
                this->cube2triIdentical( w, simplex, xRef, yRef, jacobian );
                this->tri2element( x1, x2, x3, xRef, x );
                this->tri2element( y1, y2, y3, yRef, y );
                kernel = evalHypersingularPart1( x, y ) *
                    weights[0] * weights[1] * weights[2] * weights[3] *
                    jacobian;
                kernelDDot = evalHypersingularPart2( x, y ) *
                    weights[0] * weights[1] * weights[2] * weights[3] *
                    jacobian;

                entry += kernel;
                phi1x = 1.0 - xRef[0];
                phi2x = xRef[0] - xRef[1];
                phi3x = xRef[1];
                phi1y = 1.0 - yRef[0];
                phi2y = yRef[0] - yRef[1];
                phi3y = yRef[1];

                matrix.add( 0, 0, kernelDDot * ( normalDot * phi1y * phi1x ) );
                matrix.add( 0, 1, kernelDDot * ( normalDot * phi2y * phi1x ) );
                matrix.add( 0, 2, kernelDDot * ( normalDot * phi3y * phi1x ) );
                matrix.add( 1, 0, kernelDDot * ( normalDot * phi1y * phi2x ) );
                matrix.add( 1, 1, kernelDDot * ( normalDot * phi2y * phi2x ) );
                matrix.add( 1, 2, kernelDDot * ( normalDot * phi3y * phi2x ) );
                matrix.add( 2, 0, kernelDDot * ( normalDot * phi1y * phi3x ) );
                matrix.add( 2, 1, kernelDDot * ( normalDot * phi2y * phi3x ) );
                matrix.add( 2, 2, kernelDDot * ( normalDot * phi3y * phi3x ) );
              }
              break;

              // common edge
            case 1:
              for ( int simplex = 0; simplex < 5; simplex++ ) {
                this->cube2triCommonEdge( w, simplex, xRef, yRef, jacobian );
                // <("): outer element swaps first two nodes, so that the edges agree
                this->tri2element( x1, x2, x3, xRef, outerRot, x, true );
                this->tri2element( y1, y2, y3, yRef, innerRot, y );
                kernel = jacobian * evalHypersingularPart1( x, y ) *
                    weights[0] * weights[1] * weights[2] * weights[3];
                kernelDDot = evalHypersingularPart2( x, y ) *
                    weights[0] * weights[1] * weights[2] * weights[3] *
                    jacobian;

                entry += kernel;
                phi1x = 1.0 - xRef[0];
                phi2x = xRef[0] - xRef[1];
                phi3x = xRef[1];
                phi1y = 1.0 - yRef[0];
                phi2y = yRef[0] - yRef[1];
                phi3y = yRef[1];
                // outer indices swap due to comment <(")
                matrix.add( this->mod3( 1 + outerRot ), this->mod3( innerRot ), kernelDDot * ( normalDot * phi1y * phi1x ) );
                matrix.add( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ), kernelDDot * ( normalDot * phi2y * phi1x ) );
                matrix.add( this->mod3( 1 + outerRot ), this->mod3( 2 + innerRot ), kernelDDot * ( normalDot * phi3y * phi1x ) );
                matrix.add( this->mod3( outerRot ), this->mod3( innerRot ), kernelDDot * ( normalDot * phi1y * phi2x ) );
                matrix.add( this->mod3( outerRot ), this->mod3( 1 + innerRot ), kernelDDot * ( normalDot * phi2y * phi2x ) );
                matrix.add( this->mod3( outerRot ), this->mod3( 2 + innerRot ), kernelDDot * ( normalDot * phi3y * phi2x ) );
                matrix.add( this->mod3( 2 + outerRot ), this->mod3( innerRot ), kernelDDot * ( normalDot * phi1y * phi3x ) );
                matrix.add( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ), kernelDDot * ( normalDot * phi2y * phi3x ) );
                matrix.add( this->mod3( 2 + outerRot ), this->mod3( 2 + innerRot ), kernelDDot * ( normalDot * phi3y * phi3x ) );
              }
              break;

              // common vertex
            case 2:
              for ( int simplex = 0; simplex < 2; simplex++ ) {
                this->cube2triCommonVertex( w, simplex, xRef, yRef, jacobian );
                this->tri2element( x1, x2, x3, xRef, outerRot, x );
                this->tri2element( y1, y2, y3, yRef, innerRot, y );
                kernel = evalHypersingularPart1( x, y ) *
                    weights[0] * weights[1] * weights[2] * weights[3] *
                    jacobian;
                kernelDDot = evalHypersingularPart2( x, y ) *
                    weights[0] * weights[1] * weights[2] * weights[3] *
                    jacobian;
                entry += kernel;
                phi1x = 1.0 - xRef[0];
                phi2x = xRef[0] - xRef[1];
                phi3x = xRef[1];
                phi1y = 1.0 - yRef[0];
                phi2y = yRef[0] - yRef[1];
                phi3y = yRef[1];
                matrix.add( this->mod3( outerRot ), this->mod3( innerRot ), kernelDDot * ( normalDot * phi1y * phi1x ) );
                matrix.add( this->mod3( outerRot ), this->mod3( 1 + innerRot ), kernelDDot * ( normalDot * phi2y * phi1x ) );
                matrix.add( this->mod3( outerRot ), this->mod3( 2 + innerRot ), kernelDDot * ( normalDot * phi3y * phi1x ) );
                matrix.add( this->mod3( 1 + outerRot ), this->mod3( innerRot ), kernelDDot * ( normalDot * phi1y * phi2x ) );
                matrix.add( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ), kernelDDot * ( normalDot * phi2y * phi2x ) );
                matrix.add( this->mod3( 1 + outerRot ), this->mod3( 2 + innerRot ), kernelDDot * ( normalDot * phi3y * phi2x ) );
                matrix.add( this->mod3( 2 + outerRot ), this->mod3( innerRot ), kernelDDot * ( normalDot * phi1y * phi3x ) );
                matrix.add( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ), kernelDDot * ( normalDot * phi2y * phi3x ) );
                matrix.add( this->mod3( 2 + outerRot ), this->mod3( 2 + innerRot ), kernelDDot * ( normalDot * phi3y * phi3x ) );
              }
              break;

              // disjoint triangles  
            case 3:
              this->cube2triDisjoint( w, xRef, yRef, jacobian );
              this->tri2element( x1, x2, x3, xRef, x );
              this->tri2element( y1, y2, y3, yRef, y );
              kernel = evalHypersingularPart1( x, y ) *
                  weights[0] * weights[1] * weights[2] * weights[3] *
                  jacobian;
              kernelDDot = evalHypersingularPart2( x, y ) *
                  weights[0] * weights[1] * weights[2] * weights[3] *
                  jacobian;
              entry += kernel;
              phi1x = 1.0 - xRef[0];
              phi2x = xRef[0] - xRef[1];
              phi3x = xRef[1];
              phi1y = 1.0 - yRef[0];
              phi2y = yRef[0] - yRef[1];
              phi3y = yRef[1];
              matrix.add( 0, 0, kernelDDot * ( normalDot * phi1y * phi1x ) );
              matrix.add( 0, 1, kernelDDot * ( normalDot * phi2y * phi1x ) );
              matrix.add( 0, 2, kernelDDot * ( normalDot * phi3y * phi1x ) );
              matrix.add( 1, 0, kernelDDot * ( normalDot * phi1y * phi2x ) );
              matrix.add( 1, 1, kernelDDot * ( normalDot * phi2y * phi2x ) );
              matrix.add( 1, 2, kernelDDot * ( normalDot * phi3y * phi2x ) );
              matrix.add( 2, 0, kernelDDot * ( normalDot * phi1y * phi3x ) );
              matrix.add( 2, 1, kernelDDot * ( normalDot * phi2y * phi3x ) );
              matrix.add( 2, 2, kernelDDot * ( normalDot * phi3y * phi3x ) );
              break;
          }
        }
      }
    }
  }

  matrix.add( curlMatrix, entry );
  matrix.scale( 4.0 * innerArea * outerArea );
  //std::cout << type << std::endl;
  //matrix.print();
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::computeChebyshevForHypersingular( ) {

  // at first, compute supports of psidot and psiddot depending on current test and basis functions
  SCVT minTest, maxTest, minBasis, maxBasis;
  SCVT minSupp, maxSupp;
  minTest = ( currTestF - 1 ) * dt;
  maxTest = ( currTestF + 1 ) * dt;

  minBasis = ( currBasisF - 1 ) * dt;
  maxBasis = ( currBasisF + 1 ) * dt;

  minSupp = std::max( minTest - maxBasis, 0.0 );
  maxSupp = maxTest - minBasis;

  // roots of Chebyshev polynomial
  SC *zeros = new SC[nChebPoints];
  SC *psi = new SC[nChebPoints];
  SC *psiTilde = new SC[nChebPoints];
  SC t;
  for ( int i = 0; i <= nChebIntervals; i++ ) {
    chebIntervalStarts[i] = minSupp + i * ( maxSupp - minSupp ) / ( (SC) nChebIntervals );
  }
  for ( int i = 0; i < nChebIntervals; i++ ) {
    // compute roots of Chebyshev pol
    chebyshevZeros( chebIntervalStarts[i], chebIntervalStarts[i + 1], zeros );

    // compute values of psi and psitilde in zeros of Cheb. interpolant
    for ( int j = 0; j < nChebPoints; j++ ) {
      t = zeros[j];
      SCVT intStart = std::min( std::max( minTest, minBasis + t ), ( N + nPos ) * dt );
      SCVT intEnd = std::min( maxBasis + t, maxTest );

      psi[j] = integrateHypersingularPart2( intStart, intEnd, t );
      psiTilde[j] = integrateHypersingularPart1( intStart, intEnd, t );

      if ( intEnd <= intStart ) {
        psi[j] = 0.0;
        psiTilde[j] = 0.0;
      }
    }
    chebyshevCoefficients( psi, psiChebCoeffs[i] );
    chebyshevCoefficients( psiTilde, psiTildeChebCoeffs[i] );
  }

  delete [] zeros;
  delete [] psi;
  delete [] psiTilde;
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::chebyshevZeros( SC a, SC b, SC *zeros ) const {
  SC angle;
  int i;
  for ( i = 0; i < nChebPoints; i++ ) {
    angle = (SC) ( 2 * i + 1 ) * M_PI / (SC) ( 2 * nChebPoints );
    zeros[i] = std::cos( angle );
    zeros[i] = 0.5 * ( a + b ) + zeros[i]*0.5 * ( b - a );
  }
}

template<class LO, class SC>
void BEIntegratorWave<LO, SC>::chebyshevCoefficients( SC *fx, SC *coeffs ) const {
  SC angle;
  memset( coeffs, 0, nChebPoints * sizeof (SC ) );
  for ( int i = 0; i < nChebPoints; i++ ) {
    for ( int j = 0; j < nChebPoints; j++ ) {
      angle = (SC) ( i * ( 2 * j + 1 ) ) * M_PI / (SC) ( 2 * nChebPoints );
      coeffs[i] += fx[j] * std::cos( angle );
    }
    coeffs[i] *= 2.0 / (SC) ( nChebPoints );
  }
}

template<class LO, class SC>
SC BEIntegratorWave<LO, SC>::chebyshevInterpolant( SC a, SC b, SC* coeffs, SC t ) const {

  SC di, dip1, dip2, y;

  dip1 = 0.0;
  di = 0.0;
  y = ( 2.0 * t - a - b ) / ( b - a );
  for ( int i = nChebPoints - 1; 1 <= i; i-- ) {
    dip2 = dip1;
    dip1 = di;
    di = 2.0 * y * dip1 - dip2 + coeffs[i];
  }
  return (y * di - dip1 + 0.5 * coeffs[0] );

}


}

#endif

#endif
