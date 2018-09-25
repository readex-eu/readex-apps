/*!
 * @file    BEIntegratorLame.cpp
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    April 29, 2015
 * 
 */

#ifdef BEINTEGRATORLAME_H

namespace bem4i {

template <class LO, class SC>
BEIntegratorLame<LO, SC>::BEIntegratorLame( ) {
}

template <class LO, class SC>
BEIntegratorLame<LO, SC>::BEIntegratorLame(
  const BEIntegratorLame& orig
  ) {
}

template<class LO, class SC>
BEIntegratorLame<LO, SC>::BEIntegratorLame(
  BESpace<LO, SC>* space,
  int* quadratureOrder,
  quadratureType quadrature,
  int* quadratureOrderDisjointElems
  ) {

  this->space = space;
  this->quadratureOrder = quadratureOrder;
  this->quadrature = quadrature;
  this->quadratureOrderDisjointElems = quadratureOrderDisjointElems;

  if ( quadratureOrderDisjointElems ) {
    this->initDisjointQuadratureData( quadratureOrderDisjointElems );
  }

  if ( quadrature == SauterSchwab ) {
    this->initSauterSchwabQuadratureData( quadratureOrder );
  }
}

template<class LO, class SC>
void BEIntegratorLame<LO, SC>::computeElemMatrix1Layer(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  if ( this->space->getTestFunctionType( ) == p0 &&
    this->space->getAnsatzFunctionType( ) == p0 ) {
    // use basic Gaussian quadrature for disjoint elements
    if ( this->quadratureOrderDisjointElems &&
      this->areElementsDisjoint( outerElem, innerElem ) ) {
      this->computeElemMatrix1LayerDisjointP0P0( outerElem, innerElem, matrix );
      return;
    }
    // use something a little more sophisticated for close panels
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
    // use basic Gaussian quadrature for disjoint elements
    if ( this->quadratureOrderDisjointElems &&
      this->areElementsDisjoint( outerElem, innerElem ) ) {
      this->computeElemMatrix1LayerDisjointP1P1( outerElem, innerElem, matrix );
      return;
    }
    // use something a little more sophisticated for close panels
    switch ( this->quadrature ) {
      case SauterSchwab:
        this->computeElemMatrix1LayerSauterSchwabP1P1( outerElem, innerElem,
          matrix );
        break;
      case Steinbach:
        std::cout << "Not implemented!" << std::endl;
        break;
    }
  } else {
    std::cout << "Not implemented!" << std::endl;
  }
}

template<class LO, class SC>
void BEIntegratorLame<LO, SC>::computeElemMatrix1LayerSauterSchwabP0P0(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  int outerRot;
  int innerRot;
  int type;

  this->getSauterSchwabQuadratureType( innerElem, outerElem, type, innerRot,
    outerRot );

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
  //SCVT xMinusY[ 3 ] __attribute__( ( aligned( DATA_ALIGN ) ) );
  SCVT xMinusY0, xMinusY1, xMinusY2;

  this->getSpace( )->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->getSpace( )->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );

  int nSimplex;
  SCVT ** jacobianP;

  this->getReferenceSauterSchwabData( type, nSimplex, jacobianP );

  SCVT * jacobian;

  SC entryLap = 0.0;
  SC entry11 = 0.0, entry22 = 0.0, entry33 = 0.0;
  SC entry12 = 0.0, entry13 = 0.0, entry23 = 0.0;
  SC kernel, kernel3;
  SCVT mult;

  SCVT * x1ss = this->x1ss;
  SCVT * x2ss = this->x2ss;
  SCVT * x3ss = this->x3ss;
  SCVT * y1ss = this->y1ss;
  SCVT * y2ss = this->y2ss;
  SCVT * y3ss = this->y3ss;

  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

  for ( int simplex = 0; simplex < nSimplex; ++simplex ) {
    this->updateSauterSchwabQuadratureNodes( type, simplex, x1, x2, x3,
      outerRot, y1, y2, y3, innerRot );
    jacobian = jacobianP[ simplex ];

    int i = 0;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entryLap, entry11, entry22 ) \
reduction( + : entry33, entry12, entry13, entry23 ) \
aligned( x1ss, x2ss, x3ss, y1ss, y2ss, y3ss, jacobian : align ) \
private( kernel, kernel3 ) \
private( xMinusY0, xMinusY1, xMinusY2 ) \
simdlen( width )
    for ( i = 0; i < totalSize; ++i ) {

      kernel = this->evalSingleLayerKernel( x1ss[ i ], x2ss[ i ], x3ss[ i ],
        y1ss[ i ], y2ss[ i ], y3ss[ i ] );

      kernel3 = kernel * kernel * kernel * jacobian[ i ];
      kernel *= jacobian[ i ];

      entryLap += kernel;

      xMinusY0 = x1ss[ i ] - y1ss[ i ];
      xMinusY1 = x2ss[ i ] - y2ss[ i ];
      xMinusY2 = x3ss[ i ] - y3ss[ i ];

      entry11 += xMinusY0 * xMinusY0 * kernel3;
      entry22 += xMinusY1 * xMinusY1 * kernel3;
      entry33 += xMinusY2 * xMinusY2 * kernel3;
      entry12 += xMinusY0 * xMinusY1 * kernel3;
      entry13 += xMinusY0 * xMinusY2 * kernel3;
      entry23 += xMinusY1 * xMinusY2 * kernel3;
    }
  }

  SCVT areaMult = PI_FACT * (SCVT) 4.0 *
    this->getSpace( )->getRightMesh( )->getElemArea( innerElem ) *
    this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );

  matrix.set( 0, 0, entryLap );
  matrix.set( 0, 1, entry11 );
  matrix.set( 0, 2, entry22 );
  matrix.set( 0, 3, entry33 );
  matrix.set( 0, 4, entry12 );
  matrix.set( 0, 5, entry13 );
  matrix.set( 0, 6, entry23 );
  matrix.scale( areaMult );
}

template<class LO, class SC>
void BEIntegratorLame<LO, SC>::computeElemMatrix1LayerSauterSchwabP1P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  int outerRot;
  int innerRot;
  int type;

  this->getSauterSchwabQuadratureType( innerElem, outerElem, type, innerRot,
    outerRot );

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
  //SCVT xMinusY[ 3 ] __attribute__( ( aligned( DATA_ALIGN ) ) );
  SCVT xMinusY0, xMinusY1, xMinusY2;

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

  SC entryLap11 = 0.0, entryLap12 = 0.0, entryLap13 = 0.0,
    entryLap21 = 0.0, entryLap22 = 0.0, entryLap23 = 0.0,
    entryLap31 = 0.0, entryLap32 = 0.0, entryLap33 = 0.0;
  SC entry11_11 = 0.0, entry11_12 = 0.0, entry11_13 = 0.0,
    entry11_21 = 0.0, entry11_22 = 0.0, entry11_23 = 0.0,
    entry11_31 = 0.0, entry11_32 = 0.0, entry11_33 = 0.0;
  SC entry22_11 = 0.0, entry22_12 = 0.0, entry22_13 = 0.0,
    entry22_21 = 0.0, entry22_22 = 0.0, entry22_23 = 0.0,
    entry22_31 = 0.0, entry22_32 = 0.0, entry22_33 = 0.0;
  SC entry33_11 = 0.0, entry33_12 = 0.0, entry33_13 = 0.0,
    entry33_21 = 0.0, entry33_22 = 0.0, entry33_23 = 0.0,
    entry33_31 = 0.0, entry33_32 = 0.0, entry33_33 = 0.0;
  SC entry12_11 = 0.0, entry12_12 = 0.0, entry12_13 = 0.0,
    entry12_21 = 0.0, entry12_22 = 0.0, entry12_23 = 0.0,
    entry12_31 = 0.0, entry12_32 = 0.0, entry12_33 = 0.0;
  SC entry13_11 = 0.0, entry13_12 = 0.0, entry13_13 = 0.0,
    entry13_21 = 0.0, entry13_22 = 0.0, entry13_23 = 0.0,
    entry13_31 = 0.0, entry13_32 = 0.0, entry13_33 = 0.0;
  SC entry23_11 = 0.0, entry23_12 = 0.0, entry23_13 = 0.0,
    entry23_21 = 0.0, entry23_22 = 0.0, entry23_23 = 0.0,
    entry23_31 = 0.0, entry23_32 = 0.0, entry23_33 = 0.0;
  SC kernel, kernel3;
  SCVT mult, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y;

  SCVT * x1ss = this->x1ss;
  SCVT * x2ss = this->x2ss;
  SCVT * x3ss = this->x3ss;
  SCVT * y1ss = this->y1ss;
  SCVT * y2ss = this->y2ss;
  SCVT * y3ss = this->y3ss;

  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

  for ( int simplex = 0; simplex < nSimplex; ++simplex ) {
    this->updateSauterSchwabQuadratureNodes( type, simplex, x1, x2, x3,
      outerRot, y1, y2, y3, innerRot );
    jacobian = jacobianP[ simplex ];
    x1ref = x1refP[ simplex ];
    x2ref = x2refP[ simplex ];
    y1ref = y1refP[ simplex ];
    y2ref = y2refP[ simplex ];

    int i = 0;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entryLap11, entryLap12, entryLap13 ) \
reduction( + : entryLap21, entryLap22, entryLap23 ) \
reduction( + : entryLap31, entryLap32, entryLap33 ) \
reduction( + : entry11_11, entry11_12, entry11_13 ) \
reduction( + : entry11_21, entry11_22, entry11_23 ) \
reduction( + : entry11_31, entry11_32, entry11_33 ) \
reduction( + : entry22_11, entry22_12, entry22_13 ) \
reduction( + : entry22_21, entry22_22, entry22_23 ) \
reduction( + : entry22_31, entry22_32, entry22_33 ) \
reduction( + : entry33_11, entry33_12, entry33_13 ) \
reduction( + : entry33_21, entry33_22, entry33_23 ) \
reduction( + : entry33_31, entry33_32, entry33_33 ) \
reduction( + : entry12_11, entry12_12, entry12_13 ) \
reduction( + : entry12_21, entry12_22, entry12_23 ) \
reduction( + : entry12_31, entry12_32, entry12_33 ) \
reduction( + : entry13_11, entry13_12, entry13_13 ) \
reduction( + : entry13_21, entry13_22, entry13_23 ) \
reduction( + : entry13_31, entry13_32, entry13_33 ) \
reduction( + : entry23_11, entry23_12, entry23_13 ) \
reduction( + : entry23_21, entry23_22, entry23_23 ) \
reduction( + : entry23_31, entry23_32, entry23_33 ) \
aligned( x1ss, x2ss, x3ss, y1ss, y2ss, y3ss, jacobian : align ) \
aligned( x1ref, x2ref, y1ref, y2ref : align ) \
private( kernel, kernel3, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y ) \
private( xMinusY0, xMinusY1, xMinusY2 ) \
simdlen( width )
    for ( i = 0; i < totalSize; ++i ) {

      kernel = this->evalSingleLayerKernel( x1ss[ i ], x2ss[ i ], x3ss[ i ],
        y1ss[ i ], y2ss[ i ], y3ss[ i ] );

      kernel3 = kernel * kernel * kernel * jacobian[ i ];
      kernel *= jacobian[ i ];

      phi1x = (SCVT) 1.0 - x1ref[ i ];
      phi2x = x1ref[ i ] - x2ref[ i ];
      phi3x = x2ref[ i ];
      phi1y = (SCVT) 1.0 - y1ref[ i ];
      phi2y = y1ref[ i ] - y2ref[ i ];
      phi3y = y2ref[ i ];

      entryLap11 += kernel * phi1x * phi1y;
      entryLap12 += kernel * phi1x * phi2y;
      entryLap13 += kernel * phi1x * phi3y;
      entryLap21 += kernel * phi2x * phi1y;
      entryLap22 += kernel * phi2x * phi2y;
      entryLap23 += kernel * phi2x * phi3y;
      entryLap31 += kernel * phi3x * phi1y;
      entryLap32 += kernel * phi3x * phi2y;
      entryLap33 += kernel * phi3x * phi3y;

      xMinusY0 = x1ss[ i ] - y1ss[ i ];
      xMinusY1 = x2ss[ i ] - y2ss[ i ];
      xMinusY2 = x3ss[ i ] - y3ss[ i ];

      entry11_11 += xMinusY0 * xMinusY0 * kernel3 * phi1x * phi1y;
      entry11_12 += xMinusY0 * xMinusY0 * kernel3 * phi1x * phi2y;
      entry11_13 += xMinusY0 * xMinusY0 * kernel3 * phi1x * phi3y;
      entry11_21 += xMinusY0 * xMinusY0 * kernel3 * phi2x * phi1y;
      entry11_22 += xMinusY0 * xMinusY0 * kernel3 * phi2x * phi2y;
      entry11_23 += xMinusY0 * xMinusY0 * kernel3 * phi2x * phi3y;
      entry11_31 += xMinusY0 * xMinusY0 * kernel3 * phi3x * phi1y;
      entry11_32 += xMinusY0 * xMinusY0 * kernel3 * phi3x * phi2y;
      entry11_33 += xMinusY0 * xMinusY0 * kernel3 * phi3x * phi3y;

      entry22_11 += xMinusY1 * xMinusY1 * kernel3 * phi1x * phi1y;
      entry22_12 += xMinusY1 * xMinusY1 * kernel3 * phi1x * phi2y;
      entry22_13 += xMinusY1 * xMinusY1 * kernel3 * phi1x * phi3y;
      entry22_21 += xMinusY1 * xMinusY1 * kernel3 * phi2x * phi1y;
      entry22_22 += xMinusY1 * xMinusY1 * kernel3 * phi2x * phi2y;
      entry22_23 += xMinusY1 * xMinusY1 * kernel3 * phi2x * phi3y;
      entry22_31 += xMinusY1 * xMinusY1 * kernel3 * phi3x * phi1y;
      entry22_32 += xMinusY1 * xMinusY1 * kernel3 * phi3x * phi2y;
      entry22_33 += xMinusY1 * xMinusY1 * kernel3 * phi3x * phi3y;

      entry33_11 += xMinusY2 * xMinusY2 * kernel3 * phi1x * phi1y;
      entry33_12 += xMinusY2 * xMinusY2 * kernel3 * phi1x * phi2y;
      entry33_13 += xMinusY2 * xMinusY2 * kernel3 * phi1x * phi3y;
      entry33_21 += xMinusY2 * xMinusY2 * kernel3 * phi2x * phi1y;
      entry33_22 += xMinusY2 * xMinusY2 * kernel3 * phi2x * phi2y;
      entry33_23 += xMinusY2 * xMinusY2 * kernel3 * phi2x * phi3y;
      entry33_31 += xMinusY2 * xMinusY2 * kernel3 * phi3x * phi1y;
      entry33_32 += xMinusY2 * xMinusY2 * kernel3 * phi3x * phi2y;
      entry33_33 += xMinusY2 * xMinusY2 * kernel3 * phi3x * phi3y;

      entry12_11 += xMinusY0 * xMinusY1 * kernel3 * phi1x * phi1y;
      entry12_12 += xMinusY0 * xMinusY1 * kernel3 * phi1x * phi2y;
      entry12_13 += xMinusY0 * xMinusY1 * kernel3 * phi1x * phi3y;
      entry12_21 += xMinusY0 * xMinusY1 * kernel3 * phi2x * phi1y;
      entry12_22 += xMinusY0 * xMinusY1 * kernel3 * phi2x * phi2y;
      entry12_23 += xMinusY0 * xMinusY1 * kernel3 * phi2x * phi3y;
      entry12_31 += xMinusY0 * xMinusY1 * kernel3 * phi3x * phi1y;
      entry12_32 += xMinusY0 * xMinusY1 * kernel3 * phi3x * phi2y;
      entry12_33 += xMinusY0 * xMinusY1 * kernel3 * phi3x * phi3y;

      entry13_11 += xMinusY0 * xMinusY2 * kernel3 * phi1x * phi1y;
      entry13_12 += xMinusY0 * xMinusY2 * kernel3 * phi1x * phi2y;
      entry13_13 += xMinusY0 * xMinusY2 * kernel3 * phi1x * phi3y;
      entry13_21 += xMinusY0 * xMinusY2 * kernel3 * phi2x * phi1y;
      entry13_22 += xMinusY0 * xMinusY2 * kernel3 * phi2x * phi2y;
      entry13_23 += xMinusY0 * xMinusY2 * kernel3 * phi2x * phi3y;
      entry13_31 += xMinusY0 * xMinusY2 * kernel3 * phi3x * phi1y;
      entry13_32 += xMinusY0 * xMinusY2 * kernel3 * phi3x * phi2y;
      entry13_33 += xMinusY0 * xMinusY2 * kernel3 * phi3x * phi3y;

      entry23_11 += xMinusY1 * xMinusY2 * kernel3 * phi1x * phi1y;
      entry23_12 += xMinusY1 * xMinusY2 * kernel3 * phi1x * phi2y;
      entry23_13 += xMinusY1 * xMinusY2 * kernel3 * phi1x * phi3y;
      entry23_21 += xMinusY1 * xMinusY2 * kernel3 * phi2x * phi1y;
      entry23_22 += xMinusY1 * xMinusY2 * kernel3 * phi2x * phi2y;
      entry23_23 += xMinusY1 * xMinusY2 * kernel3 * phi2x * phi3y;
      entry23_31 += xMinusY1 * xMinusY2 * kernel3 * phi3x * phi1y;
      entry23_32 += xMinusY1 * xMinusY2 * kernel3 * phi3x * phi2y;
      entry23_33 += xMinusY1 * xMinusY2 * kernel3 * phi3x * phi3y;
    }
  }

  SCVT areaMult = PI_FACT * (SCVT) 4.0 *
    this->getSpace( )->getRightMesh( )->getElemArea( innerElem ) *
    this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );

  int offset = 0;

  if ( type == 1 ) {
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ) + offset,
      entryLap11 );
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ) + offset,
      entryLap12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entryLap21 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ) + offset,
      entryLap22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entryLap31 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ) + offset,
      entryLap32 );
  } else {
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ) + offset,
      entryLap11 );
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ) + offset,
      entryLap12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ) + offset,
      entryLap21 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entryLap22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ) + offset,
      entryLap31 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entryLap32 );
  }
  matrix.set( this->mod3( outerRot ), this->mod3( 2 + innerRot ) + offset,
    entryLap13 );
  matrix.set( this->mod3( 1 + outerRot ), this->mod3( 2 + innerRot ) + offset,
    entryLap23 );
  matrix.set( this->mod3( 2 + outerRot ), this->mod3( 2 + innerRot ) + offset,
    entryLap33 );

  offset += 3;

  if ( type == 1 ) {
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry11_11 );
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ) + offset,
      entry11_12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry11_21 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ) + offset,
      entry11_22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry11_31 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ) + offset,
      entry11_32 );
  } else {
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ) + offset,
      entry11_11 );
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry11_12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ) + offset,
      entry11_21 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry11_22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ) + offset,
      entry11_31 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry11_32 );
  }
  matrix.set( this->mod3( outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry11_13 );
  matrix.set( this->mod3( 1 + outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry11_23 );
  matrix.set( this->mod3( 2 + outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry11_33 );

  offset += 3;

  if ( type == 1 ) {
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry22_11 );
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ) + offset,
      entry22_12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry22_21 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ) + offset,
      entry22_22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry22_31 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ) + offset,
      entry22_32 );
  } else {
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ) + offset,
      entry22_11 );
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry22_12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ) + offset,
      entry22_21 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry22_22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ) + offset,
      entry22_31 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry22_32 );
  }
  matrix.set( this->mod3( outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry22_13 );
  matrix.set( this->mod3( 1 + outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry22_23 );
  matrix.set( this->mod3( 2 + outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry22_33 );

  offset += 3;

  if ( type == 1 ) {
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry33_11 );
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ) + offset,
      entry33_12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry33_21 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ) + offset,
      entry33_22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry33_31 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ) + offset,
      entry33_32 );
  } else {
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ) + offset,
      entry33_11 );
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry33_12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ) + offset,
      entry33_21 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry33_22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ) + offset,
      entry33_31 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry33_32 );
  }
  matrix.set( this->mod3( outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry33_13 );
  matrix.set( this->mod3( 1 + outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry33_23 );
  matrix.set( this->mod3( 2 + outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry33_33 );

  offset += 3;

  if ( type == 1 ) {
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry12_11 );
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ) + offset,
      entry12_12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry12_21 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ) + offset,
      entry12_22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry12_31 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ) + offset,
      entry12_32 );
  } else {
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ) + offset,
      entry12_11 );
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry12_12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ) + offset,
      entry12_21 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry12_22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ) + offset,
      entry12_31 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry12_32 );
  }
  matrix.set( this->mod3( outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry12_13 );
  matrix.set( this->mod3( 1 + outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry12_23 );
  matrix.set( this->mod3( 2 + outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry12_33 );

  offset += 3;

  if ( type == 1 ) {
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry13_11 );
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ) + offset,
      entry13_12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry13_21 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ) + offset,
      entry13_22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry13_31 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ) + offset,
      entry13_32 );
  } else {
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ) + offset,
      entry13_11 );
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry13_12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ) + offset,
      entry13_21 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry13_22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ) + offset,
      entry13_31 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry13_32 );
  }
  matrix.set( this->mod3( outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry13_13 );
  matrix.set( this->mod3( 1 + outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry13_23 );
  matrix.set( this->mod3( 2 + outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry13_33 );

  offset += 3;

  if ( type == 1 ) {
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry23_11 );
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ) + offset,
      entry23_12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry23_21 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ) + offset,
      entry23_22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry23_31 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ) + offset,
      entry23_32 );
  } else {
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ) + offset,
      entry23_11 );
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry23_12 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ) + offset,
      entry23_21 );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry23_22 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ) + offset,
      entry23_31 );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ) + offset,
      entry23_32 );
  }
  matrix.set( this->mod3( outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry23_13 );
  matrix.set( this->mod3( 1 + outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry23_23 );
  matrix.set( this->mod3( 2 + outerRot ), this->mod3( 2 + innerRot ) + offset,
    entry23_33 );

  matrix.scale( areaMult );
}

template<class LO, class SC>
void BEIntegratorLame<LO, SC>::computeElemMatrix1LayerDisjointP0P0(
  LO outerElem,
  LO innerElem,
  FullMatrix< LO, SC > & matrix
  ) const {

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  //SCVT xMinusY[ 3 ] __attribute__( ( aligned( DATA_ALIGN ) ) );
  SCVT xMinusY0, xMinusY1, xMinusY2;

  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );

  this->updateDisjointQuadratureNodes( x1, x2, x3, y1, y2, y3 );

  int outerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 0 ] ];
  int innerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 1 ] ];

  LO nIters = outerPoints * innerPoints;
  LO i = 0;
  SC entryLap = 0.0;
  SC entry11 = 0.0;
  SC entry22 = 0.0;
  SC entry33 = 0.0;
  SC entry12 = 0.0;
  SC entry13 = 0.0;
  SC entry23 = 0.0;
  SC kernel, kernel3;

  SCVT * x1d = this->x1d;
  SCVT * x2d = this->x2d;
  SCVT * x3d = this->x3d;
  SCVT * y1d = this->y1d;
  SCVT * y2d = this->y2d;
  SCVT * y3d = this->y3d;
  SCVT * wxyd = this->wxyd;

  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entryLap, entry11, entry22 ) \
reduction( + : entry33, entry12, entry13, entry23 ) \
aligned( x1d, x2d, x3d, y1d, y2d, y3d : align ) \
aligned( wxyd : align ) \
private( kernel, kernel3 ) \
private( xMinusY0, xMinusY1, xMinusY2 ) \
simdlen( width )
  for ( i = 0; i < nIters; ++i ) {
    kernel = this->evalSingleLayerKernel( x1d[ i ], x2d[ i ], x3d[ i ],
      y1d[ i ], y2d[ i ], y3d[ i ] );

    kernel3 = kernel * kernel * kernel * wxyd[ i ];
    kernel *= wxyd[ i ];

    entryLap += kernel;

    xMinusY0 = x1d[ i ] - y1d[ i ];
    xMinusY1 = x2d[ i ] - y2d[ i ];
    xMinusY2 = x3d[ i ] - y3d[ i ];

    entry11 += xMinusY0 * xMinusY0 * kernel3;
    entry22 += xMinusY1 * xMinusY1 * kernel3;
    entry33 += xMinusY2 * xMinusY2 * kernel3;
    entry12 += xMinusY0 * xMinusY1 * kernel3;
    entry13 += xMinusY0 * xMinusY2 * kernel3;
    entry23 += xMinusY1 * xMinusY2 * kernel3;
  }

  SCVT areaMult = PI_FACT *
    this->space->getLeftMesh( )->getElemArea( outerElem ) *
    this->space->getRightMesh( )->getElemArea( innerElem );

  matrix.set( 0, 0, entryLap );
  matrix.set( 0, 1, entry11 );
  matrix.set( 0, 2, entry22 );
  matrix.set( 0, 3, entry33 );
  matrix.set( 0, 4, entry12 );
  matrix.set( 0, 5, entry13 );
  matrix.set( 0, 6, entry23 );
  matrix.scale( areaMult );
}

template<class LO, class SC>
void BEIntegratorLame<LO, SC>::computeElemMatrix1LayerDisjointP1P1(
  LO outerElem,
  LO innerElem,
  FullMatrix< LO, SC > & matrix
  ) const {

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  //SCVT xMinusY[ 3 ] __attribute__( ( aligned( DATA_ALIGN ) ) );
  SCVT xMinusY0, xMinusY1, xMinusY2;

  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );

  this->updateDisjointQuadratureNodes( x1, x2, x3, y1, y2, y3 );

  int outerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 0 ] ];
  int innerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 1 ] ];

  LO nIters = outerPoints * innerPoints;
  LO i = 0;
  SC entryLap11 = 0.0, entryLap12 = 0.0, entryLap13 = 0.0,
    entryLap21 = 0.0, entryLap22 = 0.0, entryLap23 = 0.0,
    entryLap31 = 0.0, entryLap32 = 0.0, entryLap33 = 0.0;
  SC entry11_11 = 0.0, entry11_12 = 0.0, entry11_13 = 0.0,
    entry11_21 = 0.0, entry11_22 = 0.0, entry11_23 = 0.0,
    entry11_31 = 0.0, entry11_32 = 0.0, entry11_33 = 0.0;
  SC entry22_11 = 0.0, entry22_12 = 0.0, entry22_13 = 0.0,
    entry22_21 = 0.0, entry22_22 = 0.0, entry22_23 = 0.0,
    entry22_31 = 0.0, entry22_32 = 0.0, entry22_33 = 0.0;
  SC entry33_11 = 0.0, entry33_12 = 0.0, entry33_13 = 0.0,
    entry33_21 = 0.0, entry33_22 = 0.0, entry33_23 = 0.0,
    entry33_31 = 0.0, entry33_32 = 0.0, entry33_33 = 0.0;
  SC entry12_11 = 0.0, entry12_12 = 0.0, entry12_13 = 0.0,
    entry12_21 = 0.0, entry12_22 = 0.0, entry12_23 = 0.0,
    entry12_31 = 0.0, entry12_32 = 0.0, entry12_33 = 0.0;
  SC entry13_11 = 0.0, entry13_12 = 0.0, entry13_13 = 0.0,
    entry13_21 = 0.0, entry13_22 = 0.0, entry13_23 = 0.0,
    entry13_31 = 0.0, entry13_32 = 0.0, entry13_33 = 0.0;
  SC entry23_11 = 0.0, entry23_12 = 0.0, entry23_13 = 0.0,
    entry23_21 = 0.0, entry23_22 = 0.0, entry23_23 = 0.0,
    entry23_31 = 0.0, entry23_32 = 0.0, entry23_33 = 0.0;
  SC kernel, kernel3;

  SCVT phi1y, phi2y, phi3y, phi1x, phi2x, phi3x;

  SCVT * x1d = this->x1d;
  SCVT * x2d = this->x2d;
  SCVT * x3d = this->x3d;
  SCVT * y1d = this->y1d;
  SCVT * y2d = this->y2d;
  SCVT * y3d = this->y3d;
  SCVT * y1dref = this->y1dref;
  SCVT * y2dref = this->y2dref;
  SCVT * x1dref = this->x1dref;
  SCVT * x2dref = this->x2dref;
  SCVT * wxyd = this->wxyd;

  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entryLap11, entryLap12, entryLap13 ) \
reduction( + : entryLap21, entryLap22, entryLap23 ) \
reduction( + : entryLap31, entryLap32, entryLap33 ) \
reduction( + : entry11_11, entry11_12, entry11_13 ) \
reduction( + : entry11_21, entry11_22, entry11_23 ) \
reduction( + : entry11_31, entry11_32, entry11_33 ) \
reduction( + : entry22_11, entry22_12, entry22_13 ) \
reduction( + : entry22_21, entry22_22, entry22_23 ) \
reduction( + : entry22_31, entry22_32, entry22_33 ) \
reduction( + : entry33_11, entry33_12, entry33_13 ) \
reduction( + : entry33_21, entry33_22, entry33_23 ) \
reduction( + : entry33_31, entry33_32, entry33_33 ) \
reduction( + : entry12_11, entry12_12, entry12_13 ) \
reduction( + : entry12_21, entry12_22, entry12_23 ) \
reduction( + : entry12_31, entry12_32, entry12_33 ) \
reduction( + : entry13_11, entry13_12, entry13_13 ) \
reduction( + : entry13_21, entry13_22, entry13_23 ) \
reduction( + : entry13_31, entry13_32, entry13_33 ) \
reduction( + : entry23_11, entry23_12, entry23_13 ) \
reduction( + : entry23_21, entry23_22, entry23_23 ) \
reduction( + : entry23_31, entry23_32, entry23_33 ) \
aligned( x1d, x2d, x3d, x1dref, x2dref : align ) \
aligned( y1d, y2d, y3d, y1dref, y2dref : align ) \
aligned( wxyd : align ) \
private( kernel, kernel3, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y ) \
private( xMinusY0, xMinusY1, xMinusY2 ) \
simdlen( width )
  for ( i = 0; i < nIters; ++i ) {
    kernel = this->evalSingleLayerKernel( x1d[ i ], x2d[ i ], x3d[ i ],
      y1d[ i ], y2d[ i ], y3d[ i ] );

    kernel3 = kernel * kernel * kernel * wxyd[ i ];
    kernel *= wxyd[ i ];

    phi1y = (SCVT) 1.0 - y1dref[ i ] - y2dref[ i ];
    phi2y = y1dref[ i ];
    phi3y = y2dref[ i ];
    phi1x = (SCVT) 1.0 - x1dref[ i ] - x2dref[ i ];
    phi2x = x1dref[ i ];
    phi3x = x2dref[ i ];

    entryLap11 += kernel * phi1x * phi1y;
    entryLap12 += kernel * phi1x * phi2y;
    entryLap13 += kernel * phi1x * phi3y;
    entryLap21 += kernel * phi2x * phi1y;
    entryLap22 += kernel * phi2x * phi2y;
    entryLap23 += kernel * phi2x * phi3y;
    entryLap31 += kernel * phi3x * phi1y;
    entryLap32 += kernel * phi3x * phi2y;
    entryLap33 += kernel * phi3x * phi3y;

    xMinusY0 = x1d[ i ] - y1d[ i ];
    xMinusY1 = x2d[ i ] - y2d[ i ];
    xMinusY2 = x3d[ i ] - y3d[ i ];

    entry11_11 += xMinusY0 * xMinusY0 * kernel3 * phi1x * phi1y;
    entry11_12 += xMinusY0 * xMinusY0 * kernel3 * phi1x * phi2y;
    entry11_13 += xMinusY0 * xMinusY0 * kernel3 * phi1x * phi3y;
    entry11_21 += xMinusY0 * xMinusY0 * kernel3 * phi2x * phi1y;
    entry11_22 += xMinusY0 * xMinusY0 * kernel3 * phi2x * phi2y;
    entry11_23 += xMinusY0 * xMinusY0 * kernel3 * phi2x * phi3y;
    entry11_31 += xMinusY0 * xMinusY0 * kernel3 * phi3x * phi1y;
    entry11_32 += xMinusY0 * xMinusY0 * kernel3 * phi3x * phi2y;
    entry11_33 += xMinusY0 * xMinusY0 * kernel3 * phi3x * phi3y;

    entry22_11 += xMinusY1 * xMinusY1 * kernel3 * phi1x * phi1y;
    entry22_12 += xMinusY1 * xMinusY1 * kernel3 * phi1x * phi2y;
    entry22_13 += xMinusY1 * xMinusY1 * kernel3 * phi1x * phi3y;
    entry22_21 += xMinusY1 * xMinusY1 * kernel3 * phi2x * phi1y;
    entry22_22 += xMinusY1 * xMinusY1 * kernel3 * phi2x * phi2y;
    entry22_23 += xMinusY1 * xMinusY1 * kernel3 * phi2x * phi3y;
    entry22_31 += xMinusY1 * xMinusY1 * kernel3 * phi3x * phi1y;
    entry22_32 += xMinusY1 * xMinusY1 * kernel3 * phi3x * phi2y;
    entry22_33 += xMinusY1 * xMinusY1 * kernel3 * phi3x * phi3y;

    entry33_11 += xMinusY2 * xMinusY2 * kernel3 * phi1x * phi1y;
    entry33_12 += xMinusY2 * xMinusY2 * kernel3 * phi1x * phi2y;
    entry33_13 += xMinusY2 * xMinusY2 * kernel3 * phi1x * phi3y;
    entry33_21 += xMinusY2 * xMinusY2 * kernel3 * phi2x * phi1y;
    entry33_22 += xMinusY2 * xMinusY2 * kernel3 * phi2x * phi2y;
    entry33_23 += xMinusY2 * xMinusY2 * kernel3 * phi2x * phi3y;
    entry33_31 += xMinusY2 * xMinusY2 * kernel3 * phi3x * phi1y;
    entry33_32 += xMinusY2 * xMinusY2 * kernel3 * phi3x * phi2y;
    entry33_33 += xMinusY2 * xMinusY2 * kernel3 * phi3x * phi3y;

    entry12_11 += xMinusY0 * xMinusY1 * kernel3 * phi1x * phi1y;
    entry12_12 += xMinusY0 * xMinusY1 * kernel3 * phi1x * phi2y;
    entry12_13 += xMinusY0 * xMinusY1 * kernel3 * phi1x * phi3y;
    entry12_21 += xMinusY0 * xMinusY1 * kernel3 * phi2x * phi1y;
    entry12_22 += xMinusY0 * xMinusY1 * kernel3 * phi2x * phi2y;
    entry12_23 += xMinusY0 * xMinusY1 * kernel3 * phi2x * phi3y;
    entry12_31 += xMinusY0 * xMinusY1 * kernel3 * phi3x * phi1y;
    entry12_32 += xMinusY0 * xMinusY1 * kernel3 * phi3x * phi2y;
    entry12_33 += xMinusY0 * xMinusY1 * kernel3 * phi3x * phi3y;

    entry13_11 += xMinusY0 * xMinusY2 * kernel3 * phi1x * phi1y;
    entry13_12 += xMinusY0 * xMinusY2 * kernel3 * phi1x * phi2y;
    entry13_13 += xMinusY0 * xMinusY2 * kernel3 * phi1x * phi3y;
    entry13_21 += xMinusY0 * xMinusY2 * kernel3 * phi2x * phi1y;
    entry13_22 += xMinusY0 * xMinusY2 * kernel3 * phi2x * phi2y;
    entry13_23 += xMinusY0 * xMinusY2 * kernel3 * phi2x * phi3y;
    entry13_31 += xMinusY0 * xMinusY2 * kernel3 * phi3x * phi1y;
    entry13_32 += xMinusY0 * xMinusY2 * kernel3 * phi3x * phi2y;
    entry13_33 += xMinusY0 * xMinusY2 * kernel3 * phi3x * phi3y;

    entry23_11 += xMinusY1 * xMinusY2 * kernel3 * phi1x * phi1y;
    entry23_12 += xMinusY1 * xMinusY2 * kernel3 * phi1x * phi2y;
    entry23_13 += xMinusY1 * xMinusY2 * kernel3 * phi1x * phi3y;
    entry23_21 += xMinusY1 * xMinusY2 * kernel3 * phi2x * phi1y;
    entry23_22 += xMinusY1 * xMinusY2 * kernel3 * phi2x * phi2y;
    entry23_23 += xMinusY1 * xMinusY2 * kernel3 * phi2x * phi3y;
    entry23_31 += xMinusY1 * xMinusY2 * kernel3 * phi3x * phi1y;
    entry23_32 += xMinusY1 * xMinusY2 * kernel3 * phi3x * phi2y;
    entry23_33 += xMinusY1 * xMinusY2 * kernel3 * phi3x * phi3y;

  }

  SCVT areaMult = PI_FACT *
    this->space->getLeftMesh( )->getElemArea( outerElem ) *
    this->space->getRightMesh( )->getElemArea( innerElem );

  int offset = 0;

  matrix.set( 0, 0 + offset, entryLap11 );
  matrix.set( 0, 1 + offset, entryLap12 );
  matrix.set( 0, 2 + offset, entryLap13 );
  matrix.set( 1, 0 + offset, entryLap21 );
  matrix.set( 1, 1 + offset, entryLap22 );
  matrix.set( 1, 2 + offset, entryLap23 );
  matrix.set( 2, 0 + offset, entryLap31 );
  matrix.set( 2, 1 + offset, entryLap32 );
  matrix.set( 2, 2 + offset, entryLap33 );

  offset += 3;

  matrix.set( 0, 0 + offset, entry11_11 );
  matrix.set( 0, 1 + offset, entry11_12 );
  matrix.set( 0, 2 + offset, entry11_13 );
  matrix.set( 1, 0 + offset, entry11_21 );
  matrix.set( 1, 1 + offset, entry11_22 );
  matrix.set( 1, 2 + offset, entry11_23 );
  matrix.set( 2, 0 + offset, entry11_31 );
  matrix.set( 2, 1 + offset, entry11_32 );
  matrix.set( 2, 2 + offset, entry11_33 );

  offset += 3;

  matrix.set( 0, 0 + offset, entry22_11 );
  matrix.set( 0, 1 + offset, entry22_12 );
  matrix.set( 0, 2 + offset, entry22_13 );
  matrix.set( 1, 0 + offset, entry22_21 );
  matrix.set( 1, 1 + offset, entry22_22 );
  matrix.set( 1, 2 + offset, entry22_23 );
  matrix.set( 2, 0 + offset, entry22_31 );
  matrix.set( 2, 1 + offset, entry22_32 );
  matrix.set( 2, 2 + offset, entry22_33 );

  offset += 3;

  matrix.set( 0, 0 + offset, entry33_11 );
  matrix.set( 0, 1 + offset, entry33_12 );
  matrix.set( 0, 2 + offset, entry33_13 );
  matrix.set( 1, 0 + offset, entry33_21 );
  matrix.set( 1, 1 + offset, entry33_22 );
  matrix.set( 1, 2 + offset, entry33_23 );
  matrix.set( 2, 0 + offset, entry33_31 );
  matrix.set( 2, 1 + offset, entry33_32 );
  matrix.set( 2, 2 + offset, entry33_33 );

  offset += 3;

  matrix.set( 0, 0 + offset, entry12_11 );
  matrix.set( 0, 1 + offset, entry12_12 );
  matrix.set( 0, 2 + offset, entry12_13 );
  matrix.set( 1, 0 + offset, entry12_21 );
  matrix.set( 1, 1 + offset, entry12_22 );
  matrix.set( 1, 2 + offset, entry12_23 );
  matrix.set( 2, 0 + offset, entry12_31 );
  matrix.set( 2, 1 + offset, entry12_32 );
  matrix.set( 2, 2 + offset, entry12_33 );

  offset += 3;

  matrix.set( 0, 0 + offset, entry13_11 );
  matrix.set( 0, 1 + offset, entry13_12 );
  matrix.set( 0, 2 + offset, entry13_13 );
  matrix.set( 1, 0 + offset, entry13_21 );
  matrix.set( 1, 1 + offset, entry13_22 );
  matrix.set( 1, 2 + offset, entry13_23 );
  matrix.set( 2, 0 + offset, entry13_31 );
  matrix.set( 2, 1 + offset, entry13_32 );
  matrix.set( 2, 2 + offset, entry13_33 );

  offset += 3;

  matrix.set( 0, 0 + offset, entry23_11 );
  matrix.set( 0, 1 + offset, entry23_12 );
  matrix.set( 0, 2 + offset, entry23_13 );
  matrix.set( 1, 0 + offset, entry23_21 );
  matrix.set( 1, 1 + offset, entry23_22 );
  matrix.set( 1, 2 + offset, entry23_23 );
  matrix.set( 2, 0 + offset, entry23_31 );
  matrix.set( 2, 1 + offset, entry23_32 );
  matrix.set( 2, 2 + offset, entry23_33 );

  matrix.scale( areaMult );
}

template <class LO, class SC >
BEIntegratorLame<LO, SC>::~BEIntegratorLame( ) {
}

template<class LO, class SC>
void BEIntegratorLame<LO, SC>::computeElemMatrixAllMIC(
  const SCVT * nodes,
  const LO * elems,
  const SCVT * areas,
  const SCVT * normals,
  LO outerElem,
  LO innerElem,
  int qOrderOuter,
  int qOrderInner,
  const SCVT * outerX1ref,
  const SCVT * outerX2ref,
  const SCVT * innerX1ref,
  const SCVT * innerX2ref,
  SCVT * outerX1,
  SCVT * outerX2,
  SCVT * outerX3,
  SCVT * innerX1,
  SCVT * innerX2,
  SCVT * innerX3,
  const SCVT * phi1Values,
  const SCVT * phi2Values,
  const SCVT * phi3Values,
  const SCVT * vOutW,
  const SCVT * vInW,
  SC * elemV11,
  SC * elemV22,
  SC * elemV33,
  SC * elemV12,
  SC * elemV13,
  SC * elemV23,
  SC * elemVLap,
  SC * elemKLap
  ) {

  elemKLap[ 0 ] = elemKLap[ 1 ] = elemKLap[ 2 ] = 0.0;
  bool evalK = true;

  if ( outerElem == innerElem ) {
    evalK = false;
  }

  SCVT x1[3] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[3] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[3] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[3] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[3] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[3] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  memcpy( x1, &nodes[elems[outerElem * 3] * 3], sizeof (SCVT ) * 3 );
  memcpy( x2, &nodes[elems[outerElem * 3 + 1] * 3], sizeof (SCVT ) * 3 );
  memcpy( x3, &nodes[elems[outerElem * 3 + 2] * 3], sizeof (SCVT ) * 3 );
  memcpy( y1, &nodes[elems[innerElem * 3] * 3], sizeof (SCVT ) * 3 );
  memcpy( y2, &nodes[elems[innerElem * 3 + 1] * 3], sizeof (SCVT ) * 3 );
  memcpy( y3, &nodes[elems[innerElem * 3 + 2] * 3], sizeof (SCVT ) * 3 );

  BEIntegrator<LO, SC, BEIntegratorLaplace<LO, SC> >::getQuadratureNodesMIC( x1,
    x2, x3, y1, y2, y3, qOrderOuter, qOrderInner, outerX1ref, outerX2ref,
    innerX1ref, innerX2ref, outerX1, outerX2, outerX3, innerX1, innerX2,
    innerX3 );

  int outerPoints = quadSizes[ qOrderOuter ];
  int innerPoints = quadSizes[ qOrderInner ];

  __assume_aligned( phi1Values, DATA_ALIGN );
  __assume_aligned( phi2Values, DATA_ALIGN );
  __assume_aligned( phi3Values, DATA_ALIGN );
  __assume_aligned( vOutW, DATA_ALIGN );
  __assume_aligned( vInW, DATA_ALIGN );
  __assume_aligned( outerX1, DATA_ALIGN );
  __assume_aligned( outerX2, DATA_ALIGN );
  __assume_aligned( outerX3, DATA_ALIGN );
  __assume_aligned( innerX1, DATA_ALIGN );
  __assume_aligned( innerX2, DATA_ALIGN );
  __assume_aligned( innerX3, DATA_ALIGN );

  LO nIters = outerPoints * innerPoints;
  LO i = 0;
  SC entryK1 = 0.0;
  SC entryK2 = 0.0;
  SC entryK3 = 0.0;
  SC entryVLap = 0.0;
  SC entryV11 = 0.0;
  SC entryV22 = 0.0;
  SC entryV33 = 0.0;
  SC entryV12 = 0.0;
  SC entryV13 = 0.0;
  SC entryV23 = 0.0;

  SCVT areasM;

  //SCVT diffxy[ 3 ] __attribute__( ( aligned( DATA_ALIGN ) ) );
  SCVT diffxy0, diffxy1, diffxy2;

  SCVT norm;
  SCVT dot = 0.0;
  SC kernelKLap = 0.0;
  SC kernelVLap;
  SC kernelCommon;

  SCVT n[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  if ( evalK ) {
    n[ 0 ] = normals[ 3 * innerElem ];
    n[ 1 ] = normals[ 3 * innerElem + 1 ];
    n[ 2 ] = normals[ 3 * innerElem + 2 ];
  }

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entryK1, entryK2, entryK3, entryVLap ) \
reduction( + : entryV11, entryV22, entryV33, entryV12, entryV13, entryV23 ) \
private( norm, kernelVLap, kernelCommon, dot, kernelKLap ) \
private( diffxy0, diffxy1, diffxy2 ) \
simdlen( DATA_WIDTH )
  for ( i = 0; i < nIters; ++i ) {

    diffxy0 = outerX1[ i ] - innerX1[ i ];
    diffxy1 = outerX2[ i ] - innerX2[ i ];
    diffxy2 = outerX3[ i ] - innerX3[ i ];

    norm = std::sqrt( diffxy0 * diffxy0 + diffxy1 * diffxy1
      + diffxy2 * diffxy2 );

    kernelVLap = vOutW[ i ] * vInW[ i ] * ( PI_FACT / norm );
    entryVLap += kernelVLap;

    kernelCommon = kernelVLap / ( norm * norm );

    entryV11 += diffxy0 * diffxy0 * kernelCommon;
    entryV22 += diffxy1 * diffxy1 * kernelCommon;
    entryV33 += diffxy2 * diffxy2 * kernelCommon;
    entryV12 += diffxy0 * diffxy1 * kernelCommon;
    entryV13 += diffxy0 * diffxy2 * kernelCommon;
    entryV23 += diffxy1 * diffxy2 * kernelCommon;

    if ( evalK ) {
      dot = diffxy0 * n[ 0 ] + diffxy1 * n[ 1 ] + diffxy2 * n[ 2 ];
      kernelKLap = kernelCommon * dot;
      entryK1 += kernelKLap * phi1Values[ i ];
      entryK2 += kernelKLap * phi2Values[ i ];
      entryK3 += kernelKLap * phi3Values[ i ];
    }
  }

  areasM = areas[ outerElem ] * areas[ innerElem ];
  elemKLap[ 0 ] = entryK1 * areasM;
  elemKLap[ 1 ] = entryK2 * areasM;
  elemKLap[ 2 ] = entryK3 * areasM;
  *elemVLap = entryVLap * areasM;
  *elemV11 = entryV11 * areasM;
  *elemV22 = entryV22 * areasM;
  *elemV33 = entryV33 * areasM;
  *elemV12 = entryV12 * areasM;
  *elemV13 = entryV13 * areasM;
  *elemV23 = entryV23 * areasM;
}

}

#endif
