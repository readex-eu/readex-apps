/*!
 * @file    BEIntegrator.cpp
 * @author  Michal Merta
 * @author  Jan Zapletal
 * @date    July 10, 2013
 *
 */

#ifdef BEINTEGRATOR_H

namespace bem4i {

template< class LO, class SC, class SpecificIntegrator >
int BEIntegrator< LO, SC, SpecificIntegrator >::mod3arr[5] = { 0, 1, 2, 0, 1 };

template<class LO, class SC, class SpecificIntegrator>
BEIntegrator<LO, SC, SpecificIntegrator>::BEIntegrator( ) {

  this->x1d = nullptr;
  this->x2d = nullptr;
  this->x3d = nullptr;
  this->y1d = nullptr;
  this->y2d = nullptr;
  this->y3d = nullptr;
  this->wxyd = nullptr;
  this->x1dref = nullptr;
  this->x2dref = nullptr;
  this->y1dref = nullptr;
  this->y2dref = nullptr;

  this->x1stref = nullptr;
  this->x2stref = nullptr;
  this->wxst = nullptr;
  this->x1st = nullptr;
  this->x2st = nullptr;
  this->x3st = nullptr;
  this->sx = nullptr;
  this->tx = nullptr;
  this->ux = nullptr;
  this->x1strefCol = nullptr;
  this->x2strefCol = nullptr;
  this->x1stCol = nullptr;
  this->x2stCol = nullptr;
  this->x3stCol = nullptr;
  this->y1strefCol = nullptr;
  this->y2strefCol = nullptr;
  this->y1stCol = nullptr;
  this->y2stCol = nullptr;
  this->y3stCol = nullptr;
  this->wxyst = nullptr;
  this->r11 = nullptr;
  this->r12 = nullptr;
  this->r13 = nullptr;
  this->r21 = nullptr;
  this->r22 = nullptr;
  this->r23 = nullptr;
  this->alpha1 = nullptr;
  this->alpha2 = nullptr;
  this->stau = nullptr;
  this->rot = nullptr;
  this->ref = nullptr;

  this->x1RefIdentical = new SCVT * [ 6 ];
  this->x2RefIdentical = new SCVT * [ 6 ];
  this->y1RefIdentical = new SCVT * [ 6 ];
  this->y2RefIdentical = new SCVT * [ 6 ];
  this->jacobianWeightIdentical = new SCVT * [ 6 ];
  for ( int i = 0; i < 6; ++i ) {
    this->x1RefIdentical[ i ] = nullptr;
    this->x2RefIdentical[ i ] = nullptr;
    this->y1RefIdentical[ i ] = nullptr;
    this->y2RefIdentical[ i ] = nullptr;
    this->jacobianWeightIdentical[ i ] = nullptr;
  }

  this->x1RefCommonEdge = new SCVT * [ 5 ];
  this->x2RefCommonEdge = new SCVT * [ 5 ];
  this->y1RefCommonEdge = new SCVT * [ 5 ];
  this->y2RefCommonEdge = new SCVT * [ 5 ];
  this->jacobianWeightCommonEdge = new SCVT * [ 5 ];
  for ( int i = 0; i < 5; ++i ) {
    this->x1RefCommonEdge[ i ] = nullptr;
    this->x2RefCommonEdge[ i ] = nullptr;
    this->y1RefCommonEdge[ i ] = nullptr;
    this->y2RefCommonEdge[ i ] = nullptr;
    this->jacobianWeightCommonEdge[ i ] = nullptr;
  }

  this->x1RefCommonVertex = new SCVT * [ 2 ];
  this->x2RefCommonVertex = new SCVT * [ 2 ];
  this->y1RefCommonVertex = new SCVT * [ 2 ];
  this->y2RefCommonVertex = new SCVT * [ 2 ];
  this->jacobianWeightCommonVertex = new SCVT * [ 2 ];
  for ( int i = 0; i < 2; ++i ) {
    this->x1RefCommonVertex[ i ] = nullptr;
    this->x2RefCommonVertex[ i ] = nullptr;
    this->y1RefCommonVertex[ i ] = nullptr;
    this->y2RefCommonVertex[ i ] = nullptr;
    this->jacobianWeightCommonVertex[ i ] = nullptr;
  }

  this->x1RefDisjoint = new SCVT * [ 1 ];
  this->x2RefDisjoint = new SCVT * [ 1 ];
  this->y1RefDisjoint = new SCVT * [ 1 ];
  this->y2RefDisjoint = new SCVT * [ 1 ];
  this->jacobianWeightDisjoint = new SCVT * [ 1 ];
  this->x1RefDisjoint[ 0 ] = nullptr;
  this->x2RefDisjoint[ 0 ] = nullptr;
  this->y1RefDisjoint[ 0 ] = nullptr;
  this->y2RefDisjoint[ 0 ] = nullptr;
  this->jacobianWeightDisjoint[ 0 ] = nullptr;

  this->x1ss = nullptr;
  this->x2ss = nullptr;
  this->x3ss = nullptr;
  this->y1ss = nullptr;
  this->y2ss = nullptr;
  this->y3ss = nullptr;

}

template<class LO, class SC, class SpecificIntegrator>
BEIntegrator<LO, SC, SpecificIntegrator>::BEIntegrator(
  const BEIntegrator & orig
  ) {
}

template<class LO, class SC, class SpecificIntegrator>
BEIntegrator<LO, SC, SpecificIntegrator>::BEIntegrator(
  BESpace<LO, SC> * space
  ) {
  this->space = space;
}

template<class LO, class SC, class SpecificIntegrator>
BEIntegrator<LO, SC, SpecificIntegrator>::~BEIntegrator( ) {

  if ( this->x1d ) _mm_free( this->x1d );
  if ( this->x2d ) _mm_free( this->x2d );
  if ( this->x3d ) _mm_free( this->x3d );
  if ( this->y1d ) _mm_free( this->y1d );
  if ( this->y2d ) _mm_free( this->y2d );
  if ( this->y3d ) _mm_free( this->y3d );
  if ( this->x1dref ) _mm_free( this->x1dref );
  if ( this->x2dref ) _mm_free( this->x2dref );
  if ( this->y1dref ) _mm_free( this->y1dref );
  if ( this->y2dref ) _mm_free( this->y2dref );
  if ( this->wxyd ) _mm_free( this->wxyd );

  if ( this->x1stref ) _mm_free( this->x1stref );
  if ( this->x2stref ) _mm_free( this->x2stref );
  if ( this->wxst ) _mm_free( this->wxst );
  if ( this->x1st ) _mm_free( this->x1st );
  if ( this->x2st ) _mm_free( this->x2st );
  if ( this->x3st ) _mm_free( this->x3st );
  if ( this->sx ) _mm_free( this->sx );
  if ( this->tx ) _mm_free( this->tx );
  if ( this->ux ) _mm_free( this->ux );
  if ( this->x1strefCol ) _mm_free( this->x1strefCol );
  if ( this->x2strefCol ) _mm_free( this->x2strefCol );
  if ( this->x1stCol ) _mm_free( this->x1stCol );
  if ( this->x2stCol ) _mm_free( this->x2stCol );
  if ( this->x3stCol ) _mm_free( this->x3stCol );
  if ( this->y1strefCol ) _mm_free( this->y1strefCol );
  if ( this->y2strefCol ) _mm_free( this->y2strefCol );
  if ( this->y1stCol ) _mm_free( this->y1stCol );
  if ( this->y2stCol ) _mm_free( this->y2stCol );
  if ( this->y3stCol ) _mm_free( this->y3stCol );
  if ( this->wxyst ) _mm_free( this->wxyst );
  if ( this->r11 ) _mm_free( this->r11 );
  if ( this->r12 ) _mm_free( this->r12 );
  if ( this->r13 ) _mm_free( this->r13 );
  if ( this->r21 ) _mm_free( this->r21 );
  if ( this->r22 ) _mm_free( this->r22 );
  if ( this->r23 ) _mm_free( this->r23 );
  if ( this->alpha1 ) _mm_free( this->alpha1 );
  if ( this->alpha2 ) _mm_free( this->alpha2 );
  if ( this->stau ) _mm_free( this->stau );
  if ( this->ref ) _mm_free( this->ref );
  if ( this->rot ) _mm_free( this->rot );

  for ( int i = 0; i < 6; ++i ) {
    if ( this->x1RefIdentical[ i ] ) _mm_free( this->x1RefIdentical[ i ] );
    if ( this->x2RefIdentical[ i ] ) _mm_free( this->x2RefIdentical[ i ] );
    if ( this->y1RefIdentical[ i ] ) _mm_free( this->y1RefIdentical[ i ] );
    if ( this->y2RefIdentical[ i ] ) _mm_free( this->y2RefIdentical[ i ] );
    if ( this->jacobianWeightIdentical[ i ] )
      _mm_free( this->jacobianWeightIdentical[ i ] );
  }
  delete [] this->x1RefIdentical;
  delete [] this->x2RefIdentical;
  delete [] this->y1RefIdentical;
  delete [] this->y2RefIdentical;
  delete [] this->jacobianWeightIdentical;

  for ( int i = 0; i < 5; ++i ) {
    if ( this->x1RefCommonEdge[ i ] ) _mm_free( this->x1RefCommonEdge[ i ] );
    if ( this->x2RefCommonEdge[ i ] ) _mm_free( this->x2RefCommonEdge[ i ] );
    if ( this->y1RefCommonEdge[ i ] ) _mm_free( this->y1RefCommonEdge[ i ] );
    if ( this->y2RefCommonEdge[ i ] ) _mm_free( this->y2RefCommonEdge[ i ] );
    if ( this->jacobianWeightCommonEdge[ i ] )
      _mm_free( this->jacobianWeightCommonEdge[ i ] );
  }
  delete [] this->x1RefCommonEdge;
  delete [] this->x2RefCommonEdge;
  delete [] this->y1RefCommonEdge;
  delete [] this->y2RefCommonEdge;
  delete [] this->jacobianWeightCommonEdge;

  for ( int i = 0; i < 2; ++i ) {
    if ( this->x1RefCommonVertex[ i ] )
      _mm_free( this->x1RefCommonVertex[ i ] );
    if ( this->x2RefCommonVertex[ i ] )
      _mm_free( this->x2RefCommonVertex[ i ] );
    if ( this->y1RefCommonVertex[ i ] )
      _mm_free( this->y1RefCommonVertex[ i ] );
    if ( this->y2RefCommonVertex[ i ] )
      _mm_free( this->y2RefCommonVertex[ i ] );
    if ( this->jacobianWeightCommonVertex[ i ] )
      _mm_free( this->jacobianWeightCommonVertex[ i ] );
  }
  delete [] this->x1RefCommonVertex;
  delete [] this->x2RefCommonVertex;
  delete [] this->y1RefCommonVertex;
  delete [] this->y2RefCommonVertex;
  delete [] this->jacobianWeightCommonVertex;

  if ( this->x1RefDisjoint[ 0 ] ) _mm_free( this->x1RefDisjoint[ 0 ] );
  if ( this->x2RefDisjoint[ 0 ] ) _mm_free( this->x2RefDisjoint[ 0 ] );
  if ( this->y1RefDisjoint[ 0 ] ) _mm_free( this->y1RefDisjoint[ 0 ] );
  if ( this->y2RefDisjoint[ 0 ] ) _mm_free( this->y2RefDisjoint[ 0 ] );
  if ( this->jacobianWeightDisjoint[ 0 ] )
    _mm_free( this->jacobianWeightDisjoint[ 0 ] );
  delete [] this->x1RefDisjoint;
  delete [] this->x2RefDisjoint;
  delete [] this->y1RefDisjoint;
  delete [] this->y2RefDisjoint;
  delete [] this->jacobianWeightDisjoint;

  if ( this->x1ss ) _mm_free( this->x1ss );
  if ( this->x2ss ) _mm_free( this->x2ss );
  if ( this->x3ss ) _mm_free( this->x3ss );
  if ( this->y1ss ) _mm_free( this->y1ss );
  if ( this->y2ss ) _mm_free( this->y2ss );
  if ( this->y3ss ) _mm_free( this->y3ss );
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::getLinValues(
  int order,
  SCVT * linValues
  ) const {

  for ( int i = 0; i < quadSizes[ order ]; i++ ) {
    linValues[ i ] = 1.0 - quadPoints[ order ][ 2 * i ] -
      quadPoints[ order ][ 2 * i + 1 ];
    linValues[ i + quadSizes[ order ] ] = quadPoints[ order ][ 2 * i ];
    linValues[ i + 2 * quadSizes[ order ] ] = quadPoints[ order ][ 2 * i + 1 ];
  }
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::getQuadratureNodes(
  const SCVT * x1,
  const SCVT * x2,
  const SCVT * x3,
  int quadratureOrder,
  SCVT * nodes
  ) const {

  int numPoints = quadSizes[ quadratureOrder ];
  for ( int i = 0; i < numPoints; i++ ) {
    nodes[ i * 3 ] = x1[ 0 ] + ( x2[ 0 ] - x1[ 0 ] ) *
      quadPoints[ quadratureOrder ][ i * 2 ] + ( x3[ 0 ] - x1[ 0 ] ) *
      quadPoints[ quadratureOrder ][ i * 2 + 1 ];
    nodes[ i * 3 + 1 ] = x1[ 1 ] + ( x2[ 1 ] - x1[ 1 ] ) *
      quadPoints[ quadratureOrder ][ i * 2 ] + ( x3[ 1 ] - x1[ 1 ] ) *
      quadPoints[ quadratureOrder ][ i * 2 + 1 ];
    nodes[ i * 3 + 2 ] = x1[ 2 ] + ( x2[ 2 ] - x1[ 2 ] ) *
      quadPoints[ quadratureOrder ][ i * 2 ] + ( x3[ 2 ] - x1[ 2 ] ) *
      quadPoints[ quadratureOrder ][ i * 2 + 1 ];
  }
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::getLocalCoordinates(
  const SCVT * x1,
  const SCVT * x2,
  const SCVT * x3,
  SCVT * r1,
  SCVT * r2,
  SCVT * n,
  SCVT & stau,
  SCVT & alpha1,
  SCVT & alpha2
  ) const {

  SCVT tk, t0;

  r2[0] = x3[0] - x2[0];
  r2[1] = x3[1] - x2[1];
  r2[2] = x3[2] - x2[2];

  tk = std::sqrt( DOT3( r2, r2 ) );

  r2[0] /= tk;
  r2[1] /= tk;
  r2[2] /= tk;

  r1[0] = x2[0] - x1[0];
  r1[1] = x2[1] - x1[1];
  r1[2] = x2[2] - x1[2];

  t0 = -DOT3( r1, r2 );

  r1[0] += t0 * r2[0];
  r1[1] += t0 * r2[1];
  r1[2] += t0 * r2[2];

  stau = std::sqrt( DOT3( r1, r1 ) );

  r1[0] /= stau;
  r1[1] /= stau;
  r1[2] /= stau;

  n[0] = r1[1] * r2[2] - r1[2] * r2[1];
  n[1] = -r1[0] * r2[2] + r1[2] * r2[0];
  n[2] = r1[0] * r2[1] - r1[1] * r2[0];

  alpha1 = -t0 / stau;
  alpha2 = ( tk - t0 ) / stau;
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::getLocalCoordinates(
  const SCVT * y1,
  const SCVT * y2,
  const SCVT * y3,
  SCVT * r1,
  SCVT * r2
  ) const {

  SCVT tk, t0, stau;

  r2[ 0 ] = y3[ 0 ] - y2[ 0 ];
  r2[ 1 ] = y3[ 1 ] - y2[ 1 ];
  r2[ 2 ] = y3[ 2 ] - y2[ 2 ];

  tk = std::sqrt( DOT3( r2, r2 ) );

  r2[ 0 ] /= tk;
  r2[ 1 ] /= tk;
  r2[ 2 ] /= tk;

  r1[ 0 ] = y2[ 0 ] - y1[ 0 ];
  r1[ 1 ] = y2[ 1 ] - y1[ 1 ];
  r1[ 2 ] = y2[ 2 ] - y1[ 2 ];

  t0 = -DOT3( r1, r2 );

  r1[ 0 ] += t0 * r2[ 0 ];
  r1[ 1 ] += t0 * r2[ 1 ];
  r1[ 2 ] += t0 * r2[ 2 ];

  stau = std::sqrt( DOT3( r1, r1 ) );

  r1[ 0 ] /= stau;
  r1[ 1 ] /= stau;
  r1[ 2 ] /= stau;
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::finalizeLocalCoordinates(
  const SCVT * y1,
  const SCVT * y2,
  const SCVT * y3
  ) const {

  SCVT tk, t0;
  int i;
  SCVT yrot11, yrot12, yrot13;
  SCVT yrot21, yrot22, yrot23;
  SCVT yrot31, yrot32, yrot33;

  SCVT * r11 = this->r11;
  SCVT * r12 = this->r12;
  SCVT * r13 = this->r13;
  SCVT * r21 = this->r21;
  SCVT * r22 = this->r22;
  SCVT * r23 = this->r23;
  SCVT * stau = this->stau;
  SCVT * alpha1 = this->alpha1;
  SCVT * alpha2 = this->alpha2;
  int * rot = this->rot;
  const int width = DATA_WIDTH;
  const int align = DATA_ALIGN;

#pragma omp simd \
private( yrot11, yrot12, yrot13 ) \
private( yrot21, yrot22, yrot23 ) \
private( yrot31, yrot32, yrot33 ) \
private( tk, t0 ) \
linear( i : 1 ) \
aligned( r11, r12, r13, r21, r22, r23, stau, alpha1, alpha2, rot : align ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXpad; ++i ) {
    if ( rot[ i ] == 1 ) {
      yrot11 = y2[ 0 ];
      yrot12 = y2[ 1 ];
      yrot13 = y2[ 2 ];
      yrot21 = y3[ 0 ];
      yrot22 = y3[ 1 ];
      yrot23 = y3[ 2 ];
      yrot31 = y1[ 0 ];
      yrot32 = y1[ 1 ];
      yrot33 = y1[ 2 ];
    } else if ( rot[ i ] == 2 ) {
      yrot11 = y3[ 0 ];
      yrot12 = y3[ 1 ];
      yrot13 = y3[ 2 ];
      yrot21 = y1[ 0 ];
      yrot22 = y1[ 1 ];
      yrot23 = y1[ 2 ];
      yrot31 = y2[ 0 ];
      yrot32 = y2[ 1 ];
      yrot33 = y2[ 2 ];
    } else { // rot[ i ] = 0
      yrot11 = y1[ 0 ];
      yrot12 = y1[ 1 ];
      yrot13 = y1[ 2 ];
      yrot21 = y2[ 0 ];
      yrot22 = y2[ 1 ];
      yrot23 = y2[ 2 ];
      yrot31 = y3[ 0 ];
      yrot32 = y3[ 1 ];
      yrot33 = y3[ 2 ];
    }

    r21[ i ] = yrot31 - yrot21;
    r22[ i ] = yrot32 - yrot22;
    r23[ i ] = yrot33 - yrot23;

    tk = std::sqrt( r21[ i ] * r21[ i ] + r22[ i ] * r22[ i ]
      + r23[ i ] * r23[ i ] );

    r21[ i ] /= tk;
    r22[ i ] /= tk;
    r23[ i ] /= tk;

    r11[ i ] = yrot21 - yrot11;
    r12[ i ] = yrot22 - yrot12;
    r13[ i ] = yrot23 - yrot13;

    t0 = -( r11[ i ] * r21[ i ] + r12[ i ] * r22[ i ]
      + r13[ i ] * r23[ i ] );

    r11[ i ] += t0 * r21[ i ];
    r12[ i ] += t0 * r22[ i ];
    r13[ i ] += t0 * r23[ i ];

    stau[ i ] = std::sqrt( r11[ i ] * r11[ i ] + r12[ i ] * r12[ i ]
      + r13[ i ] * r13[ i ] );

    r11[ i ] /= stau[ i ];
    r12[ i ] /= stau[ i ];
    r13[ i ] /= stau[ i ];

    alpha1[ i ] = -t0 / stau[ i ];
    alpha2[ i ] = ( tk - t0 ) / stau[ i ];
  }

}

template<class LO, class SC, class SpecificIntegrator>
bool BEIntegrator<LO, SC, SpecificIntegrator>::checkLocalCoordinates(
  const SCVT * r1,
  const SCVT * r2,
  const SCVT * n
  ) const {

  SCVT aux1, aux2;
  int i;
  bool ret = true;
  
  int * rot = this->rot;
  SCVT * ref = this->ref;
  SCVT * sx = this->sx;
  SCVT * tx = this->tx;
  SCVT * ux = this->ux;
  const int width = DATA_WIDTH;
  const int align = DATA_ALIGN;

//#pragma nounroll
#pragma omp simd \
private( aux1, aux2 ) \
aligned( sx, tx, ux, rot, ref : align ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXpad; ++i ) {
    aux1 = tx[ i ] * tx[ i ];
    aux2 = sx[ i ] * sx[ i ] + ux[ i ] * ux[ i ];

    ref[ i ] = aux1 / aux2;
    rot[ i ] = 0;
  }
  
//#pragma nounroll
  for( i = 0; i < this->qSizeXpad; ++i ){
    if( ref[ i ] < this->coordCoeff ){
     ret = false;
     //break;
    }
  }
  
  return ret;
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::chooseAndUpdateLocalCoordinates(
  const SCVT * y1,
  const SCVT * y2,
  const SCVT * y3,
  const SCVT * n
  ) const {

  SCVT aux1, aux2;
  SCVT r1[ 3 ], r2[ 3 ];

  int * rot = this->rot;
  SCVT * ref = this->ref;
  SCVT * sx = this->sx;
  SCVT * tx = this->tx;
  SCVT * ux = this->ux;
  const int width = DATA_WIDTH;
  const int align = DATA_ALIGN;

  int i;

  this->getLocalCoordinates( y1, y2, y3, r1, r2 );
  this->updateSteinbachLocalQuadratureNodes( y1, r1, r2, n );

#pragma omp simd \
private( aux1, aux2 ) \
aligned( sx, tx, ux, rot, ref : align ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXpad; ++i ) {
    aux1 = tx[ i ] * tx[ i ];
    aux2 = sx[ i ] * sx[ i ] + ux[ i ] * ux[ i ];

    ref[ i ] = aux1 / aux2;
    rot[ i ] = 0;
  }

  this->getLocalCoordinates( y2, y3, y1, r1, r2 );
  this->updateSteinbachLocalQuadratureNodes( y2, r1, r2, n );

#pragma omp simd \
private( aux1, aux2 ) \
aligned( sx, tx, ux, rot, ref : align ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXpad; ++i ) {
    aux1 = tx[ i ] * tx[ i ];
    aux2 = sx[ i ] * sx[ i ] + ux[ i ] * ux[ i ];

    if ( aux1 < ref[ i ] * aux2 ) {
      ref[ i ] = aux1 / aux2;
      rot[ i ] = 1;
    }
  }

  this->getLocalCoordinates( y3, y1, y2, r1, r2 );
  this->updateSteinbachLocalQuadratureNodes( y3, r1, r2, n );

#pragma omp simd \
private( aux1, aux2 ) \
aligned( sx, tx, ux, rot, ref : align ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXpad; ++i ) {
    aux1 = tx[ i ] * tx[ i ];
    aux2 = sx[ i ] * sx[ i ] + ux[ i ] * ux[ i ];

    if ( aux1 < ref[ i ] * aux2 ) {
      rot[ i ] = 2;
    }
  }

  // testing
  //  for ( i = 0; i < this->qSizeXpad; ++i )
  //    rot[ i ] = 2;

  this->finalizeLocalCoordinates( y1, y2, y3 );
  this->finalizeSteinbachLocalQuadratureNodes( y1, y2, y3, n );

}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::
chooseAndUpdateLocalCoordinatesFromRot2(
  const SCVT * y1,
  const SCVT * y2,
  const SCVT * y3,
  const SCVT * n
  ) const {

  SCVT aux1, aux2;
  SCVT r1[ 3 ], r2[ 3 ];

  int * rot = this->rot;
  SCVT * ref = this->ref;
  SCVT * sx = this->sx;
  SCVT * tx = this->tx;
  SCVT * ux = this->ux;
  const int width = DATA_WIDTH;
  const int align = DATA_ALIGN;

  int i;

  this->getLocalCoordinates( y2, y3, y1, r1, r2 );
  this->updateSteinbachLocalQuadratureNodes( y2, r1, r2, n );

#pragma omp simd \
private( aux1, aux2 ) \
aligned( sx, tx, ux, rot, ref : align ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXpad; ++i ) {
    aux1 = tx[ i ] * tx[ i ];
    aux2 = sx[ i ] * sx[ i ] + ux[ i ] * ux[ i ];

    if ( aux1 < ref[ i ] * aux2 ) {
      ref[ i ] = aux1 / aux2;
      rot[ i ] = 1;
    }
  }

  this->getLocalCoordinates( y3, y1, y2, r1, r2 );
  this->updateSteinbachLocalQuadratureNodes( y3, r1, r2, n );

#pragma omp simd \
private( aux1, aux2 ) \
aligned( sx, tx, ux, rot, ref : align ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXpad; ++i ) {
    aux1 = tx[ i ] * tx[ i ];
    aux2 = sx[ i ] * sx[ i ] + ux[ i ] * ux[ i ];

    if ( aux1 < ref[ i ] * aux2 ) {
      rot[ i ] = 2;
    }
  }

  this->finalizeLocalCoordinates( y1, y2, y3 );
  this->finalizeSteinbachLocalQuadratureNodes( y1, y2, y3, n );

}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::getLocalCoordinates(
  LO i,
  SCVT * r1,
  SCVT * r2,
  SCVT * n,
  SCVT & stau,
  SCVT & alpha1,
  SCVT & alpha2,
  int order
  ) const {

  static_cast<const SpecificIntegrator*> ( this )->
    getSpace( )->getRightMesh( )->getNormal( i, n );
  static_cast<const SpecificIntegrator*> ( this )->
    getSpace( )->getRightMesh( )->getR1( i, r1, order );
  static_cast<const SpecificIntegrator*> ( this )->
    getSpace( )->getRightMesh( )->getR2( i, r2, order );
  static_cast<const SpecificIntegrator*> ( this )->
    getSpace( )->getRightMesh( )->getAlpha1( i, alpha1, order );
  static_cast<const SpecificIntegrator*> ( this )->
    getSpace( )->getRightMesh( )->getAlpha2( i, alpha2, order );
  static_cast<const SpecificIntegrator*> ( this )->
    getSpace( )->getRightMesh( )->getSk( i, stau, order );
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::getLocalCoordinates(
  const SCVT * x1,
  const SCVT * x2,
  const SCVT * x3,
  SCVT * r1,
  SCVT * r2,
  SCVT * n,
  SCVT & stau,
  SCVT & alpha1,
  SCVT & alpha2,
  int order
  ) const {

  switch ( order ) {
    case 0:
      getLocalCoordinates( x1, x2, x3, r1, r2, n, stau, alpha1, alpha2 );
      break;
    case 1:
      getLocalCoordinates( x2, x3, x1, r1, r2, n, stau, alpha1, alpha2 );
      break;
    case 2:
      getLocalCoordinates( x3, x1, x2, r1, r2, n, stau, alpha1, alpha2 );
      break;
  }
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::getLocalQuadratureNode(
  const SCVT * x,
  const SCVT * quadratureNode,
  const SCVT * r1,
  const SCVT * r2,
  const SCVT * n,
  SCVT * xLocal
  ) const {

  SCVT tmp[ 3 ];

  tmp[ 0 ] = quadratureNode[ 0 ] - x[ 0 ];
  tmp[ 1 ] = quadratureNode[ 1 ] - x[ 1 ];
  tmp[ 2 ] = quadratureNode[ 2 ] - x[ 2 ];

  xLocal[ 0 ] = DOT3( tmp, r1 );
  xLocal[ 1 ] = DOT3( tmp, r2 );
  xLocal[ 2 ] = DOT3( tmp, n );
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::getLocalQuadratureNode(
  const SCVT * x1,
  const SCVT * x2,
  const SCVT * x3,
  const SCVT * quadratureNode,
  const SCVT * r1,
  const SCVT * r2,
  const SCVT * n,
  SCVT * xLocal,
  int order
  ) const {

  switch ( order ) {
    case 0:
      getLocalQuadratureNode( x1, quadratureNode, r1, r2, n, xLocal );
      break;
    case 1:
      getLocalQuadratureNode( x2, quadratureNode, r1, r2, n, xLocal );
      break;
    case 2:
      getLocalQuadratureNode( x3, quadratureNode, r1, r2, n, xLocal );
      break;
  }
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::getSauterSchwabQuadratureType(
  LO inner,
  LO outer,
  int & type,
  int & rotInner,
  int & rotOuter
  ) const {

  SurfaceMesh3D<LO, SC> *lMesh = static_cast<const SpecificIntegrator*>
    ( this )->getSpace( )->getLeftMesh( );
  SurfaceMesh3D<LO, SC> *rMesh = static_cast<const SpecificIntegrator*>
    ( this )->getSpace( )->getRightMesh( );

  if ( rMesh->getL2gElem( inner ) == lMesh->getL2gElem( outer ) ) {
    // we have identical triangles
    type = 0;
    rotInner = 0;
    rotOuter = 0;
    return;
  }

  LO outerElementGlobal[ 3 ];
  LO innerElementGlobal[ 3 ];

  rMesh->getElementGlobal( inner, innerElementGlobal );
  lMesh->getElementGlobal( outer, outerElementGlobal );

  // check for common edge
  // i ... inner, right
  for ( int i = 0; i < 3; ++i ) {
    // j ... outer, left
    for ( int j = 0; j < 3; ++j ) {
      if ( ( innerElementGlobal[ i ] == outerElementGlobal[ ( j + 1 ) % 3 ] ) &&
        ( innerElementGlobal[ ( i + 1 ) % 3 ] == outerElementGlobal[ j ] ) ) {
        // we have common edge
        type = 1;
        rotOuter = j;
        rotInner = i;
        return;
      }
    }
  }

  // check for common vertex
  // i ... inner, right
  for ( int i = 0; i < 3; ++i ) {
    // j ... outer, left
    for ( int j = 0; j < 3; ++j ) {
      if ( outerElementGlobal[ j ] == innerElementGlobal[ i ] ) {
        // we have common vertex
        type = 2;
        rotOuter = j;
        rotInner = i;
        return;
      }
    }
  }

  // we're clear
  type = 3;
  rotOuter = 0;
  rotInner = 0;

}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::
computeElemMatrix1LayerDisjointP0P0(
  LO outerElem,
  LO innerElem,
  FullMatrix< LO, SC > & elemMatrix
  ) const {

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );

  this->updateDisjointQuadratureNodes( x1, x2, x3, y1, y2, y3 );

  int outerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 0 ] ];
  int innerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 1 ] ];

  LO nIters = outerPoints * innerPoints;
  LO i = 0;
  SC entry = 0.0;
  SC kernel;

  const SpecificIntegrator * thisIntegrator =
    static_cast<const SpecificIntegrator *> ( this );

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
reduction( + : entry ) \
aligned( x1d, x2d, x3d, y1d, y2d, y3d : align ) \
aligned( wxyd : align ) \
private( kernel ) \
simdlen( width )
  for ( i = 0; i < nIters; ++i ) {
    kernel = thisIntegrator->evalSingleLayerKernel( x1d[ i ], x2d[ i ],
      x3d[ i ], y1d[ i ], y2d[ i ], y3d[ i ] );
    entry += wxyd[ i ] * kernel;
  }

  SCVT outerArea = this->space->getLeftMesh( )->getElemArea( outerElem );
  SCVT innerArea = this->space->getRightMesh( )->getElemArea( innerElem );

  entry *= outerArea * innerArea;
  elemMatrix.set( 0, 0, entry );
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::
computeElemMatrix2LayerDisjointP0P0(
  LO outerElem,
  LO innerElem,
  FullMatrix< LO, SC > & elemMatrix
  ) const {

  if ( outerElem == innerElem ) {
    elemMatrix.set( 0, 0, 0.0 );
    return;
  }

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT n[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->space->getRightMesh( )->getNormal( innerElem, n );

  this->updateDisjointQuadratureNodes( x1, x2, x3, y1, y2, y3 );

  int outerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 0 ] ];
  int innerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 1 ] ];

  LO nIters = outerPoints * innerPoints;
  LO i = 0;
  SC entry = 0.0;
  SC kernel;

  const SpecificIntegrator * thisIntegrator =
    static_cast<const SpecificIntegrator *> ( this );

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
reduction( + : entry ) \
aligned( x1d, x2d, x3d, y1d, y2d, y3d : align ) \
aligned( wxyd : align ) \
private( kernel ) \
simdlen( width )
  for ( i = 0; i < nIters; ++i ) {
    kernel = thisIntegrator->evalDoubleLayerKernel( x1d[ i ], x2d[ i ],
      x3d[ i ], y1d[ i ], y2d[ i ], y3d[ i ], n[ 0 ], n[ 1 ],
      n[ 2 ] );
    entry += wxyd[ i ] * kernel;
  }

  SCVT outerArea = this->space->getLeftMesh( )->getElemArea( outerElem );
  SCVT innerArea = this->space->getRightMesh( )->getElemArea( innerElem );

  entry *= outerArea * innerArea;
  elemMatrix.set( 0, 0, entry );
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::
computeElemMatrix2LayerDisjointP0P1(
  LO outerElem,
  LO innerElem,
  FullMatrix< LO, SC > & elemMatrix
  ) const {

  elemMatrix.setAll( 0.0 );
  if ( outerElem == innerElem ) {
    return;
  }

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT n[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->space->getRightMesh( )->getNormal( innerElem, n );

  this->updateDisjointQuadratureNodes( x1, x2, x3, y1, y2, y3 );

  int outerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 0 ] ];
  int innerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 1 ] ];

  LO nIters = outerPoints * innerPoints;
  LO i = 0;
  SC entry1 = 0.0;
  SC entry2 = 0.0;
  SC entry3 = 0.0;
  SC kernel;

  SCVT phi1y, phi2y, phi3y;

  const SpecificIntegrator * thisIntegrator =
    static_cast<const SpecificIntegrator *> ( this );

  SCVT * x1d = this->x1d;
  SCVT * x2d = this->x2d;
  SCVT * x3d = this->x3d;
  SCVT * y1d = this->y1d;
  SCVT * y2d = this->y2d;
  SCVT * y3d = this->y3d;
  SCVT * y1dref = this->y1dref;
  SCVT * y2dref = this->y2dref;
  SCVT * wxyd = this->wxyd;

  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry1, entry2, entry3 ) \
aligned( x1d, x2d, x3d, y1d, y2d, y3d : align ) \
aligned( y1dref, y2dref, wxyd : align ) \
private( kernel, phi1y, phi2y, phi3y ) \
simdlen( width )
  for ( i = 0; i < nIters; ++i ) {

    kernel = wxyd[ i ] *
      thisIntegrator->evalDoubleLayerKernel( x1d[ i ], x2d[ i ], x3d[ i ],
      y1d[ i ], y2d[ i ], y3d[ i ], n[ 0 ], n[ 1 ],
      n[ 2 ] );

    phi1y = (SCVT) 1.0 - y1dref[ i ] - y2dref[ i ];
    phi2y = y1dref[ i ];
    phi3y = y2dref[ i ];

    entry1 += kernel * phi1y;
    entry2 += kernel * phi2y;
    entry3 += kernel * phi3y;
  }

  SCVT areasM = this->space->getLeftMesh( )->getElemArea( outerElem ) *
    this->space->getRightMesh( )->getElemArea( innerElem );
  elemMatrix.set( 0, 0, entry1 * areasM );
  elemMatrix.set( 0, 1, entry2 * areasM );
  elemMatrix.set( 0, 2, entry3 * areasM );
}

template<class LO, class SC, class SpecificIntegrator>
bool BEIntegrator<LO, SC, SpecificIntegrator>::areElementsDisjoint(
  LO outerElem,
  LO innerElem,
  SCVT threshold
  ) const {

  SCVT centrOut[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT centrIn[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r[ 6 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->space->getLeftMesh( )->getCentroid( outerElem, centrOut );
  this->space->getRightMesh( )->getCentroid( innerElem, centrIn );
  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );

  SCVT dist = DIST3SQ( centrOut, centrIn );
  r[0] = DIST3SQ( centrOut, x1 );
  r[3] = DIST3SQ( centrIn, y1 );
  r[1] = DIST3SQ( centrOut, x2 );
  r[4] = DIST3SQ( centrIn, y2 );
  r[2] = DIST3SQ( centrOut, x3 );
  r[5] = DIST3SQ( centrIn, y3 );

  SCVT diam = 0.0;
  for ( int i = 0; i < 6; i++ ) {
    if ( r[i] > diam ) {
      diam = r[i];
    }
  }

  if ( dist / ( 4.0 * diam ) > threshold * threshold ) {
    return true;
  } else {
    return false;
  }
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::
computeElemMatrix1LayerSauterSchwabP1P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  const SpecificIntegrator * thisIntegrator =
    static_cast<const SpecificIntegrator*> ( this );

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

  SC kernel;
  SCVT mult, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y;
  int i;

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

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry11, entry21, entry31 ) \
reduction( + : entry12, entry22, entry32 ) \
reduction( + : entry13, entry23, entry33 ) \
aligned( x1ss, x2ss, x3ss, y1ss, y2ss, y3ss, jacobian : align ) \
aligned( x1ref, x2ref, y1ref, y2ref : align ) \
private( kernel, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y ) \
simdlen( width )
    for ( i = 0; i < totalSize; ++i ) {
      kernel = thisIntegrator->evalSingleLayerKernel( x1ss[ i ], x2ss[ i ],
        x3ss[ i ], y1ss[ i ], y2ss[ i ], y3ss[ i ] );
      kernel *= jacobian[ i ];

      phi1x = (SCVT) 1.0 - x1ref[ i ];
      phi2x = x1ref[ i ] - x2ref[ i ];
      phi3x = x2ref[ i ];
      phi1y = (SCVT) 1.0 - y1ref[ i ];
      phi2y = y1ref[ i ] - y2ref[ i ];
      phi3y = y2ref[ i ];

      entry11 += kernel * phi1x * phi1y;
      entry21 += kernel * phi2x * phi1y;
      entry31 += kernel * phi3x * phi1y;
      entry12 += kernel * phi1x * phi2y;
      entry22 += kernel * phi2x * phi2y;
      entry32 += kernel * phi3x * phi2y;
      entry13 += kernel * phi1x * phi3y;
      entry23 += kernel * phi2x * phi3y;
      entry33 += kernel * phi3x * phi3y;
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

  SCVT areaMult = (SCVT) 4.0 *
    this->getSpace( )->getRightMesh( )->getElemArea( innerElem ) *
    this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );
  matrix.scale( areaMult );
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::
computeElemMatrix1LayerDisjointP1P1(
  LO outerElem,
  LO innerElem,
  FullMatrix< LO, SC > & elemMatrix
  ) const {

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );

  this->updateDisjointQuadratureNodes( x1, x2, x3, y1, y2, y3 );

  int outerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 0 ] ];
  int innerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 1 ] ];

  LO nIters = outerPoints * innerPoints;
  LO i = 0;
  SC entry11 = 0.0;
  SC entry21 = 0.0;
  SC entry31 = 0.0;
  SC entry12 = 0.0;
  SC entry22 = 0.0;
  SC entry32 = 0.0;
  SC entry13 = 0.0;
  SC entry23 = 0.0;
  SC entry33 = 0.0;

  SC kernel;

  SCVT phi1y, phi2y, phi3y, phi1x, phi2x, phi3x;

  const SpecificIntegrator * thisIntegrator =
    static_cast<const SpecificIntegrator *> ( this );

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
reduction( + : entry11, entry21, entry31 ) \
reduction( + : entry12, entry22, entry32 ) \
reduction( + : entry13, entry23, entry33 ) \
aligned( x1d, x2d, x3d, x1dref, x2dref : align ) \
aligned( y1d, y2d, y3d, y1dref, y2dref : align ) \
aligned( wxyd : align ) \
private( kernel, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y ) \
simdlen( width )
  for ( i = 0; i < nIters; ++i ) {

    kernel = wxyd[ i ] *
      thisIntegrator->evalSingleLayerKernel( x1d[ i ], x2d[ i ],
      x3d[ i ], y1d[ i ], y2d[ i ], y3d[ i ] );

    phi1y = (SCVT) 1.0 - y1dref[ i ] - y2dref[ i ];
    phi2y = y1dref[ i ];
    phi3y = y2dref[ i ];
    phi1x = (SCVT) 1.0 - x1dref[ i ] - x2dref[ i ];
    phi2x = x1dref[ i ];
    phi3x = x2dref[ i ];

    entry11 += kernel * phi1x * phi1y;
    entry21 += kernel * phi2x * phi1y;
    entry31 += kernel * phi3x * phi1y;
    entry12 += kernel * phi1x * phi2y;
    entry22 += kernel * phi2x * phi2y;
    entry32 += kernel * phi3x * phi2y;
    entry13 += kernel * phi1x * phi3y;
    entry23 += kernel * phi2x * phi3y;
    entry33 += kernel * phi3x * phi3y;
  }

  SCVT areasM = this->space->getLeftMesh( )->getElemArea( outerElem ) *
    this->space->getRightMesh( )->getElemArea( innerElem );

  elemMatrix.set( 0, 0, entry11 );
  elemMatrix.set( 1, 0, entry21 );
  elemMatrix.set( 2, 0, entry31 );
  elemMatrix.set( 0, 1, entry12 );
  elemMatrix.set( 1, 1, entry22 );
  elemMatrix.set( 2, 1, entry32 );
  elemMatrix.set( 0, 2, entry13 );
  elemMatrix.set( 1, 2, entry23 );
  elemMatrix.set( 2, 2, entry33 );
  elemMatrix.scale( areasM );
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::
computeElemMatrix2LayerSauterSchwabP1P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  const SpecificIntegrator * thisIntegrator =
    static_cast<const SpecificIntegrator*> ( this );

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
  SCVT n[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->getSpace( )->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->getSpace( )->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->getSpace( )->getRightMesh( )->getNormal( innerElem, n );

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

  SC kernel;
  SCVT mult, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y;
  int i;

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

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry11, entry21, entry31 ) \
reduction( + : entry12, entry22, entry32 ) \
reduction( + : entry13, entry23, entry33 ) \
aligned( x1ss, x2ss, x3ss, y1ss, y2ss, y3ss, jacobian : align ) \
aligned( x1ref, x2ref, y1ref, y2ref : align ) \
private( kernel, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y ) \
simdlen( width )
    for ( i = 0; i < totalSize; ++i ) {
      kernel = thisIntegrator->evalDoubleLayerKernel( x1ss[ i ], x2ss[ i ],
        x3ss[ i ], y1ss[ i ], y2ss[ i ], y3ss[ i ], n[ 0 ], n[ 1 ], n[ 2 ] );
      kernel *= jacobian[ i ];

      phi1x = (SCVT) 1.0 - x1ref[ i ];
      phi2x = x1ref[ i ] - x2ref[ i ];
      phi3x = x2ref[ i ];
      phi1y = (SCVT) 1.0 - y1ref[ i ];
      phi2y = y1ref[ i ] - y2ref[ i ];
      phi3y = y2ref[ i ];

      entry11 += kernel * phi1x * phi1y;
      entry21 += kernel * phi2x * phi1y;
      entry31 += kernel * phi3x * phi1y;
      entry12 += kernel * phi1x * phi2y;
      entry22 += kernel * phi2x * phi2y;
      entry32 += kernel * phi3x * phi2y;
      entry13 += kernel * phi1x * phi3y;
      entry23 += kernel * phi2x * phi3y;
      entry33 += kernel * phi3x * phi3y;
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


  SCVT areaMult = (SCVT) 4.0 *
    this->getSpace( )->getRightMesh( )->getElemArea( innerElem ) *
    this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );
  matrix.scale( areaMult );
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::
computeElemMatrix2LayerSauterSchwabP0P0(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  const SpecificIntegrator * thisIntegrator =
    static_cast<const SpecificIntegrator *> ( this );

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
  SCVT n[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->getSpace( )->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->getSpace( )->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->getSpace( )->getRightMesh( )->getNormal( innerElem, n );

  int nSimplex;
  SCVT ** jacobianP;

  this->getReferenceSauterSchwabData( type, nSimplex, jacobianP );

  SCVT * jacobian;

  SCVT * x1ss = this->x1ss;
  SCVT * x2ss = this->x2ss;
  SCVT * x3ss = this->x3ss;
  SCVT * y1ss = this->y1ss;
  SCVT * y2ss = this->y2ss;
  SCVT * y3ss = this->y3ss;

  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

  SC entry = 0.0;
  for ( int simplex = 0; simplex < nSimplex; ++simplex ) {
    this->updateSauterSchwabQuadratureNodes( type, simplex, x1, x2, x3,
      outerRot, y1, y2, y3, innerRot );
    jacobian = jacobianP[ simplex ];

    int i = 0;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry ) \
aligned( x1ss, x2ss, x3ss, y1ss, y2ss, y3ss, jacobian : align ) \
simdlen( width )
    for ( i = 0; i < totalSize; ++i ) {
      entry += thisIntegrator->evalDoubleLayerKernel( x1ss[ i ], x2ss[ i ],
        x3ss[ i ], y1ss[ i ], y2ss[ i ], y3ss[ i ], n[ 0 ], n[ 1 ], n[ 2 ] ) *
        jacobian[ i ];
    }
  }

  SCVT innerArea =
    this->getSpace( )->getRightMesh( )->getElemArea( innerElem );
  SCVT outerArea =
    this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );
  matrix.set( 0, 0, entry * ( (SCVT) 4.0 ) * innerArea * outerArea );
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::
computeElemMatrix2LayerDisjointP1P1(
  LO outerElem,
  LO innerElem,
  FullMatrix< LO, SC > & elemMatrix
  ) const {

  if ( outerElem == innerElem ) {
    elemMatrix.setAll( 0.0 );
    return;
  }

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT n[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->space->getRightMesh( )->getNormal( innerElem, n );

  this->updateDisjointQuadratureNodes( x1, x2, x3, y1, y2, y3 );

  int outerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 0 ] ];
  int innerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 1 ] ];

  LO nIters = outerPoints * innerPoints;
  LO i = 0;
  SC entry11 = 0.0;
  SC entry21 = 0.0;
  SC entry31 = 0.0;
  SC entry12 = 0.0;
  SC entry22 = 0.0;
  SC entry32 = 0.0;
  SC entry13 = 0.0;
  SC entry23 = 0.0;
  SC entry33 = 0.0;

  SC kernel;

  SCVT phi1y, phi2y, phi3y, phi1x, phi2x, phi3x;

  const SpecificIntegrator * thisIntegrator =
    static_cast<const SpecificIntegrator *> ( this );

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
reduction( + : entry11, entry21, entry31 ) \
reduction( + : entry12, entry22, entry32 ) \
reduction( + : entry13, entry23, entry33 ) \
private( kernel, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y ) \
aligned( x1d, x2d, x3d, x1dref, x2dref : align ) \
aligned( y1d, y2d, y3d, y1dref, y2dref : align ) \
aligned( wxyd : align ) \
simdlen( width )
  for ( i = 0; i < nIters; ++i ) {

    kernel = wxyd[ i ] *
      thisIntegrator->evalDoubleLayerKernel( x1d[ i ], x2d[ i ],
      x3d[ i ], y1d[ i ], y2d[ i ], y3d[ i ], n[ 0 ], n[ 1 ],
      n[ 2 ] );

    phi1y = (SCVT) 1.0 - y1dref[ i ] - y2dref[ i ];
    phi2y = y1dref[ i ];
    phi3y = y2dref[ i ];
    phi1x = (SCVT) 1.0 - x1dref[ i ] - x2dref[ i ];
    phi2x = x1dref[ i ];
    phi3x = x2dref[ i ];

    entry11 += kernel * phi1x * phi1y;
    entry21 += kernel * phi2x * phi1y;
    entry31 += kernel * phi3x * phi1y;
    entry12 += kernel * phi1x * phi2y;
    entry22 += kernel * phi2x * phi2y;
    entry32 += kernel * phi3x * phi2y;
    entry13 += kernel * phi1x * phi3y;
    entry23 += kernel * phi2x * phi3y;
    entry33 += kernel * phi3x * phi3y;
  }

  SCVT areasM = this->space->getLeftMesh( )->getElemArea( outerElem ) *
    this->space->getRightMesh( )->getElemArea( innerElem );

  elemMatrix.set( 0, 0, entry11 );
  elemMatrix.set( 1, 0, entry21 );
  elemMatrix.set( 2, 0, entry31 );
  elemMatrix.set( 0, 1, entry12 );
  elemMatrix.set( 1, 1, entry22 );
  elemMatrix.set( 2, 1, entry32 );
  elemMatrix.set( 0, 2, entry13 );
  elemMatrix.set( 1, 2, entry23 );
  elemMatrix.set( 2, 2, entry33 );
  elemMatrix.scale( areasM );
}

// Intel Xeon Phi specific functionality
//#if N_MIC > 0

template<class LO, class SC, class SpecificIntegrator>
bool BEIntegrator<LO, SC, SpecificIntegrator>::areElementsDisjointOnMIC(
  LO outerElem,
  LO innerElem,
  const SCVT * nodes,
  const LO * elems,
  SCVT threshold
  ) {

  SCVT centrOut[3], centrIn[3];
  const SCVT * x1 = &( nodes[ 3 * elems[ 3 * outerElem ] ] );
  const SCVT * x2 = &( nodes[ 3 * elems[ 3 * outerElem + 1 ] ] );
  const SCVT * x3 = &( nodes[ 3 * elems[ 3 * outerElem + 2 ] ] );
  const SCVT * y1 = &( nodes[ 3 * elems[ 3 * innerElem ] ] );
  const SCVT * y2 = &( nodes[ 3 * elems[ 3 * innerElem + 1 ] ] );
  const SCVT * y3 = &( nodes[ 3 * elems[ 3 * innerElem + 2 ] ] );

  SCVT third = (SCVT) 1.0 / (SCVT) 3.0;
  centrOut[ 0 ] = third * ( x1[ 0 ] + x2[ 0 ] + x3[ 0 ] );
  centrOut[ 1 ] = third * ( x1[ 1 ] + x2[ 1 ] + x3[ 1 ] );
  centrOut[ 2 ] = third * ( x1[ 2 ] + x2[ 2 ] + x3[ 2 ] );
  centrIn[ 0 ] = third * ( y1[ 0 ] + y2[ 0 ] + y3[ 0 ] );
  centrIn[ 1 ] = third * ( y1[ 1 ] + y2[ 1 ] + y3[ 1 ] );
  centrIn[ 2 ] = third * ( y1[ 2 ] + y2[ 2 ] + y3[ 2 ] );

  SCVT dist = DIST3SQ( centrOut, centrIn );
  SCVT diam = 0.0;

  SCVT r[ 6 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  r[0] = DIST3SQ( centrOut, x1 );
  r[3] = DIST3SQ( centrIn, y1 );
  r[1] = DIST3SQ( centrOut, x2 );
  r[4] = DIST3SQ( centrIn, y2 );
  r[2] = DIST3SQ( centrOut, x3 );
  r[5] = DIST3SQ( centrIn, y3 );

  for ( int i = 0; i < 6; i++ ) {
    if ( r[i] > diam ) {
      diam = r[i];
    }
  }

  if ( dist / ( (SCVT) 4.0 * diam ) > threshold * threshold ) {
    return true;
  } else {
    return false;
  }
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::cube2triIdentical(
  const SCVT * w,
  int subdomain,
  SCVT * x,
  SCVT * y,
  SCVT & jacobian
  ) const {

  SCVT ksi = w[0];
  SCVT eta1 = w[1];
  SCVT eta2 = w[2];
  SCVT eta3 = w[3];

  switch ( subdomain ) {
    case 0:
      x[0] = ksi;
      x[1] = ksi * ( 1.0 - eta1 + eta1 * eta2 );
      y[0] = ksi * ( 1.0 - eta1 * eta2 * eta3 );
      y[1] = ksi * ( 1.0 - eta1 );
      break;
    case 1:
      x[0] = ksi * ( 1.0 - eta1 * eta2 * eta3 );
      x[1] = ksi * ( 1.0 - eta1 );
      y[0] = ksi;
      y[1] = ksi * ( 1 - eta1 + eta1 * eta2 );
      break;
    case 2:
      x[0] = ksi;
      x[1] = ksi * ( eta1 * ( 1.0 - eta2 + eta2 * eta3 ) );
      y[0] = ksi * ( 1.0 - eta1 * eta2 );
      y[1] = ksi * ( eta1 * ( 1.0 - eta2 ) );
      break;
    case 3:
      x[0] = ksi * ( 1.0 - eta1 * eta2 );
      x[1] = ksi * ( eta1 * ( 1.0 - eta2 ) );
      y[0] = ksi;
      y[1] = ksi * ( eta1 * ( 1.0 - eta2 + eta2 * eta3 ) );
      break;
    case 4:
      x[0] = ksi * ( 1.0 - eta1 * eta2 * eta3 );
      x[1] = ksi * ( eta1 * ( 1.0 - eta2 * eta3 ) );
      y[0] = ksi;
      y[1] = ksi * ( eta1 * ( 1.0 - eta2 ) );
      break;
    case 5:
      x[0] = ksi;
      x[1] = ksi * ( eta1 * ( 1.0 - eta2 ) );
      y[0] = ksi * ( 1.0 - eta1 * eta2 * eta3 );
      y[1] = ksi * ( eta1 * ( 1.0 - eta2 * eta3 ) );
      break;
  }
  jacobian = ksi * ksi * ksi * eta1 * eta1 * eta2;
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::cube2triCommonEdge(
  const SCVT * w,
  int subdomain,
  SCVT * x,
  SCVT * y,
  SCVT & jacobian
  ) const {

  SCVT ksi = w[0];
  SCVT eta1 = w[1];
  SCVT eta2 = w[2];
  SCVT eta3 = w[3];

  jacobian = ksi * ksi * ksi * eta1 * eta1;
  switch ( subdomain ) {
    case 0:
      x[0] = ksi;
      x[1] = ksi * eta1 * eta3;
      y[0] = ksi * ( 1.0 - eta1 * eta2 );
      y[1] = ksi * eta1 * ( 1.0 - eta2 );
      break;
    case 1:
      x[0] = ksi;
      x[1] = ksi * eta1;
      y[0] = ksi * ( 1.0 - eta1 * eta2 * eta3 );
      y[1] = ksi * eta1 * eta2 * ( 1.0 - eta3 );
      jacobian *= eta2;
      break;
    case 2:
      x[0] = ksi * ( 1.0 - eta1 * eta2 );
      x[1] = ksi * eta1 * ( 1.0 - eta2 );
      y[0] = ksi;
      y[1] = ksi * eta1 * eta2 * eta3;
      jacobian *= eta2;
      break;
    case 3:
      x[0] = ksi * ( 1.0 - eta1 * eta2 * eta3 );
      x[1] = ksi * eta1 * eta2 * ( 1.0 - eta3 );
      y[0] = ksi;
      y[1] = ksi * eta1;
      jacobian *= eta2;
      break;
    case 4:
      x[0] = ksi * ( 1.0 - eta1 * eta2 * eta3 );
      x[1] = ksi * eta1 * ( 1.0 - eta2 * eta3 );
      y[0] = ksi;
      y[1] = ksi * eta1 * eta2;
      jacobian *= eta2;
      break;
  }
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::cube2triCommonVertex(
  const SCVT * w,
  int subdomain,
  SCVT * x,
  SCVT * y,
  SCVT & jacobian
  ) const {

  SCVT ksi = w[0];
  SCVT eta1 = w[1];
  SCVT eta2 = w[2];
  SCVT eta3 = w[3];

  switch ( subdomain ) {
    case 0:
      x[0] = ksi;
      x[1] = ksi * eta1;
      y[0] = ksi * eta2;
      y[1] = ksi * eta2 * eta3;
      break;
    case 1:
      x[0] = ksi * eta2;
      x[1] = ksi * eta2 * eta3;
      y[0] = ksi;
      y[1] = ksi * eta1;
      break;
  }
  jacobian = ksi * ksi * ksi * eta2;
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::cube2triDisjoint(
  const SCVT * w,
  SCVT * x,
  SCVT * y,
  SCVT & jacobian
  ) const {

  SCVT ksi = w[0];
  SCVT eta1 = w[1];
  SCVT eta2 = w[2];
  SCVT eta3 = w[3];

  x[0] = ksi;
  x[1] = ksi * eta1;
  y[0] = eta2;
  y[1] = eta2 * eta3;

  jacobian = ksi * eta2;
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::tri2element(
  const SCVT * x1,
  const SCVT * x2,
  const SCVT * x3,
  const SCVT * x,
  SCVT * globalCoord
  ) const {
  globalCoord[ 0 ] =
    x1[ 0 ] + ( x2[ 0 ] - x1[ 0 ] ) * x[ 0 ] + ( x3[ 0 ] - x2[ 0 ] ) * x[ 1 ];
  globalCoord[ 1 ] =
    x1[ 1 ] + ( x2[ 1 ] - x1[ 1 ] ) * x[ 0 ] + ( x3[ 1 ] - x2[ 1 ] ) * x[ 1 ];
  globalCoord[ 2 ] =
    x1[ 2 ] + ( x2[ 2 ] - x1[ 2 ] ) * x[ 0 ] + ( x3[ 2 ] - x2[ 2 ] ) * x[ 1 ];
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::tri2element(
  const SCVT * x1,
  const SCVT * x2,
  const SCVT * x3,
  const SCVT * x,
  int rot,
  SCVT * globalCoord,
  bool swap
  ) const {

  switch ( rot ) {
    case 0:
      if ( swap ) {
        tri2element( x2, x1, x3, x, globalCoord );
      } else {
        tri2element( x1, x2, x3, x, globalCoord );
      }
      break;
    case 1:
      if ( swap ) {
        tri2element( x3, x2, x1, x, globalCoord );
      } else {
        tri2element( x2, x3, x1, x, globalCoord );
      }
      break;
    case 2:
      if ( swap ) {
        tri2element( x1, x3, x2, x, globalCoord );
      } else {
        tri2element( x3, x1, x2, x, globalCoord );
      }
      break;
  }
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::
computeElemMatrix1LayerSauterSchwabP0P0(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  const SpecificIntegrator * thisIntegrator =
    static_cast<const SpecificIntegrator *> ( this );

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

  this->getSpace( )->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->getSpace( )->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );

  int nSimplex;
  SCVT ** jacobianP;

  this->getReferenceSauterSchwabData( type, nSimplex, jacobianP );

  SCVT * jacobian;

  SCVT * x1ss = this->x1ss;
  SCVT * x2ss = this->x2ss;
  SCVT * x3ss = this->x3ss;
  SCVT * y1ss = this->y1ss;
  SCVT * y2ss = this->y2ss;
  SCVT * y3ss = this->y3ss;

  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

  SC entry = 0.0;
  for ( int simplex = 0; simplex < nSimplex; ++simplex ) {
    this->updateSauterSchwabQuadratureNodes( type, simplex, x1, x2, x3,
      outerRot, y1, y2, y3, innerRot );
    jacobian = jacobianP[ simplex ];

    int i = 0;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry ) \
aligned( x1ss, x2ss, x3ss, y1ss, y2ss, y3ss, jacobian : align ) \
simdlen( width )
    for ( i = 0; i < totalSize; ++i ) {
      entry += thisIntegrator->evalSingleLayerKernel( x1ss[ i ], x2ss[ i ],
        x3ss[ i ], y1ss[ i ], y2ss[ i ], y3ss[ i ] ) * jacobian[ i ];
    }
  }

  SCVT innerArea =
    this->getSpace( )->getRightMesh( )->getElemArea( innerElem );
  SCVT outerArea =
    this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );
  matrix.set( 0, 0, entry * ( (SCVT) 4.0 ) * innerArea * outerArea );
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::
computeElemMatrix1LayerSauterSchwabP1P0(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC>& matrix
  ) const {

  matrix.setAll( 0.0 );

  const SpecificIntegrator *thisIntegrator = static_cast<const SpecificIntegrator*> ( this );

  int outerRot;
  int innerRot;
  int type;

  this->getSauterSchwabQuadratureType( innerElem, outerElem, type, innerRot, outerRot );

  int qOrder0 = thisIntegrator->quadratureOrder[ 0 ];
  int qOrder1 = thisIntegrator->quadratureOrder[ 1 ];
  int qOrder2 = thisIntegrator->quadratureOrder[ 2 ];
  int qOrder3 = thisIntegrator->quadratureOrder[ 3 ];

  int qSize0 = lineQuadSizes[ qOrder0 ];
  int qSize1 = lineQuadSizes[ qOrder1 ];
  int qSize2 = lineQuadSizes[ qOrder2 ];
  int qSize3 = lineQuadSizes[ qOrder3 ];

  SCVT xRef[2], yRef[2];
  SCVT x[3], y[3];
  SCVT w[4];
  SCVT x1[3], x2[3], x3[3];
  SCVT y1[3], y2[3], y3[3];
  SCVT n[3];
  SCVT weights[4];
  SCVT jacobian;
  SC kernel = 0.0;

  thisIntegrator->getSpace( )->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  thisIntegrator->getSpace( )->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  thisIntegrator->getSpace( )->getRightMesh( )->getNormal( innerElem, n );

  SCVT innerArea = thisIntegrator->getSpace( )->getRightMesh( )->getElemArea( innerElem );
  SCVT outerArea = thisIntegrator->getSpace( )->getLeftMesh( )->getElemArea( outerElem );

  //if ( outerElem == innerElem ) {
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
                thisIntegrator->cube2triIdentical( w, simplex, xRef, yRef, jacobian );
                thisIntegrator->tri2element( x1, x2, x3, xRef, x );
                thisIntegrator->tri2element( y1, y2, y3, yRef, y );
                kernel = evalSingleLayerKernel( x, y ) *
                  weights[0] * weights[1] * weights[2] * weights[3] *
                  jacobian;
                matrix.add( 0, 0, kernel * ( 1.0 - xRef[0] ) );
                matrix.add( 1, 0, kernel * ( xRef[0] - xRef[1] ) );
                matrix.add( 2, 0, kernel * ( xRef[1] ) );
              }
              break;

              // common edge
            case 1:
              for ( int simplex = 0; simplex < 5; simplex++ ) {
                thisIntegrator->cube2triCommonEdge( w, simplex, xRef, yRef, jacobian );
                thisIntegrator->tri2element( x1, x2, x3, xRef, outerRot, x );
                // inner element swaps first two nodes, so that the edges agree (no change in matrix indices)
                thisIntegrator->tri2element( y1, y2, y3, yRef, innerRot, y, true );
                kernel = jacobian * evalSingleLayerKernel( x, y ) *
                  weights[0] * weights[1] * weights[2] * weights[3];
                matrix.add( this->mod3( innerRot ), 0, kernel * ( 1.0 - xRef[0] ) );
                matrix.add( this->mod3( 1 + innerRot ), 0, kernel * ( xRef[0] - xRef[1] ) );
                matrix.add( this->mod3( 2 + innerRot ), 0, kernel * ( xRef[1] ) );
              }
              break;

              // common vertex
            case 2:
              for ( int simplex = 0; simplex < 2; simplex++ ) {
                thisIntegrator->cube2triCommonVertex( w, simplex, xRef, yRef, jacobian );
                thisIntegrator->tri2element( x1, x2, x3, xRef, outerRot, x );
                thisIntegrator->tri2element( y1, y2, y3, yRef, innerRot, y );
                kernel = evalSingleLayerKernel( x, y ) *
                  weights[0] * weights[1] * weights[2] * weights[3] *
                  jacobian;
                matrix.add( this->mod3( innerRot ), 0, kernel * ( 1.0 - xRef[0] ) );
                matrix.add( this->mod3( 1 + innerRot ), 0, kernel * ( xRef[0] - xRef[1] ) );
                matrix.add( this->mod3( 2 + innerRot ), 0, kernel * ( xRef[1] ) );
              }
              break;

              // disjoint triangles  
            case 3:
              thisIntegrator->cube2triDisjoint( w, xRef, yRef, jacobian );
              thisIntegrator->tri2element( x1, x2, x3, xRef, x );
              thisIntegrator->tri2element( y1, y2, y3, yRef, y );
              kernel = evalSingleLayerKernel( x, y ) *
                weights[0] * weights[1] * weights[2] * weights[3] *
                jacobian;
              matrix.add( 0, 0, kernel * ( 1.0 - xRef[0] ) );
              matrix.add( 1, 0, kernel * ( xRef[0] - xRef[1] ) );
              matrix.add( 2, 0, kernel * ( xRef[1] ) );
              break;
          }
        }
      }
    }
  }

  matrix.scale( 4.0 * innerArea * outerArea );
}

template<class LO, class SC, class SpecificIntegrator>
void BEIntegrator<LO, SC, SpecificIntegrator>::
computeElemMatrix2LayerSauterSchwabP0P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  if ( outerElem == innerElem ) {
    matrix.setAll( 0.0 );
    return;
  }

  const SpecificIntegrator * thisIntegrator =
    static_cast<const SpecificIntegrator *> ( this );

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
  SCVT n[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->getSpace( )->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->getSpace( )->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->getSpace( )->getRightMesh( )->getNormal( innerElem, n );

  int nSimplex;
  SCVT ** jacobianP;
  SCVT ** y1refP;
  SCVT ** y2refP;

  this->getReferenceSauterSchwabData( type, true, nSimplex, jacobianP,
    y1refP, y2refP );

  SCVT * jacobian;
  SCVT * y1ref;
  SCVT * y2ref;

  SC entry1 = 0.0;
  SC entry2 = 0.0;
  SC entry3 = 0.0;
  SC kernel;

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
    y1ref = y1refP[ simplex ];
    y2ref = y2refP[ simplex ];

    int i = 0;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry1, entry2, entry3 ) \
aligned( x1ss, x2ss, x3ss, y1ss, y2ss, y3ss, jacobian : align ) \
aligned( y1ref, y2ref : align ) \
private( kernel ) \
simdlen( width )
    for ( i = 0; i < totalSize; ++i ) {
      kernel = thisIntegrator->evalDoubleLayerKernel( x1ss[ i ], x2ss[ i ],
        x3ss[ i ], y1ss[ i ], y2ss[ i ], y3ss[ i ], n[ 0 ], n[ 1 ], n[ 2 ] ) *
        jacobian[ i ];
      entry1 += kernel * ( (SCVT) 1.0 - y1ref[ i ] );
      entry2 += kernel * ( y1ref[ i ] - y2ref[ i ] );
      entry3 += kernel * y2ref[ i ];
    }
  }

  if ( type == 1 ) {
    matrix.set( 0, this->mod3( 1 + innerRot ), entry1 );
    matrix.set( 0, this->mod3( innerRot ), entry2 );
  } else {
    matrix.set( 0, this->mod3( innerRot ), entry1 );
    matrix.set( 0, this->mod3( 1 + innerRot ), entry2 );
  }
  matrix.set( 0, this->mod3( 2 + innerRot ), entry3 );

  SCVT areaMult = (SCVT) 4.0 *
    this->getSpace( )->getRightMesh( )->getElemArea( innerElem ) *
    this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );
  matrix.scale( areaMult );
}

}

#endif
