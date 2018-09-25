/*!
 * @file    BEIntegratorLaplace.cpp
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    July 10, 2013
 *
 */

#ifdef BEINTEGRATORLAPLACE_H

namespace bem4i {

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::computeElemMatrix1Layer(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC>& matrix
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
        this->computeElemMatrix1LayerP0P0( outerElem, innerElem, matrix );
        break;
    }
  } else if ( ( this->space->getTestFunctionType( ) == p1 ||
    this->space->getTestFunctionType( ) == p1dis ) &&
    ( this->space->getAnsatzFunctionType( ) == p1 ||
    this->space->getAnsatzFunctionType( ) == p1dis ) ) {

    if ( this->quadratureOrderDisjointElems &&
      this->areElementsDisjoint( outerElem, innerElem ) ) {
      this->computeElemMatrix1LayerDisjointP1P1( outerElem, innerElem, matrix );
      return;
    }

    switch ( this->quadrature ) {
      case SauterSchwab:
        this->computeElemMatrix1LayerSauterSchwabP1P1( outerElem, innerElem,
          matrix );
        break;
      case Steinbach:
        this->computeElemMatrix1LayerP1P1( outerElem, innerElem, matrix );
        break;
    }
  } else if ( ( this->space->getTestFunctionType( ) == p1 ||
    this->space->getTestFunctionType( ) == p1dis ) &&
    this->space->getAnsatzFunctionType( ) == p0 ) {
    switch ( this->quadrature ) {
      case SauterSchwab:
        this->computeElemMatrix1LayerSauterSchwabP1P0( outerElem, innerElem,
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
void BEIntegratorLaplace<LO, SC>::computeElemMatrix2Layer(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC>& matrix
  ) const {

  if ( this->space->getTestFunctionType( ) == p0 &&
    ( this->space->getAnsatzFunctionType( ) == p1 ||
    this->space->getAnsatzFunctionType( ) == p1dis ) ) {
    // use basic Gaussian quadrature for disjoint elements
    if ( this->quadratureOrderDisjointElems &&
      this->areElementsDisjoint( outerElem, innerElem ) ) {
      this->computeElemMatrix2LayerDisjointP0P1( outerElem, innerElem, matrix );
      return;
    }
    // use something a little more sophisticated for close panels
    switch ( this->quadrature ) {
      case SauterSchwab:
        this->computeElemMatrix2LayerSauterSchwabP0P1( outerElem, innerElem,
          matrix );
        break;
      case Steinbach:
        this->computeElemMatrix2LayerP0P1( outerElem, innerElem, matrix );
        break;
    }
  } else if ( this->space->getTestFunctionType( ) == p0 &&
    this->space->getAnsatzFunctionType( ) == p0 ) {
    if ( this->quadratureOrderDisjointElems &&
      this->areElementsDisjoint( outerElem, innerElem ) ) {
      this->computeElemMatrix2LayerDisjointP0P0( outerElem, innerElem, matrix );
      return;
    }
    switch ( this->quadrature ) {
      case SauterSchwab:
        this->computeElemMatrix2LayerSauterSchwabP0P0( outerElem, innerElem,
          matrix );
        break;
      case Steinbach:
        this->computeElemMatrix2LayerP0P0( outerElem, innerElem, matrix );
        break;
    }
  } else if ( ( this->space->getTestFunctionType( ) == p1 ||
    this->space->getTestFunctionType( ) == p1dis ) &&
    ( this->space->getAnsatzFunctionType( ) == p1 ||
    this->space->getAnsatzFunctionType( ) == p1dis ) ) {

    if ( this->quadratureOrderDisjointElems &&
      this->areElementsDisjoint( outerElem, innerElem ) ) {
      this->computeElemMatrix2LayerDisjointP1P1( outerElem, innerElem, matrix );
      return;
    }

    switch ( this->quadrature ) {
      case SauterSchwab:
        this->computeElemMatrix2LayerSauterSchwabP1P1( outerElem, innerElem,
          matrix );
        break;
      case Steinbach:
        this->computeElemMatrix2LayerP1P1( outerElem, innerElem, matrix );
        break;
    }
  } else {
    std::cout << "Not implemented!" << std::endl;
  }
}

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::computeElemMatrixHypersingular(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC>& matrix
  ) const {

  if ( ( this->space->getTestFunctionType( ) == p1 ||
    this->space->getTestFunctionType( ) == p1dis ) &&
    ( this->space->getAnsatzFunctionType( ) == p1 ||
    this->space->getAnsatzFunctionType( ) == p1dis ) ) {
    this->computeElemMatrixHypersingularP1P1( outerElem, innerElem, matrix );
  } else {
    std::cout << "Not implemented!" << std::endl;
  }
}

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::representationFormulaP1P1(
  const SCVT * xCoord,
  LO numPoints,
  const Vector<LO, SC> & dir,
  const Vector<LO, SC> & neu,
  bool interior,
  Vector<LO, SC> & values
  ) const {

  LO nElems = this->space->getMesh( )->getNElements( );
  SCVT stau, alpha1, alpha2;
  SCVT x1[3], x2[3], x3[3], r1[3], r2[3], n[3], xLocal[3];
  SCVT centroid[3];
  SCVT r[3] = { 0.0, 0.0, 0.0 };
  SCVT diam = 0.0;
  LO ind[3];
  SC val;
  SC dlkernel, slkernel;
  SC integrand;
  SC onePatch;

  SCVT * quadratureNodes =
    new SCVT[ quadSizes[ this->quadratureOrder[ 1 ] ] * 3 ];

  // points in reference triangle
  double * qPointsIn = quadPoints[ this->quadratureOrder[ 1 ] ];
  SCVT intPx, intPy;

  for ( LO i = 0; i < numPoints; ++i ) {

    val = 0.0;
    const SCVT * xSingle = xCoord + 3 * i;

    for ( LO j = 0; j < nElems; ++j ) {

      this->space->getMesh( )->getCentroid( j, centroid );
      this->space->getMesh( )->getNodes( j, x1, x2, x3 );
      r[0] = DIST3( centroid, x1 );
      r[1] = DIST3( centroid, x2 );
      r[2] = DIST3( centroid, x3 );

      diam = 2.0 * ( *std::max_element( r, r + 3 ) );

      this->space->getMesh( )->getElement( j, ind );
      this->space->getMesh( )->getNormal( j, n );
      this->getQuadratureNodes( x1, x2, x3, this->quadratureOrder[ 1 ],
        quadratureNodes );

      if ( DIST3( centroid, xSingle ) < 2.0 * diam ) {

        for ( int rot = 0; rot < 3; ++rot ) {

          this->getLocalCoordinates( x1, x2, x3, r1, r2, n, stau, alpha1,
            alpha2, rot );
          this->getLocalQuadratureNode( x1, x2, x3, xSingle, r1, r2, n,
            xLocal, rot );

          val += neu.get( ind[ rot ] ) *
            collocation1LayerP1( xLocal[ 0 ], xLocal[ 1 ], xLocal[ 2 ], stau,
            alpha1, alpha2 );

          val -= dir.get( ind[ rot ] ) *
            collocation2LayerP1( xLocal[ 0 ], xLocal[ 1 ], xLocal[ 2 ], stau,
            alpha1, alpha2 );
        }

      } else {

        onePatch = 0.0;

        for ( int k = 0; k < quadSizes[ this->quadratureOrder[ 1 ] ]; ++k ) {

          intPx = qPointsIn[ 2 * k ];
          intPy = qPointsIn[ 2 * k + 1 ];

          slkernel = evalSingleLayerKernel( xSingle,
            ( quadratureNodes + 3 * k ) );
          integrand = neu.get( ind[ 0 ] ) * slkernel *
            ( (SCVT) 1.0 - intPx - intPy );
          integrand += neu.get( ind[ 1 ] ) * slkernel * intPx;
          integrand += neu.get( ind[ 2 ] ) * slkernel * intPy;

          dlkernel = evalDoubleLayerKernel( xSingle,
            ( quadratureNodes + 3 * k ), n );
          integrand -= dir.get( ind[ 0 ] ) * dlkernel *
            ( (SCVT) 1.0 - intPx - intPy );
          integrand -= dir.get( ind[ 1 ] ) * dlkernel * intPx;
          integrand -= dir.get( ind[ 2 ] ) * dlkernel * intPy;

          onePatch += integrand *
            (SCVT) quadWeights[ this->quadratureOrder[ 1 ] ][ k ];
        }

        val += onePatch * this->space->getMesh( )->getElemArea( j );
      }
    }

    values.set( i, val );
  }

  if ( !interior ) {
    values.scale( -1.0 );
  }

  delete [] quadratureNodes;
}

//#if N_MIC > 0

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::computeElemMatrix1LayerP0P0MIC(
  const SCVT * nodes,
  const LO * elems,
  const SCVT * areas,
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
  const SCVT * vOutW,
  const SCVT * vInW,
  SC * elemMatrix
  ) {

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

  BEIntegrator<LO, SC, BEIntegratorLaplace<LO, SC> >::getQuadratureNodesMIC(
    x1, x2, x3, y1, y2, y3, qOrderOuter, qOrderInner, outerX1ref, outerX2ref,
    innerX1ref, innerX2ref, outerX1, outerX2, outerX3, innerX1, innerX2,
    innerX3 );

  int outerPoints = quadSizes[ qOrderOuter ];
  int innerPoints = quadSizes[ qOrderInner ];

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
  SC entry = 0.0;
  SCVT norm;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry ) \
private( norm ) \
simdlen( DATA_WIDTH )
  for ( i = 0; i < nIters; ++i ) {

    norm = std::sqrt(
      ( outerX1[ i ] - innerX1[ i ] ) * ( outerX1[ i ] - innerX1[ i ] ) +
      ( outerX2[ i ] - innerX2[ i ] ) * ( outerX2[ i ] - innerX2[ i ] ) +
      ( outerX3[ i ] - innerX3[ i ] ) * ( outerX3[ i ] - innerX3[ i ] ) );

    entry += vOutW[ i ] * vInW[ i ] * ( PI_FACT / norm );
  }

  entry *= areas[ outerElem ] * areas[ innerElem ];
  *elemMatrix = entry;
}

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::computeElemMatrix2LayerP0P1MIC(
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
  SC * elemMatrix
  ) {

  elemMatrix[ 0 ] = elemMatrix[ 1 ] = elemMatrix[ 2 ] = 0.0;
  if ( outerElem == innerElem ) {
    return;
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
  SC entry1 = 0.0;
  SC entry2 = 0.0;
  SC entry3 = 0.0;
  SCVT areasM;

  SCVT n[3] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  n[ 0 ] = normals[ 3 * innerElem];
  n[ 1 ] = normals[ 3 * innerElem + 1 ];
  n[ 2 ] = normals[ 3 * innerElem + 2 ];
  //SCVT diffxy[ 3 ] __attribute__( ( aligned( DATA_ALIGN ) ) );
  SCVT diffxy0, diffxy1, diffxy2;

  SCVT norm, dot;
  SC kernel;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry1, entry2, entry3 ) \
private( norm, dot, kernel ) \
private( diffxy0, diffxy1, diffxy2 ) \
simdlen( DATA_WIDTH )
  for ( i = 0; i < nIters; ++i ) {

    diffxy0 = outerX1[ i ] - innerX1[ i ];
    diffxy1 = outerX2[ i ] - innerX2[ i ];
    diffxy2 = outerX3[ i ] - innerX3[ i ];

    norm = std::sqrt( diffxy0 * diffxy0 + diffxy1 * diffxy1 +
      diffxy2 * diffxy2 );
    dot = diffxy0 * n[ 0 ] + diffxy1 * n[ 1 ] + diffxy2 * n[ 2 ];

    kernel =
      vOutW[ i ] * vInW[ i ] * ( PI_FACT * dot / ( norm * norm * norm ) );
    entry1 += kernel * phi1Values[ i ];
    entry2 += kernel * phi2Values[ i ];
    entry3 += kernel * phi3Values[ i ];
  }

  areasM = areas[ outerElem ] * areas[ innerElem ];
  elemMatrix[ 0 ] = entry1 * areasM;
  elemMatrix[ 1 ] = entry2 * areasM;
  elemMatrix[ 2 ] = entry3 * areasM;
}

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::computeElemMatrix1And2LayerP0P0P0P1MIC(
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
  SC * elemV,
  SC * elemK
  ) {

  elemK[ 0 ] = elemK[ 1 ] = elemK[ 2 ] = 0.0;
  bool evalK = true;

  if ( outerElem == innerElem ) {
    evalK = false;
  }

  SCVT x1[3] __attribute__ ( ( aligned( 64 ) ) );
  SCVT x2[3] __attribute__ ( ( aligned( 64 ) ) );
  SCVT x3[3] __attribute__ ( ( aligned( 64 ) ) );
  SCVT y1[3] __attribute__ ( ( aligned( 64 ) ) );
  SCVT y2[3] __attribute__ ( ( aligned( 64 ) ) );
  SCVT y3[3] __attribute__ ( ( aligned( 64 ) ) );
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
  SC entryV = 0.0;
  SCVT areasM;

  //SCVT diffxy[ 3 ] __attribute__( ( aligned( 64 ) ) );
  SCVT diffxy0, diffxy1, diffxy2;

  SCVT norm;
  SCVT dot = 0.0;
  SC kernelK = 0.0;
  SC kernelV;

  SCVT n[ 3 ] __attribute__ ( ( aligned( 64 ) ) );
  if ( evalK ) {
    n[ 0 ] = normals[ 3 * innerElem ];
    n[ 1 ] = normals[ 3 * innerElem + 1 ];
    n[ 2 ] = normals[ 3 * innerElem + 2 ];
  }

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entryK1, entryK2, entryK3, entryV ) \
private( norm, kernelV, dot, kernelK ) \
private( diffxy0, diffxy1, diffxy2 ) \
simdlen( DATA_WIDTH )
  for ( i = 0; i < nIters; ++i ) {

    diffxy0 = outerX1[ i ] - innerX1[ i ];
    diffxy1 = outerX2[ i ] - innerX2[ i ];
    diffxy2 = outerX3[ i ] - innerX3[ i ];

    norm = std::sqrt( diffxy0 * diffxy0 + diffxy1 * diffxy1 +
      diffxy2 * diffxy2 );

    kernelV = vOutW[ i ] * vInW[ i ] * ( PI_FACT / norm );
    entryV += kernelV;

    if ( evalK ) {
      dot = diffxy0 * n[ 0 ] + diffxy1 * n[ 1 ] + diffxy2 * n[ 2 ];
      kernelK = ( kernelV * dot ) / ( norm * norm );
      entryK1 += kernelK * phi1Values[ i ];
      entryK2 += kernelK * phi2Values[ i ];
      entryK3 += kernelK * phi3Values[ i ];
    }
  }

  areasM = areas[ outerElem ] * areas[ innerElem ];
  elemK[ 0 ] = entryK1 * areasM;
  elemK[ 1 ] = entryK2 * areasM;
  elemK[ 2 ] = entryK3 * areasM;
  *elemV = entryV * areasM;
}

template<class LO, class SC>
BEIntegratorLaplace<LO, SC>::BEIntegratorLaplace( ) {
}

template<class LO, class SC>
BEIntegratorLaplace<LO, SC>::BEIntegratorLaplace(
  const BEIntegratorLaplace& orig
  ) {

  this->space = orig.space;

}

template<class LO, class SC>
BEIntegratorLaplace<LO, SC>::BEIntegratorLaplace(
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

  if ( quadrature == Steinbach ) {
    this->initSteinbachQuadratureData( quadratureOrder );
  }
}

template<class LO, class SC>
BEIntegratorLaplace<LO, SC>::~BEIntegratorLaplace( ) {
}

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::computeElemMatrix1LayerP0P0(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  SC stau, alpha1, alpha2;

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT n[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  // getting local coordinate system
  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->getLocalCoordinates( y1, y2, y3, r1, r2, n, stau, alpha1, alpha2 );

  // getting outer quadrature points in local coordinate system
  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->updateSteinbachLocalQuadratureNodes( x1, x2, x3, y1, r1, r2, n );

  SCVT * sx = this->sx;
  SCVT * tx = this->tx;
  SCVT * ux = this->ux;
  SCVT * wxst = this->wxst;

  SC entry = 0.0;
  int i;
  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry ) \
aligned( sx, tx, ux, wxst : align ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXpad; ++i ) {
    entry += wxst[ i ] * this->collocation1LayerP0( sx[ i ], tx[ i ],
      ux[ i ], stau, alpha1, alpha2 );
  }

  entry *= this->space->getLeftMesh( )->getElemArea( outerElem );
  matrix.set( 0, 0, entry );
}

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::computeElemMatrix1LayerP1P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  SC stau, alpha1, alpha2;

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT n[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  SCVT areax = this->space->getLeftMesh( )->getElemArea( outerElem );

  SC entry1 = 0.0;
  SC entry2 = 0.0;
  SC entry3 = 0.0;
  SC kernel;
  SCVT phi1x, phi2x, phi3x;

  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

  int i;

  SCVT * x1stref = this->x1stref;
  SCVT * x2stref = this->x2stref;
  SCVT * sx = this->sx;
  SCVT * tx = this->tx;
  SCVT * ux = this->ux;
  SCVT * wxst = this->wxst;

  for ( int innerRot = 0; innerRot < 3; ++innerRot ) {
    // getting local coordinate system
    this->getLocalCoordinates( y1, y2, y3, r1, r2, n, stau, alpha1, alpha2,
      innerRot );
    // getting outer quadrature points in local coordinate system
    this->updateSteinbachLocalQuadratureNodes( x1, x2, x3, y1, y2, y3, r1, r2,
      n, innerRot );

    entry1 = entry2 = entry3 = 0.0;

    // quadrature rule
#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry1, entry2, entry3 ) \
aligned( sx, tx, ux, x1stref, x2stref, wxst : align ) \
private( kernel, phi1x, phi2x, phi3x ) \
simdlen( width )
    for ( i = 0; i < this->qSizeXpad; ++i ) {
      kernel = wxst[ i ] * this->collocation1LayerP1( sx[ i ],
        tx[ i ], ux[ i ], stau, alpha1, alpha2 );

      phi1x = (SCVT) 1.0 - x1stref[ i ] - x2stref[ i ];
      phi2x = x1stref[ i ];
      phi3x = x2stref[ i ];

      entry1 += kernel * phi1x;
      entry2 += kernel * phi2x;
      entry3 += kernel * phi3x;
    }

    matrix.set( 0, innerRot, entry1 * areax );
    matrix.set( 1, innerRot, entry2 * areax );
    matrix.set( 2, innerRot, entry3 * areax );
  }
}

template<class LO, class SC>
SC BEIntegratorLaplace<LO, SC>::collocation1LayerP0(
  SC sx,
  SC tx,
  SC ux,
  SC stau,
  SC alpha1,
  SC alpha2
  ) const {

  SC ret = 0.0;

  ret += this->f1LayerP0( sx, tx, ux, stau, alpha2 );
  ret -= this->f1LayerP0( sx, tx, ux, stau, alpha1 );

  return ( (SC) PI_FACT * ret );
}

template<class LO, class SC>
SC BEIntegratorLaplace<LO, SC>::collocation1LayerP1(
  SC sx,
  SC tx,
  SC ux,
  SC stau,
  SC alpha1,
  SC alpha2
  ) const {

  SC lap = (SCVT) 0.0;

  lap -= this->f1LayerP1( stau, alpha1, tx, sx, ux );
  lap += this->f1LayerP1( stau, alpha2, tx, sx, ux );
  lap /= stau;
  lap *= (SCVT) PI_FACT;

  return lap;
}

// Taken from GOBEM, adapted for both s=s and s=0

template<class LO, class SC>
SC BEIntegratorLaplace<LO, SC>::f1LayerP0(
  SC sx,
  SC tx,
  SC ux,
  SC s,
  SC alpha
  ) const {

  SC f; // return value
  SC beta, gamma, p, q, a2; // substitution variable
  SC tmp1, tmp2, tmp3, tmp4; // aux variable
  SC beta_sq, q_sq; // aux variable

  SC f_0, tmp1_0, tmp2_0, tmp3_0; // variables for s = 0

  SC hh1, hh2, hh3, hh3_0;

  hh2 = ux * ux;

  // first term
  f = -s;
  f_0 = (SC) 0.0;

  // aux variables for substitution
  gamma = tx - alpha * sx;
  beta_sq = (SC) 1.0 + alpha * alpha;
  beta = std::sqrt( beta_sq );
  p = ( alpha * tx + sx ) / beta_sq;
  q_sq = hh2 + gamma * gamma / beta_sq;
  q = std::sqrt( q_sq );
  a2 = beta_sq * q + gamma;

  // aux variables for backward substitution
  tmp1 = beta * ( s - p );
  tmp1_0 = beta * ( -p );
  tmp2 = std::sqrt( tmp1 * tmp1 + q_sq );
  tmp2_0 = std::sqrt( tmp1_0 * tmp1_0 + q_sq );

  // second term
  if ( std::abs( s - sx ) > (SC) EPS ) {
    // aux variables for backward substitution
    tmp3 = alpha * s - tx;

    // stable evaluation
    if ( tmp3 < (SC) 0.0 )
      hh3 = ( hh2 + ( s - sx ) * ( s - sx ) ) / ( tmp2 - tmp3 );
    else
      hh3 = tmp3 + tmp2;

  } else {
    hh3 = (SC) 1.0;
  }

  f += ( s - sx ) * std::log( hh3 );

  if ( std::abs( sx ) > (SC) EPS ) {
    // aux variables for backward substitution
    tmp3_0 = -tx;

    // stable evaluation
    if ( tmp3_0 < (SC) 0.0 )
      hh3_0 = ( hh2 + sx * sx ) / ( tmp2_0 - tmp3_0 );
    else
      hh3_0 = tmp3_0 + tmp2_0;

  } else {
    hh3_0 = (SC) 1.0;
  }

  f_0 += ( -sx ) * std::log( hh3_0 );

  hh1 = gamma / beta;

  // third term
  if ( std::abs( gamma ) > (SC) EPS ) {
    // stable evaluation
    if ( s - p < (SC) 0.0 )
      hh3 = q_sq / ( tmp2 - tmp1 );
    else
      hh3 = tmp1 + tmp2;

    if ( p > (SC) 0.0 )
      hh3_0 = q_sq / ( tmp2_0 - tmp1_0 );
    else
      hh3_0 = tmp1_0 + tmp2_0;

  } else {
    hh3 = hh3_0 = (SC) 1.0;
  }

  f -= hh1 * std::log( hh3 );
  f_0 -= hh1 * std::log( hh3_0 );

  // aux variables for backward substitution
  tmp4 = std::abs( ux );

  // fourth term
  if ( tmp4 > (SC) EPS ) {
    if ( std::abs( tmp1 ) < (SC) EPS )
      hh3 = alpha * q / tmp4;
    else
      hh3 = ( a2 * ( tmp2 - q ) / ( beta * tmp1 ) + alpha * q ) / tmp4;

    if ( std::abs( tmp1_0 ) < (SC) EPS )
      hh3_0 = alpha * q / tmp4;
    else
      hh3_0 = ( a2 * ( tmp2_0 - q ) / ( beta * tmp1_0 ) + alpha * q ) / tmp4;

  } else {
    hh3 = hh3_0 = (SC) 0.0;
  }

  f += (SC) 2.0 * tmp4 * std::atan( hh3 );
  f_0 += (SC) 2.0 * tmp4 * std::atan( hh3_0 );

  return ( f - f_0 );
}

template<class LO, class SC>
SC BEIntegratorLaplace<LO, SC>::f1LayerP1(
  SC stau,
  SC alpha,
  SC tx,
  SC sx,
  SC ux
  ) const {

  // Local variables
  SC beta, beta_2, f, p, q, v, a0, a1, a2, b1, b2, q_2;
  SC a1_0, f_0, b1_0, b2_0, v_0;
  SC hh1, hh2, hh3, hh4, hh5, hh6, hh7, hh8, hh9;

  hh3 = ux * ux;
  hh4 = tx - alpha * sx;
  hh5 = hh4 * hh4;
  hh6 = tx - alpha * stau;
  hh7 = sx - stau * (SC) 2.0;
  hh8 = stau - sx;
  hh9 = hh8 * hh8;

  beta_2 = (SC) 1.0 + alpha * alpha;
  beta = std::sqrt( beta_2 );
  p = ( alpha * tx + sx ) / beta_2;
  q_2 = hh3 + hh5 / beta_2;
  q = std::sqrt( q_2 );
  a1 = std::sqrt( hh9 + hh6 * hh6 + hh3 );
  a1_0 = std::sqrt( sx * sx + tx * tx + hh3 );

  hh1 = (SC) 0.5 / beta_2 * hh4;

  f = (SC) 0.25 + stau * stau + stau * (SC) 0.5 * hh7 + hh1 * ( a1 + q );
  f_0 = (SC) 0.25 + hh1 * ( a1_0 + q );

  if ( sx - stau != (SC) 0.0 ) {
    if ( hh6 > (SC) 0.0 )
      b1 = ( hh9 + hh3 ) / ( a1 + hh6 );
    else
      b1 = a1 - hh6;

    f += -hh8 * (SC) 0.5 * ( stau + hh7 ) * std::log( b1 );
  }

  if ( sx != (SC) 0.0 ) {
    if ( tx > (SC) 0.0 )
      b1_0 = ( sx * sx + hh3 ) / ( a1_0 + tx );
    else
      b1_0 = -tx + a1_0;

    f_0 += sx * (SC) 0.5 * hh7 * std::log( b1_0 );
  }

  b2 = stau - sx - alpha * hh6;
  b2_0 = -sx - alpha * tx;

  if ( b2 < (SC) 0.0 )
    b1 = ( hh3 + hh5 / beta_2 ) / ( a1 - b2 / beta );
  else
    b1 = b2 / beta + a1;

  if ( b2_0 < (SC) 0.0 )
    b1_0 = ( hh3 + hh5 / beta_2 )
    / ( a1_0 - b2_0 / beta );
  else
    b1_0 = b2_0 / beta + a1_0;

  hh2 = (SC) 0.5 / beta * ( alpha * q_2 + hh4 * (SC) 2.0 * ( -hh8 ) );

  if ( b1 != (SC) 0.0 )
    f += hh2 * std::log( b1 );

  if ( b1_0 != (SC) 0.0 )
    f_0 += hh2 * std::log( b1_0 );

  if ( ux != (SC) 0.0 ) {
    f -= hh3 * (SC) 0.5 * std::log( std::abs( q + a1 ) );
    f_0 -= hh3 * (SC) 0.5 * std::log( std::abs( q + a1_0 ) );
    v = beta * ( stau - p ) / ( a1 + q );
    v_0 = -beta * p / ( a1_0 + q );
    a2 = beta_2 * q + hh4;
    a0 = beta_2 * q - hh4;
    a1 = alpha * (SC) 2.0 * beta * q;
    b2 = a2 * v * v + a1 * v + a0;
    b2_0 = a2 * v_0 * v_0 + a1 * v_0 + a0;

    if ( b2 != (SC) 0.0 )
      f -= hh3 * (SC) 0.5 * std::log( std::abs( b2 ) );

    if ( b2_0 != (SC) 0.0 )
      f_0 -= hh3 * (SC) 0.5 * std::log( std::abs( b2_0 ) );

    f += ux * (SC) 2.0 * hh8
      * std::atan( ( a2 * (SC) 2.0 * v + a1 )
      / ( beta * (SC) 2.0 * ux ) );

    f_0 += ux * (SC) 2.0 * hh8
      * std::atan( ( a2 * (SC) 2.0 * v_0 + a1 )
      / ( beta * (SC) 2.0 * ux ) );
  }

  return ( f - f_0 );
}

/*
template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::computeElemMatrix2LayerP0P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  if ( outerElem == innerElem ) {
    matrix.setAll( 0.0 );
    return;
  }

  SC stau, alpha1, alpha2;

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT n[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );

  SCVT * sx = this->sx;
  SCVT * tx = this->tx;
  SCVT * ux = this->ux;
  SCVT * wxst = this->wxst;

  SC entry;
  int i;
  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

  for ( int innerRot = 0; innerRot < 3; ++innerRot ) {
    // getting local coordinate system
    this->getLocalCoordinates( y1, y2, y3, r1, r2, n, stau, alpha1, alpha2,
      innerRot );
    // getting outer quadrature points in local coordinate system
    this->updateSteinbachLocalQuadratureNodes( x1, x2, x3, y1, y2, y3, r1, r2,
      n, innerRot );

    entry = 0.0;

    // quadrature rule
#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry ) \
aligned( sx, tx, ux, wxst : align ) \
simdlen( width )
    for ( i = 0; i < this->qSizeXpad; ++i ) {
      entry += wxst[ i ] * this->collocation2LayerP1( sx[ i ], 
        tx[ i ], ux[ i ], stau, alpha1, alpha2 );
    }

    matrix.set( 0, innerRot, entry );
  }

  matrix.scale( this->space->getLeftMesh( )->getElemArea( outerElem ) );
}
 */
/*
template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::computeElemMatrix2LayerP0P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  SC stau, alpha1, alpha2;

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT n[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );

  SCVT * sx = this->sx;
  SCVT * tx = this->tx;
  SCVT * ux = this->ux;
  SCVT * wxst = this->wxst;

  SC entry1 = 0.0;
  SC entry2 = 0.0;
  SC entry3 = 0.0;
  SC ret1, ret2, ret3;
  int i;
  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

  // getting local coordinate system
  this->getLocalCoordinates( y1, y2, y3, r1, r2, n, stau, alpha1, alpha2 );
  // getting outer quadrature points in local coordinate system
  this->updateSteinbachLocalQuadratureNodes( x1, x2, x3, y1, r1, r2, n );

  // quadrature rule
#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry1, entry2, entry3 ) \
private( ret1, ret2, ret3 ) \
aligned( sx, tx, ux, wxst : align ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXpad; ++i ) {
    this->collocation2LayerP1All( 0, sx[ i ], tx[ i ], ux[ i ],
      stau, alpha1, alpha2, ret1, ret2, ret3 );
    entry1 += wxst[ i ] * ret1;
    entry2 += wxst[ i ] * ret2;
    entry3 += wxst[ i ] * ret3;
  }

  matrix.set( 0, 0, entry1 );
  matrix.set( 0, 1, entry2 );
  matrix.set( 0, 2, entry3 );
  matrix.scale( this->space->getLeftMesh( )->getElemArea( outerElem ) );
}
 */
/*
template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::computeElemMatrix2LayerP0P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT n[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->space->getRightMesh( )->getNormal( innerElem, n );

  int * rot = this->rot;
  SCVT * stau = this->stau;
  SCVT * alpha1 = this->alpha1;
  SCVT * alpha2 = this->alpha2;
  SCVT * sx = this->sx;
  SCVT * tx = this->tx;
  SCVT * ux = this->ux;
  SCVT * wxst = this->wxst;
  const int width = DATA_WIDTH;
  const int align = DATA_ALIGN;

  SC entry1 = 0.0;
  SC entry2 = 0.0;
  SC entry3 = 0.0;
  SC ret1, ret2, ret3;
  int i;

  this->updateSteinbachQuadratureNodes( x1, x2, x3 );
  this->chooseAndUpdateLocalCoordinates( y1, y2, y3, n );

  // quadrature rule
#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry1, entry2, entry3 ) \
private( ret1, ret2, ret3 ) \
aligned( rot, stau, alpha1, alpha2, sx, tx, ux, wxst : align ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXpad; ++i ) {
    this->collocation2LayerP1All( rot[ i ], sx[ i ], tx[ i ], ux[ i ], 
      stau[ i ], alpha1[ i ], alpha2[ i ], ret1, ret2, ret3 );
    entry1 += wxst[ i ] * ret1;
    entry2 += wxst[ i ] * ret2;
    entry3 += wxst[ i ] * ret3;
  }

  matrix.set( 0, 0, entry1 );
  matrix.set( 0, 1, entry2 );
  matrix.set( 0, 2, entry3 );
  matrix.scale( this->space->getLeftMesh( )->getElemArea( outerElem ) );
}
 */
///*

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::computeElemMatrix2LayerP0P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT n[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );

  SC lalpha1, lalpha2, lstau;

  int * rot = this->rot;
  SCVT * stau = this->stau;
  SCVT * alpha1 = this->alpha1;
  SCVT * alpha2 = this->alpha2;
  SCVT * sx = this->sx;
  SCVT * tx = this->tx;
  SCVT * ux = this->ux;
  SCVT * wxst = this->wxst;
  const int width = DATA_WIDTH;
  const int align = DATA_ALIGN;

  SC entry1 = 0.0;
  SC entry2 = 0.0;
  SC entry3 = 0.0;
  SC ret1, ret2, ret3;
  int i;

  // getting local coordinate system
  this->getLocalCoordinates( y1, y2, y3, r1, r2, n, lstau, lalpha1, lalpha2 );
  // getting outer quadrature points in local coordinate system
  this->updateSteinbachLocalQuadratureNodes( x1, x2, x3, y1, r1, r2, n );

  if ( this->checkLocalCoordinates( r1, r2, n ) ) {
//#pragma nounroll
#pragma omp simd linear( i : 1 ) simdlen( width ) aligned( stau, alpha1, alpha2 : align )   
    for ( i = 0; i < this->qSizeXpad; ++i ) {
      this->stau[ i ] = lstau;
      this->alpha1[ i ] = lalpha1;
      this->alpha2[ i ] = lalpha2;
    }
  } else {
    this->chooseAndUpdateLocalCoordinatesFromRot2( y1, y2, y3, n );
  }

  // quadrature rule
#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry1, entry2, entry3 ) \
private( ret1, ret2, ret3 ) \
aligned( rot, stau, alpha1, alpha2, sx, tx, ux, wxst : align ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXpad; ++i ) {
    this->collocation2LayerP1All( rot[ i ], sx[ i ], tx[ i ], ux[ i ],
      stau[ i ], alpha1[ i ], alpha2[ i ], ret1, ret2, ret3 );
    entry1 += wxst[ i ] * ret1;
    entry2 += wxst[ i ] * ret2;
    entry3 += wxst[ i ] * ret3;
  }

  matrix.set( 0, 0, entry1 );
  matrix.set( 0, 1, entry2 );
  matrix.set( 0, 2, entry3 );
  matrix.scale( this->space->getLeftMesh( )->getElemArea( outerElem ) );
}
//*/

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::computeElemMatrix2LayerP0P0(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  if ( outerElem == innerElem ) {
    matrix.setAll( 0.0 );
    return;
  }

  SC stau, alpha1, alpha2;

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT n[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  // getting local coordinate system
  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->getLocalCoordinates( y1, y2, y3, r1, r2, n, stau, alpha1, alpha2 );

  // getting outer quadrature points in local coordinate system
  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->updateSteinbachLocalQuadratureNodes( x1, x2, x3, y1, r1, r2, n );

  SCVT * sx = this->sx;
  SCVT * tx = this->tx;
  SCVT * ux = this->ux;
  SCVT * wxst = this->wxst;

  SC entry = 0.0;
  int i;
  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry ) \
aligned( sx, tx, ux, wxst : align ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXpad; ++i ) {
    entry += wxst[ i ] * this->collocation2LayerP0( sx[ i ], tx[ i ],
      ux[ i ], stau, alpha1, alpha2 );
  }

  entry *= this->space->getLeftMesh( )->getElemArea( outerElem );
  matrix.set( 0, 0, entry );
}

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::computeElemMatrix2LayerP1P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  SC stau, alpha1, alpha2;

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT n[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  SCVT areax = this->space->getLeftMesh( )->getElemArea( outerElem );

  SC entry1 = 0.0;
  SC entry2 = 0.0;
  SC entry3 = 0.0;
  SC kernel;
  SCVT phi1x, phi2x, phi3x;

  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

  int i;

  SCVT * x1stref = this->x1stref;
  SCVT * x2stref = this->x2stref;
  SCVT * sx = this->sx;
  SCVT * tx = this->tx;
  SCVT * ux = this->ux;
  SCVT * wxst = this->wxst;

  for ( int innerRot = 0; innerRot < 3; ++innerRot ) {
    // getting local coordinate system
    this->getLocalCoordinates( y1, y2, y3, r1, r2, n, stau, alpha1, alpha2,
      innerRot );
    // getting outer quadrature points in local coordinate system
    this->updateSteinbachLocalQuadratureNodes( x1, x2, x3, y1, y2, y3, r1, r2,
      n, innerRot );

    entry1 = entry2 = entry3 = 0.0;

    // quadrature rule
#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry1, entry2, entry3 ) \
aligned( sx, tx, ux, x1stref, x2stref, wxst : align ) \
private( kernel, phi1x, phi2x, phi3x ) \
simdlen( width )
    for ( i = 0; i < this->qSizeXpad; ++i ) {
      kernel = wxst[ i ] * this->collocation2LayerP1( sx[ i ],
        tx[ i ], ux[ i ], stau, alpha1, alpha2 );

      phi1x = (SCVT) 1.0 - x1stref[ i ] - x2stref[ i ];
      phi2x = x1stref[ i ];
      phi3x = x2stref[ i ];

      entry1 += kernel * phi1x;
      entry2 += kernel * phi2x;
      entry3 += kernel * phi3x;
    }

    matrix.set( 0, innerRot, entry1 * areax );
    matrix.set( 1, innerRot, entry2 * areax );
    matrix.set( 2, innerRot, entry3 * areax );
  }
}

template<class LO, class SC>
SC BEIntegratorLaplace<LO, SC>::collocation2LayerP1(
  SC sx,
  SC tx,
  SC ux,
  SC stau,
  SC alpha1,
  SC alpha2
  ) const {

  SC ret;

  // helps geometries with planes
  //if ( std::abs( ux ) < (SC) EPS ) {
  //ret = (SC) 0.0;
  //} else {
  SC f1 = this->f2LayerP1( alpha2, tx, sx, ux, stau );
  SC f2 = this->f2LayerP1( alpha1, tx, sx, ux, stau );
  ret = (SC) PI_FACT * ( f1 - f2 ) / stau;
  //}

  return ret;
}

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::collocation2LayerP1All(
  int rot,
  SC sx,
  SC tx,
  SC ux,
  SC stau,
  SC alpha1,
  SC alpha2,
  SC & ret1,
  SC & ret2,
  SC & ret3
  ) const {

  SC aux11, aux12, aux13, aux21, aux22, aux23;

  // helps geometries with planes
  //if ( std::abs( ux ) < (SC) EPS ) {
  //  ret1 = ret2 = ret3 = (SC) 0.0;
  //} else {

  this->f2LayerP1All( alpha2, tx, sx, ux, stau, alpha1, alpha2,
    aux11, aux12, aux13 );
  this->f2LayerP1All( alpha1, tx, sx, ux, stau, alpha1, alpha2,
    aux21, aux22, aux23 );

  aux11 -= aux21;
  aux12 -= aux22;
  aux13 -= aux23;
  aux11 *= (SC) PI_FACT;
  aux12 *= (SC) PI_FACT;
  aux13 *= (SC) PI_FACT;

  if ( rot == 1 ) {
    ret1 = aux12;
    ret2 = aux13;
    ret3 = aux11;
  } else if ( rot == 2 ) {
    ret1 = aux13;
    ret2 = aux11;
    ret3 = aux12;
  } else { //rot = 0
    ret1 = aux11;
    ret2 = aux12;
    ret3 = aux13;
  }
  //}

}

template<class LO, class SC>
SC BEIntegratorLaplace<LO, SC>::collocation2LayerP0(
  SC sx,
  SC tx,
  SC ux,
  SC stau,
  SC alpha1,
  SC alpha2
  ) const {

  SC ret;

  // helps geometries with planes
  //if ( std::abs( ux ) < (SC) EPS ) {
  //ret = (SC) 0.0;
  //} else {
  SC f1 = f2LayerP0( alpha2, tx, sx, ux, stau );
  SC f2 = f2LayerP0( alpha1, tx, sx, ux, stau );
  ret = (SC) PI_FACT * ( f1 - f2 );
  //}

  return ret;
}

template<class LO, class SC>
SC BEIntegratorLaplace<LO, SC>::f2LayerP1(
  SC alpha,
  SC tx,
  SC sx,
  SC ux,
  SC stau
  ) const {

  // Local variables
  SC ret_val = (SC) 0.0;
  SC beta, beta_2, a, b, f, q, f1, f2, h1, h2, h3, h4, q2, aq;
  SC v1, v2, gm, gp, v1m, v2m, v1p, v2p, ask, asl, tau_2, sl_2, tl_2, alpha_2;
  // our modifications for vectorization
  SC hh1, hh2, hh3, hh4, hh5;

  // aux variables
  beta_2 = alpha * alpha + (SC) 1.0;
  beta = std::sqrt( beta_2 );
  asl = alpha * sx;
  ask = alpha * stau;
  tau_2 = ux * ux;
  q2 = tau_2 + ( tx - asl ) * ( tx - asl ) / beta_2;
  q = std::sqrt( q2 );
  sl_2 = sx * sx;
  tl_2 = tx * tx;
  alpha_2 = alpha * alpha;

  // substitution variables v~ (S.3) for s = 0 and s = stau
  a = sl_2 + tl_2 + tau_2;
  v2 = ( stau - sx ) * ( stau - sx ) + ( ask - tx ) * ( ask - tx ) + tau_2;
  v1 = ( -( sx ) - alpha * tx ) / ( std::sqrt( a ) + q ); // S.18 v~ for s = 0
  v2 = ( stau - sx + alpha * ( ask - tx ) ) / ( std::sqrt( v2 ) + q );

  hh1 = ( tx - asl ) / beta_2;
  hh2 = ( tau_2 + alpha_2 * q2 ) / beta_2;

  // aux variables gamma+ and gamma- (S.17)
  // stable evaluation
  if ( ( tx - asl ) >= (SC) 0.0 ) {
    gp = q + hh1;
    gm = hh2 / gp;
  } else {
    gm = q - hh1;
    gp = hh2 / gm;
  }

  // arctan argument   (S.17)
  aq = alpha * q; // aux
  v2p = gp * v2 + aq;
  v2m = gm * v2 - aq;
  v1p = gp * v1 + aq;
  v1m = gm * v1 - aq;

  // stable evaluation
  if ( ( a - tl_2 ) < (SC) 1e-14 ) { //  (a - b < 1e-14)
    a = std::sqrt( a ); // sqrt(sl_2 + tl_2 + tau_2);
    b = alpha * ( sx + alpha * tx ); // onetime evaluation instead of 2 ifs

    hh3 = alpha_2 * ( sl_2 + tau_2 ) - sx * ( sx + alpha * (SC) 2.0 * tx );

    if ( b > (SC) 0.0 ) // (alpha * (sx + alpha * tx) > 0.0)
    {
      v1p = hh3;
      v1p /= alpha * a + sx + alpha * tx;
      v1p = q * v1p + alpha * tau_2 - ( tx - alpha * sx ) * sx;
      v1p /= q + a;
    }
    // if (alpha * (sx + alpha * tx) < 0.0)   // else eingefuegt
    if ( b < (SC) 0.0 ) {
      v1m = hh3;
      v1m /= alpha * a - sx - alpha * tx;
      v1m = q * v1m + alpha * tau_2 - ( tx - alpha * sx ) * sx;
      v1m = -v1m / ( q + a );
    }
  }

  //  log args
  h2 = v2m * v2m + tau_2;
  h1 = v1m * v1m + tau_2;
  h4 = v2p * v2p + tau_2;
  h3 = v1p * v1p + tau_2;

  // log and atan evaluation
  f = ux * (SC) 0.5 * std::log( h4 * h1 / ( h3 * h2 ) );

  f1 = std::atan( v2m / ux ) - std::atan( v1m / ux ); // std::atan
  f2 = std::atan( v2p / ux ) - std::atan( v1p / ux ); // std::atan

  // add together
  f += ( stau - sx ) * ( f1 - f2 );

  // -ux * alpha/beta * u
  h1 = -( sx ) - alpha * tx;
  h2 = stau - sx + alpha * ( ask - tx );

  hh4 = std::sqrt( h1 * h1 + beta_2 * q2 );

  // stable evaluation
  if ( h1 < (SC) 0.0 )
    f1 = beta_2 * q2 / ( hh4 - h1 );
  else
    f1 = h1 + hh4;

  hh5 = std::sqrt( h2 * h2 + beta_2 * q2 );

  // stable evaluation
  if ( h2 < (SC) 0.0 )
    f2 = beta_2 * q2 / ( hh5 - h2 );
  else
    f2 = h2 + hh5;

  if ( std::abs( ux ) > (SC) EPS )
    ret_val = f - ux * alpha / beta * std::log( f2 / f1 );

  return ret_val;
}

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::f2LayerP1All(
  SC alpha,
  SC tx,
  SC sx,
  SC ux,
  SC stau,
  SC alpha1,
  SC alpha2,
  SC & ret1,
  SC & ret2,
  SC & ret3
  ) const {

  ret1 = ret2 = ret3 = (SC) 0.0;
  bool abstau = std::abs( ux ) > (SC) EPS;

  SC beta, beta_2, a, b, q, h1, h2, h3, h4, q2, aq;
  SC v1, v2, gm, gp, v1m, v2m, v1p, v2p, ask, asl, tau_2, sl_2, tl_2, alpha_sq;
  SC aux1, aux2;

  SC hh1, hh2, hh3, hh4, hh5;

  beta_2 = alpha * alpha + (SC) 1.0;
  beta = std::sqrt( beta_2 );
  asl = alpha * sx;
  ask = alpha * stau;
  tau_2 = ux * ux;
  q2 = tau_2 + ( tx - asl ) * ( tx - asl ) / beta_2;
  q = std::sqrt( q2 );
  sl_2 = sx * sx;
  tl_2 = tx * tx;
  alpha_sq = alpha * alpha;

  a = sl_2 + tl_2 + tau_2;
  v2 = ( stau - sx ) * ( stau - sx ) + ( ask - tx ) * ( ask - tx ) + tau_2;
  v1 = ( -( sx ) - alpha * tx ) / ( std::sqrt( a ) + q );
  v2 = ( stau - sx + alpha * ( ask - tx ) ) / ( std::sqrt( v2 ) + q );

  hh1 = ( tx - asl ) / beta_2;
  hh2 = ( tau_2 + alpha_sq * q2 ) / beta_2;

  if ( ( tx - asl ) < (SC) 0.0 ) {
    gm = q - hh1;
    gp = hh2 / gm;
  } else {
    gp = q + hh1;
    gm = hh2 / gp;
  }

  aq = alpha * q;
  v2p = gp * v2 + aq;
  v2m = gm * v2 - aq;
  v1p = gp * v1 + aq;
  v1m = gm * v1 - aq;

  if ( ( a - tl_2 ) < (SC) 1e-14 ) {
    a = std::sqrt( a );
    b = alpha * ( sx + alpha * tx );

    hh3 = alpha_sq * ( sl_2 + tau_2 ) - sx * ( sx + alpha * (SC) 2.0 * tx );

    if ( b > (SC) 0.0 ) {
      v1p = hh3;
      v1p /= alpha * a + sx + alpha * tx;
      v1p = q * v1p + alpha * tau_2 - ( tx - alpha * sx ) * sx;
      v1p /= q + a;
    }
    if ( b < (SC) 0.0 ) {
      v1m = hh3;
      v1m /= alpha * a - sx - alpha * tx;
      v1m = q * v1m + alpha * tau_2 - ( tx - alpha * sx ) * sx;
      v1m = -v1m / ( q + a );
    }
  }

  h2 = v2m * v2m + tau_2;
  h1 = v1m * v1m + tau_2;
  h4 = v2p * v2p + tau_2;
  h3 = v1p * v1p + tau_2;

  aux1 = ux * (SC) 0.5 / stau * std::log( h4 * h1 / ( h3 * h2 ) );
  if ( abstau )
    ret1 = aux1;
  aux1 /= ( alpha2 - alpha1 );
  if ( abstau ) {
    ret2 = -alpha2 * aux1;
    ret3 = alpha1 * aux1;
  }

  aux1 = std::atan( v2m / ux ) - std::atan( v1m / ux );

  aux2 = std::atan( v2p / ux ) - std::atan( v1p / ux );
  aux1 = ( aux1 - aux2 ) / stau;
  if ( abstau )
    ret1 += ( stau - sx ) * aux1;
  aux1 /= ( alpha2 - alpha1 );
  if ( abstau ) {
    ret2 -= ( tx - alpha2 * sx ) * aux1;
    ret3 += ( tx - alpha1 * sx ) * aux1;
  }

  h1 = -( sx ) - alpha * tx;
  h2 = stau - sx + alpha * ( ask - tx );

  hh4 = std::sqrt( h1 * h1 + beta_2 * q2 );

  if ( h1 < (SC) 0.0 )
    aux1 = beta_2 * q2 / ( hh4 - h1 );
  else
    aux1 = h1 + hh4;

  hh5 = std::sqrt( h2 * h2 + beta_2 * q2 );

  if ( h2 < (SC) 0.0 )
    aux2 = beta_2 * q2 / ( hh5 - h2 );
  else
    aux2 = h2 + hh5;

  aux1 = -ux / ( beta * stau ) * std::log( aux2 / aux1 );
  if ( abstau )
    ret1 += alpha * aux1;
  aux1 /= ( alpha2 - alpha1 );
  if ( abstau ) {
    ret2 -= ( 1 + alpha * alpha2 ) * aux1;
    ret3 += ( 1 + alpha * alpha1 ) * aux1;
  }

}

template<class LO, class SC>
SC BEIntegratorLaplace<LO, SC>::f2LayerP0(
  SC alpha,
  SC tx,
  SC sx,
  SC ux,
  SC stau
  ) const {

  // Local variables
  SC beta_2, a, b, q, f1, f2, q2, aq;
  SC v1, v2, gm, gp, v1m, v2m, v1p, v2p, ask, asl, tau_2, sl_2, tl_2, alpha_2;
  SC sk_minus_sl, ask_minus_tl, tl_minus_asl, sl_plus_alpha_tl;

  SC hh1, hh2, hh3, hh4, hh5;

  // Aux variables
  beta_2 = alpha * alpha + (SC) 1.0;
  asl = alpha * sx;
  ask = alpha * stau;
  tau_2 = ux * ux;
  tl_minus_asl = tx - asl;
  q2 = tau_2 + tl_minus_asl * tl_minus_asl / beta_2;
  q = sqrt( q2 );
  sl_2 = sx * sx;
  tl_2 = tx * tx;
  alpha_2 = alpha * alpha;
  sk_minus_sl = stau - sx;
  ask_minus_tl = ask - tx;
  sl_plus_alpha_tl = sx + alpha * tx;

  // Substitution v~ (S.3) for s = 0 and s = stau
  a = sl_2 + tl_2 + tau_2;
  v2 = sk_minus_sl * sk_minus_sl + ask_minus_tl * ask_minus_tl + tau_2;
  v1 = -sl_plus_alpha_tl / ( std::sqrt( a ) + q ); // S.18 v~ for s = 0
  v2 = ( sk_minus_sl + alpha * ask_minus_tl ) / ( std::sqrt( v2 ) + q );

  // Aux variables gamma+ und gamma- (S.17)
  // Stability
  hh1 = tau_2 + alpha_2 * q2;
  hh2 = tl_minus_asl / beta_2;
  if ( tl_minus_asl >= (SC) 0.0 ) {
    gp = q + hh2;
    gm = hh1 / ( gp * beta_2 );
  } else {
    gm = q - hh2;
    gp = hh1 / ( gm * beta_2 );
  }

  // Args for the arctan   (S.17)
  aq = alpha * q; // Aux
  v2p = gp * v2 + aq;
  v2m = gm * v2 - aq;
  v1p = gp * v1 + aq;
  v1m = gm * v1 - aq;

  // Stability
  if ( a - tl_2 < (SC) 1e-14 ) //  (a - b < 1e-14)
  {
    a = std::sqrt( a ); // sqrt(sl_2 + tl_2 + tau_2);
    b = alpha * sl_plus_alpha_tl;
    hh3 = alpha_2 * ( sl_2 + tau_2 ) - sx * ( sx + alpha * (SC) 2.0 * tx );
    hh4 = alpha * tau_2 - tl_minus_asl * sx;
    hh5 = q + a;

    if ( b > (SC) 0.0 ) // (alpha * (sx + alpha * tx) > 0.0)
    {
      v1p = hh3;
      v1p /= alpha * a + sx + alpha * tx;
      v1p = q * v1p + hh4;
      v1p /= hh5;
    }
    if ( b < (SC) 0.0 ) // if (alpha * (sx + alpha * tx) < 0.0)
    {
      v1m = hh3;
      v1m /= alpha * a - sx - alpha * tx;
      v1m = q * v1m + hh4;
      v1m = -v1m / hh5;
    }
  }

  // Atans
  f1 = std::atan( v2m / ux ) - std::atan( v1m / ux );
  f2 = std::atan( v2p / ux ) - std::atan( v1p / ux );

  return ( f1 - f2 );
}

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::computeElemMatrixHypersingularP1P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  SCVT * outerCurl =
    this->space->getLeftMesh( )->getCurls( )->getData( ) + outerElem * 9;
  SCVT * innerCurl =
    this->space->getRightMesh( )->getCurls( )->getData( ) + innerElem * 9;

  SCVT curlDot1, curlDot2, curlDot3;
  SC * matrixData = matrix.getData( );

  const int width = DATA_WIDTH;

#pragma omp simd \
private( curlDot1, curlDot2, curlDot3 ) \
simdlen( width )
  for ( int innerRot = 0; innerRot < 3; ++innerRot ) {
    curlDot1 = outerCurl[ 0 ] * ( innerCurl + 3 * innerRot )[ 0 ]
      + outerCurl[ 1 ] * ( innerCurl + 3 * innerRot )[ 1 ]
      + outerCurl[ 2 ] * ( innerCurl + 3 * innerRot )[ 2 ];
    curlDot2 = ( outerCurl + 3 )[ 0 ] * ( innerCurl + 3 * innerRot )[ 0 ]
      + ( outerCurl + 3 )[ 1 ] * ( innerCurl + 3 * innerRot )[ 1 ]
      + ( outerCurl + 3 )[ 2 ] * ( innerCurl + 3 * innerRot )[ 2 ];
    curlDot3 = ( outerCurl + 6 )[ 0 ] * ( innerCurl + 3 * innerRot )[ 0 ]
      + ( outerCurl + 6 )[ 1 ] * ( innerCurl + 3 * innerRot )[ 1 ]
      + ( outerCurl + 6 )[ 2 ] * ( innerCurl + 3 * innerRot )[ 2 ];

    matrixData[ innerRot * 3 ] = curlDot1;
    matrixData[ innerRot * 3 + 1 ] = curlDot2;
    matrixData[ innerRot * 3 + 2 ] = curlDot3;
  }

  // HACK: computeElemMatrix1LayerP0P0 uses wrong BESpace, but only uses mesh, 
  // which is independent of ansatz/test types
  FullMatrix<LO, SC> localV( 1, 1 );
  // use basic Gaussian quadrature for disjoint elements
  if ( this->quadratureOrderDisjointElems &&
    this->areElementsDisjoint( outerElem, innerElem ) ) {
    this->computeElemMatrix1LayerDisjointP0P0( outerElem, innerElem, localV );
  } else {
    switch ( this->quadrature ) {
      case SauterSchwab:
        this->computeElemMatrix1LayerSauterSchwabP0P0( outerElem, innerElem,
          localV );
        break;
      case Steinbach:
        this->computeElemMatrix1LayerP0P0( outerElem, innerElem, localV );
        break;
    }
  }
  matrix.scale( localV.get( 0, 0 ) );
}

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::representationFormulaP1P0(
  const SCVT* xCoord,
  LO numPoints,
  const Vector<LO, SC> & dir,
  const Vector<LO, SC> & neu,
  bool interior,
  Vector<LO, SC> & values
  ) const {

  LO nElems = this->space->getMesh( )->getNElements( );
  SCVT stau, alpha1, alpha2;
  SCVT x1[3], x2[3], x3[3], r1[3], r2[3], n[3], xLocal[3];
  LO ind[3];
  SC val;

  for ( LO i = 0; i < numPoints; i++ ) {
    val = 0.0;
    const SCVT * xSingle = xCoord + 3 * i;

    for ( LO j = 0; j < nElems; j++ ) {
      this->space->getMesh( )->getNodes( j, x1, x2, x3 );
      this->space->getMesh( )->getElement( j, ind );

      for ( int rot = 0; rot < 3; rot++ ) {
        //this->getLocalCoordinates( x1, x2, x3, r1, r2, n, stau, alpha1, alpha2, rot );
        this->getLocalCoordinates( j, r1, r2, n, stau, alpha1, alpha2, rot );
        this->getLocalQuadratureNode( x1, x2, x3, xSingle, r1, r2, n, xLocal, rot );
        if ( rot == 0 ) {
          val += neu.get( j ) * collocation1LayerP0( xLocal[ 0 ], xLocal[ 1 ],
            xLocal[ 2 ], stau, alpha1, alpha2 );
        }
        val -= dir.get( ind[ rot ] ) * collocation2LayerP1( xLocal[ 0 ],
          xLocal[ 1 ], xLocal[ 2 ], stau, alpha1, alpha2 );
      }
    }

    values.set( i, val );
  }

  if ( !interior ) {
    values.scale( -1.0 );
  }
}

template<class LO, class SC>
void BEIntegratorLaplace<LO, SC>::doubleLayerPotentialP1(
  const SCVT* x,
  LO nPoints,
  const Vector<LO, SC> & density,
  Vector<LO, SC> & values
  ) const {

  LO nElems = this->space->getMesh( )->getNElements( );
  SCVT stau, alpha1, alpha2;
  SCVT x1[3], x2[3], x3[3], r1[3], r2[3], n[3], xLocal[3];
  LO ind[3];
  SC val;

  for ( LO i = 0; i < nPoints; i++ ) {
    val = 0.0;
    const SCVT * xSingle = x + 3 * i;

    for ( LO j = 0; j < nElems; j++ ) {
      this->space->getMesh( )->getNodes( j, x1, x2, x3 );
      this->space->getMesh( )->getElement( j, ind );

      for ( int rot = 0; rot < 3; rot++ ) {
        this->getLocalCoordinates( j, r1, r2, n, stau, alpha1, alpha2, rot );
        this->getLocalQuadratureNode( x1, x2, x3, xSingle, r1, r2, n, xLocal,
          rot );
        val += density.get( ind[ rot ] ) * collocation2LayerP1( xLocal[ 0 ],
          xLocal[ 1 ], xLocal[ 2 ], stau, alpha1, alpha2 );
      }
    }

    values.set( i, val );
  }
}

}
#endif
