/*!
 * @file    BEIntegratorHelmholtz.cpp
 * @author  Jan Zapletal
 * @date    August 12, 2013
 * 
 */

#ifdef BEINTEGRATORHELMHOLTZ_H

namespace bem4i {

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix1Layer(
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
    // use something a little more sophisticated for close panels
    switch ( this->quadrature ) {
      case SauterSchwab:
        this->computeElemMatrix1LayerSauterSchwabP1P1( outerElem, innerElem,
          matrix );
        break;
      case Steinbach:
        this->computeElemMatrix1LayerP1P1( outerElem, innerElem, matrix );
        break;
    }
  } else {
    std::cout << "Not implemented!" << std::endl;
  }
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix2Layer(
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
    switch ( this->quadrature ) {
      case SauterSchwab:
        this->computeElemMatrix2LayerSauterSchwabP0P1( outerElem, innerElem,
          matrix );
        break;
      case Steinbach:
        computeElemMatrix2LayerP0P1( outerElem, innerElem, matrix );
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
    // use something a little more sophisticated for close panels
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
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrixHypersingular(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC>& matrix
  )const {

  if ( ( this->space->getTestFunctionType( ) == p1 ||
    this->space->getTestFunctionType( ) == p1dis ) &&
    ( this->space->getAnsatzFunctionType( ) == p1 ||
    this->space->getAnsatzFunctionType( ) == p1dis ) ) {
    // use basic Gaussian quadrature for disjoint elements
    if ( this->quadratureOrderDisjointElems &&
      this->areElementsDisjoint( outerElem, innerElem ) ) {
      this->computeElemMatrixHypersingularDisjointP1P1( outerElem, innerElem,
        matrix );
      return;
    }
    // use something a little more sophisticated for close panels
    switch ( this->quadrature ) {
      case SauterSchwab:
        this->computeElemMatrixHypersingularSauterSchwabP1P1( outerElem,
          innerElem, matrix );
        break;
      case Steinbach:
        this->computeElemMatrixHypersingularP1P1( outerElem, innerElem,
          matrix );
        break;
    }
  } else {
    std::cout << "Not implemented!" << std::endl;
  }
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix1LayerSauterSchwabP0P0(
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

  this->getSpace( )->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->getSpace( )->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );

  int nSimplex;
  SCVT ** jacobianP;

  this->getReferenceSauterSchwabData( type, nSimplex, jacobianP );

  SCVT * jacobian;

  SCVT entryRe = 0.0;
  SCVT entryIm = 0.0;
  SCVT kernelRe, kernelIm, mult;

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
reduction( + : entryRe, entryIm ) \
aligned( x1ss, x2ss, x3ss, y1ss, y2ss, y3ss, jacobian : align ) \
private( kernelRe, kernelIm ) \
simdlen( width )
    for ( i = 0; i < totalSize; ++i ) {
      this->evalSingleLayerKernel( x1ss[ i ], x2ss[ i ], x3ss[ i ], y1ss[ i ],
        y2ss[ i ], y3ss[ i ], kernelRe, kernelIm );
      entryRe += jacobian[ i ] * kernelRe;
      entryIm += jacobian[ i ] * kernelIm;
    }
  }

  SCVT areaMult = (SCVT) 4.0 *
    this->getSpace( )->getRightMesh( )->getElemArea( innerElem ) *
    this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );

  matrix.set( 0, 0, areaMult * SC( entryRe, entryIm ) );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix1LayerSauterSchwabP1P1(
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

  SCVT entry11Re = 0.0;
  SCVT entry21Re = 0.0;
  SCVT entry31Re = 0.0;
  SCVT entry12Re = 0.0;
  SCVT entry22Re = 0.0;
  SCVT entry32Re = 0.0;
  SCVT entry13Re = 0.0;
  SCVT entry23Re = 0.0;
  SCVT entry33Re = 0.0;
  SCVT entry11Im = 0.0;
  SCVT entry21Im = 0.0;
  SCVT entry31Im = 0.0;
  SCVT entry12Im = 0.0;
  SCVT entry22Im = 0.0;
  SCVT entry32Im = 0.0;
  SCVT entry13Im = 0.0;
  SCVT entry23Im = 0.0;
  SCVT entry33Im = 0.0;

  SCVT kernelRe, kernelIm;
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
reduction( + : entry11Re, entry11Im ) \
reduction( + : entry12Re, entry12Im ) \
reduction( + : entry13Re, entry13Im ) \
reduction( + : entry21Re, entry21Im ) \
reduction( + : entry22Re, entry22Im ) \
reduction( + : entry23Re, entry23Im ) \
reduction( + : entry31Re, entry31Im ) \
reduction( + : entry32Re, entry32Im ) \
reduction( + : entry33Re, entry33Im ) \
aligned( x1ss, x2ss, x3ss, y1ss, y2ss, y3ss, jacobian : align ) \
aligned( x1ref, x2ref, y1ref, y2ref : align ) \
private( kernelRe, kernelIm, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y ) \
simdlen( width )
    for ( i = 0; i < totalSize; ++i ) {
      this->evalSingleLayerKernel( x1ss[ i ], x2ss[ i ], x3ss[ i ], y1ss[ i ],
        y2ss[ i ], y3ss[ i ], kernelRe, kernelIm );

      kernelRe *= jacobian[ i ];
      kernelIm *= jacobian[ i ];

      phi1x = (SCVT) 1.0 - x1ref[ i ];
      phi2x = x1ref[ i ] - x2ref[ i ];
      phi3x = x2ref[ i ];
      phi1y = (SCVT) 1.0 - y1ref[ i ];
      phi2y = y1ref[ i ] - y2ref[ i ];
      phi3y = y2ref[ i ];

      entry11Re += kernelRe * phi1x * phi1y;
      entry21Re += kernelRe * phi2x * phi1y;
      entry31Re += kernelRe * phi3x * phi1y;
      entry12Re += kernelRe * phi1x * phi2y;
      entry22Re += kernelRe * phi2x * phi2y;
      entry32Re += kernelRe * phi3x * phi2y;
      entry13Re += kernelRe * phi1x * phi3y;
      entry23Re += kernelRe * phi2x * phi3y;
      entry33Re += kernelRe * phi3x * phi3y;
      entry11Im += kernelIm * phi1x * phi1y;
      entry21Im += kernelIm * phi2x * phi1y;
      entry31Im += kernelIm * phi3x * phi1y;
      entry12Im += kernelIm * phi1x * phi2y;
      entry22Im += kernelIm * phi2x * phi2y;
      entry32Im += kernelIm * phi3x * phi2y;
      entry13Im += kernelIm * phi1x * phi3y;
      entry23Im += kernelIm * phi2x * phi3y;
      entry33Im += kernelIm * phi3x * phi3y;
    }
  }

  // first two y indices swapped for common edge!
  if ( type == 1 ) {
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ),
      SC( entry11Re, entry11Im ) );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ),
      SC( entry21Re, entry21Im ) );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ),
      SC( entry31Re, entry31Im ) );
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ),
      SC( entry12Re, entry12Im ) );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ),
      SC( entry22Re, entry22Im ) );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ),
      SC( entry32Re, entry32Im ) );
  } else {
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ),
      SC( entry11Re, entry11Im ) );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ),
      SC( entry21Re, entry21Im ) );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ),
      SC( entry31Re, entry31Im ) );
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ),
      SC( entry12Re, entry12Im ) );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ),
      SC( entry22Re, entry22Im ) );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ),
      SC( entry32Re, entry32Im ) );
  }
  matrix.set( this->mod3( outerRot ), this->mod3( 2 + innerRot ),
    SC( entry13Re, entry13Im ) );
  matrix.set( this->mod3( 1 + outerRot ), this->mod3( 2 + innerRot ),
    SC( entry23Re, entry23Im ) );
  matrix.set( this->mod3( 2 + outerRot ), this->mod3( 2 + innerRot ),
    SC( entry33Re, entry33Im ) );

  SCVT areaMult = (SCVT) 4.0 *
    this->getSpace( )->getRightMesh( )->getElemArea( innerElem ) *
    this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );
  matrix.scale( areaMult );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix1LayerDisjointP0P0(
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
  SCVT entryRe = 0.0;
  SCVT entryIm = 0.0;
  SCVT kernelRe, kernelIm;

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
reduction( + : entryRe, entryIm ) \
aligned( x1d, x2d, x3d, y1d, y2d, y3d : align ) \
aligned( wxyd : align ) \
private( kernelRe, kernelIm ) \
simdlen( width )
  for ( i = 0; i < nIters; ++i ) {
    this->evalSingleLayerKernel( x1d[ i ], x2d[ i ], x3d[ i ],
      y1d[ i ], y2d[ i ], y3d[ i ], kernelRe, kernelIm );
    entryRe += wxyd[ i ] * kernelRe;
    entryIm += wxyd[ i ] * kernelIm;
  }

  SCVT areaMult = this->space->getLeftMesh( )->getElemArea( outerElem ) *
    this->space->getRightMesh( )->getElemArea( innerElem );

  entryRe *= areaMult;
  entryIm *= areaMult;
  elemMatrix.set( 0, 0, SC( entryRe, entryIm ) );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix1LayerDisjointP1P1(
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
  SCVT entry11Re = 0.0;
  SCVT entry21Re = 0.0;
  SCVT entry31Re = 0.0;
  SCVT entry12Re = 0.0;
  SCVT entry22Re = 0.0;
  SCVT entry32Re = 0.0;
  SCVT entry13Re = 0.0;
  SCVT entry23Re = 0.0;
  SCVT entry33Re = 0.0;
  SCVT entry11Im = 0.0;
  SCVT entry21Im = 0.0;
  SCVT entry31Im = 0.0;
  SCVT entry12Im = 0.0;
  SCVT entry22Im = 0.0;
  SCVT entry32Im = 0.0;
  SCVT entry13Im = 0.0;
  SCVT entry23Im = 0.0;
  SCVT entry33Im = 0.0;

  SCVT kernelRe, kernelIm;

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
reduction( + : entry11Re, entry11Im ) \
reduction( + : entry12Re, entry12Im ) \
reduction( + : entry13Re, entry13Im ) \
reduction( + : entry21Re, entry21Im ) \
reduction( + : entry22Re, entry22Im ) \
reduction( + : entry23Re, entry23Im ) \
reduction( + : entry31Re, entry31Im ) \
reduction( + : entry32Re, entry32Im ) \
reduction( + : entry33Re, entry33Im ) \
aligned( x1d, x2d, x3d, x1dref, x2dref : align ) \
aligned( y1d, y2d, y3d, y1dref, y2dref : align ) \
aligned( wxyd : align ) \
private( kernelRe, kernelIm, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y ) \
simdlen( width )
  for ( i = 0; i < nIters; ++i ) {

    this->evalSingleLayerKernel( x1d[ i ], x2d[ i ], x3d[ i ], y1d[ i ],
      y2d[ i ], y3d[ i ], kernelRe, kernelIm );

    kernelRe *= wxyd[ i ];
    kernelIm *= wxyd[ i ];

    phi1y = (SCVT) 1.0 - y1dref[ i ] - y2dref[ i ];
    phi2y = y1dref[ i ];
    phi3y = y2dref[ i ];
    phi1x = (SCVT) 1.0 - x1dref[ i ] - x2dref[ i ];
    phi2x = x1dref[ i ];
    phi3x = x2dref[ i ];

    entry11Re += kernelRe * phi1x * phi1y;
    entry21Re += kernelRe * phi2x * phi1y;
    entry31Re += kernelRe * phi3x * phi1y;
    entry12Re += kernelRe * phi1x * phi2y;
    entry22Re += kernelRe * phi2x * phi2y;
    entry32Re += kernelRe * phi3x * phi2y;
    entry13Re += kernelRe * phi1x * phi3y;
    entry23Re += kernelRe * phi2x * phi3y;
    entry33Re += kernelRe * phi3x * phi3y;
    entry11Im += kernelIm * phi1x * phi1y;
    entry21Im += kernelIm * phi2x * phi1y;
    entry31Im += kernelIm * phi3x * phi1y;
    entry12Im += kernelIm * phi1x * phi2y;
    entry22Im += kernelIm * phi2x * phi2y;
    entry32Im += kernelIm * phi3x * phi2y;
    entry13Im += kernelIm * phi1x * phi3y;
    entry23Im += kernelIm * phi2x * phi3y;
    entry33Im += kernelIm * phi3x * phi3y;
  }

  SCVT areasM = this->space->getLeftMesh( )->getElemArea( outerElem ) *
    this->space->getRightMesh( )->getElemArea( innerElem );

  elemMatrix.set( 0, 0, SC( entry11Re, entry11Im ) );
  elemMatrix.set( 1, 0, SC( entry21Re, entry21Im ) );
  elemMatrix.set( 2, 0, SC( entry31Re, entry31Im ) );
  elemMatrix.set( 0, 1, SC( entry12Re, entry12Im ) );
  elemMatrix.set( 1, 1, SC( entry22Re, entry22Im ) );
  elemMatrix.set( 2, 1, SC( entry32Re, entry32Im ) );
  elemMatrix.set( 0, 2, SC( entry13Re, entry13Im ) );
  elemMatrix.set( 1, 2, SC( entry23Re, entry23Im ) );
  elemMatrix.set( 2, 2, SC( entry33Re, entry33Im ) );
  elemMatrix.scale( areasM );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix2LayerSauterSchwabP0P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  if ( outerElem == innerElem ) {
    matrix.setAll( 0.0 );
    return;
  }

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
  SCVT * y1ref;
  SCVT * y2ref;

  SCVT entry1Re = 0.0;
  SCVT entry2Re = 0.0;
  SCVT entry3Re = 0.0;
  SCVT entry1Im = 0.0;
  SCVT entry2Im = 0.0;
  SCVT entry3Im = 0.0;
  SCVT kernelRe, kernelIm;
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
    y1ref = y1refP[ simplex ];
    y2ref = y2refP[ simplex ];

    int i = 0;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry1Re, entry2Re, entry3Re ) \
reduction( + : entry1Im, entry2Im, entry3Im ) \
aligned( x1ss, x2ss, x3ss, y1ss, y2ss, y3ss, jacobian : align ) \
aligned( y1ref, y2ref : align ) \
private( kernelRe, kernelIm ) \
simdlen( width )
    for ( i = 0; i < totalSize; ++i ) {
      this->evalDoubleLayerKernel( x1ss[ i ], x2ss[ i ], x3ss[ i ], y1ss[ i ],
        y2ss[ i ], y3ss[ i ], n[ 0 ], n[ 1 ], n[ 2 ], kernelRe, kernelIm );
      entry1Re += jacobian[ i ] * kernelRe * ( (SCVT) 1.0 - y1ref[ i ] );
      entry2Re += jacobian[ i ] * kernelRe * ( y1ref[ i ] - y2ref[ i ] );
      entry3Re += jacobian[ i ] * kernelRe * y2ref[ i ];
      entry1Im += jacobian[ i ] * kernelIm * ( (SCVT) 1.0 - y1ref[ i ] );
      entry2Im += jacobian[ i ] * kernelIm * ( y1ref[ i ] - y2ref[ i ] );
      entry3Im += jacobian[ i ] * kernelIm * y2ref[ i ];
    }
  }

  if ( type == 1 ) {
    matrix.set( 0, this->mod3( 1 + innerRot ), SC( entry1Re, entry1Im ) );
    matrix.set( 0, this->mod3( innerRot ), SC( entry2Re, entry2Im ) );
  } else {
    matrix.set( 0, this->mod3( innerRot ), SC( entry1Re, entry1Im ) );
    matrix.set( 0, this->mod3( 1 + innerRot ), SC( entry2Re, entry2Im ) );
  }
  matrix.set( 0, this->mod3( 2 + innerRot ), SC( entry3Re, entry3Im ) );

  SCVT areaMult = (SCVT) 4.0 *
    this->getSpace( )->getRightMesh( )->getElemArea( innerElem ) *
    this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );
  matrix.scale( areaMult );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix2LayerSauterSchwabP0P0(
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

  SCVT entryRe = 0.0;
  SCVT entryIm = 0.0;
  SCVT kernelRe, kernelIm;

  for ( int simplex = 0; simplex < nSimplex; ++simplex ) {
    this->updateSauterSchwabQuadratureNodes( type, simplex, x1, x2, x3,
      outerRot, y1, y2, y3, innerRot );
    jacobian = jacobianP[ simplex ];

    int i = 0;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entryRe, entryIm ) \
aligned( x1ss, x2ss, x3ss, y1ss, y2ss, y3ss, jacobian : align ) \
private( kernelRe, kernelIm ) \
simdlen( width )
    for ( i = 0; i < totalSize; ++i ) {
      this->evalDoubleLayerKernel( x1ss[ i ], x2ss[ i ], x3ss[ i ], y1ss[ i ],
        y2ss[ i ], y3ss[ i ], n[ 0 ], n[ 1 ], n[ 2 ], kernelRe, kernelIm );

      entryRe += kernelRe * jacobian[ i ];
      entryIm += kernelIm * jacobian[ i ];
    }
  }

  SCVT innerArea = this->getSpace( )->getRightMesh( )->getElemArea( innerElem );
  SCVT outerArea = this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );
  SCVT mult = (SCVT) 4.0 * innerArea * outerArea;

  entryRe *= mult;
  entryIm *= mult;
  matrix.set( 0, 0, SC( entryRe, entryIm ) );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix2LayerSauterSchwabP1P1(
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

  SCVT entry11Re = 0.0;
  SCVT entry21Re = 0.0;
  SCVT entry31Re = 0.0;
  SCVT entry12Re = 0.0;
  SCVT entry22Re = 0.0;
  SCVT entry32Re = 0.0;
  SCVT entry13Re = 0.0;
  SCVT entry23Re = 0.0;
  SCVT entry33Re = 0.0;
  SCVT entry11Im = 0.0;
  SCVT entry21Im = 0.0;
  SCVT entry31Im = 0.0;
  SCVT entry12Im = 0.0;
  SCVT entry22Im = 0.0;
  SCVT entry32Im = 0.0;
  SCVT entry13Im = 0.0;
  SCVT entry23Im = 0.0;
  SCVT entry33Im = 0.0;

  SCVT kernelRe, kernelIm;
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
reduction( + : entry11Re, entry11Im ) \
reduction( + : entry12Re, entry12Im ) \
reduction( + : entry13Re, entry13Im ) \
reduction( + : entry21Re, entry21Im ) \
reduction( + : entry22Re, entry22Im ) \
reduction( + : entry23Re, entry23Im ) \
reduction( + : entry31Re, entry31Im ) \
reduction( + : entry32Re, entry32Im ) \
reduction( + : entry33Re, entry33Im ) \
aligned( x1ss, x2ss, x3ss, y1ss, y2ss, y3ss, jacobian : align ) \
aligned( x1ref, x2ref, y1ref, y2ref : align ) \
private( kernelRe, kernelIm, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y ) \
simdlen( width )
    for ( i = 0; i < totalSize; ++i ) {
      this->evalDoubleLayerKernel( x1ss[ i ], x2ss[ i ], x3ss[ i ], y1ss[ i ],
        y2ss[ i ], y3ss[ i ], n[ 0 ], n[ 1 ], n[ 2 ], kernelRe, kernelIm );

      kernelRe *= jacobian[ i ];
      kernelIm *= jacobian[ i ];

      phi1x = (SCVT) 1.0 - x1ref[ i ];
      phi2x = x1ref[ i ] - x2ref[ i ];
      phi3x = x2ref[ i ];
      phi1y = (SCVT) 1.0 - y1ref[ i ];
      phi2y = y1ref[ i ] - y2ref[ i ];
      phi3y = y2ref[ i ];

      entry11Re += kernelRe * phi1x * phi1y;
      entry21Re += kernelRe * phi2x * phi1y;
      entry31Re += kernelRe * phi3x * phi1y;
      entry12Re += kernelRe * phi1x * phi2y;
      entry22Re += kernelRe * phi2x * phi2y;
      entry32Re += kernelRe * phi3x * phi2y;
      entry13Re += kernelRe * phi1x * phi3y;
      entry23Re += kernelRe * phi2x * phi3y;
      entry33Re += kernelRe * phi3x * phi3y;
      entry11Im += kernelIm * phi1x * phi1y;
      entry21Im += kernelIm * phi2x * phi1y;
      entry31Im += kernelIm * phi3x * phi1y;
      entry12Im += kernelIm * phi1x * phi2y;
      entry22Im += kernelIm * phi2x * phi2y;
      entry32Im += kernelIm * phi3x * phi2y;
      entry13Im += kernelIm * phi1x * phi3y;
      entry23Im += kernelIm * phi2x * phi3y;
      entry33Im += kernelIm * phi3x * phi3y;
    }
  }

  // first two y indices swapped for common edge!
  if ( type == 1 ) {
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ),
      SC( entry11Re, entry11Im ) );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ),
      SC( entry21Re, entry21Im ) );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ),
      SC( entry31Re, entry31Im ) );
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ),
      SC( entry12Re, entry12Im ) );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ),
      SC( entry22Re, entry22Im ) );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ),
      SC( entry32Re, entry32Im ) );
  } else {
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ),
      SC( entry11Re, entry11Im ) );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ),
      SC( entry21Re, entry21Im ) );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ),
      SC( entry31Re, entry31Im ) );
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ),
      SC( entry12Re, entry12Im ) );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ),
      SC( entry22Re, entry22Im ) );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ),
      SC( entry32Re, entry32Im ) );
  }
  matrix.set( this->mod3( outerRot ), this->mod3( 2 + innerRot ),
    SC( entry13Re, entry13Im ) );
  matrix.set( this->mod3( 1 + outerRot ), this->mod3( 2 + innerRot ),
    SC( entry23Re, entry23Im ) );
  matrix.set( this->mod3( 2 + outerRot ), this->mod3( 2 + innerRot ),
    SC( entry33Re, entry33Im ) );


  SCVT areaMult = (SCVT) 4.0 *
    this->getSpace( )->getRightMesh( )->getElemArea( innerElem ) *
    this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );
  matrix.scale( areaMult );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix2LayerDisjointP0P1(
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
  SCVT entry1Re = 0.0;
  SCVT entry1Im = 0.0;
  SCVT entry2Re = 0.0;
  SCVT entry2Im = 0.0;
  SCVT entry3Re = 0.0;
  SCVT entry3Im = 0.0;
  SCVT kernelRe, kernelIm;

  SCVT phi1y, phi2y, phi3y;

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
reduction( + : entry1Re, entry2Re, entry3Re ) \
reduction( + : entry1Im, entry2Im, entry3Im ) \
aligned( x1d, x2d, x3d, y1d, y2d, y3d : align ) \
aligned( y1dref, y2dref, wxyd : align ) \
private( kernelRe, kernelIm, phi1y, phi2y, phi3y ) \
simdlen( width )
  for ( i = 0; i < nIters; ++i ) {

    this->evalDoubleLayerKernel( x1d[ i ], x2d[ i ], x3d[ i ],
      y1d[ i ], y2d[ i ], y3d[ i ], n[ 0 ], n[ 1 ], n[ 2 ],
      kernelRe, kernelIm );

    phi1y = (SCVT) 1.0 - y1dref[ i ] - y2dref[ i ];
    phi2y = y1dref[ i ];
    phi3y = y2dref[ i ];

    kernelRe *= wxyd[ i ];
    kernelIm *= wxyd[ i ];
    entry1Re += kernelRe * phi1y;
    entry1Im += kernelIm * phi1y;
    entry2Re += kernelRe * phi2y;
    entry2Im += kernelIm * phi2y;
    entry3Re += kernelRe * phi3y;
    entry3Im += kernelIm * phi3y;
  }

  SCVT areasM = this->space->getLeftMesh( )->getElemArea( outerElem ) *
    this->space->getRightMesh( )->getElemArea( innerElem );
  elemMatrix.set( 0, 0, areasM * SC( entry1Re, entry1Im ) );
  elemMatrix.set( 0, 1, areasM * SC( entry2Re, entry2Im ) );
  elemMatrix.set( 0, 2, areasM * SC( entry3Re, entry3Im ) );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix2LayerDisjointP0P0(
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
  SCVT entryRe = 0.0;
  SCVT entryIm = 0.0;
  SCVT kernelRe, kernelIm;

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
reduction( + : entryRe, entryIm ) \
aligned( x1d, x2d, x3d, y1d, y2d, y3d : align ) \
aligned( wxyd : align ) \
private( kernelRe, kernelIm ) \
simdlen( width )
  for ( i = 0; i < nIters; ++i ) {
    this->evalDoubleLayerKernel( x1d[ i ], x2d[ i ], x3d[ i ], y1d[ i ],
      y2d[ i ], y3d[ i ], n[ 0 ], n[ 1 ], n[ 2 ], kernelRe, kernelIm );

    entryRe += wxyd[ i ] * kernelRe;
    entryIm += wxyd[ i ] * kernelIm;
  }

  SCVT outerArea = this->space->getLeftMesh( )->getElemArea( outerElem );
  SCVT innerArea = this->space->getRightMesh( )->getElemArea( innerElem );
  SCVT mult = outerArea * innerArea;

  entryRe *= mult;
  entryIm *= mult;
  elemMatrix.set( 0, 0, SC( entryRe, entryIm ) );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix2LayerDisjointP1P1(
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
  SCVT entry11Re = 0.0;
  SCVT entry21Re = 0.0;
  SCVT entry31Re = 0.0;
  SCVT entry12Re = 0.0;
  SCVT entry22Re = 0.0;
  SCVT entry32Re = 0.0;
  SCVT entry13Re = 0.0;
  SCVT entry23Re = 0.0;
  SCVT entry33Re = 0.0;
  SCVT entry11Im = 0.0;
  SCVT entry21Im = 0.0;
  SCVT entry31Im = 0.0;
  SCVT entry12Im = 0.0;
  SCVT entry22Im = 0.0;
  SCVT entry32Im = 0.0;
  SCVT entry13Im = 0.0;
  SCVT entry23Im = 0.0;
  SCVT entry33Im = 0.0;

  SCVT kernelRe, kernelIm;

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
reduction( + : entry11Re, entry11Im ) \
reduction( + : entry12Re, entry12Im ) \
reduction( + : entry13Re, entry13Im ) \
reduction( + : entry21Re, entry21Im ) \
reduction( + : entry22Re, entry22Im ) \
reduction( + : entry23Re, entry23Im ) \
reduction( + : entry31Re, entry31Im ) \
reduction( + : entry32Re, entry32Im ) \
reduction( + : entry33Re, entry33Im ) \
aligned( x1d, x2d, x3d, x1dref, x2dref : align ) \
aligned( y1d, y2d, y3d, y1dref, y2dref : align ) \
aligned( wxyd : align ) \
private( kernelRe, kernelIm, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y ) \
simdlen( width )
  for ( i = 0; i < nIters; ++i ) {

    this->evalDoubleLayerKernel( x1d[ i ], x2d[ i ], x3d[ i ], y1d[ i ],
      y2d[ i ], y3d[ i ], n[ 0 ], n[ 1 ], n[ 2 ], kernelRe, kernelIm );

    kernelRe *= wxyd[ i ];
    kernelIm *= wxyd[ i ];

    phi1y = (SCVT) 1.0 - y1dref[ i ] - y2dref[ i ];
    phi2y = y1dref[ i ];
    phi3y = y2dref[ i ];
    phi1x = (SCVT) 1.0 - x1dref[ i ] - x2dref[ i ];
    phi2x = x1dref[ i ];
    phi3x = x2dref[ i ];

    entry11Re += kernelRe * phi1x * phi1y;
    entry21Re += kernelRe * phi2x * phi1y;
    entry31Re += kernelRe * phi3x * phi1y;
    entry12Re += kernelRe * phi1x * phi2y;
    entry22Re += kernelRe * phi2x * phi2y;
    entry32Re += kernelRe * phi3x * phi2y;
    entry13Re += kernelRe * phi1x * phi3y;
    entry23Re += kernelRe * phi2x * phi3y;
    entry33Re += kernelRe * phi3x * phi3y;
    entry11Im += kernelIm * phi1x * phi1y;
    entry21Im += kernelIm * phi2x * phi1y;
    entry31Im += kernelIm * phi3x * phi1y;
    entry12Im += kernelIm * phi1x * phi2y;
    entry22Im += kernelIm * phi2x * phi2y;
    entry32Im += kernelIm * phi3x * phi2y;
    entry13Im += kernelIm * phi1x * phi3y;
    entry23Im += kernelIm * phi2x * phi3y;
    entry33Im += kernelIm * phi3x * phi3y;
  }

  SCVT areasM = this->space->getLeftMesh( )->getElemArea( outerElem ) *
    this->space->getRightMesh( )->getElemArea( innerElem );

  elemMatrix.set( 0, 0, SC( entry11Re, entry11Im ) );
  elemMatrix.set( 1, 0, SC( entry21Re, entry21Im ) );
  elemMatrix.set( 2, 0, SC( entry31Re, entry31Im ) );
  elemMatrix.set( 0, 1, SC( entry12Re, entry12Im ) );
  elemMatrix.set( 1, 1, SC( entry22Re, entry22Im ) );
  elemMatrix.set( 2, 1, SC( entry32Re, entry32Im ) );
  elemMatrix.set( 0, 2, SC( entry13Re, entry13Im ) );
  elemMatrix.set( 1, 2, SC( entry23Re, entry23Im ) );
  elemMatrix.set( 2, 2, SC( entry33Re, entry33Im ) );
  elemMatrix.scale( areasM );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::
computeElemMatrixHypersingularSauterSchwabP1P1(
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
  SCVT on[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT in[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->getSpace( )->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->getSpace( )->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->getSpace( )->getLeftMesh( )->getNormal( outerElem, on );
  this->getSpace( )->getRightMesh( )->getNormal( innerElem, in );

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

  SCVT entry11Re = 0.0;
  SCVT entry21Re = 0.0;
  SCVT entry31Re = 0.0;
  SCVT entry12Re = 0.0;
  SCVT entry22Re = 0.0;
  SCVT entry32Re = 0.0;
  SCVT entry13Re = 0.0;
  SCVT entry23Re = 0.0;
  SCVT entry33Re = 0.0;
  SCVT entry11Im = 0.0;
  SCVT entry21Im = 0.0;
  SCVT entry31Im = 0.0;
  SCVT entry12Im = 0.0;
  SCVT entry22Im = 0.0;
  SCVT entry32Im = 0.0;
  SCVT entry13Im = 0.0;
  SCVT entry23Im = 0.0;
  SCVT entry33Im = 0.0;

  SCVT kappaRe = std::real( this->kappa );
  SCVT kappaIm = std::imag( this->kappa );
  SCVT ndot = DOT3( on, in );
  SCVT ndotRe = ndot * ( kappaRe * kappaRe - kappaIm * kappaIm );
  SCVT ndotIm = ndot * (SCVT) 2.0 * kappaRe * kappaIm;
  SCVT kernelRe, kernelIm;
  //SCVT aux1[ 9 ] __attribute__( ( aligned( DATA_ALIGN ) ) );
  //SCVT aux2[ 9 ] __attribute__( ( aligned( DATA_ALIGN ) ) );
  SCVT aux10, aux11, aux12, aux13, aux14, aux15, aux16, aux17, aux18;
  SCVT aux20, aux21, aux22, aux23, aux24, aux25, aux26, aux27, aux28;
  SCVT mult, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y;
  SCVT curlMultRe = 0.0;
  SCVT curlMultIm = 0.0;
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
reduction( + : entry11Re, entry21Re, entry31Re ) \
reduction( + : entry12Re, entry22Re, entry32Re ) \
reduction( + : entry13Re, entry23Re, entry33Re ) \
reduction( + : entry11Im, entry21Im, entry31Im ) \
reduction( + : entry12Im, entry22Im, entry32Im ) \
reduction( + : entry13Im, entry23Im, entry33Im ) \
reduction( + : curlMultRe, curlMultIm ) \
aligned( x1ss, x2ss, x3ss, y1ss, y2ss, y3ss, jacobian : align ) \
aligned( x1ref, x2ref, y1ref, y2ref : align ) \
private( kernelRe, kernelIm, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y ) \
private( aux10, aux11, aux12, aux13, aux14, aux15, aux16, aux17, aux18 ) \
private( aux20, aux21, aux22, aux23, aux24, aux25, aux26, aux27, aux28 ) \
simdlen( width )
    for ( i = 0; i < totalSize; ++i ) {
      this->evalSingleLayerKernel( x1ss[ i ], x2ss[ i ], x3ss[ i ], y1ss[ i ],
        y2ss[ i ], y3ss[ i ], kernelRe, kernelIm );
      kernelRe *= jacobian[ i ];
      kernelIm *= jacobian[ i ];
      curlMultRe += kernelRe;
      curlMultIm += kernelIm;

      phi1x = (SCVT) 1.0 - x1ref[ i ];
      phi2x = x1ref[ i ] - x2ref[ i ];
      phi3x = x2ref[ i ];
      phi1y = (SCVT) 1.0 - y1ref[ i ];
      phi2y = y1ref[ i ] - y2ref[ i ];
      phi3y = y2ref[ i ];

      aux10 = -ndotRe * phi1x * phi1y;
      aux11 = -ndotRe * phi2x * phi1y;
      aux12 = -ndotRe * phi3x * phi1y;
      aux13 = -ndotRe * phi1x * phi2y;
      aux14 = -ndotRe * phi2x * phi2y;
      aux15 = -ndotRe * phi3x * phi2y;
      aux16 = -ndotRe * phi1x * phi3y;
      aux17 = -ndotRe * phi2x * phi3y;
      aux18 = -ndotRe * phi3x * phi3y;

      aux20 = ndotIm * phi1x * phi1y;
      aux21 = ndotIm * phi2x * phi1y;
      aux22 = ndotIm * phi3x * phi1y;
      aux23 = ndotIm * phi1x * phi2y;
      aux24 = ndotIm * phi2x * phi2y;
      aux25 = ndotIm * phi3x * phi2y;
      aux26 = ndotIm * phi1x * phi3y;
      aux27 = ndotIm * phi2x * phi3y;
      aux28 = ndotIm * phi3x * phi3y;

      entry11Re += kernelRe * aux10 + kernelIm * aux20;
      entry21Re += kernelRe * aux11 + kernelIm * aux21;
      entry31Re += kernelRe * aux12 + kernelIm * aux22;
      entry12Re += kernelRe * aux13 + kernelIm * aux23;
      entry22Re += kernelRe * aux14 + kernelIm * aux24;
      entry32Re += kernelRe * aux15 + kernelIm * aux25;
      entry13Re += kernelRe * aux16 + kernelIm * aux26;
      entry23Re += kernelRe * aux17 + kernelIm * aux27;
      entry33Re += kernelRe * aux18 + kernelIm * aux28;
      entry11Im += kernelIm * aux10 - kernelRe * aux20;
      entry21Im += kernelIm * aux11 - kernelRe * aux21;
      entry31Im += kernelIm * aux12 - kernelRe * aux22;
      entry12Im += kernelIm * aux13 - kernelRe * aux23;
      entry22Im += kernelIm * aux14 - kernelRe * aux24;
      entry32Im += kernelIm * aux15 - kernelRe * aux25;
      entry13Im += kernelIm * aux16 - kernelRe * aux26;
      entry23Im += kernelIm * aux17 - kernelRe * aux27;
      entry33Im += kernelIm * aux18 - kernelRe * aux28;
    }
  }

  // first two x indices swapped for common edge!
  if ( type == 1 ) {
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ),
      SC( entry11Re, entry11Im ) );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ),
      SC( entry21Re, entry21Im ) );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ),
      SC( entry31Re, entry31Im ) );
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ),
      SC( entry12Re, entry12Im ) );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ),
      SC( entry22Re, entry22Im ) );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ),
      SC( entry32Re, entry32Im ) );
  } else {
    matrix.set( this->mod3( outerRot ), this->mod3( innerRot ),
      SC( entry11Re, entry11Im ) );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( innerRot ),
      SC( entry21Re, entry21Im ) );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( innerRot ),
      SC( entry31Re, entry31Im ) );
    matrix.set( this->mod3( outerRot ), this->mod3( 1 + innerRot ),
      SC( entry12Re, entry12Im ) );
    matrix.set( this->mod3( 1 + outerRot ), this->mod3( 1 + innerRot ),
      SC( entry22Re, entry22Im ) );
    matrix.set( this->mod3( 2 + outerRot ), this->mod3( 1 + innerRot ),
      SC( entry32Re, entry32Im ) );
  }
  matrix.set( this->mod3( outerRot ), this->mod3( 2 + innerRot ),
    SC( entry13Re, entry13Im ) );
  matrix.set( this->mod3( 1 + outerRot ), this->mod3( 2 + innerRot ),
    SC( entry23Re, entry23Im ) );
  matrix.set( this->mod3( 2 + outerRot ), this->mod3( 2 + innerRot ),
    SC( entry33Re, entry33Im ) );

  SC * matrixData = matrix.getData( );
  SCVT curls[ 9 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT * ocurls =
    this->space->getLeftMesh( )->getCurls( )->getData( ) + 9 * outerElem;
  SCVT * icurls =
    this->space->getRightMesh( )->getCurls( )->getData( ) + 9 * innerElem;
  for ( int oRot = 0; oRot < 3; ++oRot ) {
    for ( int iRot = 0; iRot < 3; ++iRot ) {
      curls[ oRot + 3 * iRot ] =
        DOT3( ( ocurls + 3 * oRot ), ( icurls + 3 * iRot ) );
    }
  }

  SCVT areaMult = (SCVT) 4.0 *
    this->getSpace( )->getRightMesh( )->getElemArea( innerElem ) *
    this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );

  for ( i = 0; i < 9; ++i ) {
    matrixData[ i ] += curls[ i ] * SC( curlMultRe, curlMultIm );
    matrixData[ i ] *= areaMult;
  }
}
//*/

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrixHypersingularDisjointP1P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  SCVT * ocurls =
    this->space->getLeftMesh( )->getCurls( )->getData( ) + 9 * outerElem;
  SCVT * icurls =
    this->space->getRightMesh( )->getCurls( )->getData( ) + 9 * innerElem;

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT on[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT in[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT curls[ 9 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  for ( int oRot = 0; oRot < 3; ++oRot ) {
    for ( int iRot = 0; iRot < 3; ++iRot ) {
      curls[ oRot + 3 * iRot ] =
        DOT3( ( ocurls + 3 * oRot ), ( icurls + 3 * iRot ) );
    }
  }

  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->space->getLeftMesh( )->getNormal( outerElem, on );
  this->space->getRightMesh( )->getNormal( innerElem, in );

  this->updateDisjointQuadratureNodes( x1, x2, x3, y1, y2, y3 );

  int outerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 0 ] ];
  int innerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 1 ] ];

  LO nIters = outerPoints * innerPoints;
  LO i = 0;
  SCVT entry11Re = 0.0;
  SCVT entry21Re = 0.0;
  SCVT entry31Re = 0.0;
  SCVT entry12Re = 0.0;
  SCVT entry22Re = 0.0;
  SCVT entry32Re = 0.0;
  SCVT entry13Re = 0.0;
  SCVT entry23Re = 0.0;
  SCVT entry33Re = 0.0;
  SCVT entry11Im = 0.0;
  SCVT entry21Im = 0.0;
  SCVT entry31Im = 0.0;
  SCVT entry12Im = 0.0;
  SCVT entry22Im = 0.0;
  SCVT entry32Im = 0.0;
  SCVT entry13Im = 0.0;
  SCVT entry23Im = 0.0;
  SCVT entry33Im = 0.0;

  SCVT phi1y, phi2y, phi3y, phi1x, phi2x, phi3x;

  SCVT kappaRe = std::real( this->kappa );
  SCVT kappaIm = std::imag( this->kappa );
  SCVT ndot = DOT3( on, in );
  SCVT ndotRe = ndot * ( kappaRe * kappaRe - kappaIm * kappaIm );
  SCVT ndotIm = ndot * (SCVT) 2.0 * kappaRe * kappaIm;
  SCVT kernelRe, kernelIm;
  //SCVT aux1[ 9 ] __attribute__( ( aligned( DATA_ALIGN ) ) );
  //SCVT aux2[ 9 ] __attribute__( ( aligned( DATA_ALIGN ) ) );
  SCVT aux10, aux11, aux12, aux13, aux14, aux15, aux16, aux17, aux18;
  SCVT aux20, aux21, aux22, aux23, aux24, aux25, aux26, aux27, aux28;

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
reduction( + : entry11Re, entry21Re, entry31Re ) \
reduction( + : entry12Re, entry22Re, entry32Re ) \
reduction( + : entry13Re, entry23Re, entry33Re ) \
reduction( + : entry11Im, entry21Im, entry31Im ) \
reduction( + : entry12Im, entry22Im, entry32Im ) \
reduction( + : entry13Im, entry23Im, entry33Im ) \
aligned( x1d, x2d, x3d, x1dref, x2dref : align ) \
aligned( y1d, y2d, y3d, y1dref, y2dref : align ) \
aligned( wxyd : align ) \
private( kernelRe, kernelIm, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y ) \
private( aux10, aux11, aux12, aux13, aux14, aux15, aux16, aux17, aux18 ) \
private( aux20, aux21, aux22, aux23, aux24, aux25, aux26, aux27, aux28 ) \
simdlen( width )
  for ( i = 0; i < nIters; ++i ) {

    this->evalSingleLayerKernel( x1d[ i ], x2d[ i ], x3d[ i ],
      y1d[ i ], y2d[ i ], y3d[ i ], kernelRe, kernelIm );

    kernelRe *= wxyd[ i ];
    kernelIm *= wxyd[ i ];

    phi1y = (SCVT) 1.0 - y1dref[ i ] - y2dref[ i ];
    phi2y = y1dref[ i ];
    phi3y = y2dref[ i ];
    phi1x = (SCVT) 1.0 - x1dref[ i ] - x2dref[ i ];
    phi2x = x1dref[ i ];
    phi3x = x2dref[ i ];

    aux10 = curls[ 0 ] - ndotRe * phi1x * phi1y;
    aux11 = curls[ 1 ] - ndotRe * phi2x * phi1y;
    aux12 = curls[ 2 ] - ndotRe * phi3x * phi1y;
    aux13 = curls[ 3 ] - ndotRe * phi1x * phi2y;
    aux14 = curls[ 4 ] - ndotRe * phi2x * phi2y;
    aux15 = curls[ 5 ] - ndotRe * phi3x * phi2y;
    aux16 = curls[ 6 ] - ndotRe * phi1x * phi3y;
    aux17 = curls[ 7 ] - ndotRe * phi2x * phi3y;
    aux18 = curls[ 8 ] - ndotRe * phi3x * phi3y;

    aux20 = ndotIm * phi1x * phi1y;
    aux21 = ndotIm * phi2x * phi1y;
    aux22 = ndotIm * phi3x * phi1y;
    aux23 = ndotIm * phi1x * phi2y;
    aux24 = ndotIm * phi2x * phi2y;
    aux25 = ndotIm * phi3x * phi2y;
    aux26 = ndotIm * phi1x * phi3y;
    aux27 = ndotIm * phi2x * phi3y;
    aux28 = ndotIm * phi3x * phi3y;

    entry11Re += kernelRe * aux10 + kernelIm * aux20;
    entry21Re += kernelRe * aux11 + kernelIm * aux21;
    entry31Re += kernelRe * aux12 + kernelIm * aux22;
    entry12Re += kernelRe * aux13 + kernelIm * aux23;
    entry22Re += kernelRe * aux14 + kernelIm * aux24;
    entry32Re += kernelRe * aux15 + kernelIm * aux25;
    entry13Re += kernelRe * aux16 + kernelIm * aux26;
    entry23Re += kernelRe * aux17 + kernelIm * aux27;
    entry33Re += kernelRe * aux18 + kernelIm * aux28;
    entry11Im += kernelIm * aux10 - kernelRe * aux20;
    entry21Im += kernelIm * aux11 - kernelRe * aux21;
    entry31Im += kernelIm * aux12 - kernelRe * aux22;
    entry12Im += kernelIm * aux13 - kernelRe * aux23;
    entry22Im += kernelIm * aux14 - kernelRe * aux24;
    entry32Im += kernelIm * aux15 - kernelRe * aux25;
    entry13Im += kernelIm * aux16 - kernelRe * aux26;
    entry23Im += kernelIm * aux17 - kernelRe * aux27;
    entry33Im += kernelIm * aux18 - kernelRe * aux28;
  }

  SCVT areasM = this->space->getLeftMesh( )->getElemArea( outerElem ) *
    this->space->getRightMesh( )->getElemArea( innerElem );

  matrix.set( 0, 0, SC( entry11Re, entry11Im ) );
  matrix.set( 1, 0, SC( entry21Re, entry21Im ) );
  matrix.set( 2, 0, SC( entry31Re, entry31Im ) );
  matrix.set( 0, 1, SC( entry12Re, entry12Im ) );
  matrix.set( 1, 1, SC( entry22Re, entry22Im ) );
  matrix.set( 2, 1, SC( entry32Re, entry32Im ) );
  matrix.set( 0, 2, SC( entry13Re, entry13Im ) );
  matrix.set( 1, 2, SC( entry23Re, entry23Im ) );
  matrix.set( 2, 2, SC( entry33Re, entry33Im ) );
  matrix.scale( areasM );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::representationFormulaP1P1(
  const SCVT *xCoord,
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

  SC * shiftedKernel = new SC[ quadSizes[ this->quadratureOrder[ 1 ] ] ];
  SC * shiftedKernelNeumann =
    new SC[ quadSizes[ this->quadratureOrder[ 1 ] ] ];
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
      //if ( false ) {
      if ( DIST3( centroid, xSingle ) < 2.0 * diam ) {

        for ( int k = 0; k < quadSizes[ this->quadratureOrder[ 1 ] ]; ++k ) {
          shiftedKernel[ k ] = evalShiftedKernel( xSingle,
            ( quadratureNodes + 3 * k ) );
          shiftedKernelNeumann[ k ] = evalShiftedKernelNeumann( xSingle,
            ( quadratureNodes + 3 * k ), n );
        }

        for ( int rot = 0; rot < 3; ++rot ) {

          this->getLocalCoordinates( x1, x2, x3, r1, r2, n, stau, alpha1,
            alpha2, rot );
          this->getLocalQuadratureNode( x1, x2, x3, xSingle, r1, r2, n,
            xLocal, rot );

          val += neu.get( ind[ rot ] ) *
            collocation1LayerP1( xLocal, stau, alpha1, alpha2,
            this->quadratureOrder[ 1 ], shiftedKernel, j, rot );

          val -= dir.get( ind[ rot ] ) *
            collocation2LayerP1( xLocal, stau, alpha1, alpha2,
            this->quadratureOrder[ 1 ], shiftedKernelNeumann, j, rot );
        }

      } else {

        onePatch = 0.0;

        for ( int k = 0; k < quadSizes[ this->quadratureOrder[ 1 ] ]; ++k ) {

          intPx = qPointsIn[ 2 * k ];
          intPy = qPointsIn[ 2 * k + 1 ];

          slkernel = evalSingleLayerKernel( xSingle, ( quadratureNodes + 3 * k ) );
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

  delete [] shiftedKernel;
  delete [] shiftedKernelNeumann;
  delete [] quadratureNodes;
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::representationFormulaP0P0(
  const SCVT *xCoord,
  LO numPoints,
  const Vector<LO, SC> & dir,
  const Vector<LO, SC> & neu,
  bool interior,
  Vector<LO, SC> & values
  ) const {

  LO nElems = this->space->getMesh( )->getNElements( );
  //SCVT stau, alpha1, alpha2;
  SCVT x1[3], x2[3], x3[3];
  //SCVT r1[3], r2[3], n[3], xLocal[3];
  //SCVT centroid[3];
  SCVT n[3];
  //SCVT r[3] = { 0.0, 0.0, 0.0 };
  //SCVT diam = 0.0;
  //LO ind[3];
  SC val;
  //SC dlkernel;
  SC integrand;
  SC onePatch;

  //SC * shiftedKernel = new SC[ quadSizes[ this->quadratureOrder[ 1 ] ] ];
  //SC * shiftedKernelNeumann =
  //    new SC[ quadSizes[ this->quadratureOrder[ 1 ] ] ];
  SCVT * quadratureNodes =
    new SCVT[ quadSizes[ this->quadratureOrder[ 1 ] ] * 3 ];

  // points in reference triangle
  double * qPointsIn = quadPoints[ this->quadratureOrder[ 1 ] ];
  SCVT intPx, intPy;

  for ( LO i = 0; i < numPoints; ++i ) {

    val = 0.0;
    const SCVT * xSingle = xCoord + 3 * i;

    for ( LO j = 0; j < nElems; ++j ) {

      //      this->space->getMesh( )->getCentroid( j, centroid );
      this->space->getMesh( )->getNodes( j, x1, x2, x3 );
      //      r[0] = DIST3( centroid, x1 );
      //      r[1] = DIST3( centroid, x2 );
      //      r[2] = DIST3( centroid, x3 );
      //
      //      diam = 2.0 * ( *std::max_element( r, r + 3 ) );

      //     this->space->getMesh( )->getElement( j, ind );
      this->space->getMesh( )->getNormal( j, n );
      this->getQuadratureNodes( x1, x2, x3, this->quadratureOrder[ 1 ],
        quadratureNodes );

      //      if ( DIST3( centroid, xSingle ) < 2.0 * diam ) {
      //
      //        for ( int k = 0; k < quadSizes[ this->quadratureOrder[ 1 ] ]; ++k ) {
      //          shiftedKernel[ k ] = evalShiftedKernel( xSingle,
      //              ( quadratureNodes + 3 * k ) );
      //          shiftedKernelNeumann[ k ] = evalShiftedKernelNeumann( xSingle,
      //              ( quadratureNodes + 3 * k ), n );
      //        }
      //
      //        for ( int rot = 0; rot < 3; ++rot ) {
      //
      //          this->getLocalCoordinates( x1, x2, x3, r1, r2, n, stau, alpha1,
      //              alpha2, rot );
      //          this->getLocalQuadratureNode( x1, x2, x3, xSingle, r1, r2, n,
      //              xLocal, rot );
      //
      //          if ( rot == 0 ) {
      //            val += neu.get( j ) *
      //                collocation1LayerP0( xLocal, stau, alpha1, alpha2,
      //                this->quadratureOrder[ 1 ], shiftedKernel, j );
      //          }
      //
      //          val -= dir.get( ind[ rot ] ) *
      //              collocation2LayerP1( xLocal, stau, alpha1, alpha2,
      //              this->quadratureOrder[ 1 ], shiftedKernelNeumann, j, rot );
      //        }
      //
      //      } else {

      onePatch = 0.0;

      for ( int k = 0; k < quadSizes[ this->quadratureOrder[ 1 ] ]; ++k ) {

        intPx = qPointsIn[ 2 * k ];
        intPy = qPointsIn[ 2 * k + 1 ];

        integrand = neu.get( j ) *
          evalSingleLayerKernel( xSingle, ( quadratureNodes + 3 * k ) );

        integrand -= dir.get( j ) * evalDoubleLayerKernel( xSingle,
          ( quadratureNodes + 3 * k ), n );

        onePatch += integrand *
          (SCVT) quadWeights[ this->quadratureOrder[ 1 ] ][ k ];
      }

      val += onePatch * this->space->getMesh( )->getElemArea( j );
      //}
    }

    values.set( i, val );
  }

  if ( !interior ) {
    values.scale( -1.0 );
  }

  //delete [] shiftedKernel;
  //delete [] shiftedKernelNeumann;
  delete [] quadratureNodes;
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrixH1SauterSchwabP0P0(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC>& matrix
  ) const {
  /*
  if ( this->quadratureOrderDisjointElems &&
    this->areElementsDisjoint( outerElem, innerElem ) ) {
    SCVT x1[3], x2[3], x3[3], y1[3], y2[3], y3[3];

    int qOrderOut = this->quadratureOrderDisjointElems[0];
    int qOrderIn = this->quadratureOrderDisjointElems[1];

    // getting outer quadrature points
    this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
    this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
    this->getQuadratureNodes( x1, x2, x3, qOrderOut,
      this->outerQuadratureNodesDisjoint );
    this->getQuadratureNodes( y1, y2, y3, qOrderIn,
      this->innerQuadratureNodesDisjoint );

    SC entry = 0.0;
    SCVT nx[3], ny[3];
    this->getSpace( )->getLeftMesh( )->getNormal( outerElem, nx );
    this->getSpace( )->getRightMesh( )->getNormal( innerElem, ny );
    SCVT dot = DOT3( nx, ny );

    SCVT innerArea = this->getSpace( )->getRightMesh( )->getElemArea( innerElem );
    SCVT outerArea = this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );

    int outerPoints = quadSizes[ qOrderOut ];
    int innerPoints = quadSizes[ qOrderIn ];

    for ( int i = 0; i < outerPoints; i++ ) {
      for ( int j = 0; j < innerPoints; j++ ) {
        entry += (SCVT) quadWeights[ qOrderOut ][ i ] *
          (SCVT) quadWeights[ qOrderIn ][ j ] *
          evalSingleLayerKernel(
          ( this->outerQuadratureNodesDisjoint + 3 * i ),
          ( this->innerQuadratureNodesDisjoint + 3 * j ) );
      }
    }

    matrix.set( 0, 0, kappa * kappa * dot * entry * innerArea * outerArea );
    return;
  }

  int outerRot;
  int innerRot;
  int type;

  this->getSauterSchwabQuadratureType( innerElem, outerElem, type, innerRot,
    outerRot );

  int qOrder0 = this->quadratureOrder[ 0 ];
  int qOrder1 = this->quadratureOrder[ 1 ];
  int qOrder2 = this->quadratureOrder[ 2 ];
  int qOrder3 = this->quadratureOrder[ 3 ];

  int qSize0 = lineQuadSizes[ qOrder0 ];
  int qSize1 = lineQuadSizes[ qOrder1 ];
  int qSize2 = lineQuadSizes[ qOrder2 ];
  int qSize3 = lineQuadSizes[ qOrder3 ];

  SC entry = 0.0;
  SCVT xRef[2], yRef[2];
  SCVT x[3], y[3];
  SCVT w[4];
  SCVT x1[3], x2[3], x3[3];
  SCVT y1[3], y2[3], y3[3];
  SCVT weights[4];
  SCVT jacobian;
  SCVT nx[3], ny[3];

  this->getSpace( )->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  this->getSpace( )->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );

  this->getSpace( )->getLeftMesh( )->getNormal( outerElem, nx );
  this->getSpace( )->getRightMesh( )->getNormal( innerElem, ny );
  SCVT dot = DOT3( nx, ny );

  SCVT innerArea = this->getSpace( )->getRightMesh( )->getElemArea( innerElem );
  SCVT outerArea = this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );

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
                this->cube2triIdentical( w, simplex, xRef, yRef, jacobian );
                this->tri2element( x1, x2, x3, xRef, x );
                this->tri2element( y1, y2, y3, yRef, y );
                entry += evalSingleLayerKernel( x, y ) *
                  weights[0] * weights[1] * weights[2] * weights[3] * jacobian;
              }
              break;

              // common edge
            case 1:
              for ( int simplex = 0; simplex < 5; simplex++ ) {
                this->cube2triCommonEdge( w, simplex, xRef, yRef, jacobian );
                // outer element swaps first two nodes, so that the edges agree
                this->tri2element( x1, x2, x3, xRef, outerRot, x, true );
                this->tri2element( y1, y2, y3, yRef, innerRot, y );

                entry += jacobian * evalSingleLayerKernel( x, y ) *
                  weights[0] * weights[1] * weights[2] * weights[3];
              }
              break;

              // common vertex
            case 2:
              for ( int simplex = 0; simplex < 2; simplex++ ) {
                this->cube2triCommonVertex( w, simplex, xRef, yRef, jacobian );
                this->tri2element( x1, x2, x3, xRef, outerRot, x );
                this->tri2element( y1, y2, y3, yRef, innerRot, y );
                entry += evalSingleLayerKernel( x, y ) *
                  weights[0] * weights[1] * weights[2] * weights[3] *
                  jacobian;
              }
              break;

              // disjoint triangles  
            case 3:
              this->cube2triDisjoint( w, xRef, yRef, jacobian );
              this->tri2element( x1, x2, x3, xRef, x );
              this->tri2element( y1, y2, y3, yRef, y );
              entry += evalSingleLayerKernel( x, y ) *
                weights[0] * weights[1] * weights[2] * weights[3] * jacobian;
              break;
          }
        }
      }
    }
  }

  matrix.set( 0, 0, kappa * kappa * dot * entry * ( (SCVT) 4.0 ) * innerArea * outerArea );
   */
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrixH2SauterSchwabP0P0(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC>& matrix
  ) const {
  /*
    if ( this->quadratureOrderDisjointElems &&
      this->areElementsDisjoint( outerElem, innerElem ) ) {
      SCVT x1[3], x2[3], x3[3], y1[3], y2[3], y3[3];

      int qOrderOut = this->quadratureOrderDisjointElems[0];
      int qOrderIn = this->quadratureOrderDisjointElems[1];

      // getting outer quadrature points
      this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
      this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
      this->getQuadratureNodes( x1, x2, x3, qOrderOut,
        this->outerQuadratureNodesDisjoint );
      this->getQuadratureNodes( y1, y2, y3, qOrderIn,
        this->innerQuadratureNodesDisjoint );

      SC entry = 0.0;
      SCVT nx[3], ny[3];
      this->getSpace( )->getLeftMesh( )->getNormal( outerElem, nx );
      this->getSpace( )->getRightMesh( )->getNormal( innerElem, ny );


      int outerPoints = quadSizes[ qOrderOut ];
      int innerPoints = quadSizes[ qOrderIn ];

      for ( int i = 0; i < outerPoints; i++ ) {
        for ( int j = 0; j < innerPoints; j++ ) {
          entry += (SCVT) quadWeights[ qOrderOut ][ i ] *
            (SCVT) quadWeights[ qOrderIn ][ j ] *
            evalHypersingularP0P0Kernel(
            ( this->outerQuadratureNodesDisjoint + 3 * i ),
            ( this->innerQuadratureNodesDisjoint + 3 * j ), nx, ny );
        }
      }

      entry *= this->space->getLeftMesh( )->getElemArea( outerElem ) *
        this->space->getRightMesh( )->getElemArea( innerElem );
      matrix.set( 0, 0, entry );
      return;
    }

    int outerRot;
    int innerRot;
    int type;

    this->getSauterSchwabQuadratureType( innerElem, outerElem, type, innerRot,
      outerRot );

    int qOrder0 = this->quadratureOrder[ 0 ];
    int qOrder1 = this->quadratureOrder[ 1 ];
    int qOrder2 = this->quadratureOrder[ 2 ];
    int qOrder3 = this->quadratureOrder[ 3 ];

    int qSize0 = lineQuadSizes[ qOrder0 ];
    int qSize1 = lineQuadSizes[ qOrder1 ];
    int qSize2 = lineQuadSizes[ qOrder2 ];
    int qSize3 = lineQuadSizes[ qOrder3 ];

    SC entry = 0.0;
    SCVT xRef[2], yRef[2];
    SCVT x[3], y[3];
    SCVT w[4];
    SCVT x1[3], x2[3], x3[3];
    SCVT y1[3], y2[3], y3[3];
    SCVT weights[4];
    SCVT jacobian;
    SCVT nx[3], ny[3];

    this->getSpace( )->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
    this->getSpace( )->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );

    this->getSpace( )->getLeftMesh( )->getNormal( outerElem, nx );
    this->getSpace( )->getRightMesh( )->getNormal( innerElem, ny );

    SCVT innerArea = this->getSpace( )->getRightMesh( )->getElemArea( innerElem );
    SCVT outerArea = this->getSpace( )->getLeftMesh( )->getElemArea( outerElem );

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
                  this->cube2triIdentical( w, simplex, xRef, yRef, jacobian );
                  this->tri2element( x1, x2, x3, xRef, x );
                  this->tri2element( y1, y2, y3, yRef, y );
                  entry += evalHypersingularP0P0Kernel( x, y, nx, ny ) *
                    weights[0] * weights[1] * weights[2] * weights[3] * jacobian;
                }
                break;

                // common edge
              case 1:
                for ( int simplex = 0; simplex < 5; simplex++ ) {
                  this->cube2triCommonEdge( w, simplex, xRef, yRef, jacobian );
                  // outer element swaps first two nodes, so that the edges agree
                  this->tri2element( x1, x2, x3, xRef, outerRot, x, true );
                  this->tri2element( y1, y2, y3, yRef, innerRot, y );

                  entry += jacobian * evalHypersingularP0P0Kernel( x, y, nx, ny ) *
                    weights[0] * weights[1] * weights[2] * weights[3];
                }
                break;

                // common vertex
              case 2:
                for ( int simplex = 0; simplex < 2; simplex++ ) {
                  this->cube2triCommonVertex( w, simplex, xRef, yRef, jacobian );
                  this->tri2element( x1, x2, x3, xRef, outerRot, x );
                  this->tri2element( y1, y2, y3, yRef, innerRot, y );
                  entry += evalHypersingularP0P0Kernel( x, y, nx, ny ) *
                    weights[0] * weights[1] * weights[2] * weights[3] *
                    jacobian;
                }
                break;

                // disjoint triangles  
              case 3:
                this->cube2triDisjoint( w, xRef, yRef, jacobian );
                this->tri2element( x1, x2, x3, xRef, x );
                this->tri2element( y1, y2, y3, yRef, y );
                entry += evalHypersingularP0P0Kernel( x, y, nx, ny ) *
                  weights[0] * weights[1] * weights[2] * weights[3] * jacobian;
                break;
            }
          }
        }
      }
    }

    matrix.set( 0, 0, entry * ( (SCVT) 4.0 ) * innerArea * outerArea );
   */
}

template<class LO, class SC>
BEIntegratorHelmholtz<LO, SC>::BEIntegratorHelmholtz( ) {
}

template<class LO, class SC>
BEIntegratorHelmholtz<LO, SC>::BEIntegratorHelmholtz(
  const BEIntegratorHelmholtz& orig
  ) {

  this->space = orig.space;


  // TODO copy the rest of the variables and arrays
}

template<class LO, class SC>
BEIntegratorHelmholtz<LO, SC>::BEIntegratorHelmholtz(
  BESpace<LO, SC>* space,
  int *quadratureOrder,
  SC kappa,
  quadratureType quadrature,
  int* quadratureOrderDisjointElems
  ) {

  this->space = space;
  this->quadratureOrder = quadratureOrder;
  this->quadrature = quadrature;
  this->kappa = kappa;
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
BEIntegratorHelmholtz<LO, SC>::~BEIntegratorHelmholtz( ) {
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix1LayerP0P0(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  SCVT stau, alpha1, alpha2;

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

  // getting collapsed quadrature data for regular part
  this->updateSteinbachRegularQuadratureNodes( x1, x2, x3, y1, y2, y3 );

  SCVT * sx = this->sx;
  SCVT * tx = this->tx;
  SCVT * ux = this->ux;
  SCVT * wxst = this->wxst;

  SCVT entrySing = 0.0;

  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;
  int i;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entrySing ) \
aligned( sx, tx, ux, wxst : align ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXpad; ++i ) {
    entrySing += wxst[ i ] * this->collocationSingular1LayerP0( sx[ i ],
      tx[ i ], ux[ i ], stau, alpha1, alpha2 );
  }

  SCVT * x1stCol = this->x1stCol;
  SCVT * x2stCol = this->x2stCol;
  SCVT * x3stCol = this->x3stCol;
  SCVT * y1stCol = this->y1stCol;
  SCVT * y2stCol = this->y2stCol;
  SCVT * y3stCol = this->y3stCol;
  SCVT * wxyst = this->wxyst;

  SCVT entryRe = 0.0;
  SCVT entryIm = 0.0;
  SCVT kernelRe;
  SCVT kernelIm;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entryRe, entryIm ) \
aligned( x1stCol, x2stCol, x3stCol, wxyst : align ) \
aligned( y1stCol, y2stCol, y3stCol : align ) \
private( kernelRe, kernelIm ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXYpad; ++i ) {
    this->evalShiftedKernel( x1stCol[ i ], x2stCol[ i ], x3stCol[ i ],
      y1stCol[ i ], y2stCol[ i ], y3stCol[ i ], kernelRe, kernelIm );

    entryRe += wxyst[ i ] * kernelRe;
    entryIm += wxyst[ i ] * kernelIm;
  }

  SCVT areaIn = this->space->getRightMesh( )->getElemArea( innerElem );
  SCVT areaOut = this->space->getLeftMesh( )->getElemArea( outerElem );
  SCVT areaMult = areaIn * areaOut;

  entrySing *= areaOut;
  entryRe *= areaMult;
  entryIm *= areaMult;

  matrix.set( 0, 0, SC( entrySing + entryRe, entryIm ) );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix1LayerP1P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  SCVT stau, alpha1, alpha2;

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
  SCVT areay = this->space->getRightMesh( )->getElemArea( innerElem );
  SCVT areaxy = areax * areay;

  SCVT * x1refCol = this->x1strefCol;
  SCVT * x2refCol = this->x2strefCol;
  SCVT * x1stCol = this->x1stCol;
  SCVT * x2stCol = this->x2stCol;
  SCVT * x3stCol = this->x3stCol;
  SCVT * y1refCol = this->y1strefCol;
  SCVT * y2refCol = this->y2strefCol;
  SCVT * y1stCol = this->y1stCol;
  SCVT * y2stCol = this->y2stCol;
  SCVT * y3stCol = this->y3stCol;
  SCVT * wxyst = this->wxyst;

  SCVT entry11Re = 0.0;
  SCVT entry11Im = 0.0;
  SCVT entry12Re = 0.0;
  SCVT entry12Im = 0.0;
  SCVT entry13Re = 0.0;
  SCVT entry13Im = 0.0;
  SCVT entry21Re = 0.0;
  SCVT entry21Im = 0.0;
  SCVT entry22Re = 0.0;
  SCVT entry22Im = 0.0;
  SCVT entry23Re = 0.0;
  SCVT entry23Im = 0.0;
  SCVT entry31Re = 0.0;
  SCVT entry31Im = 0.0;
  SCVT entry32Re = 0.0;
  SCVT entry32Im = 0.0;
  SCVT entry33Re = 0.0;
  SCVT entry33Im = 0.0;
  SCVT kernelRe, kernelIm;
  SCVT phi1x, phi2x, phi3x, phi1y, phi2y, phi3y;

  const int width = DATA_WIDTH;
  const int align = DATA_ALIGN;
  int i;

  // getting collapsed quadrature data for regular part
  this->updateSteinbachRegularQuadratureNodes( x1, x2, x3, y1, y2, y3 );

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry11Re, entry11Im ) \
reduction( + : entry12Re, entry12Im ) \
reduction( + : entry13Re, entry13Im ) \
reduction( + : entry21Re, entry21Im ) \
reduction( + : entry22Re, entry22Im ) \
reduction( + : entry23Re, entry23Im ) \
reduction( + : entry31Re, entry31Im ) \
reduction( + : entry32Re, entry32Im ) \
reduction( + : entry33Re, entry33Im ) \
aligned( x1stCol, x2stCol, x3stCol, x1refCol, x2refCol, wxyst : align ) \
aligned( y1stCol, y2stCol, y3stCol, y1refCol, y2refCol : align ) \
private( kernelRe, kernelIm, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXYpad; ++i ) {
    this->evalShiftedKernel( x1stCol[ i ], x2stCol[ i ], x3stCol[ i ],
      y1stCol[ i ], y2stCol[ i ], y3stCol[ i ], kernelRe, kernelIm );

    kernelRe *= wxyst[ i ];
    kernelIm *= wxyst[ i ];

    phi1y = (SCVT) 1.0 - y1refCol[ i ] - y2refCol[ i ];
    phi2y = y1refCol[ i ];
    phi3y = y2refCol[ i ];
    phi1x = (SCVT) 1.0 - x1refCol[ i ] - x2refCol[ i ];
    phi2x = x1refCol[ i ];
    phi3x = x2refCol[ i ];

    entry11Re += kernelRe * phi1x * phi1y;
    entry11Im += kernelIm * phi1x * phi1y;
    entry12Re += kernelRe * phi1x * phi2y;
    entry12Im += kernelIm * phi1x * phi2y;
    entry13Re += kernelRe * phi1x * phi3y;
    entry13Im += kernelIm * phi1x * phi3y;
    entry21Re += kernelRe * phi2x * phi1y;
    entry21Im += kernelIm * phi2x * phi1y;
    entry22Re += kernelRe * phi2x * phi2y;
    entry22Im += kernelIm * phi2x * phi2y;
    entry23Re += kernelRe * phi2x * phi3y;
    entry23Im += kernelIm * phi2x * phi3y;
    entry31Re += kernelRe * phi3x * phi1y;
    entry31Im += kernelIm * phi3x * phi1y;
    entry32Re += kernelRe * phi3x * phi2y;
    entry32Im += kernelIm * phi3x * phi2y;
    entry33Re += kernelRe * phi3x * phi3y;
    entry33Im += kernelIm * phi3x * phi3y;
  }

  matrix.set( 0, 0, SC( entry11Re, entry11Im ) * areaxy );
  matrix.set( 0, 1, SC( entry12Re, entry12Im ) * areaxy );
  matrix.set( 0, 2, SC( entry13Re, entry13Im ) * areaxy );
  matrix.set( 1, 0, SC( entry21Re, entry21Im ) * areaxy );
  matrix.set( 1, 1, SC( entry22Re, entry22Im ) * areaxy );
  matrix.set( 1, 2, SC( entry23Re, entry23Im ) * areaxy );
  matrix.set( 2, 0, SC( entry31Re, entry31Im ) * areaxy );
  matrix.set( 2, 1, SC( entry32Re, entry32Im ) * areaxy );
  matrix.set( 2, 2, SC( entry33Re, entry33Im ) * areaxy );

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

    entry11Re = entry22Re = entry33Re = 0.0;

    // quadrature rule
#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry11Re, entry22Re, entry33Re ) \
aligned( sx, tx, ux, x1stref, x2stref, wxst : align ) \
private( kernelRe, phi1x, phi2x, phi3x ) \
simdlen( width )
    for ( i = 0; i < this->qSizeXpad; ++i ) {
      kernelRe = wxst[ i ] * this->collocationSingular1LayerP1( sx[ i ],
        tx[ i ], ux[ i ], stau, alpha1, alpha2 );

      phi1x = (SCVT) 1.0 - x1stref[ i ] - x2stref[ i ];
      phi2x = x1stref[ i ];
      phi3x = x2stref[ i ];

      entry11Re += kernelRe * phi1x;
      entry22Re += kernelRe * phi2x;
      entry33Re += kernelRe * phi3x;
    }

    matrix.add( 0, innerRot, entry11Re * areax );
    matrix.add( 1, innerRot, entry22Re * areax );
    matrix.add( 2, innerRot, entry33Re * areax );
  }
}

template<class LO, class SC>
SC BEIntegratorHelmholtz<LO, SC>::collocation1LayerP0(
  const SCVT* xLocal,
  SCVT stau,
  SCVT alpha1,
  SCVT alpha2,
  int quadratureOrder,
  const SC* shiftedKernel,
  LO elem
  ) const {

  SCVT lap = 0.0;
  SC ret( 0.0, 0.0 );

  lap += f1LayerP0( xLocal[ 0 ], xLocal[ 1 ], xLocal[ 2 ], stau, alpha2 );
  lap -= f1LayerP0( xLocal[ 0 ], xLocal[ 1 ], xLocal[ 2 ], stau, alpha1 );
  lap /= ( 4.0 * M_PI );

  for ( int j = 0; j < quadSizes[ quadratureOrder ]; j++ ) {
    ret += (SCVT) quadWeights[ quadratureOrder ][ j ] * shiftedKernel[ j ];
  }
  ret *= this->space->getRightMesh( )->getElemArea( elem );
  ret += lap;

  return ret;
}

template<class LO, class SC>
typename bem4i::BEIntegratorHelmholtz<LO, SC>::SCVT
BEIntegratorHelmholtz<LO, SC>::collocationSingular1LayerP0(
  SCVT sx,
  SCVT tx,
  SCVT ux,
  SCVT stau,
  SCVT alpha1,
  SCVT alpha2
  ) const {

  SCVT ret = 0.0;

  ret += f1LayerP0( sx, tx, ux, stau, alpha2 );
  ret -= f1LayerP0( sx, tx, ux, stau, alpha1 );

  return ( (SCVT) PI_FACT * ret );
}

template<class LO, class SC>
SC BEIntegratorHelmholtz<LO, SC>::collocation1LayerP1(
  const SCVT* xLocal,
  SCVT stau,
  SCVT alpha1,
  SCVT alpha2,
  int quadratureOrder,
  const SC* shiftedKernel,
  LO elem,
  int rot
  ) const {

  SCVT lap = 0.0;
  SC ret( 0.0, 0.0 );

  lap -= f1LayerP1( stau, alpha1, xLocal[ 1 ], xLocal[ 0 ], xLocal[ 2 ] );
  lap += f1LayerP1( stau, alpha2, xLocal[ 1 ], xLocal[ 0 ], xLocal[ 2 ] );
  lap /= stau;
  lap *= (SCVT) PI_FACT;

  SCVT phi = (SCVT) 0.0;

  for ( int j = 0; j < quadSizes[ quadratureOrder ]; j++ ) {
    if ( rot == 0 ) {
      phi = (SCVT) 1.0 - (SCVT) quadPoints_X1[ quadratureOrder ][ j ]
        - (SCVT) quadPoints_X2[ quadratureOrder ][ j ];
    } else if ( rot == 1 ) {
      phi = (SCVT) quadPoints_X1[ quadratureOrder ][ j ];
    } else if ( rot == 2 ) {
      phi = (SCVT) quadPoints_X2[ quadratureOrder ][ j ];
    }

    ret += (SCVT) quadWeights[ quadratureOrder ][ j ]
      * phi * shiftedKernel[ j ];
  }
  ret *= this->space->getRightMesh( )->getElemArea( elem );
  ret += lap;

  return ret;
}

template<class LO, class SC>
typename bem4i::BEIntegratorHelmholtz<LO, SC>::SCVT
BEIntegratorHelmholtz<LO, SC>::collocationSingular1LayerP1(
  SCVT sx,
  SCVT tx,
  SCVT ux,
  SCVT stau,
  SCVT alpha1,
  SCVT alpha2
  ) const {

  SCVT lap = (SCVT) 0.0;

  lap -= f1LayerP1( stau, alpha1, tx, sx, ux );
  lap += f1LayerP1( stau, alpha2, tx, sx, ux );
  lap /= stau;
  lap *= (SCVT) PI_FACT;

  return lap;
}

// Taken from GOBEM ... modified for s = s and s = 0 at once!

template<class LO, class SC>
typename bem4i::BEIntegratorHelmholtz<LO, SC>::SCVT
BEIntegratorHelmholtz<LO, SC>::f1LayerP0(
  SCVT sx,
  SCVT tx,
  SCVT ux,
  SCVT s,
  SCVT alpha
  ) const {

  SCVT f; // return value
  SCVT beta, gamma, p, q, a2; // substitution variable
  SCVT tmp1, tmp2, tmp3, tmp4; // aux variable
  SCVT beta_sq, q_sq; // aux variable

  SCVT f_0, tmp1_0, tmp2_0, tmp3_0; // variables for s = 0

  SCVT hh1, hh2, hh3, hh3_0;

  hh2 = ux * ux;

  // first term
  f = -s;
  f_0 = (SCVT) 0.0;

  // aux variables for substitution
  gamma = tx - alpha * sx;
  beta_sq = (SCVT) 1.0 + alpha * alpha;
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
  if ( std::abs( s - sx ) > (SCVT) EPS ) {
    // aux variables for backward substitution
    tmp3 = alpha * s - tx;

    // stable evaluation
    if ( tmp3 < (SCVT) 0.0 )
      hh3 = ( hh2 + ( s - sx ) * ( s - sx ) ) / ( tmp2 - tmp3 );
    else
      hh3 = tmp3 + tmp2;

  } else {
    hh3 = (SCVT) 1.0;
  }

  f += ( s - sx ) * std::log( hh3 );

  if ( std::abs( sx ) > (SCVT) EPS ) {
    // aux variables for backward substitution
    tmp3_0 = -tx;

    // stable evaluation
    if ( tmp3_0 < (SCVT) 0.0 )
      hh3_0 = ( hh2 + sx * sx ) / ( tmp2_0 - tmp3_0 );
    else
      hh3_0 = tmp3_0 + tmp2_0;

  } else {
    hh3_0 = (SCVT) 1.0;
  }

  f_0 += ( -sx ) * std::log( hh3_0 );

  hh1 = gamma / beta;

  // third term
  if ( std::abs( gamma ) > (SCVT) EPS ) {
    // stable evaluation
    if ( s - p < (SCVT) 0.0 )
      hh3 = q_sq / ( tmp2 - tmp1 );
    else
      hh3 = tmp1 + tmp2;

    if ( p > (SCVT) 0.0 )
      hh3_0 = q_sq / ( tmp2_0 - tmp1_0 );
    else
      hh3_0 = tmp1_0 + tmp2_0;

  } else {
    hh3 = hh3_0 = (SCVT) 1.0;
  }

  f -= hh1 * std::log( hh3 );
  f_0 -= hh1 * std::log( hh3_0 );

  // aux variables for backward substitution
  tmp4 = std::abs( ux );

  // fourth term
  if ( tmp4 > (SCVT) EPS ) {
    if ( std::abs( tmp1 ) < (SCVT) EPS )
      hh3 = alpha * q / tmp4;
    else
      hh3 = ( a2 * ( tmp2 - q ) / ( beta * tmp1 ) + alpha * q ) / tmp4;

    if ( std::abs( tmp1_0 ) < (SCVT) EPS )
      hh3_0 = alpha * q / tmp4;
    else
      hh3_0 = ( a2 * ( tmp2_0 - q ) / ( beta * tmp1_0 ) + alpha * q ) / tmp4;

  } else {
    hh3 = hh3_0 = (SCVT) 0.0;
  }

  f += (SCVT) 2.0 * tmp4 * std::atan( hh3 );
  f_0 += (SCVT) 2.0 * tmp4 * std::atan( hh3_0 );

  return ( f - f_0 );
}

// Taken from GOBEM

template<class LO, class SC>
typename bem4i::BEIntegratorHelmholtz<LO, SC>::SCVT
BEIntegratorHelmholtz<LO, SC>::f1LayerP1(
  SCVT stau,
  SCVT alpha,
  SCVT tx,
  SCVT sx,
  SCVT ux
  ) const {

  // Local variables
  SCVT beta, beta_2, f, p, q, v, a0, a1, a2, b1, b2, q_2, asl;
  SCVT a1_0, f_0, b1_0, b2_0, v_0;
  SCVT hh1, hh2, hh3, hh4, hh5, hh6, hh7, hh8, hh9;

  hh3 = ux * ux;
  hh4 = tx - alpha * sx;
  hh5 = hh4 * hh4;
  hh6 = tx - alpha * stau;
  hh7 = sx - stau * (SCVT) 2.0;
  hh8 = stau - sx;
  hh9 = hh8 * hh8;

  beta_2 = (SCVT) 1.0 + alpha * alpha;
  beta = std::sqrt( beta_2 );
  p = ( alpha * tx + sx ) / beta_2;
  q_2 = hh3 + hh5 / beta_2;
  q = std::sqrt( q_2 );
  a1 = std::sqrt( hh9 + hh6 * hh6 + hh3 );
  a1_0 = std::sqrt( sx * sx + tx * tx + hh3 );

  hh1 = (SCVT) 0.5 / beta_2 * hh4;

  f = (SCVT) 0.25 + stau * stau + stau * (SCVT) 0.5 * hh7 + hh1 * ( a1 + q );
  f_0 = (SCVT) 0.25 + hh1 * ( a1_0 + q );

  if ( sx - stau != (SCVT) 0.0 ) {
    if ( hh6 > (SCVT) 0.0 )
      b1 = ( hh9 + hh3 ) / ( a1 + hh6 );
    else
      b1 = a1 - hh6;

    f += -hh8 * (SCVT) 0.5 * ( stau + hh7 ) * std::log( b1 );
  }

  if ( sx != (SCVT) 0.0 ) {
    if ( tx > (SCVT) 0.0 )
      b1_0 = ( sx * sx + hh3 ) / ( a1_0 + tx );
    else
      b1_0 = -tx + a1_0;

    f_0 += sx * (SCVT) 0.5 * hh7 * std::log( b1_0 );
  }

  b2 = stau - sx - alpha * hh6;
  b2_0 = -sx - alpha * tx;

  if ( b2 < (SCVT) 0.0 )
    b1 = ( hh3 + hh5 / beta_2 ) / ( a1 - b2 / beta );
  else
    b1 = b2 / beta + a1;

  if ( b2_0 < (SCVT) 0.0 )
    b1_0 = ( hh3 + hh5 / beta_2 )
    / ( a1_0 - b2_0 / beta );
  else
    b1_0 = b2_0 / beta + a1_0;

  hh2 = (SCVT) 0.5 / beta * ( alpha * q_2 + hh4 * (SCVT) 2.0 * ( -hh8 ) );

  if ( b1 != (SCVT) 0.0 )
    f += hh2 * std::log( b1 );

  if ( b1_0 != (SCVT) 0.0 )
    f_0 += hh2 * std::log( b1_0 );

  if ( ux != (SCVT) 0.0 ) {
    f -= hh3 * (SCVT) 0.5 * std::log( std::abs( q + a1 ) );
    f_0 -= hh3 * (SCVT) 0.5 * std::log( std::abs( q + a1_0 ) );
    v = beta * ( stau - p ) / ( a1 + q );
    v_0 = -beta * p / ( a1_0 + q );
    a2 = beta_2 * q + hh4;
    a0 = beta_2 * q - hh4;
    a1 = alpha * (SCVT) 2.0 * beta * q;
    b2 = a2 * v * v + a1 * v + a0;
    b2_0 = a2 * v_0 * v_0 + a1 * v_0 + a0;

    if ( b2 != (SCVT) 0.0 )
      f -= hh3 * (SCVT) 0.5 * std::log( std::abs( b2 ) );

    if ( b2_0 != (SCVT) 0.0 )
      f_0 -= hh3 * (SCVT) 0.5 * std::log( std::abs( b2_0 ) );

    f += ux * (SCVT) 2.0 * hh8
      * std::atan( ( a2 * (SCVT) 2.0 * v + a1 )
      / ( beta * (SCVT) 2.0 * ux ) );

    f_0 += ux * (SCVT) 2.0 * hh8
      * std::atan( ( a2 * (SCVT) 2.0 * v_0 + a1 )
      / ( beta * (SCVT) 2.0 * ux ) );
  }

  return ( f - f_0 );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix2LayerP0P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  if ( outerElem == innerElem ) {
    matrix.setAll( SC( 0.0, 0.0 ) );
    return;
  }

  SCVT stau, alpha1, alpha2;

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
  this->space->getRightMesh( )->getNormal( innerElem, n );

  SCVT * x1stCol = this->x1stCol;
  SCVT * x2stCol = this->x2stCol;
  SCVT * x3stCol = this->x3stCol;
  SCVT * y1refCol = this->y1strefCol;
  SCVT * y2refCol = this->y2strefCol;
  SCVT * y1stCol = this->y1stCol;
  SCVT * y2stCol = this->y2stCol;
  SCVT * y3stCol = this->y3stCol;
  SCVT * wxyst = this->wxyst;

  SCVT entry1Re = 0.0;
  SCVT entry1Im = 0.0;
  SCVT entry2Re = 0.0;
  SCVT entry2Im = 0.0;
  SCVT entry3Re = 0.0;
  SCVT entry3Im = 0.0;
  SCVT kernelRe;
  SCVT kernelIm;
  SCVT phi1y, phi2y, phi3y;

  int i;
  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

  // getting collapsed quadrature data for regular part
  this->updateSteinbachRegularQuadratureNodes( x1, x2, x3, y1, y2, y3 );

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry1Re, entry1Im ) \
reduction( + : entry2Re, entry2Im ) \
reduction( + : entry3Re, entry3Im ) \
aligned( x1stCol, x2stCol, x3stCol, wxyst : align ) \
aligned( y1stCol, y2stCol, y3stCol, y1refCol, y2refCol : align ) \
private( kernelRe, kernelIm, phi1y, phi2y, phi3y ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXYpad; ++i ) {
    this->evalShiftedKernelNeumann( x1stCol[ i ], x2stCol[ i ], x3stCol[ i ],
      y1stCol[ i ], y2stCol[ i ], y3stCol[ i ], n[ 0 ], n[ 1 ], n[ 2 ],
      kernelRe, kernelIm );

    kernelRe *= wxyst[ i ];
    kernelIm *= wxyst[ i ];

    phi1y = (SCVT) 1.0 - y1refCol[ i ] - y2refCol[ i ];
    phi2y = y1refCol[ i ];
    phi3y = y2refCol[ i ];

    entry1Re += kernelRe * phi1y;
    entry1Im += kernelIm * phi1y;
    entry2Re += kernelRe * phi2y;
    entry2Im += kernelIm * phi2y;
    entry3Re += kernelRe * phi3y;
    entry3Im += kernelIm * phi3y;
  }

  matrix.set( 0, 0, SC( entry1Re, entry1Im ) );
  matrix.set( 0, 1, SC( entry2Re, entry2Im ) );
  matrix.set( 0, 2, SC( entry3Re, entry3Im ) );
  matrix.scale( this->space->getRightMesh( )->getElemArea( innerElem ) );

  SCVT * sx = this->sx;
  SCVT * tx = this->tx;
  SCVT * ux = this->ux;
  SCVT * wxst = this->wxst;

  SCVT entrySing;

  for ( int innerRot = 0; innerRot < 3; ++innerRot ) {
    // getting local coordinate system
    this->getLocalCoordinates( y1, y2, y3, r1, r2, n, stau, alpha1, alpha2,
      innerRot );
    // getting outer quadrature points in local coordinate system
    this->updateSteinbachLocalQuadratureNodes( x1, x2, x3, y1, y2, y3, r1, r2,
      n, innerRot );

    entrySing = 0.0;

    // quadrature rule
#pragma omp simd \
linear( i : 1 ) \
reduction( + : entrySing ) \
aligned( sx, tx, ux, wxst : align ) \
simdlen( width )
    for ( i = 0; i < this->qSizeXpad; ++i ) {
      entrySing += wxst[ i ] * collocationSingular2LayerP1( sx[ i ],
        tx[ i ], ux[ i ], stau, alpha1, alpha2 );
    }

    matrix.add( 0, innerRot, entrySing );
  }

  matrix.scale( this->space->getLeftMesh( )->getElemArea( outerElem ) );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix2LayerP0P0(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  if ( outerElem == innerElem ) {
    matrix.setAll( SC( 0.0, 0.0 ) );
    return;
  }

  SCVT stau, alpha1, alpha2;

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

  // getting collapsed quadrature data for regular part
  this->updateSteinbachRegularQuadratureNodes( x1, x2, x3, y1, y2, y3 );

  SCVT * sx = this->sx;
  SCVT * tx = this->tx;
  SCVT * ux = this->ux;
  SCVT * wxst = this->wxst;

  SCVT entrySing = 0.0;

  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;
  int i;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entrySing ) \
aligned( sx, tx, ux, wxst : align ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXpad; ++i ) {
    entrySing += wxst[ i ] * this->collocationSingular2LayerP0( sx[ i ],
      tx[ i ], ux[ i ], stau, alpha1, alpha2 );
  }

  SCVT * x1stCol = this->x1stCol;
  SCVT * x2stCol = this->x2stCol;
  SCVT * x3stCol = this->x3stCol;
  SCVT * y1stCol = this->y1stCol;
  SCVT * y2stCol = this->y2stCol;
  SCVT * y3stCol = this->y3stCol;
  SCVT * wxyst = this->wxyst;

  SCVT entryRe = 0.0;
  SCVT entryIm = 0.0;
  SCVT kernelRe;
  SCVT kernelIm;

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entryRe, entryIm ) \
aligned( x1stCol, x2stCol, x3stCol, wxyst : align ) \
aligned( y1stCol, y2stCol, y3stCol : align ) \
private( kernelRe, kernelIm ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXYpad; ++i ) {
    this->evalShiftedKernelNeumann( x1stCol[ i ], x2stCol[ i ], x3stCol[ i ],
      y1stCol[ i ], y2stCol[ i ], y3stCol[ i ], n[ 0 ], n[ 1 ], n[ 2 ],
      kernelRe, kernelIm );

    entryRe += wxyst[ i ] * kernelRe;
    entryIm += wxyst[ i ] * kernelIm;
  }

  SCVT areaIn = this->space->getRightMesh( )->getElemArea( innerElem );
  SCVT areaOut = this->space->getLeftMesh( )->getElemArea( outerElem );
  SCVT areaMult = areaIn * areaOut;

  entrySing *= areaOut;
  entryRe *= areaMult;
  entryIm *= areaMult;

  matrix.set( 0, 0, SC( entrySing + entryRe, entryIm ) );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrix2LayerP1P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  if ( outerElem == innerElem ) {
    matrix.setAll( SC( 0.0, 0.0 ) );
    return;
  }

  SCVT stau, alpha1, alpha2;

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
  this->space->getRightMesh( )->getNormal( innerElem, n );
  SCVT areax = this->space->getLeftMesh( )->getElemArea( outerElem );
  SCVT areay = this->space->getRightMesh( )->getElemArea( innerElem );
  SCVT areaxy = areax * areay;

  SCVT * x1refCol = this->x1strefCol;
  SCVT * x2refCol = this->x2strefCol;
  SCVT * x1stCol = this->x1stCol;
  SCVT * x2stCol = this->x2stCol;
  SCVT * x3stCol = this->x3stCol;
  SCVT * y1refCol = this->y1strefCol;
  SCVT * y2refCol = this->y2strefCol;
  SCVT * y1stCol = this->y1stCol;
  SCVT * y2stCol = this->y2stCol;
  SCVT * y3stCol = this->y3stCol;
  SCVT * wxyst = this->wxyst;

  SCVT entry11Re = 0.0;
  SCVT entry11Im = 0.0;
  SCVT entry12Re = 0.0;
  SCVT entry12Im = 0.0;
  SCVT entry13Re = 0.0;
  SCVT entry13Im = 0.0;
  SCVT entry21Re = 0.0;
  SCVT entry21Im = 0.0;
  SCVT entry22Re = 0.0;
  SCVT entry22Im = 0.0;
  SCVT entry23Re = 0.0;
  SCVT entry23Im = 0.0;
  SCVT entry31Re = 0.0;
  SCVT entry31Im = 0.0;
  SCVT entry32Re = 0.0;
  SCVT entry32Im = 0.0;
  SCVT entry33Re = 0.0;
  SCVT entry33Im = 0.0;
  SCVT kernelRe, kernelIm;
  SCVT phi1x, phi2x, phi3x, phi1y, phi2y, phi3y;

  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

  int i;

  // getting collapsed quadrature data for regular part
  this->updateSteinbachRegularQuadratureNodes( x1, x2, x3, y1, y2, y3 );

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry11Re, entry11Im ) \
reduction( + : entry12Re, entry12Im ) \
reduction( + : entry13Re, entry13Im ) \
reduction( + : entry21Re, entry21Im ) \
reduction( + : entry22Re, entry22Im ) \
reduction( + : entry23Re, entry23Im ) \
reduction( + : entry31Re, entry31Im ) \
reduction( + : entry32Re, entry32Im ) \
reduction( + : entry33Re, entry33Im ) \
aligned( x1stCol, x2stCol, x3stCol, x1refCol, x2refCol, wxyst : align ) \
aligned( y1stCol, y2stCol, y3stCol, y1refCol, y2refCol : align ) \
private( kernelRe, kernelIm, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXYpad; ++i ) {
    this->evalShiftedKernelNeumann( x1stCol[ i ], x2stCol[ i ], x3stCol[ i ],
      y1stCol[ i ], y2stCol[ i ], y3stCol[ i ], n[ 0 ], n[ 1 ], n[ 2 ],
      kernelRe, kernelIm );

    kernelRe *= wxyst[ i ];
    kernelIm *= wxyst[ i ];

    phi1y = (SCVT) 1.0 - y1refCol[ i ] - y2refCol[ i ];
    phi2y = y1refCol[ i ];
    phi3y = y2refCol[ i ];
    phi1x = (SCVT) 1.0 - x1refCol[ i ] - x2refCol[ i ];
    phi2x = x1refCol[ i ];
    phi3x = x2refCol[ i ];

    entry11Re += kernelRe * phi1x * phi1y;
    entry11Im += kernelIm * phi1x * phi1y;
    entry12Re += kernelRe * phi1x * phi2y;
    entry12Im += kernelIm * phi1x * phi2y;
    entry13Re += kernelRe * phi1x * phi3y;
    entry13Im += kernelIm * phi1x * phi3y;
    entry21Re += kernelRe * phi2x * phi1y;
    entry21Im += kernelIm * phi2x * phi1y;
    entry22Re += kernelRe * phi2x * phi2y;
    entry22Im += kernelIm * phi2x * phi2y;
    entry23Re += kernelRe * phi2x * phi3y;
    entry23Im += kernelIm * phi2x * phi3y;
    entry31Re += kernelRe * phi3x * phi1y;
    entry31Im += kernelIm * phi3x * phi1y;
    entry32Re += kernelRe * phi3x * phi2y;
    entry32Im += kernelIm * phi3x * phi2y;
    entry33Re += kernelRe * phi3x * phi3y;
    entry33Im += kernelIm * phi3x * phi3y;
  }

  matrix.set( 0, 0, SC( entry11Re, entry11Im ) * areaxy );
  matrix.set( 0, 1, SC( entry12Re, entry12Im ) * areaxy );
  matrix.set( 0, 2, SC( entry13Re, entry13Im ) * areaxy );
  matrix.set( 1, 0, SC( entry21Re, entry21Im ) * areaxy );
  matrix.set( 1, 1, SC( entry22Re, entry22Im ) * areaxy );
  matrix.set( 1, 2, SC( entry23Re, entry23Im ) * areaxy );
  matrix.set( 2, 0, SC( entry31Re, entry31Im ) * areaxy );
  matrix.set( 2, 1, SC( entry32Re, entry32Im ) * areaxy );
  matrix.set( 2, 2, SC( entry33Re, entry33Im ) * areaxy );

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

    entry11Re = entry22Re = entry33Re = 0.0;

    // quadrature rule
#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry11Re, entry22Re, entry33Re ) \
aligned( sx, tx, ux, x1stref, x2stref, wxst : align ) \
private( kernelRe, phi1x, phi2x, phi3x ) \
simdlen( width )
    for ( i = 0; i < this->qSizeXpad; ++i ) {
      kernelRe = wxst[ i ] * this->collocationSingular2LayerP1( sx[ i ],
        tx[ i ], ux[ i ], stau, alpha1, alpha2 );

      phi1x = (SCVT) 1.0 - x1stref[ i ] - x2stref[ i ];
      phi2x = x1stref[ i ];
      phi3x = x2stref[ i ];

      entry11Re += kernelRe * phi1x;
      entry22Re += kernelRe * phi2x;
      entry33Re += kernelRe * phi3x;
    }

    matrix.add( 0, innerRot, entry11Re * areax );
    matrix.add( 1, innerRot, entry22Re * areax );
    matrix.add( 2, innerRot, entry33Re * areax );
  }
}

template<class LO, class SC>
SC BEIntegratorHelmholtz<LO, SC>::collocation2LayerP1(
  const SCVT* xLocal,
  SCVT stau,
  SCVT alpha1,
  SCVT alpha2,
  int quadratureOrder,
  const SC* shiftedKernelNeumann,
  LO elem,
  int rot
  ) const {

  if ( fabs( xLocal[ 2 ] ) < EPS ) {
    return 0.0;
  }
  SCVT lap = 0.0;
  SC ret( 0.0, 0.0 );

  lap += f2LayerP1( alpha2, xLocal[ 1 ], xLocal[ 0 ], xLocal[ 2 ], stau );
  lap -= f2LayerP1( alpha1, xLocal[ 1 ], xLocal[ 0 ], xLocal[ 2 ], stau );
  lap /= ( stau * 4.0 * M_PI );

  SCVT phi = (SCVT) 0.0;

  for ( int j = 0; j < quadSizes[ quadratureOrder ]; j++ ) {
    if ( rot == 0 ) {
      phi = (SCVT) 1.0 - (SCVT) quadPoints_X1[ quadratureOrder ][ j ]
        - (SCVT) quadPoints_X2[ quadratureOrder ][ j ];
    } else if ( rot == 1 ) {
      phi = (SCVT) quadPoints_X1[ quadratureOrder ][ j ];
    } else if ( rot == 2 ) {
      phi = (SCVT) quadPoints_X2[ quadratureOrder ][ j ];
    }

    ret += (SCVT) quadWeights[ quadratureOrder ][ j ]
      * phi * shiftedKernelNeumann[ j ];
  }
  ret *= this->space->getRightMesh( )->getElemArea( elem );
  ret += lap;

  return ret;
}

template<class LO, class SC>
typename bem4i::BEIntegratorHelmholtz<LO, SC>::SCVT
BEIntegratorHelmholtz<LO, SC>::collocationSingular2LayerP1(
  SCVT sx,
  SCVT tx,
  SCVT ux,
  SCVT stau,
  SCVT alpha1,
  SCVT alpha2
  ) const {

  SCVT ret;

  // helps geometries with planes
  //if ( std::abs( ux ) < (SCVT) EPS ) {
  //ret = (SCVT) 0.0;
  //} else {
  SCVT f1 = f2LayerP1( alpha2, tx, sx, ux, stau );
  SCVT f2 = f2LayerP1( alpha1, tx, sx, ux, stau );
  ret = (SCVT) PI_FACT * ( f1 - f2 ) / stau;
  //}

  return ret;
}

template<class LO, class SC>
typename bem4i::BEIntegratorHelmholtz<LO, SC>::SCVT
BEIntegratorHelmholtz<LO, SC>::collocationSingular2LayerP0(
  SCVT sx,
  SCVT tx,
  SCVT ux,
  SCVT stau,
  SCVT alpha1,
  SCVT alpha2
  ) const {

  SCVT ret;

  // helps geometries with planes
  //if ( std::abs( ux ) < (SC) EPS ) {
  //ret = (SC) 0.0;
  //} else {
  SCVT f1 = f2LayerP0( alpha2, tx, sx, ux, stau );
  SCVT f2 = f2LayerP0( alpha1, tx, sx, ux, stau );
  ret = (SCVT) PI_FACT * ( f1 - f2 );
  //}

  return ret;
}

// taken from GOBEM

template<class LO, class SC>
typename bem4i::BEIntegratorHelmholtz<LO, SC>::SCVT
BEIntegratorHelmholtz<LO, SC>::f2LayerP1(
  SCVT alpha,
  SCVT tx,
  SCVT sx,
  SCVT ux,
  SCVT stau
  ) const {

  // Local variables
  SCVT ret_val = (SCVT) 0.0;
  SCVT beta, beta_2, a, b, f, q, f1, f2, h1, h2, h3, h4, q2, aq;
  SCVT v1, v2, gm, gp, v1m, v2m, v1p, v2p, ask, asl, tau_2, sl_2, tl_2, alpha_2;
  // our modifications for vectorization
  SCVT hh1, hh2, hh3, hh4, hh5;

  // aux variables
  beta_2 = alpha * alpha + (SCVT) 1.0;
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
  if ( ( tx - asl ) >= (SCVT) 0.0 ) {
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
  if ( ( a - tl_2 ) < (SCVT) 1e-14 ) //  (a - b < 1e-14)
  {
    a = std::sqrt( a ); // sqrt(sl_2 + tl_2 + tau_2);
    b = alpha * ( sx + alpha * tx ); // onetime evaluation instead of 2 ifs

    hh3 = alpha_2 * ( sl_2 + tau_2 ) - sx * ( sx + alpha * (SCVT) 2.0 * tx );

    if ( b > (SCVT) 0.0 ) // (alpha * (sx + alpha * tx) > 0.0)
    {
      v1p = hh3;
      v1p /= alpha * a + sx + alpha * tx;
      v1p = q * v1p + alpha * tau_2 - ( tx - alpha * sx ) * sx;
      v1p /= q + a;
    }
    // if (alpha * (sx + alpha * tx) < 0.0)   // else eingefuegt
    if ( b < (SCVT) 0.0 ) {
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
  f = ux * (SCVT) 0.5 * std::log( h4 * h1 / ( h3 * h2 ) );

  f1 = std::atan( v2m / ux ) - std::atan( v1m / ux ); // std::atan
  f2 = std::atan( v2p / ux ) - std::atan( v1p / ux ); // std::atan

  // add together
  f += ( stau - sx ) * ( f1 - f2 );

  // -ux * alpha/beta * u
  h1 = -( sx ) - alpha * tx;
  h2 = stau - sx + alpha * ( ask - tx );

  hh4 = std::sqrt( h1 * h1 + beta_2 * q2 );

  // stable evaluation
  if ( h1 < (SCVT) 0.0 )
    f1 = beta_2 * q2 / ( hh4 - h1 );
  else
    f1 = h1 + hh4;

  hh5 = std::sqrt( h2 * h2 + beta_2 * q2 );

  // stable evaluation
  if ( h2 < (SCVT) 0.0 )
    f2 = beta_2 * q2 / ( hh5 - h2 );
  else
    f2 = h2 + hh5;

  if ( std::abs( ux ) > (SCVT) EPS )
    ret_val = f - ux * alpha / beta * std::log( f2 / f1 );

  return ret_val;
}

template<class LO, class SC>
typename bem4i::BEIntegratorHelmholtz<LO, SC>::SCVT
BEIntegratorHelmholtz<LO, SC>::f2LayerP0(
  SCVT alpha,
  SCVT tx,
  SCVT sx,
  SCVT ux,
  SCVT stau
  ) const {

  // Local variables
  SCVT beta_2, a, b, q, f1, f2, q2, aq;
  SCVT v1, v2, gm, gp, v1m, v2m, v1p, v2p, ask, asl, tau_2, sl_2, tl_2, alpha_2;
  SCVT sk_minus_sl, ask_minus_tl, tl_minus_asl, sl_plus_alpha_tl;

  SCVT hh1, hh2, hh3, hh4, hh5;

  // Aux variables
  beta_2 = alpha * alpha + (SCVT) 1.0;
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
  if ( tl_minus_asl >= (SCVT) 0.0 ) {
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
  if ( a - tl_2 < (SCVT) 1e-14 ) //  (a - b < 1e-14)
  {
    a = std::sqrt( a ); // sqrt(sl_2 + tl_2 + tau_2);
    b = alpha * sl_plus_alpha_tl;
    hh3 = alpha_2 * ( sl_2 + tau_2 ) - sx * ( sx + alpha * (SCVT) 2.0 * tx );
    hh4 = alpha * tau_2 - tl_minus_asl * sx;
    hh5 = q + a;

    if ( b > (SCVT) 0.0 ) // (alpha * (sx + alpha * tx) > 0.0)
    {
      v1p = hh3;
      v1p /= alpha * a + sx + alpha * tx;
      v1p = q * v1p + hh4;
      v1p /= hh5;
    }
    if ( b < (SCVT) 0.0 ) // if (alpha * (sx + alpha * tx) < 0.0)
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
void BEIntegratorHelmholtz<LO, SC>::computeElemMatrixHypersingularP1P1(
  LO outerElem,
  LO innerElem,
  FullMatrix<LO, SC> & matrix
  ) const {

  SCVT nx[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT ny[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  this->space->getRightMesh( )->getNormal( innerElem, ny );
  this->space->getLeftMesh( )->getNormal( outerElem, nx );

  SCVT normalDot = DOT3( nx, ny );

  SCVT stau, alpha1, alpha2;

  SCVT x1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT x3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT y3[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r1[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );
  SCVT r2[ 3 ] __attribute__ ( ( aligned( DATA_ALIGN ) ) );

  this->space->getRightMesh( )->getNodes( innerElem, y1, y2, y3 );
  this->space->getLeftMesh( )->getNodes( outerElem, x1, x2, x3 );
  SCVT areax = this->space->getLeftMesh( )->getElemArea( outerElem );
  SCVT areay = this->space->getRightMesh( )->getElemArea( innerElem );

  SC mult = -this->kappa * this->kappa * normalDot;
  SC areaxym = areax * areay * mult;
  SC areaxm = areax * mult;

  SCVT * x1refCol = this->x1strefCol;
  SCVT * x2refCol = this->x2strefCol;
  SCVT * x1stCol = this->x1stCol;
  SCVT * x2stCol = this->x2stCol;
  SCVT * x3stCol = this->x3stCol;
  SCVT * y1refCol = this->y1strefCol;
  SCVT * y2refCol = this->y2strefCol;
  SCVT * y1stCol = this->y1stCol;
  SCVT * y2stCol = this->y2stCol;
  SCVT * y3stCol = this->y3stCol;
  SCVT * wxyst = this->wxyst;

  SCVT entryVRe = 0.0;
  SCVT entryVIm = 0.0;
  SCVT entryVSing = 0.0;
  SCVT entry11Re = 0.0;
  SCVT entry11Im = 0.0;
  SCVT entry12Re = 0.0;
  SCVT entry12Im = 0.0;
  SCVT entry13Re = 0.0;
  SCVT entry13Im = 0.0;
  SCVT entry21Re = 0.0;
  SCVT entry21Im = 0.0;
  SCVT entry22Re = 0.0;
  SCVT entry22Im = 0.0;
  SCVT entry23Re = 0.0;
  SCVT entry23Im = 0.0;
  SCVT entry31Re = 0.0;
  SCVT entry31Im = 0.0;
  SCVT entry32Re = 0.0;
  SCVT entry32Im = 0.0;
  SCVT entry33Re = 0.0;
  SCVT entry33Im = 0.0;
  SCVT kernelRe, kernelIm;
  SCVT phi1x, phi2x, phi3x, phi1y, phi2y, phi3y;

  const int align = DATA_ALIGN;
  const int width = DATA_WIDTH;

  int i;

  // getting collapsed quadrature data for regular part
  this->updateSteinbachRegularQuadratureNodes( x1, x2, x3, y1, y2, y3 );

#pragma omp simd \
linear( i : 1 ) \
reduction( + : entryVRe, entryVIm ) \
reduction( + : entry11Re, entry11Im ) \
reduction( + : entry12Re, entry12Im ) \
reduction( + : entry13Re, entry13Im ) \
reduction( + : entry21Re, entry21Im ) \
reduction( + : entry22Re, entry22Im ) \
reduction( + : entry23Re, entry23Im ) \
reduction( + : entry31Re, entry31Im ) \
reduction( + : entry32Re, entry32Im ) \
reduction( + : entry33Re, entry33Im ) \
aligned( x1stCol, x2stCol, x3stCol, x1refCol, x2refCol, wxyst : align ) \
aligned( y1stCol, y2stCol, y3stCol, y1refCol, y2refCol : align ) \
private( kernelRe, kernelIm, phi1x, phi2x, phi3x, phi1y, phi2y, phi3y ) \
simdlen( width )
  for ( i = 0; i < this->qSizeXYpad; ++i ) {
    this->evalShiftedKernel( x1stCol[ i ], x2stCol[ i ], x3stCol[ i ],
      y1stCol[ i ], y2stCol[ i ], y3stCol[ i ], kernelRe, kernelIm );

    kernelRe *= wxyst[ i ];
    kernelIm *= wxyst[ i ];

    phi1y = (SCVT) 1.0 - y1refCol[ i ] - y2refCol[ i ];
    phi2y = y1refCol[ i ];
    phi3y = y2refCol[ i ];
    phi1x = (SCVT) 1.0 - x1refCol[ i ] - x2refCol[ i ];
    phi2x = x1refCol[ i ];
    phi3x = x2refCol[ i ];

    entryVRe += kernelRe;
    entryVIm += kernelIm;
    entry11Re += kernelRe * phi1x * phi1y;
    entry11Im += kernelIm * phi1x * phi1y;
    entry12Re += kernelRe * phi1x * phi2y;
    entry12Im += kernelIm * phi1x * phi2y;
    entry13Re += kernelRe * phi1x * phi3y;
    entry13Im += kernelIm * phi1x * phi3y;
    entry21Re += kernelRe * phi2x * phi1y;
    entry21Im += kernelIm * phi2x * phi1y;
    entry22Re += kernelRe * phi2x * phi2y;
    entry22Im += kernelIm * phi2x * phi2y;
    entry23Re += kernelRe * phi2x * phi3y;
    entry23Im += kernelIm * phi2x * phi3y;
    entry31Re += kernelRe * phi3x * phi1y;
    entry31Im += kernelIm * phi3x * phi1y;
    entry32Re += kernelRe * phi3x * phi2y;
    entry32Im += kernelIm * phi3x * phi2y;
    entry33Re += kernelRe * phi3x * phi3y;
    entry33Im += kernelIm * phi3x * phi3y;
  }

  entryVRe *= areax * areay;
  entryVIm *= areax * areay;
  matrix.set( 0, 0, SC( entry11Re, entry11Im ) * areaxym );
  matrix.set( 0, 1, SC( entry12Re, entry12Im ) * areaxym );
  matrix.set( 0, 2, SC( entry13Re, entry13Im ) * areaxym );
  matrix.set( 1, 0, SC( entry21Re, entry21Im ) * areaxym );
  matrix.set( 1, 1, SC( entry22Re, entry22Im ) * areaxym );
  matrix.set( 1, 2, SC( entry23Re, entry23Im ) * areaxym );
  matrix.set( 2, 0, SC( entry31Re, entry31Im ) * areaxym );
  matrix.set( 2, 1, SC( entry32Re, entry32Im ) * areaxym );
  matrix.set( 2, 2, SC( entry33Re, entry33Im ) * areaxym );

  SCVT * x1stref = this->x1stref;
  SCVT * x2stref = this->x2stref;
  SCVT * sx = this->sx;
  SCVT * tx = this->tx;
  SCVT * ux = this->ux;
  SCVT * wxst = this->wxst;

  for ( int innerRot = 0; innerRot < 3; ++innerRot ) {
    // getting local coordinate system
    this->getLocalCoordinates( y1, y2, y3, r1, r2, ny, stau, alpha1, alpha2,
      innerRot );
    // getting outer quadrature points in local coordinate system
    this->updateSteinbachLocalQuadratureNodes( x1, x2, x3, y1, y2, y3, r1, r2,
      ny, innerRot );

    if ( innerRot == 0 ) {
#pragma omp simd \
linear( i : 1 ) \
reduction( + : entryVSing ) \
aligned( sx, tx, ux, wxst : align ) \
simdlen( width )
      for ( i = 0; i < this->qSizeXpad; ++i ) {
        entryVSing += wxst[ i ] * this->collocationSingular1LayerP0(
          sx[ i ], tx[ i ], ux[ i ], stau, alpha1, alpha2 );
      }

      entryVSing *= areax;
    }

    entry11Re = entry22Re = entry33Re = 0.0;

    // quadrature rule
#pragma omp simd \
linear( i : 1 ) \
reduction( + : entry11Re, entry22Re, entry33Re ) \
aligned( sx, tx, ux, x1stref, x2stref, wxst : align ) \
private( kernelRe, phi1x, phi2x, phi3x ) \
simdlen( width )
    for ( i = 0; i < this->qSizeXpad; ++i ) {
      kernelRe = wxst[ i ] * this->collocationSingular1LayerP1( sx[ i ],
        tx[ i ], ux[ i ], stau, alpha1, alpha2 );

      phi1x = (SCVT) 1.0 - x1stref[ i ] - x2stref[ i ];
      phi2x = x1stref[ i ];
      phi3x = x2stref[ i ];

      entry11Re += kernelRe * phi1x;
      entry22Re += kernelRe * phi2x;
      entry33Re += kernelRe * phi3x;
    }

    matrix.add( 0, innerRot, entry11Re * areaxm );
    matrix.add( 1, innerRot, entry22Re * areaxm );
    matrix.add( 2, innerRot, entry33Re * areaxm );
  }

  SCVT * outerCurl = this->space->getLeftMesh( )->getCurls( )->getData( ) +
    outerElem * 9;
  SCVT * innerCurl = this->space->getRightMesh( )->getCurls( )->getData( ) +
    innerElem * 9;

  SC entryV( entryVRe + entryVSing, entryVIm );
  SCVT curlDot1, curlDot2, curlDot3;
  SC * matrixData = matrix.getData( );

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

    matrixData[ innerRot * 3 ] += curlDot1 * entryV;
    matrixData[ innerRot * 3 + 1 ] += curlDot2 * entryV;
    matrixData[ innerRot * 3 + 2 ] += curlDot3 * entryV;
  }

}
//*/

template<class LO, class SC>
SC BEIntegratorHelmholtz<LO, SC>::evalShiftedKernel(
  const SCVT* x,
  const SCVT* y
  ) const {

  SC i( 0.0, 1.0 );
  if ( std::abs( x[0] - y[0] ) < EPS && std::abs( x[1] - y[1] ) < EPS && std::abs( x[2] - y[2] ) < EPS ) {
    return (SCVT) PI_FACT * ( std::real( kappa ) * i - std::imag( kappa ) );
  } else {
    SCVT norm = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) + ( x[1] - y[1] )*( x[1] - y[1] ) + ( x[2] - y[2] )*( x[2] - y[2] ) );
    return (SCVT) PI_FACT * ( std::exp( i * kappa * norm ) - (SCVT) 1.0 ) / norm;
  }
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::evalShiftedKernel(
  SCVT x1,
  SCVT x2,
  SCVT x3,
  SCVT y1,
  SCVT y2,
  SCVT y3,
  SCVT & real,
  SCVT & imag
  ) const {

  SCVT rekappa = std::real( kappa );
  SCVT imkappa = std::imag( kappa );

  if ( std::abs( x1 - y1 ) < (SCVT) EPS
    && std::abs( x2 - y2 ) < (SCVT) EPS
    && std::abs( x3 - y3 ) < (SCVT) EPS ) {

    real = (SCVT) PI_FACT * ( -imkappa );
    imag = (SCVT) PI_FACT * rekappa;
  } else {

    SCVT norm = std::sqrt(
      ( x1 - y1 ) * ( x1 - y1 )
      + ( x2 - y2 ) * ( x2 - y2 )
      + ( x3 - y3 ) * ( x3 - y3 ) );
    SCVT rekappan = rekappa * norm;

    SCVT expim = std::exp( -imkappa * norm );

    real = (SCVT) PI_FACT *
      ( expim * std::cos( rekappan ) - (SCVT) 1.0 ) / norm;
    imag = (SCVT) PI_FACT *
      ( expim * std::sin( rekappan ) ) / norm;
  }
}

template<class LO, class SC>
SC BEIntegratorHelmholtz<LO, SC>::evalShiftedKernelNeumann(
  const SCVT* x,
  const SCVT* y,
  const SCVT* n
  ) const {

  SCVT norm = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) + ( x[1] - y[1] )*( x[1] - y[1] ) + ( x[2] - y[2] )*( x[2] - y[2] ) );
  SCVT dot = ( x[0] - y[0] ) * n[0] + ( x[1] - y[1] ) * n[1] + ( x[2] - y[2] ) * n[2];
  SC i( 0.0, 1.0 );
  return (SCVT) PI_FACT * ( dot / ( norm * norm * norm ) ) *
    ( ( (SCVT) 1.0 - i * kappa * norm ) * std::exp( i * kappa * norm )
    - (SCVT) 1.0 );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::evalShiftedKernelNeumann(
  SCVT x1,
  SCVT x2,
  SCVT x3,
  SCVT y1,
  SCVT y2,
  SCVT y3,
  SCVT n1,
  SCVT n2,
  SCVT n3,
  SCVT & real,
  SCVT & imag
  ) const {

  SCVT diff1 = x1 - y1;
  SCVT diff2 = x2 - y2;
  SCVT diff3 = x3 - y3;

  SCVT rekappa = std::real( kappa );
  SCVT imkappa = std::imag( kappa );
  SCVT norm = std::sqrt( diff1 * diff1 + diff2 * diff2 + diff3 * diff3 );
  SCVT dot = diff1 * n1 + diff2 * n2 + diff3 * n3;
  SCVT mult = (SCVT) PI_FACT * ( dot / ( norm * norm * norm ) );
  SCVT rekappan = rekappa * norm;
  SCVT imkappan = imkappa * norm;
  SCVT expim = std::exp( -imkappan );
  SCVT sine = std::sin( rekappan );
  SCVT cosine = std::cos( rekappan );

  real = mult * ( expim
    * ( ( (SCVT) 1.0 + imkappan ) * cosine + rekappan * sine ) - (SCVT) 1.0 );
  imag = mult * expim
    * ( ( (SCVT) 1.0 + imkappan ) * sine - rekappan * cosine );
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::representationFormulaP1P0(
  const SCVT *xCoord,
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
  SC dlkernel;
  SC integrand;
  SC onePatch;

  SC * shiftedKernel = new SC[ quadSizes[ this->quadratureOrder[ 1 ] ] ];
  SC * shiftedKernelNeumann =
    new SC[ quadSizes[ this->quadratureOrder[ 1 ] ] ];
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

        for ( int k = 0; k < quadSizes[ this->quadratureOrder[ 1 ] ]; ++k ) {
          shiftedKernel[ k ] = evalShiftedKernel( xSingle,
            ( quadratureNodes + 3 * k ) );
          shiftedKernelNeumann[ k ] = evalShiftedKernelNeumann( xSingle,
            ( quadratureNodes + 3 * k ), n );
        }

        for ( int rot = 0; rot < 3; ++rot ) {

          this->getLocalCoordinates( x1, x2, x3, r1, r2, n, stau, alpha1,
            alpha2, rot );
          this->getLocalQuadratureNode( x1, x2, x3, xSingle, r1, r2, n,
            xLocal, rot );

          if ( rot == 0 ) {
            val += neu.get( j ) *
              collocation1LayerP0( xLocal, stau, alpha1, alpha2,
              this->quadratureOrder[ 1 ], shiftedKernel, j );
          }

          val -= dir.get( ind[ rot ] ) *
            collocation2LayerP1( xLocal, stau, alpha1, alpha2,
            this->quadratureOrder[ 1 ], shiftedKernelNeumann, j, rot );
        }

      } else {

        onePatch = 0.0;

        for ( int k = 0; k < quadSizes[ this->quadratureOrder[ 1 ] ]; ++k ) {

          intPx = qPointsIn[ 2 * k ];
          intPy = qPointsIn[ 2 * k + 1 ];

          integrand = neu.get( j ) *
            evalSingleLayerKernel( xSingle, ( quadratureNodes + 3 * k ) );

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

  delete [] shiftedKernel;
  delete [] shiftedKernelNeumann;
  delete [] quadratureNodes;
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::doubleLayerPotentialP1(
  const SCVT * xCoord,
  LO numPoints,
  const Vector<LO, SC> & density,
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
  SC dlkernel;
  SC integrand;
  SC onePatch;
  SC * shiftedKernelNeumann =
    new SC[ quadSizes[ this->quadratureOrder[ 1 ] ] ];
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

        for ( int k = 0; k < quadSizes[ this->quadratureOrder[ 1 ] ]; ++k ) {
          shiftedKernelNeumann[ k ] = evalShiftedKernelNeumann( xSingle,
            ( quadratureNodes + 3 * k ), n );
        }

        for ( int rot = 0; rot < 3; ++rot ) {

          this->getLocalCoordinates( x1, x2, x3, r1, r2, n, stau, alpha1,
            alpha2, rot );
          this->getLocalQuadratureNode( x1, x2, x3, xSingle, r1, r2, n,
            xLocal, rot );

          val += density.get( ind[ rot ] ) *
            collocation2LayerP1( xLocal, stau, alpha1, alpha2,
            this->quadratureOrder[ 1 ], shiftedKernelNeumann, j, rot );
        }

      } else {

        onePatch = 0.0;

        for ( int k = 0; k < quadSizes[ this->quadratureOrder[ 1 ] ]; ++k ) {

          intPx = qPointsIn[ 2 * k ];
          intPy = qPointsIn[ 2 * k + 1 ];

          dlkernel = evalDoubleLayerKernel( xSingle,
            ( quadratureNodes + 3 * k ), n );
          integrand = density.get( ind[ 0 ] ) * dlkernel *
            ( (SCVT) 1.0 - intPx - intPy );
          integrand += density.get( ind[ 1 ] ) * dlkernel * intPx;
          integrand += density.get( ind[ 2 ] ) * dlkernel * intPy;

          onePatch += integrand * (SCVT)
            quadWeights[ this->quadratureOrder[ 1 ] ][ k ];
        }

        val += onePatch * this->space->getMesh( )->getElemArea( j );
      }
    }

    values.set( i, val );
  }

  delete [] shiftedKernelNeumann;
  delete [] quadratureNodes;
}

template<class LO, class SC>
void BEIntegratorHelmholtz<LO, SC>::singleLayerPotentialP0(
  const SCVT * xCoord,
  LO numPoints,
  const Vector<LO, SC> & density,
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
  SC integrand;
  SC onePatch;
  SC * shiftedKernel = new SC[ quadSizes[ this->quadratureOrder[ 1 ] ] ];
  SCVT * quadratureNodes =
    new SCVT[ quadSizes[ this->quadratureOrder[ 1 ] ] * 3 ];

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

        for ( int k = 0; k < quadSizes[ this->quadratureOrder[ 1 ] ]; ++k ) {
          shiftedKernel[ k ] = evalShiftedKernel( xSingle,
            ( quadratureNodes + 3 * k ) );
        }

        this->getLocalCoordinates( x1, x2, x3, r1, r2, n, stau, alpha1,
          alpha2, 0 );
        this->getLocalQuadratureNode( x1, x2, x3, xSingle, r1, r2, n,
          xLocal, 0 );

        val += density.get( j ) *
          collocation1LayerP0( xLocal, stau, alpha1, alpha2,
          this->quadratureOrder[ 1 ], shiftedKernel, j );

      } else {

        onePatch = 0.0;

        for ( int k = 0; k < quadSizes[ this->quadratureOrder[ 1 ] ]; ++k ) {

          integrand = density.get( j ) *
            evalSingleLayerKernel( xSingle, ( quadratureNodes + 3 * k ) );

          onePatch += integrand *
            (SCVT) quadWeights[ this->quadratureOrder[ 1 ] ][ k ];
        }

        val += onePatch * this->space->getMesh( )->getElemArea( j );
      }
    }

    values.set( i, val );
  }

  delete [] shiftedKernel;
  delete [] quadratureNodes;
}

}

#endif
