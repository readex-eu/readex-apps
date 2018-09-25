/*!
 * @file    BEIntegrator.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    July 10, 2013
 * @brief   Header file for pure virtual class BEIntegrator
 * 
 */

#ifndef BEINTEGRATOR_H
#define BEINTEGRATOR_H

#define EUC_SP(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

#include "BESpace.h"
#include "FullMatrix.h"
#include "Quadratures.h"
#include "Macros.h"
#include "Settings.h"
#include "BEBilinearForm.h"

#include <algorithm>

namespace bem4i {

/*!
 * abstract base class for integrators over triangular elements
 * 
 * the class uses Curiously recurring template pattern to replace virtual methods
 */
template<class LO, class SC, class SpecificIntegrator>
class BEIntegrator {
  // to get inner type of complex numbers (for Helmholtz)
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

#pragma omp declare reduction\
( + : std::complex< double >, std::complex< float > : omp_out += omp_in )\
initializer ( omp_priv = 0.0 )

  //! default constructor
  BEIntegrator( );

  //! copy constructor
  BEIntegrator( const BEIntegrator& orig );

  //! constructor taking BESpace as the argument
  BEIntegrator( BESpace<LO, SC>* space );

  //! destructor
  virtual ~BEIntegrator( );

  void initDisjointQuadratureData(
    int * quadratureOrderDisjointElems
    ) {

    int outerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 0 ] ];
    int innerPoints = quadSizes[ this->quadratureOrderDisjointElems[ 1 ] ];

    this->x1d = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
    this->x2d = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
    this->x3d = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
    this->y1d = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
    this->y2d = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
    this->y3d = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );

    this->wxyd = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
    this->x1dref = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
    this->x2dref = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
    this->y1dref = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );
    this->y2dref = (SCVT *) _mm_malloc(
      outerPoints * innerPoints * sizeof ( SCVT ), DATA_ALIGN );

    int counter = 0;
    for ( int i = 0; i < outerPoints; ++i ) {
      for ( int j = 0; j < innerPoints; ++j ) {
        this->wxyd[ counter ] =
          (SCVT) quadWeights[ this->quadratureOrderDisjointElems[ 0 ] ][ i ]
          * (SCVT) quadWeights[ this->quadratureOrderDisjointElems[ 1 ] ][ j ];
        this->x1dref[ counter ] =
          (SCVT) quadPoints_X1[ this->quadratureOrderDisjointElems[ 0 ] ][ i ];
        this->x2dref[ counter ] =
          (SCVT) quadPoints_X2[ this->quadratureOrderDisjointElems[ 0 ] ][ i ];
        this->y1dref[ counter ] =
          (SCVT) quadPoints_X1[ this->quadratureOrderDisjointElems[ 1 ] ][ j ];
        this->y2dref[ counter ] =
          (SCVT) quadPoints_X2[ this->quadratureOrderDisjointElems[ 1 ] ][ j ];
        ++counter;
      }
    }
  }

  void initSauterSchwabQuadratureData(
    int * quadratureOrder
    ) {

    int qOrder0 = quadratureOrder[ 0 ];
    int qOrder1 = quadratureOrder[ 1 ];
    int qOrder2 = quadratureOrder[ 2 ];
    int qOrder3 = quadratureOrder[ 3 ];

    int qSize0 = lineQuadSizes[ qOrder0 ];
    int qSize1 = lineQuadSizes[ qOrder1 ];
    int qSize2 = lineQuadSizes[ qOrder2 ];
    int qSize3 = lineQuadSizes[ qOrder3 ];

    int totalSize = qSize0 * qSize1 * qSize2 * qSize3;

    for ( int i = 0; i < 6; ++i ) {
      this->x1RefIdentical[ i ] =
        (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
      this->x2RefIdentical[ i ] =
        (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
      this->y1RefIdentical[ i ] =
        (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
      this->y2RefIdentical[ i ] =
        (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
      this->jacobianWeightIdentical[ i ] =
        (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
    }

    for ( int i = 0; i < 5; ++i ) {
      this->x1RefCommonEdge[ i ] =
        (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
      this->x2RefCommonEdge[ i ] =
        (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
      this->y1RefCommonEdge[ i ] =
        (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
      this->y2RefCommonEdge[ i ] =
        (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
      this->jacobianWeightCommonEdge[ i ] =
        (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
    }

    for ( int i = 0; i < 2; ++i ) {
      this->x1RefCommonVertex[ i ] =
        (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
      this->x2RefCommonVertex[ i ] =
        (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
      this->y1RefCommonVertex[ i ] =
        (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
      this->y2RefCommonVertex[ i ] =
        (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
      this->jacobianWeightCommonVertex[ i ] =
        (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
    }

    this->x1RefDisjoint[ 0 ] =
      (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
    this->x2RefDisjoint[ 0 ] =
      (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
    this->y1RefDisjoint[ 0 ] =
      (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
    this->y2RefDisjoint[ 0 ] =
      (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
    this->jacobianWeightDisjoint[ 0 ] =
      (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );

    this->x1ss = (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
    this->x2ss = (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
    this->x3ss = (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
    this->y1ss = (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
    this->y2ss = (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );
    this->y3ss = (SCVT *) _mm_malloc( totalSize * sizeof ( SCVT ), DATA_ALIGN );

    SCVT w[ 4 ];
    SCVT xRef[ 2 ];
    SCVT yRef[ 2 ];
    SCVT jacobian;
    SCVT weight;

    int counter = 0;
    for ( int dksi = 0; dksi < qSize0; ++dksi ) {
      w[ 0 ] = lineQuadPoints[ qOrder0 ][ dksi ];
      for ( int deta3 = 0; deta3 < qSize1; ++deta3 ) {
        w[ 3 ] = lineQuadPoints[ qOrder1 ][ deta3 ];
        for ( int deta2 = 0; deta2 < qSize2; ++deta2 ) {
          w[ 2 ] = lineQuadPoints[ qOrder2 ][ deta2 ];
          for ( int deta1 = 0; deta1 < qSize3; ++deta1 ) {
            w[ 1 ] = lineQuadPoints[ qOrder3 ][ deta1 ];

            weight = lineQuadWeights[ qOrder0 ][ dksi ] *
              lineQuadWeights[ qOrder3 ][ deta1 ] *
              lineQuadWeights[ qOrder2 ][ deta2 ] *
              lineQuadWeights[ qOrder1 ][ deta3 ];

            for ( int simplex = 0; simplex < 6; ++simplex ) {
              this->cube2triIdentical( w, simplex, xRef, yRef, jacobian );
              this->x1RefIdentical[ simplex ][ counter ] = xRef[ 0 ];
              this->x2RefIdentical[ simplex ][ counter ] = xRef[ 1 ];
              this->y1RefIdentical[ simplex ][ counter ] = yRef[ 0 ];
              this->y2RefIdentical[ simplex ][ counter ] = yRef[ 1 ];
              this->jacobianWeightIdentical[ simplex ][ counter ] =
                jacobian * weight;
            }

            for ( int simplex = 0; simplex < 5; ++simplex ) {
              this->cube2triCommonEdge( w, simplex, xRef, yRef, jacobian );
              this->x1RefCommonEdge[ simplex ][ counter ] = xRef[ 0 ];
              this->x2RefCommonEdge[ simplex ][ counter ] = xRef[ 1 ];
              this->y1RefCommonEdge[ simplex ][ counter ] = yRef[ 0 ];
              this->y2RefCommonEdge[ simplex ][ counter ] = yRef[ 1 ];
              this->jacobianWeightCommonEdge[ simplex ][ counter ] =
                jacobian * weight;
            }

            for ( int simplex = 0; simplex < 2; ++simplex ) {
              this->cube2triCommonVertex( w, simplex, xRef, yRef, jacobian );
              this->x1RefCommonVertex[ simplex ][ counter ] = xRef[ 0 ];
              this->x2RefCommonVertex[ simplex ][ counter ] = xRef[ 1 ];
              this->y1RefCommonVertex[ simplex ][ counter ] = yRef[ 0 ];
              this->y2RefCommonVertex[ simplex ][ counter ] = yRef[ 1 ];
              this->jacobianWeightCommonVertex[ simplex ][ counter ] =
                jacobian * weight;
            }

            this->cube2triDisjoint( w, xRef, yRef, jacobian );
            this->x1RefDisjoint[ 0 ][ counter ] = xRef[ 0 ];
            this->x2RefDisjoint[ 0 ][ counter ] = xRef[ 1 ];
            this->y1RefDisjoint[ 0 ][ counter ] = yRef[ 0 ];
            this->y2RefDisjoint[ 0 ][ counter ] = yRef[ 1 ];
            this->jacobianWeightDisjoint[ 0 ][ counter ] = jacobian * weight;

            ++counter;
          }
        }
      }
    }
  }

  void initSteinbachQuadratureData(
    int * quadratureOrder
    ) {

    this->coordCoeff = (SCVT) 4.0;

    int qOrderX = this->quadratureOrder[ 0 ];
    int qOrderY = this->quadratureOrder[ 1 ];
    int qSizeX = quadSizes[ qOrderX ];
    int qSizeY = quadSizes[ qOrderY ];

    const int width = DATA_WIDTH;
    const int align = DATA_ALIGN;
    int qSizeXY = qSizeX * qSizeY;
    this->qSizeXpad = qSizeX + width - qSizeX % width;
    this->qSizeXYpad = qSizeXY + width - qSizeXY % width;

    // outer quadrature + analytic collocation
    this->x1stref =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->x2stref =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->x1st =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->x2st =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->x3st =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->sx =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->tx =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->ux =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->wxst =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->r11 =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->r12 =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->r13 =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->r21 =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->r22 =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->r23 =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->alpha1 =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->alpha2 =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->stau =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->ref =
      (SCVT *) _mm_malloc( this->qSizeXpad * sizeof ( SCVT ), align );
    this->rot =
      (int *) _mm_malloc( this->qSizeXpad * sizeof ( int ), align );

    // outer + inner quadrature for regular part
    this->x1strefCol =
      (SCVT *) _mm_malloc( this->qSizeXYpad * sizeof ( SCVT ), align );
    this->x2strefCol =
      (SCVT *) _mm_malloc( this->qSizeXYpad * sizeof ( SCVT ), align );
    this->x1stCol =
      (SCVT *) _mm_malloc( this->qSizeXYpad * sizeof ( SCVT ), align );
    this->x2stCol =
      (SCVT *) _mm_malloc( this->qSizeXYpad * sizeof ( SCVT ), align );
    this->x3stCol =
      (SCVT *) _mm_malloc( this->qSizeXYpad * sizeof ( SCVT ), align );
    this->y1strefCol =
      (SCVT *) _mm_malloc( this->qSizeXYpad * sizeof ( SCVT ), align );
    this->y2strefCol =
      (SCVT *) _mm_malloc( this->qSizeXYpad * sizeof ( SCVT ), align );
    this->y1stCol =
      (SCVT *) _mm_malloc( this->qSizeXYpad * sizeof ( SCVT ), align );
    this->y2stCol =
      (SCVT *) _mm_malloc( this->qSizeXYpad * sizeof ( SCVT ), align );
    this->y3stCol =
      (SCVT *) _mm_malloc( this->qSizeXYpad * sizeof ( SCVT ), align );
    this->wxyst =
      (SCVT *) _mm_malloc( this->qSizeXYpad * sizeof ( SCVT ), align );

    int counter = 0;

    for ( int i = 0; i < qSizeX; ++i ) {

      this->x1stref[ i ] = (SCVT) quadPoints_X1[ qOrderX ][ i ];
      this->x2stref[ i ] = (SCVT) quadPoints_X2[ qOrderX ][ i ];
      this->wxst[ i ] = (SCVT) quadWeights[ qOrderX ][ i ];

      for ( int j = 0; j < qSizeY; ++j ) {
        this->x1strefCol[ counter ] = (SCVT) quadPoints_X1[ qOrderX ][ i ];
        this->x2strefCol[ counter ] = (SCVT) quadPoints_X2[ qOrderX ][ i ];
        this->y1strefCol[ counter ] = (SCVT) quadPoints_X1[ qOrderY ][ j ];
        this->y2strefCol[ counter ] = (SCVT) quadPoints_X2[ qOrderY ][ j ];
        this->wxyst[ counter ] =
          (SCVT) ( quadWeights[ qOrderX ][ i ] * quadWeights[ qOrderY ][ j ] );

        ++counter;
      }
    }

    for ( int i = qSizeX; i < this->qSizeXpad; ++i ) {
      this->x1stref[ i ] = this->x1stref[ qSizeX - 1 ];
      this->x2stref[ i ] = this->x2stref[ qSizeX - 1 ];
      this->wxst[ i ] = (SCVT) 0.0;
    }

    for ( int i = qSizeXY; i < this->qSizeXYpad; ++i ) {
      this->x1strefCol[ i ] = this->x1strefCol[ qSizeXY - 1 ];
      this->x2strefCol[ i ] = this->x2strefCol[ qSizeXY - 1 ];
      this->y1strefCol[ i ] = this->y1strefCol[ qSizeXY - 1 ];
      this->y2strefCol[ i ] = this->y2strefCol[ qSizeXY - 1 ];
      this->wxyst[ i ] = (SCVT) 0.0;
    }
  }

  /*! returns an element matrix of the single layer potential
   *
   * @param[in]       outerElem index of the outer integration triangle
   * @param[in]       innerElem index of the inner integration triangle
   * @param[in,out]   element matrix
   */
  void getElemMatrix1Layer(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) {
    static_cast<SpecificIntegrator*> ( this )->computeElemMatrix1Layer(
      outerElem, innerElem, matrix );
  }

  /*! returns an element matrix of the double layer potential
   *
   * @param[in]       outerElem index of the outer integration triangle
   * @param[in]       innerElem index of the inner integration triangle
   * @param[in,out]   element matrix
   */
  void getElemMatrix2Layer(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) {
    static_cast<SpecificIntegrator*> ( this )->computeElemMatrix2Layer(
      outerElem, innerElem, matrix );
  }

  /*! returns an element matrix of the hypersingular operator
   *
   * @param[in]       outerElem index of the outer integration triangle
   * @param[in]       innerElem index of the inner integration triangle
   * @param[in,out]   element matrix
   */
  void getElemMatrixHypersingular(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) {
    static_cast<SpecificIntegrator*> ( this )->computeElemMatrixHypersingular(
      outerElem, innerElem, matrix );
  }

  inline BESpace<LO, SC>* getSpace( ) const {
    return this->space;
  }


#if N_MIC > 0

  __attribute__ ( ( target( mic ) ) )
#endif
  static void getQuadratureNodesMIC(
    const SCVT* x1,
    const SCVT* x2,
    const SCVT* x3,
    int quadratureOrder,
    SCVT* nodes
    ) {
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

#if N_MIC > 0

  __attribute__ ( ( target( mic ) ) )
#endif
  static void getQuadratureNodesMIC(
    const SCVT * ox1,
    const SCVT * ox2,
    const SCVT * ox3,
    const SCVT * ix1,
    const SCVT * ix2,
    const SCVT * ix3,
    int quadratureOrderOuter,
    int quadratureOrderInner,
    SCVT * outerNodes,
    SCVT * innerNodes
    ) {

    int numPointsOuter = quadSizes[ quadratureOrderOuter ];
    int numPointsInner = quadSizes[ quadratureOrderInner ];

    __assume_aligned( outerNodes, DATA_ALIGN );
    __assume_aligned( innerNodes, DATA_ALIGN );

    int counter = 0;
    for ( int i = 0; i < numPointsOuter; ++i ) {
      for ( int j = 0; j < numPointsInner; ++j ) {
        outerNodes[ counter * 3 ] = ox1[ 0 ] + ( ox2[ 0 ] - ox1[ 0 ] ) *
          quadPoints[ quadratureOrderOuter ][ i * 2 ] +
          ( ox3[ 0 ] - ox1[ 0 ] ) *
          quadPoints[ quadratureOrderOuter ][ i * 2 + 1 ];
        outerNodes[ counter * 3 + 1 ] = ox1[ 1 ] + ( ox2[ 1 ] - ox1[ 1 ] ) *
          quadPoints[ quadratureOrderOuter ][ i * 2 ] +
          ( ox3[ 1 ] - ox1[ 1 ] ) *
          quadPoints[ quadratureOrderOuter ][ i * 2 + 1 ];
        outerNodes[ counter * 3 + 2 ] = ox1[ 2 ] + ( ox2[ 2 ] - ox1[ 2 ] ) *
          quadPoints[ quadratureOrderOuter ][ i * 2 ] +
          ( ox3[ 2 ] - ox1[ 2 ] ) *
          quadPoints[ quadratureOrderOuter ][ i * 2 + 1 ];

        innerNodes[ counter * 3 ] = ix1[ 0 ] + ( ix2[ 0 ] - ix1[ 0 ] ) *
          quadPoints[ quadratureOrderInner ][ j * 2 ] +
          ( ix3[ 0 ] - ix1[ 0 ] ) *
          quadPoints[ quadratureOrderInner ][ j * 2 + 1 ];
        innerNodes[ counter * 3 + 1 ] = ix1[ 1 ] + ( ix2[ 1 ] - ix1[ 1 ] ) *
          quadPoints[ quadratureOrderInner ][ j * 2 ] +
          ( ix3[ 1 ] - ix1[ 1 ] ) *
          quadPoints[ quadratureOrderInner ][ j * 2 + 1 ];
        innerNodes[ counter * 3 + 2 ] = ix1[ 2 ] + ( ix2[ 2 ] - ix1[ 2 ] ) *
          quadPoints[ quadratureOrderInner ][ j * 2 ] +
          ( ix3[ 2 ] - ix1[ 2 ] ) *
          quadPoints[ quadratureOrderInner ][ j * 2 + 1 ];

        ++counter;
      }
    }
  }

#if N_MIC > 0

  __attribute__ ( ( target( mic ) ) )
#endif
  static void getQuadratureNodesMIC(
    const SCVT * ox1,
    const SCVT * ox2,
    const SCVT * ox3,
    const SCVT * ix1,
    const SCVT * ix2,
    const SCVT * ix3,
    int quadratureOrderOuter,
    int quadratureOrderInner,
    const SCVT * outerX1ref,
    const SCVT * outerX2ref,
    const SCVT * innerX1ref,
    const SCVT * innerX2ref,
    SCVT * outerX1,
    SCVT * outerX2,
    SCVT * outerX3,
    SCVT * innerX1,
    SCVT * innerX2,
    SCVT * innerX3
    ) {

    __assume_aligned( outerX1, DATA_ALIGN );
    __assume_aligned( outerX2, DATA_ALIGN );
    __assume_aligned( outerX3, DATA_ALIGN );
    __assume_aligned( innerX1, DATA_ALIGN );
    __assume_aligned( innerX2, DATA_ALIGN );
    __assume_aligned( innerX3, DATA_ALIGN );
    __assume_aligned( outerX1ref, DATA_ALIGN );
    __assume_aligned( outerX2ref, DATA_ALIGN );
    __assume_aligned( innerX1ref, DATA_ALIGN );
    __assume_aligned( innerX2ref, DATA_ALIGN );

    int numPointsOuter = quadSizes[ quadratureOrderOuter ];
    int numPointsInner = quadSizes[ quadratureOrderInner ];

    int cmax = numPointsOuter * numPointsInner;
    int c;

#pragma omp simd linear( c : 1 ) simdlen( DATA_WIDTH )
    for ( c = 0; c < cmax; ++c ) {
      outerX1[ c ] = ox1[ 0 ] + ( ox2[ 0 ] - ox1[ 0 ] ) *
        outerX1ref[ c ] + ( ox3[ 0 ] - ox1[ 0 ] ) *
        outerX2ref[ c ];
      outerX2[ c ] = ox1[ 1 ] + ( ox2[ 1 ] - ox1[ 1 ] ) *
        outerX1ref[ c ] + ( ox3[ 1 ] - ox1[ 1 ] ) *
        outerX2ref[ c ];
      outerX3[ c ] = ox1[ 2 ] + ( ox2[ 2 ] - ox1[ 2 ] ) *
        outerX1ref[ c ] + ( ox3[ 2 ] - ox1[ 2 ] ) *
        outerX2ref[ c ];

      innerX1[ c ] = ix1[ 0 ] + ( ix2[ 0 ] - ix1[ 0 ] ) *
        innerX1ref[ c ] + ( ix3[ 0 ] - ix1[ 0 ] ) *
        innerX2ref[ c ];
      innerX2[ c ] = ix1[ 1 ] + ( ix2[ 1 ] - ix1[ 1 ] ) *
        innerX1ref[ c ] + ( ix3[ 1 ] - ix1[ 1 ] ) *
        innerX2ref[ c ];
      innerX3[ c ] = ix1[ 2 ] + ( ix2[ 2 ] - ix1[ 2 ] ) *
        innerX1ref[ c ] + ( ix3[ 2 ] - ix1[ 2 ] ) *
        innerX2ref[ c ];
    }

  }

#if N_MIC > 0
  __attribute__ ( ( target( mic ) ) )
#endif
  static bool areElementsDisjointOnMIC(
    LO outerElem,
    LO innerElem,
    const SCVT * nodes,
    const LO * elems,
    SCVT threshold = 2.0
    );

protected:

  BESpace<LO, SC> * space;

  static int mod3arr [ 5 ];

  //! quadrature rule
  int * quadratureOrder;
  //! quadrature order for disjoint elements
  int * quadratureOrderDisjointElems;
  //! type of quadrature (Steinbach or SauterSchwab)
  quadratureType quadrature;

  // remote quadrature data
  SCVT * x1d;
  SCVT * x2d;
  SCVT * x3d;
  SCVT * y1d;
  SCVT * y2d;
  SCVT * y3d;
  SCVT * x1dref;
  SCVT * x2dref;
  SCVT * y1dref;
  SCVT * y2dref;
  SCVT * wxyd;

  // Steinbach data
  // outer quadrature + analytic
  SCVT * x1stref;
  SCVT * x2stref;
  SCVT * wxst;
  SCVT * x1st;
  SCVT * x2st;
  SCVT * x3st;
  SCVT * sx;
  SCVT * tx;
  SCVT * ux;
  SCVT * r11;
  SCVT * r12;
  SCVT * r13;
  SCVT * r21;
  SCVT * r22;
  SCVT * r23;
  SCVT * alpha1;
  SCVT * alpha2;
  SCVT * stau;
  SCVT * ref;
  int * rot;
  int qSizeXpad;
  int qSizeXYpad;
  SCVT coordCoeff;
  // outer+inner quadrature for regular part
  SCVT * x1strefCol;
  SCVT * x2strefCol;
  SCVT * x1stCol;
  SCVT * x2stCol;
  SCVT * x3stCol;
  SCVT * y1strefCol;
  SCVT * y2strefCol;
  SCVT * y1stCol;
  SCVT * y2stCol;
  SCVT * y3stCol;
  SCVT * wxyst;

  // Sauter Schwab data
  SCVT ** x1RefIdentical;
  SCVT ** x2RefIdentical;
  SCVT ** y1RefIdentical;
  SCVT ** y2RefIdentical;
  SCVT ** jacobianWeightIdentical;
  SCVT ** x1RefCommonEdge;
  SCVT ** x2RefCommonEdge;
  SCVT ** y1RefCommonEdge;
  SCVT ** y2RefCommonEdge;
  SCVT ** jacobianWeightCommonEdge;
  SCVT ** x1RefCommonVertex;
  SCVT ** x2RefCommonVertex;
  SCVT ** y1RefCommonVertex;
  SCVT ** y2RefCommonVertex;
  SCVT ** jacobianWeightCommonVertex;
  SCVT ** x1RefDisjoint;
  SCVT ** x2RefDisjoint;
  SCVT ** y1RefDisjoint;
  SCVT ** y2RefDisjoint;
  SCVT ** jacobianWeightDisjoint;
  SCVT * x1ss;
  SCVT * x2ss;
  SCVT * x3ss;
  SCVT * y1ss;
  SCVT * y2ss;
  SCVT * y3ss;

  inline int mod3(
    int ind
    ) const {
    return this->mod3arr[ ind ];
  }

  //! returns values of linear functions in integration points
  void getLinValues(
    int order,
    SCVT* linValues
    ) const;

  //! returns integration nodes for the triangle defined by x1, x2, x3 to a user-preallocated array nodes
  void getQuadratureNodes(
    const SCVT* x1,
    const SCVT* x2,
    const SCVT* x3,
    int quadratureOrder,
    SCVT* nodes
    ) const;

  void updateDisjointQuadratureNodes(
    const SCVT * ox1,
    const SCVT * ox2,
    const SCVT * ox3,
    const SCVT * ix1,
    const SCVT * ix2,
    const SCVT * ix3
    ) const {

    int numPointsOuter = quadSizes[ this->quadratureOrderDisjointElems[ 0 ] ];
    int numPointsInner = quadSizes[ this->quadratureOrderDisjointElems[ 1 ] ];

    int cmax = numPointsOuter * numPointsInner;
    int c = 0;

    SCVT * x1d = this->x1d;
    SCVT * x2d = this->x2d;
    SCVT * x3d = this->x3d;
    SCVT * y1d = this->y1d;
    SCVT * y2d = this->y2d;
    SCVT * y3d = this->y3d;
    SCVT * x1dref = this->x1dref;
    SCVT * x2dref = this->x2dref;
    SCVT * y1dref = this->y1dref;
    SCVT * y2dref = this->y2dref;

    const int align = DATA_ALIGN;
    const int width = DATA_WIDTH;

#pragma omp simd \
linear( c : 1 ) \
aligned( x1d, x2d, x3d, x1dref, x2dref : align ) \
simdlen( width )
    for ( c = 0; c < cmax; ++c ) {
      x1d[ c ] = ox1[ 0 ] + ( ox2[ 0 ] - ox1[ 0 ] ) *
        x1dref[ c ] + ( ox3[ 0 ] - ox1[ 0 ] ) *
        x2dref[ c ];

      x2d[ c ] = ox1[ 1 ] + ( ox2[ 1 ] - ox1[ 1 ] ) *
        x1dref[ c ] + ( ox3[ 1 ] - ox1[ 1 ] ) *
        x2dref[ c ];

      x3d[ c ] = ox1[ 2 ] + ( ox2[ 2 ] - ox1[ 2 ] ) *
        x1dref[ c ] + ( ox3[ 2 ] - ox1[ 2 ] ) *
        x2dref[ c ];
    }

#pragma omp simd \
linear( c : 1 ) \
aligned( y1d, y2d, y3d, y1dref, y2dref : align ) \
simdlen( width )
    for ( c = 0; c < cmax; ++c ) {
      y1d[ c ] = ix1[ 0 ] + ( ix2[ 0 ] - ix1[ 0 ] ) *
        y1dref[ c ] + ( ix3[ 0 ] - ix1[ 0 ] ) *
        y2dref[ c ];

      y2d[ c ] = ix1[ 1 ] + ( ix2[ 1 ] - ix1[ 1 ] ) *
        y1dref[ c ] + ( ix3[ 1 ] - ix1[ 1 ] ) *
        y2dref[ c ];

      y3d[ c ] = ix1[ 2 ] + ( ix2[ 2 ] - ix1[ 2 ] ) *
        y1dref[ c ] + ( ix3[ 2 ] - ix1[ 2 ] ) *
        y2dref[ c ];
    }
  }

  void updateSteinbachRegularQuadratureNodes(
    const SCVT * x1,
    const SCVT * x2,
    const SCVT * x3,
    const SCVT * y1,
    const SCVT * y2,
    const SCVT * y3
    ) const {

    SCVT * x1stCol = this->x1stCol;
    SCVT * x2stCol = this->x2stCol;
    SCVT * x3stCol = this->x3stCol;
    SCVT * x1refCol = this->x1strefCol;
    SCVT * x2refCol = this->x2strefCol;
    SCVT * y1stCol = this->y1stCol;
    SCVT * y2stCol = this->y2stCol;
    SCVT * y3stCol = this->y3stCol;
    SCVT * y1refCol = this->y1strefCol;
    SCVT * y2refCol = this->y2strefCol;

    int c;
    const int width = DATA_WIDTH;
    const int align = DATA_ALIGN;

#pragma omp simd \
linear( c : 1 )   \
aligned( x1stCol, x2stCol, x3stCol, x1refCol, x2refCol : align ) \
simdlen( width )
    for ( c = 0; c < this->qSizeXYpad; ++c ) {
      x1stCol[ c ] = x1[ 0 ]
        + ( x2[ 0 ] - x1[ 0 ] ) * x1refCol[ c ]
        + ( x3[ 0 ] - x1[ 0 ] ) * x2refCol[ c ];

      x2stCol[ c ] = x1[ 1 ]
        + ( x2[ 1 ] - x1[ 1 ] ) * x1refCol[ c ]
        + ( x3[ 1 ] - x1[ 1 ] ) * x2refCol[ c ];

      x3stCol[ c ] = x1[ 2 ]
        + ( x2[ 2 ] - x1[ 2 ] ) * x1refCol[ c ]
        + ( x3[ 2 ] - x1[ 2 ] ) * x2refCol[ c ];
    }

#pragma omp simd \
linear( c : 1 )   \
aligned( y1stCol, y2stCol, y3stCol, y1refCol, y2refCol : align ) \
simdlen( width )
    for ( c = 0; c < this->qSizeXYpad; ++c ) {
      y1stCol[ c ] = y1[ 0 ]
        + ( y2[ 0 ] - y1[ 0 ] ) * y1refCol[ c ]
        + ( y3[ 0 ] - y1[ 0 ] ) * y2refCol[ c ];

      y2stCol[ c ] = y1[ 1 ]
        + ( y2[ 1 ] - y1[ 1 ] ) * y1refCol[ c ]
        + ( y3[ 1 ] - y1[ 1 ] ) * y2refCol[ c ];

      y3stCol[ c ] = y1[ 2 ]
        + ( y2[ 2 ] - y1[ 2 ] ) * y1refCol[ c ]
        + ( y3[ 2 ] - y1[ 2 ] ) * y2refCol[ c ];
    }
  }

  void updateSteinbachLocalQuadratureNodes(
    const SCVT * x1,
    const SCVT * x2,
    const SCVT * x3,
    const SCVT * y1,
    const SCVT * r1,
    const SCVT * r2,
    const SCVT * n
    ) const {

    SCVT * x1st = this->x1st;
    SCVT * x2st = this->x2st;
    SCVT * x3st = this->x3st;
    SCVT * x1stref = this->x1stref;
    SCVT * x2stref = this->x2stref;
    SCVT * sx = this->sx;
    SCVT * tx = this->tx;
    SCVT * ux = this->ux;

    int c;
    const int align = DATA_ALIGN;
    const int width = DATA_WIDTH;

#pragma omp simd \
aligned( x1st, x2st, x3st, x1stref, x2stref : align ) \
simdlen( width )
    for ( c = 0; c < this->qSizeXpad; ++c ) {
      x1st[ c ] = x1[ 0 ]
        + ( x2[ 0 ] - x1[ 0 ] ) * x1stref[ c ]
        + ( x3[ 0 ] - x1[ 0 ] ) * x2stref[ c ];

      x2st[ c ] = x1[ 1 ]
        + ( x2[ 1 ] - x1[ 1 ] ) * x1stref[ c ]
        + ( x3[ 1 ] - x1[ 1 ] ) * x2stref[ c ];

      x3st[ c ] = x1[ 2 ]
        + ( x2[ 2 ] - x1[ 2 ] ) * x1stref[ c ]
        + ( x3[ 2 ] - x1[ 2 ] ) * x2stref[ c ];
    }

#pragma omp simd \
aligned( x1st, x2st, x3st, sx, tx, ux : align ) \
simdlen( width )
    for ( c = 0; c < this->qSizeXpad; ++c ) {
      sx[ c ] = ( x1st[ c ] - y1[ 0 ] ) * r1[ 0 ]
        + ( x2st[ c ] - y1[ 1 ] ) * r1[ 1 ]
        + ( x3st[ c ] - y1[ 2 ] ) * r1[ 2 ];

      tx[ c ] = ( x1st[ c ] - y1[ 0 ] ) * r2[ 0 ]
        + ( x2st[ c ] - y1[ 1 ] ) * r2[ 1 ]
        + ( x3st[ c ] - y1[ 2 ] ) * r2[ 2 ];

      ux[ c ] = ( x1st[ c ] - y1[ 0 ] ) * n[ 0 ]
        + ( x2st[ c ] - y1[ 1 ] ) * n[ 1 ]
        + ( x3st[ c ] - y1[ 2 ] ) * n[ 2 ];
    }
  }

  void updateSteinbachLocalQuadratureNodes(
    const SCVT * x1,
    const SCVT * x2,
    const SCVT * x3,
    const SCVT * y1,
    const SCVT * y2,
    const SCVT * y3,
    const SCVT * r1,
    const SCVT * r2,
    const SCVT * n,
    int rot
    ) const {

    switch ( rot ) {
      case 0:
        this->updateSteinbachLocalQuadratureNodes( x1, x2, x3, y1, r1, r2, n );
        break;
      case 1:
        this->updateSteinbachLocalQuadratureNodes( x1, x2, x3, y2, r1, r2, n );
        break;
      case 2:
        this->updateSteinbachLocalQuadratureNodes( x1, x2, x3, y3, r1, r2, n );
        break;
    }
  }

  void updateSteinbachLocalQuadratureNodes(
    const SCVT * y1,
    const SCVT * r1,
    const SCVT * r2,
    const SCVT * n
    ) const {

    SCVT * x1st = this->x1st;
    SCVT * x2st = this->x2st;
    SCVT * x3st = this->x3st;
    SCVT * sx = this->sx;
    SCVT * tx = this->tx;
    SCVT * ux = this->ux;
    const int width = DATA_WIDTH;
    const int align = DATA_ALIGN;

    int c;

#pragma omp simd \
aligned( x1st, x2st, x3st, sx, tx, ux : align ) \
simdlen( width )
    for ( c = 0; c < this->qSizeXpad; ++c ) {
      sx[ c ] = ( x1st[ c ] - y1[ 0 ] ) * r1[ 0 ]
        + ( x2st[ c ] - y1[ 1 ] ) * r1[ 1 ]
        + ( x3st[ c ] - y1[ 2 ] ) * r1[ 2 ];

      tx[ c ] = ( x1st[ c ] - y1[ 0 ] ) * r2[ 0 ]
        + ( x2st[ c ] - y1[ 1 ] ) * r2[ 1 ]
        + ( x3st[ c ] - y1[ 2 ] ) * r2[ 2 ];

      ux[ c ] = ( x1st[ c ] - y1[ 0 ] ) * n[ 0 ]
        + ( x2st[ c ] - y1[ 1 ] ) * n[ 1 ]
        + ( x3st[ c ] - y1[ 2 ] ) * n[ 2 ];
    }
  }

  void finalizeSteinbachLocalQuadratureNodes(
    const SCVT * y1,
    const SCVT * y2,
    const SCVT * y3,
    const SCVT * n
    ) const {

    SCVT * x1st = this->x1st;
    SCVT * x2st = this->x2st;
    SCVT * x3st = this->x3st;
    SCVT * sx = this->sx;
    SCVT * tx = this->tx;
    SCVT * ux = this->ux;
    SCVT * r11 = this->r11;
    SCVT * r12 = this->r12;
    SCVT * r13 = this->r13;
    SCVT * r21 = this->r21;
    SCVT * r22 = this->r22;
    SCVT * r23 = this->r23;
    int * rot = this->rot;
    const int width = DATA_WIDTH;
    const int align = DATA_ALIGN;

    SCVT yrot11, yrot12, yrot13;
    int c;

#pragma omp simd \
private( yrot11, yrot12, yrot13 ) \
linear( c : 1 ) \
aligned( x1st, x2st, x3st, sx, tx, ux : align ) \
aligned( r11, r12, r13, r21, r22, r23, rot : align ) \
simdlen( width )
    for ( c = 0; c < this->qSizeXpad; ++c ) {
      if ( rot[ c ] == 1 ) {
        yrot11 = y2[ 0 ];
        yrot12 = y2[ 1 ];
        yrot13 = y2[ 2 ];
      } else if ( rot[ c ] == 2 ) {
        yrot11 = y3[ 0 ];
        yrot12 = y3[ 1 ];
        yrot13 = y3[ 2 ];
      } else { // rot[ c ] = 0
        yrot11 = y1[ 0 ];
        yrot12 = y1[ 1 ];
        yrot13 = y1[ 2 ];
      }

      sx[ c ] = ( x1st[ c ] - yrot11 ) * r11[ c ]
        + ( x2st[ c ] - yrot12 ) * r12[ c ]
        + ( x3st[ c ] - yrot13 ) * r13[ c ];

      tx[ c ] = ( x1st[ c ] - yrot11 ) * r21[ c ]
        + ( x2st[ c ] - yrot12 ) * r22[ c ]
        + ( x3st[ c ] - yrot13 ) * r23[ c ];

      ux[ c ] = ( x1st[ c ] - yrot11 ) * n[ 0 ]
        + ( x2st[ c ] - yrot12 ) * n[ 1 ]
        + ( x3st[ c ] - yrot13 ) * n[ 2 ];
    }
  }

  void updateSteinbachQuadratureNodes(
    const SCVT * x1,
    const SCVT * x2,
    const SCVT * x3
    ) const {

    SCVT * x1st = this->x1st;
    SCVT * x2st = this->x2st;
    SCVT * x3st = this->x3st;
    SCVT * x1stref = this->x1stref;
    SCVT * x2stref = this->x2stref;
    const int width = DATA_WIDTH;
    const int align = DATA_ALIGN;

    int c;

#pragma omp simd \
aligned( x1st, x2st, x3st, x1stref, x2stref : align ) \
simdlen( width )
    for ( c = 0; c < this->qSizeXpad; ++c ) {
      x1st[ c ] = x1[ 0 ]
        + ( x2[ 0 ] - x1[ 0 ] ) * x1stref[ c ]
        + ( x3[ 0 ] - x1[ 0 ] ) * x2stref[ c ];

      x2st[ c ] = x1[ 1 ]
        + ( x2[ 1 ] - x1[ 1 ] ) * x1stref[ c ]
        + ( x3[ 1 ] - x1[ 1 ] ) * x2stref[ c ];

      x3st[ c ] = x1[ 2 ]
        + ( x2[ 2 ] - x1[ 2 ] ) * x1stref[ c ]
        + ( x3[ 2 ] - x1[ 2 ] ) * x2stref[ c ];
    }

  }

  void updateSauterSchwabQuadratureNodes(
    int type,
    int simplex,
    const SCVT * x1,
    const SCVT * x2,
    const SCVT * x3,
    int xRot,
    const SCVT * y1,
    const SCVT * y2,
    const SCVT * y3,
    int yRot
    ) const {

    int qSize0 = lineQuadSizes[ this->quadratureOrder[ 0 ] ];
    int qSize1 = lineQuadSizes[ this->quadratureOrder[ 1 ] ];
    int qSize2 = lineQuadSizes[ this->quadratureOrder[ 2 ] ];
    int qSize3 = lineQuadSizes[ this->quadratureOrder[ 3 ] ];
    int totalSize = qSize0 * qSize1 * qSize2 * qSize3;

    const SCVT * x1rot = nullptr;
    const SCVT * x2rot = nullptr;
    const SCVT * x3rot = nullptr;
    const SCVT * y1rot = nullptr;
    const SCVT * y2rot = nullptr;
    const SCVT * y3rot = nullptr;

    switch ( yRot ) {
      case 0:
        if ( type == 1 ) {
          y1rot = y2;
          y2rot = y1;
          y3rot = y3;
        } else {
          y1rot = y1;
          y2rot = y2;
          y3rot = y3;
        }
        break;
      case 1:
        if ( type == 1 ) {
          y1rot = y3;
          y2rot = y2;
          y3rot = y1;
        } else {
          y1rot = y2;
          y2rot = y3;
          y3rot = y1;
        }
        break;
      case 2:
        if ( type == 1 ) {
          y1rot = y1;
          y2rot = y3;
          y3rot = y2;
        } else {
          y1rot = y3;
          y2rot = y1;
          y3rot = y2;
        }
        break;
    }

    switch ( xRot ) {
      case 0:
        x1rot = x1;
        x2rot = x2;
        x3rot = x3;
        break;
      case 1:
        x1rot = x2;
        x2rot = x3;
        x3rot = x1;
        break;
      case 2:
        x1rot = x3;
        x2rot = x1;
        x3rot = x2;
        break;
    }

    const SCVT * x1ref = nullptr;
    const SCVT * x2ref = nullptr;
    const SCVT * y1ref = nullptr;
    const SCVT * y2ref = nullptr;

    switch ( type ) {
      case 0:
        x1ref = this->x1RefIdentical[ simplex ];
        x2ref = this->x2RefIdentical[ simplex ];
        y1ref = this->y1RefIdentical[ simplex ];
        y2ref = this->y2RefIdentical[ simplex ];
        break;
      case 1:
        x1ref = this->x1RefCommonEdge[ simplex ];
        x2ref = this->x2RefCommonEdge[ simplex ];
        y1ref = this->y1RefCommonEdge[ simplex ];
        y2ref = this->y2RefCommonEdge[ simplex ];
        break;
      case 2:
        x1ref = this->x1RefCommonVertex[ simplex ];
        x2ref = this->x2RefCommonVertex[ simplex ];
        y1ref = this->y1RefCommonVertex[ simplex ];
        y2ref = this->y2RefCommonVertex[ simplex ];
        break;
      case 3:
        x1ref = this->x1RefDisjoint[ 0 ];
        x2ref = this->x2RefDisjoint[ 0 ];
        y1ref = this->y1RefDisjoint[ 0 ];
        y2ref = this->y2RefDisjoint[ 0 ];
        break;
    }

    int c = 0;

    SCVT * x1ss = this->x1ss;
    SCVT * x2ss = this->x2ss;
    SCVT * x3ss = this->x3ss;
    SCVT * y1ss = this->y1ss;
    SCVT * y2ss = this->y2ss;
    SCVT * y3ss = this->y3ss;

    const int align = DATA_ALIGN;
    const int width = DATA_WIDTH;

#pragma omp simd \
linear( c : 1 ) \
aligned( x1ss, x2ss, x3ss, x1ref, x2ref : align ) \
simdlen( width )
    for ( c = 0; c < totalSize; ++c ) {
      x1ss[ c ] = x1rot[ 0 ] + ( x2rot[ 0 ] - x1rot[ 0 ] ) * x1ref[ c ] +
        ( x3rot[ 0 ] - x2rot[ 0 ] ) * x2ref[ c ];

      x2ss[ c ] = x1rot[ 1 ] + ( x2rot[ 1 ] - x1rot[ 1 ] ) * x1ref[ c ] +
        ( x3rot[ 1 ] - x2rot[ 1 ] ) * x2ref[ c ];

      x3ss[ c ] = x1rot[ 2 ] + ( x2rot[ 2 ] - x1rot[ 2 ] ) * x1ref[ c ] +
        ( x3rot[ 2 ] - x2rot[ 2 ] ) * x2ref[ c ];
    }

#pragma omp simd \
linear( c : 1 ) \
aligned( y1ss, y2ss, y3ss, y1ref, y2ref : align ) \
simdlen( width )
    for ( c = 0; c < totalSize; ++c ) {
      y1ss[ c ] = y1rot[ 0 ] + ( y2rot[ 0 ] - y1rot[ 0 ] ) * y1ref[ c ] +
        ( y3rot[ 0 ] - y2rot[ 0 ] ) * y2ref[ c ];

      y2ss[ c ] = y1rot[ 1 ] + ( y2rot[ 1 ] - y1rot[ 1 ] ) * y1ref[ c ] +
        ( y3rot[ 1 ] - y2rot[ 1 ] ) * y2ref[ c ];

      y3ss[ c ] = y1rot[ 2 ] + ( y2rot[ 2 ] - y1rot[ 2 ] ) * y1ref[ c ] +
        ( y3rot[ 2 ] - y2rot[ 2 ] ) * y2ref[ c ];
    }
  }

  void getReferenceSauterSchwabData(
    int type,
    int & nSimplex,
    SCVT ** & jacobianP,
    SCVT ** & x1refP,
    SCVT ** & x2refP,
    SCVT ** & y1refP,
    SCVT ** & y2refP
    ) const {

    switch ( type ) {
      case 0:
        nSimplex = 6;
        jacobianP = this->jacobianWeightIdentical;
        x1refP = this->x1RefIdentical;
        x2refP = this->x2RefIdentical;
        y1refP = this->y1RefIdentical;
        y2refP = this->y2RefIdentical;
        break;
      case 1:
        nSimplex = 5;
        jacobianP = this->jacobianWeightCommonEdge;
        x1refP = this->x1RefCommonEdge;
        x2refP = this->x2RefCommonEdge;
        y1refP = this->y1RefCommonEdge;
        y2refP = this->y2RefCommonEdge;
        break;
      case 2:
        nSimplex = 2;
        jacobianP = this->jacobianWeightCommonVertex;
        x1refP = this->x1RefCommonVertex;
        x2refP = this->x2RefCommonVertex;
        y1refP = this->y1RefCommonVertex;
        y2refP = this->y2RefCommonVertex;
        break;
      case 3:
        nSimplex = 1;
        jacobianP = this->jacobianWeightDisjoint;
        x1refP = this->x1RefDisjoint;
        x2refP = this->x2RefDisjoint;
        y1refP = this->y1RefDisjoint;
        y2refP = this->y2RefDisjoint;
        break;
      default:
        nSimplex = 0;
        jacobianP = nullptr;
        x1refP = nullptr;
        x2refP = nullptr;
        y1refP = nullptr;
        y2refP = nullptr;
        break;
    }
  }

  void getReferenceSauterSchwabData(
    int type,
    int & nSimplex,
    SCVT ** & jacobianP
    ) const {

    switch ( type ) {
      case 0:
        nSimplex = 6;
        jacobianP = this->jacobianWeightIdentical;
        break;
      case 1:
        nSimplex = 5;
        jacobianP = this->jacobianWeightCommonEdge;
        break;
      case 2:
        nSimplex = 2;
        jacobianP = this->jacobianWeightCommonVertex;
        break;
      case 3:
        nSimplex = 1;
        jacobianP = this->jacobianWeightDisjoint;
        break;
      default:
        nSimplex = 0;
        jacobianP = nullptr;
        break;
    }
  }

  void getReferenceSauterSchwabData(
    int type,
    bool inner,
    int & nSimplex,
    SCVT ** & jacobianP,
    SCVT ** & z1refP,
    SCVT ** & z2refP
    ) const {

    switch ( type ) {
      case 0:
        nSimplex = 6;
        jacobianP = this->jacobianWeightIdentical;
        if ( !inner ) {
          z1refP = this->x1RefIdentical;
          z2refP = this->x2RefIdentical;
        } else {
          z1refP = this->y1RefIdentical;
          z2refP = this->y2RefIdentical;
        }
        break;
      case 1:
        nSimplex = 5;
        jacobianP = this->jacobianWeightCommonEdge;
        if ( !inner ) {
          z1refP = this->x1RefCommonEdge;
          z2refP = this->x2RefCommonEdge;
        } else {
          z1refP = this->y1RefCommonEdge;
          z2refP = this->y2RefCommonEdge;
        }
        break;
      case 2:
        nSimplex = 2;
        jacobianP = this->jacobianWeightCommonVertex;
        if ( !inner ) {
          z1refP = this->x1RefCommonVertex;
          z2refP = this->x2RefCommonVertex;
        } else {
          z1refP = this->y1RefCommonVertex;
          z2refP = this->y2RefCommonVertex;
        }
        break;
      case 3:
        nSimplex = 1;
        jacobianP = this->jacobianWeightDisjoint;
        if ( !inner ) {
          z1refP = this->x1RefDisjoint;
          z2refP = this->x2RefDisjoint;
        } else {
          z1refP = this->y1RefDisjoint;
          z2refP = this->y2RefDisjoint;
        }
        break;
      default:
        nSimplex = 0;
        jacobianP = nullptr;
        z1refP = nullptr;
        z2refP = nullptr;
        break;
    }
  }

  //! returns local coordinate system and triangle properties
  void getLocalCoordinates(
    const SCVT * x1,
    const SCVT * x2,
    const SCVT * x3,
    SCVT * r1,
    SCVT * r2,
    SCVT * n,
    SCVT & stau,
    SCVT & alpha1,
    SCVT & alpha2
    ) const;

  //! returns local coordinate system and triangle properties
  void getLocalCoordinates(
    const SCVT * y1,
    const SCVT * y2,
    const SCVT * y3,
    SCVT * r1,
    SCVT * r2
    ) const;

  bool checkLocalCoordinates(
    const SCVT * r1,
    const SCVT * r2,
    const SCVT * n
    ) const;

  void chooseAndUpdateLocalCoordinates(
    const SCVT * y1,
    const SCVT * y2,
    const SCVT * y3,
    const SCVT * n
    ) const;
  
  void chooseAndUpdateLocalCoordinatesFromRot2(
    const SCVT * y1,
    const SCVT * y2,
    const SCVT * y3,
    const SCVT * n
    ) const;

  void finalizeLocalCoordinates(
    const SCVT * y1,
    const SCVT * y2,
    const SCVT * y3
    ) const;

  //! returns local coordinate system and triangle properties for different order of x1, x2, x3
  void getLocalCoordinates(
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
    ) const;

  //! returns PRECOMPUTED local coordinate system and triangle properties for different rotations
  void getLocalCoordinates(
    LO i,
    SCVT * r1,
    SCVT * r2,
    SCVT * n,
    SCVT & stau,
    SCVT & alpha1,
    SCVT & alpha2,
    int order
    ) const;

  //! returns local coordinates of the quadrature node
  void getLocalQuadratureNode(
    const SCVT * x,
    const SCVT * quadratureNode,
    const SCVT * r1,
    const SCVT * r2,
    const SCVT * n,
    SCVT * xLocal
    ) const;

  //! returns local coordinates of the quadrature node for different order of x1, x2, x3
  void getLocalQuadratureNode(
    const SCVT * x1,
    const SCVT * x2,
    const SCVT * x3,
    const SCVT * quadratureNode,
    const SCVT * r1,
    const SCVT * r2,
    const SCVT * n,
    SCVT * xLocal,
    int order
    ) const;

  //! returns element matrix of single layer potential by Sauter-Schwab quadrature
  void computeElemMatrix1LayerSauterSchwabP0P0(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  //! returns element matrix of single layer potential by Sauter-Schwab quadrature
  void computeElemMatrix1LayerSauterSchwabP1P0(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  /*
   * returns local matrix for  single layer operator with p1p1 approximation
   * 
   * @param[in] outerElem index of outer element
   * @param[in] innerElem index of inner element
   * @param[out] matrix preallocated local matrix
   */
  void computeElemMatrix1LayerSauterSchwabP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  void computeElemMatrix2LayerSauterSchwabP0P0(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  //! returns element matrix of double layer potential by Sauter-Schwab quadrature
  void computeElemMatrix2LayerSauterSchwabP0P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  //! returns element matrix of double layer potential by Sauter-Schwab quadrature
  void computeElemMatrix2LayerSauterSchwabP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;
  //public:
  //! determines the type of Sauter Schwab quadrature
  void getSauterSchwabQuadratureType(
    LO inner,
    LO outer,
    int &type,
    int &rotInner,
    int &rotOuter
    ) const;
  //protected:
  //! returns element matrix of single layer potential for regular pairs
  void computeElemMatrix1LayerDisjointP0P0(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  /*
   * returns local matrix for single layer operator with p1p1 approximation
   * 
   * @param[in] outerElem index of outer element
   * @param[in] innerElem index of inner element
   * @param[out] matrix preallocated local matrix
   */
  void computeElemMatrix1LayerDisjointP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  //! returns element matrix of double layer potential for regular pairs
  void computeElemMatrix2LayerDisjointP0P0(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  //! returns element matrix of single layer potential for regular pairs
  void computeElemMatrix2LayerDisjointP0P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  void computeElemMatrix2LayerDisjointP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  //! returns whether two panels are disjoint
  bool areElementsDisjoint(
    LO outerElem,
    LO innerElem,
    SCVT threshold = 2.0
    ) const;

  /*! maps coordinates from 4D unit cube to the reference triangle
   *  for two identical faces
   *
   * @param[in]       w = [ksi, eta1, eta2, et3] pointer to array of coordinates on the 4D unit cube
   * @param[in]       subdomain number of 4D unit cube subdomain
   * @param[in,out]   x coordinates on the reference triangle
   * @param[in,out]   y coordinates on the reference triangle
   * @param[in,out]   jacobian  jacobian of the transformation
   */
  void cube2triIdentical(
    const SCVT * w,
    int subdomain,
    SCVT * x,
    SCVT * y,
    SCVT & jacobian
    ) const;

  /*! maps coordinates from 4D unit cube to the reference triangle
   *  for two faces with common edge
   *
   * @param[in]       w = [ksi, eta1, eta2, et3] pointer to array of coordinates on the 4D unit cube
   * @param[in]       subdomain number of 4D unit cube subdomain
   * @param[in,out]   x coordinates on the reference triangle
   * @param[in,out]   y coordinates on the reference triangle
   * @param[in,out]   jacobian  jacobian of the transformation
   */
  void cube2triCommonEdge(
    const SCVT * w,
    int subdomain,
    SCVT * x,
    SCVT * y,
    SCVT & jacobian
    ) const;

  /*! maps coordinates from 4D unit cube to the reference triangle
   *  for two faces with common vertex
   *
   * @param[in]       w = [ksi, eta1, eta2, et3] pointer to array of coordinates on the 4D unit cube
   * @param[in]       subdomain number of 4D unit cube subdomain
   * @param[in,out]   x coordinates in the reference triangle
   * @param[in,out]   y coordinates in the reference triangle
   * @param[in,out]   jacobian  jacobian of the transformation
   */
  void cube2triCommonVertex(
    const SCVT * w,
    int subdomain,
    SCVT * x,
    SCVT * y,
    SCVT & jacobian
    ) const;

  /*! maps coordinates from 4D unit cube to the reference triangle
   *  for two disjoint faces
   *
   * @param[in]       w = [ksi, eta1, eta2, et3] pointer to array of coordinates on the 4D unit cube
   * @param[in]       subdomain number of 4D unit cube subdomain
   * @param[in,out]   x coordinates in the reference triangle
   * @param[in,out]   y coordinates in the reference triangle
   * @param[in,out]   jacobian  jacobian of the transformation
   */
  void cube2triDisjoint(
    const SCVT * w,
    SCVT * x,
    SCVT * y,
    SCVT & jacobian
    ) const;

  /*! maps from reference triangle (0,0), (1,0), (0,1) to element triangle given by its vertices
   *
   * @param[in]       x1  vertex of triangular element
   * @param[in]       x2  vertex of triangular element
   * @param[in]       x3  vertex of triangular element
   * @param[in]       x   coordinates in the reference element
   * @param[in,out]   globalCoord  coordinates in a global 3D space
   */
  void tri2element(
    const SCVT * x1,
    const SCVT * x2,
    const SCVT * x3,
    const SCVT * x,
    SCVT * globalCoord
    ) const;

  /*! maps from reference triangle (0,0), (1,0), (0,1) to element triangle given by its vertices
   *
   * @param[in]       x1  vertex of triangular element
   * @param[in]       x2  vertex of triangular element
   * @param[in]       x3  vertex of triangular element
   * @param[in]       x   coordinates in the reference element
   * @param[in]       rot rotates the output triangle
   * @param[in,out]   globalCoord  coordinates in a global 3D space
   * @param[in]       swap swaps first two nodes
   */
  void tri2element(
    const SCVT * x1,
    const SCVT * x2,
    const SCVT * x3,
    const SCVT * x,
    int rot,
    SCVT * globalCoord,
    bool swap = false
    ) const;

  //! returns specific kernel evaluated in given points (x, y)

  SC evalSingleLayerKernel(
    const SCVT * x,
    const SCVT * y
    ) const {
    return static_cast<const SpecificIntegrator*> ( this )->
      evalSingleLayerKernel( x, y );
  };

  SC evalSingleLayerKernel(
    SCVT x1,
    SCVT x2,
    SCVT x3,
    SCVT y1,
    SCVT y2,
    SCVT y3
    ) const {
    return static_cast<const SpecificIntegrator*> ( this )->
      evalSingleLayerKernel( x1, x2, x3, y1, y2, y3 );
  };

  /*! returns specific kernel evaluated in given points (x, y)
   *
   * @param[in]       x 
   * @param[in]       y
   * @param[in,out]   n unit outer normal to the given triangle
   */
  SC evalDoubleLayerKernel(
    const SCVT * x,
    const SCVT * y,
    const SCVT * n
    ) const {
    return static_cast<const SpecificIntegrator*> ( this )->
      evalDoubleLayerKernel( x, y, n );
  };

  SC evalDoubleLayerKernel(
    SCVT x1,
    SCVT x2,
    SCVT x3,
    SCVT y1,
    SCVT y2,
    SCVT y3,
    SCVT n1,
    SCVT n2,
    SCVT n3
    ) const {
    return static_cast<const SpecificIntegrator*> ( this )->
      evalDoubleLayerKernel( x1, x2, x3, y1, y2, y3, n1, n2, n3 );
  };

};

}

// include .cpp file to overcome linking problems due to templates
#include "BEIntegrator.cpp"

#endif /* BEINTEGRATOR_H */
