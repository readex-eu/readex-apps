/*!
 * @file    BEIntegratorLaplace.h
 * @author  Michal Merta
 * @author  Jan Zapletal
 * @date    July 10, 2013
 * @brief   Header file for class BEIntegratorLaplace
 *
 */

#ifndef BEINTEGRATORLAPLACE_H
#define BEINTEGRATORLAPLACE_H

#include "BEIntegrator.h"
//#include "BESpace.h"


namespace bem4i {

/*!
 * concrete class for Laplace integrators over triangular elements
 *
 * the class uses Curiously recurring template pattern to replace virtual methods
 */
template<class LO, class SC>
class BEIntegratorLaplace : public BEIntegrator<LO, SC, BEIntegratorLaplace<LO, SC> > {
  // to get inner type of complex numbers (for Helmholtz)
  typedef typename GetType<LO, SC>::SCVT SCVT;

  // we have to enable BEIntegrator to use kernel evaluation private methods
  friend class BEIntegrator<LO, SC, BEIntegratorLaplace<LO, SC> >;

public:

  //! default constructor
  BEIntegratorLaplace( );

  //! copy constructor
  BEIntegratorLaplace(
    const BEIntegratorLaplace& orig
    );

  //! constructor taking BESpace as the argument
  BEIntegratorLaplace(
    BESpace<LO, SC>* space,
    int* quadratureOrder,
    quadratureType quadrature = SauterSchwab,
    int* quadratureOrderDisjointElems = nullptr
    );

  //! destructor
  virtual ~BEIntegratorLaplace( );

  //! returns element matrix of single layer potential
  void computeElemMatrix1Layer(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  //! returns element matrix of double layer potential
  void computeElemMatrix2Layer(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  //! returns element matrix of hypersingular operator
  void computeElemMatrixHypersingular(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  //! returns element matrix of hypersingular operator
  void computeElemMatrixHypersingular(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& V,
    FullMatrix<LO, SC>& matrix
    ) const;

  /*!
   * evaluates the representation formula in points x, stores values in
   * preallocated vector values
   *
   * @param[in] x pointer to array with evaluation points
   * @param[in] n number of evaluation points
   * @param[in] dir Dirichlet data
   * @param[in] neu Neumann data
   * @param[in] interior flag interior/exterior
   * @param[out] values preallocated vector for storing results
   */
  void representationFormula(
    const SC * x,
    LO nPoints,
    const Vector<LO, SC> & dir,
    const Vector<LO, SC> & neu,
    bool interior,
    Vector<LO, SC> & values
    ) const {

    if ( this->space->getTestFunctionType( ) == p1 &&
      this->space->getAnsatzFunctionType( ) == p1 ) {
      this->representationFormulaP1P1( x, nPoints, dir, neu, interior, values );
    } else if ( this->space->getTestFunctionType( ) == p0 &&
      this->space->getAnsatzFunctionType( ) == p0 ) {
      std::cout << "Not implemented!" << std::endl;
      //this->representationFormulaP0P0( x, nPoints, dir, neu, interior, values );
    } else {
      this->representationFormulaP1P0( x, nPoints, dir, neu, interior, values );
    }

  }

  /*!
   * evaluates the representation formula in points x
   * for p1 Dirchlet, p0 Neumann data,
   * stores values in preallocated vector values
   *
   * @param[in] x pointer to array with evaluation points
   * @param[in] n number of evaluation points
   * @param[in] dir Dirichlet data
   * @param[in] neu Neumann data
   * @param[in] interior flag interior/exterior
   * @param[out] values preallocated vector for storing results
   */
  void representationFormulaP1P0(
    const SCVT * x,
    LO n,
    const Vector<LO, SC> & dir,
    const Vector<LO, SC> & neu,
    bool interior,
    Vector<LO, SC> & values
    ) const;

  void representationFormulaP1P1(
    const SCVT * x,
    LO n,
    const Vector<LO, SC> & dir,
    const Vector<LO, SC> & neu,
    bool interior,
    Vector<LO, SC> & values
    ) const;

  /*
   * evaluates double layer potential in points x, stores values in
   * preallocated vector values
   *
   * @param[in] x pointer to array with evaluation points
   * @param[in] nPoints number of evaluation points
   * @param[in] density density function
   * @param[out] values preallocated vector for storing results
   */
  void doubleLayerPotential(
    const SCVT * x,
    LO nPoints,
    const Vector<LO, SC> & density,
    Vector<LO, SC> & values
    ) const {
    doubleLayerPotentialP1( x, nPoints, density, values );
  }

  /*
   * evaluates double layer potential in points x
   * for p1 density,
   * stores values in preallocated vector values
   *
   * @param[in] x pointer to array with evaluation points
   * @param[in] nPoints number of evaluation points
   * @param[in] density density function
   * @param[out] values preallocated vector for storing results
   */
  void doubleLayerPotentialP1(
    const SCVT * x,
    LO nPoints,
    const Vector<LO, SC> & density,
    Vector<LO, SC> & values
    ) const;


  //#if N_MIC > 0
  //! computes element matrix on MIC
#if N_MIC > 0
  __attribute__ ( ( target( mic ) ) )
#endif
  static void computeElemMatrix1LayerP0P0MIC(
    const SCVT * nodes,
    const LO * elems,
    const SCVT * areas,
    LO outerElem,
    LO innerElem,
    int qOrdOuter,
    int qOrdInner,
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
    );

#if N_MIC > 0
  __attribute__ ( ( target( mic ) ) )
#endif
  static void computeElemMatrix2LayerP0P1MIC(
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
    );

#if N_MIC > 0
  __attribute__ ( ( target( mic ) ) )
#endif
  static void computeElemMatrix1And2LayerP0P0P0P1MIC(
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
    );

  //#endif

private:

  //! returns element matrix of single layer potential with p0p0 approximation
  void computeElemMatrix1LayerP0P0(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  //! returns element matrix of single layer potential with p0p0 approximation
  void computeElemMatrix1LayerP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  //! returns element matrix of double layer potential with p0p1 approximation
  void computeElemMatrix2LayerP0P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  //! returns element matrix of double layer potential with p0p0 approximation
  void computeElemMatrix2LayerP0P0(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC> & matrix
    ) const;

  //! returns element matrix of double layer potential with p0p0 approximation
  void computeElemMatrix2LayerP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC> & matrix
    ) const;

  //! returns element matrix of hypersingular operator with p1p1 approximation
  void computeElemMatrixHypersingularP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  //! returns specific kernel evaluated in given points (x, y)

  SC evalSingleLayerKernel(
    const SCVT *x,
    const SCVT *y
    ) const {
    SCVT norm = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) +
      ( x[1] - y[1] )*( x[1] - y[1] ) +
      ( x[2] - y[2] )*( x[2] - y[2] ) );
    return ( PI_FACT / norm );
  };

#pragma omp declare simd simdlen( DATA_WIDTH )

  SC evalSingleLayerKernel(
    SCVT x1,
    SCVT x2,
    SCVT x3,
    SCVT y1,
    SCVT y2,
    SCVT y3
    ) const {
    SCVT norm = std::sqrt( ( x1 - y1 ) * ( x1 - y1 ) + ( x2 - y2 ) * ( x2 - y2 ) +
      ( x3 - y3 ) * ( x3 - y3 ) );
    return ( PI_FACT / norm );
  };

  //! returns specific kernel evaluated in given points (x, y)

#pragma omp declare simd uniform( n1, n2, n3 ) simdlen( DATA_WIDTH )

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

    SCVT diff1 = x1 - y1;
    SCVT diff2 = x2 - y2;
    SCVT diff3 = x3 - y3;
    SCVT norm = std::sqrt( diff1 * diff1 + diff2 * diff2 + diff3 * diff3 );
    SCVT dot = diff1 * n1 + diff2 * n2 + diff3 * n3;

    return ( PI_FACT * dot / ( norm * norm * norm ) );
  };

  SC evalDoubleLayerKernel( const SCVT *x, const SCVT *y, const SCVT* n ) const {
    SCVT norm = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) +
      ( x[1] - y[1] )*( x[1] - y[1] ) +
      ( x[2] - y[2] )*( x[2] - y[2] ) );
    SCVT dot = ( x[0] - y[0] ) * n[0] +
      ( x[1] - y[1] ) * n[1] +
      ( x[2] - y[2] ) * n[2];
    return ( PI_FACT * dot / ( norm * norm * norm ) );
  };

#pragma omp declare simd uniform( stau, alpha1, alpha2 ) simdlen( DATA_WIDTH )
  SC collocation1LayerP0(
    SC sx,
    SC tx,
    SC ux,
    SC stau,
    SC alpha1,
    SC alpha2
    ) const;

  //! evaluates p0 Laplace single layer operator in xLocal
#pragma omp declare simd uniform( stau, alpha1, alpha2 ) simdlen( DATA_WIDTH )
  SC collocation1LayerP1(
    SC sx,
    SC tx,
    SC ux,
    SC stau,
    SC alpha1,
    SC alpha2
    ) const;

#pragma omp declare simd uniform( s, alpha ) simdlen( DATA_WIDTH )
  SC f1LayerP0(
    SC sx,
    SC tx,
    SC ux,
    SC s,
    SC alpha
    ) const;

  //! help function for collocation of single layer operator
#pragma omp declare simd uniform( stau, alpha ) simdlen( DATA_WIDTH )
  SC f1LayerP1(
    SC stau,
    SC alpha,
    SC tx,
    SC sx,
    SC ux
    ) const;

#pragma omp declare simd uniform( stau, alpha1, alpha2 ) simdlen( DATA_WIDTH )
  SC collocation2LayerP1(
    SC sx,
    SC tx,
    SC ux,
    SC stau,
    SC alpha1,
    SC alpha2
    ) const;

#pragma omp declare simd simdlen( DATA_WIDTH )
  void collocation2LayerP1All(
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
    ) const;

  //! help function for collocation of double layer operator
#pragma omp declare simd uniform( alpha, stau ) simdlen( DATA_WIDTH )
  SC f2LayerP1(
    SC alpha,
    SC tx,
    SC sx,
    SC ux,
    SC stau
    ) const;

#pragma omp declare simd simdlen( DATA_WIDTH )
  void f2LayerP1All(
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
    ) const;

  //! evaluates p0 Laplace double layer operator in xLocal
#pragma omp declare simd uniform( stau, alpha1, alpha2 ) simdlen( DATA_WIDTH )
  SC collocation2LayerP0(
    SC sx,
    SC tx,
    SC ux,
    SC stau,
    SC alpha1,
    SC alpha2
    ) const;

  //! help function for collocation of double layer operator
#pragma omp declare simd uniform( alpha, stau ) simdlen( DATA_WIDTH )
  SC f2LayerP0(
    SC alpha,
    SC tx,
    SC sx,
    SC ux,
    SC stau
    ) const;

};

}
// include .cpp file to overcome linking problems due to templates
#include "BEIntegratorLaplace.cpp"

#endif /* BEINTEGRATORLAPLACE_H */
