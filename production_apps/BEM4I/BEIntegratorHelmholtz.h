/*!
 * @file    BEIntegratorHelmholtz.h
 * @author  Jan Zapletal
 * @date    August 12, 2013
 * @brief   Header file for class BEIntegratorHelmholtz
 * 
 */

#ifndef BEINTEGRATORHELMHOLTZ_H
#define BEINTEGRATORHELMHOLTZ_H

#include "BEIntegrator.h"
//#include "BESpace.h"

namespace bem4i {

/*!
 * concrete class for Helmholtz integrators over triangular elements
 * 
 * the class uses Curiously recurring template pattern to replace virtual methods
 */
template<class LO, class SC>
class BEIntegratorHelmholtz : public BEIntegrator<LO, SC, BEIntegratorHelmholtz<LO, SC> > {
  typedef typename GetType<LO, SC>::SCVT SCVT;

  friend class BEIntegrator<LO, SC, BEIntegratorHelmholtz<LO, SC> >;

public:

  //! default constructor
  BEIntegratorHelmholtz( );

  //! copy constructor
  BEIntegratorHelmholtz( const BEIntegratorHelmholtz& orig );

  //! constructor taking BESpace as the argument
  BEIntegratorHelmholtz(
    BESpace<LO, SC>* space,
    int* quadratureOrder,
    SC kappa,
    quadratureType quadrature = SauterSchwab,
    int* quadratureOrderDisjointElems = nullptr
    );

  //! destructor
  virtual ~BEIntegratorHelmholtz( );

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

  void computeElemMatrixH1SauterSchwabP0P0(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  void computeElemMatrixH2SauterSchwabP0P0(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  /*
   * evaluates the representation formula in points x, stores values in
   * preallocated vector values
   * 
   * @param[in] xCoord pointer to array of evaluation points
   * @param[in] n number of evaluation points
   * @param[in] dir Dirichlet data
   * @param[in] neu Neumann data
   * @param[in] interior flag interior/exterior
   * @param[out] values preallocated vector for storing results
   */
  void representationFormula(
    const SCVT *xCoord,
    LO n,
    const Vector<LO, SC> & dir,
    const Vector<LO, SC> & neu,
    bool interior,
    Vector<LO, SC> & values
    ) const {

    if ( this->space->getTestFunctionType( ) == p1 &&
      this->space->getAnsatzFunctionType( ) == p1 ) {
      representationFormulaP1P1( xCoord, n, dir, neu, interior, values );
    } else if ( this->space->getTestFunctionType( ) == p0 &&
      this->space->getAnsatzFunctionType( ) == p0 ) {
      representationFormulaP0P0( xCoord, n, dir, neu, interior, values );
    } else {
      representationFormulaP1P0( xCoord, n, dir, neu, interior, values );
    }
  };

  /*
   * evaluates the representation formula in points x 
   * for p1 Dirchlet, p0 Neumann data,
   * stores values in preallocated vector values
   * 
   * @param[in] xCoord pointer to array of evaluation points
   * @param[in] n number of evaluation points
   * @param[in] dir Dirichlet data
   * @param[in] neu Neumann data
   * @param[in] interior flag interior/exterior
   * @param[out] values preallocated vector for storing results
   */
  void representationFormulaP1P0(
    const SCVT *x,
    LO n,
    const Vector<LO, SC> & dir,
    const Vector<LO, SC> & neu,
    bool interior,
    Vector<LO, SC> & values
    ) const;

  /*
   * evaluates the representation formula in points x 
   * for p1 Dirchlet, p1 Neumann data,
   * stores values in preallocated vector values
   * 
   * @param[in] xCoord pointer to array of evaluation points
   * @param[in] n number of evaluation points
   * @param[in] dir Dirichlet data
   * @param[in] neu Neumann data
   * @param[in] interior flag interior/exterior
   * @param[out] values preallocated vector for storing results
   */
  void representationFormulaP1P1(
    const SCVT *x,
    LO n,
    const Vector<LO, SC> & dir,
    const Vector<LO, SC> & neu,
    bool interior,
    Vector<LO, SC> & values
    ) const;

  /*
   * evaluates the representation formula in points x 
   * for p0 Dirchlet, p0 Neumann data,
   * stores values in preallocated vector values
   * 
   * @param[in] xCoord pointer to array of evaluation points
   * @param[in] n number of evaluation points
   * @param[in] dir Dirichlet data
   * @param[in] neu Neumann data
   * @param[in] interior flag interior/exterior
   * @param[out] values preallocated vector for storing results
   */
  void representationFormulaP0P0(
    const SCVT *x,
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
    this->doubleLayerPotentialP1( x, nPoints, density, values );
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

  /*
   * evaluates single layer potential in points x, stores values in
   * preallocated vector values
   * 
   * @param[in] x pointer to array with evaluation points
   * @param[in] nPoints number of evaluation points
   * @param[in] density density function
   * @param[out] values preallocated vector for storing results
   */
  void singleLayerPotential(
    const SCVT * x,
    LO nPoints,
    const Vector<LO, SC> & density,
    Vector<LO, SC> & values
    ) const {
    this->singleLayerPotentialP0( x, nPoints, density, values );
  }

  /*
   * evaluates single layer potential in points x 
   * for p1 density,
   * stores values in preallocated vector values
   * 
   * @param[in] x pointer to array with evaluation points
   * @param[in] nPoints number of evaluation points
   * @param[in] density density function
   * @param[out] values preallocated vector for storing results
   */
  void singleLayerPotentialP0(
    const SCVT * x,
    LO nPoints,
    const Vector<LO, SC> & density,
    Vector<LO, SC> & values
    ) const;

protected:

  //! returns element matrix of single layer potential for regular pairs
  void computeElemMatrixHypersingularDisjointP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

private:

  /*
   * returns local matrix for Helmholtz single layer operator with p0p0 approximation
   * 
   * @param[in] outerElem index of outer element
   * @param[in] innerElem index of inner element
   * @param[out] matrix preallocated local matrix
   */
  void computeElemMatrix1LayerP0P0(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  //! returns element matrix of single layer potential with p1p1 approximation
  void computeElemMatrix1LayerP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC> & matrix
    ) const;


  void computeElemMatrix1LayerSauterSchwabP0P0(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  void computeElemMatrix1LayerSauterSchwabP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  void computeElemMatrix1LayerDisjointP0P0(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  void computeElemMatrix1LayerDisjointP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  /*
   * returns local matrix for Helmholtz double layer operator with p0p1 approximation
   * 
   * @param[in] outerElem index of outer element
   * @param[in] innerElem index of inner element
   * @param[out] matrix preallocated local matrix
   */
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

  //! returns element matrix of double layer potential with p1p1 approximation
  void computeElemMatrix2LayerP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC> & matrix
    ) const;

  void computeElemMatrix2LayerSauterSchwabP0P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  void computeElemMatrix2LayerSauterSchwabP0P0(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  void computeElemMatrix2LayerSauterSchwabP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  void computeElemMatrix2LayerDisjointP0P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  void computeElemMatrix2LayerDisjointP0P0(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  void computeElemMatrix2LayerDisjointP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;


  /*
   * returns local matrix for Helmholtz hypersingular operator with p0p0 approximation
   * 
   * @param[in] outerElem index of outer element
   * @param[in] innerElem index of inner element
   * @param[in] precomputed surface curls of test function
   * @param[in] precomputed surface curls of ansatz functions
   * @param[out] matrix preallocated local matrix
   */
  void computeElemMatrixHypersingularP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

  /*
   * returns local matrix for Helmholtz hypersingular operator with p0p0 approximation
   * 
   * @param[in] outerElem index of outer element
   * @param[in] innerElem index of inner element
   * @param[out] matrix preallocated local matrix
   */
  void computeElemMatrixHypersingularSauterSchwabP1P1(
    LO outerElem,
    LO innerElem,
    FullMatrix<LO, SC>& matrix
    ) const;

private:

  /*! returns specific kernel evaluated in given points (x, y)
   *
   * @param[in]       x 
   * @param[in]       y
   */
  SC evalSingleLayerKernel( const SCVT *x, const SCVT *y ) const {
    SCVT norm = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) +
      ( x[1] - y[1] )*( x[1] - y[1] ) +
      ( x[2] - y[2] )*( x[2] - y[2] ) );
    SC i( 0.0, 1.0 );
    return ( PI_FACT * std::exp( i * kappa * norm ) / norm );
  };

  SC evalSingleLayerKernel(
    SCVT x1,
    SCVT x2,
    SCVT x3,
    SCVT y1,
    SCVT y2,
    SCVT y3
    ) const {

    SCVT real, imag;

    SCVT rekappa = std::real( kappa );
    SCVT imkappa = std::imag( kappa );

    SCVT norm = std::sqrt(
      ( x1 - y1 ) * ( x1 - y1 )
      + ( x2 - y2 ) * ( x2 - y2 )
      + ( x3 - y3 ) * ( x3 - y3 ) );
    SCVT rekappan = rekappa * norm;

    SCVT expim = std::exp( -imkappa * norm );

    real = (SCVT) PI_FACT * expim * std::cos( rekappan ) / norm;
    imag = (SCVT) PI_FACT * expim * std::sin( rekappan ) / norm;

    return SC( real, imag );
  };

#pragma omp declare simd simdlen( DATA_WIDTH )

  void evalSingleLayerKernel(
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

    SCVT norm = std::sqrt(
      ( x1 - y1 ) * ( x1 - y1 )
      + ( x2 - y2 ) * ( x2 - y2 )
      + ( x3 - y3 ) * ( x3 - y3 ) );
    SCVT rekappan = rekappa * norm;

    SCVT expim = std::exp( -imkappa * norm );

    real = (SCVT) PI_FACT * expim * std::cos( rekappan ) / norm;
    imag = (SCVT) PI_FACT * expim * std::sin( rekappan ) / norm;
  };

  /*! returns specific kernel evaluated in given points (x, y)
   *
   * @param[in]       x 
   * @param[in]       y
   * @param[in,out]   n unit outer normal to the given triangle
   */
  SC evalDoubleLayerKernel( const SCVT *x, const SCVT *y, const SCVT *n ) const {
    SCVT norm = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) + ( x[1] - y[1] )*
      ( x[1] - y[1] ) + ( x[2] - y[2] )*( x[2] - y[2] ) );
    SCVT dot = ( x[0] - y[0] ) * n[0] + ( x[1] - y[1] ) * n[1] + ( x[2] - y[2] )
      * n[2];
    SC i( 0.0, 1.0 );
    return PI_FACT * ( dot / ( norm * norm * norm ) ) *
      ( (SCVT) 1.0 - i * kappa * norm ) * std::exp( i * kappa * norm );
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

    SCVT real, imag;

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
      * ( ( (SCVT) 1.0 + imkappan ) * cosine + rekappan * sine ) );
    imag = mult * expim
      * ( ( (SCVT) 1.0 + imkappan ) * sine - rekappan * cosine );

    return SC( real, imag );
  };

#pragma omp declare simd uniform( n1, n2, n3 ) simdlen( DATA_WIDTH )

  void evalDoubleLayerKernel(
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
      * ( ( (SCVT) 1.0 + imkappan ) * cosine + rekappan * sine ) );
    imag = mult * expim
      * ( ( (SCVT) 1.0 + imkappan ) * sine - rekappan * cosine );
  };

  SC evalHypersingularP0P0Kernel(
    const SCVT *x,
    const SCVT *y,
    const SCVT *nx,
    const SCVT *ny ) const {
    SC i( 0.0, 1.0 );
    SCVT norm = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) + ( x[1] - y[1] )*
      ( x[1] - y[1] ) + ( x[2] - y[2] )*( x[2] - y[2] ) );
    SCVT dotnx = ( x[0] - y[0] ) * nx[0] + ( x[1] - y[1] ) * nx[1] +
      ( x[2] - y[2] ) * nx[2];
    SCVT dotny = ( x[0] - y[0] ) * ny[0] + ( x[1] - y[1] ) * ny[1] +
      ( x[2] - y[2] ) * ny[2];

    SC e = std::exp( i * kappa * norm );
    SC norm3 = norm * norm*norm;
    SC tmp = 0.0;
    for ( int j = 0; j < 3; j++ ) {
      tmp += ( nx[j] * ny[j]*( e * ( i * kappa * norm - (SCVT) 1.0 ) ) / ( norm3 ) );
    }
    return (SCVT) PI_FACT * ( ( e * dotnx * dotny * ( kappa * kappa * norm *
      norm + (SCVT) 3.0 *
      ( i * kappa * norm - (SCVT) 1.0 ) ) ) / ( norm3 * norm * norm ) - tmp );
  };

  //! wave number
  SC kappa;



  /*
   * evaluates Helmholtz p0 single layer operator in xLocal
   * 
   * @param[in] x outer quadrature point in cartesian coordinates
   * @param[in] xLocal outer quadrature point in local coordinates
   * @param[in] stau inner triangle parameter
   * @param[in] alpha1 inner triangle parameter
   * @param[in] alpha2 inner triangle parameter
   * @param[in] quadratureOrder quadrature order for inner quadrature
   * @param[in] quadratureNodes inner quadrature nodes
   * 
   * @return value in xLocal
   * 
   * @todo send shifted kernel values directly to be consistent with double layer
   *       and hypersingular operators
   */
  SC collocation1LayerP0(
    const SCVT* xLocal,
    SCVT stau,
    SCVT alpha1,
    SCVT alpha2,
    int quadratureOrder,
    const SC* shiftedKernel,
    LO elem ) const;

#pragma omp declare simd uniform( stau, alpha1, alpha2 ) simdlen( DATA_WIDTH )
  SCVT collocationSingular1LayerP0(
    SCVT sx,
    SCVT tx,
    SCVT ux,
    SCVT stau,
    SCVT alpha1,
    SCVT alpha2
    ) const;

  /*
   * evaluates p1 Helmholtz single layer operator in xLocal
   * (used for p1 hypersingular operator)
   * 
   * @param[in] xLocal outer quadrature point in local coordinates
   * @param[in] stau inner triangle parameter
   * @param[in] alpha1 inner triangle parameter
   * @param[in] alpha2 inner triangle parameter
   * @param[in] quadratureOrder quadrature order for inner quadrature
   * @param[in] shiftedKernel shifted kernel evaluated in x and inner quadrature points
   * @param[in] elem index of inner element
   * @param[in] rot rotation index of inner element
   * 
   * @return value in xLocal
   */
  SC collocation1LayerP1(
    const SCVT* xLocal,
    SCVT stau,
    SCVT alpha1,
    SCVT alpha2,
    int quadratureOrder,
    const SC* shiftedKernel,
    LO elem,
    int rot
    ) const;

#pragma omp declare simd uniform( stau, alpha1, alpha2 ) simdlen( DATA_WIDTH )
  SCVT collocationSingular1LayerP1(
    SCVT sx,
    SCVT tx,
    SCVT ux,
    SCVT stau,
    SCVT alpha1,
    SCVT alpha2
    ) const;

  /*
   * auxiliary function for p0 single layer operator
   * 
   * @param[in] xLocal outer quadrature point in local coordinates
   * @param[in] s inner triangle parameter
   * @param[in] alpha inner triangle parameter
   * 
   * @return value in xLocal
   */
#pragma omp declare simd uniform( s, alpha ) simdlen( DATA_WIDTH )
  SCVT f1LayerP0(
    SCVT sx,
    SCVT tx,
    SCVT ux,
    SCVT s,
    SCVT alpha
    ) const;

  /*
   * auxiliary function for p1 single layer operator
   * (used for the hypersingular operator)
   * 
   * @param[in] s inner triangle parameter
   * @param[in] stau inner triangle parameter
   * @param[in] alpha inner triangle parameter
   * @param[in] tx xLocal coordinate
   * @param[in] sx xLocal coordinate
   * @param[in] ux xLocal coordinate
   * 
   * @return value in xLocal
   */
#pragma omp declare simd uniform( stau, alpha ) simdlen( DATA_WIDTH )
  SCVT f1LayerP1(
    SCVT stau,
    SCVT alpha,
    SCVT tx,
    SCVT sx,
    SCVT ux
    ) const;

  /*
   * evaluates p1 Helmholtz double layer operator in xLocal
   * 
   * @param[in] xLocal outer quadrature point in local coordinates
   * @param[in] stau inner triangle parameter
   * @param[in] alpha1 inner triangle parameter
   * @param[in] alpha2 inner triangle parameter
   * @param[in] quadratureOrder quadrature order for inner quadrature
   * @param[in] shiftedKernelNeumann normal derivative of the shifted kernel 
   *            evaluated in x and inner quadrature points
   * @param[in] elem index of inner element
   * @param[in] rot rotation index of inner element
   * 
   * @return value in xLocal
   */
  SC collocation2LayerP1(
    const SCVT* xLocal,
    SCVT stau,
    SCVT alpha1,
    SCVT alpha2,
    int quadratureOrder,
    const SC* shiftedKernelNeumann,
    LO elem,
    int rot
    ) const;

#pragma omp declare simd uniform( stau, alpha1, alpha2 ) simdlen( DATA_WIDTH )
  SCVT collocationSingular2LayerP1(
    SCVT sx,
    SCVT tx,
    SCVT ux,
    SCVT stau,
    SCVT alpha1,
    SCVT alpha2
    ) const;

#pragma omp declare simd uniform( stau, alpha1, alpha2 ) simdlen( DATA_WIDTH )
  SCVT collocationSingular2LayerP0(
    SCVT sx,
    SCVT tx,
    SCVT ux,
    SCVT stau,
    SCVT alpha1,
    SCVT alpha2
    ) const;

  /*
   * auxiliary function for p1 double layer operator
   * 
   * @param[in] alpha inner triangle parameter
   * @param[in] tx xLocal coordinate
   * @param[in] sx xLocal coordinate
   * @param[in] ux xLocal coordinate
   * @param[in] stau inner triangle parameter
   * 
   * @return value in xLocal
   */
#pragma omp declare simd uniform( alpha, stau ) simdlen( DATA_WIDTH )
  SCVT f2LayerP1(
    SCVT alpha,
    SCVT tx,
    SCVT sx,
    SCVT ux,
    SCVT stau
    ) const;

  //! help function for collocation of double layer operator
#pragma omp declare simd uniform( alpha, stau ) simdlen( DATA_WIDTH )
  SCVT f2LayerP0(
    SCVT alpha,
    SCVT tx,
    SCVT sx,
    SCVT ux,
    SCVT stau
    ) const;

  /*
   * evaluates the shifted Helmholtz kernel in x, y
   * (incl. the limit case x=y)
   * 
   * @param[in] x outer point
   * @param[in] y inner point
   * 
   * @return value in xLocal
   */
  SC evalShiftedKernel(
    const SCVT* x,
    const SCVT* y
    ) const;

#pragma omp declare simd simdlen( DATA_WIDTH )
  void evalShiftedKernel(
    SCVT x1,
    SCVT x2,
    SCVT x3,
    SCVT y1,
    SCVT y2,
    SCVT y3,
    SCVT & real,
    SCVT & imag
    ) const;

  /*
   * evaluates normal derivative of the shifted Helmholtz kernel in x, y
   * 
   * @param[in] x outer point
   * @param[in] y inner point
   * @param[in] n inner normal
   * 
   * @return value in xLocal
   */
  SC evalShiftedKernelNeumann(
    const SCVT* x,
    const SCVT* y,
    const SCVT* n
    ) const;

#pragma omp declare simd uniform( n1, n2, n3 ) simdlen( DATA_WIDTH )
  void evalShiftedKernelNeumann(
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
    ) const;

};

}
// include .cpp file to overcome linking problems due to templates
#include "BEIntegratorHelmholtz.cpp"

#endif /* BEINTEGRATORHELMHOLTZ_H */

