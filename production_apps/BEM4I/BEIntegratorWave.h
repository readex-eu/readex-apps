/*!
 * @file    BEIntegratorWave.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    December 12, 2013
 * @brief   Header file for class BEIntegratorWave
 * 
 */

#ifndef BEINTEGRATORWAVE_H
#define BEINTEGRATORWAVE_H

#include <algorithm>

#include "BEIntegrator.h"

#ifndef EXPERIMENTAL_WAVE

// use ordinary basis functions described in Sauter, Veit: Retarded time-domain...

namespace bem4i {

/*!
 * class for integrators for Sauter-Veit time domain BEM for wave equation
 * 
 */
template<class LO, class SC>
class BEIntegratorWave : public BEIntegrator<LO, SC, BEIntegratorWave<LO, SC> > {
  // to get inner type of complex numbers (for Helmholtz)
  typedef typename GetType<LO, SC>::SCVT SCVT;

  // we have to enable BEIntegrator to use kernel evaluation private methods
  friend class BEIntegrator<LO, SC, BEIntegratorWave<LO, SC> >;

public:

  //! default constructor
  BEIntegratorWave( );

  //! copy constructor
  BEIntegratorWave( const BEIntegratorWave& orig );

  //! constructor taking BESpace as the argument
  BEIntegratorWave(
      BESpace<LO, SC>* space,
      int* quadratureOrder,
      int timeQuadOrder,
      int* quadratureOrderDisjointElems = nullptr,
      int nChebIntervals = 20,
      int nChebPoints = 20
      );

  //! destructor
  virtual ~BEIntegratorWave( );

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

  //! returns element matrix of single layer potential for regular pairs
  void computeElemMatrixHypersingularDisjointP1P1(
      LO outerElem,
      LO innerElem,
      FullMatrix<LO, SC>& matrix
      ) const;

  //! sets current time

  inline void setCurrentTime(
      SC time
      ) {
    currentTime = time;
  }

  //! returns current time

  inline SC getCurrentTime( ) {
    return currentTime;
  }

  //! sets current temporal basis (i) and test (j) functions

  inline void setCurrentFunctions(
      LO i,
      LO j
      ) {
    this->currBasisF = i;
    this->currTestF = j;

  }

  inline void setCurrentLegendreOrder(
      LO orderBasis,
      LO orderTest
      ) {
    this->currLegOrderBasis = orderBasis;
    this->currLegOrderTest = orderTest;
    computeChebyshevForHypersingular( );
    computeChebyshevFor1Layer( );
  }

  //! auxiliary function which generates right hand side
  void getDirichletRHS(
      Vector<LO, SC> &rhs
      ) const;

  //! auxiliary function which generates right hand side
  void getNeumannRHS(
      Vector<LO, SC> &rhs
      ) const;

  //! auxiliary function which generates right hand side
  void getNeumannRHS(
      Vector<LO, SC> &rhs,
      SC( *incWave ) ( SCVT, SCVT*, SCVT* )
      ) const;


  //! evaluates i-th temporal basis function in time t multiplied by the Legendre polynomial of order order
  SCVT evalB(
      SCVT t,
      LO i,
      LO m = 0
      ) const;

  /*
   * evaluates double layer potential in points x, stores values in
   * preallocated vector values
   * 
   * @param[in] x pointer to array with evaluation points
   * @param[in] n number of evaluation points
   * @param[in] density Neumann data
   * @param[out] values preallocated vector for storing results
   */
  void doubleLayerPotential(
      const SCVT *x,
      LO nPoints,
      SCVT t,
      const Vector<LO, SC> & density,
      Vector<LO, SC> & values
      ) const {
    doubleLayerPotentialP1( x, nPoints, t, density, values );
  }



  /*
   * evaluates double layer potential in points x 
   * for p1 density,
   * stores values in preallocated vector values
   * 
   * @param[in] x pointer to array with evaluation points
   * @param[in] n number of evaluation points
   * @param[in] density discretized density function
   * @param[out] values preallocated vector for storing results
   */
  void doubleLayerPotentialP1(
      const SCVT *x,
      LO n,
      SCVT t,
      const Vector<LO, SC> & density,
      Vector<LO, SC> & values
      ) const;

protected:

  //! returns specific kernel evaluated in given points (x, y)

  SC evalDoubleLayerKernel(
      const SCVT *x,
      const SCVT *y,
      const SCVT* n
      ) const {

    // here we have to do some fancy integration!
    std::cout << "Not implemented!" << std::endl;
    return 0.0;
  };

  //! evaluates Legendre polynomial
  SCVT evalLegendrePolynomial(
      SCVT t,
      int order = 0
      ) const;

  //! evaluates the first derivative of Legendre polynomial
  SCVT evalLegendrePolynomialDot(
      SCVT t,
      int order = 0
      ) const;

  //! evaluates the second derivative of Legendre polynomial
  SCVT evalLegendrePolynomialDDot(
      SCVT t,
      int order = 0
      ) const;

  //! evaluates h_a,b function 
  SCVT evalErf(
      SCVT t,
      SCVT a,
      SCVT b
      ) const;

  //! evaluates the i-th partition of unity bump function in time t
  SCVT evalPUM(
      SCVT t,
      LO i
      ) const;

  //! evaluates first time derivative of the h_a,b function
  SCVT evalErfDot(
      SCVT t,
      SCVT a,
      SCVT b
      ) const;

  //! evaluates second time derivative of the h_a,b function
  SCVT evalErfDDot(
      SCVT t,
      SCVT a,
      SCVT b
      ) const;

  // evaluates first time derivative of the i-th PUM function
  SCVT evalPUMDot(
      SCVT t,
      LO i
      ) const;

  // evaluates second time derivative of the i-th PUM function
  SCVT evalPUMDDot(
      SCVT t,
      LO i
      ) const;

  // evaluates time derivative of the i-th temporal basis function
  SCVT evalBDot(
      SCVT t,
      LO i,
      LO m = 0
      ) const;

  // evaluates second time derivative of the i-th temporal basis function
  SCVT evalBDDot(
      SCVT t,
      LO i,
      LO m = 0
      ) const;

  //! current time 
  SCVT currentTime;

  //! current temporal test function
  LO currTestF;

  //! current temporal basis function
  LO currBasisF;

  //! current order of Legendre polynomial
  LO currLegOrderBasis;

  //! current order of Legendre polynomial
  LO currLegOrderTest;

  //! time step
  SCVT dt;

  //! number of time steps
  LO N;

  //! order of temporal Gauss integration
  int timeQuadOrder;

  //! number of Chebyshev intervals for psi interpolation
  int nChebIntervals;

  //! number of Chebyshev points for psi interpolation
  int nChebPoints;
  
  //! chebyshev interval
  double chebDelta;

  //! beginnings of the Chebyshev interpolation intervals
  SCVT *chebIntervalStarts;

  //! coefficients for Cheb. interpolation - pointers to arrays - one for each Cheb. subinterval
  SCVT **psiChebCoeffs;

  SCVT **psiTildeChebCoeffs;

  SCVT **psi1LayerChebCoeffs;

  //! integrates the product of current test and basis function over given interval
  SCVT integrate1Layer(
      SCVT start,
      SCVT end,
      SCVT r
      ) const;

  //! integrates the product of current test and basis function over given interval
  SCVT integrateHypersingularPart2(
      SCVT start,
      SCVT end,
      SCVT r
      ) const;

  //! integrates the product of current test and basis function over given interval
  SCVT integrateHypersingularPart1(
      SCVT start,
      SCVT end,
      SCVT r
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

  //! computes coefficients of Chebyshev interpolation for efficient approximation of function psiddot, psidot
  void computeChebyshevForHypersingular( );

  //! computes coefficients of Chebyshev interpolation for efficient approximation of function psiddot, psidot
  void computeChebyshevFor1Layer( );

  //! computes zeros of Chebyshev interpolant over interval a, b
  void chebyshevZeros(
      SC a,
      SC b,
      SC *zeros
      ) const;

  /*
   * computes coefficients of the Chebyshev interpolant
   * 
   * @param[in]  fx   values of a given function in Cheb. zeros
   * @param[out] coeffs coefficients of the Chebyshev interpolation
   */
  void chebyshevCoefficients(
      SC * fx,
      SC *coeffs
      ) const;

private:

  SC uIncNeu(
      SCVT t,
      SCVT * x,
      SCVT * n,
      int type
      ) const;

  //! computes Chebyshev interpolation in point t, inte
  SC chebyshevInterpolant(
      const SC a,
      const SC b,
      const SC * const coeffs,
      const SC t
      ) const;

  //! returns specific kernel evaluated in given points (x, y)

  SC evalSingleLayerKernel(
      const SCVT *x,
      const SCVT *y
      ) const {

    // here we have to do some fancy integration!

    SCVT norm = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) +
        ( x[1] - y[1] )*( x[1] - y[1] ) +
        ( x[2] - y[2] )*( x[2] - y[2] ) );

    // check whether the supports are overlapping
    SCVT testStart, testEnd;
    SCVT basisStart, basisEnd;
    if ( currTestF == 0 ) {
      testStart = 0.0;
      testEnd = dt;
    } else if ( currTestF == N - 1 ) {
      testStart = ( N - 2 ) * dt;
      testEnd = ( N - 1 ) * dt;
    } else {
      testStart = ( currTestF - 1 ) * dt;
      testEnd = ( currTestF + 1 ) * dt;
    }
    if ( currBasisF == 0 ) {
      basisStart = norm;
      basisEnd = dt + norm;
    } else if ( currBasisF == N - 1 ) {
      basisStart = ( N - 2 ) * dt + norm;
      basisEnd = ( N - 1 ) * dt + norm;
    } else {
      basisStart = ( currBasisF - 1 ) * dt + norm;
      basisEnd = ( currBasisF + 1 ) * dt + norm;
    }

    SCVT intStart = std::min( std::max( testStart, basisStart ),
        ( N - 1 ) * dt );
    SCVT intEnd = std::min( basisEnd, testEnd );

    if ( intStart >= intEnd ) {
      return 0.0;
    }

    SC a = 0;
    SC b = 0;
    int interval;

    for ( interval = 0; interval < nChebIntervals; interval++ ) {
      if ( norm > chebIntervalStarts[interval] && norm <=
          chebIntervalStarts[interval + 1] ) {
        a = chebIntervalStarts[interval];
        b = chebIntervalStarts[interval + 1];
        break;
      }
    }

    ///std::cout << chebyshevInterpolant(a, b, psiTildeChebCoeffs[interval-1], norm) << " , " <<integrateHypersingularPart1( intStart, intEnd, norm ) << std::endl;
    return ( ( PI_FACT / norm ) * chebyshevInterpolant(
        a, b, psi1LayerChebCoeffs[interval], norm ) );

    //return ( ( PI_FACT / norm ) * integrate1Layer( intStart, intEnd, norm ) );
  };

  SC evalSingleLayerKernel(
      SCVT x1,
      SCVT x2,
      SCVT x3,
      SCVT y1,
      SCVT y2,
      SCVT y3
      ) const {

    // here we have to do some fancy integration!

    SCVT norm = std::sqrt( ( x1 - y1 )*( x1 - y1 ) +
        ( x2 - y2 )*( x2 - y2 ) +
        ( x3 - y3 )*( x3 - y3 ) );

    // check whether the supports are overlapping
    SCVT testStart, testEnd;
    SCVT basisStart, basisEnd;
    if ( currTestF == 0 ) {
      testStart = 0.0;
      testEnd = dt;
    } else if ( currTestF == N - 1 ) {
      testStart = ( N - 2 ) * dt;
      testEnd = ( N - 1 ) * dt;
    } else {
      testStart = ( currTestF - 1 ) * dt;
      testEnd = ( currTestF + 1 ) * dt;
    }
    if ( currBasisF == 0 ) {
      basisStart = norm;
      basisEnd = dt + norm;
    } else if ( currBasisF == N - 1 ) {
      basisStart = ( N - 2 ) * dt + norm;
      basisEnd = ( N - 1 ) * dt + norm;
    } else {
      basisStart = ( currBasisF - 1 ) * dt + norm;
      basisEnd = ( currBasisF + 1 ) * dt + norm;
    }

    SCVT intStart = std::min( std::max( testStart, basisStart ),
        ( N - 1 ) * dt );
    SCVT intEnd = std::min( basisEnd, testEnd );

    if ( intStart >= intEnd ) {
      return 0.0;
    }

    SC a = 0;
    SC b = 0;
    int interval;

    for ( interval = 0; interval < nChebIntervals; interval++ ) {
      if ( norm > chebIntervalStarts[interval] && norm <=
          chebIntervalStarts[interval + 1] ) {
        a = chebIntervalStarts[interval];
        b = chebIntervalStarts[interval + 1];
        break;
      }
    }

    ///std::cout << chebyshevInterpolant(a, b, psiTildeChebCoeffs[interval-1], norm) << " , " <<integrateHypersingularPart1( intStart, intEnd, norm ) << std::endl;
    return ( ( PI_FACT / norm ) * chebyshevInterpolant(
        a, b, psi1LayerChebCoeffs[interval], norm ) );

    //return ( ( PI_FACT / norm ) * integrate1Layer( intStart, intEnd, norm ) );
  };


  //! returns specific kernel evaluated in given points (x, y)
  // TODO evaluate laplace kernel separately and multiply later

  SC evalHypersingularPart1(
      const SCVT *x,
      const SCVT *y
      ) const {

    // here we have to do some fancy integration!

    SCVT norm = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) +
        ( x[1] - y[1] )*( x[1] - y[1] ) +
        ( x[2] - y[2] )*( x[2] - y[2] ) );

    // check whether the supports are overlapping
    SCVT testStart, testEnd;
    SCVT basisStart, basisEnd;
    if ( currTestF == 0 ) {
      testStart = 0.0;
      testEnd = dt;
    } else if ( currTestF == N - 1 ) {
      testStart = ( N - 2 ) * dt;
      testEnd = ( N - 1 ) * dt;
    } else {
      testStart = ( currTestF - 1 ) * dt;
      testEnd = ( currTestF + 1 ) * dt;
    }
    if ( currBasisF == 0 ) {
      basisStart = norm;
      basisEnd = dt + norm;
    } else if ( currBasisF == N - 1 ) {
      basisStart = ( N - 2 ) * dt + norm;
      basisEnd = ( N - 1 ) * dt + norm;
    } else {
      basisStart = ( currBasisF - 1 ) * dt + norm;
      basisEnd = ( currBasisF + 1 ) * dt + norm;
    }

    SCVT intStart = std::min( std::max( testStart, basisStart ),
        ( N - 1 ) * dt );
    SCVT intEnd = std::min( basisEnd, testEnd );

    if ( intStart >= intEnd ) {
      return 0.0;
    }

    SC a = 0;
    SC b = 0;
    int interval;

    for ( interval = 0; interval < nChebIntervals; interval++ ) {
      if ( norm > chebIntervalStarts[interval] && norm <=
          chebIntervalStarts[interval + 1] ) {
        a = chebIntervalStarts[interval];
        b = chebIntervalStarts[interval + 1];
        break;
      }
    }
    
    //std::cout << norm << std::endl;
    ///std::cout << chebyshevInterpolant(a, b, psiTildeChebCoeffs[interval-1], norm) << " , " <<integrateHypersingularPart1( intStart, intEnd, norm ) << std::endl;
    return ( ( PI_FACT / norm ) * chebyshevInterpolant(
        a, b, psiTildeChebCoeffs[interval], norm ) );

    //return( integrateHypersingularPart1( intStart, intEnd, norm ) );
    //return ( ( PI_FACT / norm ) * integrateHypersingularPart1( intStart, intEnd, norm ) );
  };
  
//#pragma omp declare simd simdlen( DATA_WIDTH )
  SC evalHypersingularPart1(
      SCVT x1,
      SCVT x2, 
      SCVT x3, 
      SCVT y1,
      SCVT y2, 
      SCVT y3
      ) const {

    SCVT diff1 = x1 - y1;
    SCVT diff2 = x2 - y2;
    SCVT diff3 = x3 - y3;
    SCVT norm = std::sqrt( diff1 * diff1 + diff2 * diff2 + diff3 * diff3 );

    // check whether the supports are overlapping
    SCVT testStart, testEnd;
    SCVT basisStart, basisEnd;
    if ( currTestF == 0 ) {
      testStart = 0.0;
      testEnd = dt;
    } else if ( currTestF == N - 1 ) {
      testStart = ( N - 2 ) * dt;
      testEnd = ( N - 1 ) * dt;
    } else {
      testStart = ( currTestF - 1 ) * dt;
      testEnd = ( currTestF + 1 ) * dt;
    }
    if ( currBasisF == 0 ) {
      basisStart = norm;
      basisEnd = dt + norm;
    } else if ( currBasisF == N - 1 ) {
      basisStart = ( N - 2 ) * dt + norm;
      basisEnd = ( N - 1 ) * dt + norm;
    } else {
      basisStart = ( currBasisF - 1 ) * dt + norm;
      basisEnd = ( currBasisF + 1 ) * dt + norm;
    }

    SCVT intStart = std::min( std::max( testStart, basisStart ),
        ( N - 1 ) * dt );
    SCVT intEnd = std::min( basisEnd, testEnd );

    if ( intStart >= intEnd ) {
      return 0.0;
    }

    SC a = 0;
    SC b = 0;
        
    int interval = (int) ( ( norm - chebIntervalStarts[0] ) / chebDelta );
    a = chebIntervalStarts[ interval ];
    b = chebIntervalStarts[ interval + 1];

    return ( ( PI_FACT / norm ) * chebyshevInterpolant(
        a, b, psiTildeChebCoeffs[interval], norm ) );

  };

  //! returns specific kernel evaluated in given points (x, y)
  // TODO evaluate laplace kernel separately and multiply later

  SC evalHypersingularPart2( 
  const SCVT *x, 
      const SCVT *y 
  ) const {

    // here we have to do some fancy integration!

    SCVT norm = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) +
        ( x[1] - y[1] )*( x[1] - y[1] ) +
        ( x[2] - y[2] )*( x[2] - y[2] ) );

    // check whether the supports are overlapping
    SCVT testStart, testEnd;
    SCVT basisStart, basisEnd;
    if ( currTestF == 0 ) {
      testStart = 0.0;
      testEnd = dt;
    } else if ( currTestF == N - 1 ) {
      testStart = ( N - 2 ) * dt;
      testEnd = ( N - 1 ) * dt;
    } else {
      testStart = ( currTestF - 1 ) * dt;
      testEnd = ( currTestF + 1 ) * dt;
    }
    if ( currBasisF == 0 ) {
      basisStart = norm;
      basisEnd = dt + norm;
    } else if ( currBasisF == N - 1 ) {
      basisStart = ( N - 2 ) * dt + norm;
      basisEnd = ( N - 1 ) * dt + norm;
    } else {
      basisStart = ( currBasisF - 1 ) * dt + norm;
      basisEnd = ( currBasisF + 1 ) * dt + norm;
    }

    SCVT intStart = std::min( std::max( testStart, basisStart ), 
        ( N - 1 ) * dt );
    SCVT intEnd = std::min( basisEnd, testEnd );

    if ( intStart >= intEnd ) {
      return 0.0;
    }

    SC a = 0.0, b = 0.0;
    int interval;

    for ( interval = 0; interval < nChebIntervals; interval++ ) {
      if ( norm > chebIntervalStarts[interval] && norm <=
          chebIntervalStarts[interval + 1] ) {
        a = chebIntervalStarts[interval];
        b = chebIntervalStarts[interval + 1];
        break;
      }
    }
    //std::cout<< interval<<std::endl;
    return ( ( PI_FACT / norm ) * chebyshevInterpolant(
        a, b, psiChebCoeffs[interval], norm ) );


    //return ( ( PI_FACT / norm ) * integrateHypersingularPart2( intStart, intEnd, norm ) );
  };

  //#pragma omp declare simd simdlen( DATA_WIDTH )
  SC evalHypersingularPart2( 
      SCVT x1,
      SCVT x2, 
      SCVT x3, 
      SCVT y1,
      SCVT y2, 
      SCVT y3
  ) const {

    SCVT diff1 = x1 - y1;
    SCVT diff2 = x2 - y2;
    SCVT diff3 = x3 - y3;
    SCVT norm = std::sqrt( diff1 * diff1 + diff2 * diff2 + diff3 * diff3 );
    
    // check whether the supports are overlapping
    SCVT testStart, testEnd;
    SCVT basisStart, basisEnd;
    if ( currTestF == 0 ) {
      testStart = 0.0;
      testEnd = dt;
    } else if ( currTestF == N - 1 ) {
      testStart = ( N - 2 ) * dt;
      testEnd = ( N - 1 ) * dt;
    } else {
      testStart = ( currTestF - 1 ) * dt;
      testEnd = ( currTestF + 1 ) * dt;
    }
    if ( currBasisF == 0 ) {
      basisStart = norm;
      basisEnd = dt + norm;
    } else if ( currBasisF == N - 1 ) {
      basisStart = ( N - 2 ) * dt + norm;
      basisEnd = ( N - 1 ) * dt + norm;
    } else {
      basisStart = ( currBasisF - 1 ) * dt + norm;
      basisEnd = ( currBasisF + 1 ) * dt + norm;
    }

    SCVT intStart = std::min( std::max( testStart, basisStart ), 
        ( N - 1 ) * dt );
    SCVT intEnd = std::min( basisEnd, testEnd );

    if ( intStart >= intEnd ) {
      return 0.0;
    }

    SC a = 0;
    SC b = 0;
        
    int interval = (int) ( ( norm - chebIntervalStarts[0] ) / chebDelta );
    a = chebIntervalStarts[ interval ];
    b = chebIntervalStarts[ interval + 1];

    return ( ( PI_FACT / norm ) * chebyshevInterpolant(
        a, b, psiChebCoeffs[interval], norm ) );

  };
  
  SC evalHypersingular( 
      SCVT x1,
      SCVT x2, 
      SCVT x3, 
      SCVT y1,
      SCVT y2, 
      SCVT y3,
      SCVT & kernel,
      SCVT & kernelDDot
  ) const {

    SCVT diff1 = x1 - y1;
    SCVT diff2 = x2 - y2;
    SCVT diff3 = x3 - y3;
    SCVT norm = std::sqrt( diff1 * diff1 + diff2 * diff2 + diff3 * diff3 );

    // check whether the supports are overlapping
    SCVT testStart, testEnd;
    SCVT basisStart, basisEnd;
    if ( currTestF == 0 ) {
      testStart = 0.0;
      testEnd = dt;
    } else if ( currTestF == N - 1 ) {
      testStart = ( N - 2 ) * dt;
      testEnd = ( N - 1 ) * dt;
    } else {
      testStart = ( currTestF - 1 ) * dt;
      testEnd = ( currTestF + 1 ) * dt;
    }
    if ( currBasisF == 0 ) {
      basisStart = norm;
      basisEnd = dt + norm;
    } else if ( currBasisF == N - 1 ) {
      basisStart = ( N - 2 ) * dt + norm;
      basisEnd = ( N - 1 ) * dt + norm;
    } else {
      basisStart = ( currBasisF - 1 ) * dt + norm;
      basisEnd = ( currBasisF + 1 ) * dt + norm;
    }

    SCVT intStart = std::min( std::max( testStart, basisStart ),
        ( N - 1 ) * dt );
    SCVT intEnd = std::min( basisEnd, testEnd );

    if ( intStart >= intEnd ) {
      return 0.0;
    }

    SC a = 0;
    SC b = 0;
        
    int interval = (int) ( ( norm - chebIntervalStarts[0] ) / chebDelta );
    a = chebIntervalStarts[ interval ];
    b = chebIntervalStarts[ interval + 1];

    kernel = ( ( PI_FACT / norm ) * chebyshevInterpolant(
        a, b, psiTildeChebCoeffs[interval], norm ) );
    
    kernelDDot = ( ( PI_FACT / norm ) * chebyshevInterpolant(
        a, b, psiChebCoeffs[interval], norm ) );

  };
  
};

}

#else

// use experimental basis functions which are not a partition of unity

namespace bem4i {

/*!
 * class for integrators for Sauter-Veit time domain BEM for wave equation
 * 
 */
template<class LO, class SC>
class BEIntegratorWave : public BEIntegrator<LO, SC, BEIntegratorWave<LO, SC> > {
  // to get inner type of complex numbers (for Helmholtz)
  typedef typename GetType<LO, SC>::SCVT SCVT;

  // we have to enable BEIntegrator to use kernel evaluation private methods
  friend class BEIntegrator<LO, SC, BEIntegratorWave<LO, SC> >;

public:

  //! default constructor
  BEIntegratorWave( );

  //! copy constructor
  BEIntegratorWave( const BEIntegratorWave& orig );

  //! constructor taking BESpace as the argument
  BEIntegratorWave( BESpace<LO, SC>* space, int* quadratureOrder, int timeQuadOrder, int nPre = 0, int nPos = 3, int nChebIntervals = 20, int nChebPoints = 20 );

  //! destructor
  virtual ~BEIntegratorWave( );

  //! returns element matrix of single layer potential
  void computeElemMatrix1Layer( LO outerElem, LO innerElem, FullMatrix<LO, SC>& matrix ) const;

  //! returns element matrix of double layer potential
  void computeElemMatrix2Layer( LO outerElem, LO innerElem, FullMatrix<LO, SC>& matrix ) const;

  //! returns element matrix of hypersingular operator
  void computeElemMatrixHypersingular( LO outerElem, LO innerElem, FullMatrix<LO, SC>& matrix ) const;

  //! returns element matrix of hypersingular operator
  void computeElemMatrixHypersingular( LO outerElem, LO innerElem, FullMatrix<LO, SC>& V, FullMatrix<LO, SC>& matrix ) const;

  //! sets current time

  inline void setCurrentTime( SC time ) {
    currentTime = time;
  }

  //! returns current time

  inline SC getCurrentTime( ) {
    return currentTime;
  }

  //! sets current temporal basis (i) and test (j) functions

  inline void setCurrentFunctions( LO i, LO j ) {
    this->currBasisF = i;
    this->currTestF = j;
    computeChebyshevForHypersingular( );
    //computeChebyshevFor1Layer( );
  }

  //! auxiliary function which generates right hand side
  void getDirichletRHS( Vector<LO, SC> &rhs ) const;

  //! auxiliary function which generates right hand side
  void getNeumannRHS( Vector<LO, SC> &rhs ) const;

protected:


  //! returns specific kernel evaluated in given points (x, y)

  SC evalSingleLayerKernel( const SCVT *x, const SCVT *y ) const {

    // here we have to do some fancy integration!

    SCVT norm = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) +
        ( x[1] - y[1] )*( x[1] - y[1] ) +
        ( x[2] - y[2] )*( x[2] - y[2] ) );

    // check whether the supports are overlapping
    SCVT testStart, testEnd;
    SCVT basisStart, basisEnd;
    if ( currTestF == 0 ) {
      testStart = 0.0;
      testEnd = dt;
    } else if ( currTestF == N - 1 ) {
      testStart = ( N - 2 ) * dt;
      testEnd = ( N - 1 ) * dt;
    } else {
      testStart = ( currTestF - 1 ) * dt;
      testEnd = ( currTestF + 1 ) * dt;
    }
    if ( currBasisF == 0 ) {
      basisStart = norm;
      basisEnd = dt + norm;
    } else if ( currBasisF == N - 1 ) {
      basisStart = ( N - 2 ) * dt + norm;
      basisEnd = ( N - 1 ) * dt + norm;
    } else {
      basisStart = ( currBasisF - 1 ) * dt + norm;
      basisEnd = ( currBasisF + 1 ) * dt + norm;
    }

    SCVT intStart = std::min( std::max( testStart, basisStart ), ( N - 1 ) * dt );
    SCVT intEnd = std::min( basisEnd, testEnd );

    if ( intStart >= intEnd ) {
      return 0.0;
    }


    return ( ( PI_FACT / norm ) * integrate1Layer( intStart, intEnd, norm ) );
  };


  //! returns specific kernel evaluated in given points (x, y)
  // TODO evaluate laplace kernel separately and multiply later

  SC evalHypersingularPart1( const SCVT *x, const SCVT *y ) const {

    // here we have to do some fancy integration!

    SCVT norm = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) +
        ( x[1] - y[1] )*( x[1] - y[1] ) +
        ( x[2] - y[2] )*( x[2] - y[2] ) );

    // check whether the supports are overlapping
    SCVT testStart, testEnd;
    SCVT basisStart, basisEnd;

    testStart = ( currTestF - 1 ) * dt;
    testEnd = ( currTestF + 1 ) * dt;


    basisStart = ( currBasisF - 1 ) * dt + norm;
    basisEnd = ( currBasisF + 1 ) * dt + norm;


    SCVT intStart = std::min( std::max( testStart, basisStart ), ( N + nPos ) * dt );
    SCVT intEnd = std::min( basisEnd, testEnd );

    if ( intStart >= intEnd ) {
      return 0.0;
    }

    SC a, b;
    int interval;

    for ( interval = 0; interval < nChebIntervals; interval++ ) {
      if ( norm > chebIntervalStarts[interval] && norm < chebIntervalStarts[interval + 1] ) {
        a = chebIntervalStarts[interval];
        b = chebIntervalStarts[interval + 1];
        break;
      }
    }
    ///std::cout << chebyshevInterpolant(a, b, psiTildeChebCoeffs[interval-1], norm) << " , " <<integrateHypersingularPart1( intStart, intEnd, norm ) << std::endl;
    return ( ( PI_FACT / norm ) * chebyshevInterpolant( a, b, psiTildeChebCoeffs[interval], norm ) );

    //return( integrateHypersingularPart1( intStart, intEnd, norm ) );
    //return ( ( PI_FACT / norm ) * integrateHypersingularPart1( intStart, intEnd, norm ) );
  };

  //! returns specific kernel evaluated in given points (x, y)
  // TODO evaluate laplace kernel separately and multiply later

  SC evalHypersingularPart2( const SCVT *x, const SCVT *y ) const {

    // here we have to do some fancy integration!
    SCVT norm = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) +
        ( x[1] - y[1] )*( x[1] - y[1] ) +
        ( x[2] - y[2] )*( x[2] - y[2] ) );

    // check whether the supports are overlapping
    SCVT testStart, testEnd;
    SCVT basisStart, basisEnd;

    testStart = ( currTestF - 1 ) * dt;
    testEnd = ( currTestF + 1 ) * dt;


    basisStart = ( currBasisF - 1 ) * dt + norm;
    basisEnd = ( currBasisF + 1 ) * dt + norm;


    SCVT intStart = std::min( std::max( testStart, basisStart ), ( N + nPos ) * dt );
    SCVT intEnd = std::min( basisEnd, testEnd );

    if ( intStart >= intEnd ) {
      return 0.0;
    }

    SC a = 0.0, b = 0.0;
    int interval;

    for ( interval = 0; interval < nChebIntervals; interval++ ) {
      if ( norm > chebIntervalStarts[interval] && norm < chebIntervalStarts[interval + 1] ) {
        a = chebIntervalStarts[interval];
        b = chebIntervalStarts[interval + 1];
        break;
      }
    }
    //std::cout<< interval<<std::endl;
    return ( ( PI_FACT / norm ) * chebyshevInterpolant( a, b, psiChebCoeffs[interval], norm ) );


    //return ( ( PI_FACT / norm ) * integrateHypersingularPart2( intStart, intEnd, norm ) );
  };


  //! returns specific kernel evaluated in given points (x, y)

  SC evalDoubleLayerKernel( const SCVT *x, const SCVT *y, const SCVT* n ) const {

    // here we have to do some fancy integration!
    std::cout << "Not implemented!" << std::endl;
    return 0.0;
  };

  SCVT evalLegendrePolynomial( SCVT t, int order = 0 ) const {
    // for now we just have polynomial of 0-th order
    return 1.0;
  }

  //! evaluates h_a,b function 
  SCVT evalErf( SCVT t, SCVT a, SCVT b ) const;

  //! evaluates the i-th partition of unity bump function in time t
  SCVT evalPUM( SCVT t, LO i ) const;

  //! evaluates i-th temporal basis function in time t multiplied by the Legendre polynomial of order order
  SCVT evalB( SCVT t, LO i, LO m = 0 ) const;

  //! evaluates first time derivative of the h_a,b function
  SCVT evalErfDot( SCVT t, SCVT a, SCVT b ) const;

  //! evaluates second time derivative of the h_a,b function
  SCVT evalErfDDot( SCVT t, SCVT a, SCVT b ) const;

  // evaluates first time derivative of the i-th PUM function
  SCVT evalPUMDot( SCVT t, LO i ) const;

  // evaluates second time derivative of the i-th PUM function
  SCVT evalPUMDDot( SCVT t, LO i ) const;

  // evaluates time derivative of the i-th temporal basis function
  SCVT evalBDot( SCVT t, LO i, LO m = 0 ) const;

  // evaluates second time derivative of the i-th temporal basis function
  SCVT evalBDDot( SCVT t, LO i, LO m = 0 ) const;

  //! current time 
  SCVT currentTime;

  //! current temporal test function
  LO currTestF;

  //! current temporal basis function
  LO currBasisF;

  //! time step
  SCVT dt;

  //! number of time steps
  LO N;

  //! order of temporal Gauss integration
  int timeQuadOrder;

  //! number of Chebyshev intervals for psi interpolation
  int nChebIntervals;

  //! number of Chebyshev points for psi interpolation
  int nChebPoints;

  //! number of pre points
  int nPre;

  //! number of post points
  int nPos;

  //! beginnings of the Chebyshev interpolation intervals
  SCVT *chebIntervalStarts;

  //! coefficients for Cheb. interpolation - pointers to arrays - one for each Cheb. subinterval
  SCVT **psiChebCoeffs;

  SCVT **psiTildeChebCoeffs;

  //! integrates the product of current test and basis function over given interval
  SCVT integrate1Layer( SCVT start, SCVT end, SCVT r ) const;

  //! integrates the product of current test and basis function over given interval
  SCVT integrateHypersingularPart2( SCVT start, SCVT end, SCVT r ) const;

  //! integrates the product of current test and basis function over given interval
  SCVT integrateHypersingularPart1( SCVT start, SCVT end, SCVT r ) const;

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

  //! computes coefficients of Chebyshev interpolation for efficient approximation of function psiddot, psidot
  void computeChebyshevForHypersingular( );

  //! computes zeros of Chebyshev interpolant over interval a, b
  void chebyshevZeros( SC a, SC b, SC *zeros ) const;

  /*
   * computes coefficients of the Chebyshev interpolant
   * 
   * @param[in]  fx   values of a given function in Cheb. zeros
   * @param[out] coeffs coefficients of the Chebyshev interpolation
   */
  void chebyshevCoefficients( SC *fx, SC *coeffs ) const;

  //! computes Chebyshev interpolation in point t, inte
  SC chebyshevInterpolant( SC a, SC b, SC *coeffs, SC t ) const;
};

}

#endif
// include .cpp file to overcome linking problems due to templates
#include "BEIntegratorWave.cpp"

#endif /* BEINTEGRATORWAVE_H */

