/*!
 * @file    MathFun.h
 * @author  Jan Zapletal
 * @date    February 20, 2015
 * @brief   Various mathematical functions
 * 
 */

#ifndef MATHFUN_H
#define	MATHFUN_H

#include <vector>
#include <math.h>
#include "omp.h"

namespace bem4i {

template<class LO, class SC>
class MathFun {

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  MathFun( ) = delete;

  MathFun( const MathFun & orig ) = delete;

  /*!
   * Evaluates associated Legendre polynomials in x (P_l^m)
   * for maxL = l creating values P_0^0, 
   *                              P_1^0, P_1^1, 
   *                              P_2^0, P_2^1, P_2^2, 
   *                              ...
   *                              P_l^0 ... P_l^l
   * I.e., P_l^m is stored at position 0.5*l*(l+1)+m
   * 
   * @param   x      evaluation point in (-1,1)
   * @param   maxL   maximum degree
   * @param   res    vector of results (0.5*(maxL+1)*(maxL+2))
   */
  static void legendre(
      SCVT x,
      LO maxL,
      std::vector< SCVT > & res
      );

  //  /*!
  //   * Evaluates the real basis of spherical harmonics
  //   * for maxL = l creating values Y_0^0, 
  //   *                              Y_1^0, Y_1^1, 
  //   *                              Y_2^0, Y_2^1, Y_2^2, 
  //   *                              ...
  //   *                              Y_l^0 ... Y_l^l
  //   * I.e., Y_l^m is stored at position 0.5*l*(l+1)+m
  //   * 
  //   * @param   theta  ranges from 0 at the North Pole to pi at the South Pole
  //   * @param   phi    longitude ranging from 0 to 2pi 
  //   * @param   maxL   maximum degree
  //   * @param   res    pre-allocated array for results (0.5*(maxL+1)*(maxL+2))
  //   */
  //  static void sphericalHarmonicsRealBasis(
  //      SCVT theta,
  //      SCVT phi,
  //      LO maxL,
  //      SCVT * res
  //      );

  /*!
   * Evaluates the real basis of spherical harmonics
   * for maxL = l creating values Y_0^0, 
   *                              Y_1^-1, Y_1^0,  Y_1^1, 
   *                              ...
   *                              Y_l^-l ... Y_l^l
   * I.e., Y_l^m is stored at position l*(l+1)+m
   * 
   * @param   x      point on a unit sphere
   * @param   maxL   maximum degree
   * @param   res    array of results (maxL+1)^2
   */
  static void sphericalHarmonicsRealBasis(
      const SCVT * x,
      LO maxL,
      std::vector< SCVT > & res
      );

  /*!
   * Evaluates the real basis of spherical harmonics
   * for maxL = l creating values Y_0^0, 
   *                              Y_1^-1, Y_1^0,  Y_1^1, 
   *                              ...
   *                              Y_l^-l ... Y_l^l
   * I.e., Y_l^m is stored at position l*(l+1)+m
   * 
   * @param   x        points on a unit sphere
   * @param   nPoints  number of points in x
   * @param   maxL     maximum degree
   * @param   res      entry in res corresponds to Y_l^m evaluated in every x
   */
  static void sphericalHarmonicsRealBasis(
      const SCVT * x,
      LO nPoints,
      LO maxL,
      std::vector< Vector< LO, SCVT > * > & res
      );

private:

  static SCVT legendreIncrementM(
      LO l,
      LO m,
      SCVT x,
      SCVT p_l_m,
      SCVT p_lm_m,
      SCVT sqroot
      );

  static SCVT legendreIncrementL(
      LO l,
      LO m,
      SCVT x,
      SCVT p_l_m,
      SCVT p_lm_m
      );

  static LO lmtoiP(
      LO l,
      LO m
      ) {
    return (LO) ( 0.5 * l * ( l + 1 ) + m );
  }

  static LO lmtoiY(
      LO l,
      LO m
      ) {
    return l * ( l + 1 ) + m;
  }

};

}

#include "MathFun.cpp"

#endif	/* MATHFUN_H */
