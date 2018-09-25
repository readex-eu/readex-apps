/*!
 * @file    BEIntegratorLame.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    April 29, 2015
 * @brief   Header file for class BEIntegratorLame
 * 
 */

#ifndef BEINTEGRATORLAME_H
#define	BEINTEGRATORLAME_H

#include "BEIntegrator.h"

namespace bem4i {

template <class LO, class SC>
class BEIntegratorLame : public BEIntegrator<LO, SC, BEIntegratorLame<LO, SC> > {
  // to get inner type of complex numbers (for Helmholtz)
  typedef typename GetType<LO, SC>::SCVT SCVT;

  // we have to enable BEIntegrator to use kernel evaluation private methods  
  friend class BEIntegrator<LO, SC, BEIntegratorLame<LO, SC> >;

public:

  //! constructor taking BESpace as the argument
  BEIntegratorLame(
      BESpace<LO, SC>* space,
      int* quadratureOrder,
      quadratureType quadrature = SauterSchwab,
      int* quadratureOrderDisjointElems = nullptr
      );

  virtual ~BEIntegratorLame( );

  //! sets the current row in the aux. matrix V 

  inline void setI( int i ) {
    this->subi = i;
  }

  //! sets the current column in the aux. matrix V 

  inline void setJ( int j ) {
    this->subj = j;
  }

  //! returns element matrix of single layer potential
  void computeElemMatrix1Layer(
      LO outerElem,
      LO innerElem,
      FullMatrix<LO, SC>& matrix
      ) const;

#if N_MIC > 0
  __attribute__( ( target( mic ) ) )
#endif
  static void computeElemMatrixAllMIC(
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
      );

protected:

  int subi, subj;

  //! computes members of the element matrix and stores them in matrix in order:
  //! V_laplace, V11, V22, V33, V12, V13, V23
  void computeElemMatrix1LayerSauterSchwabP0P0(
      LO outerElem,
      LO innerElem,
      FullMatrix<LO, SC>& matrix
      ) const;

  void computeElemMatrix1LayerDisjointP0P0(
      LO outerElem,
      LO innerElem,
      FullMatrix<LO, SC>& matrix
      ) const;
  
  void computeElemMatrix1LayerSauterSchwabP1P1(
      LO outerElem,
      LO innerElem,
      FullMatrix<LO, SC>& matrix
      ) const;

  void computeElemMatrix1LayerDisjointP1P1(
      LO outerElem,
      LO innerElem,
      FullMatrix<LO, SC>& matrix
      ) const;

private:

  BEIntegratorLame( );

  BEIntegratorLame(
      const BEIntegratorLame& orig
      );

  // just common part!!!

  SC evalSingleLayerKernel(
      const SCVT *x,
      const SCVT *y
      ) const {
    SCVT norm = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) +
        ( x[1] - y[1] )*( x[1] - y[1] ) +
        ( x[2] - y[2] )*( x[2] - y[2] ) );

    //    return ( PI_FACT * ( ( ( x[subi] - y[subi] )*( x[subj] - y[subj] ) ) /
    //        ( norm * norm * norm ) ) );
    return (SCVT) 1.0 / norm;
  };

  // just common part!!!

  SC evalSingleLayerKernel(
      SCVT x1,
      SCVT x2,
      SCVT x3,
      SCVT y1,
      SCVT y2,
      SCVT y3
      ) const {

    SCVT norm = std::sqrt( ( x1 - y1 ) * ( x1 - y1 ) +
        ( x2 - y2 ) * ( x2 - y2 ) + ( x3 - y3 ) * ( x3 - y3 ) );

    return (SCVT) 1.0 / norm;
  };

};

}

#include "BEIntegratorLame.cpp"

#endif	/* BEINTEGRATORLAME_H */

