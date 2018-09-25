/*!
 * @file    BEBilinearFormLaplaceHypersingular.h
 * @author  Jan Zapletal
 * @date    July 30, 2013
 * @brief   Header file for class BEBilinearFormLaplaceHypersingular
 * 
 */

#ifndef BEBILINEARFORMLAPLACEHYPERSINGULAR_H
#define	BEBILINEARFORMLAPLACEHYPERSINGULAR_H

#include "BEBilinearForm.h"
#include "LaplaceHypersingularOperator.h"

namespace bem4i {

/*! 
 * Class representing the bilinear form for the Laplace hypersingular operator
 * 
 * Provides methods for system matrix assembly.
 * 
 */
template<class LO, class SC>
class BEBilinearFormLaplaceHypersingular : public BEBilinearForm<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  using BEBilinearForm<LO, SC>::assemble;

  //! default constructor
  BEBilinearFormLaplaceHypersingular( );

  //! copy constructor
  BEBilinearFormLaplaceHypersingular(
      const BEBilinearFormLaplaceHypersingular & orig
      );

  /*! 
   * Constructor taking the boundary element space on which the form acts
   */
  BEBilinearFormLaplaceHypersingular(
      BESpace<LO, SC> * space,
      int * quadratureOrder = nullptr,
      quadratureType quadrature = SauterSchwab,
      int * quadratureOrderDisjointElems = nullptr
      );

  //! destructor
  virtual ~BEBilinearFormLaplaceHypersingular( );

  /*!
   * The method assembles the Galerkin matrix for the Laplace hypersingular operator
   */
  virtual void assemble(
      FullMatrix<LO, SC> & matrix
      ) const;

  virtual void assembleP1P1(
      FullMatrix<LO, SC> & matrix
      ) const;

  /*!
   * The method assembles the operator T (for T^T*V*T)
   */
  virtual void assemble(
      LaplaceHypersingularOperator<LO, SC> & op
      ) const;
/*
  virtual void assemble(
      FullMatrix<LO, SC> & matrix,
      FullMatrix<LO, SC> & V
      ) const;
*/
protected:

  //! method assembles idx-th row of a local block

  virtual void assembleRow(
      const BEBlockCluster< LO, SC > & block,
      LO idx,
      const std::vector< LO > & cols,
      Vector<LO, SC>& row,
      void * voidIntegrator
      ) const;

  void assembleRowP1DisP1Dis(
      const BEBlockCluster< LO, SC > & block,
      LO idx,
      const std::vector< LO > & cols,
      Vector<LO, SC>& row,
      void * voidIntegrator
      ) const;

  //! method assembles idx-th row of a local block

  virtual void assembleColumn(
      const BEBlockCluster< LO, SC > & block,
      LO idx,
      const std::vector< LO > & rows,
      Vector<LO, SC>& col,
      void * voidIntegrator
      ) const;

  void assembleColumnP1DisP1Dis(
      const BEBlockCluster< LO, SC > & block,
      LO idx,
      const std::vector< LO > & rows,
      Vector<LO, SC>& col,
      void * voidIntegrator
      ) const;

  /*!
   * Submatrix of a full Galerkin matrix for single layer Laplace operator
   * 
   * @param[in]   rowElems  row indices
   * @param[in]   colElems  column indices
   * 
   */
  virtual void assemble(
      BECluster<LO, SC> const *leftCluster,
      BECluster<LO, SC> const *rightCluster,
      FullMatrix<LO, SC>& matrix,
      void * voidIntegrator
      ) const;

  void assembleP1DisP1Dis(
      const BECluster< LO, SC > * leftCluster,
      const BECluster< LO, SC > * rightCluster,
      FullMatrix< LO, SC > & matrix,
      void * voidIntegrator
      ) const;

  virtual void * createIntegrator( ) const {
    return static_cast<void *> ( new BEIntegratorLaplace< LO, SC >(
        this->space, this->quadratureOrder, this->quadrature,
        this->quadratureOrderDisjointElems ) );
  }

  virtual void destroyIntegrator(
      void * voidIntegrator
      ) const {
    delete static_cast<BEIntegratorLaplace< LO, SC > *> ( voidIntegrator );
  }

private:

  SC x1[3], x2[3], x3[3];

  Vector<LO, SC> *auxCurl;

  inline void initCurl( ) {
    this->auxCurl = this->space->getMesh( )->getCurls( );
  }

  //! default quadrature rule for this class (specific values set in .cpp file)
  static int defaultQuadratureOrder[2];

};

}

// include .cpp file to overcome linking problems due to templates
#include "BEBilinearFormLaplaceHypersingular.cpp"

#endif	/* BEBILINEARFORMLAPLACEHYPERSINGULAR_H */

