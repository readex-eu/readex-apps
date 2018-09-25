/*!
 * @file    BEBilinearFormLaplace1Layer.h
 * @author  Jan Zapletal
 * @date    August 12, 2013
 * @brief   Header file the class BEBilinearFormHelmholtz1Layer
 * 
 */

#ifndef BEBILINEARFORMHELMHOLTZ1LAYER_H
#define	BEBILINEARFORMHELMHOLTZ1LAYER_H

#include "BEBilinearForm.h"

namespace bem4i {

/*! 
 * Class representing the bilinear form for the Helmholtz single layer operator
 * 
 * Provides methods for system matrix assembly.
 * 
 */
template<class LO, class SC>
class BEBilinearFormHelmholtz1Layer : public BEBilinearForm<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  using BEBilinearForm<LO, SC>::assemble;

  //! default constructor
  BEBilinearFormHelmholtz1Layer( );

  //! copy constructor
  BEBilinearFormHelmholtz1Layer( const BEBilinearFormHelmholtz1Layer& orig );

  /*! 
   * Constructor taking the boundary element space on which the form acts
   */
  BEBilinearFormHelmholtz1Layer(
      BESpace<LO, SC>* space,
      int* quadratureOrder = nullptr,
      SC kappa = 2.0,
      quadratureType quadrature = SauterSchwab,
      int* quadratureOrderDisjointElems = nullptr
      );

  //! destructor
  virtual ~BEBilinearFormHelmholtz1Layer( );

  /*!
   * The method assembles the Galerkin matrix for the Helmholtz single layer operator
   */
  virtual void assemble(
      FullMatrix<LO, SC> & matrix
      ) const;

  virtual void assembleP0P0(
      FullMatrix<LO, SC> & matrix
      ) const;

  virtual void assembleP1P1(
      FullMatrix<LO, SC> & matrix
      ) const;

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

  void assembleRowP0P0(
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

  void assembleColumnP0P0(
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
   * Submatrix of a full Galerkin matrix for single layer Helmholtz operator
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
  
  void assembleP0P0(
      const BECluster< LO, SC > * leftCluster,
      const BECluster< LO, SC > * rightCluster,
      FullMatrix< LO, SC > & matrix,
      void * voidIntegrator
      ) const;

  void assembleP1DisP1Dis(
      const BECluster< LO, SC > * leftCluster,
      const BECluster< LO, SC > * rightCluster,
      FullMatrix< LO, SC > & matrix,
      void * voidIntegrator
      ) const;

  virtual void * createIntegrator( ) const {
    return static_cast<void *> ( new BEIntegratorHelmholtz< LO, SC >(
        this->space, this->quadratureOrder, this->kappa, this->quadrature,
        this->quadratureOrderDisjointElems ) );
  }

  virtual void destroyIntegrator(
      void * voidIntegrator
      ) const {
    delete static_cast<BEIntegratorHelmholtz< LO, SC > *> ( voidIntegrator );
  }

private:

  //! default quadrature rule for this class (specific values set in .cpp file)
  static int defaultQuadratureOrder[2];

  //! kappa
  SC kappa;

};

}

// include .cpp file to overcome linking problems due to templates
#include "BEBilinearFormHelmholtz1Layer.cpp"

#endif	/* BEBILINEARFORMHELMHOLTZ1LAYER_H */
