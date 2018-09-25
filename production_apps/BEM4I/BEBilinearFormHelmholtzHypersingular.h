/*!
 * @file    BEBilinearFormHelmholtzHypersingular.h
 * @author  Jan Zapletal
 * @date    August 19, 2013
 * @brief   Header file for class BEBilinearFormHelmholtzHypersingular
 * 
 */

#ifndef BEBILINEARFORMHELMHOLTZHYPERSINGULAR_H
#define	BEBILINEARFORMHELMHOLTZHYPERSINGULAR_H

#include "BEBilinearForm.h"

namespace bem4i {

/*! 
 * Class representing the bilinear form for the Helmholtz hypersingular operator
 * 
 * Provides methods for system matrix assembly.
 * 
 */
template<class LO, class SC>
class BEBilinearFormHelmholtzHypersingular : public BEBilinearForm<LO, SC> {
public:

  using BEBilinearForm<LO, SC>::assemble;

  typedef typename GetType<LO, SC>::SCVT SCVT;
  //! default constructor
  BEBilinearFormHelmholtzHypersingular( );

  //! copy constructor
  BEBilinearFormHelmholtzHypersingular( const BEBilinearFormHelmholtzHypersingular& orig );

  /*! 
   * Constructor taking the boundary element space on which the form acts
   */
  BEBilinearFormHelmholtzHypersingular(
      BESpace<LO, SC> * space,
      int * quadratureRule = nullptr,
      SC kappa = 2.0,
      quadratureType quadrature = SauterSchwab,
      int* quadratureOrderDisjointElems = nullptr
      );

  //! destructor
  virtual ~BEBilinearFormHelmholtzHypersingular( );

  /*!
   * The method assembles the Galerkin matrix for the Helmholtz hypersingular operator
   */
  virtual void assemble(
      FullMatrix<LO, SC> &
      matrix ) const;

  virtual void assembleP1P1(
      FullMatrix<LO, SC> &
      matrix ) const;

  /*!
   * The method assembles the Galerkin matrix for the Helmholtz hypersingular operator
   */
  //virtual void assemble( FullMatrix<LO, SC>& matrix, FullMatrix<LO, SC>& V );

  //! assembles 4th matrix from p0p0 formulation
  virtual void assembleP0P0(
      FullMatrix<LO, SC>& matrix
      ) const;

  //! assembles 4th matrix from p0p0 formulation
  virtual void assembleH1P0P0(
      FullMatrix<LO, SC>& matrix
      ) const;
  //! assembles 4th matrix from p0p0 formulation
  virtual void assembleH2P0P0(
      FullMatrix<LO, SC>& matrix
      ) const;

  //! assembles 4th matrix from p0p0 formulation
  virtual void assembleH1P0P0(
      ACAMatrix<LO, SC>& matrix
      ) const;
  //! assembles 4th matrix from p0p0 formulation

  virtual void assembleH2P0P0(
      ACAMatrix<LO, SC>& matrix
      ) const;

  //! assembles 4th matrix from p0p0 formulation
  virtual void assembleH1P0P0(
      MPIACAMatrix<LO, SC>& matrix
      ) const;
  //! assembles 4th matrix from p0p0 formulation

  virtual void assembleH2P0P0(
      MPIACAMatrix<LO, SC>& matrix
      ) const;

  virtual void assembleACABlockH1P0P0(
      BEBlockCluster<LO, SC>* block,
      FullMatrix<LO, SC>& U,
      FullMatrix<LO, SC>& V,
      void * voidIntegrator
      ) const;

  virtual void assembleACABlockH2P0P0(
      BEBlockCluster<LO, SC>* block,
      FullMatrix<LO, SC>& U,
      FullMatrix<LO, SC>& V,
      void * voidIntegrator
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
   * Submatrix of a full Galerkin matrix for hypersingular Helmholtz operator
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


  // methods for assembling p0p0 matrices

  virtual void assembleH1P0P0(
      BECluster<LO, SC> const *leftCluster,
      BECluster<LO, SC> const *rightCluster,
      FullMatrix<LO, SC>& matrix
      ) const;

  virtual void assembleH2P0P0(
      BECluster<LO, SC> const *leftCluster,
      BECluster<LO, SC> const *rightCluster,
      FullMatrix<LO, SC>& matrix
      ) const;

  virtual void assembleRowH1P0P0(
      const BEBlockCluster< LO, SC > & block,
      LO idx,
      const std::vector< LO > & cols,
      Vector<LO, SC>& row,
      void * voidIntegrator
      ) const;

  virtual void assembleColumnH1P0P0(
      const BEBlockCluster< LO, SC > & block,
      LO idx,
      const std::vector< LO > & rows,
      Vector<LO, SC>& col,
      void * voidIntegrator
      ) const;

  virtual void assembleRowH2P0P0(
      const BEBlockCluster< LO, SC > & block,
      LO idx,
      const std::vector< LO > & cols,
      Vector<LO, SC>& row,
      void * voidIntegrator
      ) const;

  virtual void assembleColumnH2P0P0(
      const BEBlockCluster< LO, SC > & block,
      LO idx,
      const std::vector< LO > & rows,
      Vector<LO, SC>& col,
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

  Vector<LO, SCVT> *auxCurl;

  inline void initCurl( ) {
    this->auxCurl = this->space->getMesh( )->getCurls( );
  }

  //! default quadrature rule for this class (specific values set in .cpp file)
  static int defaultQuadratureOrder[2];

  SC kappa;

};

}

// include .cpp file to overcome linking problems due to templates
#include "BEBilinearFormHelmholtzHypersingular.cpp"

#endif	/* BEBILINEARFORMHELMHOLTZHYPERSINGULAR_H */

