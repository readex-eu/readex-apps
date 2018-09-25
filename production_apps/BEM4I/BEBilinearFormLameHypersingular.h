/*!
 * @file    BEBilinearFormLameHypersingular.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    April 30, 2015
 * @brief   Header file for class BEBilinearFormLameHypersingular
 * 
 */

#ifndef BEBILINEARFORMLAMEHYPERSINGULAR_H
#define	BEBILINEARFORMLAMEHYPERSINGULAR_H

#include "BEBilinearForm.h"

namespace bem4i {

template <class LO, class SC>
class BEBilinearFormLameHypersingular : public BEBilinearForm<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:
  
  using BEBilinearForm<LO, SC>::assemble;

  virtual ~BEBilinearFormLameHypersingular( );

  BEBilinearFormLameHypersingular(
      BESpace<LO, SC>* space,
      int* quadratureOrder = nullptr,
      quadratureType quadrature = SauterSchwab,
      int* quadratureOrderDisjointElems = nullptr
      );

  /*!
   * The method assembles the Galerkin matrix for the Lame hypersingular 
   * operator
   */
  virtual void assemble(
      FullMatrix<LO, SC>& matrix
      ) const {
  };

  /*!
   * The method assembles the Galerkin matrix for the Lame hypersingular 
   * operator
   */
  virtual void assemble(
      FullMatrix<LO, SC>& matrix,
      FullMatrix<LO, SC> & Vlaplace,
      FullMatrix<LO, SC> & Vlame,
      const std::vector< SparseMatrix<LO, SC>* > & T
      ) const;

  virtual void assembleP1P1(
      FullMatrix<LO, SC>& matrix,
      FullMatrix<LO, SC> & Vlaplace,
      FullMatrix<LO, SC> & Vlame,
      const std::vector< SparseMatrix<LO, SC>* > & T
      ) const;

  inline void setNu( SCVT nu ) {
    this->nu = nu;
    this->mu = this->E / ( 2 * ( 1 + this->nu ) );
  }

  inline void setE( SCVT E ) {
    this->E = E;
    this->mu = this->E / ( 2 * ( 1 + this->nu ) );
  }

  inline SCVT getNu( ) {
    return this->nu;
  }

  inline SCVT getE( ) {
    return this->E;
  }

protected:
  //! method assembles idx-th row of a local block

  virtual void assembleRow(
      const BEBlockCluster< LO, SC > & block,
      LO idx,
      const std::vector< LO > & cols,
      Vector<LO, SC>& row,
      void * voidIntegrator
      ) const {
  };

  //! method assembles idx-th row of a local block

  virtual void assembleColumn(
      const BEBlockCluster< LO, SC > & block,
      LO idx,
      const std::vector< LO > & rows,
      Vector<LO, SC>& col,
      void * voidIntegrator
      ) const {
  };

  /*!
   * Submatrix of a full Galerkin matrix for single layer Lame operator
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
      ) const {
  };

  virtual void * createIntegrator( ) const {
    return static_cast<void *> ( new BEIntegratorLame< LO, SC >(
        this->space, this->quadratureOrder, this->quadrature,
        this->quadratureOrderDisjointElems ) );
  }

  virtual void destroyIntegrator(
      void * voidIntegrator
      ) const {
    delete static_cast<BEIntegratorLame< LO, SC > *> ( voidIntegrator );
  }

private:

  BEBilinearFormLameHypersingular( );

  BEBilinearFormLameHypersingular(
      const BEBilinearFormLameHypersingular& orig
      );

  SCVT nu;

  SCVT E;

  SCVT mu;

};

}

#include "BEBilinearFormLameHypersingular.cpp"

#endif	/* BEBILINEARFORMLAMEHYPERSINGULAR_H */

