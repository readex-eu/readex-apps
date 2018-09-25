/*!
 * @file    BEBilinearFormLame1Layer.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    April 29, 2015
 * @brief   Header file for class BEBilinearFormLame1Layer
 * 
 */

#ifndef BEBILINEARFORMLAME1LAYER_H
#define	BEBILINEARFORMLAME1LAYER_H

#include "BEBilinearForm.h"
#include "BEIntegratorLame.h"

namespace bem4i {

template<class LO, class SC>
class BEBilinearFormLame1Layer : public BEBilinearForm<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  using BEBilinearForm<LO, SC>::assemble;

  BEBilinearFormLame1Layer(
      BESpace<LO, SC>* space,
      int* quadratureOrder = nullptr,
      quadratureType quadrature = SauterSchwab,
      int* quadratureOrderDisjointElems = nullptr,
      bool employSymmetricity = true
      );

  virtual ~BEBilinearFormLame1Layer( );

  /*!
   * The method assembles the Galerkin matrix for the Lame single layer operator
   */
  virtual void assemble(
      FullMatrix<LO, SC>& matrix
      ) const {
  };


  /*!
   * The method assembles the Galerkin matrix for the Lame single layer operator
   */
  virtual void assemble(
      FullMatrix<LO, SC>& matrix,
      FullMatrix<LO, SC>& V
      ) const;

  virtual void assembleP0P0(
      FullMatrix<LO, SC>& matrix,
      FullMatrix<LO, SC>& V
      ) const;

  virtual void assembleP1P1(
      FullMatrix<LO, SC>& matrix,
      FullMatrix<LO, SC>& V
      ) const;

  inline void setNu( SCVT nu ) {
    this->nu = nu;
  }

  inline void setE( SCVT E ) {
    this->E = E;
  }

  inline SCVT getNu( ) {
    return this->nu;
  }

  inline SCVT getE( ) {
    return this->E;
  }

  /*!
   * The method assembles the Galerkin matrices for the Lame operator
   * using the Intel Xeon Phi coprocessors
   */
  virtual void assembleAllMIC(
      FullMatrix<LO, SC>& VLame,
      FullMatrix<LO, SC>& VLaplace,
      FullMatrix<LO, SC>& KLaplace,
      BESpace< LO, SC > & bespaceK
      ) const;

  /*!
   * The method assembles the Galerkin matrices for the Lame operator
   * using the Intel Xeon Phi coprocessors
   */
  virtual void assembleAllMIC(
      BlockMatrix<LO, SC>& VLame,
      FullMatrix<LO, SC>& VLaplace,
      FullMatrix<LO, SC>& KLaplace,
      BESpace< LO, SC > & bespaceK
      ) const;

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

  SCVT nu;

  SCVT E;

  bool employSymmetricity;

private:
  BEBilinearFormLame1Layer( );
  BEBilinearFormLame1Layer(
      const BEBilinearFormLame1Layer& orig
      );


  void getGlobalIndices(
      LO locRows,
      LO locCols,
      const FullMatrix<LO, SC> & values,
      std::vector<LO> & glRows,
      std::vector<LO> & glCols,
      std::vector<SC> & glValues
      ) const;

  void getGlobalIndices(
      LO locRows,
      LO locCols,
      const FullMatrix<LO, SC> & values,
      LO * glRows,
      LO * glCols,
      SC * glValues
      ) const;

};

}

// include .cpp file to overcome linking problems due to templates
#include "BEBilinearFormLame1Layer.cpp"

#endif	/* BEBILINEARFORMLAME1LAYER_H */

