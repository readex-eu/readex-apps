/*!
 * @file    BEBilinearFormLame2Layer.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    April 29, 2015
 * @brief   Header file for class BEBilinearFormLame2Layer
 * 
 */

#ifndef BEBILINEARFORMLAME2LAYER_H
#define	BEBILINEARFORMLAME2LAYER_H


#include "BEBilinearForm.h"
#include "BEIntegratorLame.h"

namespace bem4i {

template<class LO, class SC>
class BEBilinearFormLame2Layer : public BEBilinearForm<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  using BEBilinearForm<LO, SC>::assemble;

  BEBilinearFormLame2Layer(
      BESpace<LO, SC>* space,
      int* quadratureOrder = nullptr,
      quadratureType quadrature = SauterSchwab,
      int* quadratureOrderDisjointElems = nullptr
      );


  virtual ~BEBilinearFormLame2Layer( );

  /*!
   * The method assembles the Galerkin matrix for the Lame double layer operator
   */
  virtual void assemble(
      FullMatrix<LO, SC>& matrix
      ) const {
  };

  /*!
   * The method assembles the Galerkin matrix for the Lame single layer operator
   */
  virtual void assemble(
      FullMatrix<LO, SC> & matrix,
      FullMatrix<LO, SC> & Vlaplace,
      FullMatrix<LO, SC> & Vlame,
      FullMatrix<LO, SC> & Klaplace,
      const std::vector< SparseMatrix<LO, SC> * > & T
      ) const;

  virtual void assembleP0P1(
      FullMatrix<LO, SC> & matrix,
      FullMatrix<LO, SC> & Vlaplace,
      FullMatrix<LO, SC> & Vlame,
      FullMatrix<LO, SC> & Klaplace,
      const std::vector< SparseMatrix<LO, SC> * > & T
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

  BEBilinearFormLame2Layer( );

  BEBilinearFormLame2Layer(
      const BEBilinearFormLame2Layer& orig
      );


  SCVT nu;

  SCVT E;

};

}

#include "BEBilinearFormLame2Layer.cpp"

#endif	/* BEBILINEARFORMLAME2LAYER_H */

