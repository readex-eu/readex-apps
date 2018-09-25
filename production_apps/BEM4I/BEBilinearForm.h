/*!
 * @file    BEBilinearForm.h
 * @author  Michal Merta 
 * @date    July 12, 2013
 * @brief   Header file for pure virtual class BEBilinearForm
 * 
 */

#ifndef BEBILINEARFORM_H
#define	BEBILINEARFORM_H

#include <algorithm>
#include <omp.h>
#include <new>

#include "FullMatrix.h"
#include "FastBESpace.h"
#include "Settings.h"
#include "Macros.h"
#include "Quadratures.h"
#include "MPIACAMatrix.h"
#include "MPIBlockACAMatrix.h"
#include "ProgressMonitor.h"


namespace bem4i {

//! type of quadrature

enum quadratureType {
  SauterSchwab,
  Steinbach
};

/*! 
 * Abstract class representing a bilinear form over some BEM space
 * 
 */
template<class LO, class SC>
class BEBilinearForm {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:


  //! default constructor
  BEBilinearForm( );

  //! copy constructor
  BEBilinearForm(
      const BEBilinearForm& orig
      );

  //! destructor
  virtual ~BEBilinearForm( );

  //! method assembles a full matrix of a given bilinear form
  virtual void assemble(
      FullMatrix<LO, SC>& matrix
      ) const = 0;

  virtual void assembleACABlock(
      BEBlockCluster<LO, SC>* block,
      FullMatrix<LO, SC>& U,
      FullMatrix<LO, SC>& V,
      void * voidIntegrator
      ) const;

  //! methods assembles an ACA matrix
  virtual void assemble(
      ACAMatrix<LO, SC>& matrix
      ) const;

  //! methods assembles a distributed ACA matrix
  virtual void assemble(
      MPIACAMatrix<LO, SC>& matrix
      ) const;

  inline LO getRowCount( ) {
    return this->counterRow;
  }

  inline LO getColCount( ) {
    return this->counterCol;
  }

protected:

  //! method assembles idx-th row of a local block
  virtual void assembleRow(
      const BEBlockCluster< LO, SC > & block,
      LO idx,
      const std::vector< LO > & cols,
      Vector<LO, SC>& row,
      void * voidIntegrator
      ) const = 0;

  //! method assembles idx-th row of a local block
  virtual void assembleColumn(
      const BEBlockCluster< LO, SC > & block,
      LO idx,
      const std::vector< LO > & rows,
      Vector<LO, SC>& col,
      void * voidIntegrator
      ) const = 0;

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
      ) const = 0;

  /*!
   * The method assembles matrix transforming linear basis functions to 
   * linear discontinuous basis functions.
   */
  void assembleP12p1disMat(
      ACAMatrix<LO, SC> & matrix
      ) const;

  /*!
   * The method assembles matrix transforming linear basis functions to 
   * linear discontinuous basis functions.
   */
  void assembleP1dis2p1Mat(
      ACAMatrix<LO, SC> & matrix
      ) const;

  virtual void * createIntegrator( ) const = 0;

  virtual void destroyIntegrator(
      void * voidIntegrator
      ) const = 0;

  void leavesInnerP12p1dis( ) const;

  void leavesInnerP1dis2p1( ) const;

  void leavesOuterP12p1dis( ) const;

  void leavesOuterP1dis2p1( ) const;

  //! boundary element space for full BEM
  BESpace<LO, SC>* space;

  //! quadrature rule
  int* quadratureOrder;

  //! quadrature type
  quadratureType quadrature;

  //! disjoint elements Gaussian quadrature order
  int* quadratureOrderDisjointElems;

  mutable LO counterRow;

  mutable LO counterCol;

};

}
// include .cpp file to overcome linking problems due to templates
#include "BEBilinearForm.cpp"


#endif	/* BEBILINEARFORM_H */

