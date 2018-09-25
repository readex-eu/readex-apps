/*!
 * @file    BEBilinearFormLaplace2Layer.h
 * @author  Michal Merta
 * @date    July 18, 2013
 * @brief   Header file for class BEBilinearFormLaplace2Layer
 *
 */

#ifndef BEBILINEARFORMLAPLACE2LAYER_H
#define	BEBILINEARFORMLAPLACE2LAYER_H

#include "BEBilinearForm.h"
#include "FMMMatrix.h"
#include "ACAMatrix.h"
#include "FMMKernelLaplace2Layer.h"
#include "Quadratures.h"

namespace bem4i {

/*!
 * Class representing the bilinear form for the Laplace double layer operator
 *
 * Provides methods for system matrix assembly.
 *
 */
template<class LO, class SC>
class BEBilinearFormLaplace2Layer : public BEBilinearForm<LO, SC> {
  // to get inner type of complex numbers (for Helmholtz)
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:
  using BEBilinearForm<LO, SC>::assemble;

  //! default constructor
  BEBilinearFormLaplace2Layer( );

  //! copy constructor
  BEBilinearFormLaplace2Layer(
      const BEBilinearFormLaplace2Layer& orig
      );

  /*!
   * Constructor taking the boundary element space on which the form acts
   */
  BEBilinearFormLaplace2Layer(
      BESpace<LO, SC>* space,
      int* quadratureOrder = nullptr,
      quadratureType quadrature = SauterSchwab,
      int* quadratureOrderDisjointElems = nullptr
      );

  //! destructor
  virtual ~BEBilinearFormLaplace2Layer( );

  /*!
   * The method assembles the Galerkin matrix for the Laplace double layer operator
   */
  virtual void assemble(
      FullMatrix<LO, SC>& matrix
      ) const;

  /*!
   * The method assembles the Galerkin matrix for the Laplace double layer
   * operator with p0p1 elements
   */
  void assembleP0P1(
      FullMatrix<LO, SC> & matrix
      ) const;

  void assembleP1P1(
      FullMatrix<LO, SC> & matrix
      ) const;

  void assembleP0P0(
      FullMatrix<LO, SC> & matrix
      ) const;

  /*!
   * Fast BEM - The method assembles the sparsified
   * Galerkin matrix for the Laplace single layer operator
   *
   * The Fast Multipole Method is used for the sparsification
   */
  virtual void assemble(
      FMMMatrix<LO, SC>& matrix
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

  void assembleRowP0P1Dis(
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

  void assembleColumnP0P1Dis(
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
      const BECluster< LO, SC > * leftCluster,
      const BECluster< LO, SC > * rightCluster,
      FullMatrix< LO, SC > & matrix,
      void * voidIntegrator
      ) const;

  void assembleP0P0(
      const BECluster< LO, SC > * leftCluster,
      const BECluster< LO, SC > * rightCluster,
      FullMatrix< LO, SC > & matrix,
      void * voidIntegrator
      ) const;

  void assembleP0P1Dis(
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

  /*!
   * The method assembles the Galerkin matrix for the Laplace single layer
   * operator using Intel Xeon Phi accelerator
   */
  void assembleMIC(
      FullMatrix<LO, SC> & matrix
      ) const {

    if ( this->space->getAnsatzFunctionType( ) == p1 &&
        this->space->getTestFunctionType( ) == p0 ) {
      this->assembleP0P1MIC( matrix );
    } else {
      std::cout << "Not implemented!" << std::endl;
    }
  };

  /*!
   * The method assembles the Galerkin matrix for the Laplace double layer
   * operator with p0p1 elements using Intel Xeon Phi accelerator
   */
  void assembleP0P1MIC(
      FullMatrix<LO, SC> & matrix
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

  //! default quadrature rule for this class (specific values set in .cpp file)
  static int defaultQuadratureOrder[2];

};

}

// include .cpp file to overcome linking problems due to templates
#include "BEBilinearFormLaplace2Layer.cpp"

#endif	/* BEBILINEARFORMLAPLACE2LAYER_H */
