/*!
 * @file    BEBilinearFormLaplace1Layer.h
 * @author  Michal Merta 
 * @date    July 12, 2013
 * @brief   Header file for class BEBilinearFormLaplace1Layer
 * 
 */

#ifndef BEBILINEARFORMLAPLACE1LAYER_H
#define	BEBILINEARFORMLAPLACE1LAYER_H

#include "BEBilinearForm.h"
#include "FMMMatrix.h"
#include "ACAMatrix.h"
#include "FMMKernelLaplace1Layer.h"

namespace bem4i {

/*! 
 * Class representing the bilinear form for the Laplace single layer operator
 * 
 * Provides methods for system matrix assembly.
 * 
 */
template<class LO, class SC>
class BEBilinearFormLaplace1Layer : public BEBilinearForm<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;



public:
  using BEBilinearForm<LO, SC>::assemble;

  //! default constructor
  BEBilinearFormLaplace1Layer( );

  //! copy constructor
  BEBilinearFormLaplace1Layer(
      const BEBilinearFormLaplace1Layer& orig
      );

  /*! 
   * Constructor taking the boundary element space on which the form acts
   */
  BEBilinearFormLaplace1Layer(
      BESpace<LO, SC>* space,
      int* quadratureOrder = nullptr,
      quadratureType quadrature = SauterSchwab,
      int* quadratureOrderDisjointElems = nullptr,
      bool employSymmetricity = false
      );

  //! destructor
  virtual ~BEBilinearFormLaplace1Layer( );

  /*!
   * The method assembles the Galerkin matrix for the Laplace single layer operator
   */
  virtual void assemble(
      FullMatrix<LO, SC> & matrix
      ) const;

  //! methods assembles a distributed ACA matrix
  virtual void assemble(
      MPIBlockACAMatrix<LO, SC>& matrix
      ) const;

  /*!
   * The method assembles the Galerkin matrix for the Laplace single layer operator
   */
  virtual void assembleP0P0(
      FullMatrix<LO, SC> & matrix
      ) const;

  virtual void assembleP1P1(
      FullMatrix<LO, SC> & matrix
      ) const;

  /*!
   * The method assembles the Galerkin matrix for the Laplace single layer 
   * operator employing the symmetry of the matrix
   */
  virtual void assembleSymmetric(
      FullMatrix<LO, SC>& matrix
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

  //! methods assembles an ACA matrix

  /*virtual void assemble(
      ACAMatrix<LO, SC>& matrix
      ) const;
   */
  void assembleWith2Layer(
      FullMatrix<LO, SC> & V,
      FullMatrix<LO, SC> & K,
      BESpace< LO, SC > & bespaceK
      ) const {

#if N_MIC > 0
    if ( this->space->getAnsatzFunctionType( ) == p0 &&
        this->space->getTestFunctionType( ) == p0 &&
        bespaceK.getAnsatzFunctionType( ) == p1 &&
        bespaceK.getTestFunctionType( ) == p0
        ) {
      this->assembleWith2LayerP0P0P0P1MIC( V, K, bespaceK );
    } else {
      std::cout << "Not implemented!" << std::endl;
    }

#else
    std::cout << "Not implemented!" << std::endl;
#endif
  }


protected:
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

  //! use symmetricity of the matrix to speed up computation
  bool employSymmetricity;

  /*!
   * The method assembles the Galerkin matrix for the Laplace single layer 
   * operator using Intel Xeon Phi accelerator
   */
  void assembleMIC(
      FullMatrix<LO, SC> & matrix
      ) const {

    if ( this->space->getAnsatzFunctionType( ) == p0 &&
        this->space->getTestFunctionType( ) == p0 ) {
      this->assembleP0P0MIC( matrix );
    } else {
      std::cout << "Not implemented!" << std::endl;
    }
  };

  /*!
   * The method assembles the Galerkin matrix for the Laplace single layer 
   * operator with p0p0 elements using Intel Xeon Phi accelerator
   */
  void assembleP0P0MIC(
      FullMatrix<LO, SC> & matrix
      ) const;

  /*!
   * The method assembles the Galerkin matrix for the Laplace single and double
   * layer operator with p0p0 and p1p0 elements using Intel Xeon Phi accelerator
   */
  void assembleWith2LayerP0P0P0P1MIC(
      FullMatrix<LO, SC> & V,
      FullMatrix<LO, SC> & K,
      BESpace< LO, SC > & bespaceK
      ) const;

  //! methods assembles an ACA matrix
  void assembleMIC(
      ACAMatrix<LO, SC>& matrix
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
#include "BEBilinearFormLaplace1Layer.cpp"

#endif	/* BEBILINEARFORMLAPLACE1LAYER_H */

