/*!
 * @file    BEBilinearFormWave1Layer.h
 * @author  Michal Merta 
 * @date    December 10, 2013
 * @brief   Header file for class BEBilinearFormWave1Layer
 * 
 */

#ifndef BEBILINEARFORMWAVE1LAYER_H
#define	BEBILINEARFORMWAVE1LAYER_H

#include "BEBilinearForm.h"
#include "BESpaceTime.h"
#include "BEIntegratorWave.h"
#include "SparseMatrix.h"

namespace bem4i {

/*! 
 * Class representing the bilinear form for the single layer operator for wave equation
 * 
 * Provides methods for system matrix assembly.
 * 
 */
template <class LO, class SC>
class BEBilinearFormWave1Layer : public BEBilinearForm<LO, SC> {
public:

  using BEBilinearForm<LO, SC>::assemble;

  //! default constructor
  BEBilinearFormWave1Layer( );

  //! copy constructor
  BEBilinearFormWave1Layer(
      const BEBilinearFormWave1Layer& orig
      );

  /*! 
   * Constructor taking the boundary element space on which the form acts
   */
  BEBilinearFormWave1Layer(
      BESpace<LO, SC>* space,
      int* quadratureOrder = nullptr,
      int timeQuadratureOrder = defaultTimeQuadrature,
      int* quadratureOrderDisjointElems = nullptr
      );

  //! destructor
  virtual ~BEBilinearFormWave1Layer( );

  //! method assembles a full matrix of a given bilinear form
  virtual void assemble(
      FullMatrix<LO, SC>& matrix
      ) const;

  //! method assembles a sparse matrix of a given bilinear form
  virtual void assemble(
      SparseMatrix<LO, SC>& matrix
      ) const;

  /*! Method assembles a block matrix of sparse matricesof a given bilinear form
   * 
   * matrix must have M*M initiated blocks of appropriate sizes
   * where M is the number of time-steps
   * 
   */
  virtual void assemble(
      BlockMatrix<LO, SC>& matrix
      ) const;

  //! methods assembles sparse matrix corresponding to the n-th time basis and m-th time test function
  void assembleBlock(
      LO row,
      LO column,
      SparseMatrix<LO, SC> &block
      ) const;

private:

protected:

  //! method assembles idx-th row of a local block

  virtual void assembleRow(
      const BEBlockCluster< LO, SC > & block,
      LO idx,
      const std::vector< LO > & cols,
      Vector<LO, SC>& row,
      void * voidIntegrator
      ) const {
  }

  //! method assembles idx-th row of a local block

  virtual void assembleColumn(
      const BEBlockCluster< LO, SC > & block,
      LO idx,
      const std::vector< LO > & rows,
      Vector<LO, SC>& col,
      void * voidIntegrator
      ) const {
  }

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
      ) const {
  }

  virtual void * createIntegrator( ) const {
    return static_cast<void *> ( new BEIntegratorWave< LO, SC >(
        this->space, this->quadratureOrder, this->timeQuadOrder ) );
  }
  
  virtual void destroyIntegrator(
      void * voidIntegrator
      ) const {
    delete static_cast<BEIntegratorWave< LO, SC > *> ( voidIntegrator );
  }

  int timeQuadOrder;
};

}

// include .cpp file to overcome linking problems due to templates
#include "BEBilinearFormWave1Layer.cpp"

#endif	/* BEBILINEARFORMWAVE1LAYER_H */

