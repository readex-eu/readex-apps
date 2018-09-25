/*!
 * @file    BEBilinearFormWaveHypersingular.h
 * @author  Michal Merta 
 * @date    December 19, 2013
 * @brief   Header file for class BEBilinearFormWaveHypersingular
 * 
 */

#ifndef BEBILINEARFORMWAVEHYPERSINGULAR_H
#define	BEBILINEARFORMWAVEHYPERSINGULAR_H

#include "BEBilinearForm.h"
#include "BESpaceTime.h"
#include "BEIntegratorWave.h"
#include "SparseMatrix.h"
#include "MPIBlockMatrix.h"
#include "Eigen/Dense" 

namespace bem4i {

/*! 
 * Class representing the bilinear form for the single layer operator for wave equation
 * 
 * Provides methods for system matrix assembly.
 * 
 */
template <class LO, class SC>
class BEBilinearFormWaveHypersingular : public BEBilinearForm<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:
  using BEBilinearForm<LO, SC>::assemble;

  //! default constructor
  BEBilinearFormWaveHypersingular( );

  //! copy constructor
  BEBilinearFormWaveHypersingular(
      const BEBilinearFormWaveHypersingular& orig
      );

  /*! 
   * Constructor taking the boundary element space on which the form acts
   */
  BEBilinearFormWaveHypersingular(
      BESpace<LO, SC>* space,
      int* quadratureOrder = nullptr,
      int timeQuadratureOrder = defaultTimeQuadrature,
      int* quadratureOrderDisjointElems = nullptr
      );

#ifdef EXPERIMENTAL_WAVE

  /*!
   * constructor taking number of pre and post basis functions
   */
  BEBilinearFormWaveHypersingular(
      BESpace<LO, SC>* space,
      int nPre, int nPos,
      int* quadratureOrder = nullptr,
      int timeQuadratureOrder = defaultTimeQuadrature
      ) {
    this->space = space;

    if ( quadratureOrder ) {
      this->quadratureOrder = quadratureOrder;
    } else {
      this->quadratureOrder = defaultQuadraturesSauterSchwab;
    }

    this->timeQuadOrder = timeQuadratureOrder;
    this->nPre = nPre;
    this->nPos = nPos;
  };
#endif

  //! destructor
  virtual ~BEBilinearFormWaveHypersingular( );

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
#ifndef EXPERIMENTAL_WAVE  
  virtual void assemble(
      BlockMatrix<LO, SC>& matrix
      ) const;
#endif
  /*! Method assembles a distributed block matrix of sparse matrices of a given bilinear form
   * 
   * matrix must have M*M initiated blocks of appropriate sizes
   * where M is the number of time-steps
   * 
   */
#ifndef EXPERIMENTAL_WAVE
  virtual void assemble(
      MPIBlockMatrix<LO, SC>& matrix
      ) const;
#endif

private:

  //! methods assembles triplet list (indexed globally in system matrix) corresponding to the n-th time basis and m-th time test function
  void assembleBlock(
      LO row,
      LO column,
      std::vector<Eigen::Triplet<SC, LO> > &tripletList,
      LO &nnz
      ) const;

  //! methods assembles sparse matrix corresponding to the n-th time basis and m-th time test function
  void assembleBlock(
      LO row,
      LO column,
      SparseMatrix<LO, SC> &block
      ) const;

  //! methods assembles sparse matrix corresponding to the n-th time basis and m-th time test function
  void assembleMPIBlock(
      LO row,
      LO column,
      SparseMatrix<LO, SC> &block,
      int owner ) const;

  /*! Method estimates fill of the matrix blocks
   * 
   * @param[in,out]   fill user preallocated array of estimated fill percentage
   * 
   */
  void estimateNnz( LO *nnz ) const;

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

  /*!
   * Sends a matrix block from one process to another
   */
  void sendBlock(
      MPIBlockMatrix<LO, SC>& matrix,
      int sourceRow,
      int sourceColumn,
      int destColumn,
      int destRow,
      int source,
      int destination ) const;
  
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

  int nPre, nPos;
};

}

// include .cpp file to overcome linking problems due to templates
#include "BEBilinearFormWaveHypersingular.cpp"

#endif	/* BEBILINEARFORMWAVEhYPERSINGULAR_H */

