/*!
 * @file    STFPreconditioner.h
 * @author  Michal Merta, Jan Zapletal 
 * @date    June 3, 2015
 * @brief   Header file for abstract class for preconditioner for STF transm.
 * 
 */

#ifndef STFPRECONDITIONER_H
#define	STFPRECONDITIONER_H

#include "LeftPreconditioner.h"
#include "BlockMatrix.h"
#include "SurfaceMesh3D.h"

namespace bem4i {

template<class LO, class SC>
class STFPreconditioner : public LeftPreconditioner<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;
  
public:
    /*!
   * @brief Constructor
   *
   * @param Vi     user supplied interior Helmholtz single-layer operator
   * @param Ve     user supplied exterior Helmholtz single-layer operator
   * @param Ki     user supplied interior Helmholtz double-layer operator
   * @param Ke     user supplied exterior Helmholtz double-layer operator
   * @param Di     user supplied interior Helmholtz hypersingular operator
   * @param De     user supplied interior Helmholtz hypersingular operator
   */
  STFPreconditioner(
      SurfaceMesh3D<LO, SC> * mesh, 
      LinearOperator< LO, SC > * Vi,
      LinearOperator< LO, SC > * Ve,
      LinearOperator< LO, SC > * Ki,
      LinearOperator< LO, SC > * Ke,
      LinearOperator< LO, SC > * Di,
      LinearOperator< LO, SC > * De,
      SparseMatrix< LO, SC > * M,
      bool p0p0 = false
      );

  /*!
   * @brief Applies operator on a vector ( y = beta * y + alpha * this * x )
   *
   * @param x      operator argument
   * @param y      result (user pre-allocated)
   * @param alpha  y = beta * y + alpha * this * x
   * @param beta   y = beta * y + alpha * this * x
   */
  void apply(
      Vector<LO, SC> const & x,
      Vector<LO, SC> & y,
      bool transA = false,
      SC alpha = 1.0,
      SC beta = 0.0
      );
  
  virtual ~STFPreconditioner( ) {};
private:
  
  //! whether we apply P0P0 approximation
  bool p0p0;
  
    //! Surface mesh
  SurfaceMesh3D< LO, SC > * mesh;
  
   //! interior single-layer operator
  LinearOperator< LO, SC > * Vi;

  //! exterior single-layer operator
  LinearOperator< LO, SC > * Ve;

  //! interior double-layer operator
  LinearOperator< LO, SC > * Ki;

  //! exterior double-layer operator
  LinearOperator< LO, SC > * Ke;

  //! interior hypersingular operator
  LinearOperator< LO, SC > * Di;

  //! exterior hypersingular operator
  LinearOperator< LO, SC > * De;
  
  //! identity operator
  SparseMatrix< LO, SC > * M;
  
  STFPreconditioner( );
  STFPreconditioner( const STFPreconditioner& orig );

};

} // namespace bem4i

// include .cpp file to overcome linking problems due to templates
#include "STFPreconditioner.cpp"

#endif	/* STFPRECONDITIONER_H */

