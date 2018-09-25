/*!
 * @file    STFOperator.h
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    June 26, 2015
 * @brief   Header file for class STFOperator
 */

#ifndef STFOPERATOR_H
#define	STFOPERATOR_H

#include <vector>

#include "LinearOperator.h"
#include "Vector.h"
#include "BESpace.h"
#include "SparseMatrix.h"
#include "IterativeSolver.h"

namespace bem4i {

/*!
 * Class representing a weak operator on the given discretized space
 */
template<class LO, class SC>
class STFOperator :
public LinearOperator<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! todo: assembly without V, K, D, M
  /*!
   * @brief Constructor with SurfaceMesh3D
   *
   * @param mesh   surface mesh
   * @param ki     interior wave number
   * @param ke     exterior wave number
   * @param Vi     user supplied interior Helmholtz single-layer operator
   * @param Ve     user supplied exterior Helmholtz single-layer operator
   * @param Ki     user supplied interior Helmholtz double-layer operator
   * @param Ke     user supplied exterior Helmholtz double-layer operator
   * @param Di     user supplied interior Helmholtz hypersingular operator
   * @param De     user supplied interior Helmholtz hypersingular operator
   * @param M      user supplied identity operator
   */
  STFOperator(
      SurfaceMesh3D< LO, SC > * mesh,
      SC ki,
      SC ke,
      LinearOperator< LO, SC > * Vi,
      LinearOperator< LO, SC > * Ve,
      LinearOperator< LO, SC > * Ki,
      LinearOperator< LO, SC > * Ke,
      LinearOperator< LO, SC > * Di,
      LinearOperator< LO, SC > * De,
      LinearOperator< LO, SC > * M,
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

  bool LUSolve(
      Vector<LO, SC> & rhs
      );

  void getRHS(
      const Vector< LO, SC > & g,
      const Vector< LO, SC > & f,
      Vector< LO, SC > & rhs
      ) const;

  void getTraces(
      const Vector< LO, SC > & res,
      const Vector< LO, SC > & f,
      const Vector< LO, SC > & g,
      Vector< LO, SC > & ui,
      Vector< LO, SC > & ue,
      Vector< LO, SC > & ti,
      Vector< LO, SC > & te
      ) const;

  //! Destructor
  virtual ~STFOperator( );


private:
  
  //! whether we apply P0P0 approximation
  bool p0p0;

  //! Surface mesh
  SurfaceMesh3D< LO, SC > * mesh;

  //! interior wave number
  SC ki;

  //! exterior wave number
  SC ke;

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

  //! Identity operator
  LinearOperator< LO, SC > * M;

  //! Default constructor (empty)
  STFOperator( );

  //! Copy constructor (empty)
  STFOperator(
      const STFOperator & orig
      );

};

}
// include .cpp file to overcome linking problems due to templates
#include "STFOperator.cpp"

#endif	/* STFOPERATOR_H */
