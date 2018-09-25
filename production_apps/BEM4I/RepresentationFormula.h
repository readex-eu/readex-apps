/*!
 * @file    RepresentationFormula.h
 * @author  Michal Merta 
 * @date    November 2, 2013
 * @brief   Header file for pure virtual class RepresentationFormula
 * 
 */

#ifndef REPRESENTATIONFORMULA_H
#define	REPRESENTATIONFORMULA_H

#include "Vector.h"
#include "BESpace.h"
#include "Macros.h"
#include "Mesh.h"

namespace bem4i {

/*! 
 * Abstract class for representation formula
 * 
 */
template<class LO, class SC>
class RepresentationFormula {

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! default constructor
  RepresentationFormula( );

  //! copy constructor
  RepresentationFormula(
      const RepresentationFormula& orig
      );

  //! destructor
  virtual ~RepresentationFormula( );


  /*!
   * function evaluates formula in a given point x
   *  
   * @param[in]   x
   */
  virtual SC evaluate(
      const SCVT *x
      ) const = 0;

  /*!
   * function evaluates formula in given array of points
   *  
   * @param[in]       x
   * @param[in,out]   result
   */
  virtual void evaluate(
      const SCVT *x,
      LO n,
      bool interior,
      Vector<LO, SC> & values
      ) const = 0;

  /*!
   * function evaluates formula in nodes of given mesh
   *  
   * @param[in]       mesh input mesh
   * @param[in]       interior interior/exterior flag
   * @param[in,out]   values
   */
  virtual void evaluate(
      Mesh<LO, SC> &mesh,
      bool interior,
      Vector<LO, SC> & values
      ) const = 0;

protected:

  //! boundary element space
  BESpace<LO, SC>* space;

  //! Dirichlet data
  Vector<LO, SC> *dirichlet;

  //! Neumann data
  Vector<LO, SC> *neumann;



};

} // end of namespace bem4i

// include .cpp file to overcome linking problems due to templates
#include "RepresentationFormula.cpp"

#endif	/* REPRESENTATIONFORMULA_H */

