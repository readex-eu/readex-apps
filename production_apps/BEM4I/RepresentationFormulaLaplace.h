/*!
 * @file    RepresentationFormula.h
 * @author  Michal Merta 
 * @date    November 5, 2013
 * @brief   Header file for class RepresentationFormulaLaplace
 * 
 */

#ifndef REPRESENTATIONFORMULALAPLACE_H
#define	REPRESENTATIONFORMULALAPLACE_H

#include "RepresentationFormula.h"

namespace bem4i {

/*! 
 * Class for representation formula for Laplace operator
 * 
 */
template<class LO, class SC>
class RepresentationFormulaLaplace : public RepresentationFormula<LO, SC> {

  typedef typename GetType<LO, SC>::SCVT SCVT;
public:

  //! default constructor
  RepresentationFormulaLaplace( );

  /*!
   * constructor taking Dirichlet and Neumann data as arguments
   * 
   * @param[in]     space
   * @param[in]     dirichlet - vector of dirichlet data
   * @param[in]     neumann - vector of neumann data
   * 
   */
  RepresentationFormulaLaplace(
      BESpace<LO, SC>* space,
      Vector<LO, SCVT> *dirichlet,
      Vector<LO, SCVT>* neumann,
      int quadOrder = 4
      );

  //! copy constructor
  RepresentationFormulaLaplace(
      const RepresentationFormulaLaplace& orig
      );

  // destructor
  virtual ~RepresentationFormulaLaplace( );


  /*!
   * function evaluates formula in a given point x
   *  
   * @param[in]   x
   */
  virtual SC evaluate(
      const SCVT *x
      ) const;

  /*!
   * function evaluates formula in given array of points
   *  
   * @param[in]       x input array of point coordinates
   * @param[in]       n number of points
   * @param[in,out]   result
   * @param[in]       interior flag interior/exterior
   */
  virtual void evaluate(
      const SCVT *x,
      LO n,
      bool interior,
      Vector<LO, SC> & values
      ) const;

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
      ) const;

  /*!
   * functions saves points to paraview point cloud format
   */
  void printParaviewVtu(
      const string& pointFile,
      SC* points,
      LO nPoints,
      int nNodal,
      string* nodeNames,
      Vector< LO, SC >** nodalData
      ) const;

private:
  
  int quadOrder;

};

} // end of namespace bem4i

// include .cpp file to overcome linking problems due to templates
#include "RepresentationFormulaLaplace.cpp"

#endif	/* REPRESENTATIONFORMULALAPLACE_H */

