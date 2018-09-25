/*!
 * @file    VolumeMesh3D.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    November 6, 2013
 * @brief   Header file for class VolumeMesh3D
 * 
 */

#ifndef VOLUMEMESH3D_H
#define	VOLUMEMESH3D_H

#include <iostream>
#include <fstream>
#include <set>
#include "Mesh.h"
#include "Tree.h"
#include "Quadratures.h"
#include "Macros.h"

namespace bem4i {

/*! 
 * Class representing a volume mesh of a 3D body 
 * 
 */
template<class LO, class SC>
class VolumeMesh3D : public Mesh<LO, SC> {

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:
  //! default constructor
  VolumeMesh3D( );

  //! copy constructor
  VolumeMesh3D( const VolumeMesh3D& orig );

  //! destructor
  virtual ~VolumeMesh3D( );

  /*!
   * Loads mesh from a default file format
   * 
   * @param   meshFile a string with a path to the mesh file
   */
  virtual void load( const string& meshFile, SCVT scaleFactor = 1.0 );

  /*!
   * Prints mesh info to stdout
   */
  virtual void printInfo( );

  /*!
   * Prints mesh to the xml paraview file format
   * 
   * @param[in]   meshFile string with the target file name
   */
  virtual void printParaviewVtu( const string& meshFile );

  /*!
   * Prints mesh to the xml paraview file format
   * 
   * @param[in]   meshFile string with the target file name
   */
  virtual void printParaviewVtu( const string& meshFile, int nNodal, string* nodeNames, Vector< LO, SC >** nodalData, int nElem, string* elemNames, Vector< LO, SC >** elemData );

  virtual void loadFromNetgen( const string& meshFile ) {
  };

  /*!
   * Loads mesh from a paraview file format
   * 
   * @param   meshFile a string with a path to the mesh file
   */
  virtual void loadFromParaview( const string& meshFile ) {
  };

  /*!
   * Prints mesh to the legacy paraview file format
   * 
   * @param[in]   meshFile string with the target file name
   */
  virtual void printParaviewVtk( const string& meshFile ) {
  };

  /*!
   * computes areas of all elements
   */
  virtual void initArea( ) {
  };
protected:

private:

};

} // end of namespace bem4i

// include .cpp file to overcome linking problems due to templates
#include "VolumeMesh3D.cpp"

#endif	/* VOLUMEMESH3D_H */

