/*!
 * @file    Mesh.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    July 4, 2013
 * @brief   Header file for abstract class Mesh
 * 
 */

#ifndef MESH_H
#define	MESH_H

#include <string>
#include <vector>
#include <set>
#include <cmath>
#include "Vector.h"
#include "Macros.h"

using std::string;
using std::vector;

namespace bem4i {

#ifdef OPTIMIZATION
template< class LO, class SC >
class OpenMeshWrapper;
#endif

//! struct of multipole coefficients of cluster - for fast BEM

template<class LO, class SC>
struct MultipoleMoments {
  std::complex<SC>* momentsSL; // multipole moments of a cell - filled by the 
  // upward method of appropriate kernel 
  std::complex<SC>* momentsDL; // moments for double layer potential

  std::complex<SC>* momTilde;
  std::complex<SC>* momTildeDL;

  std::complex<SC>* localExpCoefSL;
  std::complex<SC>* localExpCoefDL;

  ~MultipoleMoments( ) {
    if ( momentsSL )
      delete [] momentsSL;
    if ( momentsDL )
      delete [] momentsDL;
    if ( momTilde )
      delete [] momTilde;
    if ( momTildeDL )
      delete [] momTildeDL;
    if ( localExpCoefSL )
      delete [] localExpCoefSL;
    if ( localExpCoefDL )
      delete [] localExpCoefDL;
  }
};

//! cluster of boundary elements (for creation of cluster tree)  

template<class LO, class SC>
struct BECluster {
  typedef typename GetType<LO, SC>::SCVT SCVT;

  //! coordinates of the centroid
  SCVT *centroid;

  //! radius of the cluster
  SCVT radius;

  //! mass (length/area) of the cluster
  SCVT mass;

  //! number of elements                          
  LO nelems;

  //! number of nodes
  LO nnodes;

  //! element indices
  std::vector<LO>* elems;

  //! global nodes indices
  std::set<LO>* nodes;

  // The following data are used in fast multipole method

  MultipoleMoments<LO, SC>* moments;

  BECluster<LO, SC> *parent;


  std::vector<BECluster<LO, SC>* > *admClusters;

  // M2L and L2L coefficients
  std::complex<SCVT> **M2L;
  std::complex<SCVT> **L2L;

  BECluster( ) {
    centroid = nullptr;
    elems = nullptr;
    moments = nullptr;
    admClusters = nullptr;
    M2L = nullptr;
    L2L = nullptr;
    nodes = nullptr;
  }

  ~BECluster( ) {
    if ( centroid ) {
      delete [] centroid;
    }
    if ( elems ) {
      delete elems;
    }
    if ( nodes ) {
      delete nodes;
    }
    if ( moments )
      delete moments;
    if ( admClusters )
      delete admClusters;
  }
};

/*! 
 * Abstract class representing an arbitrary mesh in 2D or 3D 
 * 
 */
template<class LO, class SC>
class Mesh {
  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! default constructor
  Mesh( );

  //! copy constructor
  Mesh( const Mesh& orig );

  //! destructor
  virtual ~Mesh( );

  //! returns a dimension of a mesh

  inline int getDim( ) const {
    return dim;
  }

  //! returns number of nodes

  inline LO getNNodes( ) const {
    return nNodes;
  }

  //! returns number of edges

  inline LO getNEdges( ) const {
    return nEdges;
  }

  //! returns number of elements

  inline LO getNElements( ) const {
    return nElems;
  }

  /*!
   * Returns element of a mesh
   * 
   * @param   idx index of element to return
   * @param   element pointer to <em>preallocated</em> array of LOs where the element indices will be copied
   */
  inline void getElement(
      LO idx,
      LO * element
      ) const {

    memcpy( element, &elems[idx * nNodesPerElem], sizeof (LO ) * nNodesPerElem );
  }

  /*!
   * Returns node of a mesh
   * 
   * @param   idx index of node to return
   * @param   node pointer to <em>preallocated</em> array of doubles where the node coordinates will be copied
   */
  inline void getNode( LO idx, SCVT* node ) const {
    memcpy( node, &nodes[idx * dim], sizeof (SCVT ) * dim );
  }

  /*!
   * Returns three nodes of a triangular mesh associated with a given element
   * 
   * @param[in]       elemIdx   index of element
   * @param[in,out]   x1        pointer to preallocated array for first node
   * @param[in,out]   x2        pointer to preallocated array for second node
   * @param[in,out]   x3        pointer to preallocated array for third node
   */
  inline void getNodes( LO elemIdx, SCVT * x1, SCVT * x2, SCVT * x3 ) const {
    memcpy( x1, &nodes[elems[elemIdx * 3] * dim], sizeof (SCVT ) * dim );
    memcpy( x2, &nodes[elems[elemIdx * 3 + 1] * dim], sizeof (SCVT ) * dim );
    memcpy( x3, &nodes[elems[elemIdx * 3 + 2] * dim], sizeof (SCVT ) * dim );
  }

  inline vector<SCVT>* getAreas( ) {
    return & this->area;
  }




  /*!
   * Loads mesh from a default file format
   * 
   * @param   meshFile a string with a path to the mesh file
   */
  virtual bool load( const string& meshFile, SCVT scaleFactor = 1.0 ) = 0;

  /*!
   * Loads mesh from a netgen file format
   * 
   * @param   meshFile a string with a path to the mesh file
   */
  virtual bool loadFromNetgen( const string& meshFile ) = 0;

  /*!
   * Loads mesh from a paraview file format
   * 
   * @param   meshFile a string with a path to the mesh file
   */
  virtual bool loadFromParaview( const string& meshFile ) = 0;

  /*!
   * Prints mesh info to stdout
   */
  virtual void printInfo( ) = 0;

  /*!
   * Prints mesh to the legacy paraview file format
   * 
   * @param[in]   meshFile string with the target file name
   */
  virtual void printParaviewVtk( const string& meshFile ) = 0;

  /*!
   * Prints mesh to the xml paraview file format
   * 
   * @param[in]   meshFile string with the target file name
   */
  virtual void printParaviewVtu( const string& meshFile ) = 0;

  /*!
   * computes areas of all elements
   */
  virtual void initArea( ) = 0;

  //! returns area of given element (must be precomputed by initArea())

  SCVT getElemArea( LO elemIdx ) const {
    return (SCVT) area[elemIdx];
  }

  //! returns number of nodes per one element

  inline int getNNodesPerElem( ) {
    return nNodesPerElem;
  }

  inline vector<LO>* getElements( ) {
    return & this->elems;
  }

  inline vector<SCVT>* getNodes( ) {
    return & this->nodes;
  }

protected:

  //! dimension of the mesh
  int dim;

  //! number of mesh nodes
  LO nNodes;

  //! number of mesh elements
  LO nElems;

  //! number of edges
  LO nEdges;

  //! number of nodes per one element
  int nNodesPerElem;

  //! 1D vector of node coordinates
  /*! Nodes stored as [x1_1, x1_2, x1_3, x2_1, x2_2, x2_3, ...]. */
  vector<SCVT> nodes;

  //! 1D vector of elements mapping to the vector of nodes
  /*! Elements stored as [node1_1, node1_2, node1_3, node2_1, node2_2, node2_3, ...]. */
  vector<LO> elems;

  //! indices of edges
  vector<LO> edges;

  //! vector of element areas
  vector<SCVT> area;

  //! total area of a mesh
  SCVT totalArea;

private:

};

}

// include .cpp file to overcome linking problems due to templates
#include "Mesh.cpp"

#endif	/* MESH_H */
