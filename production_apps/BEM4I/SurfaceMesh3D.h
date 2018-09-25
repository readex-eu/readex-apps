/*!
 * @file    SurfaceMesh3D.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    July 5, 2013
 * @brief   Header file for class SurfaceMesh3D
 * 
 */

#ifndef SURFACEMESH3D_H
#define	SURFACEMESH3D_H

#include <iostream>
#include <fstream>
#include <set>
#include <algorithm>
#include <iterator>
#include "Mesh.h"
#include "Tree.h"
#include "Quadratures.h"
#include "Macros.h"
#include "BESpace.h"
#include "PotentialsLaplace.h"
#include "PotentialsHelmholtz.h"
#include "BEIntegratorLaplace.h"
#include "SparseMatrix.h"

namespace bem4i {

#ifdef OPTIMIZATION
template< class LO, class SC >
class OpenMeshWrapper;
#endif

//template< class LO, class SC >
//class BESpace;

/*! 
 * Class representing a surface mesh of a 3D body 
 * 
 */
template<class LO, class SC>
class SurfaceMesh3D : public Mesh<LO, SC> {
  typedef typename GetType<LO, SC>::SCVT SCVT;

  friend class SurfaceMesh3D<LO, SCVT>;
  friend class SurfaceMesh3D<LO, std::complex<SCVT> >;

#ifdef OPTIMIZATION
  friend class OpenMeshWrapper< LO, SC >;
#endif

public:
  //! default constructor
  SurfaceMesh3D( );

  //! read from file
  SurfaceMesh3D(
      const string& meshFile
      );

  //! construct from vector
  SurfaceMesh3D(
      const std::vector< SCVT > & nodes,
      const std::vector< LO > & elems
      );

  //! copy constructor
  SurfaceMesh3D(
      const SurfaceMesh3D& orig
      );

  //! destructor
  virtual ~SurfaceMesh3D( );

  void copyToComplex(
      SurfaceMesh3D< LO, std::complex< SCVT > > & copy
      ) const;

  void copyToReal(
      SurfaceMesh3D< LO, SCVT > & copy
      ) const;

  /*!
   * Loads mesh from a default file format
   * 
   * @param   meshFile a string with a path to the mesh file
   */
  virtual bool load(
      const string& meshFile,
      SCVT scaleFactor = 1.0
      );

  /*!
   * Loads mesh from a default file format
   * 
   * @param   meshFile a string with a path to the mesh file
   */
  virtual bool loadSmf(
      const string& meshFile
      );

  virtual bool loadSmfOpenFTLOptimize(
      const string& meshFile,
      Vector< LO, SCVT > & fixity,
      Vector< LO, SCVT > & curvature
      );

  /*!
   * Loads mesh from a netgen file format
   * 
   * @param   meshFile a string with a path to the mesh file
   */
  virtual bool loadFromNetgen(
      const string& meshFile
      );

  /*!
   * Loads mesh from a paraview file format
   * 
   * @param   meshFile a string with a path to the mesh file
   */
  virtual bool loadFromParaview(
      const string& meshFile
      );

  /*!
   * Prints mesh info to stdout
   */
  virtual void printInfo( );

  /*!
   * Prints mesh to the legacy paraview file format
   * 
   * @param[in]   meshFile string with the target file name
   */
  virtual void printParaviewVtk(
      const string & meshFile
      );

  /*!
   * Prints mesh to the paraview file format
   * 
   * @param[in]   meshFile string with the target file name
   */
  virtual void printParaviewVtu(
      const string & meshFile
      ) {

    this->printParaviewVtu( meshFile, nullptr, nullptr, nullptr, nullptr );
  }

  /*!
   * Prints mesh to the input format
   * 
   * @param[in]   meshFile string with the target file name
   */
  virtual void print(
      const string & meshFile
      );

  /*!
   * Prints mesh to the smf format
   * 
   * @param[in]   meshFile string with the target file name
   */
  virtual void printSmf(
      const string & meshFile
      );

  //  /*!
  //   * Prints mesh to the xml paraview file format
  //   * 
  //   * @param[in]   meshFile string with the target file name
  //   */
  //  virtual void printParaviewVtu( const string& meshFile );

  /*!
   * Prints mesh to the xml paraview file format
   * 
   * @param[in]   meshFile string with the target file name
   */
  virtual void printParaviewVtu(
      const string& meshFile,
      int nNodal,
      string* nodeNames,
      Vector< LO, SC >** nodalData,
      int nElem,
      string* elemNames,
      Vector< LO, SC >** elemData
      );

  /*!
   * Prints mesh to the xml paraview file format
   * 
   * @param[in] meshFile string with the target file name
   */
  void printParaviewVtu(
      const string & meshFile,
      const std::vector< string > * nodeNames,
      const std::vector< Vector< LO, SC >* > * nodalData,
      const std::vector< string > * elemNames,
      const std::vector< Vector< LO, SC >* > * elemData,
      const std::vector< string > * nodeVNames = nullptr,
      const std::vector< Vector< LO, SC >* > * nodalVData = nullptr,
      const std::vector< string > * elemVNames = nullptr,
      const std::vector< Vector< LO, SC >* > * elemVData = nullptr
      ) const;

  /*!
   * Updates node position
   * 
   * @param[in] idx node index
   * @param[in] node new coordinates 
   */
  void setNode(
      LO idx,
      const SCVT * node
      ) {
    this->nodes[ 3 * idx ] = node[ 0 ];
    this->nodes[ 3 * idx + 1 ] = node[ 1 ];
    this->nodes[ 3 * idx + 2 ] = node[ 2 ];
  }

  //! computes area of one element

  inline SCVT computeElemArea(
      LO elemIdx
      ) const {

    const SCVT *x1 = &this->nodes[this->elems[elemIdx * 3] * this->dim];
    const SCVT *x2 = &this->nodes[this->elems[elemIdx * 3 + 1] * this->dim];
    const SCVT *x3 = &this->nodes[this->elems[elemIdx * 3 + 2] * this->dim];

    SCVT u[3], v[3];
    for ( int i = 0; i < 3; i++ ) {
      u[i] = x1[i] - x2[i];
      v[i] = x3[i] - x2[i];
    }

    SCVT s[3];
    s[0] = u[1] * v[2] - u[2] * v[1];
    s[1] = u[2] * v[0] - u[0] * v[2];
    s[2] = u[0] * v[1] - u[1] * v[0];
    return 0.5 * std::sqrt( s[0] * s[0] + s[1] * s[1] + s[2] * s[2] );
  }

  //! returns mesh discretization parameter

  inline SCVT getDiscretizationParameter( ) {
    return std::sqrt(
        *std::max_element( this->area.begin( ), this->area.end( ) ) );
  }

  //! returns the centroid of the given element

  inline void getCentroid(
      LO i,
      SCVT centroid[ 3 ]
      ) const {
    SCVT x1[ 3 ], x2[ 3 ], x3[ 3 ];
    this->getNodes( i, x1, x2, x3 );
    centroid[ 0 ] = 1.0 / 3.0 * ( x1[ 0 ] + x2[ 0 ] + x3[ 0 ] );
    centroid[ 1 ] = 1.0 / 3.0 * ( x1[ 1 ] + x2[ 1 ] + x3[ 1 ] );
    centroid[ 2 ] = 1.0 / 3.0 * ( x1[ 2 ] + x2[ 2 ] + x3[ 2 ] );
  }

  //! returns the centroid of the mesh

  inline void getCentroid(
      SCVT centroid[ 3 ]
      ) const {
    SCVT x[ 3 ];
    centroid[ 0 ] = centroid[ 1 ] = centroid[ 2 ] = 0.0;
    for ( LO i = 0; i < this->nNodes; ++i ) {
      this->getNode( i, x );
      centroid[ 0 ] += x[ 0 ];
      centroid[ 1 ] += x[ 1 ];
      centroid[ 2 ] += x[ 2 ];
    }
    centroid[ 0 ] /= this->nNodes;
    centroid[ 1 ] /= this->nNodes;
    centroid[ 2 ] /= this->nNodes;
  }

  bool isInside(
      const SCVT * x
      ) const;

  void isInside(
      const SCVT * x,
      LO nPoints,
      bool * ret
      ) const;

  void trimEvaluationGrid(
      SurfaceMesh3D< LO, SC > & grid,
      SurfaceMesh3D< LO, SC > & gridTrimmed,
      bool trimInterior = true
      ) const;

  void trimEvaluationGrid(
      SurfaceMesh3D< LO, SC > & grid,
      SurfaceMesh3D< LO, SC > & interior,
      SurfaceMesh3D< LO, SC > & exterior
      ) const;

  /*!
   * computes area of all elements
   */
  virtual void initArea( );

  /*!
   * Initializes mapping from nodes to adjacent elements
   */
  void initNode2Elems( );

  /*!
   * Get elems connected to node idx (need to initNode2Elems first!)
   * 
   * @param[in] idx index of the node
   * @param[in,out] elems vector of elems indices connected to the node idx
   */
  inline void getElements(
      LO idx,
      vector< LO > & elems
      ) const {
    if ( idx < this->node2elems.size( ) )
      elems = this->node2elems[ idx ];
  };

  inline void getElementGlobal(
      LO idx,
      LO * element
      ) const {

    this->getElement( idx, element );
    element[ 0 ] = this->l2gNodes[ element[ 0 ] ];
    element[ 1 ] = this->l2gNodes[ element[ 1 ] ];
    element[ 2 ] = this->l2gNodes[ element[ 2 ] ];
  }

  using Mesh< LO, SC >::getElements;

  /*!
   * Initializes mapping from nodes to adjacent elements
   */
  void initAdditiveCurvature( );

  /*!
   * Returns volume of closed surfaces
   */
  inline SCVT getVolume( ) const;

  /*!
   * Returns total surface area
   */
  inline SCVT getArea( ) const;

  /*!
   * Returns the curvature in node idx
   */
  inline SCVT getCurvature(
      LO idx
      ) const;

  /*!
   * Returns the curvature vector
   */
  inline void getCurvatureVector(
      Vector< LO, SC > & curvature
      ) const;

  /*!
   * Checks the size of curvature vector, returns true if equal to nNodes.
   */
  inline bool curvatureReady( ) const {
    return ( this->additiveCurvature.size( ) == this->nNodes ) ?
        true : false;
  }

  /*!
   * computes outer normal to each element
   */
  void initNormals( );

  /*!
   * changes orientation of surface elements
   */
  void flipNormals( );

  /*!
   * computes outer normal to each element
   */
  void initLocalCoordinates( );

  /*!
   * initializes default mapping from local to global indexing for elements
   */
  void initL2g( );

  /*!
   * Method sets the vector of local to global mapping of elements
   */
  void setL2gElems( std::vector<LO> & l2g ) {
    this->l2gElems = l2g;
  }

  /*!
   * Method sets the vector of local to global mapping of nodes
   */
  void setL2gNodes( std::vector<LO> & l2g ) {
    this->l2gNodes = l2g;
  }

  /*!
   * returns normal to idx-th element
   * 
   * [in]     idx     index of given element
   * [in,out] normal  user allocated 3D array of type SCVT
   */
  inline void getNormal(
      LO idx,
      SCVT* normal
      ) const {
    //std::cout << "Size n: " << this->normals.size() << std::endl;
    normal[0] = this->normals[3 * idx];
    normal[1] = this->normals[3 * idx + 1];
    normal[2] = this->normals[3 * idx + 2];
  }

  /*!
   * returns weighted normal to idx-th node
   * 
   * [in]     idx     index of given node
   * [in,out] normal  user allocated 3D array of type SCVT
   */
  inline void getNormalNodal(
      LO nodeIdx,
      SCVT* normal
      ) const {
    //TODO: need to init first. some flag?
    //this->initNode2Elems( );

    normal[ 0 ] = normal[ 1 ] = normal[ 2 ] = 0.0;
    SCVT area;
    SCVT totalArea = 0;
    LO elemIdx;
    SCVT elemNormal[ 3 ];

    for ( int i = 0; i < this->node2elems[ nodeIdx ].size( ); ++i ) {
      elemIdx = this->node2elems[ nodeIdx ][ i ];
      area = this->getElemArea( elemIdx );
      totalArea += area;
      this->getNormal( elemIdx, elemNormal );
      normal[ 0 ] += area * elemNormal[ 0 ];
      normal[ 1 ] += area * elemNormal[ 1 ];
      normal[ 2 ] += area * elemNormal[ 2 ];
    }

    SCVT norm = std::sqrt( normal[0] * normal[0] + normal[1] * normal[1] +
        normal[2] * normal[2] );
    normal[ 0 ] /= norm;
    normal[ 1 ] /= norm;
    normal[ 2 ] /= norm;
  }

  inline void getR1( LO idx, SCVT* r1, int rotate ) const {
    //std::cout << "Size r1: " << this->r1.size() << std::endl;
    r1[0] = this->r1[3 * ( 3 * idx + rotate )];
    r1[1] = this->r1[3 * ( 3 * idx + rotate ) + 1];
    r1[2] = this->r1[3 * ( 3 * idx + rotate ) + 2];
  }

  inline void getR2( LO idx, SCVT* r2, int rotate ) const {
    r2[0] = this->r2[3 * ( 3 * idx + rotate )];
    r2[1] = this->r2[3 * ( 3 * idx + rotate ) + 1];
    r2[2] = this->r2[3 * ( 3 * idx + rotate ) + 2];
  }

  inline void getAlpha1( LO idx, SCVT &alpha1, int rotate ) const {
    alpha1 = this->alpha1[ 3 * idx + rotate ];
  }

  inline void getAlpha2( LO idx, SCVT &alpha2, int rotate ) const {
    alpha2 = this->alpha2[ 3 * idx + rotate ];
  }

  inline void getSk( LO idx, SCVT &stau, int rotate ) const {
    stau = this->stau[ 3 * idx + rotate ];
  }

  /*!
   * fills given (user preallocated) array with indices of element edges
   * 
   * [in]     idx     index of given element
   * [in,out] edges   user allocated array for 3 LO indices of edges
   */
  inline void getEdges(
      LO idx,
      LO *edges
      ) {
    edges[0] = this->elem2edges[idx * 3];
    edges[1] = this->elem2edges[idx * 3 + 1];
    edges[2] = this->elem2edges[idx * 3 + 2];
  }

  inline void getEdge(
      LO idx,
      LO *edge
      ) {
    edge[0] = this->edges[2 * idx];
    edge[1] = this->edges[2 * idx + 1];
  }

  /*!
   * computes L2 norm for piecewise constant approximation
   */
  //SCVT l2RelativeErrorConst( Vector< LO, SC > const & anal, Vector< LO, SC > const & res ) const;

  /*!
   * computes L2 norm for piecewise linear approximation
   */
  //SCVT l2RelativeErrorLin( Vector< LO, SC > const & anal, Vector< LO, SC > const & res ) const;

  void getQuadratureNodes(
      SCVT const *x1,
      SCVT const *x2,
      SCVT const *x3,
      SCVT* nodes,
      int order = 5
      ) const;

  //! testing
  SCVT l2RelativeErrorConst(
      Vector< LO, SC > const & res,
      int order = 5,
      SC kappa = 2.0
      ) const;

  //! testing 
  SCVT l2RelativeErrorLin(
      Vector< LO, SC > const & res,
      int order = 5,
      SC kappa = 2.0
      ) const;

  /*!
   * Creates a tree of boundary element clusters using the nested dissection algorithm
   * 
   * @param[in,out] clusterTree           initialized tree of clusters
   * @param[in]     maxElemsPerCluseter   max. number of elements per one cluster
   */
  void nestedDissection(
      Tree<BECluster<LO, SC>*, LO > & clusterTree,
      LO maxElemsPerCluster
      );

  //! computes the ration between masses at each level, 
  //  - data has be a vector of n-level vectors!
  void getNestedProjectionMassRatio(
      Tree<BECluster<LO, SC>*, LO > & clusterTree,
      std::vector< std::vector< SC > > & data
      ) const;

  void getNestedProjectionMassRatio(
      Tree<BECluster<LO, SC>*, LO > & clusterTree,
      std::vector< Vector< LO, SC > * > & data
      ) const;

  void getNestedProjectionDiam(
      Tree<BECluster<LO, SC>*, LO > & clusterTree,
      std::vector< std::vector< SC > > & data
      ) const;


  //! Maps elements onto lower level of nestedDissection clusters
  void elemToNestedDissectionMapping(
      Tree<BECluster<LO, SC>*, LO > & clusterTree,
      std::vector< LO > & data
      ) const;

  //! returns pointer to the auxiliary array of curls

  Vector<LO, SCVT>* getCurls( ) const {
    return this->auxCurl;
  }

  /*!
   * Uniformly refines the mesh. 
   */
  void refine(
      int level = 1,
      int type = 2
      );

  /*!
   * Uniformly refines the mesh by quadrisection. 
   */
  void refine_2sect(
      int level
      );

  /*!
   * Uniformly refines the mesh by dividing each side to thirds. 
   */
  void refine_3sect(
      int level
      );

  /*!
   * Maps the mesh nodes to unit ball. 
   */
  void mapToUnitBall( );

  /*!
   * Appends this mesh by mesh. Does not perform any check if the resulting
   * mesh is computationally regular!
   * 
   * @param[ in ] mesh mesh to be appended the this mesh
   */
  void append(
      const SurfaceMesh3D< LO, SC >& mesh
      );

  /*!
   * Scales the mesh nodes by factor. 
   * 
   * @param[ in ] factor scaling factor
   */
  void scale(
      SCVT factor
      ) {
    this->scale( factor, factor, factor );
  }

  /*!
   * Scales the mesh nodes by different factors in each direction. 
   * 
   * @param[ in ] f1 x1 scaling factor
   * @param[ in ] f2 x2 scaling factor
   * @param[ in ] f3 x3 scaling factor
   */
  void scale(
      SCVT f1,
      SCVT f2,
      SCVT f3
      );

  /*!
   * Move centroid to [0,0,0], scale, move back
   * 
   * @param[ in ] factor scaling factor
   */
  void scaleAroundCentroid(
      SCVT factor
      ) {
    this->scaleAroundCentroid( factor, factor, factor );
  }

  /*!
   * Move centroid to [0,0,0], scale, move back
   * 
   * @param[ in ] f1 x1 scaling factor
   * @param[ in ] f2 x2 scaling factor
   * @param[ in ] f3 x3 scaling factor
   */
  void scaleAroundCentroid(
      SCVT f1,
      SCVT f2,
      SCVT f3
      );

  /*!
   * Rotates mesh about the axis u by the angle alpha
   * 
   * @param[ in ] u rotation vector
   * @param[ in ] alpha rotation angle
   */
  void rotate(
      SCVT * v,
      SCVT alpha
      );

  //! Moves centroid to the origin
  void moveCentroidToOrigin( );

  /*!
   * Moves the mesh nodes in direction [x1,x2,x3]. 
   * 
   * @param[ in ] x1 shift in x1 direction
   * @param[ in ] x2 shift in x2 direction
   * @param[ in ] x2 shift in x3 direction
   */
  void move(
      SCVT x1,
      SCVT x2,
      SCVT x3
      );

  /*!
   * Moves the node idx in direction [x1,x2,x3]. 
   * 
   * @param[ in ] idx 
   * @param[ in ] x1 shift in x1 direction
   * @param[ in ] x2 shift in x2 direction
   * @param[ in ] x2 shift in x3 direction
   */
  void move(
      LO idx,
      SCVT x1,
      SCVT x2,
      SCVT x3
      );

  /*!
   * Expands the mesh nodes in normal direction [x1,x2,x3]. 
   * 
   * @param[ in ] alpha moves in alpha*normal
   */
  void expand(
      SCVT alpha = 1.0
      );

  /*!
   * Computes and returns a radius of a mesh
   */
  SCVT getRadius( );

  //! initilizes matrices T for Lame equation
  void assembleT(
      SparseMatrix< LO, SCVT > & T12,
      SparseMatrix< LO, SCVT > & T13,
      SparseMatrix< LO, SCVT > & T23,
      SparseMatrix< LO, SCVT > & T
      );

  inline vector<SCVT>* getNormals( ) {
    return & this->normals;
  }

  LO getL2gElem( LO idx ) {
    return l2gElems[idx];
  }

  LO getL2gNode( LO idx ) {
    return l2gNodes[idx];
  }

protected:

  //! recursive method for creation of cluster tree
  void doNestedDissection(
      Tree<BECluster<LO, SC>*,
      LO >& clusterTree,
      LO maxElemsPerCluster
      ) const;

  void doNestedProjectionMassRatio(
      Tree<BECluster<LO, SC>*, LO > & clusterTree,
      std::vector< std::vector< SC > > & data
      ) const;

  void doNestedProjectionMassRatio(
      Tree<BECluster<LO, SC>*, LO > & clusterTree,
      std::vector< Vector< LO, SC > * > & data
      ) const;

  void doGetNestedProjectionDiam(
      Tree<BECluster<LO, SC>*, LO > & clusterTree,
      std::vector< std::vector< SC > > & data
      ) const;

  //! center of mass of each element
  vector<SCVT>* centroids;

  //! vector of outer normals to each element [n1_x, n1_y, n1_z, n2_x, n2_y, ...]
  vector<SCVT> normals;

  vector< SCVT > additiveCurvature;

  //! mapping from element to its edges
  vector<LO> elem2edges;

  //! mapping from nodes to elements
  vector< vector< LO > > node2elems;

  //! local coordinate vectors for all orientations
  vector<SCVT> r1, r2, alpha1, alpha2, stau;

  //! vector of curls of p1 functions on all elements
  Vector<LO, SCVT> *auxCurl;

  //! mapping from local indexing to global indexing (for MPI distributed)
  vector<LO> l2gElems;

  //! mapping from local indexing to global indexing (for MPI distributed)
  vector<LO> l2gNodes;

  //! initializes array of edges
  void initEdges( );

private:

  //! get local coordinates for a specific rotation of a triangle
  void getLocalCoordinates(
      const SCVT* x1,
      const SCVT* x2,
      const SCVT* x3,
      SCVT* r1,
      SCVT* r2,
      SCVT* stau,
      SCVT* alpha1,
      SCVT* alpha2
      ) const;

  //! initializes vector of curls
  void initCurl( );

  //! compute curvature in one node
  SCVT computeNodalAdditiveCurvature( LO i );

};

}

// include .cpp file to overcome linking problems due to templates
#include "SurfaceMesh3D.cpp"

#endif	/* SURFACEMESH3D_H */

