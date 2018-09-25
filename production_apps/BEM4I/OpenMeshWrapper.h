/*!
 * @file    OpenMeshWrapper.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    May 13, 2014
 * @brief   Header file for class OpenMeshWrapper
 * 
 */

#ifndef OPENMESHWRAPPER_H
#define	OPENMESHWRAPPER_H

#include <iostream>
#include <fstream>
#include <vector>

#include "SurfaceMesh3D.h"
#include "Macros.h"
#include "MeshTags.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "LoopT.h"

namespace bem4i {

template< class LO, class SC >
struct MyTraits : public OpenMesh::DefaultTraits {

  typedef typename GetType<LO, SC>::SCVT SCVT;

  typedef OpenMesh::VectorT< SCVT, 3 > Point;
  typedef OpenMesh::VectorT< SCVT, 3 > Normal;

  VertexAttributes( OpenMesh::Attributes::Normal );
  FaceAttributes( OpenMesh::Attributes::Normal );

  VertexTraits{
    public:

    VertexT( ) : vtag( VERTEX_SMOOTH ) {
    };

    void tag( VertexTag vtag ) {
      this->vtag = vtag;
    }

    VertexTag tag( ) {
      return this->vtag;
    }

    private:

    VertexTag vtag;
  };

  EdgeTraits{
    public:

    EdgeT( ) : etag( EDGE_SMOOTH ) {
    };

    void tag( EdgeTag etag ) {
      this->etag = etag;
    }

    EdgeTag tag( ) {
      return this->etag;
    }

    private:

    EdgeTag etag;
  };
};

/*! 
 * Class representing a surface mesh of a 3D body 
 * 
 */
template<class LO, class SC>
class OpenMeshWrapper {

  typedef typename GetType<LO, SC>::SCVT SCVT;
  typedef OpenMesh::TriMesh_ArrayKernelT< MyTraits< LO, SC > > MyTriMesh;

  typedef typename MyTriMesh::Point myPoint;
  typedef typename MyTriMesh::VertexHandle myVertexHandle;
  typedef typename MyTriMesh::EdgeHandle myEdgeHandle;
  typedef typename MyTriMesh::HalfedgeHandle myHalfedgeHandle;

public:
  //! default constructor
  OpenMeshWrapper( );

  //! read from  file
  OpenMeshWrapper(
      const string& meshFile
      );

  //! copy constructor (only creates nodes and elems structure
  // (similar to construction from file or SurfaceMesh3D)
  OpenMeshWrapper(
      const OpenMeshWrapper< LO, SC > & orig
      );

  //! constructor from SurfaceMesh3D
  OpenMeshWrapper(
      const SurfaceMesh3D< LO, SC > & orig
      );

  //! destructor
  ~OpenMeshWrapper( );

  void copy(
      OpenMeshWrapper< LO, SC > & copy
      ) const;

  //! returns number of nodes

  inline LO getNNodes( ) const {
    return mesh->n_vertices( );
  }

  //! returns number of edges

  inline LO getNEdges( ) const {
    return mesh->n_edges( );
  }

  //! returns number of elements

  inline LO getNElements( ) const {
    return mesh->n_faces( );
  }

  /*!
   * Returns element of a mesh
   * 
   * @param   idx index of element to return
   * @param   element pointer to <em>preallocated</em> array of LOs where the element indices will be copied
   */
  inline void getElement(
      LO idx,
      LO* element
      ) const {

    std::vector<typename MyTriMesh::VertexHandle> vhandles;
    for ( typename MyTriMesh::CFVIter fv_it =
        mesh->cfv_iter( typename MyTriMesh::FaceHandle( idx ) );
        fv_it.is_valid( ); ++fv_it ) {
      vhandles.push_back( *fv_it );
    }
    element[ 0 ] = vhandles[ 0 ].idx( );
    element[ 1 ] = vhandles[ 1 ].idx( );
    element[ 2 ] = vhandles[ 2 ].idx( );
  }

  /*!
   * Returns node of a mesh
   * 
   * @param   idx index of node to return
   * @param   node pointer to <em>preallocated</em> array of doubles where
   *  the node coordinates will be copied
   */
  inline void getNode(
      LO idx,
      SCVT * node
      ) const {

    typename MyTriMesh::Point v;
    v = mesh->point( typename MyTriMesh::VertexHandle( idx ) );
    node[ 0 ] = v[ 0 ];
    node[ 1 ] = v[ 1 ];
    node[ 2 ] = v[ 2 ];
  }

  /*!
   * Returns nodes of a mesh
   * 
   * @param   nodes pointer to <em>preallocated</em> array of doubles where
   *  the node coordinates will be copied
   */
  inline void getNodes(
      SCVT * nodes
      ) const {

    typename MyTriMesh::Point v;
    LO nNodes = this->getNNodes( );

    for ( LO i = 0; i < nNodes; ++i ) {
      v = mesh->point( typename MyTriMesh::VertexHandle( i ) );
      nodes[ 3 * i ] = v[ 0 ];
      nodes[ 3 * i + 1 ] = v[ 1 ];
      nodes[ 3 * i + 2 ] = v[ 2 ];
    }
  }
  
  /*!
   * Returns nodes of a mesh
   * 
   * @param   nodes vector of doubles where the node coordinates will be copied
   */
  inline void getNodes(
      std::vector< SCVT > & nodes
      ) const {

    typename MyTriMesh::Point v;
    LO nNodes = this->getNNodes( );
    
    nodes.clear();
    nodes.reserve( 3 * nNodes );

    for ( LO i = 0; i < nNodes; ++i ) {
      v = mesh->point( typename MyTriMesh::VertexHandle( i ) );
      nodes.push_back( v[ 0 ] );
      nodes.push_back( v[ 1 ] );
      nodes.push_back( v[ 2 ] );
    }
  }

  /*!
   * Returns index of the edge vertexStart to vertexEnd (either direction)
   * 
   * @param vertexStart start point of the edge
   * @param vertexEnd end point of the edge
   */
  inline LO getEdge(
      LO vertexStart,
      LO vertexEnd
      ) const {

    myVertexHandle vStart( vertexStart );
    myVertexHandle vEnd( vertexEnd );
    myHalfedgeHandle heHandle = this->mesh->find_halfedge( vStart, vEnd );
    myEdgeHandle eHandle = this->mesh->edge_handle( heHandle );

    return eHandle.idx( );
  }

  /*!
   * Returns indices of the vertices building the edge
   * 
   * @param eInd index of the edge
   * @param vertexStart start point of the edge
   * @param vertexEnd end point of the edge
   */
  inline void getEdge(
      LO eInd,
      LO & vertexStart,
      LO & vertexEnd
      ) const {

    myHalfedgeHandle heHandle = this->mesh->halfedge_handle(
        myEdgeHandle( eInd ), 0 );
    vertexStart = this->mesh->to_vertex_handle( heHandle ).idx( );
    heHandle = this->mesh->halfedge_handle(
        myEdgeHandle( eInd ), 1 );
    vertexEnd = this->mesh->to_vertex_handle( heHandle ).idx( );
  }

  inline void getNodalNormal(
      LO idx,
      SCVT * n
      ) const {

    typename MyTriMesh::Normal v;
    v = this->mesh->normal( typename MyTriMesh::VertexHandle( idx ) );
    n[ 0 ] = v[ 0 ];
    n[ 1 ] = v[ 1 ];
    n[ 2 ] = v[ 2 ];
  };

  SCVT getVertexAverageEdgeLength(
      LO idx
      ) const;

  /*!
   * Loads mesh from a default file format
   * 
   * @param   meshFile a string with a path to the mesh file
   */
  bool load(
      const string& meshFile,
      SCVT scaleFactor = 1.0
      );

  /*!
   * Loads mesh from a default file format
   * 
   * @param   meshFile a string with a path to the mesh file
   */
  bool loadSmf(
      const string& meshFile
      );

  /*!
   * Reads tag file for edges and vertices
   * 
   * @param   tagFile a string with a path to the tag file
   */
  bool readTags(
      const string& meshFile
      );

  /*!
   * Prints mesh info to stdout
   */
  void printInfo( ) const;

  /*!
   * Prints mesh to the input format
   * 
   * @param[in]   meshFile string with the target file name
   */
  void print(
      const string& meshFile
      ) const;

  /*!
   * Prints mesh to the off format
   * 
   * @param[in]   meshFile string with the target file name
   */
  void printOff(
      const string& meshFile
      ) const;

  /*!
   * Prints mesh to the xml paraview file format
   * 
   * @param[in]   meshFile string with the target file name
   */
  void printParaviewVtu(
      const string& meshFile
      ) const;

  /*!
   * Prints mesh to the xml paraview file format
   * 
   * @param[in] meshFile string with the target file name
   */
  void printParaviewVtu(
      const string& meshFile,
      int nNodal,
      string* nodeNames,
      Vector< LO, SC >** nodalData,
      int nElem,
      string* elemNames,
      Vector< LO, SC >** elemData
      ) const;

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
      const std::vector< Vector< LO, SC >* > * nodalVData = nullptr
      ) const;

  /*!
   * Creates a SurfaceMesh3D instance
   * 
   * @param[in,out] copy reference to a SurfaceMesh3D instance
   */
  void createSurfaceMesh3D(
      SurfaceMesh3D< LO, SC > & copy
      ) const;

  //  /*!
  //   * Subdivides the mesh n times using the Loop subdivision scheme
  //   * 
  //   * @param[in,out] n number of subdivisions
  //   */
  //  void subdivideLoop(
  //      int n,
  //      bool update = true
  //      );

  /*!
   * Subdivides the mesh n times using the Loop subdivision scheme,
   * get subdivision matrices
   * 
   * @param[in,out] n number of subdivisions
   */
  void subdivideLoop(
      int n = 1,
      bool update = true,
      std::vector< SparseMatrix< LO, SCVT >* > * matrices = nullptr
      );

  //  /*!
  //   * Subdivides the mesh n times using the Loop subdivision scheme,
  //   * get subdivision matrix (one for levels 0 -> n)
  //   * 
  //   * @param[in,out] n number of subdivisions
  //   */
  //  void subdivideLoop(
  //      int n,
  //      SparseMatrix< int, SCVT > * matrix,
  //      bool update = true
  //      );

  /*!
   * Updates node position
   * 
   * @param[in] idx node index
   * @param[in] node new coordinates 
   */
  void setNode(
      LO idx,
      const SCVT * node
      );

  /*!
   * Updates position of all nodes
   * 
   * @param[in] node new coordinates 
   */
  void setNodes(
      const SCVT * nodes
      );

  VertexTag getNodeTag(
      LO idx
      ) const {

    typename MyTriMesh::VertexHandle vh( idx );

    return this->mesh->data( vh ).tag( );
  }

  void setNodeTag(
      LO idx,
      VertexTag vtag
      ) {

    typename MyTriMesh::VertexHandle vh( idx );
    this->mesh->data( vh ).tag( vtag );
  }

  EdgeTag getEdgeTag(
      LO idx
      ) const {

    typename MyTriMesh::EdgeHandle eh( idx );

    return this->mesh->data( eh ).tag( );
  }

  void setEdgeTag(
      LO idx,
      EdgeTag etag
      ) {

    typename MyTriMesh::EdgeHandle eh( idx );
    this->mesh->data( eh ).tag( etag );
  }

protected:

private:

  MyTriMesh * mesh;

  //! updates face and vertex normals (must be called after mesh modification)

  void updateNormals( ) {
    this->mesh->update_face_normals( );
    this->mesh->update_vertex_normals( );
  }

};

}

// include .cpp file to overcome linking problems due to templates
#include "OpenMeshWrapper.cpp"

#endif	/* OPENMESHWRAPPER_H */

