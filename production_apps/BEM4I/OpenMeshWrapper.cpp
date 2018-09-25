/*!
 * @file    OpenMeshWrapper.cpp
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    May 13, 2014
 * 
 */

#ifdef OPENMESHWRAPPER_H

namespace bem4i {

template<class LO, class SC>
OpenMeshWrapper<LO, SC>::OpenMeshWrapper( ) {
  this->mesh = nullptr;
}

template<class LO, class SC>
OpenMeshWrapper<LO, SC>::OpenMeshWrapper(
    const string & meshFile
    ) {

  this->mesh = new MyTriMesh( );
  this->load( meshFile );

  this->updateNormals( );
}

template<class LO, class SC>
OpenMeshWrapper<LO, SC>::OpenMeshWrapper(
    const OpenMeshWrapper< LO, SC > & orig
    ) {

  this->mesh = new MyTriMesh( );

  SCVT node[ 3 ];
  LO elem[ 3 ];
  LO nNodes = orig.mesh->n_vertices( );
  LO nElems = orig.mesh->n_faces( );
  LO nEdges = orig.mesh->n_edges( );

  // copy nodes
  typename MyTriMesh::VertexHandle * vhandle =
      new typename MyTriMesh::VertexHandle[ nNodes ];
  for ( LO i = 0; i < nNodes; ++i ) {
    orig.getNode( i, node );
    vhandle[ i ] = mesh->add_vertex(
        typename MyTriMesh::Point( node[ 0 ], node[ 1 ], node[ 2 ] ) );
    this->setNodeTag( i, orig.getNodeTag( i ) );
  }

  // copy elements
  std::vector<typename MyTriMesh::VertexHandle> face_vhandles;
  for ( LO i = 0; i < nElems; ++i ) {
    orig.getElement( i, elem );
    face_vhandles.clear( );
    face_vhandles.push_back( vhandle[ elem[ 0 ] ] );
    face_vhandles.push_back( vhandle[ elem[ 1 ] ] );
    face_vhandles.push_back( vhandle[ elem[ 2 ] ] );
    mesh->add_face( face_vhandles );
  }

  delete [] vhandle;

  // tag edges
  LO start;
  LO end;
  EdgeTag etag;
  for ( LO i = 0; i < nEdges; ++i ) {
    orig.getEdge( i, start, end );
    etag = orig.getEdgeTag( i );
    this->setEdgeTag( this->getEdge( start, end ), etag );
  }

  this->updateNormals( );
}

template<class LO, class SC>
OpenMeshWrapper<LO, SC>::OpenMeshWrapper(
    const SurfaceMesh3D< LO, SC > & orig
    ) {

  this->mesh = new MyTriMesh( );

  SCVT node[ 3 ];
  LO elem[ 3 ];
  LO nNodes = orig.getNNodes( );
  LO nElems = orig.getNElements( );

  typename MyTriMesh::VertexHandle * vhandle =
      new typename MyTriMesh::VertexHandle[ nNodes ];
  for ( LO i = 0; i < nNodes; ++i ) {
    orig.getNode( i, node );
    vhandle[ i ] = mesh->add_vertex(
        typename MyTriMesh::Point( node[ 0 ], node[ 1 ], node[ 2 ] ) );
  }

  std::vector<typename MyTriMesh::VertexHandle> face_vhandles;
  for ( LO i = 0; i < nElems; ++i ) {
    orig.getElement( i, elem );
    face_vhandles.clear( );
    face_vhandles.push_back( vhandle[ elem[ 0 ] ] );
    face_vhandles.push_back( vhandle[ elem[ 1 ] ] );
    face_vhandles.push_back( vhandle[ elem[ 2 ] ] );
    mesh->add_face( face_vhandles );
  }

  delete [] vhandle;

  this->updateNormals( );
}

template<class LO, class SC>
OpenMeshWrapper<LO, SC>::~OpenMeshWrapper( ) {
  if ( this->mesh ) delete mesh;
}

template<class LO, class SC>
void OpenMeshWrapper<LO, SC>::copy(
    OpenMeshWrapper<LO, SC> & copy
    ) const {

  copy.mesh = new MyTriMesh( );

  SCVT node[ 3 ];
  LO elem[ 3 ];
  LO nNodes = this->mesh->n_vertices( );
  LO nElems = this->mesh->n_faces( );
  LO nEdges = this->mesh->n_edges( );

  // copy nodes
  typename MyTriMesh::VertexHandle * vhandle =
      new typename MyTriMesh::VertexHandle[ nNodes ];
  for ( LO i = 0; i < nNodes; ++i ) {
    this->getNode( i, node );
    vhandle[ i ] = copy.mesh->add_vertex(
        typename MyTriMesh::Point( node[ 0 ], node[ 1 ], node[ 2 ] ) );
    copy.setNodeTag( i, this->getNodeTag( i ) );
  }

  // copy elements
  std::vector<typename MyTriMesh::VertexHandle> face_vhandles;
  for ( LO i = 0; i < nElems; ++i ) {
    this->getElement( i, elem );
    face_vhandles.clear( );
    face_vhandles.push_back( vhandle[ elem[ 0 ] ] );
    face_vhandles.push_back( vhandle[ elem[ 1 ] ] );
    face_vhandles.push_back( vhandle[ elem[ 2 ] ] );
    copy.mesh->add_face( face_vhandles );
  }

  delete [] vhandle;

  // tag edges
  LO start;
  LO end;
  EdgeTag etag;
  for ( LO i = 0; i < nEdges; ++i ) {
    this->getEdge( i, start, end );
    etag = this->getEdgeTag( i );
    copy.setEdgeTag( this->getEdge( start, end ), etag );
  }

  copy.updateNormals( );
}

template<class LO, class SC>
typename bem4i::OpenMeshWrapper<LO, SC>::SCVT
OpenMeshWrapper<LO, SC>::getVertexAverageEdgeLength( LO idx ) const {

  SCVT ret = 0.0;
  int nEdges = 0;

  typename MyTriMesh::VertexHandle vh = this->mesh->vertex_handle( idx );
  typename MyTriMesh::VertexEdgeIter veit = this->mesh->ve_iter( vh );

  for (; veit.is_valid( ); ++veit ) {
    ret += this->mesh->calc_edge_length( *veit );
    ++nEdges;
  }

  ret /= nEdges;

  return ret;
}

template<class LO, class SC>
bool OpenMeshWrapper<LO, SC>::load(
    const string & meshFile,
    SCVT scale
    ) {

  if ( mesh ) delete mesh;
  mesh = new MyTriMesh( );

  std::cout << "Reading file '" << meshFile << "' ... ";
  std::ifstream file( meshFile.c_str( ) );

  if ( !file.is_open( ) ) {
    std::cout << "File '" << meshFile << "' could not be opened!" << std::endl;
    return false;
  }

  LO dummy;
  SCVT node[ 3 ];
  LO elem[ 3 ];
  LO nNodes;
  LO nElems;

  file >> dummy;
  file >> dummy;
  file >> nNodes;

  typename MyTriMesh::VertexHandle * vhandle =
      new typename MyTriMesh::VertexHandle[ nNodes ];
  for ( LO i = 0; i < nNodes; ++i ) {
    file >> node[ 0 ];
    file >> node[ 1 ];
    file >> node[ 2 ];
    vhandle[ i ] = mesh->add_vertex(
        typename MyTriMesh::Point(
        scale * node[ 0 ], scale * node[ 1 ], scale * node[ 2 ] ) );
  }

  file >> nElems;

  std::vector<typename MyTriMesh::VertexHandle> face_vhandles;
  for ( LO i = 0; i < nElems; ++i ) {
    file >> elem[ 0 ];
    file >> elem[ 1 ];
    file >> elem[ 2 ];
    face_vhandles.clear( );
    face_vhandles.push_back( vhandle[ elem[ 0 ] ] );
    face_vhandles.push_back( vhandle[ elem[ 1 ] ] );
    face_vhandles.push_back( vhandle[ elem[ 2 ] ] );
    mesh->add_face( face_vhandles );
  }

  delete [] vhandle;

  file.close( );

  std::cout << "done." << std::endl;

  this->updateNormals( );
  return true;
}

template<class LO, class SC>
bool OpenMeshWrapper<LO, SC>::loadSmf(
    const string& meshFile
    ) {

  if ( mesh ) delete mesh;
  mesh = new MyTriMesh( );

  std::cout << "Reading file '" << meshFile << "' ... ";
  std::ifstream file( meshFile.c_str( ) );

  if ( !file.is_open( ) ) {
    std::cout << "File '" << meshFile << "' could not be opened!" << std::endl;
    return false;
  }

  LO nNodes;
  LO nElems;

  // skip first two lines
  file.ignore( 50, '\n' );
  file.ignore( 50, '\n' );
  file >> nNodes;
  file >> nElems;

  SCVT node[ 3 ];
  typename MyTriMesh::VertexHandle * vhandle =
      new typename MyTriMesh::VertexHandle[ this->nNodes ];
  for ( LO i = 0; i < nNodes; ++i ) {
    file >> node[ 0 ];
    file >> node[ 1 ];
    file >> node[ 2 ];
    vhandle[ i ] = mesh->add_vertex(
        typename MyTriMesh::Point( node[ 0 ], node[ 1 ], node[ 2 ] ) );
  }

  LO elem[ 3 ];
  std::vector<typename MyTriMesh::VertexHandle> face_vhandles;
  for ( LO i = 0; i < nElems; ++i ) {
    file >> elem[ 0 ];
    file >> elem[ 1 ];
    file >> elem[ 2 ];
    face_vhandles.clear( );
    face_vhandles.push_back( vhandle[ elem[ 0 ] ] );
    face_vhandles.push_back( vhandle[ elem[ 1 ] ] );
    face_vhandles.push_back( vhandle[ elem[ 2 ] ] );
    mesh->add_face( face_vhandles );
  }

  delete [] vhandle;

  file.close( );

  std::cout << "done." << std::endl;

  this->updateNormals( );
  return true;
}

template<class LO, class SC>
bool OpenMeshWrapper<LO, SC>::readTags(
    const string& meshFile
    ) {

  std::cout << "Reading file '" << meshFile << "' ... ";
  std::ifstream file( meshFile.c_str( ) );

  if ( !file.is_open( ) ) {
    std::cout << "File '" << meshFile << "' could not be opened!" << std::endl;
    return false;
  }

  LO nVTags;
  LO nETags;
  LO vIndStart;
  LO vIndEnd;
  LO eInd;
  int tag;

  file >> nVTags;

  for ( LO i = 0; i < nVTags; ++i ) {
    file >> vIndStart;
    file >> tag;
    if ( vIndStart >= this->getNNodes( ) || tag >= VertexTag::NO_VERTEX_TAGS )
      return false;
    this->setNodeTag( vIndStart, (VertexTag) tag );
  }

  file >> nETags;

  for ( LO i = 0; i < nETags; ++i ) {
    file >> vIndStart;
    file >> vIndEnd;
    file >> tag;
    eInd = this->getEdge( vIndStart, vIndEnd );
    if ( eInd >= this->getNEdges( ) || tag >= EdgeTag::NO_EDGE_TAGS )
      return false;
    this->setEdgeTag( eInd, (EdgeTag) tag );
  }
  // TODO: return to original state if unsuccessful

  file.close( );

  std::cout << "done." << std::endl;

  return true;
}

template<class LO, class SC>
void OpenMeshWrapper<LO, SC>::printInfo( ) const {
  std::cout << "Mesh type: OpenMeshWrapper, ";
  std::cout << "dim: " << 3;
  std::cout << ", nodes: " << mesh->n_vertices( );
  std::cout << ", elements: " << mesh->n_faces( );
  std::cout << ", edges: " << mesh->n_edges( );
  std::cout << "." << std::endl;
}

template<class LO, class SC>
void OpenMeshWrapper<LO, SC>::printParaviewVtu(
    const string& meshFile
    ) const {

  LO nNodes = mesh->n_vertices( );
  LO nElems = mesh->n_faces( );

  std::cout << "Printing '" << meshFile << "' ... ";

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 10 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File '" << meshFile << "' could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"
      << std::endl;
  file_vtu << "  <UnstructuredGrid>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << nNodes <<
      "\" NumberOfCells=\"" << nElems << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\" "
      "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  SCVT node[ 3 ];
  for ( LO i = 0; i < nNodes; i++ ) {
    this->getNode( i, node );
    file_vtu << "          "
        << node[ 0 ] << " "
        << node[ 1 ] << " "
        << node[ 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Cells>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\" "
      "format=\"ascii\">" << std::endl;

  LO elem[ 3 ];
  for ( LO i = 0; i < nElems; i++ ) {
    this->getElement( i, elem );

    file_vtu << "          "
        << elem[ 0 ] << " "
        << elem[ 1 ] << " "
        << elem[ 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"offsets\" "
      "format=\"ascii\">" << std::endl;

  for ( LO offset = 1; offset <= nElems; offset++ ) {
    file_vtu << "          " << offset * 3 << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"UInt8\" Name=\"types\" "
      "format=\"ascii\">" << std::endl;
  for ( LO i = 1; i <= nElems; i++ ) {
    file_vtu << "          5" << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Cells>" << std::endl;

  file_vtu << "    </Piece>" << std::endl;
  file_vtu << "  </UnstructuredGrid>" << std::endl;
  file_vtu << "</VTKFile>" << std::endl;
  file_vtu.close( );

  std::cout << "done." << std::endl;
}

template<class LO, class SC>
void OpenMeshWrapper<LO, SC>::printParaviewVtu(
    const string& meshFile,
    int nNodal,
    string* nodeNames,
    Vector< LO, SC >** nodalData,
    int nElem, string* elemNames,
    Vector< LO, SC >** elemData
    ) const {

  LO nNodes = mesh->n_vertices( );
  LO nElems = mesh->n_faces( );

  std::cout << "Printing '" << meshFile << "' ... ";

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 10 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File '" << meshFile << "' could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"
      << std::endl;
  file_vtu << "  <UnstructuredGrid>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << nNodes <<
      "\" NumberOfCells=\"" << nElems << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\" "
      "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  SCVT node[ 3 ];
  for ( LO i = 0; i < nNodes; i++ ) {
    this->getNode( i, node );
    file_vtu << "          "
        << node[ 0 ] << " "
        << node[ 1 ] << " "
        << node[ 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Cells>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\" "
      "format=\"ascii\">" << std::endl;

  LO elem[ 3 ];
  for ( LO i = 0; i < nElems; i++ ) {
    this->getElement( i, elem );

    file_vtu << "          "
        << elem[ 0 ] << " "
        << elem[ 1 ] << " "
        << elem[ 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"offsets\" "
      "format=\"ascii\">" << std::endl;

  for ( LO offset = 1; offset <= nElems; offset++ ) {
    file_vtu << "          " << offset * 3 << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"UInt8\" Name=\"types\" "
      "format=\"ascii\">" << std::endl;
  for ( LO i = 1; i <= nElems; i++ ) {
    file_vtu << "          5" << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Cells>" << std::endl;

  string header;
  if ( nNodal > 0 ) {
    string header = nodeNames[ 0 ];
    for ( int j = 1; j < nNodal; j++ ) {
      header += "," + nodeNames[ j ];
    }
    file_vtu << "      <PointData Scalars=\"" + header + "\">" << std::endl;
    for ( int j = 0; j < nNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" + nodeNames[ j ]
          + "\" format=\"ascii\">" << std::endl;
      for ( LO i = 0; i < nNodes; i++ ) {
        file_vtu << "          " << nodalData[ j ]->get( i ) << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }
    file_vtu << "      </PointData>" << std::endl;
  }

  if ( nElem > 0 ) {
    string header = elemNames[ 0 ];
    for ( int j = 1; j < nElem; j++ ) {
      header += "," + elemNames[ j ];
    }
    file_vtu << "      <CellData Scalars=\"" + header + "\">" << std::endl;
    for ( int j = 0; j < nElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" + elemNames[ j ]
          + "\" format=\"ascii\">" << std::endl;
      for ( LO i = 0; i < nElems; i++ ) {
        file_vtu << "          " << elemData[ j ]->get( i ) << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }
    file_vtu << "      </CellData>" << std::endl;
  }

  file_vtu << "    </Piece>" << std::endl;
  file_vtu << "  </UnstructuredGrid>" << std::endl;
  file_vtu << "</VTKFile>" << std::endl;
  file_vtu.close( );

  std::cout << "done." << std::endl;
}

template<class LO, class SC>
void OpenMeshWrapper<LO, SC>::printParaviewVtu(
    const string & meshFile,
    const std::vector< string > * nodeNames,
    const std::vector< Vector< LO, SC >* > * nodalData,
    const std::vector< string > * elemNames,
    const std::vector< Vector< LO, SC >* > * elemData,
    const std::vector< string > * nodeVNames,
    const std::vector< Vector< LO, SC >* > * nodalVData
    ) const {

  LO nNodes = mesh->n_vertices( );
  LO nElems = mesh->n_faces( );

  std::cout << "Printing '" << meshFile << "' ... ";

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 10 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File '" << meshFile << "' could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"
      << std::endl;
  file_vtu << "  <UnstructuredGrid>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << nNodes <<
      "\" NumberOfCells=\"" << nElems << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\" "
      "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  SCVT node[ 3 ];
  for ( LO i = 0; i < nNodes; i++ ) {
    this->getNode( i, node );
    file_vtu << "          "
        << node[ 0 ] << " "
        << node[ 1 ] << " "
        << node[ 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Cells>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\" "
      "format=\"ascii\">" << std::endl;

  LO elem[ 3 ];
  for ( LO i = 0; i < nElems; i++ ) {
    this->getElement( i, elem );

    file_vtu << "          "
        << elem[ 0 ] << " "
        << elem[ 1 ] << " "
        << elem[ 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"offsets\" "
      "format=\"ascii\">" << std::endl;

  for ( LO offset = 1; offset <= nElems; offset++ ) {
    file_vtu << "          " << offset * 3 << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"UInt8\" Name=\"types\" "
      "format=\"ascii\">" << std::endl;
  for ( LO i = 1; i <= nElems; i++ ) {
    file_vtu << "          5" << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Cells>" << std::endl;

  int nNodal = 0;
  if ( nodalData ) nNodal = nodalData->size( );
  int nVNodal = 0;
  if ( nodalVData ) nVNodal = nodalVData->size( );

  string header, vheader;
  if ( nNodal > 0 || nVNodal > 0 ) {
    file_vtu << "      <PointData ";

    if ( nNodal > 0 ) {
      header = ( *nodeNames )[ 0 ];
      for ( int j = 1; j < nNodal; j++ ) {
        header += "," + ( *nodeNames )[ j ];
      }
      file_vtu << "Scalars=\"" + header;
    }

    if ( nVNodal > 0 ) {
      vheader = ( *nodeVNames )[ 0 ];
      for ( int j = 1; j < nVNodal; j++ ) {
        vheader += "," + ( *nodeVNames )[ j ];
      }
      file_vtu << " Vectors=\"" + vheader;
    }

    file_vtu << "\">" << std::endl;
    for ( int j = 0; j < nNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
          ( *nodeNames )[ j ] + "\" format=\"ascii\">" << std::endl;
      for ( LO i = 0; i < nNodes; i++ ) {
        file_vtu << "          " << ( *nodalData )[ j ]->get( i ) << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    for ( int j = 0; j < nVNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
          ( *nodeVNames )[ j ] + "\" NumberOfComponents=\"3" + "\" "
          "format=\"ascii\">" << std::endl;
      for ( LO i = 0; i < 3 * nNodes; i++ ) {
        file_vtu << "          " << ( *nodalVData )[ j ]->get( i ) << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }
    file_vtu << "      </PointData>" << std::endl;
  }

  int nElem = 0;
  if ( elemData ) nElem = elemData->size( );
  if ( nElem > 0 ) {
    header.clear( );
    header = ( *elemNames )[ 0 ];
    for ( int j = 1; j < nElem; j++ ) {
      header += "," + ( *elemNames )[ j ];
    }
    file_vtu << "      <CellData Scalars=\"" + header + "\">" << std::endl;
    for ( int j = 0; j < nElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
          ( *elemNames )[ j ] + "\" format=\"ascii\">" << std::endl;
      for ( LO i = 0; i < nElems; i++ ) {
        file_vtu << "          " << ( *elemData )[ j ]->get( i ) << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }
    file_vtu << "      </CellData>" << std::endl;
  }

  file_vtu << "    </Piece>" << std::endl;
  file_vtu << "  </UnstructuredGrid>" << std::endl;
  file_vtu << "</VTKFile>" << std::endl;
  file_vtu.close( );

  std::cout << "done." << std::endl;
}

template<class LO, class SC>
void OpenMeshWrapper<LO, SC>::print(
    const string & meshFile
    ) const {

  LO nNodes = mesh->n_vertices( );
  LO nElems = mesh->n_faces( );

  std::cout << "Printing '" << meshFile << "' ... ";

  std::ofstream file( meshFile.c_str( ) );

  file.setf( std::ios::showpoint | std::ios::scientific );
  file.precision( 10 );

  if ( !file.is_open( ) ) {
    std::cout << "File '" << meshFile << "' could not be opened!" << std::endl;
    return;
  }

  file << 3 << std::endl << 3 << std::endl << std::endl;
  file << nNodes << std::endl;

  SCVT node[ 3 ];
  for ( LO i = 0; i < nNodes; ++i ) {
    this->getNode( i, node );
    file << node[ 0 ] << " "
        << node[ 1 ] << " "
        << node[ 2 ] << std::endl;
  }

  file << std::endl << nElems << std::endl;

  LO elem[ 3 ];
  for ( LO i = 0; i < nElems; ++i ) {
    this->getElement( i, elem );
    file << elem[ 0 ] << " "
        << elem[ 1 ] << " "
        << elem[ 2 ] << std::endl;
  }

  file.close( );

  std::cout << "done." << std::endl;
}

template<class LO, class SC>
void OpenMeshWrapper<LO, SC>::printOff(
    const string& meshFile
    ) const {

  try {
    if ( !OpenMesh::IO::write_mesh( *mesh, meshFile.c_str( ) ) ) {
      std::cout << "Cannot write mesh to file '" << meshFile << "'."
          << std::endl;
    }
  } catch ( std::exception& x ) {
    std::cout << x.what( ) << std::endl;
  }
}

template<class LO, class SC>
void OpenMeshWrapper<LO, SC>::createSurfaceMesh3D(
    SurfaceMesh3D< LO, SC > & copy
    ) const {

  LO nNodes = mesh->n_vertices( );
  LO nElems = mesh->n_faces( );

  copy.dim = 3;
  copy.nNodesPerElem = 3;
  copy.nNodes = nNodes;
  copy.nElems = nElems;

  copy.nodes.clear( );
  copy.nodes.reserve( copy.dim * copy.nNodes );
  copy.elems.clear( );
  copy.elems.reserve( copy.nNodesPerElem * copy.nElems );

  LO elem[ 3 ];
  SCVT node[ 3 ];

  for ( LO i = 0; i < nNodes; i++ ) {
    this->getNode( i, node );
    copy.nodes.push_back( node[ 0 ] );
    copy.nodes.push_back( node[ 1 ] );
    copy.nodes.push_back( node[ 2 ] );
  }

  for ( LO i = 0; i < nElems; i++ ) {
    this->getElement( i, elem );
    copy.elems.push_back( elem[ 0 ] );
    copy.elems.push_back( elem[ 1 ] );
    copy.elems.push_back( elem[ 2 ] );
  }

  copy.auxCurl = nullptr;
  copy.initArea( );
  copy.initNormals( );
  copy.initLocalCoordinates( );
  copy.initEdges( );
  copy.initCurl( );
}

//template<class LO, class SC>
//void OpenMeshWrapper<LO, SC>::subdivideLoop(
//    int n,
//    bool update
//    ) {
//
//  this->subdivideLoop( n, nullptr, update );
//}

template<class LO, class SC>
void OpenMeshWrapper<LO, SC>::subdivideLoop(
    int n,
    bool update,
    std::vector< SparseMatrix< LO, SCVT >* > * matrices
    ) {

  OpenMesh::Subdivider::Uniform::LoopT< MyTriMesh, SCVT, LO > loop;
  loop.attach( *mesh );
  loop( n, matrices, update );
  loop.detach( );

  this->updateNormals( );
}

template<class LO, class SC>
void OpenMeshWrapper<LO, SC>::setNode(
    LO idx,
    const SCVT * node
    ) {
  this->mesh->set_point( typename MyTriMesh::VertexHandle( idx ),
      typename MyTriMesh::Point( node[ 0 ], node[ 1 ], node[ 2 ] ) );

  // TODO: not necessary after each setNode?
  this->updateNormals( );
}

template<class LO, class SC>
void OpenMeshWrapper<LO, SC>::setNodes(
    const SCVT * nodes
    ) {

  LO nNodes = this->mesh->n_vertices( );
  for ( LO i = 0; i < nNodes; ++i ) {
    this->mesh->set_point( typename MyTriMesh::VertexHandle( i ),
        typename MyTriMesh::Point( nodes[ 3 * i ], nodes[ 3 * i + 1 ],
        nodes[ 3 * i + 2 ] ) );
  }

  this->updateNormals( );
}

} // end namespace bem4i

#endif /*OPENMESHWRAPPER_H*/
