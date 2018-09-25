/*!
 * @file    SurfaceMesh3D.cpp
 * @author  Michal Merta
 * @author  Jan Zapletal
 * @date    July 5, 2013
 *
 */

#ifdef SURFACEMESH3D_H

namespace bem4i {

template<class LO, class SC>
SurfaceMesh3D<LO, SC>::SurfaceMesh3D( ) {

  this->dim = 3;
  this->nNodes = 0;
  this->nElems = 0;
  this->auxCurl = nullptr;
}

template<class LO, class SC>
SurfaceMesh3D<LO, SC>::SurfaceMesh3D(
  const string & meshFile
  ) {

  this->auxCurl = nullptr;
  this->load( meshFile );
}

template<class LO, class SC>
SurfaceMesh3D<LO, SC>::SurfaceMesh3D(
  const SurfaceMesh3D& orig
  ) {

  this->dim = orig.dim;
  this->nNodes = orig.nNodes;
  this->nElems = orig.nElems;
  this->nNodesPerElem = orig.nNodesPerElem;
  this->auxCurl = nullptr;

  this->nodes = orig.nodes;
  this->elems = orig.elems;

  initArea( );
  initNormals( );
  initLocalCoordinates( );
  initEdges( );
  initCurl( );
  initNode2Elems( );
  initL2g( );
}

template<class LO, class SC>
SurfaceMesh3D<LO, SC>::SurfaceMesh3D(
  const std::vector< SCVT > & nodes,
  const std::vector< LO > & elems
  ) {

  this->dim = 3;
  this->nNodes = nodes.size( ) / 3;
  this->nElems = elems.size( ) / 3;
  this->nNodesPerElem = 3;
  this->auxCurl = nullptr;

  this->nodes = nodes;
  this->elems = elems;

  initArea( );
  initNormals( );
  initLocalCoordinates( );
  initEdges( );
  initCurl( );
  initNode2Elems( );
  initL2g( );
}

template<class LO, class SC>
SurfaceMesh3D<LO, SC>::~SurfaceMesh3D( ) {

  if ( this->auxCurl ) {
    delete auxCurl;
  }
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::copyToComplex(
  SurfaceMesh3D< LO, std::complex< SCVT > > & copy
  ) const {

  copy.dim = this->dim;
  copy.nNodes = this->nNodes;
  copy.nElems = this->nElems;
  copy.nNodesPerElem = this->nNodesPerElem;
  copy.auxCurl = nullptr;

  copy.nodes = this->nodes;
  copy.elems = this->elems;

  copy.initArea( );
  copy.initNormals( );
  copy.initLocalCoordinates( );
  copy.initEdges( );
  copy.initCurl( );
  copy.initNode2Elems( );
  copy.initL2g( );
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::copyToReal(
  SurfaceMesh3D< LO, SCVT > & copy
  ) const {

  copy.dim = this->dim;
  copy.nNodes = this->nNodes;
  copy.nElems = this->nElems;
  copy.nNodesPerElem = this->nNodesPerElem;
  copy.auxCurl = nullptr;

  copy.nodes = this->nodes;
  copy.elems = this->elems;

  copy.initArea( );
  copy.initNormals( );
  copy.initLocalCoordinates( );
  copy.initEdges( );
  copy.initCurl( );
  copy.initNode2Elems( );
}

template<class LO, class SC>
bool SurfaceMesh3D<LO, SC>::load(
    const string & meshFile,
    SCVT scaleFactor
    ) {
#ifdef VERBOSE
  std::cout << "Reading file " << meshFile << " ... ";
#endif
  std::ifstream file( meshFile.c_str( ) );

  if ( !file.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return false;
  }

  file >> this->dim;
  file >> this->nNodesPerElem;
  file >> this->nNodes;
  this->nodes.resize( 3 * this->nNodes );

  for ( LO i = 0; i < this->nNodes; i++ ) {
    file >> this->nodes[3 * i];
    this->nodes[3 * i] *= scaleFactor;
    file >> this->nodes[3 * i + 1];
    this->nodes[3 * i + 1] *= scaleFactor;
    file >> this->nodes[3 * i + 2];
    this->nodes[3 * i + 2] *= scaleFactor;
  }

  file >> this->nElems;
  this->elems.resize( this->nElems * this->nNodesPerElem );

  for ( LO i = 0; i < this->nElems; i++ ) {
    for ( int j = 0; j < this->nNodesPerElem; j++ ) {
      file >> this->elems[ this->nNodesPerElem * i + j ];
    }
  }

  file.close( );

  //  initialize all auxiliary variables
  initArea( );
  initNormals( );
  initLocalCoordinates( );
  initEdges( );
  initCurl( );
  initNode2Elems( );
  initL2g( );

#ifdef VERBOSE
  std::cout << "done." << std::endl;
#endif

  return true;
}

template<class LO, class SC>
bool SurfaceMesh3D<LO, SC>::loadSmfOpenFTLOptimize(
  const string & meshFile,
  Vector< LO, SCVT > & fixity,
  Vector< LO, SCVT > & curvature
  ) {

#ifdef VERBOSE
  std::cout << "Loading mesh, fixity, and curvature from " << meshFile
    << " ... ";
#endif

  std::ifstream file( meshFile.c_str( ) );

  if ( !file.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return false;
  }

  this->dim = 3;
  this->nNodesPerElem = 3;

  // skip first two lines
  file.ignore( 50, '\n' );
  file.ignore( 50, '\n' );
  file >> this->nNodes;
  file >> this->nElems;
  this->nodes.resize( 3 * this->nNodes );

  for ( LO i = 0; i < this->nNodes; i++ ) {
    file >> this->nodes[3 * i];
    file >> this->nodes[3 * i + 1];
    file >> this->nodes[3 * i + 2];
  }

  this->elems.resize( this->nElems * this->nNodesPerElem );

  for ( LO i = 0; i < this->nElems; i++ ) {
    for ( int j = 0; j < this->nNodesPerElem; j++ ) {
      file >> this->elems[ this->nNodesPerElem * i + j ];
    }
  }

  fixity.resize( this->nNodes );
  curvature.resize( this->nNodes );
  SCVT val;

  for ( LO i = 0; i < this->nNodes; ++i ) {
    file >> val;
    fixity.set( i, val );
  }

  for ( LO i = 0; i < this->nNodes; ++i ) {
    file >> val;
    // openFTL gives mean curvature, we need additive (in 3D 2*mean)
    curvature.set( i, 2.0 * val );
  }

  file.close( );

  //  initialize all auxiliary variables
  initArea( );
  initNormals( );
  initLocalCoordinates( );
  initEdges( );
  initCurl( );
  initNode2Elems( );
  initL2g( );

#ifdef VERBOSE
  std::cout << "done." << std::endl;
#endif

  return true;
}

template<class LO, class SC>
bool SurfaceMesh3D<LO, SC>::loadSmf(
    const string & meshFile
    ) {
#ifdef VERBOSE
  std::cout << "Reading file " << meshFile << " ... ";
#endif
  std::ifstream file( meshFile.c_str( ) );

  if ( !file.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return false;
  }

  this->dim = 3;
  this->nNodesPerElem = 3;

  // skip first two lines
  file.ignore( 50, '\n' );
  file.ignore( 50, '\n' );
  file >> this->nNodes;
  file >> this->nElems;
  this->nodes.resize( 3 * this->nNodes );

  for ( LO i = 0; i < this->nNodes; i++ ) {
    file >> this->nodes[3 * i];
    file >> this->nodes[3 * i + 1];
    file >> this->nodes[3 * i + 2];
  }

  this->elems.resize( this->nElems * this->nNodesPerElem );

  for ( LO i = 0; i < this->nElems; i++ ) {
    for ( int j = 0; j < this->nNodesPerElem; j++ ) {
      file >> this->elems[ this->nNodesPerElem * i + j ];
    }
  }

  file.close( );

  //  initialize all auxiliary variables
  initArea( );
  initNormals( );
  initLocalCoordinates( );
  initEdges( );
  initCurl( );
  initNode2Elems( );
  initL2g( );

#ifdef VERBOSE
  std::cout << "done." << std::endl;
#endif

  return true;
}

template<class LO, class SC>
bool SurfaceMesh3D<LO, SC>::loadFromNetgen(
  const string & meshFile
  ) {

  return false;
}

template<class LO, class SC>
bool SurfaceMesh3D<LO, SC>::loadFromParaview(
  const string & meshFile
  ) {

  return false;
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::printInfo( ) {
  std::cout << "Mesh type: SurfaceMesh3D, ";
  std::cout << "dim: " << this->dim;
  std::cout << ", nodes: " << this->nNodes;
  std::cout << ", elements: " << this->nElems;
  std::cout << ", edges: " << this->nEdges;
  std::cout << ", h: " << this->getDiscretizationParameter( );
  std::cout << "." << std::endl;
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::printParaviewVtk(
  const string & meshFile
  ) {

  std::ofstream file( meshFile.c_str( ) );

  file.setf( std::ios::showpoint | std::ios::scientific );
  file.precision( 6 );

#ifdef VERBOSE
  std::cout << "Printing '" << meshFile << "' ... ";
#endif

  if ( !file.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }
  file << "# vtk DataFile Version 3.0\n";
  file << "3D scalar data\n";
  file << "ASCII\n\n";
  file << "DATASET UNSTRUCTURED_GRID\n";
  file << "POINTS " << this->nNodes << " float\n";

  for ( LO i = 0; i < this->nNodes; i++ ) {
    file << this->nodes[ 3 * i ] << " " <<
      this->nodes[ 3 * i + 1 ] << " " <<
      this->nodes[ 3 * i + 2 ] << std::endl;
  }

  file << std::endl << "CELLS " << this->nElems << " " << 4 * this->nElems
    << std::endl;

  for ( LO i = 0; i < this->nElems; i++ ) {
    file << "3 " <<
      this->elems[ 3 * i ] << " " <<
      this->elems[ 3 * i + 1 ] << " " <<
      this->elems[ 3 * i + 2 ] << std::endl;
  }

  file << std::endl << "CELL_TYPES " << this->nElems << std::endl;

  for ( LO i = 0; i < this->nElems; i++ ) {
    file << "5" << std::endl;
  }

#ifdef VERBOSE
  std::cout << "done.\n";
#endif

  file.close( );
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::printParaviewVtu(
  const string& meshFile,
  int nNodal,
  string* nodeNames,
  Vector< LO, SC >** nodalData,
  int nElem, string* elemNames,
  Vector< LO, SC >** elemData
  ) {

#ifdef VERBOSE
  std::cout << "Printing '" << meshFile << "' ... ";
  std::cout.flush( );
#endif

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 6 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"
    << std::endl;
  file_vtu << "  <UnstructuredGrid>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << this->nNodes <<
    "\" NumberOfCells=\"" << this->nElems << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\""
    " NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  for ( LO i = 0; i < this->nNodes; i++ ) {
    file_vtu << "          "
      << this->nodes[ 3 * i ] << " "
      << this->nodes[ 3 * i + 1 ] << " "
      << this->nodes[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Cells>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\""
    " format=\"ascii\">" << std::endl;

  for ( LO i = 0; i < this->nElems; i++ ) {
    file_vtu << "          "
      << this->elems[ 3 * i ] << " "
      << this->elems[ 3 * i + 1 ] << " "
      << this->elems[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"offsets\""
    " format=\"ascii\">" << std::endl;

  for ( LO offset = 1; offset <= this->nElems; offset++ ) {
    file_vtu << "          " << offset * 3 << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"UInt8\" Name=\"types\""
    " format=\"ascii\">" << std::endl;
  for ( LO i = 1; i <= this->nElems; i++ ) {
    file_vtu << "          5" << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Cells>" << std::endl;

  if ( nNodal > 0 ) {
    string header = nodeNames[ 0 ];
    for ( int j = 1; j < nNodal; j++ ) {
      header += "," + nodeNames[ j ];
    }
    file_vtu << "      <PointData Scalars=\"" + header + "\">" << std::endl;
    for ( int j = 0; j < nNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + nodeNames[ j ] + "\" format=\"ascii\">" << std::endl;
      for ( LO i = 0; i < this->nNodes; i++ ) {
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
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + elemNames[ j ] + "\" format=\"ascii\">" << std::endl;
      for ( LO i = 0; i < this->nElems; i++ ) {
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

#ifdef VERBOSE
  std::cout << "done." << std::endl;
#endif

}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::printParaviewVtu(
  const string & meshFile,
  const std::vector< string > * nodeNames,
  const std::vector< Vector< LO, SC >* > * nodalData,
  const std::vector< string > * elemNames,
  const std::vector< Vector< LO, SC >* > * elemData,
  const std::vector< string > * nodeVNames,
  const std::vector< Vector< LO, SC >* > * nodalVData,
  const std::vector< string > * elemVNames,
  const std::vector< Vector< LO, SC >* > * elemVData
  ) const {

#ifdef VERBOSE
  std::cout << "Printing '" << meshFile << "' ... ";
  std::cout.flush( );
#endif

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 6 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"
    << std::endl;
  file_vtu << "  <UnstructuredGrid>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << this->nNodes <<
    "\" NumberOfCells=\"" << this->nElems << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\""
    " NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  for ( LO i = 0; i < this->nNodes; i++ ) {
    file_vtu << "          "
      << this->nodes[ 3 * i ] << " "
      << this->nodes[ 3 * i + 1 ] << " "
      << this->nodes[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Cells>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\""
    " format=\"ascii\">" << std::endl;

  for ( LO i = 0; i < this->nElems; i++ ) {
    file_vtu << "          "
      << this->elems[ 3 * i ] << " "
      << this->elems[ 3 * i + 1 ] << " "
      << this->elems[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu <<
    "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
    << std::endl;

  for ( LO offset = 1; offset <= this->nElems; offset++ ) {
    file_vtu << "          " << offset * 3 << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu <<
    "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
    << std::endl;
  for ( LO i = 1; i <= this->nElems; i++ ) {
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
      file_vtu << "Scalars=\"" + header + "\"";
    }

    if ( nVNodal > 0 ) {
      vheader = ( *nodeVNames )[ 0 ];
      for ( int j = 1; j < nVNodal; j++ ) {
        vheader += "," + ( *nodeVNames )[ j ];
      }
      file_vtu << " Vectors=\"" + vheader + "\"";
    }

    file_vtu << ">" << std::endl;
    for ( int j = 0; j < nNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *nodeNames )[ j ] + "\" format=\"ascii\">" << std::endl;
      for ( LO i = 0; i < this->nNodes; i++ ) {
        file_vtu << "          " << ( *nodalData )[ j ]->get( i ) << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    for ( int j = 0; j < nVNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *nodeVNames )[ j ] + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( LO i = 0; i < 3 * this->nNodes; i++ ) {
        file_vtu << "          " << ( *nodalVData )[ j ]->get( i ) << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }
    file_vtu << "      </PointData>" << std::endl;
  }

  int nElem = 0;
  if ( elemData ) nElem = elemData->size( );
  int nVElem = 0;
  if ( elemVData ) nVElem = elemVData->size( );

  if ( nElem > 0 || nVElem > 0 ) {
    file_vtu << "      <CellData ";

    if ( nElem > 0 ) {
      header = ( *elemNames )[ 0 ];
      for ( int j = 1; j < nElem; j++ ) {
        header += "," + ( *elemNames )[ j ];
      }
      file_vtu << "Scalars=\"" + header + "\"";
    }

    if ( nVElem > 0 ) {
      vheader = ( *elemVNames )[ 0 ];
      for ( int j = 1; j < nVElem; j++ ) {
        vheader += "," + ( *elemVNames )[ j ];
      }
      file_vtu << " Vectors=\"" + vheader + "\"";
    }

    file_vtu << ">" << std::endl;
    for ( int j = 0; j < nElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *elemNames )[ j ] + "\" format=\"ascii\">" << std::endl;
      for ( LO i = 0; i < this->nElems; i++ ) {
        file_vtu << "          " << ( *elemData )[ j ]->get( i ) << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    for ( int j = 0; j < nVElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *elemVNames )[ j ] + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( LO i = 0; i < 3 * this->nElems; i++ ) {
        file_vtu << "          " << ( *elemVData )[ j ]->get( i ) << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }
    file_vtu << "      </CellData>" << std::endl;
  }

  file_vtu << "    </Piece>" << std::endl;
  file_vtu << "  </UnstructuredGrid>" << std::endl;
  file_vtu << "</VTKFile>" << std::endl;
  file_vtu.close( );

#ifdef VERBOSE
  std::cout << "done." << std::endl;
#endif

}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::print(
  const string & meshFile
  ) {

#ifdef VERBOSE
  std::cout << "Printing  '" << meshFile << "' ... ";
#endif

  std::ofstream file( meshFile.c_str( ) );

  file.setf( std::ios::showpoint | std::ios::scientific );
  file.precision( 6 );

  if ( !file.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file << 3 << std::endl << 3 << std::endl << std::endl;
  file << this->nNodes << std::endl;

  for ( LO i = 0; i < this->nNodes; i++ ) {
    file << this->nodes[ 3 * i ] << " "
      << this->nodes[ 3 * i + 1 ] << " "
      << this->nodes[ 3 * i + 2 ] << std::endl;
  }

  file << std::endl << this->nElems << std::endl;

  for ( LO i = 0; i < this->nElems; i++ ) {
    file << this->elems[ 3 * i ] << " "
      << this->elems[ 3 * i + 1 ] << " "
      << this->elems[ 3 * i + 2 ] << std::endl;
  }

  file.close( );

#ifdef VERBOSE
  std::cout << "done." << std::endl;
#endif

}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::printSmf(
  const string & meshFile
  ) {

#ifdef VERBOSE
  std::cout << "Printing '" << meshFile << "' ... ";
#endif

  std::ofstream file( meshFile.c_str( ) );

  file.setf( std::ios::showpoint | std::ios::scientific );
  file.precision( 6 );

  if ( !file.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file << "! elementShape triangle" << std::endl << "! elementNumPoints 3"
    << std::endl << std::endl;
  file << this->nNodes << " " << this->nElems << std::endl;

  for ( LO i = 0; i < this->nNodes; i++ ) {
    file << this->nodes[ 3 * i ] << " "
      << this->nodes[ 3 * i + 1 ] << " "
      << this->nodes[ 3 * i + 2 ] << std::endl;
  }

  for ( LO i = 0; i < this->nElems; i++ ) {
    file << this->elems[ 3 * i ] << " "
      << this->elems[ 3 * i + 1 ] << " "
      << this->elems[ 3 * i + 2 ] << std::endl;
  }

  file.close( );

#ifdef VERBOSE
  std::cout << "done." << std::endl;
#endif

}

template<class LO, class SC>
bool SurfaceMesh3D<LO, SC>::isInside(
  const SCVT * x
  ) const {

  Vector< LO, SC > density( this->nNodes );
  density.setAll( -1.0 );

  Vector< LO, SC > result( 1 );
  // todo: is the copy necessary?
  SurfaceMesh3D< LO, SC > copyOfThis( *this );
  BESpace< LO, SC > bespace( &copyOfThis, p0, p0 );
  PotentialsLaplace< LO, SC > pot( &bespace, &density );
  pot.doubleLayerPotential( x, 1, result );

  if ( result.get( 0 ) > 1.0 - EPS ) return true;
  return false;
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::isInside(
  const SCVT * x,
  LO nPoints,
  bool * ret
  ) const {

  if ( typeid ( SC ) == typeid ( float ) ||
    typeid ( SC ) == typeid ( double ) ) {

    Vector< LO, SCVT > density( this->nNodes );
    density.setAll( -1.0 );
    Vector< LO, SCVT > result( nPoints );
    // todo: is the copy necessary?
    SurfaceMesh3D< LO, SCVT > copyOfThis;
    this->copyToReal( copyOfThis );
    BESpace< LO, SCVT > bespace( &copyOfThis, p0, p0 );

    PotentialsLaplace< LO, SCVT > potl( &bespace, &density );
    potl.doubleLayerPotential( x, nPoints, result );

    for ( LO i = 0; i < nPoints; ++i ) {
      if ( std::real( result.get( i ) ) > 0.9 ) {
        ret[ i ] = true;
      } else {
        ret[ i ] = false;
      }
    }

  } else if ( typeid ( SC ) == typeid ( std::complex< float > ) ||
    typeid ( SC ) == typeid ( std::complex< double > ) ) {
    Vector< LO, std::complex< SCVT > > density( this->nNodes );
    density.setAll( -1.0 );
    Vector< LO, std::complex< SCVT > > result( nPoints );
    // todo: is the copy necessary?
    SurfaceMesh3D< LO, std::complex< SCVT > > copyOfThis;
    this->copyToComplex( copyOfThis );
    BESpace< LO, std::complex< SCVT > > bespace( &copyOfThis, p0, p0 );

    PotentialsHelmholtz< LO, std::complex< SCVT > >
      poth( &bespace, 0.0, &density );
    poth.doubleLayerPotential( x, nPoints, result );

    for ( LO i = 0; i < nPoints; ++i ) {
      if ( std::real( result.get( i ) ) > 0.9 ) { //1.0 - EPS ) {
        ret[ i ] = true;
      } else {
        ret[ i ] = false;
      }
    }
  }
}

// todo: grid should be const but getNodes returns non-const pointer

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::trimEvaluationGrid(
  SurfaceMesh3D<LO, SC> & grid,
  SurfaceMesh3D<LO, SC> & gridTrimmed,
  bool trimInterior
  ) const {

  LO nNodes = grid.getNNodes( );
  LO nElems = grid.getNElements( );

  const SCVT * nodes = grid.getNodes( )->data( );
  bool * interiorNodes = new bool[ nNodes ];
  std::cout << "TEST" << std::endl;
  this->isInside( nodes, nNodes, interiorNodes );

  bool * keepNode = new bool[ nNodes ];
  LO nNodesNew = 0;
  bool * keepElem = new bool[ nElems ];
  LO nElemsNew = 0;
  LO * nodesMap = new LO[ nNodes ];

  // finding out which nodes to keep
  for ( LO i = 0; i < nNodes; ++i ) {
    if ( trimInterior == interiorNodes[ i ] ) {
      keepNode[ i ] = false;
      nodesMap[ i ] = 0; // not used
    } else {
      keepNode[ i ] = true;
      nodesMap[ i ] = nNodesNew++;
    }
  }

  // finding out which elements to keep
  LO elem[ 3 ];
  bool keep;
  for ( LO i = 0; i < nElems; ++i ) {
    keep = true;
    grid.getElement( i, elem );

    for ( int j = 0; j < 3; ++j ) {
      if ( keepNode[ elem[ j ] ] == false ) {
        keep = false;
        break;
      }
    }

    if ( keep ) ++nElemsNew;
    keepElem[ i ] = keep;
  }

  gridTrimmed.dim = 3;
  gridTrimmed.nNodesPerElem = 3;
  gridTrimmed.nNodes = nNodesNew;
  gridTrimmed.nElems = nElemsNew;

  gridTrimmed.nodes.clear( );
  gridTrimmed.nodes.reserve( gridTrimmed.dim * gridTrimmed.nNodes );
  gridTrimmed.elems.clear( );
  gridTrimmed.elems.reserve( gridTrimmed.nNodesPerElem * gridTrimmed.nElems );

  for ( LO i = 0; i < nNodes; ++i ) {
    if ( keepNode[ i ] ) {
      gridTrimmed.nodes.push_back( grid.nodes[ 3 * i ] );
      gridTrimmed.nodes.push_back( grid.nodes[ 3 * i + 1 ] );
      gridTrimmed.nodes.push_back( grid.nodes[ 3 * i + 2 ] );
    }
  }

  for ( LO i = 0; i < nElems; ++i ) {
    if ( keepElem[ i ] ) {
      grid.getElement( i, elem );
      gridTrimmed.elems.push_back( nodesMap[ elem[ 0 ] ] );
      gridTrimmed.elems.push_back( nodesMap[ elem[ 1 ] ] );
      gridTrimmed.elems.push_back( nodesMap[ elem[ 2 ] ] );
    }
  }

  gridTrimmed.auxCurl = nullptr;
  gridTrimmed.initArea( );
  gridTrimmed.initNormals( );
  gridTrimmed.initLocalCoordinates( );
  gridTrimmed.initEdges( );
  gridTrimmed.initCurl( );
  gridTrimmed.initL2g( );

  delete [] keepNode;
  delete [] keepElem;
  delete [] nodesMap;
  delete [] interiorNodes;
}

// todo: if node is not referenced by any interior and exterior element, it 
//       should be deleted (only important if only one layer of elements is
//       in/out )

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::trimEvaluationGrid(
  SurfaceMesh3D<LO, SC> & grid,
  SurfaceMesh3D<LO, SC> & interior,
  SurfaceMesh3D<LO, SC> & exterior
  ) const {

  LO nNodes = grid.getNNodes( );
  LO nElems = grid.getNElements( );

  const SCVT * nodes = grid.getNodes( )->data( );
  bool * interiorNodes = new bool[ nNodes ];
  this->isInside( nodes, nNodes, interiorNodes );

  LO nNodesNewIn = 0;
  LO nNodesNewEx = 0;
  bool * interiorElem = new bool[ nElems ];
  bool * exteriorElem = new bool[ nElems ];
  LO nElemsNewIn = 0;
  LO nElemsNewEx = 0;
  LO * nodesMapIn = new LO[ nNodes ];
  LO * nodesMapEx = new LO[ nNodes ];

  // finding out which nodes to keep
  for ( LO i = 0; i < nNodes; ++i ) {
    if ( interiorNodes[i] ) {
      nodesMapIn[i] = nNodesNewIn++;
      nodesMapEx[i] = 0;
    } else {
      nodesMapEx[i] = nNodesNewEx++;
      nodesMapIn[i] = 0;
    }
  }

  // finding out which elements lie in/outside
  LO elem[ 3 ];
  bool in, out;
  for ( LO i = 0; i < nElems; ++i ) {
    in = true;
    out = true;
    grid.getElement( i, elem );

    for ( int j = 0; j < 3; ++j ) {
      if ( interiorNodes[ elem[ j ] ] == true ) {
        out = false;
      }
      if ( interiorNodes[ elem[ j ] ] == false ) {
        in = false;
      }
    }

    if ( in ) ++nElemsNewIn;
    if ( out ) ++nElemsNewEx;
    interiorElem[ i ] = in;
    exteriorElem[ i ] = out;
  }

  interior.dim = 3;
  interior.nNodesPerElem = 3;
  interior.nNodes = nNodesNewIn;
  interior.nElems = nElemsNewIn;

  exterior.dim = 3;
  exterior.nNodesPerElem = 3;
  exterior.nNodes = nNodesNewEx;
  exterior.nElems = nElemsNewEx;

  interior.nodes.clear( );
  interior.nodes.reserve( interior.dim * interior.nNodes );
  interior.elems.clear( );
  interior.elems.reserve( interior.nNodesPerElem * interior.nElems );

  exterior.nodes.clear( );
  exterior.nodes.reserve( exterior.dim * exterior.nNodes );
  exterior.elems.clear( );
  exterior.elems.reserve( exterior.nNodesPerElem * exterior.nElems );

  for ( LO i = 0; i < nNodes; ++i ) {
    if ( interiorNodes[ i ] ) {
      interior.nodes.push_back( grid.nodes[ 3 * i ] );
      interior.nodes.push_back( grid.nodes[ 3 * i + 1 ] );
      interior.nodes.push_back( grid.nodes[ 3 * i + 2 ] );
    } else {
      exterior.nodes.push_back( grid.nodes[ 3 * i ] );
      exterior.nodes.push_back( grid.nodes[ 3 * i + 1 ] );
      exterior.nodes.push_back( grid.nodes[ 3 * i + 2 ] );
    }
  }

  for ( LO i = 0; i < nElems; ++i ) {
    if ( interiorElem[ i ] ) {
      grid.getElement( i, elem );
      interior.elems.push_back( nodesMapIn[ elem[ 0 ] ] );
      interior.elems.push_back( nodesMapIn[ elem[ 1 ] ] );
      interior.elems.push_back( nodesMapIn[ elem[ 2 ] ] );
    } else if ( exteriorElem[ i ] ) {
      grid.getElement( i, elem );
      exterior.elems.push_back( nodesMapEx[ elem[ 0 ] ] );
      exterior.elems.push_back( nodesMapEx[ elem[ 1 ] ] );
      exterior.elems.push_back( nodesMapEx[ elem[ 2 ] ] );
    }
  }

  interior.auxCurl = nullptr;
  interior.initArea( );
  interior.initNormals( );
  interior.initLocalCoordinates( );
  interior.initEdges( );
  interior.initCurl( );
  interior.initL2g( );

  exterior.auxCurl = nullptr;
  exterior.initArea( );
  exterior.initNormals( );
  exterior.initLocalCoordinates( );
  exterior.initEdges( );
  exterior.initCurl( );
  exterior.initL2g( );

  delete [] exteriorElem;
  delete [] interiorElem;
  delete [] nodesMapEx;
  delete [] nodesMapIn;
  delete [] interiorNodes;
}

//template<class LO, class SC>
//void SurfaceMesh3D<LO, SC>::initArea( ) {
//  this->totalArea = 0.0;
//  this->area.resize( this->nElems );
//  for ( LO i = 0; i < this->nElems; i++ ) {
//    this->area[i] = computeElemArea( i );
//    this->totalArea += this->area[i];
//  }
//}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::initArea( ) {

  SCVT total = 0.0;
  this->area.resize( this->nElems );

#pragma omp parallel for reduction( + : total )
  for ( LO i = 0; i < this->nElems; i++ ) {
    this->area[ i ] = this->computeElemArea( i );
    total += this->area[ i ];
  }

  this->totalArea = total;
}

//template<class LO, class SC>
//void SurfaceMesh3D<LO, SC>::initNormals( ) {
//  this->normals.resize( 3 * this->nElems );
//  Vector<LO, SCVT> edge1( 3 );
//  Vector<LO, SCVT> edge2( 3 );
//  Vector<LO, SCVT> normal( 3 );
//
//  SCVT x1[3], x2[3], x3[3];
//
//  // compute outer normal of each element using cross product
//  for ( LO i = 0; i < this->nElems; i++ ) {
//    this->getNodes( i, x1, x2, x3 );
//    edge1.set( 0, x2[0] - x1[0] );
//    edge1.set( 1, x2[1] - x1[1] );
//    edge1.set( 2, x2[2] - x1[2] );
//    edge2.set( 0, x3[0] - x1[0] );
//    edge2.set( 1, x3[1] - x1[1] );
//    edge2.set( 2, x3[2] - x1[2] );
//    edge1.cross( edge2, normal );
//
//    SCVT norm = normal.norm2( );
//
//    this->normals[3 * i] = normal.get( 0 ) / norm;
//    this->normals[3 * i + 1] = normal.get( 1 ) / norm;
//    this->normals[3 * i + 2] = normal.get( 2 ) / norm;
//
//  }
//}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::initNormals( ) {

  this->normals.resize( 3 * this->nElems );

#pragma omp parallel
  {
    Vector<LO, SCVT> edge1( 3 );
    Vector<LO, SCVT> edge2( 3 );
    Vector<LO, SCVT> normal( 3 );
    SCVT x1[3], x2[3], x3[3];

    // compute outer normal of each element using cross product
#pragma omp for
    for ( LO i = 0; i < this->nElems; ++i ) {
      this->getNodes( i, x1, x2, x3 );
      edge1.set( 0, x2[0] - x1[0] );
      edge1.set( 1, x2[1] - x1[1] );
      edge1.set( 2, x2[2] - x1[2] );
      edge2.set( 0, x3[0] - x1[0] );
      edge2.set( 1, x3[1] - x1[1] );
      edge2.set( 2, x3[2] - x1[2] );
      edge1.cross( edge2, normal );

      SCVT norm = normal.norm2( );

      this->normals[ 3 * i ] = normal.get( 0 ) / norm;
      this->normals[ 3 * i + 1 ] = normal.get( 1 ) / norm;
      this->normals[ 3 * i + 2 ] = normal.get( 2 ) / norm;
    }
  }
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::initL2g( ) {
  this->l2gElems.clear( );
  this->l2gElems.reserve( this->nElems );
  for ( LO i = 0; i < this->nElems; ++i ) {
    l2gElems.push_back( i );
  }

  this->l2gNodes.clear( );
  this->l2gNodes.reserve( this->nNodes );
  for ( LO i = 0; i < this->nNodes; ++i ) {
    l2gNodes.push_back( i );
  }
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::flipNormals( ) {

  LO aux;
  typename std::vector< LO >::iterator it;

#pragma omp parallel for private( aux )
  for ( it = this->elems.begin( ); it < this->elems.end( ); it += 3 ) {
    aux = *( it + 1 );
    *( it + 1 ) = *( it + 2 );
    *( it + 2 ) = aux;
  }

  initNormals( );
  initLocalCoordinates( );
  initEdges( );
  initCurl( );
  initNode2Elems( );
  initL2g( );
}

//template<class LO, class SC>
//void SurfaceMesh3D<LO, SC>::initLocalCoordinates( ) {
//  this->r1.resize( 9 * this->nElems );
//  this->r2.resize( 9 * this->nElems );
//  this->alpha1.resize( 3 * this->nElems );
//  this->alpha2.resize( 3 * this->nElems );
//  this->stau.resize( 3 * this->nElems );
//  SCVT x1[ 3 ], x2[ 3 ], x3[ 3 ];
//
//  for ( LO i = 0; i < this->nElems; i++ ) {
//    this->getNodes( i, x1, x2, x3 );
//    this->getLocalCoordinates( x1, x2, x3, &r1[ i * 9 ], &r2[ i * 9 ], &stau[ i * 3 ], &alpha1[ i * 3 ], &alpha2[ i * 3 ] );
//    this->getLocalCoordinates( x2, x3, x1, &r1[ i * 9 + 3 ], &r2[ i * 9 + 3 ], &stau[ i * 3 + 1 ], &alpha1[ i * 3 + 1 ], &alpha2[ i * 3 + 1 ] );
//    this->getLocalCoordinates( x3, x1, x2, &r1[ i * 9 + 6 ], &r2[ i * 9 + 6 ], &stau[ i * 3 + 2 ], &alpha1[ i * 3 + 2 ], &alpha2[ i * 3 + 2 ] );
//  }
//}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::initLocalCoordinates( ) {

  this->r1.resize( 9 * this->nElems );
  this->r2.resize( 9 * this->nElems );
  this->alpha1.resize( 3 * this->nElems );
  this->alpha2.resize( 3 * this->nElems );
  this->stau.resize( 3 * this->nElems );
  SCVT x1[ 3 ], x2[ 3 ], x3[ 3 ];

#pragma omp parallel for private( x1, x2, x3 )
  for ( LO i = 0; i < this->nElems; ++i ) {
    this->getNodes( i, x1, x2, x3 );
    this->getLocalCoordinates( x1, x2, x3, &r1[ i * 9 ], &r2[ i * 9 ],
      &stau[ i * 3 ], &alpha1[ i * 3 ], &alpha2[ i * 3 ] );
    this->getLocalCoordinates( x2, x3, x1, &r1[ i * 9 + 3 ], &r2[ i * 9 + 3 ],
      &stau[ i * 3 + 1 ], &alpha1[ i * 3 + 1 ], &alpha2[ i * 3 + 1 ] );
    this->getLocalCoordinates( x3, x1, x2, &r1[ i * 9 + 6 ], &r2[ i * 9 + 6 ],
      &stau[ i * 3 + 2 ], &alpha1[ i * 3 + 2 ], &alpha2[ i * 3 + 2 ] );
  }
}

//template<class LO, class SC>
//void SurfaceMesh3D<LO, SC>::initEdges( ) {
//
//  // allocate class variables
//  this->elem2edges.clear( );
//  this->elem2edges.resize( 3 * this->nElems );
//
//  std::vector< bool > edgeNotMapped( 3 * this->nElems, true );
//
//  // allocate aux. variables
//  LO element[ 3 ];
//  std::vector< std::pair< LO, LO > > localEdges;
//
//  // resize the container to appropriate size 
//  // (Euler formula - only valid for closed meshes)
//  localEdges.reserve( this->nElems + this->nNodes - 2 );
//
//  // iterate through elements and insert edges with index of starting point 
//  // lower than ending point
//  for ( LO i = 0; i < this->nElems; ++i ) {
//
//    this->getElement( i, element );
//
//    if ( element[ 0 ] < element[ 1 ] ) {
//      localEdges.push_back( std::pair< LO, LO >( element[ 0 ], element[ 1 ] ) );
//      this->elem2edges[ 3 * i ] = localEdges.size( ) - 1;
//      edgeNotMapped[ 3 * i ] = false;
//    }
//    if ( element[ 1 ] < element[ 2 ] ) {
//      localEdges.push_back( std::pair< LO, LO >( element[ 1 ], element[ 2 ] ) );
//      this->elem2edges[ 3 * i + 1 ] = localEdges.size( ) - 1;
//      edgeNotMapped[ 3 * i + 1 ] = false;
//    }
//    if ( element[ 2 ] < element[ 0 ] ) {
//      localEdges.push_back( std::pair< LO, LO >( element[ 2 ], element[ 0 ] ) );
//      this->elem2edges[ 3 * i + 2 ] = localEdges.size( ) - 1;
//      edgeNotMapped[ 3 * i + 2 ] = false;
//    }
//  }
//
//  typedef typename std::vector< std::pair< LO, LO > >::iterator pairIterator;
//  std::pair< LO, LO > pair;
//  pairIterator edgePosition;
//
//  // iterate through elements, find its edges in vector edges and 
//  // fill the mapping vector from element to edges
//  for ( LO i = 0; i < this->nElems; ++i ) {
//
//    this->getElement( i, element );
//
//    if ( edgeNotMapped[ 3 * i ] ) {
//      pair.first = element[ 1 ];
//      pair.second = element[ 0 ];
//
//      edgePosition = std::find < pairIterator, std::pair<LO, LO> > (
//          localEdges.begin( ), localEdges.end( ), pair );
//      if ( edgePosition != localEdges.end( ) && *edgePosition == pair ) {
//        this->elem2edges[ 3 * i ] =
//            std::distance( localEdges.begin( ), edgePosition );
//      } else { // for mesh not being a boundary of a domain (square, etc.)
//        localEdges.push_back( pair );
//        this->elem2edges[ 3 * i ] = localEdges.size( ) - 1;
//      }
//    }
//
//    if ( edgeNotMapped[ 3 * i + 1 ] ) {
//      pair.first = element[ 2 ];
//      pair.second = element[ 1 ];
//
//      edgePosition = std::find < pairIterator, std::pair< LO, LO > > (
//          localEdges.begin( ), localEdges.end( ), pair );
//      if ( edgePosition != localEdges.end( ) && *edgePosition == pair ) {
//        this->elem2edges[ 3 * i + 1 ] =
//            std::distance( localEdges.begin( ), edgePosition );
//      } else {
//        localEdges.push_back( pair );
//        this->elem2edges[ 3 * i + 1 ] = localEdges.size( ) - 1;
//      }
//    }
//
//    if ( edgeNotMapped[ 3 * i + 2 ] ) {
//      pair.first = element[ 0 ];
//      pair.second = element[ 2 ];
//
//      edgePosition = std::find < pairIterator, std::pair< LO, LO > > (
//          localEdges.begin( ), localEdges.end( ), pair );
//      if ( edgePosition != localEdges.end( ) && *edgePosition == pair ) {
//        this->elem2edges[ 3 * i + 2 ] =
//            std::distance( localEdges.begin( ), edgePosition );
//      } else {
//        localEdges.push_back( pair );
//        this->elem2edges[ 3 * i + 2 ] = localEdges.size( ) - 1;
//      }
//    }
//  }
//
//  // finally convert edges from vector of std::pairs to 1D vector of LO
//  LO size = localEdges.size( );
//  this->edges.clear( );
//  this->edges.resize( 2 * size );
//
//#pragma omp parallel for
//  for ( LO i = 0; i < size; ++i ) {
//    this->edges[ 2 * i ] = localEdges[ i ].first;
//    this->edges[ 2 * i + 1 ] = localEdges[ i ].second;
//  }
//
//  this->nEdges = size;
//}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::initEdges( ) {

  // allocate class variables
  this->elem2edges.clear( );
  this->elem2edges.resize( 3 * this->nElems );

  // allocate aux. variables
  LO element[ 3 ];

  std::vector< std::vector< LO > > localEdges;
  localEdges.resize( this->nNodes );
  std::vector< std::pair< LO, LO > > elem2EdgesTemp;
  elem2EdgesTemp.resize( 3 * this->nElems );

  // iterate through elements and insert edges with index of starting point 
  // lower than ending point
  for ( LO i = 0; i < this->nElems; ++i ) {

    this->getElement( i, element );

    if ( element[ 0 ] < element[ 1 ] ) {
      localEdges[ element[ 0 ] ].push_back( element[ 1 ] );
    }
    if ( element[ 1 ] < element[ 2 ] ) {
      localEdges[ element[ 1 ] ].push_back( element[ 2 ] );
    }
    if ( element[ 2 ] < element[ 0 ] ) {
      localEdges[ element[ 2 ] ].push_back( element[ 0 ] );
    }
  }

  typedef typename std::vector< LO >::iterator it;
  it edgePosition;
  LO f, s;

  // iterate through elements, find its edges in vector edges and 
  // fill the mapping vector from element to edges
  for ( LO i = 0; i < this->nElems; ++i ) {

    this->getElement( i, element );

    if ( element[ 0 ] < element[ 1 ] ) {
      f = element[ 0 ];
      s = element[ 1 ];
    } else {
      f = element[ 1 ];
      s = element[ 0 ];
    }

    edgePosition = std::find < it, LO > (
      localEdges[ f ].begin( ), localEdges[ f ].end( ), s );
    if ( edgePosition != localEdges[ f ].end( ) && *edgePosition == s ) {
      elem2EdgesTemp[ 3 * i ] = std::pair< LO, LO >( f,
        edgePosition - localEdges[ f ].begin( ) );
    } else { // for mesh not being a boundary of a domain (square, etc.)
      localEdges[ f ].push_back( s );
      elem2EdgesTemp[ 3 * i ] =
        std::pair< LO, LO >( f, localEdges[ f ].size( ) - 1 );
    }

    if ( element[ 1 ] < element[ 2 ] ) {
      f = element[ 1 ];
      s = element[ 2 ];
    } else {
      f = element[ 2 ];
      s = element[ 1 ];
    }

    edgePosition = std::find < it, LO > (
      localEdges[ f ].begin( ), localEdges[ f ].end( ), s );
    if ( edgePosition != localEdges[ f ].end( ) && *edgePosition == s ) {
      elem2EdgesTemp[ 3 * i + 1 ] = std::pair< LO, LO >( f,
        edgePosition - localEdges[ f ].begin( ) );
    } else {
      localEdges[ f ].push_back( s );
      elem2EdgesTemp[ 3 * i + 1 ] =
        std::pair< LO, LO >( f, localEdges[ f ].size( ) - 1 );
    }

    if ( element[ 2 ] < element[ 0 ] ) {
      f = element[ 2 ];
      s = element[ 0 ];
    } else {
      f = element[ 0 ];
      s = element[ 2 ];
    }

    edgePosition = std::find < it, LO > (
      localEdges[ f ].begin( ), localEdges[ f ].end( ), s );
    if ( edgePosition != localEdges[ f ].end( ) && *edgePosition == s ) {
      elem2EdgesTemp[ 3 * i + 2 ] = std::pair< LO, LO >( f,
        edgePosition - localEdges[ f ].begin( ) );
    } else {
      localEdges[ f ].push_back( s );
      elem2EdgesTemp[ 3 * i + 2 ] =
        std::pair< LO, LO >( f, localEdges[ f ].size( ) - 1 );
    }

  }

  std::vector< LO > offsets( this->nNodes );
  LO offset = 0;
  for ( LO i = 0; i < this->nNodes; ++i ) {
    offsets[ i ] = offset;
    offset += localEdges[ i ].size( );
  }
  this->nEdges = offset;

  for ( LO i = 0; i < this->nElems; ++i ) {
    for ( int j = 0; j < 3; ++j ) {
      this->elem2edges[ 3 * i + j ] = elem2EdgesTemp[ 3 * i + j ].second +
        offsets[ elem2EdgesTemp[ 3 * i + j ].first ];
    }
  }

  // finally convert edges from vector of std::pairs to 1D vector of LO
  this->edges.clear( );
  this->edges.resize( 2 * this->nEdges );

#pragma omp parallel for
  for ( LO i = 0; i < this->nNodes; ++i ) {
    for ( LO j = 0; j < localEdges[ i ].size( ); ++j ) {
      this->edges[ 2 * ( offsets[ i ] + j ) ] = i;
      this->edges[ 2 * ( offsets[ i ] + j ) + 1 ] = localEdges[ i ][ j ];
    }
  }

}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::getLocalCoordinates(
  const SCVT * x1,
  const SCVT * x2,
  const SCVT * x3,
  SCVT * r1,
  SCVT * r2,
  SCVT * stau,
  SCVT * alpha1,
  SCVT * alpha2
  ) const {

  SCVT tk, t0;

  r2[0] = x3[0] - x2[0];
  r2[1] = x3[1] - x2[1];
  r2[2] = x3[2] - x2[2];

  tk = std::sqrt( DOT3( r2, r2 ) );

  r2[0] /= tk;
  r2[1] /= tk;
  r2[2] /= tk;

  r1[0] = x2[0] - x1[0];
  r1[1] = x2[1] - x1[1];
  r1[2] = x2[2] - x1[2];

  t0 = -DOT3( r1, r2 );

  r1[0] += t0 * r2[0];
  r1[1] += t0 * r2[1];
  r1[2] += t0 * r2[2];

  stau[0] = std::sqrt( DOT3( r1, r1 ) );

  r1[0] /= stau[0];
  r1[1] /= stau[0];
  r1[2] /= stau[0];

  alpha1[0] = -t0 / stau[0];
  alpha2[0] = ( tk - t0 ) / stau[0];
}

//template<class LO, class SC>
//typename bem4i::SurfaceMesh3D<LO, SC>::SCVT SurfaceMesh3D<LO, SC>::l2RelativeErrorConst( Vector<LO, SC> const & anal, Vector<LO, SC> const & res ) const {
//  SCVT ret = 0.0;
//  SCVT div = 0.0;
//  for ( LO i = 0; i < this->nElems; i++ ) {
//    ret += pow( std::abs( anal.get( i ) - res.get( i ) ), 2 ) * this->getElemArea( i );
//    div += pow( std::abs( anal.get( i ) ), 2 ) * this->getElemArea( i );
//  }
//
//  ret = std::sqrt( ret / div );
//  return ret;
//}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::getQuadratureNodes(
  const SCVT * x1,
  const SCVT * x2,
  const SCVT * x3,
  SCVT * nodes,
  int order
  ) const {

  int numPoints = quadSizes[ order ];
  for ( int i = 0; i < numPoints; i++ ) {
    nodes[ i * 3 ] = x1[ 0 ] + ( x2[ 0 ] - x1[ 0 ] )
      * quadPoints[ order ][ i * 2 ] + ( x3[ 0 ] - x1[ 0 ] )
      * quadPoints[ order ][ i * 2 + 1 ];
    nodes[ i * 3 + 1 ] = x1[ 1 ] + ( x2[ 1 ] - x1[ 1 ] )
      * quadPoints[ order ][ i * 2 ] + ( x3[ 1 ] - x1[ 1 ] )
      * quadPoints[ order ][ i * 2 + 1 ];
    nodes[ i * 3 + 2 ] = x1[ 2 ] + ( x2[ 2 ] - x1[ 2 ] )
      * quadPoints[ order ][ i * 2 ] + ( x3[ 2 ] - x1[ 2 ] )
      * quadPoints[ order ][ i * 2 + 1 ];
  }
}

template<class LO, class SC>
typename bem4i::SurfaceMesh3D<LO, SC>::SCVT SurfaceMesh3D<LO, SC>::
l2RelativeErrorConst(
  Vector<LO, SC> const & res,
  int order,
  SC kappa
  ) const {

  SCVT ret = 0.0;
  SCVT div = 0.0;
  SCVT aux1, aux2;
  SC anal;
  SCVT x1[3], x2[3], x3[3], x[3], n[3];
  int size = quadSizes[ order ];
  SCVT* quad = new SCVT[ 3 * size ];
  for ( LO i = 0; i < this->nElems; i++ ) {
    this->getNodes( i, x1, x2, x3 );
    this->getQuadratureNodes( x1, x2, x3, quad, order );
    aux1 = 0.0;
    aux2 = 0.0;
    this->getNormal( i, n );

    for ( int j = 0; j < size; j++ ) {
      x[0] = quad[3 * j];
      x[1] = quad[3 * j + 1];
      x[2] = quad[3 * j + 2];
      ///*
      anal = std::exp( 2.0 * M_PI * x[1] ) * ( n[0] * std::cos( 2.0 * M_PI * x[2] )
        + 2.0 * M_PI * ( 1.0 + x[0] ) * n[1] * std::cos( 2.0 * M_PI * x[2] )
        - 2.0 * M_PI * ( 1.0 + x[0] ) * n[2] * std::sin( 2.0 * M_PI * x[2] ) );
      //*/
      //anal = n[ 0 ] + n[ 1 ] + n[ 2 ];
      aux1 += ( std::abs( anal - res.get( i ) )
        * std::abs( anal - res.get( i ) ) * quadWeights[ order ][ j ] );
      aux2 += ( std::abs( anal ) * std::abs( anal )
        * quadWeights[ order ][ j ] );
    }
    ret += aux1 * this->getElemArea( i );
    div += aux2 * this->getElemArea( i );

  }

  delete [] quad;

  ret = std::sqrt( ret / div );
  return ret;
}

template<>
double SurfaceMesh3D< int, std::complex< double > >::l2RelativeErrorConst(
  Vector< int, std::complex< double > > const & res,
  int order,
  std::complex< double > kappa
  ) const {

  double ret = 0.0;
  double div = 0.0;
  double aux1, aux2;
  std::complex<double> anal;
  std::complex<double> iUnit( 0.0, 1.0 );
  double x1[3], x2[3], x3[3], x[3], n[3];
  int size = quadSizes[ order ];
  double* quad = new double[ 3 * size ];

  std::complex<double> b = 2.0 * kappa;
  std::complex<double> a = -std::sqrt( 3.0 ) * kappa;

  for ( int i = 0; i < this->nElems; i++ ) {
    this->getNodes( i, x1, x2, x3 );
    this->getQuadratureNodes( x1, x2, x3, quad, order );
    aux1 = 0.0;
    aux2 = 0.0;
    this->getNormal( i, n );
    for ( int j = 0; j < size; j++ ) {
      x[0] = quad[3 * j];
      x[1] = quad[3 * j + 1];
      x[2] = quad[3 * j + 2];
      //anal = std::exp(2.0*sqrt(3.0)*x[1]+4.0*iUnit*x[2]) * 4.0 * x[0] 
      //   * ( 1.0 + 2.0*sqrt(3.0)*x[1] + 4.0*iUnit*x[2] ); // Helmholtz
      //anal = std::exp( 2.0 * std::sqrt( 3.0 ) * x[1] + 4.0 * iUnit * x[2] ) 
      //   * 4.0 * ( n[0] + 2.0 * std::sqrt( 3.0 ) * n[1] * x[0] 
      //   + 4.0 * iUnit * x[0] * n[2] ); //Helmholtz n

      anal = b * std::exp( -a * x[ 1 ] + iUnit * b * x[ 2 ] ) * ( n[ 0 ] -
        a * x[ 0 ] * n[ 1 ] + iUnit * b * x[ 0 ] * n[ 2 ] );
      //anal = x[0] * x[1] * x[2];
      aux1 += pow( std::abs( anal - res.get( i ) ), 2 )
        * quadWeights[ order ][ j ];
      aux2 += pow( std::abs( anal ), 2 ) * quadWeights[ order ][ j ];
    }
    ret += aux1 * this->getElemArea( i );
    div += aux2 * this->getElemArea( i );
  }

  delete [] quad;

  ret = std::sqrt( ret / div );
  return ret;
}

template<>
double SurfaceMesh3D< long, std::complex< double > >::l2RelativeErrorConst(
  Vector< long, std::complex< double > > const & res,
  int order,
  std::complex< double > kappa
  ) const {

  double ret = 0.0;
  double div = 0.0;
  double aux1, aux2;
  std::complex<double> anal;
  std::complex<double> iUnit( 0.0, 1.0 );
  double x1[3], x2[3], x3[3], x[3], n[3];
  int size = quadSizes[ order ];
  double* quad = new double[ 3 * size ];
  for ( long i = 0; i < this->nElems; i++ ) {
    this->getNodes( i, x1, x2, x3 );
    this->getQuadratureNodes( x1, x2, x3, quad, order );
    aux1 = 0.0;
    aux2 = 0.0;
    this->getNormal( i, n );
    for ( int j = 0; j < size; j++ ) {
      x[0] = quad[3 * j];
      x[1] = quad[3 * j + 1];
      x[2] = quad[3 * j + 2];
      //anal = std::exp(2.0*sqrt(3.0)*x[1]+4.0*iUnit*x[2]) * 4.0 * x[0] 
      //   * ( 1.0 + 2.0*sqrt(3.0)*x[1] + 4.0*iUnit*x[2] ); // Helmholtz
      anal = std::exp( 2.0 * std::sqrt( 3.0 ) * x[1] + 4.0 * iUnit * x[2] ) * 4.0
        * ( n[0] + 2.0 * std::sqrt( 3.0 ) * n[1] * x[0]
        + 4.0 * iUnit * x[0] * n[2] ); //Helmholtz n

      aux1 += pow( std::abs( anal - res.get( i ) ), 2 )
        * quadWeights[ order ][ j ];
      aux2 += pow( std::abs( anal ), 2 ) * quadWeights[ order ][ j ];
    }
    ret += aux1 * this->getElemArea( i );
    div += aux2 * this->getElemArea( i );
  }

  delete [] quad;

  ret = std::sqrt( ret / div );
  return ret;
}

template<>
float SurfaceMesh3D< int, std::complex< float > >::l2RelativeErrorConst(
  Vector< int, std::complex< float > > const & res,
  int order,
  std::complex< float > kappa
  ) const {

  float ret = 0.0;
  float div = 0.0;
  float aux1, aux2;
  float two = 2.0;
  float three = 3.0;
  float four = 4.0;
  std::complex<float> anal;
  std::complex<float> iUnit( 0.0, 1.0 );
  float x1[3], x2[3], x3[3], x[3], n[3];
  int size = quadSizes[ order ];
  float* quad = new float[ 3 * size ];
  for ( int i = 0; i < this->nElems; i++ ) {
    this->getNodes( i, x1, x2, x3 );
    this->getQuadratureNodes( x1, x2, x3, quad, order );
    aux1 = 0.0;
    aux2 = 0.0;
    this->getNormal( i, n );
    for ( int j = 0; j < size; j++ ) {
      x[0] = quad[3 * j];
      x[1] = quad[3 * j + 1];
      x[2] = quad[3 * j + 2];
      //anal = std::exp(2.0*sqrt(3.0)*x[1]+4.0*iUnit*x[2]) * 4.0 * x[0] 
      //   * ( 1.0 + 2.0 * std::sqrt( 3.0 ) * x[1] + 4.0 * iUnit * x[2] );
      anal = std::exp( two * std::sqrt( three ) * x[1] + four * iUnit * x[2] )
        * four * ( n[0] + two * std::sqrt( three ) * n[1] * x[0]
        + four * iUnit * x[0] * n[2] ); //Helmholtz n

      aux1 += pow( std::abs( anal - res.get( i ) ), 2 )
        * quadWeights[ order ][ j ];
      aux2 += pow( std::abs( anal ), 2 ) * quadWeights[ order ][ j ];
    }
    ret += aux1 * this->getElemArea( i );
    div += aux2 * this->getElemArea( i );
  }

  delete [] quad;

  ret = std::sqrt( ret / div );
  return ret;
}

//template<class LO, class SC>
//typename bem4i::SurfaceMesh3D<LO, SC>::SCVT SurfaceMesh3D<LO, SC>::l2RelativeErrorLin( Vector<LO, SC> const & anal, Vector<LO, SC> const & res ) const {
//  SCVT ret = 0.0;
//  SCVT div = 0.0;
//  SCVT aux1, aux2;
//  LO ind[ 3 ];
//
//  for ( LO i = 0; i < this->nElems; i++ ) {
//    this->getElement( i, ind );
//    aux1 = 0.0;
//    aux2 = 0.0;
//    for ( int j = 0; j < 7; j++ ) {
//      aux1 += quadWeights5[ j ] * pow( std::abs(
//          ( anal.get( ind[0] ) - res.get( ind[0] ) ) * ( 1.0 - quadPoints5[ 2 * j ] - quadPoints5[ 2 * j + 1 ] ) +
//          ( anal.get( ind[1] ) - res.get( ind[1] ) ) * quadPoints5[ 2 * j ] +
//          ( anal.get( ind[2] ) - res.get( ind[2] ) ) * quadPoints5[ 2 * j + 1 ]
//          ), 2 );
//      aux2 += quadWeights5[ j ] * pow( std::abs(
//          anal.get( ind[0] ) * ( 1.0 - quadPoints5[ 2 * j ] - quadPoints5[ 2 * j + 1 ] ) +
//          anal.get( ind[1] ) * quadPoints5[ 2 * j ] +
//          anal.get( ind[2] ) * quadPoints5[ 2 * j + 1 ]
//          ), 2 );
//    }
//    ret += aux1 * this->getElemArea( i );
//    div += aux2 * this->getElemArea( i );
//  }
//
//  ret = std::sqrt( ret / div );
//  return ret;
//}

template<class LO, class SC>
typename bem4i::SurfaceMesh3D<LO, SC>::SCVT SurfaceMesh3D<LO, SC>::
l2RelativeErrorLin(
  Vector< LO, SC > const & res,
  int order,
  SC kappa
  ) const {

  SCVT ret = 0.0;
  SCVT div = 0.0;
  SCVT aux1, aux2;
  LO ind[ 3 ];
  SCVT x1[3], x2[3], x3[3], x[3], n[3];
  int size = quadSizes[ order ];
  SCVT* quad = new SCVT[ 3 * size ];
  SC anal;

  for ( LO i = 0; i < this->nElems; i++ ) {
    this->getElement( i, ind );
    this->getNodes( i, x1, x2, x3 );
    this->getNormal( i, n );
    this->getQuadratureNodes( x1, x2, x3, quad, order );
    aux1 = 0.0;
    aux2 = 0.0;
    for ( int j = 0; j < size; j++ ) {
      x[0] = quad[3 * j];
      x[1] = quad[3 * j + 1];
      x[2] = quad[3 * j + 2];
      anal = ( 1.0 + x[0] ) * std::exp( 2.0 * M_PI * x[1] ) *
          std::cos( 2.0 * M_PI * x[2] ); // Laplace Dirichlet
      //anal = std::exp( 2.0 * M_PI * x[1] )
      //  * ( n[0] * std::cos( 2.0 * M_PI * x[2] )
      //  + 2.0 * M_PI * ( 1.0 + x[0] ) * n[1] * std::cos( 2.0 * M_PI * x[2] )
      //  - 2.0 * M_PI * ( 1.0 + x[0] ) * n[2] * std::sin( 2.0 * M_PI * x[2] ) ); // Laplace Neumann
      aux1 += ( (SCVT) quadWeights[ order ][ j ] ) * pow( std::abs( anal -
        ( res.get( ind[0] ) * ( (SCVT) ( 1.0 - quadPoints[ order ][ 2 * j ]
        - quadPoints[ order ][ 2 * j + 1 ] ) ) +
        res.get( ind[1] ) * ( (SCVT) quadPoints[ order ][ 2 * j ] ) +
        res.get( ind[2] ) * ( (SCVT) quadPoints[ order ][ 2 * j + 1 ] ) )
        ), 2 );
      aux2 += quadWeights[ order ][ j ] * pow( std::abs( anal ), 2 );
    }
    ret += aux1 * this->getElemArea( i );
    div += aux2 * this->getElemArea( i );
  }

  delete [] quad;

  ret = std::sqrt( ret / div );
  return ret;
}

template<>
double SurfaceMesh3D< int, std::complex< double > >::l2RelativeErrorLin(
  Vector< int, std::complex< double > > const & res,
  int order,
  std::complex< double > kappa
  ) const {

  double ret = 0.0;
  double div = 0.0;
  double aux1, aux2;
  int ind[ 3 ];
  double x1[3], x2[3], x3[3], x[3], n[3];
  int size = quadSizes[ order ];
  double* quad = new double[ 3 * size ];
  std::complex<double> anal;
  std::complex<double> iUnit( 0.0, 1.0 );

  std::complex<double> b = 2.0 * kappa;
  std::complex<double> a = -std::sqrt( 3.0 ) * kappa;

  for ( int i = 0; i < this->nElems; i++ ) {
    this->getElement( i, ind );
    this->getNodes( i, x1, x2, x3 );
    this->getQuadratureNodes( x1, x2, x3, quad, order );
    this->getNormal( i, n );
    aux1 = 0.0;
    aux2 = 0.0;
    for ( int j = 0; j < size; j++ ) {
      x[0] = quad[3 * j];
      x[1] = quad[3 * j + 1];
      x[2] = quad[3 * j + 2];
      anal = b * x[ 0 ] * std::exp( -a * x[ 1 ] + iUnit * b * x[ 2 ] ); // Helmholtz Dirichlet
      //anal = b * std::exp( -a * x[ 1 ] + iUnit * b * x[ 2 ] ) * ( n[ 0 ] -
      //  a * x[ 0 ] * n[ 1 ] + iUnit * b * x[ 0 ] * n[ 2 ] ); // Helmholtz Neumann

      aux1 += quadWeights[ order ][ j ] * pow( std::abs( anal -
        ( res.get( ind[0] ) * ( 1.0 - quadPoints[ order ][ 2 * j ]
        - quadPoints[ order ][ 2 * j + 1 ] ) +
        res.get( ind[1] ) * quadPoints[ order ][ 2 * j ] +
        res.get( ind[2] ) * quadPoints[ order ][ 2 * j + 1 ] )
        ), 2 );
      aux2 += quadWeights[ order ][ j ] * pow( std::abs( anal ), 2 );
    }
    ret += aux1 * this->getElemArea( i );
    div += aux2 * this->getElemArea( i );
  }

  delete [] quad;

  ret = std::sqrt( ret / div );
  return ret;
}

template<>
double SurfaceMesh3D< long, std::complex< double > >::l2RelativeErrorLin(
  Vector< long, std::complex< double > > const & res,
  int order,
  std::complex< double > kappa
  ) const {

  double ret = 0.0;
  double div = 0.0;
  double aux1, aux2;
  long ind[ 3 ];
  double x1[3], x2[3], x3[3], x[3], n[3];
  long size = quadSizes[ order ];
  double* quad = new double[ 3 * size ];
  std::complex<double> anal;
  std::complex<double> iUnit( 0.0, 1.0 );

  std::complex<double> b = 2.0 * kappa;
  std::complex<double> a = -std::sqrt( 3.0 ) * kappa;

  for ( int i = 0; i < this->nElems; i++ ) {
    this->getElement( i, ind );
    this->getNodes( i, x1, x2, x3 );
    this->getQuadratureNodes( x1, x2, x3, quad, order );
    this->getNormal( i, n );
    aux1 = 0.0;
    aux2 = 0.0;
    for ( long j = 0; j < size; j++ ) {
      x[0] = quad[3 * j];
      x[1] = quad[3 * j + 1];
      x[2] = quad[3 * j + 2];
      anal = b * x[ 0 ] * std::exp( -a * x[ 1 ] + iUnit * b * x[ 2 ] );
      //anal = b * std::exp( -a * x[ 1 ] + iUnit * b * x[ 2 ] ) * ( n[ 0 ] -
      //  a * x[ 0 ] * n[ 1 ] + iUnit * b * x[ 0 ] * n[ 2 ] );
      aux1 += quadWeights[ order ][ j ] * pow( std::abs( anal -
        ( res.get( ind[0] ) * ( 1.0 - quadPoints[ order ][ 2 * j ]
        - quadPoints[ order ][ 2 * j + 1 ] ) +
        res.get( ind[1] ) * quadPoints[ order ][ 2 * j ] +
        res.get( ind[2] ) * quadPoints[ order ][ 2 * j + 1 ] )
        ), 2 );
      aux2 += quadWeights[ order ][ j ] * pow( std::abs( anal ), 2 );
    }
    ret += aux1 * this->getElemArea( i );
    div += aux2 * this->getElemArea( i );
  }

  delete [] quad;

  ret = std::sqrt( ret / div );
  return ret;
}

template<>
float SurfaceMesh3D< int, std::complex< float > >::l2RelativeErrorLin(
  Vector< int, std::complex< float > > const & res,
  int order,
  std::complex< float > kappa
  ) const {

  float ret = 0.0;
  float div = 0.0;
  float aux1, aux2;
  int ind[ 3 ];
  float x1[3], x2[3], x3[3], x[3];
  int size = quadSizes[ order ];
  float* quad = new float[ 3 * size ];
  std::complex<float> anal;
  std::complex<float> iUnit( 0.0, 1.0 );

  std::complex<float> b = (float) 2.0 * kappa;
  std::complex<float> a = -std::sqrt( (float) 3.0 ) * kappa;

  for ( int i = 0; i < this->nElems; i++ ) {
    this->getElement( i, ind );
    this->getNodes( i, x1, x2, x3 );
    this->getQuadratureNodes( x1, x2, x3, quad, order );
    aux1 = 0.0;
    aux2 = 0.0;
    for ( int j = 0; j < size; j++ ) {
      x[0] = quad[3 * j];
      x[1] = quad[3 * j + 1];
      x[2] = quad[3 * j + 2];
      //anal = four * x[0] 
      //    * std::exp( two * std::sqrt( three ) * x[1] + four * iUnit * x[2] );
      anal = b * x[ 0 ] * std::exp( -a * x[ 1 ] + iUnit * b * x[ 2 ] );
      aux1 += quadWeights[ order ][ j ] * pow( std::abs( anal -
        ( res.get( ind[0] ) * ( (float) 1.0
        - ( (float) quadPoints[ order ][ 2 * j ] )
        - ( (float) quadPoints[ order ][ 2 * j + 1 ] ) ) +
        res.get( ind[1] ) * ( (float) quadPoints[ order ][ 2 * j ] ) +
        res.get( ind[2] ) * ( (float) quadPoints[ order ][ 2 * j + 1 ] ) )
        ), 2 );
      aux2 += quadWeights[ order ][ j ] * pow( std::abs( anal ), 2 );
    }
    ret += aux1 * this->getElemArea( i );
    div += aux2 * this->getElemArea( i );
  }

  delete [] quad;

  ret = std::sqrt( ret / div );
  return ret;
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::nestedDissection(
  Tree<BECluster<LO, SC>*, LO >& clusterTree,
  LO const maxElemsPerCluster
  ) {

  // NOTE: all dynamically created objects will be destroyed by BECluster
  //       destructor

  // temporary variable to store centroids of all elements
  centroids = new std::vector< SCVT >( this->dim * this->nElems );

  // create the root cluster
  BECluster<LO, SC> *root = new BECluster<LO, SC>( );

  root->elems = new std::vector< LO >( this->nElems );
  for ( LO i = 0; i < this->nElems; i++ ) {
    // the root cluster contains all elements
    ( *root->elems )[ i ] = i;
    // we will also need centroids of all elements
    ( *centroids )[ 3 * i ] = ( this->nodes[ 3 * this->elems[ i * 3 ] ] +
      this->nodes[ 3 * this->elems[ i * 3 + 1] ] +
      this->nodes[ 3 * this->elems[ i * 3 + 2] ] ) / 3.0;
    ( *centroids )[ 3 * i + 1 ] = ( this->nodes[ 3 * this->elems[ i * 3 ] + 1 ]
      + this->nodes[ 3 * this->elems[ i * 3 + 1] + 1 ]
      + this->nodes[ 3 * this->elems[ i * 3 + 2] + 1 ] ) / 3.0;
    ( *centroids )[ 3 * i + 2 ] = ( this->nodes[ 3 * this->elems[ i * 3 ] + 2 ]
      + this->nodes[ 3 * this->elems[ i * 3 + 1] + 2 ]
      + this->nodes[ 3 * this->elems[ i * 3 + 2] + 2 ] ) / 3.0;
  }

  // the root cluster contains all nodes
  root->nodes = new std::set< LO >( this->elems.begin( ), this->elems.end( ) );

  root->nelems = this->nElems;
  root->nnodes = this->nNodes;

  if ( this->area.size( ) != this->nElems ) {
    initArea( );
  }

  clusterTree.createRoot( root );

  doNestedDissection( clusterTree, maxElemsPerCluster );

  delete centroids;

}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::doNestedDissection(
  Tree<BECluster<LO, SC>*, LO>& clusterTree,
  LO const maxElemsPerCluster
  ) const {

  // for algorithm see Steinbach, Rjasanow: Fast BEM, p. 109

  BECluster<LO, SC>* currentCluster = clusterTree.get( );
  std::vector<LO>* globalElems = currentCluster->elems;
  currentCluster->centroid = new SCVT[3];
  currentCluster->centroid[0] = 0.0;
  currentCluster->centroid[1] = 0.0;
  currentCluster->centroid[2] = 0.0;
  SCVT totArea = 0.0;

  // compute the centroid of the current cluster
  for ( LO i = 0; i < currentCluster->nelems; i++ ) {
    currentCluster->centroid[ 0 ] +=
      ( *centroids )[ 3 * ( *globalElems )[ i ] ]
      * this->area[ ( *globalElems )[ i ] ];
    currentCluster->centroid[ 1 ] +=
      ( *centroids )[ 3 * ( *globalElems )[ i ] + 1 ]
      * this->area[ ( *globalElems )[ i ] ];
    currentCluster->centroid[ 2 ] +=
      ( *centroids )[ 3 * ( *globalElems )[ i ] + 2 ]
      * this->area[ ( *globalElems )[ i ] ];
    totArea += this->area[ ( *globalElems )[ i ] ];
  }
  currentCluster->centroid[ 0 ] /= totArea;
  currentCluster->centroid[ 1 ] /= totArea;
  currentCluster->centroid[ 2 ] /= totArea;



  currentCluster->mass = totArea;

  // compute the radius of the cluster
  SCVT radius = 0.0;
  currentCluster->radius = 0.0;
  for ( LO i = 0; i < currentCluster->nelems; i++ ) {
    radius = std::sqrt(
      ( currentCluster->centroid[0]
      - ( *centroids )[ 3 * ( *globalElems )[ i ] ] ) *
      ( currentCluster->centroid[0]
      - ( *centroids )[ 3 * ( *globalElems )[ i ] ] ) +
      ( currentCluster->centroid[1]
      - ( *centroids )[ 3 * ( *globalElems )[ i ] + 1 ] ) *
      ( currentCluster->centroid[1]
      - ( *centroids )[ 3 * ( *globalElems )[ i ] + 1 ] ) +
      ( currentCluster->centroid[2]
      - ( *centroids )[ 3 * ( *globalElems )[ i ] + 2 ] ) *
      ( currentCluster->centroid[2]
      - ( *centroids )[ 3 * ( *globalElems )[ i ] + 2 ] ) );
    if ( radius > currentCluster->radius ) {
      currentCluster->radius = radius;
    }
  }

  // stop splitting if we are deep enough
  if ( currentCluster->nelems < maxElemsPerCluster ) {
    return;
  }

  // create the covariance matrix of the cluster
  SCVT *variance = new SCVT[9];
  for ( int i = 0; i < 9; i++ ) {
    variance[i] = (SCVT) 0.0;
  }

  for ( LO i = 0; i < currentCluster->nelems; i++ ) {
    for ( int j = 0; j < 3; j++ ) {
      for ( int k = 0; k < 3; k++ ) {
        variance[ 3 * k + j ] += this->area[ ( *globalElems )[ i ] ] *
          ( ( *centroids )[ 3 * ( *globalElems )[ i ] + j ]
          - currentCluster->centroid[j] ) *
          ( ( *centroids )[ 3 * ( *globalElems )[ i ] + k ]
          - currentCluster->centroid[k] );
      }
    }
  }

  FullMatrix<LO, SCVT> varMatrix( 3, 3, variance, false );

  SCVT eigenvalues[3];
  SCVT eigenvectors[9];
  varMatrix.eigs( eigenvectors, eigenvalues );

  // find the largest eigenvalue
  SCVT maxEig = 0;
  int idx = 0;
  for ( int i = 0; i < 3; i++ ) {
    if ( fabs( maxEig ) < fabs( eigenvalues[i] ) ) {
      maxEig = eigenvalues[i];
      idx = i;
    }
  }

  SCVT normal[3]; //normal to the cutting plane
  normal[0] = eigenvectors[3 * idx];
  normal[1] = eigenvectors[3 * idx + 1];
  normal[2] = eigenvectors[3 * idx + 2];
  delete [] variance;

  // split elements into two groups by the plane defined by the normal vector
  std::vector<LO> *leftElems = new std::vector<LO>( );
  std::vector<LO> *rightElems = new std::vector<LO>( );

  // auxiliary sets to count number of nodes in left and right clusters
  std::set<LO>* leftNodes = new std::set<LO>;
  std::set<LO>* rightNodes = new std::set<LO>;
  LO currentNodes[3];

  SCVT s;
  LO nLeftElems = 0, nRightElems = 0;
  for ( LO i = 0; i < currentCluster->nelems; i++ ) {
    this->getElement( ( *globalElems )[ i ], currentNodes );
    s = ( ( *centroids )[ 3 * ( *globalElems )[ i ] ]
      - currentCluster->centroid[0] ) * normal[0] +
      ( ( *centroids )[ 3 * ( *globalElems )[ i ] + 1 ]
      - currentCluster->centroid[1] ) * normal[1] +
      ( ( *centroids )[ 3 * ( *globalElems )[ i ] + 2 ]
      - currentCluster->centroid[2] ) * normal[2];
    if ( s >= 0 ) {
      leftElems->push_back( ( *globalElems )[ i ] );
      nLeftElems++;
      leftNodes->insert( currentNodes[0] );
      leftNodes->insert( currentNodes[1] );
      leftNodes->insert( currentNodes[2] );
    } else {
      rightElems->push_back( ( *globalElems )[ i ] );
      nRightElems++;
      rightNodes->insert( currentNodes[0] );
      rightNodes->insert( currentNodes[1] );
      rightNodes->insert( currentNodes[2] );
    }
  }

  // create new clusters and split them recursively
  if ( nLeftElems > 0 ) {
    BECluster<LO, SC> *leftCluster = new BECluster<LO, SC>( );
    leftCluster->elems = leftElems;
    leftCluster->nelems = nLeftElems;
    leftCluster->nnodes = leftNodes->size( );
    leftCluster->nodes = leftNodes;
    clusterTree.addSon( leftCluster );
  }
  if ( nRightElems > 0 ) {
    BECluster<LO, SC> *rightCluster = new BECluster<LO, SC>( );
    rightCluster->elems = rightElems;
    rightCluster->nelems = nRightElems;
    rightCluster->nnodes = rightNodes->size( );
    rightCluster->nodes = rightNodes;
    clusterTree.addSon( rightCluster );
  }

  for ( int i = 0; i < clusterTree.getNSons( ); i++ ) {
    clusterTree.setPointerToSon( i );
    doNestedDissection( clusterTree, maxElemsPerCluster );
    clusterTree.setPointerToParent( );
  }

}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::getNestedProjectionDiam(
  Tree<BECluster<LO, SC>*, LO > & clusterTree,
  std::vector< std::vector< SC > > & data
  ) const {
  clusterTree.setPointerToRoot( );
  this->doGetNestedProjectionDiam( clusterTree, data );
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::doGetNestedProjectionDiam(
  Tree<BECluster<LO, SC>*, LO > & clusterTree,
  std::vector< std::vector< SC > > & data
  ) const {

  BECluster<LO, SC>* currentCluster = clusterTree.get( );
  LO nSons = clusterTree.getNSons( );
  //! depends on the element dimension: 1D-mass, 2D-sqrt(mass)
  //data[clusterTree.getLevel( )].push_back( std::sqrt( currentCluster->mass ) );
  data[clusterTree.getLevel( )].push_back( 2.0 * currentCluster->radius );
  if ( nSons != 0 ) {
    for ( LO i = 0; i < nSons; ++i ) {
      clusterTree.setPointerToSon( i );
      this->doGetNestedProjectionDiam( clusterTree, data );
      clusterTree.setPointerToParent( );
    }
  }
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::getNestedProjectionMassRatio(
  Tree<BECluster<LO, SC>*, LO > & clusterTree,
  std::vector< std::vector< SC > > & data
  ) const {
  clusterTree.setPointerToRoot( );
  this->doNestedProjectionMassRatio( clusterTree, data );
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::doNestedProjectionMassRatio(
  Tree<BECluster<LO, SC>*, LO > & clusterTree,
  std::vector< std::vector< SC> > & data
  ) const {

  BECluster<LO, SC>* currentCluster = clusterTree.get( );
  BECluster<LO, SC>* childCluster;
  std::vector<LO>* clusterElems;
  LO elemID;
  SC elemArea;
  SCVT massPar, massSon;
  LO nSons = clusterTree.getNSons( );
  massPar = currentCluster->mass;
  if ( nSons != 0 ) {
    clusterTree.setPointerToSon( 0 );
    this->doNestedProjectionMassRatio( clusterTree, data );
    clusterTree.setPointerToParent( );

    clusterTree.setPointerToSon( 1 );
    this->doNestedProjectionMassRatio( clusterTree, data );
    clusterTree.setPointerToParent( );

    clusterElems = currentCluster->elems;
    for ( LO i = 0; i < clusterElems->size( ); i++ ) {
      elemID = clusterElems->at( i );
      elemArea = this->getElemArea( elemID );
      // maps phi_elemID to next-level cluster
      data[elemID].push_back( elemArea / massPar );
    }
  } else {
    clusterElems = currentCluster->elems;
    for ( LO i = 0; i < clusterElems->size( ); i++ ) {
      elemID = clusterElems->at( i );
      elemArea = this->getElemArea( elemID );
      // maps phi_elemID to higher-level cluster
      data[elemID].push_back( elemArea / massPar );
    }
  }
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::getNestedProjectionMassRatio(
  Tree<BECluster<LO, SC>*, LO > & clusterTree,
  std::vector< Vector< LO, SC > * > & data
  ) const {
  clusterTree.setPointerToRoot( );
  this->doNestedProjectionMassRatio( clusterTree, data );
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::doNestedProjectionMassRatio(
  Tree<BECluster<LO, SC>*, LO > & clusterTree,
  std::vector< Vector< LO, SC > * > & data
  ) const {

  BECluster<LO, SC>* currentCluster = clusterTree.get( );
  std::vector<LO>* clusterElems;
  LO elemID;
  LO currentLevel;
  SC elemArea;
  SCVT massPar, massSon;
  LO nSons = clusterTree.getNSons( );
  massPar = currentCluster->mass;
  if ( nSons != 0 ) {
    clusterTree.setPointerToSon( 0 );
    this->doNestedProjectionMassRatio( clusterTree, data );
    clusterTree.setPointerToParent( );

    clusterTree.setPointerToSon( 1 );
    this->doNestedProjectionMassRatio( clusterTree, data );
    clusterTree.setPointerToParent( );

    currentLevel = clusterTree.getLevel( );
    clusterElems = currentCluster->elems;
    for ( LO i = 0; i < clusterElems->size( ); i++ ) {
      elemID = clusterElems->at( i );
      elemArea = this->getElemArea( elemID );
      // maps phi_elemID to next-level cluster
      //data[elemID].push_back( elemArea/massPar );
      data[currentLevel]->set( elemID, elemArea / massPar );
    }
  } else {
    clusterElems = currentCluster->elems;
    currentLevel = clusterTree.getLevel( );
    for ( LO i = 0; i < clusterElems->size( ); i++ ) {
      elemID = clusterElems->at( i );
      elemArea = this->getElemArea( elemID );
      // maps phi_elemID to higher-level cluster
      //data[elemID].push_back( elemArea/massPar ); 
      data[currentLevel]->set( elemID, elemArea / massPar );
    }
  }
}

//! Specialization for LO = int, SC = complex<double>

template<>
void SurfaceMesh3D< int, std::complex<double> >::printParaviewVtu(
  const string& meshFile,
  int nNodal,
  string * nodeNames,
  Vector< int, std::complex<double> > ** nodalData,
  int nElem,
  string * elemNames,
  Vector< int, std::complex<double> > ** elemData
  ) {

#ifdef VERBOSE
  std::cout << "Printing '" << meshFile << "' ... ";
#endif

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 6 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"
    << std::endl;
  file_vtu << "  <UnstructuredGrid>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << this->nNodes <<
    "\" NumberOfCells=\"" << this->nElems << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\""
    " NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  for ( int i = 0; i < this->nNodes; i++ ) {
    file_vtu << "          "
      << this->nodes[ 3 * i ] << " "
      << this->nodes[ 3 * i + 1 ] << " "
      << this->nodes[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Cells>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\""
    " format=\"ascii\">" << std::endl;

  for ( int i = 0; i < this->nElems; i++ ) {
    file_vtu << "          "
      << this->elems[ 3 * i ] << " "
      << this->elems[ 3 * i + 1 ] << " "
      << this->elems[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"offsets\""
    " format=\"ascii\">" << std::endl;

  for ( int offset = 1; offset <= this->nElems; offset++ ) {
    file_vtu << "          " << offset * 3 << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"UInt8\" Name=\"types\""
    " format=\"ascii\">" << std::endl;
  for ( int i = 1; i <= this->nElems; i++ ) {
    file_vtu << "          5" << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Cells>" << std::endl;

  string header;
  if ( nNodal > 0 ) {
    header = nodeNames[ 0 ] + "_real," + nodeNames[ 0 ] + "_imag";
    for ( int j = 1; j < nNodal; j++ ) {
      header += "," + nodeNames[ j ] + "_real";
      header += "," + nodeNames[ j ] + "_imag";
    }
    file_vtu << "      <PointData Scalars=\"" + header + "\">" << std::endl;
    for ( int j = 0; j < nNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + nodeNames[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nNodes; i++ ) {
        file_vtu << "          " << nodalData[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + nodeNames[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nNodes; i++ ) {
        file_vtu << "          " << nodalData[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }
    file_vtu << "      </PointData>" << std::endl;
  }

  if ( nElem > 0 ) {
    header = elemNames[ 0 ] + "_real," + elemNames[ 0 ] + "_imag";
    for ( int j = 1; j < nElem; j++ ) {
      header += "," + elemNames[ j ] + "_real";
      header += "," + elemNames[ j ] + "_imag";
    }
    file_vtu << "      <CellData Scalars=\"" + header + "\">" << std::endl;
    for ( int j = 0; j < nElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + elemNames[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nElems; i++ ) {
        file_vtu << "          " << elemData[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + elemNames[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nElems; i++ ) {
        file_vtu << "          " << elemData[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }
    file_vtu << "      </CellData>" << std::endl;
  }

  file_vtu << "    </Piece>" << std::endl;
  file_vtu << "  </UnstructuredGrid>" << std::endl;
  file_vtu << "</VTKFile>" << std::endl;
  file_vtu.close( );

#ifdef VERBOSE
  std::cout << "done." << std::endl;
#endif

}

template<>
void SurfaceMesh3D< long, std::complex<double> >::printParaviewVtu(
  const string& meshFile,
  int nNodal,
  string* nodeNames,
  Vector< long, std::complex<double> >** nodalData,
  int nElem,
  string* elemNames,
  Vector< long, std::complex<double> >** elemData
  ) {

#ifdef VERBOSE
  std::cout << "Printing '" << meshFile << "' ... ";
#endif

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 6 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"
    << std::endl;
  file_vtu << "  <UnstructuredGrid>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << this->nNodes <<
    "\" NumberOfCells=\"" << this->nElems << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\""
    " NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  for ( long i = 0; i < this->nNodes; i++ ) {
    file_vtu << "          "
      << this->nodes[ 3 * i ] << " "
      << this->nodes[ 3 * i + 1 ] << " "
      << this->nodes[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Cells>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\""
    " format=\"ascii\">" << std::endl;

  for ( long i = 0; i < this->nElems; i++ ) {
    file_vtu << "          "
      << this->elems[ 3 * i ] << " "
      << this->elems[ 3 * i + 1 ] << " "
      << this->elems[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"offsets\""
    " format=\"ascii\">" << std::endl;

  for ( long offset = 1; offset <= this->nElems; offset++ ) {
    file_vtu << "          " << offset * 3 << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"UInt8\" Name=\"types\""
    " format=\"ascii\">" << std::endl;
  for ( long i = 1; i <= this->nElems; i++ ) {
    file_vtu << "          5" << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Cells>" << std::endl;

  string header;
  if ( nNodal > 0 ) {
    header = nodeNames[ 0 ] + "_real," + nodeNames[ 0 ] + "_imag";
    for ( long j = 1; j < nNodal; j++ ) {
      header += "," + nodeNames[ j ] + "_real";
      header += "," + nodeNames[ j ] + "_imag";
    }
    file_vtu << "      <PointData Scalars=\"" + header + "\">" << std::endl;
    for ( int j = 0; j < nNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + nodeNames[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nNodes; i++ ) {
        file_vtu << "          " << nodalData[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + nodeNames[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nNodes; i++ ) {
        file_vtu << "          " << nodalData[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }
    file_vtu << "      </PointData>" << std::endl;
  }

  if ( nElem > 0 ) {
    header = elemNames[ 0 ] + "_real," + elemNames[ 0 ] + "_imag";
    for ( int j = 1; j < nElem; j++ ) {
      header += "," + elemNames[ j ] + "_real";
      header += "," + elemNames[ j ] + "_imag";
    }
    file_vtu << "      <CellData Scalars=\"" + header + "\">" << std::endl;
    for ( int j = 0; j < nElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + elemNames[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nElems; i++ ) {
        file_vtu << "          " << elemData[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + elemNames[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nElems; i++ ) {
        file_vtu << "          " << elemData[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }
    file_vtu << "      </CellData>" << std::endl;
  }

  file_vtu << "    </Piece>" << std::endl;
  file_vtu << "  </UnstructuredGrid>" << std::endl;
  file_vtu << "</VTKFile>" << std::endl;
  file_vtu.close( );

#ifdef VERBOSE
  std::cout << "done." << std::endl;
#endif

}

template<>
void SurfaceMesh3D< int, std::complex<float> >::printParaviewVtu(
  const string& meshFile,
  int nNodal,
  string * nodeNames,
  Vector< int, std::complex<float> >* * nodalData,
  int nElem,
  string * elemNames,
  Vector< int, std::complex<float> > ** elemData
  ) {

#ifdef VERBOSE
  std::cout << "Printing '" << meshFile << "' ... ";
#endif

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 6 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"
    << std::endl;
  file_vtu << "  <UnstructuredGrid>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << this->nNodes <<
    "\" NumberOfCells=\"" << this->nElems << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\""
    " NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  for ( int i = 0; i < this->nNodes; i++ ) {
    file_vtu << "          "
      << this->nodes[ 3 * i ] << " "
      << this->nodes[ 3 * i + 1 ] << " "
      << this->nodes[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Cells>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\""
    " format=\"ascii\">" << std::endl;

  for ( int i = 0; i < this->nElems; i++ ) {
    file_vtu << "          "
      << this->elems[ 3 * i ] << " "
      << this->elems[ 3 * i + 1 ] << " "
      << this->elems[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"offsets\""
    " format=\"ascii\">" << std::endl;

  for ( int offset = 1; offset <= this->nElems; offset++ ) {
    file_vtu << "          " << offset * 3 << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"UInt8\" Name=\"types\""
    " format=\"ascii\">" << std::endl;
  for ( int i = 1; i <= this->nElems; i++ ) {
    file_vtu << "          5" << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Cells>" << std::endl;

  string header;
  if ( nNodal > 0 ) {
    header = nodeNames[ 0 ] + "_real" + nodeNames[ 0 ] + "_imag";
    for ( int j = 1; j < nNodal; j++ ) {
      header += "," + nodeNames[ j ] + "_real";
      header += "," + nodeNames[ j ] + "_imag";
    }
    file_vtu << "      <PointData Scalars=\"" + header + "\">" << std::endl;
    for ( int j = 0; j < nNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + nodeNames[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nNodes; i++ ) {
        file_vtu << "          " << nodalData[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + nodeNames[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nNodes; i++ ) {
        file_vtu << "          " << nodalData[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }
    file_vtu << "      </PointData>" << std::endl;
  }

  if ( nElem > 0 ) {
    header = elemNames[ 0 ] + "_real," + elemNames[ 0 ] + "_imag";
    for ( int j = 1; j < nElem; j++ ) {
      header += "," + elemNames[ j ] + "_real";
      header += "," + elemNames[ j ] + "_imag";
    }
    file_vtu << "      <CellData Scalars=\"" + header + "\">" << std::endl;
    for ( int j = 0; j < nElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + elemNames[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nElems; i++ ) {
        file_vtu << "          " << elemData[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + elemNames[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nElems; i++ ) {
        file_vtu << "          " << elemData[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }
    file_vtu << "      </CellData>" << std::endl;
  }

  file_vtu << "    </Piece>" << std::endl;
  file_vtu << "  </UnstructuredGrid>" << std::endl;
  file_vtu << "</VTKFile>" << std::endl;
  file_vtu.close( );

#ifdef VERBOSE
  std::cout << "done." << std::endl;
#endif

}

template<>
void SurfaceMesh3D< long, std::complex<float> >::printParaviewVtu(
  const string& meshFile,
  int nNodal,
  string* nodeNames,
  Vector< long, std::complex<float> >** nodalData,
  int nElem,
  string* elemNames,
  Vector< long, std::complex<float> >** elemData
  ) {

#ifdef VERBOSE
  std::cout << "Printing '" << meshFile << "' ... ";
#endif

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 6 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"
    << std::endl;
  file_vtu << "  <UnstructuredGrid>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << this->nNodes <<
    "\" NumberOfCells=\"" << this->nElems << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\""
    " NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  for ( long i = 0; i < this->nNodes; i++ ) {
    file_vtu << "          "
      << this->nodes[ 3 * i ] << " "
      << this->nodes[ 3 * i + 1 ] << " "
      << this->nodes[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Cells>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\""
    " format=\"ascii\">" << std::endl;

  for ( long i = 0; i < this->nElems; i++ ) {
    file_vtu << "          "
      << this->elems[ 3 * i ] << " "
      << this->elems[ 3 * i + 1 ] << " "
      << this->elems[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"offsets\""
    " format=\"ascii\">" << std::endl;

  for ( long offset = 1; offset <= this->nElems; offset++ ) {
    file_vtu << "          " << offset * 3 << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"UInt8\" Name=\"types\""
    " format=\"ascii\">" << std::endl;
  for ( long i = 1; i <= this->nElems; i++ ) {
    file_vtu << "          5" << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Cells>" << std::endl;

  string header;
  if ( nNodal > 0 ) {
    header = nodeNames[ 0 ] + "_real," + nodeNames[ 0 ] + "_imag";
    for ( long j = 1; j < nNodal; j++ ) {
      header += "," + nodeNames[ j ] + "_real";
      header += "," + nodeNames[ j ] + "_imag";
    }
    file_vtu << "      <PointData Scalars=\"" + header + "\">" << std::endl;
    for ( int j = 0; j < nNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + nodeNames[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nNodes; i++ ) {
        file_vtu << "          " << nodalData[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + nodeNames[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nNodes; i++ ) {
        file_vtu << "          " << nodalData[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }
    file_vtu << "      </PointData>" << std::endl;
  }

  if ( nElem > 0 ) {
    header = elemNames[ 0 ] + "_real," + elemNames[ 0 ] + "_imag";
    for ( int j = 1; j < nElem; j++ ) {
      header += "," + elemNames[ j ] + "_real";
      header += "," + elemNames[ j ] + "_imag";
    }
    file_vtu << "      <CellData Scalars=\"" + header + "\">" << std::endl;
    for ( int j = 0; j < nElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + elemNames[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nElems; i++ ) {
        file_vtu << "          " << elemData[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + elemNames[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nElems; i++ ) {
        file_vtu << "          " << elemData[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }
    file_vtu << "      </CellData>" << std::endl;
  }

  file_vtu << "    </Piece>" << std::endl;
  file_vtu << "  </UnstructuredGrid>" << std::endl;
  file_vtu << "</VTKFile>" << std::endl;
  file_vtu.close( );

#ifdef VERBOSE
  std::cout << "done." << std::endl;
#endif

}

template<>
void SurfaceMesh3D< int, std::complex< double > >::printParaviewVtu(
  const string & meshFile,
  const std::vector< string > * nodeNames,
  const std::vector< Vector< int, std::complex< double > > * > * nodalData,
  const std::vector< string > * elemNames,
  const std::vector< Vector< int, std::complex< double > > * > * elemData,
  const std::vector< string > * nodeVNames,
  const std::vector< Vector< int, std::complex< double > > * > * nodalVData,
  const std::vector< string > * elemVNames,
  const std::vector< Vector< int, std::complex< double > > * > * elemVData
  ) const {

#ifdef VERBOSE
  std::cout << "Printing '" << meshFile << "' ... ";
#endif

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 6 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"
    << std::endl;
  file_vtu << "  <UnstructuredGrid>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << this->nNodes <<
    "\" NumberOfCells=\"" << this->nElems << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\""
    " NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  for ( int i = 0; i < this->nNodes; i++ ) {
    file_vtu << "          "
      << this->nodes[ 3 * i ] << " "
      << this->nodes[ 3 * i + 1 ] << " "
      << this->nodes[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Cells>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\""
    " format=\"ascii\">" << std::endl;

  for ( int i = 0; i < this->nElems; i++ ) {
    file_vtu << "          "
      << this->elems[ 3 * i ] << " "
      << this->elems[ 3 * i + 1 ] << " "
      << this->elems[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"offsets\""
    " format=\"ascii\">" << std::endl;

  for ( int offset = 1; offset <= this->nElems; offset++ ) {
    file_vtu << "          " << offset * 3 << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"UInt8\" Name=\"types\""
    " format=\"ascii\">" << std::endl;
  for ( int i = 1; i <= this->nElems; i++ ) {
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
      header = ( *nodeNames )[ 0 ] + "_real," + ( *nodeNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nNodal; ++j ) {
        header += "," + ( *nodeNames )[ j ] + "_real";
        header += "," + ( *nodeNames )[ j ] + "_imag";
      }
      file_vtu << "Scalars=\"" + header + "\"";
    }

    if ( nVNodal > 0 ) {
      vheader =
        ( *nodeVNames )[ 0 ] + "_real," + ( *nodeVNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nVNodal; ++j ) {
        vheader += "," + ( *nodeVNames )[ j ] + "_real";
        vheader += "," + ( *nodeVNames )[ j ] + "_imag";
      }
      file_vtu << " Vectors=\"" + vheader + "\"";
    }

    file_vtu << ">" << std::endl;

    for ( int j = 0; j < nNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *nodeNames )[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *nodeNames )[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    for ( int j = 0; j < nVNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *nodeVNames )[ j ] + "_real" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < 3 * this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalVData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *nodeVNames )[ j ] + "_imag" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < 3 * this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalVData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    file_vtu << "      </PointData>" << std::endl;
  }

  int nElem = 0;
  if ( elemData ) nElem = elemData->size( );
  int nVElem = 0;
  if ( elemVData ) nVElem = elemVData->size( );

  if ( nElem > 0 || nVElem > 0 ) {
    file_vtu << "      <CellData ";

    if ( nElem > 0 ) {
      header = ( *elemNames )[ 0 ] + "_real," + ( *elemNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nElem; ++j ) {
        header += "," + ( *elemNames )[ j ] + "_real";
        header += "," + ( *elemNames )[ j ] + "_imag";
      }
      file_vtu << "Scalars=\"" + header + "\"";
    }

    if ( nVElem > 0 ) {
      vheader =
        ( *elemVNames )[ 0 ] + "_real," + ( *elemVNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nVElem; ++j ) {
        vheader += "," + ( *elemVNames )[ j ] + "_real";
        vheader += "," + ( *elemVNames )[ j ] + "_imag";
      }
      file_vtu << " Vectors=\"" + vheader + "\"";
    }

    file_vtu << ">" << std::endl;

    for ( int j = 0; j < nElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *elemNames )[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nElems; ++i ) {
        file_vtu << "          " << ( *elemData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *elemNames )[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nElems; ++i ) {
        file_vtu << "          " << ( *elemData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    for ( int j = 0; j < nVElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *elemVNames )[ j ] + "_real" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < 3 * this->nElems; ++i ) {
        file_vtu << "          " << ( *elemVData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *elemVNames )[ j ] + "_imag" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < 3 * this->nElems; ++i ) {
        file_vtu << "          " << ( *elemVData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    file_vtu << "      </CellData>" << std::endl;
  }

  file_vtu << "    </Piece>" << std::endl;
  file_vtu << "  </UnstructuredGrid>" << std::endl;
  file_vtu << "</VTKFile>" << std::endl;
  file_vtu.close( );

#ifdef VERBOSE
  std::cout << "done." << std::endl;
#endif

}

template<>
void SurfaceMesh3D< int, std::complex< float > >::printParaviewVtu(
  const string & meshFile,
  const std::vector< string > * nodeNames,
  const std::vector< Vector< int, std::complex< float > > * > * nodalData,
  const std::vector< string > * elemNames,
  const std::vector< Vector< int, std::complex< float > > * > * elemData,
  const std::vector< string > * nodeVNames,
  const std::vector< Vector< int, std::complex< float > > * > * nodalVData,
  const std::vector< string > * elemVNames,
  const std::vector< Vector< int, std::complex< float > > * > * elemVData
  ) const {

#ifdef VERBOSE
  std::cout << "Printing '" << meshFile << "' ... ";
#endif

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 6 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"
    << std::endl;
  file_vtu << "  <UnstructuredGrid>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << this->nNodes <<
    "\" NumberOfCells=\"" << this->nElems << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\""
    " NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  for ( int i = 0; i < this->nNodes; i++ ) {
    file_vtu << "          "
      << this->nodes[ 3 * i ] << " "
      << this->nodes[ 3 * i + 1 ] << " "
      << this->nodes[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Cells>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\""
    " format=\"ascii\">" << std::endl;

  for ( int i = 0; i < this->nElems; i++ ) {
    file_vtu << "          "
      << this->elems[ 3 * i ] << " "
      << this->elems[ 3 * i + 1 ] << " "
      << this->elems[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"offsets\""
    " format=\"ascii\">" << std::endl;

  for ( int offset = 1; offset <= this->nElems; offset++ ) {
    file_vtu << "          " << offset * 3 << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"UInt8\" Name=\"types\""
    " format=\"ascii\">" << std::endl;
  for ( int i = 1; i <= this->nElems; i++ ) {
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
      header = ( *nodeNames )[ 0 ] + "_real," + ( *nodeNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nNodal; ++j ) {
        header += "," + ( *nodeNames )[ j ] + "_real";
        header += "," + ( *nodeNames )[ j ] + "_imag";
      }
      file_vtu << "Scalars=\"" + header + "\"";
    }

    if ( nVNodal > 0 ) {
      vheader =
        ( *nodeVNames )[ 0 ] + "_real," + ( *nodeVNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nVNodal; ++j ) {
        vheader += "," + ( *nodeVNames )[ j ] + "_real";
        vheader += "," + ( *nodeVNames )[ j ] + "_imag";
      }
      file_vtu << " Vectors=\"" + vheader + "\"";
    }

    file_vtu << ">" << std::endl;

    for ( int j = 0; j < nNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *nodeNames )[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *nodeNames )[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    for ( int j = 0; j < nVNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *nodeVNames )[ j ] + "_real" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < 3 * this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalVData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *nodeVNames )[ j ] + "_imag" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < 3 * this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalVData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    file_vtu << "      </PointData>" << std::endl;
  }

  int nElem = 0;
  if ( elemData ) nElem = elemData->size( );
  int nVElem = 0;
  if ( elemVData ) nVElem = elemVData->size( );

  if ( nElem > 0 || nVElem > 0 ) {
    file_vtu << "      <CellData ";

    if ( nElem > 0 ) {
      header = ( *elemNames )[ 0 ] + "_real," + ( *elemNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nElem; ++j ) {
        header += "," + ( *elemNames )[ j ] + "_real";
        header += "," + ( *elemNames )[ j ] + "_imag";
      }
      file_vtu << "Scalars=\"" + header + "\"";
    }

    if ( nVElem > 0 ) {
      vheader =
        ( *elemVNames )[ 0 ] + "_real," + ( *elemVNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nVElem; ++j ) {
        vheader += "," + ( *elemVNames )[ j ] + "_real";
        vheader += "," + ( *elemVNames )[ j ] + "_imag";
      }
      file_vtu << " Vectors=\"" + vheader + "\"";
    }

    file_vtu << ">" << std::endl;

    for ( int j = 0; j < nElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *elemNames )[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nElems; ++i ) {
        file_vtu << "          " << ( *elemData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *elemNames )[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < this->nElems; ++i ) {
        file_vtu << "          " << ( *elemData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    for ( int j = 0; j < nVElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *elemVNames )[ j ] + "_real" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < 3 * this->nElems; ++i ) {
        file_vtu << "          " << ( *elemVData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *elemVNames )[ j ] + "_imag" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( int i = 0; i < 3 * this->nElems; ++i ) {
        file_vtu << "          " << ( *elemVData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    file_vtu << "      </CellData>" << std::endl;
  }

  file_vtu << "    </Piece>" << std::endl;
  file_vtu << "  </UnstructuredGrid>" << std::endl;
  file_vtu << "</VTKFile>" << std::endl;
  file_vtu.close( );

#ifdef VERBOSE
  std::cout << "done." << std::endl;
#endif

}

template<>
void SurfaceMesh3D< long, std::complex< double > >::printParaviewVtu(
  const string & meshFile,
  const std::vector< string > * nodeNames,
  const std::vector< Vector< long, std::complex< double > > * > * nodalData,
  const std::vector< string > * elemNames,
  const std::vector< Vector< long, std::complex< double > > * > * elemData,
  const std::vector< string > * nodeVNames,
  const std::vector< Vector< long, std::complex< double > > * > * nodalVData,
  const std::vector< string > * elemVNames,
  const std::vector< Vector< long, std::complex< double > > * > * elemVData
  ) const {

#ifdef VERBOSE
  std::cout << "Printing '" << meshFile << "' ... ";
#endif

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 6 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"
    << std::endl;
  file_vtu << "  <UnstructuredGrid>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << this->nNodes <<
    "\" NumberOfCells=\"" << this->nElems << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\""
    " NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  for ( long i = 0; i < this->nNodes; i++ ) {
    file_vtu << "          "
      << this->nodes[ 3 * i ] << " "
      << this->nodes[ 3 * i + 1 ] << " "
      << this->nodes[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Cells>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\""
    " format=\"ascii\">" << std::endl;

  for ( long i = 0; i < this->nElems; i++ ) {
    file_vtu << "          "
      << this->elems[ 3 * i ] << " "
      << this->elems[ 3 * i + 1 ] << " "
      << this->elems[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"offsets\""
    " format=\"ascii\">" << std::endl;

  for ( long offset = 1; offset <= this->nElems; offset++ ) {
    file_vtu << "          " << offset * 3 << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"UInt8\" Name=\"types\""
    " format=\"ascii\">" << std::endl;
  for ( long i = 1; i <= this->nElems; i++ ) {
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
      header = ( *nodeNames )[ 0 ] + "_real," + ( *nodeNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nNodal; ++j ) {
        header += "," + ( *nodeNames )[ j ] + "_real";
        header += "," + ( *nodeNames )[ j ] + "_imag";
      }
      file_vtu << "Scalars=\"" + header + "\"";
    }

    if ( nVNodal > 0 ) {
      vheader =
        ( *nodeVNames )[ 0 ] + "_real," + ( *nodeVNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nVNodal; ++j ) {
        vheader += "," + ( *nodeVNames )[ j ] + "_real";
        vheader += "," + ( *nodeVNames )[ j ] + "_imag";
      }
      file_vtu << " Vectors=\"" + vheader + "\"";
    }

    file_vtu << ">" << std::endl;

    for ( int j = 0; j < nNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *nodeNames )[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *nodeNames )[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    for ( int j = 0; j < nVNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *nodeVNames )[ j ] + "_real" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < 3 * this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalVData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *nodeVNames )[ j ] + "_imag" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < 3 * this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalVData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    file_vtu << "      </PointData>" << std::endl;
  }

  int nElem = 0;
  if ( elemData ) nElem = elemData->size( );
  int nVElem = 0;
  if ( elemVData ) nVElem = elemVData->size( );

  if ( nElem > 0 || nVElem > 0 ) {
    file_vtu << "      <CellData ";

    if ( nElem > 0 ) {
      header = ( *elemNames )[ 0 ] + "_real," + ( *elemNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nElem; ++j ) {
        header += "," + ( *elemNames )[ j ] + "_real";
        header += "," + ( *elemNames )[ j ] + "_imag";
      }
      file_vtu << "Scalars=\"" + header + "\"";
    }

    if ( nVElem > 0 ) {
      vheader =
        ( *elemVNames )[ 0 ] + "_real," + ( *elemVNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nVElem; ++j ) {
        vheader += "," + ( *elemVNames )[ j ] + "_real";
        vheader += "," + ( *elemVNames )[ j ] + "_imag";
      }
      file_vtu << " Vectors=\"" + vheader + "\"";
    }

    file_vtu << ">" << std::endl;

    for ( int j = 0; j < nElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *elemNames )[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nElems; ++i ) {
        file_vtu << "          " << ( *elemData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *elemNames )[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nElems; ++i ) {
        file_vtu << "          " << ( *elemData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    for ( int j = 0; j < nVElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *elemVNames )[ j ] + "_real" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < 3 * this->nElems; ++i ) {
        file_vtu << "          " << ( *elemVData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *elemVNames )[ j ] + "_imag" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < 3 * this->nElems; ++i ) {
        file_vtu << "          " << ( *elemVData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    file_vtu << "      </CellData>" << std::endl;
  }

  file_vtu << "    </Piece>" << std::endl;
  file_vtu << "  </UnstructuredGrid>" << std::endl;
  file_vtu << "</VTKFile>" << std::endl;
  file_vtu.close( );

#ifdef VERBOSE
  std::cout << "done." << std::endl;
#endif

}

template<>
void SurfaceMesh3D< long, std::complex< float > >::printParaviewVtu(
  const string & meshFile,
  const std::vector< string > * nodeNames,
  const std::vector< Vector< long, std::complex< float > > * > * nodalData,
  const std::vector< string > * elemNames,
  const std::vector< Vector< long, std::complex< float > > * > * elemData,
  const std::vector< string > * nodeVNames,
  const std::vector< Vector< long, std::complex< float > > * > * nodalVData,
  const std::vector< string > * elemVNames,
  const std::vector< Vector< long, std::complex< float > > * > * elemVData
  ) const {

#ifdef VERBOSE
  std::cout << "Printing '" << meshFile << "' ... ";
#endif

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 6 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">"
    << std::endl;
  file_vtu << "  <UnstructuredGrid>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << this->nNodes <<
    "\" NumberOfCells=\"" << this->nElems << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\""
    " NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  for ( long i = 0; i < this->nNodes; i++ ) {
    file_vtu << "          "
      << this->nodes[ 3 * i ] << " "
      << this->nodes[ 3 * i + 1 ] << " "
      << this->nodes[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Cells>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\""
    " format=\"ascii\">" << std::endl;

  for ( long i = 0; i < this->nElems; i++ ) {
    file_vtu << "          "
      << this->elems[ 3 * i ] << " "
      << this->elems[ 3 * i + 1 ] << " "
      << this->elems[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"offsets\""
    " format=\"ascii\">" << std::endl;

  for ( long offset = 1; offset <= this->nElems; offset++ ) {
    file_vtu << "          " << offset * 3 << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"UInt8\" Name=\"types\""
    " format=\"ascii\">" << std::endl;
  for ( long i = 1; i <= this->nElems; i++ ) {
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
      header = ( *nodeNames )[ 0 ] + "_real," + ( *nodeNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nNodal; ++j ) {
        header += "," + ( *nodeNames )[ j ] + "_real";
        header += "," + ( *nodeNames )[ j ] + "_imag";
      }
      file_vtu << "Scalars=\"" + header + "\"";
    }

    if ( nVNodal > 0 ) {
      vheader =
        ( *nodeVNames )[ 0 ] + "_real," + ( *nodeVNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nVNodal; ++j ) {
        vheader += "," + ( *nodeVNames )[ j ] + "_real";
        vheader += "," + ( *nodeVNames )[ j ] + "_imag";
      }
      file_vtu << " Vectors=\"" + vheader + "\"";
    }

    file_vtu << ">" << std::endl;

    for ( int j = 0; j < nNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *nodeNames )[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *nodeNames )[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    for ( int j = 0; j < nVNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *nodeVNames )[ j ] + "_real" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < 3 * this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalVData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *nodeVNames )[ j ] + "_imag" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < 3 * this->nNodes; ++i ) {
        file_vtu << "          " << ( *nodalVData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    file_vtu << "      </PointData>" << std::endl;
  }

  int nElem = 0;
  if ( elemData ) nElem = elemData->size( );
  int nVElem = 0;
  if ( elemVData ) nVElem = elemVData->size( );

  if ( nElem > 0 || nVElem > 0 ) {
    file_vtu << "      <CellData ";

    if ( nElem > 0 ) {
      header = ( *elemNames )[ 0 ] + "_real," + ( *elemNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nElem; ++j ) {
        header += "," + ( *elemNames )[ j ] + "_real";
        header += "," + ( *elemNames )[ j ] + "_imag";
      }
      file_vtu << "Scalars=\"" + header + "\"";
    }

    if ( nVElem > 0 ) {
      vheader =
        ( *elemVNames )[ 0 ] + "_real," + ( *elemVNames )[ 0 ] + "_imag";
      for ( int j = 1; j < nVElem; ++j ) {
        vheader += "," + ( *elemVNames )[ j ] + "_real";
        vheader += "," + ( *elemVNames )[ j ] + "_imag";
      }
      file_vtu << " Vectors=\"" + vheader + "\"";
    }

    file_vtu << ">" << std::endl;

    for ( int j = 0; j < nElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *elemNames )[ j ] + "_real" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nElems; ++i ) {
        file_vtu << "          " << ( *elemData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\""
        + ( *elemNames )[ j ] + "_imag" + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < this->nElems; ++i ) {
        file_vtu << "          " << ( *elemData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    for ( int j = 0; j < nVElem; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *elemVNames )[ j ] + "_real" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < 3 * this->nElems; ++i ) {
        file_vtu << "          " << ( *elemVData )[ j ]->get( i ).real( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" +
        ( *elemVNames )[ j ] + "_imag" + "\" NumberOfComponents=\"3"
        + "\" format=\"ascii\">" << std::endl;
      for ( long i = 0; i < 3 * this->nElems; ++i ) {
        file_vtu << "          " << ( *elemVData )[ j ]->get( i ).imag( )
          << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }

    file_vtu << "      </CellData>" << std::endl;
  }

  file_vtu << "    </Piece>" << std::endl;
  file_vtu << "  </UnstructuredGrid>" << std::endl;
  file_vtu << "</VTKFile>" << std::endl;
  file_vtu.close( );

#ifdef VERBOSE
  std::cout << "done." << std::endl;
#endif

}

//template<class LO, class SC>
//void SurfaceMesh3D<LO, SC>::initCurl( ) {
//  LO nElems = this->getNElements( );
//  if ( auxCurl ) {
//    delete auxCurl;
//  }
//  auxCurl = new Vector<LO, SCVT>( 9 * nElems );
//  SCVT x1[3], x2[3], x3[3];
//  Vector<LO, SCVT> vx1( 3, x1 ), vx2( 3, x2 ), vx3( 3, x3 );
//  Vector<LO, SCVT> vn( 3 ), vx21( 3 ), vx31( 3 ), res( 9 );
//  FullMatrix<LO, SCVT> R( 3, 3 );
//  SCVT grad[ 9 ] = { -1.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
//  Vector<LO, SCVT> *vgrad;
//  for ( LO i = 0; i < nElems; i++ ) {
//    this->getNodes( i, x1, x2, x3 );
//    vx2.add( vx1, vx21, -1.0 );
//    vx3.add( vx1, vx31, -1.0 );
//    vx21.cross( vx31, vn );
//    vn.scale( 1.0 / vn.norm2( ) );
//
//    R.set( 0, 0, vx21.get( 0 ) );
//    R.set( 0, 1, vx21.get( 1 ) );
//    R.set( 0, 2, vx21.get( 2 ) );
//    R.set( 1, 0, vx31.get( 0 ) );
//    R.set( 1, 1, vx31.get( 1 ) );
//    R.set( 1, 2, vx31.get( 2 ) );
//    R.set( 2, 0, vn.get( 0 ) );
//    R.set( 2, 1, vn.get( 1 ) );
//    R.set( 2, 2, vn.get( 2 ) );
//
//    vgrad = new Vector<LO, SCVT>( 9, grad, true );
//    R.LUSolve( *vgrad, 3 );
//    vn.cross( *vgrad, res, 3 );
//
//    for ( int j = 0; j < 9; j++ ) {
//      auxCurl->set( j + 9 * i, (SCVT) res.get( j ) );
//    }
//    delete vgrad;
//  }
//  return;
//}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::initCurl( ) {

  LO nElems = this->getNElements( );

  if ( auxCurl ) {
    delete auxCurl;
  }

  auxCurl = new Vector<LO, SCVT>( 9 * nElems );

#pragma omp parallel
  {
    SCVT x1[3], x2[3], x3[3];
    Vector<LO, SCVT> vx1( 3, x1 ), vx2( 3, x2 ), vx3( 3, x3 );
    Vector<LO, SCVT> vn( 3 ), vx21( 3 ), vx31( 3 ), res( 9 );
    FullMatrix<LO, SCVT> R( 3, 3 );
    SCVT grad[ 9 ] = { -1.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
    Vector<LO, SCVT> *vgrad;
#pragma omp for
    for ( LO i = 0; i < nElems; i++ ) {
      this->getNodes( i, x1, x2, x3 );
      vx2.add( vx1, vx21, -1.0 );
      vx3.add( vx1, vx31, -1.0 );
      vx21.cross( vx31, vn );
      vn.scale( 1.0 / vn.norm2( ) );

      R.set( 0, 0, vx21.get( 0 ) );
      R.set( 0, 1, vx21.get( 1 ) );
      R.set( 0, 2, vx21.get( 2 ) );
      R.set( 1, 0, vx31.get( 0 ) );
      R.set( 1, 1, vx31.get( 1 ) );
      R.set( 1, 2, vx31.get( 2 ) );
      R.set( 2, 0, vn.get( 0 ) );
      R.set( 2, 1, vn.get( 1 ) );
      R.set( 2, 2, vn.get( 2 ) );

      vgrad = new Vector<LO, SCVT>( 9, grad, true );
      R.LUSolve( *vgrad, 3 );
      vn.cross( *vgrad, res, 3 );

      for ( int j = 0; j < 9; j++ ) {
        auxCurl->set( j + 9 * i, (SCVT) res.get( j ) );
      }
      delete vgrad;
    }
  }
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::assembleT(
  SparseMatrix< LO, SCVT > & T12,
  SparseMatrix< LO, SCVT > & T13,
  SparseMatrix< LO, SCVT > & T23,
  SparseMatrix< LO, SCVT > & T
  ) {

  LO nElems = this->getNElements( );
  LO nNodes = this->getNNodes( );

  std::vector< LO > rowInd;
  rowInd.reserve( 3 * nElems );
  std::vector< LO > colInd;
  colInd.reserve( 3 * nElems );
  std::vector< SCVT > values12;
  values12.reserve( 3 * nElems );
  std::vector< SCVT > values13;
  values13.reserve( 3 * nElems );
  std::vector< SCVT > values23;
  values23.reserve( 3 * nElems );
  LO elem[ 3 ];
  Vector< LO, SCVT > * curls = this->getCurls( );

  for ( LO i = 0; i < nElems; ++i ) {

    this->getElement( i, elem );

    for ( int node = 0; node < 3; ++node ) {
      rowInd.push_back( i );
      colInd.push_back( elem[ node ] );
      values23.push_back( -curls->get( 9 * i + 3 * node ) );
      values13.push_back( curls->get( 9 * i + 3 * node + 1 ) );
      values12.push_back( -curls->get( 9 * i + 3 * node + 2 ) );
    }
  }

  T12.setFromTriplets( nElems, nNodes, rowInd, colInd, values12 );
  T13.setFromTriplets( nElems, nNodes, rowInd, colInd, values13 );
  T23.setFromTriplets( nElems, nNodes, rowInd, colInd, values23 );

  rowInd.clear( );
  rowInd.reserve( 9 * nElems );
  colInd.clear( );
  colInd.reserve( 9 * nElems );
  std::vector<SCVT> values;
  values.reserve( 9 * nElems );

  for ( LO i = 0; i < nElems; ++i ) {

    this->getElement( i, elem );

    for ( int node = 0; node < 3; ++node ) {

      // write to matrix T12
      rowInd.push_back( i );
      colInd.push_back( elem[ node ] + nNodes );
      values.push_back( -curls->get( 9 * i + 3 * node + 2 ) );
      // write to matrix -T12
      rowInd.push_back( i + nElems );
      colInd.push_back( elem[ node ] );
      values.push_back( curls->get( 9 * i + 3 * node + 2 ) );

      // write to matrix T13
      rowInd.push_back( i );
      colInd.push_back( elem[ node ] + 2 * nNodes );
      values.push_back( curls->get( 9 * i + 3 * node + 1 ) );
      // write to matrix -T13
      rowInd.push_back( i + 2 * nElems );
      colInd.push_back( elem[ node ] );
      values.push_back( -curls->get( 9 * i + 3 * node + 1 ) );

      // write to matrix T23
      rowInd.push_back( i + nElems );
      colInd.push_back( elem[ node ] + 2 * nNodes );
      values.push_back( -curls->get( 9 * i + 3 * node ) );
      // write to matrix -T23
      rowInd.push_back( i + 2 * nElems );
      colInd.push_back( elem[ node ] + nNodes );
      values.push_back( curls->get( 9 * i + 3 * node ) );
    }
  }

  T.setFromTriplets( 3 * nElems, 3 * nNodes, rowInd, colInd, values );

}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::refine(
  int level,
  int type
  ) {

  switch ( type ) {
    case 2:
      this->refine_2sect( level );
      break;
    case 3:
      this->refine_3sect( level );
      break;
  }
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::refine_2sect(
  int level
  ) {

  LO newNNodes, newNElems, newNEdges;
  SCVT x1[ 3 ], x2[ 3 ];
  LO edge[ 2 ], edges[ 3 ], element[ 3 ];
  LO node1, node2, node3, node4, node5, node6;

  for ( int l = 0; l < level; l++ ) {

    // allocate new arrays
    newNNodes = this->nNodes + this->nEdges;
    newNElems = 4 * this->nElems;
    vector<SCVT> newNodes( this->nodes );
    newNodes.resize( 3 * newNNodes );
    vector<LO> newElems;
    newElems.resize( 3 * newNElems );

    newNEdges = 2 * this->nEdges + 3 * this->nElems;
    vector< LO > newEdges;
    newEdges.resize( 2 * newNEdges );

    vector< LO > newElems2Edges;
    newElems2Edges.resize( 3 * newNElems );

    // loop through edges to insert new nodes
#pragma omp parallel for private( x1, x2, edge )
    for ( LO i = 0; i < this->nEdges; i++ ) {
      this->getEdge( i, edge );
      this->getNode( edge[ 0 ], x1 );
      this->getNode( edge[ 1 ], x2 );
      newNodes[ 3 * ( this->nNodes + i ) ] = ( x1[ 0 ] + x2[ 0 ] ) / 2.0;
      newNodes[ 3 * ( this->nNodes + i ) + 1 ] = ( x1[ 1 ] + x2[ 1 ] ) / 2.0;
      newNodes[ 3 * ( this->nNodes + i ) + 2 ] = ( x1[ 2 ] + x2[ 2 ] ) / 2.0;
      newEdges[ 4 * i ] = edge[ 0 ];
      newEdges[ 4 * i + 1 ] = this->nNodes + i;
      newEdges[ 4 * i + 2 ] = edge[ 1 ];
      newEdges[ 4 * i + 3 ] = this->nNodes + i;
    }

    /* loop through old elements to define new elements
                      node1
                       /\
                      /  \
                     / 1  \
               node4 ______ node6
                   / \ 4  / \
                  / 2 \  / 3 \
                 /_____\/_____\
               node2  node5  node3
     */
#pragma omp parallel for \
private( element, edges, node1, node2, node3, node4, node5, node6 )
    for ( LO i = 0; i < this->nElems; i++ ) {
      this->getElement( i, element );
      this->getEdges( i, edges );

      node1 = element[ 0 ];
      node2 = element[ 1 ];
      node3 = element[ 2 ];
      node4 = this->nNodes + edges[ 0 ];
      node5 = this->nNodes + edges[ 1 ];
      node6 = this->nNodes + edges[ 2 ];

      // first element
      newElems[ 12 * i ] = node1;
      newElems[ 12 * i + 1 ] = node4;
      newElems[ 12 * i + 2 ] = node6;
      // second element
      newElems[ 12 * i + 3 ] = node4;
      newElems[ 12 * i + 4 ] = node2;
      newElems[ 12 * i + 5 ] = node5;
      // third element
      newElems[ 12 * i + 6 ] = node5;
      newElems[ 12 * i + 7 ] = node3;
      newElems[ 12 * i + 8 ] = node6;
      // fourth element
      newElems[ 12 * i + 9 ] = node4;
      newElems[ 12 * i + 10 ] = node5;
      newElems[ 12 * i + 11 ] = node6;

      if ( node4 < node5 ) {
        newEdges[ 4 * this->nEdges + 6 * i ] = node4;
        newEdges[ 4 * this->nEdges + 6 * i + 1 ] = node5;
      } else {
        newEdges[ 4 * this->nEdges + 6 * i ] = node5;
        newEdges[ 4 * this->nEdges + 6 * i + 1 ] = node4;
      }

      if ( node5 < node6 ) {
        newEdges[ 4 * this->nEdges + 6 * i + 2 ] = node5;
        newEdges[ 4 * this->nEdges + 6 * i + 3 ] = node6;
      } else {
        newEdges[ 4 * this->nEdges + 6 * i + 2 ] = node6;
        newEdges[ 4 * this->nEdges + 6 * i + 3 ] = node5;
      }

      if ( node4 < node6 ) {
        newEdges[ 4 * this->nEdges + 6 * i + 4 ] = node4;
        newEdges[ 4 * this->nEdges + 6 * i + 5 ] = node6;
      } else {
        newEdges[ 4 * this->nEdges + 6 * i + 4 ] = node6;
        newEdges[ 4 * this->nEdges + 6 * i + 5 ] = node4;
      }

      newElems2Edges[ 12 * i ] = 2 * edges[ 0 ];
      newElems2Edges[ 12 * i + 1 ] = 2 * this->nEdges + 3 * i + 2;
      newElems2Edges[ 12 * i + 2 ] = 2 * edges[ 2 ];

      newElems2Edges[ 12 * i + 3 ] = 2 * edges[ 0 ];
      newElems2Edges[ 12 * i + 4 ] = 2 * edges[ 1 ];
      newElems2Edges[ 12 * i + 5 ] = 2 * this->nEdges + 3 * i;

      newElems2Edges[ 12 * i + 6 ] = 2 * edges[ 1 ];
      newElems2Edges[ 12 * i + 7 ] = 2 * edges[ 2 ];
      newElems2Edges[ 12 * i + 8 ] = 2 * this->nEdges + 3 * i + 1;

      if ( node1 > node2 ) {
        ++newElems2Edges[ 12 * i ];
      } else {
        ++newElems2Edges[ 12 * i + 3 ];
      }

      if ( node1 > node3 ) {
        ++newElems2Edges[ 12 * i + 2 ];
      } else {
        ++newElems2Edges[ 12 * i + 7 ];
      }

      if ( node2 > node3 ) {
        ++newElems2Edges[ 12 * i + 4 ];
      } else {
        ++newElems2Edges[ 12 * i + 6 ];
      }

      newElems2Edges[ 12 * i + 9 ] = 2 * this->nEdges + 3 * i;
      newElems2Edges[ 12 * i + 10 ] = 2 * this->nEdges + 3 * i + 1;
      newElems2Edges[ 12 * i + 11 ] = 2 * this->nEdges + 3 * i + 2;
    }

    // update the mesh
    this->nNodes = newNNodes;
    this->nElems = newNElems;
    this->nodes.swap( newNodes );
    this->elems.swap( newElems );
    //initEdges( );
    this->nEdges = newNEdges;
    this->edges.swap( newEdges );
    this->elem2edges.swap( newElems2Edges );
  }

  initArea( );
  initNormals( );
  initLocalCoordinates( );
  initCurl( );
  initNode2Elems( );
  initL2g( );
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::refine_3sect(
  int level
  ) {

  LO newNNodes, newNElems, newNEdges;
  SCVT x1[ 3 ], x2[ 3 ];
  LO edge[ 2 ], edges[ 3 ], element[ 3 ];
  LO n1, n2, n3, n4, n5, n6, n7, n8, n9, n10;
  SCVT onethird = 1.0 / 3.0;
  SCVT twothirds = 2.0 / 3.0;

  for ( int l = 0; l < level; l++ ) {

    // allocate new arrays
    newNNodes = this->nNodes + 2 * this->nEdges + this->nElems;
    newNElems = 9 * this->nElems;
    vector<SCVT> newNodes( this->nodes );
    newNodes.resize( 3 * newNNodes );
    vector<LO> newElems;
    newElems.resize( 3 * newNElems );

    newNEdges = 3 * this->nEdges + 9 * this->nElems;
    vector< LO > newEdges;
    newEdges.resize( 2 * newNEdges );
    vector< LO > newElem2Edges;
    newElem2Edges.resize( 3 * newNElems );

    // loop through elements to insert centroids
#pragma omp parallel for private( x1 )
    for ( LO i = 0; i < this->nElems; ++i ) {
      this->getCentroid( i, x1 );
      newNodes[ 3 * ( this->nNodes + i ) ] = x1[ 0 ];
      newNodes[ 3 * ( this->nNodes + i ) + 1 ] = x1[ 1 ];
      newNodes[ 3 * ( this->nNodes + i ) + 2 ] = x1[ 2 ];
    }

    LO offset = 3 * ( this->nNodes + this->nElems );

    // loop through edges to insert new nodes
#pragma omp parallel for private( edge, x1, x2 )
    for ( LO i = 0; i < this->nEdges; ++i ) {
      this->getEdge( i, edge );
      this->getNode( edge[ 0 ], x1 );
      this->getNode( edge[ 1 ], x2 );
      newNodes[ offset + 6 * i ] = twothirds * x1[ 0 ] + onethird * x2[ 0 ];
      newNodes[ offset + 6 * i + 1 ] = twothirds * x1[ 1 ] + onethird * x2[ 1 ];
      newNodes[ offset + 6 * i + 2 ] = twothirds * x1[ 2 ] + onethird * x2[ 2 ];
      newNodes[ offset + 6 * i + 3 ] = onethird * x1[ 0 ] + twothirds * x2[ 0 ];
      newNodes[ offset + 6 * i + 4 ] = onethird * x1[ 1 ] + twothirds * x2[ 1 ];
      newNodes[ offset + 6 * i + 5 ] = onethird * x1[ 2 ] + twothirds * x2[ 2 ];
      newEdges[ 6 * i ] = edge[ 0 ];
      newEdges[ 6 * i + 1 ] = this->nNodes + this->nElems + 2 * i;
      newEdges[ 6 * i + 2 ] = this->nNodes + this->nElems + 2 * i;
      newEdges[ 6 * i + 3 ] = this->nNodes + this->nElems + 2 * i + 1;
      newEdges[ 6 * i + 4 ] = edge[ 1 ];
      newEdges[ 6 * i + 5 ] = this->nNodes + this->nElems + 2 * i + 1;
    }

    /* loop through old elements to define new elements
                       n1
                       /\
                      /  \
                     / 1  \
                   n2______n9
                   / \ 3  / \
                  / 2 \  / 4 \
                n3_____\/n10__n8
                /\     /\    / \
               /  \ 6 /  \ 8/   \
              /_5__\ /_7__\/__9__\
             n4    n5     n6     n7
     */
    offset = this->nNodes + this->nElems;
#pragma omp parallel for \
private( element, edges, edge, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10 )
    for ( LO i = 0; i < this->nElems; ++i ) {
      this->getElement( i, element );
      this->getEdges( i, edges );

      n1 = element[ 0 ];
      n2 = n3 = offset + 2 * edges[ 0 ];
      this->getEdge( edges[ 0 ], edge );
      ( edge[ 0 ] == n1 ) ? ( ++n3 ) : ( ++n2 );

      n4 = element[ 1 ];
      n5 = n6 = offset + 2 * edges[ 1 ];
      this->getEdge( edges[ 1 ], edge );
      ( edge[ 0 ] == n4 ) ? ( ++n6 ) : ( ++n5 );

      n7 = element[ 2 ];
      n8 = n9 = offset + 2 * edges[ 2 ];
      this->getEdge( edges[ 2 ], edge );
      ( edge[ 0 ] == n7 ) ? ( ++n9 ) : ( ++n8 );

      n10 = this->nNodes + i;

      // element 1
      newElems[ 27 * i ] = n1;
      newElems[ 27 * i + 1 ] = n2;
      newElems[ 27 * i + 2 ] = n9;
      // element 2
      newElems[ 27 * i + 3 ] = n3;
      newElems[ 27 * i + 4 ] = n10;
      newElems[ 27 * i + 5 ] = n2;
      // element 3
      newElems[ 27 * i + 6 ] = n2;
      newElems[ 27 * i + 7 ] = n10;
      newElems[ 27 * i + 8 ] = n9;
      // element 4
      newElems[ 27 * i + 9 ] = n10;
      newElems[ 27 * i + 10 ] = n8;
      newElems[ 27 * i + 11 ] = n9;
      // element 5
      newElems[ 27 * i + 12 ] = n4;
      newElems[ 27 * i + 13 ] = n5;
      newElems[ 27 * i + 14 ] = n3;
      // element 6
      newElems[ 27 * i + 15 ] = n3;
      newElems[ 27 * i + 16 ] = n5;
      newElems[ 27 * i + 17 ] = n10;
      // element 7
      newElems[ 27 * i + 18 ] = n5;
      newElems[ 27 * i + 19 ] = n6;
      newElems[ 27 * i + 20 ] = n10;
      // element 8
      newElems[ 27 * i + 21 ] = n10;
      newElems[ 27 * i + 22 ] = n6;
      newElems[ 27 * i + 23 ] = n8;
      // element 9
      newElems[ 27 * i + 24 ] = n6;
      newElems[ 27 * i + 25 ] = n7;
      newElems[ 27 * i + 26 ] = n8;

      if ( n2 < n9 ) {
        newEdges[ 6 * this->nEdges + 18 * i ] = n2;
        newEdges[ 6 * this->nEdges + 18 * i + 1 ] = n9;
      } else {
        newEdges[ 6 * this->nEdges + 18 * i ] = n9;
        newEdges[ 6 * this->nEdges + 18 * i + 1 ] = n2;
      }

      newEdges[ 6 * this->nEdges + 18 * i + 2 ] = n10;
      newEdges[ 6 * this->nEdges + 18 * i + 3 ] = n2;

      newEdges[ 6 * this->nEdges + 18 * i + 4 ] = n10;
      newEdges[ 6 * this->nEdges + 18 * i + 5 ] = n9;

      newEdges[ 6 * this->nEdges + 18 * i + 6 ] = n10;
      newEdges[ 6 * this->nEdges + 18 * i + 7 ] = n3;

      newEdges[ 6 * this->nEdges + 18 * i + 8 ] = n10;
      newEdges[ 6 * this->nEdges + 18 * i + 9 ] = n8;

      if ( n3 < n5 ) {
        newEdges[ 6 * this->nEdges + 18 * i + 10 ] = n3;
        newEdges[ 6 * this->nEdges + 18 * i + 11 ] = n5;
      } else {
        newEdges[ 6 * this->nEdges + 18 * i + 10 ] = n5;
        newEdges[ 6 * this->nEdges + 18 * i + 11 ] = n3;
      }

      newEdges[ 6 * this->nEdges + 18 * i + 12 ] = n10;
      newEdges[ 6 * this->nEdges + 18 * i + 13 ] = n5;

      newEdges[ 6 * this->nEdges + 18 * i + 14 ] = n10;
      newEdges[ 6 * this->nEdges + 18 * i + 15 ] = n6;

      if ( n6 < n8 ) {
        newEdges[ 6 * this->nEdges + 18 * i + 16 ] = n6;
        newEdges[ 6 * this->nEdges + 18 * i + 17 ] = n8;
      } else {
        newEdges[ 6 * this->nEdges + 18 * i + 16 ] = n8;
        newEdges[ 6 * this->nEdges + 18 * i + 17 ] = n6;
      }
      // element1
      newElem2Edges[ 27 * i ] = 3 * edges[ 0 ];
      newElem2Edges[ 27 * i + 1 ] = 3 * this->nEdges + 9 * i;
      newElem2Edges[ 27 * i + 2 ] = 3 * edges[ 2 ];
      // element2
      newElem2Edges[ 27 * i + 3 ] = 3 * this->nEdges + 9 * i + 3;
      newElem2Edges[ 27 * i + 4 ] = 3 * this->nEdges + 9 * i + 1;
      newElem2Edges[ 27 * i + 5 ] = 3 * edges[ 0 ] + 1;
      // element3
      newElem2Edges[ 27 * i + 6 ] = 3 * this->nEdges + 9 * i + 1;
      newElem2Edges[ 27 * i + 7 ] = 3 * this->nEdges + 9 * i + 2;
      newElem2Edges[ 27 * i + 8 ] = 3 * this->nEdges + 9 * i;
      // element4
      newElem2Edges[ 27 * i + 9 ] = 3 * this->nEdges + 9 * i + 4;
      newElem2Edges[ 27 * i + 10 ] = 3 * edges[ 2 ] + 1;
      newElem2Edges[ 27 * i + 11 ] = 3 * this->nEdges + 9 * i + 2;
      // element5
      newElem2Edges[ 27 * i + 12 ] = 3 * edges[ 1 ];
      newElem2Edges[ 27 * i + 13 ] = 3 * this->nEdges + 9 * i + 5;
      newElem2Edges[ 27 * i + 14 ] = 3 * edges[ 0 ];
      // element6
      newElem2Edges[ 27 * i + 15 ] = 3 * this->nEdges + 9 * i + 5;
      newElem2Edges[ 27 * i + 16 ] = 3 * this->nEdges + 9 * i + 6;
      newElem2Edges[ 27 * i + 17 ] = 3 * this->nEdges + 9 * i + 3;
      // element7
      newElem2Edges[ 27 * i + 18 ] = 3 * edges[ 1 ] + 1;
      newElem2Edges[ 27 * i + 19 ] = 3 * this->nEdges + 9 * i + 7;
      newElem2Edges[ 27 * i + 20 ] = 3 * this->nEdges + 9 * i + 6;
      // element8
      newElem2Edges[ 27 * i + 21 ] = 3 * this->nEdges + 9 * i + 7;
      newElem2Edges[ 27 * i + 22 ] = 3 * this->nEdges + 9 * i + 8;
      newElem2Edges[ 27 * i + 23 ] = 3 * this->nEdges + 9 * i + 4;
      // element9
      newElem2Edges[ 27 * i + 24 ] = 3 * edges[ 1 ];
      newElem2Edges[ 27 * i + 25 ] = 3 * edges[ 2 ];
      newElem2Edges[ 27 * i + 26 ] = 3 * this->nEdges + 9 * i + 8;

      if ( n1 > n4 ) {
        newElem2Edges[ 27 * i ] += 2;
      } else {
        newElem2Edges[ 27 * i + 14 ] += 2;
      }

      if ( n4 > n7 ) {
        newElem2Edges[ 27 * i + 12 ] += 2;
      } else {
        newElem2Edges[ 27 * i + 24 ] += 2;
      }

      if ( n1 > n7 ) {
        newElem2Edges[ 27 * i + 2 ] += 2;
      } else {
        newElem2Edges[ 27 * i + 25 ] += 2;
      }

    }

    // update the mesh
    this->nNodes = newNNodes;
    this->nElems = newNElems;
    this->nodes.swap( newNodes );
    this->elems.swap( newElems );
    //initEdges( );
    this->nEdges = newNEdges;
    this->edges.swap( newEdges );
    this->elem2edges.swap( newElem2Edges );
  }

  initArea( );
  initNormals( );
  initLocalCoordinates( );
  initCurl( );
  initNode2Elems( );
  initL2g( );
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::mapToUnitBall( ) {

  SCVT x[ 3 ];
  SCVT norm;


  for ( LO i = 0; i < this->nNodes; i++ ) {
    this->getNode( i, x );
    norm = std::sqrt( x[ 0 ] * x[ 0 ] + x[ 1 ] * x[ 1 ] + x[ 2 ] * x[ 2 ] );
    this->nodes[ 3 * i ] = x[ 0 ] / norm;
    this->nodes[ 3 * i + 1 ] = x[ 1 ] / norm;
    this->nodes[ 3 * i + 2 ] = x[ 2 ] / norm;
  }

  initArea( );
  initNormals( );
  initLocalCoordinates( );
  initCurl( );
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::append(
  const SurfaceMesh3D<LO, SC>& mesh
  ) {

  this->elems.reserve( 3 * ( this->nElems + mesh.nElems ) );
  this->nodes.reserve( 3 * ( this->nNodes + mesh.nNodes ) );

  for ( LO i = 0; i < mesh.nNodes; i++ ) {
    this->nodes.push_back( mesh.nodes[ 3 * i ] );
    this->nodes.push_back( mesh.nodes[ 3 * i + 1 ] );
    this->nodes.push_back( mesh.nodes[ 3 * i + 2 ] );
  }

  for ( LO i = 0; i < mesh.nElems; i++ ) {
    this->elems.push_back( mesh.elems[ 3 * i ] + this->nNodes );
    this->elems.push_back( mesh.elems[ 3 * i + 1 ] + this->nNodes );
    this->elems.push_back( mesh.elems[ 3 * i + 2 ] + this->nNodes );
  }

  this->nElems += mesh.nElems;
  this->nNodes += mesh.nNodes;

  //! todo can be only appended
  initEdges( );
  initArea( );
  initNormals( );
  initLocalCoordinates( );
  initCurl( );
  initNode2Elems( );
  initL2g( );
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::scale(
  SCVT f1,
  SCVT f2,
  SCVT f3
  ) {

  typename std::vector< SCVT >::iterator it;

#pragma omp parallel for
  for ( it = this->nodes.begin( ); it < this->nodes.end( ); it += 3 ) {
    *it *= f1;
    *( it + 1 ) *= f2;
    *( it + 2 ) *= f3;
  }

  this->initArea( );
  this->initNormals( );
  this->initLocalCoordinates( );
  this->initCurl( );
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::scaleAroundCentroid(
  SCVT f1,
  SCVT f2,
  SCVT f3
  ) {

  SCVT centroid[ 3 ];
  this->getCentroid( centroid );

  typename std::vector< SCVT >::iterator it;

#pragma omp parallel for
  for ( it = this->nodes.begin( ); it < this->nodes.end( ); it += 3 ) {
    *it = ( *it - centroid[ 0 ] ) * f1 + centroid[ 0 ];
    *( it + 1 ) = ( *( it + 1 ) - centroid[ 1 ] ) * f2 + centroid[ 1 ];
    *( it + 2 ) = ( *( it + 2 ) - centroid[ 2 ] ) * f3 + centroid[ 2 ];
  }

  this->initArea( );
  this->initNormals( );
  this->initLocalCoordinates( );
  this->initCurl( );
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::rotate(
  SCVT* v,
  SCVT alpha
  ) {

  SCVT norm = std::sqrt( v[0] * v[0] + v[1] * v[1] + v[2] * v[2] );
  SCVT u[ 3 ];
  u[ 0 ] = v[ 0 ] / norm;
  u[ 1 ] = v[ 1 ] / norm;
  u[ 2 ] = v[ 2 ] / norm;

  SCVT ca = std::cos( alpha );
  SCVT sa = std::sin( alpha );

  SCVT x[ 3 ];
  SCVT y[ 3 ];
  for ( LO i = 0; i < this->nNodes; ++i ) {
    this->getNode( i, x );
    y[ 0 ] = ( ca + u[0] * u[0] * ( 1.0 - ca ) ) * x[0] +
      ( u[0] * u[1] * ( 1.0 - ca ) - u[2] * sa ) * x[1] +
      ( u[0] * u[2] * ( 1.0 - ca ) + u[1] * sa ) * x[2];
    y[ 1 ] = ( u[0] * u[1] * ( 1.0 - ca ) + u[2] * sa ) * x[0] +
      ( ca + u[1] * u[1] * ( 1.0 - ca ) ) * x[1] +
      ( u[1] * u[2] * ( 1.0 - ca ) - u[0] * sa ) * x[2];
    y[ 2 ] = ( u[1] * u[2] * ( 1.0 - ca ) - u[1] * sa ) * x[0] +
      ( u[1] * u[2] * ( 1.0 - ca ) + u[0] * sa ) * x[1] +
      ( ca + u[2] * u[2] * ( 1.0 - ca ) ) * x[2];
    this->nodes[ 3 * i ] = y[ 0 ];
    this->nodes[ 3 * i + 1 ] = y[ 1 ];
    this->nodes[ 3 * i + 2 ] = y[ 2 ];
  }

  initNormals( );
  // todo: is it necessary? (data may change slightly after rotation)
  initArea( );
  initLocalCoordinates( );
  initCurl( );
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::moveCentroidToOrigin( ) {

  SCVT centroid[ 3 ];
  this->getCentroid( centroid );

  for ( LO i = 0; i < this->nNodes; i++ ) {
    this->nodes[ 3 * i ] -= centroid[ 0 ];
    this->nodes[ 3 * i + 1 ] -= centroid[ 1 ];
    this->nodes[ 3 * i + 2 ] -= centroid[ 2 ];
  }
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::move(
  SCVT x1,
  SCVT x2,
  SCVT x3
  ) {

  typename std::vector< SCVT >::iterator it;

#pragma omp parallel for
  for ( it = this->nodes.begin( ); it < this->nodes.end( ); it += 3 ) {
    *it += x1;
    *( it + 1 ) += x2;
    *( it + 2 ) += x3;
  }

}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::move( LO idx, SCVT x1, SCVT x2, SCVT x3 ) {
  this->nodes[ 3 * idx ] += x1;
  this->nodes[ 3 * idx + 1 ] += x2;
  this->nodes[ 3 * idx + 2 ] += x3;
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::expand( SCVT alpha ) {
  this->initNode2Elems( );
  SCVT n[ 3 ];

  for ( LO i = 0; i < this->nNodes; i++ ) {
    this->getNormalNodal( i, n );
    this->nodes[ 3 * i ] += alpha * n[ 0 ];
    this->nodes[ 3 * i + 1 ] += alpha * n[ 1 ];
    this->nodes[ 3 * i + 2 ] += alpha * n[ 2 ];
  }
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::initNode2Elems( ) {

  this->node2elems.clear( );
  this->node2elems.resize( this->nNodes );
  LO elem[ 3 ];

  for ( LO i = 0; i < this->nElems; ++i ) {
    this->getElement( i, elem );
    this->node2elems[ elem[ 0 ] ].push_back( i );
    this->node2elems[ elem[ 1 ] ].push_back( i );
    this->node2elems[ elem[ 2 ] ].push_back( i );
  }
}

template<class LO, class SC>
typename SurfaceMesh3D< LO, SC >::SCVT SurfaceMesh3D<LO, SC>::getArea( ) const {

  SCVT ret = 0.0;
  typename std::vector< SCVT >::const_iterator it;

#pragma omp parallel for reduction( + : ret )
  for ( it = this->area.begin( ); it < this->area.end( ); ++it ) {
    ret += *it;
  }

  return ret;
}

template<class LO, class SC>
typename SurfaceMesh3D< LO, SC >::SCVT
SurfaceMesh3D<LO, SC>::getVolume( ) const {

  SCVT ret = 0.0;
  SCVT xc[ 3 ], n[ 3 ];
  SCVT third = 1.0 / 3.0;

#pragma omp parallel for reduction( + : ret ) private( xc, n )
  for ( LO i = 0; i < this->nElems; ++i ) {
    this->getCentroid( i, xc );
    this->getNormal( i, n );
    ret += this->getElemArea( i ) * DOT3( xc, n ) * third;
  }

  return ret;
}

template<class LO, class SC>
typename bem4i::SurfaceMesh3D<LO, SC>::SCVT SurfaceMesh3D<LO, SC>::getCurvature(
  LO idx
  ) const {

  if ( this->additiveCurvature.size( ) > idx )
    return this->additiveCurvature[ idx ];
  return -1.0;
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::getCurvatureVector(
  Vector< LO, SC > & curvature
  ) const {

  LO size = this->additiveCurvature.size( );
  curvature.resize( size );
  for ( LO i = 0; i < size; ++i ) {
    curvature.set( i, this->additiveCurvature[ i ] );
  }
}

template<class LO, class SC>
void SurfaceMesh3D<LO, SC>::initAdditiveCurvature( ) {
  //this->initNode2Elems( ); // done at construction hopefully
  this->additiveCurvature.clear( );
  this->additiveCurvature.resize( this->nNodes );

#pragma omp parallel for
  for ( int i = 0; i < this->nNodes; ++i ) {
    this->additiveCurvature[ i ] = this->computeNodalAdditiveCurvature( i );
  }
}

template<class LO, class SC>
typename bem4i::SurfaceMesh3D<LO, SC>::SCVT SurfaceMesh3D<LO, SC>::
computeNodalAdditiveCurvature(
  LO nodeIdx
  ) {

  SCVT curvOp[ 3 ] = { 0.0, 0.0, 0.0 };
  SCVT meanNormal[ 3 ] = { 0.0, 0.0, 0.0 };
  SCVT mixedArea = 0.0;
  SCVT totalArea = 0.0;
  SCVT piHalf = M_PI / 2.0;

  LO elem[ 3 ];
  SCVT A[ 3 ], B[ 3 ], C[ 3 ];
  SCVT normal[ 3 ];

  SCVT a, b, c;
  SCVT beta, gamma;
  SCVT cotBeta, cotGamma;
  LO elemIdx;
  int aux;
  int obtuse; // 0 ... non-obtuse, 1 ... obtuse at A, 2 ... obtuse elsewhere

  SCVT area;

  for ( int i = 0; i < this->node2elems[ nodeIdx ].size( ); ++i ) {
    elemIdx = this->node2elems[ nodeIdx ][ i ];
    area = this->getElemArea( elemIdx );
    totalArea += area;
    this->getNormal( elemIdx, normal );
    meanNormal[ 0 ] += area * normal[ 0 ];
    meanNormal[ 1 ] += area * normal[ 1 ];
    meanNormal[ 2 ] += area * normal[ 2 ];
    this->getElement( elemIdx, elem );
    for ( aux = 0; aux < 3; ++aux ) {
      if ( nodeIdx == elem[ aux ] ) {
        break;
      }
    }
    this->getNode( elem[ aux % 3 ], A );
    this->getNode( elem[ ( aux + 1 ) % 3 ], B );
    this->getNode( elem[ ( aux + 2 ) % 3 ], C );
    a = std::sqrt( ( C[ 0 ] - B[ 0 ] )*( C[ 0 ] - B[ 0 ] )
      + ( C[ 1 ] - B[ 1 ] )*( C[ 1 ] - B[ 1 ] )
      + ( C[ 2 ] - B[ 2 ] )*( C[ 2 ] - B[ 2 ] ) );
    b = std::sqrt( ( C[ 0 ] - A[ 0 ] )*( C[ 0 ] - A[ 0 ] )
      + ( C[ 1 ] - A[ 1 ] )*( C[ 1 ] - A[ 1 ] )
      + ( C[ 2 ] - A[ 2 ] )*( C[ 2 ] - A[ 2 ] ) );
    c = std::sqrt( ( A[ 0 ] - B[ 0 ] )*( A[ 0 ] - B[ 0 ] )
      + ( A[ 1 ] - B[ 1 ] )*( A[ 1 ] - B[ 1 ] )
      + ( A[ 2 ] - B[ 2 ] )*( A[ 2 ] - B[ 2 ] ) );
    beta = std::acos( ( a * a + c * c - b * b ) / ( 2.0 * a * c ) );
    gamma = std::acos( ( b * b + a * a - c * c ) / ( 2.0 * b * a ) );

    obtuse = 0;
    if ( beta + gamma <= piHalf ) obtuse = 1;
    if ( beta > piHalf || gamma > piHalf ) obtuse = 2;
    cotBeta = std::tan( piHalf - beta );
    cotGamma = std::tan( piHalf - gamma );

    switch ( obtuse ) {
      case 0:
        mixedArea += ( 1.0 / 8.0 ) * ( c * c * cotGamma + b * b * cotBeta );
        break;
      case 1:
        mixedArea += area / 2.0;
        break;
      case 2:
        mixedArea += area / 4.0;
        break;
    }

    curvOp[ 0 ] += cotBeta * ( A[ 0 ] - C[ 0 ] )
      + cotGamma * ( A[ 0 ] - B[ 0 ] );
    curvOp[ 1 ] += cotBeta * ( A[ 1 ] - C[ 1 ] )
      + cotGamma * ( A[ 1 ] - B[ 1 ] );
    curvOp[ 2 ] += cotBeta * ( A[ 2 ] - C[ 2 ] )
      + cotGamma * ( A[ 2 ] - B[ 2 ] );
  }

  curvOp[ 0 ] /= 2.0 * mixedArea;
  curvOp[ 1 ] /= 2.0 * mixedArea;
  curvOp[ 2 ] /= 2.0 * mixedArea;
  meanNormal[ 0 ] /= totalArea;
  meanNormal[ 1 ] /= totalArea;
  meanNormal[ 2 ] /= totalArea; // need to normalize?

  SCVT sign = 1.0;
  SCVT dot = DOT3( curvOp, meanNormal );
  if ( dot < 0 ) sign = -1.0;

  SCVT curv = sign * std::sqrt( curvOp[ 0 ] * curvOp[ 0 ]
    + curvOp[ 1 ] * curvOp[ 1 ] + curvOp[ 2 ] * curvOp[ 2 ] );

  return curv;
}

template<class LO, class SC>
typename bem4i::SurfaceMesh3D<LO, SC>::SCVT SurfaceMesh3D<LO, SC>::
getRadius( ) {
  // compute centroid of a mesh

  SCVT cX = 0.0;
  SCVT cY = 0.0;
  SCVT cZ = 0.0;

  for ( int i = 0; i < this->getNNodes( ); i++ ) {
    cX += this->nodes[ 3 * i ];
    cY += this->nodes[ 3 * i + 1 ];
    cZ += this->nodes[ 3 * i + 2 ];
  }

  cX /= this->getNElements( );
  cY /= this->getNElements( );
  cZ /= this->getNElements( );

  // compute distance of each node to the centroid
  SCVT radius = 0.0;
  SCVT dist = 0.0;
  for ( int i = 0; i < this->getNNodes( ); i++ ) {
    dist = std::sqrt(
      ( cX - this->nodes[ 3 * i ] ) * ( cX - this->nodes[ 3 * i ] ) +
      ( cY - this->nodes[ 3 * i + 1 ] ) * ( cY - this->nodes[ 3 * i + 1 ] ) +
      ( cZ - this->nodes[ 3 * i + 2 ] ) * ( cZ - this->nodes[ 3 * i + 2] ) );
    if ( dist > radius ) {
      radius = dist;
    }
  }
  return radius;
}

} // end namespace bem4i

#endif
