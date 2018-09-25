/*!
 * @file    VolumeMesh3D.cpp
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    November 6, 2013
 * 
 */

#ifdef VOLUMEMESH3D_H

namespace bem4i {

template<class LO, class SC>
VolumeMesh3D<LO, SC>::VolumeMesh3D( ) {
  this->dim = 3;
  this->nNodes = 0;
  this->nElems = 0;
}

template<class LO, class SC>
VolumeMesh3D<LO, SC>::VolumeMesh3D( const VolumeMesh3D& orig ) {
}

template<class LO, class SC>
VolumeMesh3D<LO, SC>::~VolumeMesh3D( ) {
}

template<class LO, class SC>
void VolumeMesh3D<LO, SC>::load( const string& meshFile, SCVT scaleFactor ) {

  std::cout << "Reading file " << meshFile << " ... ";
  std::ifstream file( meshFile.c_str( ) );

  if ( !file.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
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

  std::cout << "done." << std::endl;
}

template<class LO, class SC>
void VolumeMesh3D<LO, SC>::printInfo( ) {
  std::cout << "Mesh info:" << std::endl;
  std::cout << "  class: VolumeMesh3D,";
  std::cout << " dim: " << this->dim;
  std::cout << ", nodes: " << this->nNodes;
  std::cout << ", elements: " << this->nElems;
  std::cout << "." << std::endl;
}

template<class LO, class SC>
void VolumeMesh3D<LO, SC>::printParaviewVtu( const string& meshFile ) {
  std::cout << "Printing  " << meshFile << " ... ";

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 10 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">" << std::endl;
  file_vtu << "  <UnstructuredGrid>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << this->nNodes <<
      "\" NumberOfCells=\"" << this->nElems << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  for ( LO i = 0; i < this->nNodes; i++ ) {
    file_vtu << "          "
        << this->nodes[ 3 * i ] << " "
        << this->nodes[ 3 * i + 1 ] << " "
        << this->nodes[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Cells>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;

  for ( LO i = 0; i < this->nElems; i++ ) {
    file_vtu << "          "
        << this->elems[ 4 * i ] << " "
        << this->elems[ 4 * i + 1 ] << " "
        << this->elems[ 4 * i + 2 ] << " "
        << this->elems[ 4 * i + 3 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;

  for ( LO offset = 1; offset <= this->nElems; offset++ ) {
    file_vtu << "          " << offset * 4 << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
  for ( LO i = 1; i <= this->nElems; i++ ) {
    file_vtu << "          10" << std::endl;
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
void VolumeMesh3D<LO, SC>::printParaviewVtu( const string& meshFile, int nNodal, string* nodeNames, Vector< LO, SC >** nodalData, int nElem, string* elemNames, Vector< LO, SC >** elemData ) {
  std::cout << "Printing  " << meshFile << " ... ";

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 10 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">" << std::endl;
  file_vtu << "  <UnstructuredGrid>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << this->nNodes <<
      "\" NumberOfCells=\"" << this->nElems << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  for ( LO i = 0; i < this->nNodes; i++ ) {
    file_vtu << "          "
        << this->nodes[ 3 * i ] << " "
        << this->nodes[ 3 * i + 1 ] << " "
        << this->nodes[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Cells>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;

  for ( LO i = 0; i < this->nElems; i++ ) {
    file_vtu << "          "
        << this->elems[ 4 * i ] << " "
        << this->elems[ 4 * i + 1 ] << " "
        << this->elems[ 4 * i + 2 ] << " "
        << this->elems[ 4 * i + 3 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;

  for ( LO offset = 1; offset <= this->nElems; offset++ ) {
    file_vtu << "          " << offset * 4 << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
  for ( LO i = 1; i <= this->nElems; i++ ) {
    file_vtu << "          10" << std::endl;
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
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" + nodeNames[ j ] + "\" format=\"ascii\">" << std::endl;
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
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" + elemNames[ j ] + "\" format=\"ascii\">" << std::endl;
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

  std::cout << "done." << std::endl;
}

} // end of namespace bem4i

#endif
