/*!
 * @file    RepresentationFormulaLaplace.cpp
 * @author  Michal Merta 
 * @date    November 5, 2013
 * 
 */

#ifdef REPRESENTATIONFORMULALAPLACE_H

namespace bem4i {

template<class LO, class SC>
RepresentationFormulaLaplace<LO, SC>::RepresentationFormulaLaplace( ) {
}

template<class LO, class SC>
RepresentationFormulaLaplace<LO, SC>::RepresentationFormulaLaplace(
    const RepresentationFormulaLaplace& orig
    ) {
}

template<class LO, class SC>
RepresentationFormulaLaplace<LO, SC>::RepresentationFormulaLaplace(
    BESpace<LO, SC>* space,
    Vector<LO, SCVT> *dirichlet,
    Vector<LO, SCVT>* neumann,
    int quadOrder
    ) {

  this->space = space;
  this->dirichlet = dirichlet;
  this->neumann = neumann;
  this->quadOrder = quadOrder;
}

template<class LO, class SC>
RepresentationFormulaLaplace<LO, SC>::~RepresentationFormulaLaplace( ) {
}

template<class LO, class SC>
SC RepresentationFormulaLaplace<LO, SC>::evaluate(
    const SCVT *x
    ) const {

  return 0.0;
}

template<class LO, class SC>
void RepresentationFormulaLaplace<LO, SC>::evaluate(
    const SCVT *x,
    LO n,
    bool interior,
    Vector<LO, SC> & values
    ) const {

  int qOrder[ 2 ] = { quadOrder, quadOrder };

  int mpiSize = 1;
  int mpiRank = 0;
  int mpiInitialized = 0;
  MPI_Initialized( &mpiInitialized );
  if ( mpiInitialized != 0 ) {
    MPI_Comm_size( MPI_COMM_WORLD, &mpiSize );
    MPI_Comm_rank( MPI_COMM_WORLD, &mpiRank );
  }

  int *numElems = new int[mpiSize];
  int *displs = new int[mpiSize];
  displs[0] = 0;
  int nElemsPerMPI = n / mpiSize;
  for ( int i = 0; i < mpiSize - 1; i++ ) {
    numElems[i] = nElemsPerMPI;
    displs[i + 1] = displs[i] + numElems[i];
  }
  numElems[mpiSize - 1] = nElemsPerMPI + ( n % mpiSize );

  SC *valuesPerMPI = new SC[numElems[mpiRank]];

#pragma omp parallel
  {
    int ompSize = omp_get_num_threads( );
    int ompRank = omp_get_thread_num( );
    LO nLocalElems = numElems[mpiRank] / ompSize;
    LO myElems = nLocalElems;
    if ( ompRank == ompSize - 1 ) {
      myElems = nLocalElems + numElems[mpiRank] % ompSize;
    }

    BEIntegratorLaplace<LO, SC> integrator( this->space, qOrder, Steinbach );
    Vector<LO, SC> localValues( myElems );
    integrator.representationFormula( x + 3 * ( displs[mpiRank] +
        ompRank * nLocalElems ), myElems, *this->dirichlet, *this->neumann,
        interior, localValues );
    for ( LO i = 0; i < myElems; ++i ) {
      valuesPerMPI[ompRank * nLocalElems + i] = localValues.get( i );
    }
  }

  if ( mpiInitialized != 0 ) {
    MPI_Datatype MPI_SC = GetType<LO, SC>::MPI_SC( );
    MPI_Allgatherv( valuesPerMPI, numElems[mpiRank],
        MPI_SC, values.getData( ), numElems, displs,
        MPI_SC, MPI_COMM_WORLD );
  } else {
    memcpy( values.getData( ), valuesPerMPI, sizeof (SC ) * numElems[0] );
  }

  delete [] numElems;
  delete [] displs;
  delete [] valuesPerMPI;
}

template<class LO, class SC>
void RepresentationFormulaLaplace<LO, SC>::evaluate(
    Mesh<LO, SC> &mesh,
    bool interior,
    Vector<LO, SC> & values
    ) const {

  SCVT *nodes = &( ( *( mesh.getNodes( ) ) )[0] );
  LO nNodes = mesh.getNNodes( );

  this->evaluate( nodes, nNodes, interior, values );
}

template<class LO, class SC>
void RepresentationFormulaLaplace<LO, SC>::printParaviewVtu(
    const string& meshFile,
    SC* points,
    LO nPoints,
    int nNodal,
    string* nodeNames,
    Vector< LO, SC >** nodalData
    ) const {

  std::cout << "Printing  " << meshFile << " ... ";

  std::ofstream file_vtu( meshFile.c_str( ) );

  file_vtu.setf( std::ios::showpoint | std::ios::scientific );
  file_vtu.precision( 10 );

  if ( !file_vtu.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file_vtu << "<?xml version=\"1.0\"?>" << std::endl;
  file_vtu << "<VTKFile type=\"PolyData\" version=\"0.1\">" << std::endl;
  file_vtu << "  <PolyData>" << std::endl;

  file_vtu << "    <Piece NumberOfPoints=\"" << nPoints <<
      "\" NumberOfVerts=\"" << nPoints << "\">" << std::endl;

  file_vtu << "      <Points>" << std::endl;
  file_vtu << "        <DataArray type=\"Float32\" Name=\"points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

  for ( LO i = 0; i < nPoints; i++ ) {
    file_vtu << "          "
        << points[ 3 * i ] << " "
        << points[ 3 * i + 1 ] << " "
        << points[ 3 * i + 2 ] << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Points>" << std::endl;
  file_vtu << "      <Verts>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;

  for ( LO i = 0; i < nPoints; i++ ) {
    file_vtu << "          " << i << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;

  for ( LO offset = 1; offset <= nPoints; offset++ ) {
    file_vtu << "          " << offset << std::endl;
  }

  file_vtu << "        </DataArray>" << std::endl;
  file_vtu << "      </Verts>" << std::endl;

  string header;
  if ( nNodal > 0 ) {
    string header = nodeNames[ 0 ];
    for ( int j = 1; j < nNodal; j++ ) {
      header += "," + nodeNames[ j ];
    }
    file_vtu << "      <PointData Scalars=\"" + header + "\">" << std::endl;
    for ( int j = 0; j < nNodal; j++ ) {
      file_vtu << "        <DataArray type=\"Float32\" Name=\"" + nodeNames[ j ] + "\" format=\"ascii\">" << std::endl;
      for ( LO i = 0; i < nPoints; i++ ) {
        file_vtu << "          " << nodalData[ j ]->get( i ) << std::endl;
      }
      file_vtu << "        </DataArray>" << std::endl;
    }
    file_vtu << "      </PointData>" << std::endl;
  }

  file_vtu << "    </Piece>" << std::endl;
  file_vtu << "  </PolyData>" << std::endl;
  file_vtu << "</VTKFile>" << std::endl;
  file_vtu.close( );

  std::cout << "done." << std::endl;
}

}

#endif
