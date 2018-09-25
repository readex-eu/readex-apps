/*!
 * @file    HomogenizationProblem.cpp 
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    December 2, 2014
 * 
 */

#ifdef HOMOGENIZATIONPROBLEM_H

namespace bem4i {

template<class LO, class SC>
HomogenizationProblem<LO, SC>::HomogenizationProblem( ) {

  std::fill( parameters, parameters + SETTINGS_SIZE, false );

  this->mesh = nullptr;
  this->meshCell = nullptr;
  this->meshInclusion = nullptr;
  this->P = nullptr;
  this->homoMatrix = new FullMatrix< LO, SC >( 3, 3, true );
  this->chi1 = nullptr;
  this->chi2 = nullptr;

  setDefaultParameters( );
}

template<class LO, class SC>
HomogenizationProblem<LO, SC>::HomogenizationProblem(
    const HomogenizationProblem & orig
    ) {
}

template<class LO, class SC>
HomogenizationProblem<LO, SC>::~HomogenizationProblem( ) {

  if ( this->mesh ) delete this->mesh;
  if ( this->meshCell ) delete this->meshCell;
  if ( this->P ) delete this->P;
  if ( this->homoMatrix ) delete this->homoMatrix;
  if ( this->chi1 ) delete this->chi1;
  if ( this->chi2 ) delete this->chi2;
}

//template<class LO, class SC>
//template< class argType >
//void HomogenizationProblem<LO, SC>::set(
//    Settings name,
//    argType & value
//    ) {
//
//  switch ( name ) {
//
//    case MESH_INCLUSION:
//      if ( typeid ( value ) != typeid ( SurfaceMesh3D< LO, SC > ) )
//        break;
//      this->meshInclusion = &value;
//      this->parameters[ MESH_INCLUSION ] = true;
//      break;
//
//    case N_CELL_REFINE:
//      if( value < 0 )
//        break;
//      this->nCellRefine = value;
//      this->parameters[ N_CELL_REFINE ] = true;
//      break;
//      break;
//
//    case SYSTEM_TYPE:
//      if( value < 0 || value >= SYSTEM_TYPE_SIZE )
//        break;
//      this->systemType = value;
//      this->parameters[ SYSTEM_TYPE ] = true;
//      break;
//      
//    default:
//      break;
//  }
//
//}

template<class LO, class SC>
void HomogenizationProblem<LO, SC>::setDefaultParameters( ) {

  this->systemType = HomogenizationProblem<LO, SC>::FULL_STEKLOV_POINCARE;
  this->nSegsPerEdge = 3;
  this->a1 = 1;
  this->a2 = 10;
  this->quadType = SauterSchwab;
  this->quadOrder = 3;
}

template<class LO, class SC>
bool HomogenizationProblem<LO, SC>::checkParameters( ) const {
  return true;
}

template<class LO, class SC>
bool HomogenizationProblem<LO, SC>::solve( ) {

  bool ret = false;
  this->prepareMesh( );

  if ( this->chi1 ) delete this->chi1;
  if ( this->chi2 ) delete this->chi2;
  this->chi1 = new FullMatrix<LO, SC>( this->meshInclusion->getNNodes( ), 3 );
  this->chi2 = new FullMatrix<LO, SC>( this->meshCell->getNNodes( ), 3 );
  
  switch ( this->systemType ) {
    case HomogenizationProblem<LO, SC>::FULL_STEKLOV_POINCARE:
      ret = this->solveFullSteklovPoincare( );
      break;
    default:
      break;
  }

  ProgressMonitor::init( "Assembling homogenized matrix" );
  this->assembleHomoMatrix( );
  ProgressMonitor::step( );

  return ret;
}

template<class LO, class SC>
void HomogenizationProblem<LO, SC>::assembleHomoMatrix( ) {

  this->homoMatrix->setAll( 0.0 );

  LO nElemsInclusion = this->meshInclusion->getNElements( );
  SCVT a;
  LO elem[ 3 ];
  SCVT n[ 3 ];
  SCVT centroid[ 3 ];

  for ( LO i = 0; i < nElemsInclusion; ++i ) {
    this->meshInclusion->getElement( i, elem );
    a = this->meshInclusion->getElemArea( i );
    this->meshInclusion->getNormal( i, n );
    this->meshInclusion->getCentroid( i, centroid );
    // invert normals!!!
    n[ 0 ] *= -1.0;
    n[ 1 ] *= -1.0;
    n[ 2 ] *= -1.0;

    SC chi1X1 = this->chi1->get( elem[ 0 ], 0 );
    SC chi1X2 = this->chi1->get( elem[ 1 ], 0 );
    SC chi1X3 = this->chi1->get( elem[ 2 ], 0 );
    SC chi1Mid = ( chi1X1 + chi1X2 + chi1X3 ) / 3.0;
    SC chi2X1 = this->chi1->get( elem[ 0 ], 1 );
    SC chi2X2 = this->chi1->get( elem[ 1 ], 1 );
    SC chi2X3 = this->chi1->get( elem[ 2 ], 1 );
    SC chi2Mid = ( chi2X1 + chi2X2 + chi2X3 ) / 3.0;
    SC chi3X1 = this->chi1->get( elem[ 0 ], 2 );
    SC chi3X2 = this->chi1->get( elem[ 1 ], 2 );
    SC chi3X3 = this->chi1->get( elem[ 2 ], 2 );
    SC chi3Mid = ( chi3X1 + chi3X2 + chi3X3 ) / 3.0;

    this->homoMatrix->add( 0, 0, n[ 0 ] * ( centroid[ 0 ] - chi1Mid ) * a );
    this->homoMatrix->add( 0, 1, n[ 0 ] * ( -chi2Mid ) * a );
    this->homoMatrix->add( 0, 2, n[ 0 ] * ( -chi3Mid ) * a );
    this->homoMatrix->add( 1, 0, n[ 1 ] * ( -chi1Mid ) * a );
    this->homoMatrix->add( 1, 1, n[ 1 ] * ( centroid[ 1 ] - chi2Mid ) * a );
    this->homoMatrix->add( 1, 2, n[ 1 ] * ( -chi3Mid ) * a );
    this->homoMatrix->add( 2, 0, n[ 2 ] * ( -chi1Mid ) * a );
    this->homoMatrix->add( 2, 1, n[ 2 ] * ( -chi2Mid ) * a );
    this->homoMatrix->add( 2, 2, n[ 2 ] * ( centroid[ 2 ] - chi3Mid ) * a );
  }


  homoMatrix->scale( this->a1 - this->a2 );

  homoMatrix->add( 0, 0, this->a2 );
  homoMatrix->add( 1, 1, this->a2 );
  homoMatrix->add( 2, 2, this->a2 );
}

template<class LO, class SC>
bool HomogenizationProblem<LO, SC>::solveFullSteklovPoincare( ) {

  LO nNodesInclusion = this->meshInclusion->getNNodes( );
  LO nReferenceNodes = this->P->getNCols( );
  LO dim = nNodesInclusion + nReferenceNodes;

  FullMatrix< LO, SC > * A = new FullMatrix< LO, SC >( dim, dim );
  FullMatrix< LO, SC > chi( dim, 3 );
  this->prepareFullSteklovPoincareMatrix( *A );
  this->prepareFullSteklovPoincareRHS( chi );
  ProgressMonitor::init( "Solving the system by LU" );
  A->LUSolve( chi, 3 );
  ProgressMonitor::step( );
  delete A;

  for ( LO i = 0; i < nNodesInclusion; ++i ) {
    this->chi1->set( i, 0, chi.get( i, 0 ) );
    this->chi1->set( i, 1, chi.get( i, 1 ) );
    this->chi1->set( i, 2, chi.get( i, 2 ) );
  }

  FullMatrix< LO, SC > chi2aux( nReferenceNodes, 3 );

  for ( LO i = 0; i < nReferenceNodes; ++i ) {
    chi2aux.set( i, 0, chi.get( nNodesInclusion + i, 0 ) );
    chi2aux.set( i, 1, chi.get( nNodesInclusion + i, 1 ) );
    chi2aux.set( i, 2, chi.get( nNodesInclusion + i, 2 ) );
  }

  this->chi2->multiply( *this->P, chi2aux, false, false, 1.0, 0.0 );

  return true;
}

template<class LO, class SC>
void HomogenizationProblem<LO, SC>::prepareMesh( ) {

  std::stringstream file;
  file << "input/homo/homocube_" << this->nSegsPerEdge << ".txt";
  if ( this->meshCell ) delete this->meshCell;
  this->meshCell = new SurfaceMesh3D<LO, SC>( file.str( ) );

  file.str( std::string( ) );
  file.clear( );
  file << "input/homo/Psparse_" << this->nSegsPerEdge << ".txt";
  if ( this->P ) delete this->P;
  this->P = new SparseMatrix<LO, SC>( 0, 0 );
  this->P->loadTriplets( file.str( ) );

  if ( this->mesh ) delete this->mesh;
  this->mesh = new SurfaceMesh3D< LO, SC >( *this->meshInclusion );
  this->mesh->append( *this->meshCell );
  
  std::cout << "Inclusion, cell, total:" << std::endl;
  this->meshInclusion->printInfo();
  this->meshCell->printInfo();
  this->mesh->printInfo();
}

template<class LO, class SC>
void HomogenizationProblem<LO, SC>::assembleIdentity(
    SparseMatrix<LO, SC> & M,
    Vector<LO, SC> & a
    ) {

  LO nElems = this->mesh->getNElements( );
  LO nNodes = this->mesh->getNNodes( );

  std::vector< LO > rowInd;
  rowInd.reserve( 3 * nElems );
  std::vector< LO > colInd;
  colInd.reserve( 3 * nElems );
  std::vector< SC > val;
  val.reserve( 3 * nElems );

  a.resize( nNodes, true );

  LO elem[ 3 ];
  SCVT integral;

  for ( LO i = 0; i < nElems; ++i ) {
    this->mesh->getElement( i, elem );
    integral = this->mesh->getElemArea( i ) / 3.0;
    a.add( elem[ 0 ], integral );
    a.add( elem[ 1 ], integral );
    a.add( elem[ 2 ], integral );
    rowInd.push_back( i );
    rowInd.push_back( i );
    rowInd.push_back( i );
    colInd.push_back( elem[ 0 ] );
    colInd.push_back( elem[ 1 ] );
    colInd.push_back( elem[ 2 ] );
    val.push_back( integral );
    val.push_back( integral );
    val.push_back( integral );
  }

  M.setFromTriplets( nElems, nNodes, rowInd, colInd, val );
}

template<class LO, class SC>
void HomogenizationProblem<LO, SC>::assembleIdentityG1G1(
    SparseMatrix<LO, SC> & MG1G1
    ) {

  LO nElems = this->meshInclusion->getNElements( );
  LO nNodes = this->meshInclusion->getNNodes( );

  std::vector< LO > rowInd;
  rowInd.reserve( 3 * nElems );
  std::vector< LO > colInd;
  colInd.reserve( 3 * nElems );
  std::vector< SC > val;
  val.reserve( 3 * nElems );

  LO elem[ 3 ];
  SCVT integral;

  for ( LO i = 0; i < nElems; ++i ) {
    this->mesh->getElement( i, elem );
    integral = this->mesh->getElemArea( i ) / 3.0;
    rowInd.push_back( i );
    rowInd.push_back( i );
    rowInd.push_back( i );
    colInd.push_back( elem[ 0 ] );
    colInd.push_back( elem[ 1 ] );
    colInd.push_back( elem[ 2 ] );
    val.push_back( integral );
    val.push_back( integral );
    val.push_back( integral );
  }

  MG1G1.setFromTriplets( nElems, nNodes, rowInd, colInd, val );
}

template<class LO, class SC>
void HomogenizationProblem<LO, SC>::prepareFullSteklovPoincareMatrix(
    FullMatrix<LO, SC> & A
    ) {

  int * quadOrderArray = nullptr;
  switch ( this->quadType ) {
    case SauterSchwab:
      quadOrderArray = new int[ 4 ];
      quadOrderArray[ 0 ] = quadOrderArray[ 1 ] = quadOrderArray[ 2 ] =
          quadOrderArray[ 3 ] = this->quadOrder;
      break;
    case Steinbach:
      quadOrderArray = new int[ 2 ];
      quadOrderArray[ 0 ] = quadOrderArray[ 1 ] = this->quadOrder;
      break;
  }

  LO nNodes = this->mesh->getNNodes( );
  LO nElemsInclusion = this->meshInclusion->getNElements( );
  LO nNodesInclusion = this->meshInclusion->getNNodes( );
  LO nNodesCell = this->meshCell->getNNodes( );
  LO nNodesReference = this->P->getNCols( );

  BESpace< LO, SC > bespaceP0P0( this->mesh, p0, p0 );
  BESpace< LO, SC > bespaceP1P0( this->mesh, p1, p0 );
  BESpace< LO, SC > bespaceP1P1( this->mesh, p1, p1 );

  // setting up BEM matrices for the whole mesh
  FullMatrix< LO, SC > * V = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLaplace1Layer< LO, SC > formV( &bespaceP0P0, quadOrderArray,
      this->quadType );
  formV.assemble( *V );

  FullMatrix< LO, SC > * K = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLaplace2Layer< LO, SC > formK( &bespaceP1P0, quadOrderArray,
      this->quadType );
  formK.assemble( *K );

  ProgressMonitor::init( "Setting up identity." );
  SparseMatrix< LO, SC > M;
  Vector< LO, SC > a;
  this->assembleIdentity( M, a );
  ProgressMonitor::step( );

  // assembling S1
  ProgressMonitor::init( "Assembling S1" );
  SparseMatrix< LO, SC > MG1G1;
  this->assembleIdentityG1G1( MG1G1 );

  FullMatrix< LO, SC > * KG1G1 = new FullMatrix< LO, SC >( nElemsInclusion,
      nNodesInclusion );
  this->copyIntoSubmatrix( *K, 0, nElemsInclusion, 0, nNodesInclusion, *KG1G1 );
  KG1G1->scale( -1.0 );
  KG1G1->add( MG1G1, 0.5 );

  FullMatrix< LO, SC > * VG1G1 = new FullMatrix< LO, SC >( nElemsInclusion,
      nElemsInclusion );
  this->copyIntoSubmatrix( *V, 0, nElemsInclusion, 0, nElemsInclusion, *VG1G1 );

  FullMatrix< LO, SC > * VinvKG1G1 = new FullMatrix< LO, SC >( *KG1G1 );
  VG1G1->LUSolve( *VinvKG1G1, nNodesInclusion );
  delete VG1G1;

  FullMatrix< LO, SC > * S1 = new FullMatrix< LO, SC >( nNodesInclusion,
      nNodesInclusion );
  S1->multiply( *KG1G1, *VinvKG1G1, true, false, 1.0, 0.0 );
  delete VinvKG1G1;
  delete KG1G1;
  ProgressMonitor::step( );

  // assembling S2
  FullMatrix< LO, SC > * D = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLaplaceHypersingular< LO, SC > formD( &bespaceP1P1,
      quadOrderArray, this->quadType );
  formD.assemble( *D, *V );

  ProgressMonitor::init( "Assembling S2" );
  K->add( M, 0.5 );
  FullMatrix< LO, SC > * VinvK = new FullMatrix< LO, SC >( *K );
  V->LUSolve( *VinvK, nNodes );
  delete V;

  FullMatrix< LO, SC > * S2 = new FullMatrix< LO, SC >( nNodes, nNodes );
  S2->multiply( *K, *VinvK, true, false, 1.0, 0.0 );
  delete VinvK;
  delete K;

  // add hypersingular parts
  this->addLeftUpper( *D, nNodesInclusion, nNodesInclusion, *S1 );
  S2->add( *D, 1.0 );
  delete D;

  // modify S1, S2
  S1->scale( this->a1 );
  S2->scale( this->a2 );
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nNodes; ++i ) {
      S2->add( i, j, a.get( i ) * a.get( j ) );
    }
  }
  ProgressMonitor::step( );

  // assemble G1G1
  ProgressMonitor::init( "Assembling system matrix G1G1" );
  for ( LO j = 0; j < nNodesInclusion; ++j ) {
    for ( LO i = 0; i < nNodesInclusion; ++i ) {
      A.set( i, j, S1->get( i, j ) + S2->get( i, j ) );
    }
  }
  ProgressMonitor::step( );

  // assemble G1G
  ProgressMonitor::init( "Assembling system matrix G1G" );
  FullMatrix< LO, SC > * S2G1G = new FullMatrix< LO, SC >( nNodesInclusion,
      nNodesCell );
  this->copyIntoSubmatrix( *S2, 0, nNodesInclusion, nNodesInclusion, nNodes,
      *S2G1G );
  FullMatrix< LO, SC > * S2G1GP = new FullMatrix< LO, SC >( nNodesInclusion,
      nNodesReference, true );
  S2G1GP->multiply( *S2G1G, *this->P, false, false, 1.0, 0.0 );
  this->copyFromSubmatrix( A, 0, nNodesInclusion, *S2G1GP );
  delete S2G1G;
  delete S2G1GP;
  ProgressMonitor::step( );

  // assemble GG1
  ProgressMonitor::init( "Assembling system matrix GG1" );
  FullMatrix< LO, SC > * S2GG1 = new FullMatrix< LO, SC >( nNodesCell,
      nNodesInclusion );
  this->copyIntoSubmatrix( *S2, nNodesInclusion, nNodes, 0, nNodesInclusion,
      *S2GG1 );
  FullMatrix< LO, SC > * PS2GG1 = new FullMatrix< LO, SC >( nNodesReference,
      nNodesInclusion, true );
  PS2GG1->multiply( *this->P, *S2GG1, true, false, 1.0, 0.0 );
  this->copyFromSubmatrix( A, nNodesInclusion, 0, *PS2GG1 );
  delete S2GG1;
  delete PS2GG1;
  ProgressMonitor::step( );

  // assemble GG
  ProgressMonitor::init( "Assembling system matrix GG" );
  FullMatrix< LO, SC > * S2GG = new FullMatrix< LO, SC >( nNodesCell,
      nNodesCell );
  this->copyIntoSubmatrix( *S2, nNodesInclusion, nNodes, nNodesInclusion,
      nNodes, *S2GG );
  FullMatrix< LO, SC > * PS2GG = new FullMatrix< LO, SC >( nNodesReference,
      nNodesCell, true );
  PS2GG->multiply( *this->P, *S2GG, true, false, 1.0, 0.0 );
  delete S2GG;
  FullMatrix< LO, SC > * PS2GGP = new FullMatrix< LO, SC >( nNodesReference,
      nNodesReference, true );
  PS2GGP->multiply( *PS2GG, *this->P, false, false, 1.0, 0.0 );
  this->copyFromSubmatrix( A, nNodesInclusion, nNodesInclusion, *PS2GGP );
  delete PS2GG;
  delete PS2GGP;
  ProgressMonitor::step( );

  delete S1;
  delete S2;
  if ( quadOrderArray ) delete [] quadOrderArray;
}

template<class LO, class SC>
void HomogenizationProblem<LO, SC>::prepareFullSteklovPoincareRHS(
    FullMatrix<LO, SC>& rhs
    ) {

  ProgressMonitor::init( "Assembling rhs" );
  rhs.setAll( 0.0 );

  LO nElemsInclusion = this->meshInclusion->getNElements( );
  LO elem[ 3 ];
  SCVT n[ 3 ];
  SCVT a;

  for ( LO i = 0; i < nElemsInclusion; ++i ) {
    this->meshInclusion->getElement( i, elem );
    this->meshInclusion->getNormal( i, n );
    a = this->meshInclusion->getElemArea( i ) / 3.0;
    // invert normal!!!
    n[ 0 ] *= -a;
    n[ 1 ] *= -a;
    n[ 2 ] *= -a;

    rhs.add( elem[ 0 ], 0, n[ 0 ] );
    rhs.add( elem[ 1 ], 0, n[ 0 ] );
    rhs.add( elem[ 2 ], 0, n[ 0 ] );
    rhs.add( elem[ 0 ], 1, n[ 1 ] );
    rhs.add( elem[ 1 ], 1, n[ 1 ] );
    rhs.add( elem[ 2 ], 1, n[ 1 ] );
    rhs.add( elem[ 0 ], 2, n[ 2 ] );
    rhs.add( elem[ 1 ], 2, n[ 2 ] );
    rhs.add( elem[ 2 ], 2, n[ 2 ] );
  }

  rhs.scale( this->a1 - this->a2 );
  ProgressMonitor::step( );
}

template<class LO, class SC>
void HomogenizationProblem<LO, SC>::copyIntoSubmatrix(
    const FullMatrix<LO, SC>& M,
    LO rowsStart,
    LO rowsEnd,
    LO colsStart,
    LO colsEnd,
    FullMatrix<LO, SC>& MG1G1
    ) {

  for ( LO j = 0; j < colsEnd - colsStart; ++j ) {
    for ( LO i = 0; i < rowsEnd - rowsStart; ++i ) {
      MG1G1.set( i, j, M.get( rowsStart + i, colsStart + j ) );
    }
  }

}

template<class LO, class SC>
void HomogenizationProblem<LO, SC>::copyFromSubmatrix(
    FullMatrix<LO, SC>& M,
    LO rowsStart,
    LO colsStart,
    const FullMatrix<LO, SC>& MSub
    ) {

  for ( LO j = 0; j < MSub.getNCols( ); ++j ) {
    for ( LO i = 0; i < MSub.getNRows( ); ++i ) {
      M.set( rowsStart + i, colsStart + j, MSub.get( i, j ) );
    }
  }

}

template<class LO, class SC>
void HomogenizationProblem<LO, SC>::addLeftUpper(
    const FullMatrix<LO, SC>& M,
    LO nRows,
    LO nCols,
    FullMatrix<LO, SC>& MG1G1
    ) {

  for ( LO j = 0; j < nCols; ++j ) {
    for ( LO i = 0; i < nRows; ++i ) {
      MG1G1.add( i, j, M.get( i, j ) );
    }
  }

}

template<class LO, class SC>
void HomogenizationProblem<LO, SC>::printVtu(
    const std::string& fileName
    ) {

  if ( this->mesh && this->meshInclusion && this->meshCell ) {

    std::vector< string > nodeNames;
    std::vector< string > elemNames;
    std::vector< Vector< LO, SC > * > nodalData;
    std::vector< Vector< LO, SC > * > elemData;

    LO nNodesInclusion = this->meshInclusion->getNNodes( );
    LO nNodesCell = this->meshCell->getNNodes( );
    Vector<LO, SC> chi_1( nNodesInclusion + nNodesCell );
    Vector<LO, SC> chi_2( nNodesInclusion + nNodesCell );
    Vector<LO, SC> chi_3( nNodesInclusion + nNodesCell );

    if ( this->chi1 && this->chi2 ) {
      for ( LO i = 0; i < nNodesInclusion; ++i ) {
        chi_1.set( i, this->chi1->get( i, 0 ) );
        chi_2.set( i, this->chi1->get( i, 1 ) );
        chi_3.set( i, this->chi1->get( i, 2 ) );
      }

      for ( LO i = 0; i < nNodesCell; ++i ) {
        chi_1.set( i + nNodesInclusion, this->chi2->get( i, 0 ) );
        chi_2.set( i + nNodesInclusion, this->chi2->get( i, 1 ) );
        chi_3.set( i + nNodesInclusion, this->chi2->get( i, 2 ) );
      }

      nodeNames.push_back( "chi1" );
      nodalData.push_back( &chi_1 );
      nodeNames.push_back( "chi2" );
      nodalData.push_back( &chi_2 );
      nodeNames.push_back( "chi3" );
      nodalData.push_back( &chi_3 );
    }

    this->mesh->printParaviewVtu( fileName, &nodeNames, &nodalData,
        &elemNames, &elemData );
  }
}

}
#endif
