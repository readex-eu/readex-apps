/*!
 * @file    WaveScatteringProblem.cpp
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    June 10, 2014
 * 
 */

#ifdef WAVESCATTERINGPROBLEM_H

namespace bem4i {

template<class LO, class SC>
WaveScatteringProblem<LO, SC>::WaveScatteringProblem( ) {
  parameters = new bool[this->maxParameters];
  std::fill( parameters, parameters + this->maxParameters, 0 );
  spatialQuad = new int[4];
  setDefaultParameters( );
}

template<class LO, class SC>
WaveScatteringProblem<LO, SC>::WaveScatteringProblem( const WaveScatteringProblem& orig ) {
}

template<class LO, class SC>
WaveScatteringProblem<LO, SC>::~WaveScatteringProblem( ) {
  delete [] parameters;
}

template<class LO, class SC>
void WaveScatteringProblem<LO, SC>::set(
    int name,
    SCVT value
    ) {
  switch ( name ) {
    case END_TIME:
      this->endTime = value;
      parameters[END_TIME] = true;
      break;
    case SOLVER_PRECISION:
      this->solverPrecision = value;
      parameters[SOLVER_PRECISION] = true;
      break;
    case SCALE_EVAL:
      this->evalScaleFactor = value;
      parameters[SCALE_EVAL] = true;
      break;
    case SCALE_INPUT:
      this->meshScaleFactor = value;
      parameters[SCALE_INPUT] = true;
      break;
  }
}

template<class LO, class SC>
void WaveScatteringProblem<LO, SC>::set(
    int name,
    const string & value
    ) {
  switch ( name ) {
    case INPUT_MESH_FILE:
      this->inputMeshFile = value;
      parameters[INPUT_MESH_FILE] = true;
      break;
    case EVAL_MESH_FILE:
      this->evalMeshFile = value;
      parameters[EVAL_MESH_FILE] = true;
      break;
    case OUTPUT_FILE:
      this->outputFile = value;
      parameters[OUTPUT_FILE] = true;
      break;
  }
}

template<class LO, class SC>
void WaveScatteringProblem<LO, SC>::set(
    int name,
    SC( *value )( SCVT t, SCVT* x, SCVT *n )
    ) {
  switch ( name ) {
    case INCIDENT_WAVE_DU_DN:
      this->incWaveDuDn = value;
      parameters[INCIDENT_WAVE_DU_DN] = true;
      break;
  }
}

template<class LO, class SC>
void WaveScatteringProblem<LO, SC>::set(
    int name,
    SC( *value )( SCVT t, SCVT* x )
    ) {
  switch ( name ) {
    case INCIDENT_WAVE:
      this->incWave = value;
      parameters[INCIDENT_WAVE] = true;
      break;
  }
}

template<class LO, class SC>
void WaveScatteringProblem<LO, SC>::set(
    int name,
    int value
    ) {
  switch ( name ) {
    case PROBLEM_TYPE:
      this->problemType = (ProblemType) value;
      parameters[PROBLEM_TYPE] = true;
      break;
    case SOLVER_TYPE:
      this->solverType = (SolverType) value;
      parameters[SOLVER_TYPE] = true;
      break;
    case ITERATIVE_SOLVER:
      this->iterativeSolver = (IterativeSolverType) value;
      parameters[ITERATIVE_SOLVER] = true;
      break;
    case N_REFINE_INPUT:
      this->nInputRefines = value;
      parameters[N_REFINE_INPUT] = true;
      break;
    case N_TIME_STEPS:
      this->nTimeSteps = value;
      parameters[N_TIME_STEPS] = true;
      break;
    case LEGENDRE_ORDER:
      this->legendreOrder = value;
      parameters[LEGENDRE_ORDER] = true;
      break;
    case TEMP_QUAD_ORDER:
      this->tempQuadOrder = value;
      parameters[TEMP_QUAD_ORDER] = true;
      break;
    case N_REFINE_EVAL:
      parameters[N_REFINE_EVAL] = true;
      this->nEvalRefines = value;
      break;
    case SOLVER_MAX_IT:
      this->maxIt = value;
      parameters[SOLVER_MAX_IT] = true;
      break;
    case SOLVER_RESTARTS:
      this->gmresRestarts = value;
      parameters[SOLVER_RESTARTS] = true;
      break;
    case DGMRES_MAX_DIM:
      this->dgmresMaxDim = value;
      parameters[DGMRES_MAX_DIM] = true;
      break;
    case DGMRES_MAX_DIM_IT:
      this->dgmresMaxDimIt = value;
      parameters[DGMRES_MAX_DIM_IT] = true;
      break;
  }
}

template<class LO, class SC>
void WaveScatteringProblem<LO, SC>::set(
    int name,
    void * value
    ) {
  switch ( name ) {
    case INPUT_MESH:
      this->mesh = ( SurfaceMesh3D<LO, SC>* ) value;
      parameters[INPUT_MESH] = true;
      break;
    case EVAL_MESH:
      this->evalMesh = ( SurfaceMesh3D<LO, SC>* ) value;
      parameters[EVAL_MESH] = true;
      break;
    case SPACE_QUAD_ORDER:
      memcpy( this->spatialQuad, (int*) value, 4 * sizeof (int ) );
      parameters[SPACE_QUAD_ORDER] = true;
      break;
  }
}

template<class LO, class SC>
void WaveScatteringProblem<LO, SC>::setDefaultParameters( ) {
  solverType = ITERATIVE;
  nInputRefines = 0;
  nEvalRefines = 0;
  legendreOrder = 0;
  memcpy( this->spatialQuad, defaultQuadraturesSauterSchwab, 4 * sizeof (int ) );
  tempQuadOrder = defaultTimeQuadrature;
  solverPrecision = 1e-6;
  maxIt = 10000;
  gmresRestarts = 1000;
  evalScaleFactor = 1.0;
  meshScaleFactor = 1.0;
  iterativeSolver = GMRES;
  dgmresMaxDim = 1;
  dgmresMaxDimIt = 1;
}

template<class LO, class SC>
bool WaveScatteringProblem<LO, SC>::checkParameters( ) const {
  bool ret = true;
  if ( !parameters[INPUT_MESH_FILE] && !parameters[INPUT_MESH] ) {
    ret = false;
    std::cout << "Input mesh not set!" << std::endl;
  }
  if ( !parameters[INCIDENT_WAVE_DU_DN] &&
      parameters[PROBLEM_TYPE] == NEUMANN ) {
    ret = false;
    std::cout << "du/dn of incident wave not set!" << std::endl;
  }
  if ( !parameters[INCIDENT_WAVE] ) {
    ret = false;
    std::cout << "Incident wave not set!" << std::endl;
  }
  if ( !parameters[PROBLEM_TYPE] ) {
    ret = false;
    std::cout << "Problem type not set!" << std::endl;
  }
  if ( !parameters[END_TIME] || !parameters[N_TIME_STEPS] ) {
    ret = false;
    std::cout << "End time or number of time-steps not set!" << std::endl;
  }
  return ret;
}

template<class LO, class SC>
bool WaveScatteringProblem<LO, SC>::solve( ) {
  if ( !checkParameters( ) ) {
    return false;
  }
  bool deleteMesh = false;
  LO localLength;

  // if not provided directly, load mesh from file
  if ( !parameters[INPUT_MESH] ) {
    mesh = new SurfaceMesh3D<LO, SC>;
    mesh->load( inputMeshFile.c_str( ) );
    deleteMesh = true;
  }
  // mesh->printParaviewVtk("sitka.vtk");
  if ( meshScaleFactor != 1.0 ) {
    mesh->scale( meshScaleFactor );
  }
  if ( nInputRefines > 0 ) {
    mesh->refine( nInputRefines );
  }

#ifdef VERBOSE
  mesh->printInfo( );
#endif

  if ( problemType == DIRICHLET ) {
    localLength = mesh->getNElements( );
  } else if ( problemType == NEUMANN ) {
    localLength = mesh->getNNodes( );
  } else {
    localLength = 0;
  }

  MPIBlockMatrix<LO, SC> sysMatrix;
  Vector<LO, SC> rhs;
  Vector<LO, SC> solution( ( legendreOrder + 1 ) *( nTimeSteps )
      * localLength );

  switch ( problemType ) {
    case DIRICHLET:
      prepareDirichletSystem( sysMatrix, rhs );
      break;
    case NEUMANN:
      prepareNeumannSystem( sysMatrix, rhs );
      break;
    case MIXED:
      std::cout << "Not yet implemented!" << std::endl;
      break;
  }

  // based on chosen solver type solve the system on the boundary
  switch ( solverType ) {
    case ITERATIVE:
      solveIteratively( sysMatrix, rhs, solution );
      break;
    case DIRECT:
      solveDirectly( sysMatrix, rhs, solution );
      break;
  }

  if ( parameters[OUTPUT_FILE] ) {
    saveVtu( solution );
  }

  if ( deleteMesh ) {
    delete mesh;
  }
  return true;
}

template<class LO, class SC>
bool WaveScatteringProblem<LO, SC>::solveDirectly(
    MPIBlockMatrix<LO, SC> &matrix,
    Vector<LO, SC> &rhs,
    Vector<LO, SC> &solution
    ) {

  return true;

}

template<class LO, class SC>
bool WaveScatteringProblem<LO, SC>::solveIteratively(
    MPIBlockMatrix<LO, SC> &matrix,
    Vector<LO, SC> &rhs,
    Vector<LO, SC> &solution
    ) {

  LeftIdentityPreconditioner<LO, SC>* M = new LeftIdentityPreconditioner<LO, SC>;
  //WavePreconditioner<LO, SC>* M = new WavePreconditioner<LO, SC>( &matrix, 0 );
  SCVT precision = solverPrecision;
  LO iterations = maxIt;
  bool converged = false;

  switch ( iterativeSolver ) {
    case GMRES:
      converged = matrix.GMRESSolve( rhs, solution, precision, iterations,
          gmresRestarts, M );
      break;
    case DGMRES:
      converged = matrix.DGMRESSolve( rhs, solution, precision, iterations,
          gmresRestarts, dgmresMaxDim, dgmresMaxDimIt, M );
      break;
  }



#ifdef VERBOSE
  std::cout << "Solved in " << iterations << " iterations. Relative erorr: " <<
      precision << std::endl;
#endif

  return converged;
}

template<class LO, class SC>
void WaveScatteringProblem<LO, SC>::prepareDirichletSystem(
    MPIBlockMatrix<LO, SC> &A,
    Vector<LO, SC> &rhs
    ) {

}

template<class LO, class SC>
void WaveScatteringProblem<LO, SC>::prepareNeumannSystem(
    MPIBlockMatrix<LO, SC> &A,
    Vector<LO, SC> &rhs
    ) {
  BESpaceTime< LO, SC > spaceTime( mesh, p1, p1, legendreOrder,
      endTime, nTimeSteps );
  BEBilinearFormWaveHypersingular<LO, SC> formD( &spaceTime,
      spatialQuad, tempQuadOrder );

  // assemble matrix
  formD.assemble( A );

  // assemble RHS
  BEIntegratorWave<LO, SC> integrator( &spaceTime, spatialQuad, tempQuadOrder );

  integrator.getNeumannRHS( rhs, incWaveDuDn );
}

template<class LO, class SC>
void WaveScatteringProblem<LO, SC>::saveVtu(
    Vector<LO, SC> &solution
    ) {

  bool deleteMesh = false;
  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  // if ( rank == 0 ) {
  // if not provided directly, load mesh from file
  if ( !parameters[EVAL_MESH] ) {
    evalMesh = new SurfaceMesh3D<LO, SC>;
    evalMesh->load( evalMeshFile.c_str( ) );
    deleteMesh = true;
  }
  if ( nEvalRefines > 0 ) {
    evalMesh->refine( nEvalRefines );
  }
  if ( evalScaleFactor != 1.0 ) {
    evalMesh->scale( evalScaleFactor );
  }

  LO nPoints = evalMesh->getNNodes( );
  SC *nodes = &( ( *( evalMesh->getNodes( ) ) )[0] );

  // compute on the grid
  std::stringstream file;
  file.fill( '0' );
  file.width( 4 );

  BESpaceTime< LO, SC > spaceTime( mesh, p1, p1, legendreOrder,
      endTime, nTimeSteps );

  SCVT dt = endTime / ( nTimeSteps - 1 );

  int myStart = ( nTimeSteps / size ) * rank;
  int myEnd;
  if ( nTimeSteps % size > rank ) {
    myStart += rank;
    myEnd = myStart + ( nTimeSteps / size ) + 1;
  } else {
    myStart += nTimeSteps % size;
    myEnd = myStart + ( nTimeSteps / size );
  }

  for ( int i = myStart; i < myEnd; i++ ) {
    file.str( std::string( ) );
    file.clear( );
    file << outputFile << "_" << i << ".vtu";

    Vector< LO, SC > scatter( nPoints );
    Vector< LO, SC > incident( nPoints );

    SCVT t = i * dt;

    for ( LO j = 0; j < nPoints; j++ ) {
      SC *x = ( nodes + j * 3 );
      incident.set( j, incWave( t, x ) );
    }

#ifdef VERBOSE
    std::cout << "Evaluating double layer potential (timestep " << i
        << " of " << nTimeSteps << "). " << std::endl;
#endif


    PotentialsWave<LO, SC> pw( &spaceTime, &solution );
    pw.doubleLayerPotential( nodes, nPoints, scatter, i * dt );


    Vector< LO, SC > total( incident );
    total.add( scatter );
    string nodeNamesSol[] = { "scatter", "incident", "total" };
    Vector< LO, SC >* nodalDataSol[] = { &scatter, &incident, &total };
    evalMesh->printParaviewVtu( file.str( ).c_str( ), 3, nodeNamesSol,
        nodalDataSol, 0, nullptr, nullptr );
  }
  //}

  if ( deleteMesh ) {
    delete evalMesh;
  }
}

}
#endif
