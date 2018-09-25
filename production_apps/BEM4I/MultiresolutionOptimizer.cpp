/*!
 * @file    MultiresolutionOptimizer.cpp
 * @author  Michal Merta
 * @author  Jan Zapletal
 * @date    May 19, 2014
 *
 */

#ifdef MULTIRESOLUTIONOPTIMIZER_H

#define BVP_NOT_SOLVED_J (1000000000.0)
#define BVP_NOT_SOLVED_G (100.0)

namespace bem4i {

template<class LO, class SC>
MultiresolutionOptimizer<LO, SC>::MultiresolutionOptimizer(
MOType type
) {

  this->currOptimLevel = 0;
  this->analysisLevel = 0;
  this->maxOptimLevel = 0;
  this->outputPath = "./";
  this->fixedPart = nullptr;
  this->freePartAnal = nullptr;
  this->freePartOptim = nullptr;
  this->freePartOptimRef = nullptr;
  this->subdivMatrix = nullptr;
  this->bvp = nullptr;
  this->optimizer = nullptr;
  this->iter = 0;
  this->iterLevel = 0;
  this->cost = 0.0;
  this->costOpt = std::numeric_limits< double >::max( );
  this->costFirst = 0.0;
  this->gradNormFirst = 0.0;
  this->normalizeCost = false;
  this->normalizeGrad = false;
  this->normalizeOnEachLevel = false;
  this->jump = 0;
  this->nloptJumpTol = 1;
  this->motype = type;
  this->nloptFRel = 1e-2;

#ifdef IPOPT
  this->iAmShallowCopy = false;
  this->ipopt_tol = 1e-4;
  this->ipopt_acceptable_tol = 1e-2;
  this->ipopt_acceptable_iter = 2;
  this->ipopt_max_iter = 30;
  this->ipopt_acceptable_obj_change_tol = 1e-2;
#endif
}

//template<class LO, class SC>
//MultiresolutionOptimizer<LO, SC>::MultiresolutionOptimizer(
//    const SurfaceMesh3D< LO, SC > * fixedPart,
//    const SurfaceMesh3D< LO, SC > * freePart,
//    int analysisLevel,
//    int maxOptimLevel,
//    OptimizationSubproblem< LO, SC > & bvp,
//    const string & outputPath,
//    const string & tagFile,
//    SCVT nloptFRel,
//    int jumpTol
//    ) {
//
//  this->currOptimLevel = 0;
//  this->analysisLevel = analysisLevel;
//  this->maxOptimLevel = maxOptimLevel;
//  this->fixedPart = nullptr;
//  if ( fixedPart ) {
//    this->fixedPart = new SurfaceMesh3D< LO, SC >( *fixedPart );
//  }
//  this->freePartOptim = new OpenMeshWrapper< LO, SC >( *freePart );
//  if ( !tagFile.empty( ) ) {
//    this->freePartOptim->readTags( tagFile );
//  }
//  this->freePartOptimRef = nullptr;
//
//  this->subdivMatrix = nullptr;
//
//  this->setFreePartAnalAndInitMatrices( );
//
//  this->bvp = &bvp;
//  this->iter = 0;
//  this->iterLevel = 0;
//  this->cost = 0.0;
//  //  this->costPrev = 0;
//  this->costOpt = std::numeric_limits< double >::max( );
//  this->costFirst = 0.0;
//  this->gradNormFirst = 0.0;
//  this->normalizeCost = false;
//  this->normalizeGrad = false;
//  this->normalizeOnEachLevel = false;
//  this->jump = 0;
//  this->motype = freeForm;
//
//  const std::string createDir = "mkdir -p " + outputPath;
//  if ( system( createDir.c_str( ) ) ) exit( -1 );
//  this->outputPath = outputPath + "/";
//
//  this->jumpTol = jumpTol;
//  this->nloptFRel = nloptFRel;
//
//  SCVT node[ 3 ];
//  for ( LO i = 0; i < this->freePartOptim->getNNodes( ); ++i ) {
//    this->freePartOptim->getNode( i, node );
//    this->xOpt.push_back( node[ 0 ] );
//    this->xOpt.push_back( node[ 1 ] );
//    this->xOpt.push_back( node[ 2 ] );
//  }
//
//  this->resetBounds( ); //
//
//}

template<class LO, class SC>
MultiresolutionOptimizer<LO, SC>::~MultiresolutionOptimizer( ) {
#ifdef IPOPT
  if ( !this->iAmShallowCopy ) {
#endif
    if ( this->fixedPart ) delete this->fixedPart;
    if ( this->freePartOptim ) delete this->freePartOptim;
    if ( this->freePartOptimRef ) delete this->freePartOptimRef;
    if ( this->freePartAnal ) delete this->freePartAnal;
    if ( this->subdivMatrix ) delete this->subdivMatrix;
    for ( int i = 0; i < this->subdivMatrices.size( ); ++i ) {
      if ( this->subdivMatrices[ i ] ) delete this->subdivMatrices[ i ];
    }
#ifdef IPOPT
  }
#endif
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::coarseToFinePropagate(
    const Vector<LO, SC> & coarse,
    Vector<LO, SC> & fine
    ) const {

  if ( this->currOptimLevel == this->analysisLevel ) {
    coarse.copy( fine );
  } else {
    Vector< LO, SC > tmp( coarse );
    for ( int i = this->currOptimLevel; i < this->analysisLevel; ++i ) {
      fine.resize( this->subdivMatrices[ i ]->getNRows( ) );
      this->subdivMatrices[ i ]->apply( tmp, fine );
      if ( i < this->analysisLevel - 1 ) fine.copy( tmp );
    }
  }
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::fineToCoarsePropagate(
    const Vector<LO, SC> & fine,
    Vector<LO, SC> & coarse
    ) {

  if ( this->currOptimLevel == this->analysisLevel ) {
    fine.copy( coarse );
  } else {
    if ( !this->subdivMatrix ||
        this->subdivMatrix->getNRows( ) != this->freePartAnal->getNNodes( ) ||
        this->subdivMatrix->getNCols( ) != this->freePartOptim->getNNodes( ) ) {
      this->computeSubdivMatrix( );
    }

    coarse.resize( this->freePartOptim->getNNodes( ) );
    this->subdivMatrix->QRSolve( fine, coarse );
  }
}

template<class LO, class SC>
bool MultiresolutionOptimizer<LO, SC>::readRelativeBounds(
    const std::string& boundsFile
    ) {

  std::cout << "Reading file '" << boundsFile << "' ... ";
  std::ifstream file( boundsFile.c_str( ) );

  if ( !file.is_open( ) ) {
    std::cout << "File '" << boundsFile << "' could not be opened!"
        << std::endl;
    return false;
  }

  LO nNodes = this->freePartOptim->getNNodes( );
  LO nConstraints;
  LO nodeIdx;
  SCVT lower, upper;
  LO nConstraintsPerNode = 0;

  switch ( this->motype ) {
    case freeForm:
      nConstraintsPerNode = 3;
      break;
    case fixedRef:
      nConstraintsPerNode = 1;
      break;
    default:
      std::cout << "Wrong optimization type!" << std::endl;
      break;
  }

  file >> nConstraints;

  for ( LO i = 0; i < nConstraints; ++i ) {
    file >> nodeIdx;

    if ( nodeIdx >= nNodes ) {
      this->resetBounds( );
      return false;
    }

    for ( int j = 0; j < nConstraintsPerNode; ++j ) {
      file >> lower;
      file >> upper;
      if ( lower < 0.0 ) {
        lower = std::numeric_limits<double>::infinity( );
      }
      if ( upper < 0.0 ) {
        upper = std::numeric_limits<double>::infinity( );
      }
      this->relativeLowerBounds[ nConstraintsPerNode * nodeIdx + j ] = lower;
      this->relativeUpperBounds[ nConstraintsPerNode * nodeIdx + j ] = upper;
    }
  }

  file.close( );

  std::cout << "done." << std::endl;

  return true;
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::resetBounds( ) {

  LO dim = 0;

  switch ( this->motype ) {
    case freeForm:
      dim = 3 * this->freePartOptim->getNNodes( );
      break;
    case fixedRef:
      dim = this->freePartOptim->getNNodes( );
      break;
    default:
      std::cout << "Wrong optimization type!" << std::endl;
      break;
  }

  this->relativeLowerBounds.clear( );
  this->relativeLowerBounds.resize( dim,
      std::numeric_limits<double>::infinity( ) );
  this->relativeUpperBounds.clear( );
  this->relativeUpperBounds.resize( dim,
      std::numeric_limits<double>::infinity( ) );
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::setRelativeBounds(
    SCVT bound
    ) {

  LO dim = 0;

  switch ( this->motype ) {
    case freeForm:
      dim = 3 * this->freePartOptim->getNNodes( );
      break;
    case fixedRef:
      dim = this->freePartOptim->getNNodes( );
      break;
    default:
      std::cout << "Wrong optimization type!" << std::endl;
      break;
  }

  this->relativeLowerBounds.clear( );
  this->relativeLowerBounds.resize( dim, (double) bound );
  this->relativeUpperBounds.clear( );
  this->relativeUpperBounds.resize( dim, (double) bound );

}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::setUpAbsoluteBounds(
    std::vector< double > & lb,
    std::vector< double > & ub
    ) const {

  LO nNodes = this->freePartOptim->getNNodes( );
  lb.clear( );
  ub.clear( );

  switch ( this->motype ) {
    case freeForm:
      lb.reserve( 3 * nNodes );
      ub.reserve( 3 * nNodes );
      SCVT node[ 3 ];
      for ( LO i = 0; i < nNodes; ++i ) {
        this->freePartOptim->getNode( i, node );
        lb.push_back( node[ 0 ] - this->relativeLowerBounds[ 3 * i ] );
        lb.push_back( node[ 1 ] - this->relativeLowerBounds[ 3 * i + 1 ] );
        lb.push_back( node[ 2 ] - this->relativeLowerBounds[ 3 * i + 2 ] );
        ub.push_back( node[ 0 ] + this->relativeUpperBounds[ 3 * i ] );
        ub.push_back( node[ 1 ] + this->relativeUpperBounds[ 3 * i + 1 ] );
        ub.push_back( node[ 2 ] + this->relativeUpperBounds[ 3 * i + 2 ] );
      }
      break;

    case fixedRef:
      lb.reserve( nNodes );
      ub.reserve( nNodes );
      for ( LO i = 0; i < nNodes; ++i ) {
        lb.push_back( this->alphaOpt[ i ] -
            this->relativeLowerBounds[ i ] );
        ub.push_back( this->alphaOpt[ i ] +
            this->relativeUpperBounds[ i ] );
      }
      break;

    default:
      std::cout << "Wrong optimization type!" << std::endl;
      break;
  }


}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::printInfo( ) const {
  std::cout << "Class MultiresolutionOptimizer, " << std::endl;
  std::cout << "\tcurrent optimization level: "
      << this->currOptimLevel << "," << std::endl;
  std::cout << "\tmax optimization level: "
      << this->maxOptimLevel << "," << std::endl;
  std::cout << "\tanalysis level: "
      << this->analysisLevel << "," << std::endl;
  if ( fixedPart ) {
    std::cout << "\tfixed mesh: ";
    this->fixedPart->printInfo( );
  }
  std::cout << "\tfree mesh on current optimization level: ";
  this->freePartOptim->printInfo( );
  std::cout << "\tfree mesh on analysis level: ";
  this->freePartAnal->printInfo( );
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::printFixedPartVtu(
    const string & filename
    ) const {
  this->fixedPart->printParaviewVtu( filename );
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::printFreePartAnalVtu(
    const string & filename
    ) const {
  this->freePartAnal->printParaviewVtu( filename );
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::printFreePartAnalVtu(
    const string & meshFile,
    const std::vector< string > * nodeNames,
    const std::vector< Vector< LO, SC >* > * nodalData,
    const std::vector< string > * elemNames,
    const std::vector< Vector< LO, SC >* > * elemData
    ) const {
  this->freePartAnal->printParaviewVtu( meshFile, nodeNames, nodalData,
      elemNames, elemData );
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::printFreePartOptimVtu(
    const string & meshFile,
    const std::vector< string > * nodeNames,
    const std::vector< Vector< LO, SC >* > * nodalData,
    const std::vector< string > * elemNames,
    const std::vector< Vector< LO, SC >* > * elemData
    ) const {
  this->freePartOptim->printParaviewVtu( meshFile, nodeNames, nodalData,
      elemNames, elemData );
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::printFreePartOptimVtu(
    const string & filename
    ) const {
  this->freePartOptim->printParaviewVtu( filename );
}

//---NLOPT--------------------------------------------------------------------//

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::optimize( ) {
  switch ( this->motype ) {
    case freeForm:
      this->optimizeFreeForm( );
      break;
    case fixedRef:
#ifdef IPOPT
      this->optimizeFixedReferenceIpopt( );
#else
      this->optimizeFixedReference( );
#endif
      break;
    default:
      std::cout << "Wrong optimization type!" << std::endl;
      break;
  }
}

template<class LO, class SC>
double MultiresolutionOptimizer<LO, SC>::nloptCostFunctionFreeForm(
    const std::vector< double > & x,
    std::vector< double > & grad,
    void * data
    ) {

  // cast to obtain 'this' pointer
  MultiresolutionOptimizer< LO, SC > * myself =
      static_cast<MultiresolutionOptimizer< LO, SC > *> ( data );

  // updating free mesh on optim level
  myself->freePartOptim->setNodes( x.data( ) );

  // subdivide to obtain current analysis mesh
  myself->updateFreePartAnal( );

  // set current mesh for bvp
  myself->bvp->setProblemData( myself->freePartAnal, myself->fixedPart );

  // solve the bvp
  if ( !myself->bvp->solve( ) ) {
    //    std::cout << "Force stop: bvp not solved!" << std::endl;
    //    myself->optimizer->force_stop( );
    std::cout << "BVP not solved!" << std::endl;
    std::fill( grad.begin( ), grad.end( ), BVP_NOT_SOLVED_G );
    return BVP_NOT_SOLVED_J;
  }

  bool renormalized;
  if ( myself->normalizeOnEachLevel ) {
    renormalized = ( myself->iterLevel == 0 );
  } else {
    renormalized = ( myself->iter == 0 );
  }

  // get gradient
  Vector< LO, SC > gradV;
  //myself->getGradientFree( gradV );
  myself->getGradientNormal( gradV );
  //myself->getGradientKernel( gradV );
  if ( renormalized ) {
    myself->gradNormFirst = gradV.norm2( );
  }
  if ( myself->normalizeGrad ) {
    gradV.scale( 1.0 / myself->gradNormFirst );
  }
  if ( !grad.empty( ) ) {
    for ( LO i = 0; i < gradV.getLength( ); ++i ) {
      grad[ i ] = gradV.get( i );
    }
  }

  // get cost
  SCVT cost;
  myself->bvp->getCost( cost );
  if ( renormalized ) {
    myself->costFirst = cost;
  }
  if ( myself->normalizeCost ) {
    cost /= myself->costFirst;
  }

  //  myself->costPrev = myself->cost;
  myself->cost = cost;
  if ( cost < myself->costOpt || renormalized ) {
    myself->xOpt = x;
  }

  // print info and files
  myself->bvp->printInfo( );
  std::cout << "\tJ = " << cost << "." << std::endl;
  std::cout << "\t|g|_2 = " << gradV.norm2( ) << "." << std::endl;

  // print analysis level mesh (only print not worse result)
  // when changing levels = also prints same result on finer optimization mesh
  if ( cost <= myself->costOpt || renormalized ) {
    std::stringstream file;
    std::stringstream iter2string;
    iter2string.width( 4 );
    iter2string.fill( '0' );
    iter2string << myself->iter;
    file << myself->outputPath << "optim_a_" << iter2string.str( ) << ".vtu";
    myself->bvp->printVtu( file.str( ) );

    // print optimization level mesh (only free part) with gradient data
    std::vector< Vector< LO, SC > * > nodeVData;
    gradV.scale( -1.0 );
    nodeVData.push_back( &gradV );
    std::vector< string > nodeVNames;
    nodeVNames.push_back( "neg_shape_gradient" );
    file.str( std::string( ) );
    file.clear( );
    file << myself->outputPath << "optim_o_" << iter2string.str( ) << ".vtu";
    myself->freePartOptim->printParaviewVtu( file.str( ), nullptr, nullptr,
        nullptr, nullptr, &nodeVNames, &nodeVData );

    // print analysis level of free mesh
    file.str( std::string( ) );
    file.clear( );
    file << myself->outputPath << "optim_s_" << iter2string.str( ) << ".vtu";
    myself->freePartAnal->printParaviewVtu( file.str( ) );
  }

  if ( myself->iterLevel > 0 )
    myself->convergenceMonitor( );

  if ( cost < myself->costOpt || renormalized ) {
    myself->costOpt = cost;
  }

  ++myself->iter;
  ++myself->iterLevel;

  return cost;
}

template<class LO, class SC>
double MultiresolutionOptimizer<LO, SC>::nloptCostFunctionFixedReference(
    const std::vector< double > & alpha,
    std::vector< double > & grad,
    void * data
    ) {

  // cast to obtain 'this' pointer
  MultiresolutionOptimizer< LO, SC > * myself =
      static_cast<MultiresolutionOptimizer< LO, SC > *> ( data );

  // updating free mesh on optim level
  myself->updateFreePartOptimFromReference( alpha );

  // subdivide to obtain current analysis mesh
  myself->updateFreePartAnal( );

  // set current mesh for bvp
  myself->bvp->setProblemData( myself->freePartAnal, myself->fixedPart );

  // solve the bvp
  if ( !myself->bvp->solve( ) ) {
    //std::cout << "Force stop: bvp not solved!" << std::endl;
    //myself->optimizer->force_stop( );
    std::cout << "BVP not solved!" << std::endl;
    std::fill( grad.begin( ), grad.end( ), BVP_NOT_SOLVED_G );
    return BVP_NOT_SOLVED_J;
  }

  bool renormalized;
  if ( myself->normalizeOnEachLevel ) {
    renormalized = ( myself->iterLevel == 0 );
  } else {
    renormalized = ( myself->iter == 0 );
  }

  // get gradient
  Vector< LO, SC > gradV;
  myself->getGradientNormalFixedReference( gradV );
  if ( renormalized ) {
    myself->gradNormFirst = gradV.norm2( );
  }
  if ( myself->normalizeGrad ) {
    gradV.scale( 1.0 / myself->gradNormFirst );
  }
  if ( !grad.empty( ) ) {
    for ( LO i = 0; i < gradV.getLength( ); ++i ) {
      grad[ i ] = gradV.get( i );
    }
  }

  // get cost
  SCVT cost;
  myself->bvp->getCost( cost );
  if ( renormalized ) {
    myself->costFirst = cost;
  }
  if ( myself->normalizeCost ) {
    cost /= myself->costFirst;
  }

  //  myself->costPrev = myself->cost;
  myself->cost = cost;
  if ( cost < myself->costOpt || renormalized ) {
    myself->freePartOptim->getNodes( myself->xOpt );
    myself->alphaOpt = alpha;
  }

  // print info and files
  myself->bvp->printInfo( );
  std::cout << "\tJ = " << cost << "." << std::endl;
  std::cout << "\t|g|_2 = " << gradV.norm2( ) << "." << std::endl;

  // print analysis level mesh (only print not worse result)
  // when changing levels = also prints same result on finer optimization mesh
  if ( cost <= myself->costOpt || renormalized ) {
    std::stringstream file;
    std::stringstream iter2string;
    iter2string.width( 4 );
    iter2string.fill( '0' );
    iter2string << myself->iter;
    file << myself->outputPath << "optim_a_" << iter2string.str( ) << ".vtu";
    myself->bvp->printVtu( file.str( ) );

    // print optimization level mesh (only free part) with gradient data
    LO dim = gradV.getLength( );
    std::vector< Vector< LO, SC > * > nodeVData;
    gradV.scale( -1.0 );
    Vector< LO, SC > gradVN( 3 * dim );
    SCVT n[ 3 ];
    for ( LO i = 0; i < dim; ++i ) {
      myself->freePartOptim->getNodalNormal( i, n );
      gradVN.set( 3 * i, gradV.get( i ) * n[ 0 ] );
      gradVN.set( 3 * i + 1, gradV.get( i ) * n[ 1 ] );
      gradVN.set( 3 * i + 2, gradV.get( i ) * n[ 2 ] );
    }
    nodeVData.push_back( &gradVN );
    std::vector< string > nodeVNames;
    nodeVNames.push_back( "neg_shape_gradient" );
    file.str( std::string( ) );
    file.clear( );
    file << myself->outputPath << "optim_o_" << iter2string.str( ) << ".vtu";
    myself->freePartOptim->printParaviewVtu( file.str( ), nullptr, nullptr,
        nullptr, nullptr, &nodeVNames, &nodeVData );

    // print analysis level of free mesh
    file.str( std::string( ) );
    file.clear( );
    file << myself->outputPath << "optim_s_" << iter2string.str( ) << ".vtu";
    myself->freePartAnal->printParaviewVtu( file.str( ) );
  }

  if ( myself->iterLevel > 0 )
    myself->convergenceMonitor( );

  if ( cost < myself->costOpt || renormalized ) {
    myself->costOpt = cost;
  }

  ++myself->iter;
  ++myself->iterLevel;

  return cost;
}

//template<class LO, class SC>
//void MultiresolutionOptimizer<LO, SC>::getBoundsOnCoarseMesh( ) {
//
//  LO nNodes = this->freePartOptim->getNNodes( );
//  LO size = nNodes * 3;
//  SCVT * x = new SCVT[ size ];
//  this->freePartOptim->getNodes( x );
//
//  this->relativeLowerBounds.clear( );
//  this->relativeUpperBounds.clear( );
//  this->relativeLowerBounds.reserve( size );
//  this->relativeUpperBounds.reserve( size );
//
//  SCVT length;
//  SCVT mult = 0.6;
//  SCVT shift;
//  for ( LO i = 0; i < nNodes; ++i ) {
//    length = this->freePartOptim->getVertexAverageEdgeLength( i );
//    shift = mult * length;
//    this->relativeLowerBounds.push_back( x[ 3 * i ] - shift );
//    this->relativeLowerBounds.push_back( x[ 3 * i + 1 ] - shift );
//    this->relativeLowerBounds.push_back( x[ 3 * i + 2 ] - shift );
//    this->relativeUpperBounds.push_back( x[ 3 * i ] + shift );
//    this->relativeUpperBounds.push_back( x[ 3 * i + 1 ] + shift );
//    this->relativeUpperBounds.push_back( x[ 3 * i + 2 ] + shift );
//  }
//}

//---PRIVATE-METHODS----------------------------------------------------------//

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::getGradientFree(
    Vector<LO, SC> & grad
    ) const {

  LO nNodes = this->freePartOptim->getNNodes( );
  grad.resize( nNodes * 3 );
  Vector< LO, SCVT > pertCoarse( nNodes );
  Vector< LO, SCVT > pertFine;
  SCVT dx1, dx2, dx3;

  for ( LO i = 0; i < nNodes; ++i ) {
    pertCoarse.setAll( 0.0 );
    pertCoarse.set( i, 1.0 );
    this->coarseToFinePropagate( pertCoarse, pertFine );
    this->bvp->getShapeGradient( pertFine, dx1, dx2, dx3 );
    grad.set( 3 * i, dx1 );
    grad.set( 3 * i + 1, dx2 );
    grad.set( 3 * i + 2, dx3 );
  }
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::getGradientNormal(
    Vector<LO, SC> & grad
    ) const {

  LO nNodes = this->freePartOptim->getNNodes( );
  grad.resize( nNodes * 3 );
  Vector< LO, SCVT > pertCoarseX1( nNodes );
  Vector< LO, SCVT > pertCoarseX2( nNodes );
  Vector< LO, SCVT > pertCoarseX3( nNodes );
  Vector< LO, SCVT > pertFineX1;
  Vector< LO, SCVT > pertFineX2;
  Vector< LO, SCVT > pertFineX3;
  SCVT dx;
  SCVT n[ 3 ];

  for ( LO i = 0; i < nNodes; ++i ) {

    this->freePartOptim->getNodalNormal( i, n );

    pertCoarseX1.setAll( 0.0 );
    pertCoarseX1.set( i, n[ 0 ] );
    pertCoarseX2.setAll( 0.0 );
    pertCoarseX2.set( i, n[ 1 ] );
    pertCoarseX3.setAll( 0.0 );
    pertCoarseX3.set( i, n[ 2 ] );

    this->coarseToFinePropagate( pertCoarseX1, pertFineX1 );
    this->coarseToFinePropagate( pertCoarseX2, pertFineX2 );
    this->coarseToFinePropagate( pertCoarseX3, pertFineX3 );
    this->bvp->getShapeGradient( pertFineX1, pertFineX2, pertFineX3, dx );
    grad.set( 3 * i, dx * n[ 0 ] );
    grad.set( 3 * i + 1, dx * n[ 1 ] );
    grad.set( 3 * i + 2, dx * n[ 2 ] );
  }

}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::getGradientNormalFixedReference(
    Vector<LO, SC> & grad
    ) const {

  LO nNodes = this->freePartOptim->getNNodes( );
  grad.resize( nNodes );
  Vector< LO, SCVT > pertCoarseX1( nNodes );
  Vector< LO, SCVT > pertCoarseX2( nNodes );
  Vector< LO, SCVT > pertCoarseX3( nNodes );
  Vector< LO, SCVT > pertFineX1;
  Vector< LO, SCVT > pertFineX2;
  Vector< LO, SCVT > pertFineX3;
  SCVT dx;
  SCVT n[ 3 ];

  for ( LO i = 0; i < nNodes; ++i ) {
// todo: should be freePartOptimRef??? (the perturbation goes in this direction)
    //this->freePartOptim->getNodalNormal( i, n );
    this->freePartOptimRef->getNodalNormal( i, n );

    pertCoarseX1.setAll( 0.0 );
    pertCoarseX1.set( i, n[ 0 ] );
    pertCoarseX2.setAll( 0.0 );
    pertCoarseX2.set( i, n[ 1 ] );
    pertCoarseX3.setAll( 0.0 );
    pertCoarseX3.set( i, n[ 2 ] );

    this->coarseToFinePropagate( pertCoarseX1, pertFineX1 );
    this->coarseToFinePropagate( pertCoarseX2, pertFineX2 );
    this->coarseToFinePropagate( pertCoarseX3, pertFineX3 );
    this->bvp->getShapeGradient( pertFineX1, pertFineX2, pertFineX3, dx );
    grad.set( i, dx );
  }

}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::getGradientKernel(
    Vector<LO, SC>& grad
    ) {

  LO nNodes = this->freePartOptim->getNNodes( );
  grad.resize( nNodes * 3 );
  SCVT n[ 3 ];

  Vector< LO, SC > kernelFine;
  this->bvp->getShapeGradient( kernelFine );
  Vector< LO, SC > kernelCoarse;
  this->fineToCoarsePropagate( kernelFine, kernelCoarse );
  SC kernel;

  for ( LO i = 0; i < nNodes; ++i ) {
    kernel = kernelCoarse.get( i );
    this->freePartOptim->getNodalNormal( i, n );
    grad.set( 3 * i, kernel * n[ 0 ] );
    grad.set( 3 * i + 1, kernel * n[ 1 ] );
    grad.set( 3 * i + 2, kernel * n[ 2 ] );
  }

}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::optimizeFreeForm( ) {

  typedef bem4i::MultiresolutionOptimizer< LO, SC > myClass;

  LO nNodes;
  LO dim;
  std::vector< double > x, lb, ub;
  double minf;
  nlopt::result result;

  while ( true ) {
    nNodes = this->freePartOptim->getNNodes( );
    dim = nNodes * 3;

    // choosing algorithm, set optimization dimension (3 * nNodes))
    // MMA
    this->optimizer = new nlopt::opt( nlopt::LD_MMA, dim );
    // StoGo
    //    this->optimizer = new nlopt::opt( nlopt::GD_STOGO, dim );
    // MLSL
    //    this->optimizer = new nlopt::opt( nlopt::GD_MLSL_LDS, dim );
    //    nlopt::opt localOptimizer( nlopt::LD_MMA, dim );
    //    localOptimizer.set_ftol_rel( this->nloptFRel );
    //    this->optimizer->set_local_optimizer( localOptimizer );
    // BFGS
    //    this->optimizer = new nlopt::opt( nlopt::LD_LBFGS, dim );

    // set cost function
    this->optimizer->set_min_objective( &myClass::nloptCostFunctionFreeForm,
        this );

    // set bounds on coarse mesh
    this->setUpAbsoluteBounds( lb, ub );
    this->optimizer->set_lower_bounds( lb );
    this->optimizer->set_upper_bounds( ub );

    // set stopping criterion for relative cost decrease
    this->optimizer->set_ftol_rel( this->nloptFRel );

    // set current best as starting point
    x = this->xOpt;

    // optimize
    try {
      result = this->optimizer->optimize( x, minf );
    } catch ( nlopt::forced_stop e ) {
      result = nlopt::FORCED_STOP;
    }

    delete this->optimizer;

    this->iterLevel = 0;

    std::cout << "Optimization return code: " << result << "." << std::endl;

    if ( this->currOptimLevel < this->maxOptimLevel ) {
      // increase optimization level

      // fixme: doesnt have to be current best?
      std::stringstream file;
      file << this->outputPath << "optim_level_" << this->currOptimLevel
          << "_iter_" << this->iter - 1 << ".txt";
      this->freePartOptim->print( file.str( ) );

      this->increaseOptimizationLevel( );

      std::cout << "------Optimization level increased------" << std::endl;
      std::cout << "\tglobal iteration: "
          << this->iter - 1 << "," << std::endl;
      std::cout << "\tcurrent optimization level: "
          << this->currOptimLevel << "," << std::endl;
      std::cout << "\tmax optimization level: "
          << this->maxOptimLevel << "," << std::endl;
      std::cout << "\tanalysis level: "
          << this->analysisLevel << "." << std::endl;
    } else {
      // end of optimization

      // fixme: doesnt have to be current best?
      std::stringstream file;
      file << this->outputPath << "optim_level_" << this->currOptimLevel
          << "_iter_" << this->iter - 1 << ".txt";
      this->freePartOptim->print( file.str( ) );

      break;
    }
  }
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::optimizeFixedReference( ) {

  typedef bem4i::MultiresolutionOptimizer< LO, SC > myClass;

  LO nNodes = this->freePartOptim->getNNodes( );
  LO dim = nNodes;
  std::vector< double > alpha, lb, ub;
  double minf;
  nlopt::result result;

  while ( true ) {
    nNodes = this->freePartOptim->getNNodes( );
    dim = nNodes;

    // choosing algorithm, set optimization dimension (3 * nNodes))
    // MMA
    this->optimizer = new nlopt::opt( nlopt::LD_MMA, dim );
    // StoGo
    //    this->optimizer = new nlopt::opt( nlopt::GD_STOGO, dim );
    // MLSL
    //    this->optimizer = new nlopt::opt( nlopt::GD_MLSL_LDS, dim );
    //    nlopt::opt localOptimizer( nlopt::LD_MMA, dim );
    //    localOptimizer.set_ftol_rel( this->nloptFRel );
    //    this->optimizer->set_local_optimizer( localOptimizer );
    // BFGS
    //this->optimizer = new nlopt::opt( nlopt::LD_LBFGS, dim );

    // set cost function
    this->optimizer->set_min_objective(
        &myClass::nloptCostFunctionFixedReference, this );

    // set bounds on coarse mesh
    this->setUpAbsoluteBounds( lb, ub );
    this->optimizer->set_lower_bounds( lb );
    this->optimizer->set_upper_bounds( ub );

    // set stopping criterion for relative cost decrease
    this->optimizer->set_ftol_rel( this->nloptFRel );

    alpha = this->alphaOpt;

    // optimize
    try {
      result = this->optimizer->optimize( alpha, minf );
    } catch ( nlopt::forced_stop e ) {
      result = nlopt::FORCED_STOP;
    }

    delete this->optimizer;

    this->iterLevel = 0;

    std::cout << "Optimization return code: " << result << "." << std::endl;

    if ( this->currOptimLevel < this->maxOptimLevel ) {
      // increase optimization level

      std::stringstream file;
      file << this->outputPath << "optim_level_" << this->currOptimLevel
          << "_iter_" << this->iter - 1 << ".txt";
      this->freePartOptim->print( file.str( ) );

      this->increaseOptimizationLevel( );

      std::cout << "------Optimization level increased------" << std::endl;
      std::cout << "\tglobal iteration: "
          << this->iter - 1 << "," << std::endl;
      std::cout << "\tcurrent optimization level: "
          << this->currOptimLevel << "," << std::endl;
      std::cout << "\tmax optimization level: "
          << this->maxOptimLevel << "," << std::endl;
      std::cout << "\tanalysis level: "
          << this->analysisLevel << "." << std::endl;
    } else {
      // end of optimization

      std::stringstream file;
      file << this->outputPath << "optim_level_" << this->currOptimLevel
          << "_iter_" << this->iter - 1 << ".txt";
      this->freePartOptim->print( file.str( ) );

      break;
    }
  }
}

#ifdef IPOPT

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::optimizeFixedReferenceIpopt( ) {

  Ipopt::SmartPtr< MultiresolutionOptimizer< LO, SC > > copy =
      this->getShallowCopy( );

  while ( true ) {

    Ipopt::SmartPtr< Ipopt::IpoptApplication > app = IpoptApplicationFactory( );

    app->RethrowNonIpoptException( true );
    app->Options( )->SetNumericValue( "tol", copy->ipopt_tol );
    app->Options( )->SetNumericValue( "acceptable_tol",
        this->ipopt_acceptable_tol );
    app->Options( )->SetIntegerValue( "acceptable_iter",
        copy->ipopt_acceptable_iter );
    app->Options( )->SetIntegerValue( "max_iter",
        copy->ipopt_max_iter );
    app->Options( )->SetNumericValue( "acceptable_obj_change_tol",
        copy->ipopt_acceptable_obj_change_tol );

    app->Options( )->SetStringValue( "hessian_approximation",
        "limited-memory" );
    std::stringstream file;
    file << copy->outputPath << "ipopt_" << copy->currOptimLevel << ".txt";
    app->Options( )->SetStringValue( "output_file", file.str( ) );

    Ipopt::ApplicationReturnStatus status;
    status = app->Initialize( );
    if ( status != Ipopt::Solve_Succeeded ) {
      std::cout << std::endl << std::endl <<
          "*** Error during initialization!" << std::endl;
      std::cout << "Status: " << (int) status << "." << std::endl;
    }

    status = app->OptimizeTNLP( copy );

    std::cout << "Status: " << (int) status << "." << std::endl;

    copy->iterLevel = 0;

    if ( copy->currOptimLevel < copy->maxOptimLevel ) {

      copy->increaseOptimizationLevel( );

      std::cout << "------Optimization level increased------" << std::endl;
      std::cout << "\tglobal iteration: "
          << copy->iter - 1 << "," << std::endl;
      std::cout << "\tcurrent optimization level: "
          << copy->currOptimLevel << "," << std::endl;
      std::cout << "\tmax optimization level: "
          << copy->maxOptimLevel << "," << std::endl;
      std::cout << "\tanalysis level: "
          << copy->analysisLevel << "." << std::endl;
    } else {
      // end of optimization
      break;
    }
  }

  copy->shallowCopy( *this );
}

#endif

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::increaseOptimizationLevel( ) {
  switch ( this->motype ) {
    case freeForm:
      this->increaseOptimizationLevelFreeForm( );
      break;
    case fixedRef:
      this->increaseOptimizationLevelFixedRef( );
      break;
    default:
      std::cout << "Wrong optimization type!" << std::endl;
      break;
  }
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::increaseOptimizationLevelFreeForm( ) {
  if ( this->currOptimLevel >= this->maxOptimLevel )
    return;

  this->setBoundsForNextLevel( );

  this->freePartOptim->setNodes( this->xOpt.data( ) );
  this->freePartOptim->subdivideLoop( 1 );
  this->currOptimLevel++;
  //this->nloptFRel /= 10;

  // setting up subdivided optimum
  SCVT node[ 3 ];
  LO nNodes = this->freePartOptim->getNNodes( );
  this->xOpt.clear( );
  this->xOpt.reserve( 3 * nNodes );
  for ( LO i = 0; i < nNodes; ++i ) {
    this->freePartOptim->getNode( i, node );
    this->xOpt.push_back( node[ 0 ] );
    this->xOpt.push_back( node[ 1 ] );
    this->xOpt.push_back( node[ 2 ] );
  }

}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::increaseOptimizationLevelFixedRef( ) {
  if ( this->currOptimLevel >= this->maxOptimLevel )
    return;

  this->setBoundsForNextLevel( );

  this->freePartOptim->setNodes( this->xOpt.data( ) );
  this->freePartOptim->subdivideLoop( 1 );
  this->currOptimLevel++;
  //this->nloptFRel /= 10;

  // setting up subdivided optimum
  SCVT node[ 3 ];
  LO nNodes = this->freePartOptim->getNNodes( );
  this->xOpt.clear( );
  this->xOpt.reserve( 3 * nNodes );
  for ( LO i = 0; i < nNodes; ++i ) {
    this->freePartOptim->getNode( i, node );
    this->xOpt.push_back( node[ 0 ] );
    this->xOpt.push_back( node[ 1 ] );
    this->xOpt.push_back( node[ 2 ] );
  }

  this->alphaOpt.clear( );
  this->alphaOpt.resize( nNodes, 0.0 );

  if ( this->freePartOptimRef ) delete this->freePartOptimRef;
  this->freePartOptimRef =
      new OpenMeshWrapper< LO, SC >( *this->freePartOptim );
  this->freePartOptimRef->setNodes( this->xOpt.data( ) );

}

//template<class LO, class SC>
//void MultiresolutionOptimizer<LO, SC>::increaseOptimizationLevelFixedRef( ) {
//  if ( this->currOptimLevel >= this->maxOptimLevel )
//    return;
//
//  // this->setBoundsForNextLevelFixedRef( );
//
//  // current optimum
//  this->freePartOptim->setNodes( this->xOpt.data( ) );
//
//  // setting ups rhs for transferring alpha
//  LO nNodes = this->freePartOptim->getNNodes( );
//  LO nNodesNew = nNodes + this->freePartOptim->getNEdges( );
//  Vector< LO, SC > an1( nNodes );
//  Vector< LO, SC > an2( nNodes );
//  Vector< LO, SC > an3( nNodes );
//  Vector< LO, SC > san1( nNodesNew );
//  Vector< LO, SC > san2( nNodesNew );
//  Vector< LO, SC > san3( nNodesNew );
//  Vector< LO, SC > b( 3 * nNodesNew );
//  Vector< LO, SC > alpha( nNodesNew );
//
//  SCVT n[ 3 ];
//  for ( LO i = 0; i < nNodes; ++i ) {
//    this->freePartOptim->getNodalNormal( i, n );
//    an1.set( i, this->alphaOpt[ i ] * n[ 0 ] );
//    an2.set( i, this->alphaOpt[ i ] * n[ 1 ] );
//    an3.set( i, this->alphaOpt[ i ] * n[ 2 ] );
//  }
//
//  this->subdivMatrices[ this->currOptimLevel ]->apply( an1, san1 );
//  this->subdivMatrices[ this->currOptimLevel ]->apply( an2, san2 );
//  this->subdivMatrices[ this->currOptimLevel ]->apply( an3, san3 );
//
//  for ( LO i = 0; i < nNodesNew; ++i ) {
//    b.set( i, san1.get( i ) );
//    b.set( i + nNodesNew, san2.get( i ) );
//    b.set( i + 2 * nNodesNew, san3.get( i ) );
//  }
//
//  this->freePartOptim->subdivideLoop( 1 );
//  this->currOptimLevel++;
//  //this->nloptFRel /= 10;
//
//  // setting up system matrix for transferring alpha
//  std::vector< LO > rowInd;
//  std::vector< LO > colInd;
//  std::vector< SCVT > val;
//  rowInd.reserve( 3 * nNodesNew );
//  colInd.reserve( 3 * nNodesNew );
//  val.reserve( 3 * nNodesNew );
//  for ( LO i = 0; i < nNodesNew; ++i ) {
//    this->freePartOptim->getNodalNormal( i, n );
//    rowInd.push_back( i );
//    rowInd.push_back( i + nNodesNew );
//    rowInd.push_back( i + 2 * nNodesNew );
//    colInd.push_back( i );
//    colInd.push_back( i );
//    colInd.push_back( i );
//    val.push_back( n[ 0 ] );
//    val.push_back( n[ 1 ] );
//    val.push_back( n[ 2 ] );
//  }
//  SparseMatrix< LO, SCVT > A( 3 * nNodesNew, nNodesNew, rowInd, colInd, val );
//  A.QRSolve( b, alpha );
//
//  // copy result to this->alphaOpt
//  this->alphaOpt.clear();
//  this->alphaOpt.reserve( nNodesNew );
//  for ( LO i = 0; i < nNodesNew; ++i ) {
//    this->alphaOpt.push_back( alpha.get( i ) );
//  }
//
//  // setting up subdivided optimum
//  SCVT node[ 3 ];
//  this->xOpt.clear( );
//  this->xOpt.reserve( 3 * nNodes );
//  for ( LO i = 0; i < nNodes; ++i ) {
//    this->freePartOptim->getNode( i, node );
//    this->xOpt.push_back( node[ 0 ] );
//    this->xOpt.push_back( node[ 1 ] );
//    this->xOpt.push_back( node[ 2 ] );
//  }
//
//}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::setBoundsForNextLevel( ) {

  SCVT infty = std::numeric_limits<double>::infinity( );

  LO nConstraintsPerNode = 0;

  switch ( this->motype ) {
    case freeForm:
      nConstraintsPerNode = 3;
      break;
    case fixedRef:
      nConstraintsPerNode = 1;
      break;
    default:
      std::cout << "Wrong optimization type!" << std::endl;
      break;
  }

  LO nEdges = this->freePartOptim->getNEdges( );
  this->relativeLowerBounds.reserve( this->relativeLowerBounds.size( ) +
      nConstraintsPerNode * nEdges );
  this->relativeUpperBounds.reserve( this->relativeUpperBounds.size( ) +
      nConstraintsPerNode * nEdges );

  LO start;
  LO end;
  SCVT lNew;
  SCVT uNew;
  for ( LO i = 0; i < nEdges; ++i ) {
    this->freePartOptim->getEdge( i, start, end );

    for ( int j = 0; j < nConstraintsPerNode; ++j ) {
      lNew = infty;
      uNew = infty;

      if ( this->relativeLowerBounds[ nConstraintsPerNode * start + j ] < infty
          &&
          this->relativeLowerBounds[ nConstraintsPerNode * end + j ] < infty ) {
        lNew = 0.5 *
            ( this->relativeLowerBounds[ nConstraintsPerNode * start + j ] +
            this->relativeLowerBounds[ nConstraintsPerNode * end + j ] );
      }

      if ( this->relativeUpperBounds[ nConstraintsPerNode * start + j ] < infty
          &&
          this->relativeUpperBounds[ nConstraintsPerNode * end + j ] < infty ) {
        uNew = 0.5 *
            ( this->relativeUpperBounds[ nConstraintsPerNode * start + j ] +
            this->relativeUpperBounds[ nConstraintsPerNode * end + j ] );
      }

      this->relativeLowerBounds.push_back( lNew );
      this->relativeUpperBounds.push_back( uNew );
    }
  }
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::convergenceMonitor( ) {

  SCVT fRel =
      std::abs( ( this->cost - this->costOpt ) / this->cost );

  SCVT stopFRel = this->optimizer->get_ftol_rel( );

  std::cout << "|(f(x_new)-f(x_best))/f(x_new)| = " << fRel <<
      ", stopping criterion: " << stopFRel << "." << std::endl;

  if ( ( this->jump >= this->nloptJumpTol ) && ( fRel < stopFRel ) ) {
    std::cout << "STOP on f_tol_rel." << std::endl;
    this->optimizer->force_stop( );
    this->jump = 0;
  } else if ( fRel < stopFRel ) {
    ++this->jump;
  } else {
    this->jump = 0;
  }
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::setFreePartAnalAndInitMatrices( ) {
  this->subdivMatrices.reserve( analysisLevel );
  for ( int i = 0; i < analysisLevel; ++i ) {
    subdivMatrices.push_back( new SparseMatrix< LO, SCVT > );
  }
  OpenMeshWrapper< LO, SC > temp( *this->freePartOptim );
  temp.subdivideLoop( analysisLevel, true, &this->subdivMatrices );
  this->freePartAnal = new SurfaceMesh3D< LO, SC >( );
  temp.createSurfaceMesh3D( *this->freePartAnal );
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::updateFreePartAnal( ) {
  if ( this->freePartAnal ) delete this->freePartAnal;
  OpenMeshWrapper< LO, SC > temp( *this->freePartOptim );
  temp.subdivideLoop( this->analysisLevel - this->currOptimLevel );
  this->freePartAnal = new SurfaceMesh3D< LO, SC >( );
  temp.createSurfaceMesh3D( *this->freePartAnal );
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::updateFreePartOptimFromReference(
    const std::vector<double> & alpha
    ) {

  SCVT n[ 3 ];
  SCVT x[ 3 ];

  for ( LO i = 0; i < this->freePartOptimRef->getNNodes( ); ++i ) {
    this->freePartOptimRef->getNodalNormal( i, n );
    this->freePartOptimRef->getNode( i, x );
    x[ 0 ] += alpha[ i ] * n[ 0 ];
    x[ 1 ] += alpha[ i ] * n[ 1 ];
    x[ 2 ] += alpha[ i ] * n[ 2 ];
    this->freePartOptim->setNode( i, x );
  }
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::computeSubdivMatrix( ) {

  if ( this->subdivMatrix ) delete this->subdivMatrix;

  this->subdivMatrix = new SparseMatrix< int, SCVT >(
      *( this->subdivMatrices[ this->currOptimLevel ] ) );
  for ( int i = this->currOptimLevel + 1; i < this->analysisLevel; ++i ) {
    SparseMatrix< int, SCVT > tmp( *this->subdivMatrix );
    this->subdivMatrix->resize(
        this->subdivMatrices[ i ]->getNRows( ), tmp.getNCols( ) );
    this->subdivMatrix->multiply( *( this->subdivMatrices[ i ] ), tmp );
  }
}


//---IPOPT-METHODS------------------------------------------------------------//
#ifdef IPOPT

template<class LO, class SC>
bool MultiresolutionOptimizer<LO, SC>::get_nlp_info(
    Ipopt::Index & n,
    Ipopt::Index & m,
    Ipopt::Index & nnz_jac_g,
    Ipopt::Index & nnz_h_lag,
    IndexStyleEnum & index_style
    ) {

  n = this->freePartOptim->getNNodes( );
  m = 0;
  nnz_jac_g = 0;
  nnz_h_lag = 0;
  index_style = Ipopt::TNLP::C_STYLE;

  return true;
}

template<class LO, class SC>
bool MultiresolutionOptimizer<LO, SC>::get_bounds_info(
    Ipopt::Index n,
    Ipopt::Number * x_l,
    Ipopt::Number * x_u,
    Ipopt::Index m,
    Ipopt::Number * g_l,
    Ipopt::Number * g_u
    ) {

  assert( n == this->freePartOptim->getNNodes( ) );
  assert( m == 0 );

  std::vector< Ipopt::Number > lb, ub;
  this->setUpAbsoluteBounds( lb, ub );

  assert( x_l );
  assert( x_u );
  for ( Ipopt::Index i = 0; i < n; ++i ) {
    x_l[ i ] = lb[ i ];
    x_u[ i ] = ub[ i ];
  }

  return true;
}

template<class LO, class SC>
bool MultiresolutionOptimizer<LO, SC>::get_starting_point(
    Ipopt::Index n,
    bool init_x,
    Ipopt::Number * x,
    bool init_z,
    Ipopt::Number * z_L,
    Ipopt::Number * z_U,
    Ipopt::Index m,
    bool init_lambda,
    Ipopt::Number * lambda
    ) {

  assert( init_x == true );
  assert( init_z == false );
  assert( init_lambda == false );

  assert( x );
  for ( Ipopt::Index i = 0; i < n; ++i ) {
    x[ i ] = this->alphaOpt[ i ];
  }

  return true;
}

template<class LO, class SC>
bool MultiresolutionOptimizer<LO, SC>::eval_f(
    Ipopt::Index n,
    const Ipopt::Number * x,
    bool new_x,
    Ipopt::Number & obj_value
    ) {

  assert( n == this->freePartOptim->getNNodes( ) );

  if ( new_x ) {
    // updating free mesh on optim level
    std::vector< Ipopt::Number > alpha;
    alpha.assign( x, x + n );
    this->updateFreePartOptimFromReference( alpha );

    // subdivide to obtain current analysis mesh
    this->updateFreePartAnal( );

    // set current mesh for bvp
    this->bvp->setProblemData( this->freePartAnal, this->fixedPart );

    // solve the bvp
    if ( !this->bvp->solve( ) ) {
      std::cout << "BVP not solved!" << std::endl;
      obj_value = BVP_NOT_SOLVED_J;
      return true;
      //return false;
    }
  }

  bool renormalized;
  if ( this->normalizeOnEachLevel ) {
    renormalized = ( this->iterLevel == 0 );
  } else {
    renormalized = ( this->iter == 0 );
  }

  this->bvp->getCost( obj_value );
  if ( renormalized ) {
    this->costFirst = obj_value;
  }
  if ( this->normalizeCost ) {
    obj_value /= this->costFirst;
  }

  this->cost = obj_value;
  if ( obj_value < this->costOpt || renormalized ) {
    this->freePartOptim->getNodes( this->xOpt );
    this->alphaOpt.assign( x, x + n );
  }

  // print info and files
  this->bvp->printInfo( );
  std::cout << "\tJ = " << obj_value << "." << std::endl;

  if ( obj_value < this->costOpt || renormalized ) {
    this->costOpt = obj_value;
  }

  return true;
}

template<class LO, class SC>
bool MultiresolutionOptimizer<LO, SC>::eval_grad_f(
    Ipopt::Index n,
    const Ipopt::Number * x,
    bool new_x,
    Ipopt::Number * grad_f
    ) {

  assert( n == this->freePartOptim->getNNodes( ) );
  assert( grad_f );

  if ( new_x ) {
    // updating free mesh on optim level
    std::vector< Ipopt::Number > alpha;
    alpha.assign( x, x + n );
    this->updateFreePartOptimFromReference( alpha );

    // subdivide to obtain current analysis mesh
    this->updateFreePartAnal( );

    // set current mesh for bvp
    this->bvp->setProblemData( this->freePartAnal, this->fixedPart );

    // solve the bvp
    if ( !this->bvp->solve( ) ) {
      std::cout << "BVP not solved!" << std::endl;
      return false;
    }
  }

  bool renormalized;
  if ( this->normalizeOnEachLevel ) {
    renormalized = ( this->iterLevel == 0 );
  } else {
    renormalized = ( this->iter == 0 );
  }

  // get gradient
  Vector< LO, SC > gradV;
  this->getGradientNormalFixedReference( gradV );
  if ( renormalized ) {
    this->gradNormFirst = gradV.norm2( );
  }
  if ( this->normalizeGrad ) {
    gradV.scale( 1.0 / this->gradNormFirst );
  }

  for ( Ipopt::Index i = 0; i < gradV.getLength( ); ++i ) {
    grad_f[ i ] = gradV.get( i );
  }

  std::cout << "\t|g|_2 = " << gradV.norm2( ) << "." << std::endl;

  return true;
}

template<class LO, class SC>
bool MultiresolutionOptimizer<LO, SC>::eval_g(
    Ipopt::Index n,
    const Ipopt::Number * x,
    bool new_x,
    Ipopt::Index m,
    Ipopt::Number * g
    ) {

  assert( n == this->freePartOptim->getNNodes( ) );
  assert( m == 0 );

  return true;
}

template<class LO, class SC>
bool MultiresolutionOptimizer<LO, SC>::eval_jac_g(
    Ipopt::Index n,
    const Ipopt::Number * x,
    bool new_x,
    Ipopt::Index m,
    Ipopt::Index nele_jac,
    Ipopt::Index * iRow,
    Ipopt::Index * jCol,
    Ipopt::Number * values
    ) {

  assert( n == this->freePartOptim->getNNodes( ) );
  assert( m == 0 );

  return true;
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::finalize_solution(
    Ipopt::SolverReturn status,
    Ipopt::Index n,
    const Ipopt::Number * x,
    const Ipopt::Number * z_L,
    const Ipopt::Number * z_U,
    Ipopt::Index m,
    const Ipopt::Number * g,
    const Ipopt::Number * lambda,
    Ipopt::Number obj_value,
    const Ipopt::IpoptData * ip_data,
    Ipopt::IpoptCalculatedQuantities * ip_cq
    ) {

  // save xOpt
  SCVT normal[ 3 ];
  SCVT node[ 3 ];
  for ( LO i = 0; i < n; ++i ) {
    this->freePartOptimRef->getNodalNormal( i, normal );
    this->freePartOptimRef->getNode( i, node );
    for ( LO j = 0; j < 3; ++j ) {
      this->xOpt[ 3 * i + j ] = node[ j ] + x[ i ] * normal[ j ];
    }
  }

}

template<class LO, class SC>
bool MultiresolutionOptimizer<LO, SC>::intermediate_callback(
    Ipopt::AlgorithmMode mode,
    Ipopt::Index iter,
    Ipopt::Number obj_value,
    Ipopt::Number inf_pr,
    Ipopt::Number inf_du,
    Ipopt::Number mu,
    Ipopt::Number d_norm,
    Ipopt::Number regularization_size,
    Ipopt::Number alpha_du,
    Ipopt::Number alpha_pr,
    Ipopt::Index ls_trials,
    const Ipopt::IpoptData * ip_data,
    Ipopt::IpoptCalculatedQuantities * ip_cq
    ) {

  Ipopt::TNLPAdapter * tnlp_adapter = nullptr;
  if ( ip_cq ) {
    Ipopt::OrigIpoptNLP * orignlp;
    orignlp = dynamic_cast<Ipopt::OrigIpoptNLP*> (
        Ipopt::GetRawPtr( ip_cq->GetIpoptNLP( ) ) );
    if ( orignlp )
      tnlp_adapter = dynamic_cast<Ipopt::TNLPAdapter*> (
        Ipopt::GetRawPtr( orignlp->nlp( ) ) );
  }

  LO dim = this->freePartOptim->getNNodes( );
  Ipopt::Number * alpha = new Ipopt::Number[ dim ];
  tnlp_adapter->ResortX( *ip_data->curr( )->x( ), alpha );

  OpenMeshWrapper< LO, SC > mesh( *this->freePartOptimRef );
  SCVT n[ 3 ];
  SCVT x[ 3 ];
  for ( LO i = 0; i < dim; ++i ) {
    this->freePartOptimRef->getNodalNormal( i, n );
    this->freePartOptimRef->getNode( i, x );
    x[ 0 ] += alpha[ i ] * n[ 0 ];
    x[ 1 ] += alpha[ i ] * n[ 1 ];
    x[ 2 ] += alpha[ i ] * n[ 2 ];
    mesh.setNode( i, x );
  }

  std::stringstream file;
  std::stringstream iter2string;
  iter2string.width( 4 );
  iter2string.fill( '0' );
  iter2string << this->iter;
  file << this->outputPath << "optim_a_" << iter2string.str( ) << ".vtu";
  this->bvp->printVtu( file.str( ) );

  file.str( std::string( ) );
  file.clear( );
  file << this->outputPath << "optim_o_" << iter2string.str( ) << ".vtu";
  mesh.printParaviewVtu( file.str( ) );

  mesh.subdivideLoop( this->analysisLevel - this->currOptimLevel );
  file.str( std::string( ) );
  file.clear( );
  file << this->outputPath << "optim_s_" << iter2string.str( ) << ".vtu";
  mesh.printParaviewVtu( file.str( ) );

  ++this->iter;
  ++this->iterLevel;

  delete [] alpha;

  return true;
}

template<class LO, class SC>
MultiresolutionOptimizer<LO, SC> *
MultiresolutionOptimizer<LO, SC>::getShallowCopy( ) {

  MultiresolutionOptimizer< LO, SC > * copy =
      new MultiresolutionOptimizer< LO, SC >( );

  copy->iAmShallowCopy = true;

  copy->currOptimLevel = this->currOptimLevel;
  copy->analysisLevel = this->analysisLevel;
  copy->maxOptimLevel = this->maxOptimLevel;
  copy->outputPath = this->outputPath;
  copy->fixedPart = this->fixedPart;
  copy->freePartAnal = this->freePartAnal;
  copy->freePartOptim = this->freePartOptim;
  copy->freePartOptimRef = this->freePartOptimRef;
  copy->subdivMatrices = this->subdivMatrices;
  copy->subdivMatrix = this->subdivMatrix;
  copy->bvp = this->bvp;
  copy->relativeLowerBounds = this->relativeLowerBounds;
  copy->relativeUpperBounds = this->relativeUpperBounds;
  copy->motype = this->motype;
  copy->optimizer = this->optimizer;
  copy->nloptFRel = this->nloptFRel;
  copy->cost = this->cost;
  copy->costOpt = this->costOpt;
  copy->costFirst = this->costFirst;
  copy->gradNormFirst = this->gradNormFirst;
  copy->normalizeCost = this->normalizeCost;
  copy->normalizeGrad = this->normalizeGrad;
  copy->normalizeOnEachLevel = this->normalizeOnEachLevel;
  copy->xOpt = this->xOpt;
  copy->alphaOpt = this->alphaOpt;
  copy->iter = this->iter;
  copy->iterLevel = this->iterLevel;
  copy->jump = this->jump;
  copy->nloptJumpTol = this->nloptJumpTol;

  copy->ipopt_tol = this->ipopt_tol;
  copy->ipopt_acceptable_iter = this->ipopt_acceptable_iter;
  copy->ipopt_acceptable_obj_change_tol = this->ipopt_acceptable_obj_change_tol;
  copy->ipopt_max_iter = this->ipopt_max_iter;
  copy->ipopt_acceptable_tol = this->ipopt_acceptable_tol;

  return copy;
}

template<class LO, class SC>
void MultiresolutionOptimizer<LO, SC>::shallowCopy(
    MultiresolutionOptimizer<LO, SC> & copy
    ) const {

  copy.currOptimLevel = this->currOptimLevel;
  copy.analysisLevel = this->analysisLevel;
  copy.maxOptimLevel = this->maxOptimLevel;
  copy.outputPath = this->outputPath;
  copy.fixedPart = this->fixedPart;
  copy.freePartAnal = this->freePartAnal;
  copy.freePartOptim = this->freePartOptim;
  copy.freePartOptimRef = this->freePartOptimRef;
  copy.subdivMatrices = this->subdivMatrices;
  copy.subdivMatrix = this->subdivMatrix;
  copy.bvp = this->bvp;
  copy.relativeLowerBounds = this->relativeLowerBounds;
  copy.relativeUpperBounds = this->relativeUpperBounds;
  copy.motype = this->motype;
  copy.optimizer = this->optimizer;
  copy.nloptFRel = this->nloptFRel;
  copy.cost = this->cost;
  copy.costOpt = this->costOpt;
  copy.costFirst = this->costFirst;
  copy.gradNormFirst = this->gradNormFirst;
  copy.normalizeCost = this->normalizeCost;
  copy.normalizeGrad = this->normalizeGrad;
  copy.normalizeOnEachLevel = this->normalizeOnEachLevel;
  copy.xOpt = this->xOpt;
  copy.alphaOpt = this->alphaOpt;
  copy.iter = this->iter;
  copy.iterLevel = this->iterLevel;
  copy.jump = this->jump;
  copy.nloptJumpTol = this->nloptJumpTol;

  copy.ipopt_tol = this->ipopt_tol;
  copy.ipopt_acceptable_iter = this->ipopt_acceptable_iter;
  copy.ipopt_acceptable_obj_change_tol = this->ipopt_acceptable_obj_change_tol;
  copy.ipopt_max_iter = this->ipopt_max_iter;
  copy.ipopt_acceptable_tol = this->ipopt_acceptable_tol;

}

#endif

} // end namespace bem4i

#endif /*MULTIRESOLUTIONOPTIMIZER_H*/
