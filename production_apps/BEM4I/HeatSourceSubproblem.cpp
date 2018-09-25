/*!
 * @file    HeatSourceSubproblem.cpp
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    September 1, 2014
 * 
 */

#ifdef HEATSOURCESUBPROBLEM_H

namespace bem4i {

template< class LO, class SC >
HeatSourceSubproblem< LO, SC >::HeatSourceSubproblem( ) {
}

template< class LO, class SC >
HeatSourceSubproblem< LO, SC >::HeatSourceSubproblem(
    const Vector< LO, SC > & targetNeumannDataSensor
    ) {

  this->cost = 0.0;
  this->costMultiplier = 1.0;
  this->regularizationPar = 0.0;

  this->targetNeumannDataSensor =
      new Vector< LO, SC >( targetNeumannDataSensor );

  this->meshSensor = nullptr;
  this->meshSource = nullptr;
  this->primalDirichletData = nullptr;
  this->primalNeumannData = nullptr;
  this->adjointDirichletData = nullptr;
  this->adjointNeumannData = nullptr;
  this->adjointOnSource = nullptr;
}

template< class LO, class SC >
HeatSourceSubproblem< LO, SC >::~HeatSourceSubproblem( ) {
  if ( this->primalDirichletData ) delete this->primalDirichletData;
  if ( this->primalNeumannData ) delete this->primalNeumannData;
  if ( this->adjointDirichletData ) delete this->adjointDirichletData;
  if ( this->adjointNeumannData ) delete this->adjointNeumannData;
  if ( this->targetNeumannDataSensor ) delete this->targetNeumannDataSensor;
  if ( this->meshSensor ) delete this->meshSensor;
  if ( this->meshSource ) delete this->meshSource;
  if ( this->adjointOnSource ) delete this->adjointOnSource;
}

template< class LO, class SC >
void HeatSourceSubproblem<LO, SC>::printInfo( ) const {
  std::cout << "Class HeatSourceSubproblem, " << std::endl;
  if ( this->meshSource ) {
    std::cout << "\t Mesh source ";
    this->meshSource->printInfo( );
  }
  if ( this->meshSensor ) {
    std::cout << "\t Mesh sensor ";
    this->meshSensor->printInfo( );
  }
//  std::cout << "\tJ = " << this->cost << "." << std::endl;
}

template<class LO, class SC>
void HeatSourceSubproblem<LO, SC>::printVtu(
    const string & filename
    ) const {

  std::vector< string > nodeNames;
  std::vector< string > elemNames;
  std::vector< Vector< LO, SC > * > nodalData;
  std::vector< Vector< LO, SC > * > elemData;

  if ( this->primalDirichletData ) {
    nodeNames.push_back( "Dirichlet" );
    nodalData.push_back( this->primalDirichletData );
  }

  if ( this->adjointDirichletData ) {
    nodeNames.push_back( "Dirichlet_adjoint" );
    nodalData.push_back( this->adjointDirichletData );
  }

  if ( this->primalNeumannData ) {
    elemNames.push_back( "Neumann" );
    elemData.push_back( this->primalNeumannData );
  }

  if ( this->adjointNeumannData ) {
    elemNames.push_back( "Neumann_adjoint" );
    elemData.push_back( this->adjointNeumannData );
  }

  if ( this->targetNeumannDataSensor ) {
    elemNames.push_back( "Neumann_target" );
    elemData.push_back( this->targetNeumannDataSensor );
  }

  this->meshSensor->printParaviewVtu( filename, &nodeNames, &nodalData,
      &elemNames, &elemData );

}

template< class LO, class SC >
bool HeatSourceSubproblem<LO, SC>::solve( ) {

  const SCVT CGeps = 1e-12;
  const LO maxIter = 5000;
  const int order = 3;
  int quadOrder [4] = { order, order, order, order };
  const quadratureType quadrature = SauterSchwab;
  const int maxElemPerCluster = 40;
  const SCVT eta = 1.0;
  const SCVT epsilon = 1e-5;

  Tree<BECluster<LO, SC>*, LO> tree, tree2;
  ProgressMonitor::init( "Nested dissection" );
  this->meshSensor->nestedDissection( tree, maxElemPerCluster );
  this->meshSensor->nestedDissection( tree2, maxElemPerCluster );
  ProgressMonitor::step( );

  FastBESpace< LO, SC > bespace00( this->meshSensor, p0, p0, eta );
  bespace00.setEpsilonACA( epsilon );
  FastBESpace< LO, SC > bespace10( this->meshSensor, p1, p0, eta );
  bespace10.setEpsilonACA( epsilon );
  BESpace< LO, SC > bespace11( this->meshSensor, p1, p1 );

  std::cout << "Assembling V" << std::endl;
  FullMatrix< LO, SC > * V = new FullMatrix< LO, SC > ( 0, 0 );
  BEBilinearFormLaplace1Layer< LO, SC > formV( &bespace00, quadOrder,
      quadrature );
  formV.assemble( *V );

  std::cout << "Assembling K" << std::endl;
  FullMatrix< LO, SC > * K = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLaplace2Layer< LO, SC > formK( &bespace10, quadOrder,
      quadrature );
  formK.assemble( *K );

  IdentityOperator< LO, SC > id01( &bespace10 );

  Vector< LO, SC > rhs;
  Vector< LO, SC > aux( this->nElemsFixed );

  this->setUpPrimalRHS( rhs );

  ProgressMonitor::init( "Solving the primal system" );
  if ( !V->CGSolve( rhs, *this->primalNeumannData, CGeps, maxIter ) )
    return false;
  ProgressMonitor::step( );

  Vector< LO, SC > primalNeumannMinusTarget( *this->primalNeumannData );
  primalNeumannMinusTarget.add( *this->targetNeumannDataSensor, -1.0 );

  ProgressMonitor::init( "Computing cost" );
  this->cost = 0.0;
  // add L2 norm
  for ( LO i = 0; i < this->nElemsFixed; ++i ) {
    this->cost += 0.5 * ( primalNeumannMinusTarget.get( i ) *
        primalNeumannMinusTarget.get( i ) ) *
        this->meshSensor->getElemArea( i );
  }
  // add free surface
  for ( LO i = 0; i < this->nElemsFree; ++i ) {
    this->cost += this->regularizationPar * this->meshSource->getElemArea( i );
  }

  this->cost *= this->costMultiplier;
  ProgressMonitor::step( );

  ProgressMonitor::init(
      "Computing nodal Dirichlet BC [du/dn] for the adjoint" );
  this->elem2NodalAreaWeighted( primalNeumannMinusTarget,
      *this->adjointDirichletData );
  ProgressMonitor::step( );

  Vector< LO, SC > rhsAdjoint( this->nElemsFixed );

  ProgressMonitor::init( "Setting up the rhs for the adjoint problem" );
  id01.apply( *this->adjointDirichletData, aux, false, 0.5, 0.0 );
  K->apply( *this->adjointDirichletData, rhsAdjoint );
  rhsAdjoint.add( aux );
  ProgressMonitor::step( );

  ProgressMonitor::init( "Solving the adjoint system" );
  if ( !V->CGSolve( rhsAdjoint, *this->adjointNeumannData, CGeps, maxIter ) )
    return false;
  ProgressMonitor::step( );

  RepresentationFormulaLaplace<LO, SC> formula( &bespace10,
      this->adjointDirichletData, this->adjointNeumannData );
  ProgressMonitor::init( "Evaluating adjoint on source" );
  formula.evaluate( *this->meshSource, true, *this->adjointOnSource );
  ProgressMonitor::step( );

  delete V;
  delete K;

  return true;
}

template<class LO, class SC>
void HeatSourceSubproblem<LO, SC>::setUpPrimalRHS(
    Vector<LO, SC> & rhs
    ) const {
  rhs.resize( this->nElemsFixed );

  int quadOrder = 5;
  int quadSize = quadSizes[ quadOrder ];

#ifdef VERBOSE
  ProgressMonitor::init( "Setting up the rhs for the primal problem",
      this->nElemsFixed );
#endif

#pragma omp parallel shared( rhs )
  {
    SC val;
    SC valInner;
    SCVT x1[ 3 ], x2[ 3 ], x3[ 3 ];
    SCVT y1[ 3 ], y2[ 3 ], y3[ 3 ];
    SCVT * x = new SCVT[ quadSize * 3 ];
    SCVT * y = new SCVT[ quadSize * 3 ];
    SCVT ny[ 3 ];

#pragma omp for
    for ( LO outerEl = 0; outerEl < this->nElemsFixed; ++outerEl ) {
      val = 0.0;
      this->meshSensor->getNodes( outerEl, x1, x2, x3 );
      this->getQuadratureNodes( x1, x2, x3, quadOrder, x );

      // go over integration nodes in tau_outerEl
      for ( int outerPoint = 0; outerPoint < quadSize; ++outerPoint ) {

        // go over inner elements
        for ( LO innerEl = 0; innerEl < this->nElemsFree; ++innerEl ) {
          valInner = 0.0;
          this->meshSource->getNodes( innerEl, y1, y2, y3 );
          this->getQuadratureNodes( y1, y2, y3, quadOrder, y );
          this->meshSource->getNormal( innerEl, ny );

          // go over integration nodes in in tau_innerEl
          for ( int innerPoint = 0; innerPoint < quadSize; ++innerPoint ) {
            valInner += quadWeights[ quadOrder ][ innerPoint ] *
                this->evalNewtonKernel( x + 3 * outerPoint, y + 3 * innerPoint,
                ny );
          }

          val += quadWeights[ quadOrder ][ outerPoint ] *
              valInner * this->meshSource->getElemArea( innerEl );
        }
      }

      val *= this->meshSensor->getElemArea( outerEl );
      rhs.set( outerEl, val );

#ifdef VERBOSE
      ProgressMonitor::step( );
#endif
    }

    delete x;
    delete y;
  }
}

template<class LO, class SC>
void HeatSourceSubproblem<LO, SC>::getQuadratureNodes(
    const SCVT * x1,
    const SCVT * x2,
    const SCVT * x3,
    int quadratureOrder,
    SCVT * nodes )
const {
  int numPoints = quadSizes[ quadratureOrder ];
  for ( int i = 0; i < numPoints; i++ ) {
    nodes[ i * 3 ] = x1[ 0 ] + ( x2[ 0 ] - x1[ 0 ] ) *
        quadPoints[ quadratureOrder ][ i * 2 ] + ( x3[ 0 ] - x1[ 0 ] ) *
        quadPoints[ quadratureOrder ][ i * 2 + 1 ];
    nodes[ i * 3 + 1 ] = x1[ 1 ] + ( x2[ 1 ] - x1[ 1 ] ) *
        quadPoints[ quadratureOrder ][ i * 2 ] + ( x3[ 1 ] - x1[ 1 ] ) *
        quadPoints[ quadratureOrder ][ i * 2 + 1 ];
    nodes[ i * 3 + 2 ] = x1[ 2 ] + ( x2[ 2 ] - x1[ 2 ] ) *
        quadPoints[ quadratureOrder ][ i * 2 ] + ( x3[ 2 ] - x1[ 2 ] ) *
        quadPoints[ quadratureOrder ][ i * 2 + 1 ];
  }
}

template<class LO, class SC>
SC HeatSourceSubproblem<LO, SC>::evalNewtonKernel(
    const SCVT * x,
    const SCVT * y,
    const SCVT * n
    ) const {
  SCVT norm = std::sqrt( ( x[0] - y[0] )*( x[0] - y[0] ) +
      ( x[1] - y[1] )*( x[1] - y[1] ) +
      ( x[2] - y[2] )*( x[2] - y[2] ) );
  SCVT dot = ( x[0] - y[0] ) * n[0] +
      ( x[1] - y[1] ) * n[1] +
      ( x[2] - y[2] ) * n[2];
  return ( ( PI_FACT / 2.0 ) * dot / norm );
}

template< class LO, class SC >
void HeatSourceSubproblem<LO, SC>::getCost( SCVT & cost ) const {
  cost = this->cost;
}

template< class LO, class SC >
void HeatSourceSubproblem<LO, SC>::getShapeGradient(
    const Vector< LO, SCVT > & perturbation,
    SCVT & dx1,
    SCVT & dx2,
    SCVT & dx3
    ) const {

  dx1 = dx2 = dx3 = 0.0;
  LO elem[ 3 ];
  SCVT n[ 3 ];
  SCVT area;
  SCVT elemPert;
  SC g;
  if ( this->regularizationPar > 0 && !this->meshSource->curvatureReady( ) )
    this->meshSource->initAdditiveCurvature( );

  for ( LO i = 0; i < this->nElemsFree; ++i ) {
    this->meshSource->getElement( i, elem );

    elemPert = ( perturbation.get( elem[ 0 ] ) +
        perturbation.get( elem[ 1 ] ) +
        perturbation.get( elem[ 2 ] ) ) / 3.0;

    if ( std::abs( elemPert ) < EPS ) continue;

    area = this->meshSource->getElemArea( i );

    this->meshSource->getNormal( i, n );

    // add -p
    g = -( this->adjointOnSource->get( elem[ 0 ] ) +
        this->adjointOnSource->get( elem[ 1 ] ) +
        this->adjointOnSource->get( elem[ 2 ] ) ) / 3.0;

    // add eps*H 
    g += this->regularizationPar * (
        this->meshSource->getCurvature( elem[ 0 ] ) +
        this->meshSource->getCurvature( elem[ 1 ] ) +
        this->meshSource->getCurvature( elem[ 2 ] )
        ) / 3.0;

    g *= this->costMultiplier;

    dx1 += area * g * elemPert * n[ 0 ];
    dx2 += area * g * elemPert * n[ 1 ];
    dx3 += area * g * elemPert * n[ 2 ];
  }

}

template< class LO, class SC >
void HeatSourceSubproblem<LO, SC>::getShapeGradient(
    const Vector< LO, SCVT > & perturbationX1,
    const Vector< LO, SCVT > & perturbationX2,
    const Vector< LO, SCVT > & perturbationX3,
    SCVT & dx
    ) const {

  dx = 0.0;
  LO elem[ 3 ];
  SCVT n[ 3 ];
  SCVT area;
  SCVT elemPertX1, elemPertX2, elemPertX3;
  SC g;
  if ( this->regularizationPar > 0 && !this->meshSource->curvatureReady( ) )
    this->meshSource->initAdditiveCurvature( );

  for ( LO i = 0; i < this->nElemsFree; ++i ) {

    this->meshSource->getElement( i, elem );

    elemPertX1 = ( perturbationX1.get( elem[ 0 ] ) +
        perturbationX1.get( elem[ 1 ] ) +
        perturbationX1.get( elem[ 2 ] ) ) / 3.0;
    elemPertX2 = ( perturbationX2.get( elem[ 0 ] ) +
        perturbationX2.get( elem[ 1 ] ) +
        perturbationX2.get( elem[ 2 ] ) ) / 3.0;
    elemPertX3 = ( perturbationX3.get( elem[ 0 ] ) +
        perturbationX3.get( elem[ 1 ] ) +
        perturbationX3.get( elem[ 2 ] ) ) / 3.0;

    if ( std::sqrt( elemPertX1 * elemPertX1 + elemPertX2 * elemPertX2 + elemPertX3 *
        elemPertX3 ) < EPS ) continue;

    area = this->meshSource->getElemArea( i );

    this->meshSource->getNormal( i, n );

    // add p
    g = -( this->adjointOnSource->get( elem[ 0 ] ) +
        this->adjointOnSource->get( elem[ 1 ] ) +
        this->adjointOnSource->get( elem[ 2 ] ) ) / 3.0;

    // add eps*H 
    g += this->regularizationPar * (
        this->meshSource->getCurvature( elem[ 0 ] ) +
        this->meshSource->getCurvature( elem[ 1 ] ) +
        this->meshSource->getCurvature( elem[ 2 ] )
        ) / 3.0;

    g *= this->costMultiplier;

    dx += area * g * ( elemPertX1 * n[ 0 ] +
        elemPertX2 * n[ 1 ] + elemPertX3 * n[ 2 ] );
  }

}

template<class LO, class SC>
void HeatSourceSubproblem<LO, SC>::getShapeGradient(
    Vector< LO, SC > & grad
    ) const {

  grad.resize( this->nNodesFree );
  if ( this->regularizationPar > 0 && !this->meshSource->curvatureReady( ) )
    this->meshSource->initAdditiveCurvature( );

  for ( LO i = 0; i < this->nNodesFree; ++i ) {
    grad.set( i, ( -this->adjointOnSource->get( i ) +
        this->regularizationPar * this->meshSource->getCurvature( i ) ) *
        this->costMultiplier );
  }
}

template<class LO, class SC>
void HeatSourceSubproblem<LO, SC>::setProblemData(
    const SurfaceMesh3D< LO, SC > * freeMesh,
    const SurfaceMesh3D< LO, SC > * fixedMesh
    ) {

  if ( fixedMesh ) {
    this->nNodesFixed = fixedMesh->getNNodes( );
    this->nElemsFixed = fixedMesh->getNElements( );
  } else {
    this->nNodesFixed = 0;
    this->nElemsFixed = 0;
  }

  this->nNodesFree = freeMesh->getNNodes( );
  this->nElemsFree = freeMesh->getNElements( );

  if ( this->meshSensor ) delete this->meshSensor;
  if ( this->meshSource ) delete this->meshSource;

  this->meshSource = new SurfaceMesh3D< LO, SC >( *freeMesh );
  if ( fixedMesh ) {
    this->meshSensor = new SurfaceMesh3D< LO, SC >( *fixedMesh );
  }

  // set primal Dirichlet data to zero
  if ( this->primalDirichletData ) delete this->primalDirichletData;
  this->primalDirichletData = new Vector< LO, SC >( nNodesFixed, true );

  if ( this->primalNeumannData ) delete this->primalNeumannData;
  this->primalNeumannData = new Vector< LO, SC >( nElemsFixed );

  if ( this->adjointDirichletData ) delete this->adjointDirichletData;
  this->adjointDirichletData = new Vector< LO, SC >( nNodesFixed );

  if ( this->adjointNeumannData ) delete this->adjointNeumannData;
  this->adjointNeumannData = new Vector< LO, SC >( nElemsFixed );

  if ( this->adjointOnSource ) delete this->adjointOnSource;
  this->adjointOnSource = new Vector< LO, SC >( nNodesFree );
}

//---PRIVATE-METHODS----------------------------------------------------------//
// todo: move to SurfaceMesh3d

template<class LO, class SC>
void HeatSourceSubproblem<LO, SC>::elem2NodalAreaWeighted(
    const Vector<LO, SC> & elem,
    Vector<LO, SC> & nodal
    ) const {

  nodal.resize( this->nNodesFixed, true );
  SCVT areaSum;
  SCVT sum;
  SCVT area;
  vector< LO > elems;

  for ( LO i = 0; i < this->nNodesFixed; ++i ) {
    this->meshSensor->getElements( i, elems );
    areaSum = 0;
    sum = 0;

    for ( int j = 0; j < elems.size( ); ++j ) {
      area = this->meshSensor->getElemArea( elems[ j ] );
      sum += area * elem.get( elems[ j ] );
      areaSum += area;
    }

    nodal.set( i, sum / areaSum );
  }

}


} // end namespace bem4i

#endif /* HEATSOURCESUBPROBLEM_H */
