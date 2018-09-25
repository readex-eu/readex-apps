#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../ProgressMonitor.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../FastBESpace.h"
#include "../../Vector.h"
#include "../../BEIntegratorHelmholtz.h"
#include "../../BEBilinearFormHelmholtz1Layer.h"
#include "../../BEBilinearFormHelmholtz2Layer.h"
#include "../../IdentityOperator.h"
#include "../../RepresentationFormulaHelmholtz.h"

using namespace std;
using namespace bem4i;

void testHelmholtzSoftScatterACA(
    string const &filename
    );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  string filename = "input/cube_12.txt";
  testHelmholtzSoftScatterACA( filename );

  return 0;
}

void testHelmholtzSoftScatterACA(
    string const &filename
    ) {
  typedef complex<double> SC;
  typedef int LO;
  typedef typename SC::value_type SCVT;

  std::cout.precision( 8 );
  std::cout.setf( std::ios::scientific );

  const SCVT GMRESPrec = 1e-6;
  const int GMRESIter = 2000;
  const SCVT ACAPrec = 1e-6;
  const int GMRESRestart = 1000;
  const SC kappa = 10.0;
  const int order = 4;
  const quadratureType qType = SauterSchwab;

  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  //mesh.scale(2.0);
  //mesh.scale(1.5);
  //mesh.move( 1.5, 1.0, 0.0 );
  mesh.refine( 3 );
  //mesh.mapToUnitBall( );
  mesh.printInfo( );

  Tree<BECluster<LO, SC>*, LO> tree, tree2;
  ProgressMonitor::init( "Nested dissection" );
  mesh.nestedDissection( tree, 20 );
  mesh.nestedDissection( tree2, 20 );
  ProgressMonitor::step( );

  int quadrature[4] = { order, order, order, order };

  FastBESpace< LO, SC > bespaceK( &mesh, p1, p0, &tree, 1.0, 8, 5, 0 );
  bespaceK.setEpsilonACA( ACAPrec );
  ACAMatrix< LO, SC > *K = new ACAMatrix<LO, SC>( 0, 0 );
  BEBilinearFormHelmholtz2Layer< LO, SC > formK( &bespaceK, quadrature, kappa,
      qType );
  formK.assemble( *K );
  std::cout << "compression rate: " << K->getCompressionRatio( ) << std::endl;

  IdentityOperator< LO, SC > id( &bespaceK );

  SC iUnit( 0.0, 1.0 );
  LO nNodes = mesh.getNNodes( );
  LO nElems = mesh.getNElements( );
  Vector< LO, SC > dir( nNodes );
  Vector< LO, SC > neu( nElems );
  Vector< LO, SC > sol( nElems );
  Vector< LO, SC > aux( nElems );
  SCVT x[ 3 ], d[ 3 ];
  d[ 0 ] = //-1.0; ///std::sqrt(2.0);
      d[ 1 ] = -1.0; //1.0 / std::sqrt( 3.0 );
  d[ 2 ] = 0.0; //1.0 / std::sqrt( 3.0 );

  for ( LO i = 0; i < nNodes; i++ ) {
    mesh.getNode( i, x );
    dir.set( i, -std::exp( iUnit * kappa * DOT3( x, d ) ) );
  }

  ProgressMonitor::init( "Setting up rhs" );
  id.apply( dir, aux, false, -0.5, 0.0 );
  K->apply( dir, neu );
  delete K;
  neu.add( aux );
  ProgressMonitor::step( );

  FastBESpace< LO, SC > bespaceV( &mesh, p0, p0, &tree2, 1.0, 8, 5, 0 );
  bespaceV.setEpsilonACA( ACAPrec );
  ACAMatrix< LO, SC > V( 0, 0 );
  BEBilinearFormHelmholtz1Layer< LO, SC > formV( &bespaceV, quadrature, kappa,
      qType );
  formV.assemble( V );
  std::cout << "compression rate: " << V.getCompressionRatio( ) << std::endl;

  ProgressMonitor::init( "Solving the system" );
  LeftIdentityPreconditioner<LO, SC>* M =
      new LeftIdentityPreconditioner<LO, SC>;
  V.GMRESSolve( neu, sol, GMRESPrec, GMRESIter, GMRESRestart, M );
  ProgressMonitor::step( );

  Vector< LO, SC > incident( dir );
  incident.scale( -1.0 );
  Vector< LO, SC > total( nNodes );
  total.setAll( 0.0 );
  string nodeNames[] = { "dirichlet_scatter", "incident", "total" };
  string elemNames[] = { "neumann_scatter" };
  Vector< LO, SC >* nodalData[] = { &dir, &incident, &total };
  Vector< LO, SC >* elemData[] = { &sol };
  mesh.printParaviewVtu( "output/scatter/output_soft_scatter.vtu", 3, nodeNames,
      nodalData, 1, elemNames, elemData );

  string file_scatter = "input/scatter_grid_x.txt";
  SurfaceMesh3D< LO, SC > mesh_scatter;
  mesh_scatter.load( file_scatter.c_str( ) );
  mesh_scatter.scale( 3.0 );
  mesh_scatter.refine( 6 );
  mesh_scatter.printInfo( );
  LO nPoints = mesh_scatter.getNNodes( );
  Vector< LO, SC > pointSolutionScatter( nPoints );
  ProgressMonitor::init( "Evaluating representation formula" );
  RepresentationFormulaHelmholtz<LO, SC> formula( &bespaceK, &dir, &sol,
      kappa );
  formula.evaluate( mesh_scatter, false, pointSolutionScatter );
  ProgressMonitor::step( );
  Vector< LO, SC > pointSolutionInc( nPoints );
  Vector< LO, SC > pointSolutionTotal( nPoints );
  for ( LO i = 0; i < nPoints; i++ ) {
    mesh_scatter.getNode( i, x );
    pointSolutionInc.set( i, std::exp( iUnit * kappa * DOT3( x, d ) ) );
    pointSolutionTotal.set( i, pointSolutionScatter.get( i ) +
        pointSolutionInc.get( i ) );
  }
  string nodeNamesSol[] = { "incident", "scatter", "total" };
  Vector< LO, SC >* nodalDataSol[] = { &pointSolutionInc, &pointSolutionScatter,
    &pointSolutionTotal };
  mesh_scatter.printParaviewVtu( "output/scatter/solution_soft_scatter.vtu", 3,
      nodeNamesSol, nodalDataSol, 0, nullptr, nullptr );

  SCVT velocity = 350.0;
  SCVT endTime = 2.0 * M_PI / ( velocity * kappa.real( ) );
  int nSteps = 50;
  SCVT timeStep = endTime / nSteps;
  Vector< LO, SC > *pointSolutionTotalTime;
  std::stringstream file;
  file.fill( '0' );
  file.width( 4 );
  for ( int i = 0; i < nSteps; i++ ) {

    pointSolutionTotalTime = new Vector< LO, SC >( pointSolutionTotal );
    pointSolutionTotalTime->scale( std::exp( -iUnit * velocity * kappa * (SCVT) i *
        timeStep ) );
    nodalDataSol[ 0 ] = pointSolutionTotalTime;
    nodeNamesSol[ 0 ] = "total";
    file.str( std::string( ) );
    file.clear( );
    file << "output/scatter/solution_soft_scatter_2_" << i << ".vtu";
    mesh_scatter.printParaviewVtu( file.str( ).c_str( ), 1, nodeNamesSol,
        nodalDataSol, 0, nullptr, nullptr );
    delete pointSolutionTotalTime;
  }

}

