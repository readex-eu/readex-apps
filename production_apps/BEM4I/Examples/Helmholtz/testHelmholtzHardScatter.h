#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../ProgressMonitor.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../FastBESpace.h"
#include "../../Vector.h"
#include "../../BEIntegratorHelmholtz.h"
#include "../../BEBilinearFormHelmholtz2Layer.h"
#include "../../BEBilinearFormHelmholtzHypersingular.h"
#include "../../IdentityOperator.h"
#include "../../RepresentationFormulaHelmholtz.h"

using namespace std;
using namespace bem4i;

void testHelmholtzHardScatter(
    string const &filename,
    int printPrecision
    );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  string filename = "input/cube_12.txt";
  int printPrecision = 15;
  testHelmholtzHardScatter( filename, printPrecision );

  return 0;
}

void testHelmholtzHardScatter(
    string const &filename,
    int printPrecision
    ) {

  typedef complex<double> SC;
  typedef int LO;
  typedef typename SC::value_type SCVT;
  timeval start, stop;
  std::cout.setf( std::ios::showpoint | std::ios::fixed );
  std::cout.precision( printPrecision );

  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  //mesh.scale( 2.0 );
  //mesh.move( 1.5, 1.0, 0.0 );
  mesh.refine( 3 );
  //mesh.mapToUnitBall( );

  //  SurfaceMesh3D< LO, SC > mesh2;
  //  mesh2.load( "input/ball_icosahedron_20.txt" );
  //  mesh2.refine( 3 );
  //  mesh2.mapToUnitBall();
  //  mesh2.printInfo();
  //  mesh2.move( -2.0, -1.0, 0.0 );
  //
  //  mesh.append( mesh2 );
  mesh.printInfo( );

  int quadDisj[2] = { 3, 3 };

  SCVT lambda = 0.1;
  SCVT kappa = 2 * M_PI / lambda;
  kappa = 10.0;

  BESpace< LO, SC > bespaceK( &mesh, p1, p0 );
  FullMatrix< LO, SC > *K = new FullMatrix<LO, SC>( 0, 0 );
  BEBilinearFormHelmholtz2Layer< LO, SC > formK( &bespaceK, nullptr, kappa, SauterSchwab, quadDisj );
  std::cout << "Assembling K ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  formK.assemble( *K );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  IdentityOperator< LO, SC > id( &bespaceK );

  SC iUnit( 0.0, 1.0 );
  LO nNodes = mesh.getNNodes( );
  LO nElems = mesh.getNElements( );
  Vector< LO, SC > dir( nNodes );
  Vector< LO, SC > neu( nElems );
  Vector< LO, SC > aux( nElems );
  SCVT x[ 3 ], d[ 3 ], n[ 3 ], x1[3], x2[3], x3[3];
  d[ 0 ] = -1.0 / std::sqrt( 2.0 ); ///std::sqrt(2.0);//0.0;
  d[ 1 ] = 0.0; //1.0 / std::sqrt( 3.0 );
  d[ 2 ] = -1.0 / std::sqrt( 2.0 ); // / std::sqrt( 2.0 );

  for ( LO i = 0; i < nElems; i++ ) {
    mesh.getNormal( i, n );
    mesh.getNodes( i, x1, x2, x3 );
    x[0] = 1.0 / 3.0 * ( x1[0] + x2[0] + x3[0] );
    x[1] = 1.0 / 3.0 * ( x1[1] + x2[1] + x3[1] );
    x[2] = 1.0 / 3.0 * ( x1[2] + x2[2] + x3[2] );
    neu.set( i, -iUnit * kappa * DOT3( x, n ) * std::exp( iUnit * kappa * DOT3( x, d ) ) );
  }

  std::cout << "Setting up the rhs ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  id.apply( neu, aux, true, -0.5, 0.0 );
  K->apply( neu, dir, true, -1.0, 0.0 );
  dir.add( aux );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  delete K;

  BESpace< LO, SC > bespaceD( &mesh, p1, p1 );
  FullMatrix< LO, SC > * D = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtzHypersingular< LO, SC > formD( &bespaceD, nullptr, kappa, SauterSchwab, quadDisj );
  std::cout << "Assembling D ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  formD.assemble( *D );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  std::cout << "Solving the system ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  D->LUSolve( dir );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  delete D;

  Vector< LO, SC > incident( nNodes );
  for ( LO i = 0; i < nNodes; i++ ) {
    mesh.getNode( i, x );
    incident.set( i, std::exp( iUnit * kappa * DOT3( x, d ) ) );
  }
  Vector< LO, SC > total( nNodes );
  incident.add( dir, total );
  string nodeNames[] = { "dirichlet_scatter", "incident", "total" };
  string elemNames[] = { "neumann_scatter" };
  Vector< LO, SC >* nodalData[] = { &dir, &incident, &total };
  Vector< LO, SC >* elemData[] = { &neu };
  mesh.printParaviewVtu( "output/surface_hard_scatter.vtu", 3, nodeNames, nodalData, 1, elemNames, elemData );

  string file_scatter = "input/scatter_grid_y.txt";
  SurfaceMesh3D< LO, SC > mesh_scatter;
  mesh_scatter.load( file_scatter.c_str( ) );
  mesh_scatter.printInfo( );
  mesh_scatter.scale( 2.0 );
  mesh_scatter.refine( 6 );
  mesh_scatter.printInfo( );
  LO nPoints = mesh_scatter.getNNodes( );
  Vector< LO, SC > pointSolutionScatter( nPoints );
  std::cout << "Evaluating representation formula ... ";
  std::cout.flush( );
  RepresentationFormulaHelmholtz<LO, SC> formula( &bespaceK, &dir, &neu, kappa );
  gettimeofday( &start, nullptr );
  formula.evaluate( mesh_scatter, false, pointSolutionScatter );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;
  Vector< LO, SC > pointSolutionInc( nPoints );
  Vector< LO, SC > pointSolutionTotal( nPoints );
  for ( LO i = 0; i < nPoints; i++ ) {
    mesh_scatter.getNode( i, x );
    pointSolutionInc.set( i, std::exp( iUnit * kappa * DOT3( x, d ) ) );
    pointSolutionTotal.set( i, pointSolutionScatter.get( i ) + pointSolutionInc.get( i ) );
  }
  string nodeNamesSol[] = { "incident", "scatter", "total" };
  Vector< LO, SC >* nodalDataSol[] = { &pointSolutionInc, &pointSolutionScatter, &pointSolutionTotal };
  mesh_scatter.printParaviewVtu( "output/domain_hard_scatter.vtu", 3, nodeNamesSol, nodalDataSol, 0, nullptr, nullptr );

  //return;

  SCVT velocity = 350.0;
  SCVT endTime = 2.0 * M_PI / ( velocity * kappa );
  int nSteps = 50;
  SCVT timeStep = endTime / nSteps;
  Vector< LO, SC > *pointSolutionTotalTime;
  Vector< LO, SC > *surfaceSolutionTotalTime;
  std::stringstream file;
  file.fill( '0' );
  file.width( 4 );
  for ( int i = 0; i < nSteps; i++ ) {

    pointSolutionTotalTime = new Vector< LO, SC >( pointSolutionTotal );
    surfaceSolutionTotalTime = new Vector< LO, SC >( total );
    pointSolutionTotalTime->scale( std::exp( -iUnit * velocity * kappa * (SCVT) i * timeStep ) );
    surfaceSolutionTotalTime->scale( std::exp( -iUnit * velocity * kappa * (SCVT) i * timeStep ) );
    nodalDataSol[ 0 ] = pointSolutionTotalTime;
    nodeNamesSol[ 0 ] = "total";
    file.str( std::string( ) );
    file.clear( );
    file << "output/domain_hard_scatter_" << i << ".vtu";
    mesh_scatter.printParaviewVtu( file.str( ).c_str( ), 1, nodeNamesSol, nodalDataSol, 0, nullptr, nullptr );
    nodalDataSol[ 0 ] = surfaceSolutionTotalTime;
    file.str( std::string( ) );
    file.clear( );
    file << "output/surface_hard_scatter_" << i << ".vtu";
    mesh.printParaviewVtu( file.str( ).c_str( ), 1, nodeNamesSol, nodalDataSol, 0, nullptr, nullptr );
    delete pointSolutionTotalTime;
    delete surfaceSolutionTotalTime;
  }
}

