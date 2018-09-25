#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../BESpace.h"
#include "../../Vector.h"
#include "../../BEIntegratorHelmholtz.h"
#include "../../BEBilinearFormHelmholtz1Layer.h"
#include "../../BEBilinearFormHelmholtz2Layer.h"
#include "../../IdentityOperator.h"
#include "../../RepresentationFormulaHelmholtz.h"

using namespace std;
using namespace bem4i;

void testHelmholtzSoftScatter(
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
  testHelmholtzSoftScatter( filename, printPrecision );

  return 0;
}

void testHelmholtzSoftScatter( string const &filename, int printPrecision ) {
  typedef complex<double> SC;
  typedef int LO;
  typedef typename SC::value_type SCVT;
  timeval start, stop;
  std::cout.setf( std::ios::showpoint | std::ios::fixed );
  std::cout.precision( printPrecision );

  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  mesh.scale( 0.5 );
  //mesh.move( 1.5, 1.0, 0.0 );
  mesh.refine( 4 );
  //mesh.mapToUnitBall( );

  //SurfaceMesh3D< LO, SC > mesh2;
  //mesh2.load( "input/ball_icosahedron_20.txt" );
  //mesh2.refine( 3 );
  //mesh2.mapToUnitBall();
  //mesh2.printInfo();
  //mesh2.move( -2.0, -1.0, 0.0 );

  //mesh.append( mesh2 );
  mesh.printInfo( );

  SC kappa = 10.0;

  int quadDisj[2] = { 4, 4 };

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
  Vector< LO, SC > sol( nElems );
  Vector< LO, SC > aux( nElems );
  SCVT x[ 3 ], d[ 3 ];
  d[ 0 ] = 0.0; //-1.0; ///std::sqrt(2.0);
  d[ 1 ] = -1.0; // / sqrt(2.0); //1.0 / std::sqrt( 3.0 );
  d[ 2 ] = 0.0; //-1.0 / sqrt(2.0); //1.0 / std::sqrt( 3.0 );

  for ( LO i = 0; i < nNodes; i++ ) {
    mesh.getNode( i, x );
    dir.set( i, -std::exp( iUnit * kappa * DOT3( x, d ) ) );
  }

  std::cout << "Setting up the rhs ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  id.apply( dir, aux, false, -0.5, 0.0 );
  K->apply( dir, neu );
  delete K;
  //neu.scale( -1.0 );
  neu.add( aux );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  BESpace< LO, SC > bespaceV( &mesh, p0, p0 );
  FullMatrix< LO, SC > V( 0, 0 );
  BEBilinearFormHelmholtz1Layer< LO, SC > formV( &bespaceV, nullptr, kappa, SauterSchwab, quadDisj );
  std::cout << "Assembling V ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  formV.assemble( V );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  std::cout << "Solving the system ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  //V.solve( neu );
  LeftIdentityPreconditioner<LO, SC>* M = new LeftIdentityPreconditioner<LO, SC>;
  SCVT prec = 1e-12;
  LO maxIter = 3000;

  V.GMRESSolve( neu, sol, prec, maxIter, 100, nullptr );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  Vector< LO, SC > incident( dir );
  incident.scale( -1.0 );
  Vector< LO, SC > total( nNodes );
  total.setAll( 0.0 );
  string nodeNames[] = { "dirichlet_scatter", "incident", "total" };
  string elemNames[] = { "neumann_scatter" };
  Vector< LO, SC >* nodalData[] = { &dir, &incident, &total };
  Vector< LO, SC >* elemData[] = { &neu };
  mesh.printParaviewVtu( "output/output_soft_scatter.vtu", 3, nodeNames, nodalData, 1, elemNames, elemData );

  //  int aCount = 30;
  //  int bCount = aCount;
  //  int cCount = aCount;
  //  SCVT aSize = 10;
  //  SCVT bSize = aSize;
  //  SCVT cSize = aSize;
  //  SCVT aStep = aSize / (SCVT) aCount;
  //  SCVT bStep = bSize / (SCVT) bCount;
  //  SCVT cStep = cSize / (SCVT) cCount;
  //  int nPoints = ( aCount + 1 ) * ( bCount + 1 ) * ( cCount + 1 );
  //  SCVT *pointCloud = new SCVT[ 3 * nPoints ];
  //  int counter = 0;
  //  for ( int a = 0; a <= aCount; a++ ) {
  //    for ( int b = 0; b <= bCount; b++ ) {
  //      for ( int c = 0; c <= cCount; c++ ) {
  //        pointCloud[ counter++ ] = a * aStep - aSize / 2.0;
  //        pointCloud[ counter++ ] = b * bStep - bSize / 2.0;
  //        pointCloud[ counter++ ] = c * cStep - cSize / 2.0;
  //      }
  //    }
  //  }
  //
  //  Vector< LO, SC > pointSolutionScatter( nPoints );
  //  std::cout << "Evaluating representation formula ... ";
  //  std::cout.flush( );
  //  RepresentationFormulaHelmholtz<LO, SC> formula( &bespaceK, &dir, &neu, kappa, order );
  //  gettimeofday( &start, nullptr );
  //  formula.evaluate( pointCloud, nPoints, false, pointSolutionScatter );
  //  gettimeofday( &stop, nullptr );
  //  Vector< LO, SCVT > pointSolutionScatterReal( nPoints );
  //  Vector< LO, SCVT > pointSolutionScatterImag( nPoints );
  //  Vector< LO, SCVT > pointSolutionIncReal( nPoints );
  //  Vector< LO, SCVT > pointSolutionIncImag( nPoints );
  //  Vector< LO, SCVT > pointSolutionTotalReal( nPoints );
  //  Vector< LO, SCVT > pointSolutionTotalImag( nPoints );
  //  for ( LO i = 0; i < nPoints; i++ ) {
  //    pointSolutionScatterReal.set( i, pointSolutionScatter.get( i ).real( ) );
  //    pointSolutionScatterImag.set( i, pointSolutionScatter.get( i ).imag( ) );
  //    x[ 0 ] = pointCloud[ 3 * i ];
  //    x[ 1 ] = pointCloud[ 3 * i + 1 ];
  //    x[ 2 ] = pointCloud[ 3 * i + 2 ];
  //    pointSolutionIncReal.set( i, std::cos( kappa.real( ) * DOT3( x, d ) ) );
  //    pointSolutionIncImag.set( i, std::sin( kappa.real( ) * DOT3( x, d ) ) );
  //    pointSolutionTotalReal.set( i, pointSolutionScatterReal.get( i ) + pointSolutionIncReal.get( i ) );
  //    pointSolutionTotalImag.set( i, pointSolutionScatterImag.get( i ) + pointSolutionIncImag.get( i ) );
  //  }
  //  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;
  //  string nodeNamesSol[] = { "incident_real", "incident_imag", "scatter_real",
  //    "scatter_imag", "total_real", "total_imag" };
  //  Vector< LO, SCVT >* nodalDataSol[] = { &pointSolutionIncReal, &pointSolutionIncImag,
  //    &pointSolutionScatterReal, &pointSolutionScatterImag,
  //    &pointSolutionTotalReal, &pointSolutionTotalImag };
  //  formula.printParaviewVtu( "output/solution_soft_scatter.vtp", pointCloud, nPoints, 6, nodeNamesSol, nodalDataSol );
  //
  //  delete [] pointCloud;

  string file_scatter = "input/scatter_grid_x.txt";
  SurfaceMesh3D< LO, SCVT > mesh_scvt, scatterer_scvt, scatterer_scvt_t;
  SurfaceMesh3D< LO, SC > mesh_scatterer;
  mesh.copyToReal( mesh_scvt );

  scatterer_scvt.load( file_scatter.c_str( ) );
  scatterer_scvt.scale( 2 );
  scatterer_scvt.refine( 6 );
  scatterer_scvt.printInfo( );
  //scatterer_scvt.scaleAroundCentroid( scaling );
  //scatterer_scvt.refine( scatterRefine4, 2 );
  //scatterer_scvt.refine( scatterRefine9, 3 );
  mesh_scvt.trimEvaluationGrid( scatterer_scvt, scatterer_scvt_t, true );
  scatterer_scvt_t.copyToComplex( mesh_scatterer );

  //mesh_scatterer.printInfo( );
  //return;
  LO nPoints = mesh_scatterer.getNNodes( );
  Vector< LO, SC > pointSolutionScatter( nPoints );
  std::cout << "Evaluating representation formula ... ";
  std::cout.flush( );
  RepresentationFormulaHelmholtz<LO, SC> formula( &bespaceK, &dir, &neu, kappa );
  gettimeofday( &start, nullptr );
  formula.evaluate( mesh_scatterer, false, pointSolutionScatter );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;
  Vector< LO, SC > pointSolutionInc( nPoints );
  Vector< LO, SC > pointSolutionTotal( nPoints );
  for ( LO i = 0; i < nPoints; i++ ) {
    mesh_scatterer.getNode( i, x );
    pointSolutionInc.set( i, std::exp( iUnit * kappa * DOT3( x, d ) ) );
    pointSolutionTotal.set( i, pointSolutionScatter.get( i ) + pointSolutionInc.get( i ) );
  }
  string nodeNamesSol[] = { "incident", "scatter", "total" };
  Vector< LO, SC >* nodalDataSol[] = { &pointSolutionInc, &pointSolutionScatter, &pointSolutionTotal };
  mesh_scatterer.printParaviewVtu( "output/solution_soft_scatter.vtu", 3, nodeNamesSol, nodalDataSol, 0, nullptr, nullptr );

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
    pointSolutionTotalTime->scale( std::exp( -iUnit * velocity * kappa * (SCVT) i * timeStep ) );
    nodalDataSol[ 0 ] = pointSolutionTotalTime;
    nodeNamesSol[ 0 ] = "total";
    file.str( std::string( ) );
    file.clear( );
    file << "output/solution_soft_scatter_2_" << i << ".vtu";
    mesh_scatterer.printParaviewVtu( file.str( ).c_str( ), 1, nodeNamesSol, nodalDataSol, 0, nullptr, nullptr );
    delete pointSolutionTotalTime;
  }

}