#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../ProgressMonitor.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../BESpace.h"
#include "../../Vector.h"
#include "../../BEIntegratorHelmholtz.h"
#include "../../BEBilinearFormHelmholtz1Layer.h"
#include "../../BEBilinearFormHelmholtz2Layer.h"
#include "../../BEBilinearFormHelmholtzHypersingular.h"
#include "../../HelmholtzRegularizedExteriorNeumannOperator.h"
#include "../../IdentityOperator.h"

using namespace std;
using namespace bem4i;

void testHelmholtzHardScatterReg( );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  testHelmholtzHardScatterReg( );

  return 0;
}

void testHelmholtzHardScatterReg( ) {

  typedef complex<double> SC;
  typedef int LO;
  typedef typename SC::value_type SCVT;

  std::string filename = "input/icosahedron.txt";
  SCVT scale = 1.0;
  int refine4 = 3;
  int refine9 = 0;
  bool mapToUnitBall = true;

  SCVT eta = 1.0;
  SCVT epsilonACA = 1e-5;
  LO maxElems = 30;

  int orderNear = 3;
  int orderFar = 4;
  int orderRhs = 4;
  int orderRep = 4;

  SCVT d[ 3 ] = { 0.0, 0.0, 1.0 };

  string file_scatter = "input/scatter_grid_y.txt";
  SCVT scatterScale = 2.0;
  int scatterRefine4 = 4;
  int scatterRefine9 = 0;

  SCVT lambda = 0.7;
  //  SCVT velocity = 350.0;
  SCVT kappa = 2 * M_PI / lambda;
  //  int nSteps = 50;

  SCVT scaling = kappa;
  SCVT kappaScaled = kappa / scaling;

  //SCVT cgprec = 1e-12;
  //LO cgmaxit = 500;
  SCVT gmresprec = 1e-8;
  LO gmresmaxit = 5000;

  ProgressMonitor::init( "Initializing mesh SC" );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  mesh.refine( refine4, 2 );
  mesh.refine( refine9, 3 );
  if ( mapToUnitBall ) mesh.mapToUnitBall( );
  mesh.scaleAroundCentroid( scale );
  ProgressMonitor::step( );
  mesh.printInfo( );

  mesh.scaleAroundCentroid( scaling );
  mesh.printInfo( );

  int orderNearA [] = { orderNear, orderNear, orderNear, orderNear };
  int orderFarA [] = { orderFar, orderFar };
  std::cout << "Sauter-Schwab near-field order " << orderNear <<
      ", Gauss far-field order " << orderFar << std::endl;

  Tree< BECluster< LO, SC > *, LO > tree1, tree2, tree3;
  mesh.nestedDissection( tree1, maxElems );
  mesh.nestedDissection( tree2, maxElems );
  mesh.nestedDissection( tree3, maxElems );

  FastBESpace< LO, SC > bespace00( &mesh, p0, p0, &tree1, eta );
  FastBESpace< LO, SC > bespace10( &mesh, p1, p0, &tree2, eta );
  FastBESpace< LO, SC > bespace11( &mesh, p1, p1, &tree3, eta );
  bespace00.setEpsilonACA( epsilonACA );
  bespace10.setEpsilonACA( epsilonACA );
  bespace11.setEpsilonACA( epsilonACA );

  ACAMatrix< LO, SC > * Vlap = new ACAMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtz1Layer< LO, SC > formV( &bespace00, orderNearA, 0.0,
      SauterSchwab, orderFarA );
  formV.assemble( *Vlap );

  ACAMatrix< LO, SC > * K = new ACAMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtz2Layer< LO, SC > formK( &bespace10, orderNearA,
      kappaScaled, SauterSchwab, orderFarA );
  formK.assemble( *K );

  ACAMatrix< LO, SC > * D = new ACAMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtzHypersingular< LO, SC > formD( &bespace11, orderNearA,
      kappaScaled, SauterSchwab, orderFarA );
  formD.assemble( *D );

  IdentityOperator< LO, SC > M( &bespace10 );

  SCVT etaCoupling = 1.0;
  HelmholtzRegularizedExteriorNeumannOperator< LO, SC > op( &bespace11, Vlap,
      K, D, &M );
  op.setEta( etaCoupling );
  //op.setEpsSingleLayer( cgprec );
  //op.setMaxItSingleLayer( cgmaxit );

  LO nElems = mesh.getNElements( );
  LO nNodes = mesh.getNNodes( );

  Vector< LO, SC > rhs( nNodes, true );

  SCVT x1[3], x2[3], x3[3], n[3];
  SCVT * xc;
  int qsize = quadSizes[ orderRhs ];
  SCVT * qNodes = new SCVT[ 3 * qsize ];
  SC iUnit( 0.0, 1.0 );
  SCVT area;
  SC val, val1, val2, val3;
  LO elem[ 3 ];
  for ( LO i = 0; i < nElems; ++i ) {
    val1 = val2 = val3 = 0.0;
    mesh.getElement( i, elem );
    mesh.getNodes( i, x1, x2, x3 );
    mesh.getNormal( i, n );
    mesh.getQuadratureNodes( x1, x2, x3, qNodes, orderRhs );
    area = mesh.getElemArea( i );
    for ( int j = 0; j < qsize; ++j ) {
      xc = qNodes + 3 * j;
      val = (SCVT) quadWeights[ orderRhs ][ j ] * ( -iUnit * kappaScaled *
          DOT3( xc, n ) * std::exp( iUnit * kappaScaled * DOT3( xc, d ) ) ) /
          scaling;
      val1 += ( (SCVT) 1.0 - (SCVT) quadPoints[ orderRhs ][ 2 * j ] -
          (SCVT) quadPoints[ orderRhs ][ 2 * j + 1 ] ) * val;
      val2 += (SCVT) quadPoints[ orderRhs ][ 2 * j ] * val;
      val3 += (SCVT) quadPoints[ orderRhs ][ 2 * j + 1 ] * val;
    }
    val1 *= area;
    val2 *= area;
    val3 *= area;
    rhs.add( elem[ 0 ], val1 );
    rhs.add( elem[ 1 ], val2 );
    rhs.add( elem[ 2 ], val3 );
  }
  delete [] qNodes;

  /*
  Vector< LO, SC > v( nNodes );
  ProgressMonitor::init( "Solving the system" );
  op.GMRESSolve( rhs, v, gmresprec, gmresmaxit, gmresmaxit );
  ProgressMonitor::step( );
  Vector< LO, SC > w;
  op.getW( v, w );
   */
  ///*
  Vector< LO, SC > wv( nElems + nNodes );
  Vector< LO, SC > rhsex( nElems + nNodes );
  for ( LO i = 0; i < nNodes; ++i ) {
    rhsex.set( i + nElems, rhs.get( i ) );
  }
  ProgressMonitor::init( "Solving the system" );
  op.GMRESSolve( rhsex, wv, gmresprec, gmresmaxit, gmresmaxit );
  ProgressMonitor::step( );
  Vector< LO, SC > v( nNodes, false );
  Vector< LO, SC > w( nElems, false );
  for ( LO i = 0; i < nElems; ++i ) {
    w.set( i, wv.get( i ) );
  }
  for ( LO i = 0; i < nNodes; ++i ) {
    v.set( i, wv.get( i + nElems ) );
  }
  //*/

  delete Vlap;
  delete K;
  delete D;

  SurfaceMesh3D< LO, SCVT > mesh_scvt, scatterer_scvt, scatterer_scvt_t;
  SurfaceMesh3D< LO, SC > scatterer;
  mesh.copyToReal( mesh_scvt );
  scatterer_scvt.load( file_scatter.c_str( ) );
  scatterer_scvt.scaleAroundCentroid( scatterScale );
  scatterer_scvt.scaleAroundCentroid( scaling );
  scatterer_scvt.refine( scatterRefine4, 2 );
  scatterer_scvt.refine( scatterRefine9, 3 );
  mesh_scvt.trimEvaluationGrid( scatterer_scvt, scatterer_scvt_t, true );
  scatterer_scvt_t.copyToComplex( scatterer );
  scatterer.printInfo( );
  LO nPoints = scatterer.getNNodes( );

  ProgressMonitor::init( "Evaluating potentials" );
  Vector< LO, SC > res_slp( nPoints );
  Vector< LO, SC > res_dlp( nPoints );
  PotentialsHelmholtz< LO, SC > slp( &bespace11, kappaScaled, &w, orderRep );
  PotentialsHelmholtz< LO, SC > dlp( &bespace11, kappaScaled, &v, orderRep );
  slp.singleLayerPotential( scatterer, res_slp );
  dlp.doubleLayerPotential( scatterer, res_dlp );
  Vector< LO, SC > pointSolutionScatter( res_dlp );
  pointSolutionScatter.scale( -1.0 );
  pointSolutionScatter.add( res_slp, iUnit * etaCoupling );
  ProgressMonitor::step( );

  SCVT x[ 3 ];
  Vector< LO, SC > pointSolutionInc( nPoints );
  Vector< LO, SC > pointSolutionTotal( nPoints );
  for ( LO i = 0; i < nPoints; ++i ) {
    scatterer.getNode( i, x );
    pointSolutionInc.set( i, std::exp( iUnit * kappaScaled * DOT3( x, d ) ) );
    pointSolutionTotal.set( i,
        pointSolutionScatter.get( i ) + pointSolutionInc.get( i ) );
  }

  // careful, bespaces hold pointers!!!
  mesh.scaleAroundCentroid( 1.0 / scaling );
  scatterer.scaleAroundCentroid( 1.0 / scaling );

  mesh.printParaviewVtu( "output/scatter/surface_hard_scatter.vtu" );
  string nodeNamesSol[] = { "incident", "scatter", "total" };
  Vector< LO, SC > * nodalDataSol[] = { &pointSolutionInc,
    &pointSolutionScatter, &pointSolutionTotal };
  scatterer.printParaviewVtu( "output/scatter/domain_hard_scatter.vtu", 3,
      nodeNamesSol, nodalDataSol, 0, nullptr, nullptr );
  /*
  SCVT endTime = 2.0 * M_PI / ( velocity * kappa.real( ) );
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
   */
}