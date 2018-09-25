#include "../../Settings.h"

#include <iostream>
#include <sys/time.h>

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

void testHelmholtzNeumannReg(
    string const &filename,
    int refine4,
    int refine9,
    bool mapToUnitBall,
    int quadTypeInt,
    int order,
    int nPoints
    );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  string filename = "input/cube_12.txt";
  testHelmholtzNeumannReg( filename, 4, 0, true, 1, 3, 200 );

  return 0;
}

void testHelmholtzNeumannReg(
    string const & filename,
    int refine4,
    int refine9,
    bool mapToUnitBall,
    int quadTypeInt,
    int order,
    int nPoints
    ) {

  typedef complex<double> SC;
  typedef int LO;
  typedef typename SC::value_type SCVT;

  //SCVT cgprec = 1e-12;
  //LO cgmaxit = 5000;
  SCVT gmresprec = 1e-8;
  LO gmresmaxit = 5000;

  ProgressMonitor::init( "Initializing mesh SC" );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  mesh.refine( refine4, 2 );
  mesh.refine( refine9, 3 );
  if ( mapToUnitBall ) mesh.mapToUnitBall( );
  ProgressMonitor::step( );
  mesh.printInfo( );

  //  ProgressMonitor::init( "Initializing mesh SCVT" );
  //  SurfaceMesh3D< LO, SCVT > mesh_scvt;
  //  mesh_scvt.load( filename.c_str( ) );
  //  mesh_scvt.refine( refine4, 2 );
  //  mesh_scvt.refine( refine9, 3 );
  //  if ( mapToUnitBall ) mesh_scvt.mapToUnitBall( );
  //  ProgressMonitor::step( );
  //  mesh_scvt.printInfo( );

  int quadOrder[ 4 ];
  quadratureType quadType;
  if ( quadTypeInt == 0 ) {
    quadType = Steinbach;
    quadOrder[ 0 ] = quadOrder[ 1 ] = order;
    std::cout << "Using Steinbach quadrature, order " << order << "." <<
        std::endl;
  } else {
    quadType = SauterSchwab;
    quadOrder[ 0 ] = quadOrder[ 1 ] = quadOrder[ 2 ] = quadOrder[ 3 ] = order;
    std::cout << "Using Sauter-Schwab quadrature, order " << order << "." <<
        std::endl;
  }

  int quadDisjoint[] = { 4, 4 };
  //int * quadDisjoint = nullptr;

  SC kappa = 4.0;
  //SC kappa = 4.493;

  //  BESpace< LO, SCVT > bespace00( &mesh_scvt, p0, p0 );
  BESpace< LO, SC > bespace00( &mesh, p0, p0 );
  BESpace< LO, SC > bespace10( &mesh, p1, p0 );
  BESpace< LO, SC > bespace11( &mesh, p1, p1 );

  //  FullMatrix< LO, SCVT > * VlapSCVT = new FullMatrix< LO, SCVT >( 0, 0 );
  //  BEBilinearFormLaplace1Layer< LO, SCVT > formV( &bespace00, quadOrder,
  //      quadType, quadDisjoint );
  //  formV.assemble( *VlapSCVT );
  //
  //  ProgressMonitor::init( "Copying to complex V (dirty, dirty hack)" );
  //  FullMatrix< LO, SC > * Vlap = new FullMatrix< LO, SC >( 0, 0 );
  //  VlapSCVT->copyToComplex( *Vlap );
  //  ProgressMonitor::step( );
  //  delete VlapSCVT;

  FullMatrix< LO, SC > * Vlap = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtz1Layer< LO, SC > formV( &bespace00, quadOrder, 0.0,
      quadType, quadDisjoint );
  formV.assemble( *Vlap );

  FullMatrix< LO, SC > * K = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtz2Layer< LO, SC > formK( &bespace10, quadOrder, kappa,
      quadType, quadDisjoint );
  formK.assemble( *K );

  FullMatrix< LO, SC > * D = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtzHypersingular< LO, SC > formD( &bespace11, quadOrder,
      kappa, quadType, quadDisjoint );
  formD.assemble( *D );

  IdentityOperator< LO, SC > M( &bespace10 );

  SCVT eta = 1.0;
  HelmholtzRegularizedExteriorNeumannOperator< LO, SC > op( &bespace11, Vlap,
      K, D, &M );
  op.setEta( eta );
  //op.setEpsSingleLayer( cgprec );
  //op.setMaxItSingleLayer( cgmaxit );

  LO nElems = mesh.getNElements( );
  LO nNodes = mesh.getNNodes( );

  Vector< LO, SC > rhs( nNodes, true );

  SCVT x1[3], x2[3], x3[3], n[3];
  SCVT * xc;
  int rhsOrder = 5;
  int qsize = quadSizes[ rhsOrder ];
  SCVT * qNodes = new SCVT[ 3 * qsize ];
  SC iUnit( 0.0, 1.0 );
  SCVT xs[ 3 ] = { 0.0, 0.0, 0.9 };
  SCVT area, norm, dot;
  SC val, val1, val2, val3;
  LO elem[ 3 ];
  for ( LO i = 0; i < nElems; ++i ) {
    val1 = val2 = val3 = 0.0;
    mesh.getElement( i, elem );
    mesh.getNodes( i, x1, x2, x3 );
    mesh.getNormal( i, n );
    mesh.getQuadratureNodes( x1, x2, x3, qNodes, rhsOrder );
    area = mesh.getElemArea( i );
    for ( int j = 0; j < qsize; ++j ) {
      xc = qNodes + 3 * j;
      norm = std::sqrt( ( xc[0] - xs[0] )*( xc[0] - xs[0] ) +
          ( xc[1] - xs[1] )*( xc[1] - xs[1] ) +
          ( xc[2] - xs[2] )*( xc[2] - xs[2] ) );
      dot = n[0] * ( xc[0] - xs[0] ) + n[1] * ( xc[1] - xs[1] ) +
          n[2] * ( xc[2] - xs[2] );
      val = (SCVT) quadWeights[ rhsOrder ][ j ] * ( iUnit * kappa * norm -
          (SCVT) 1.0 ) *
          std::exp( iUnit * kappa * norm ) * dot / ( norm * norm * norm );
      val1 += ( (SCVT) 1.0 - (SCVT) quadPoints[ rhsOrder ][ 2 * j ] -
          (SCVT) quadPoints[ rhsOrder ][ 2 * j + 1 ] ) * val;
      val2 += (SCVT) quadPoints[ rhsOrder ][ 2 * j ] * val;
      val3 += (SCVT) quadPoints[ rhsOrder ][ 2 * j + 1 ] * val;
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

  int repOrder = 5;
  ProgressMonitor::init( "Evaluating potentials" );
  SCVT * xe = new SCVT[ 3 * nPoints ];
  for ( LO i = 0; i < nPoints; i++ ) {
    xe[3 * i] = 1.5;
    xe[3 * i + 1] = 0.0;
    xe[3 * i + 2] = 0.0;
  }
  norm = std::sqrt( ( xe[0] - xs[0] )*( xe[0] - xs[0] ) +
      ( xe[1] - xs[1] )*( xe[1] - xs[1] ) +
      ( xe[2] - xs[2] )*( xe[2] - xs[2] ) );
  SC exact = std::exp( iUnit * kappa * norm ) / norm;
  Vector< LO, SC > res_slp( nPoints );
  Vector< LO, SC > res_dlp( nPoints );
  PotentialsHelmholtz< LO, SC > slp( &bespace11, kappa, &w, repOrder );
  PotentialsHelmholtz< LO, SC > dlp( &bespace11, kappa, &v, repOrder );
  slp.singleLayerPotential( xe, nPoints, res_slp );
  dlp.doubleLayerPotential( xe, nPoints, res_dlp );
  ProgressMonitor::step( );

  Vector< LO, SC > res( res_dlp );
  res.scale( -1.0 );
  res.add( res_slp, iUnit * eta );

  std::cout << "Point evaluation in ["
      << xe[ 0 ] << ", "
      << xe[ 1 ] << ", "
      << xe[ 2 ] << "]: "
      << res.get( 0 ) << "." << std::endl;

  if ( nPoints > 1 ) {
    std::cout << "Point evaluation in ["
        << xe[ 0 ] << ", "
        << xe[ 1 ] << ", "
        << xe[ 2 ] << "]: "
        << res.get( 0 ) << "." << std::endl;
  }

  std::cout << "Exact solution: " << exact << std::endl;
  std::cout << "Absolute error: "
      << std::abs( res.get( 0 ) - exact ) << std::endl;

  std::cout << "Relative error: "
      << std::abs( ( res.get( 0 ) - exact ) / exact ) << std::endl;

  delete [] xe;

  string nodeNames[] = { "rhs", "v" };
  string elemNames[] = { "w" };
  Vector< LO, SC >* nodalData[] = { &rhs, &v };
  Vector< LO, SC >* elemData[] = { &w };
  mesh.printParaviewVtu( "output/output.vtu", 2, nodeNames, nodalData, 1,
      elemNames, elemData );
}