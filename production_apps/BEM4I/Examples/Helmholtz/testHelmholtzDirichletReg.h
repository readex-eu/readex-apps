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
#include "../../HelmholtzRegularizedExteriorDirichletOperator.h"
#include "../../BEBilinearFormLaplaceHypersingular.h"
#include "../../IdentityOperator.h"

using namespace std;
using namespace bem4i;

void testHelmholtzDirichletReg(
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
  testHelmholtzDirichletReg( filename, 3, 0, true, 1, 3, 200 );

  return 0;
}

void testHelmholtzDirichletReg(
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

  SCVT cgprec = 1e-12;
  LO cgmaxit = 500;
  SCVT gmresprec = 1e-8;
  LO gmresmaxit = cgmaxit;

  ProgressMonitor::init( "Initializing mesh SC" );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  mesh.refine( refine4, 2 );
  mesh.refine( refine9, 3 );
  if ( mapToUnitBall ) mesh.mapToUnitBall( );
  ProgressMonitor::step( );
  mesh.printInfo( );

  ProgressMonitor::init( "Initializing mesh SCVT" );
  SurfaceMesh3D< LO, SCVT > mesh_scvt;
  mesh_scvt.load( filename.c_str( ) );
  mesh_scvt.refine( refine4, 2 );
  mesh_scvt.refine( refine9, 3 );
  if ( mapToUnitBall ) mesh_scvt.mapToUnitBall( );
  ProgressMonitor::step( );
  mesh_scvt.printInfo( );

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

  BESpace< LO, SC > bespace00( &mesh, p0, p0 );
  BESpace< LO, SC > bespace10( &mesh, p1, p0 );
  BESpace< LO, SCVT > bespace11( &mesh_scvt, p1, p1 );

  FullMatrix< LO, SCVT > * DlapSCVT = new FullMatrix< LO, SCVT >( 0, 0 );
  BEBilinearFormLaplaceHypersingular< LO, SCVT > formD( &bespace11, quadOrder,
      quadType, quadDisjoint );
  formD.assemble( *DlapSCVT );

  LO nElems = mesh.getNElements( );
  LO nNodes = mesh.getNNodes( );

  ProgressMonitor::init( "Stabilizing Dlap" );
  Vector<LO, SCVT> ar( nNodes, true );
  LO idx[ 3 ];
  SCVT areaW;
  for ( LO i = 0; i < nElems; ++i ) {
    mesh.getElement( i, idx );
    areaW = mesh.getElemArea( i ) / 3.0;
    ar.add( idx[ 0 ], areaW );
    ar.add( idx[ 1 ], areaW );
    ar.add( idx[ 2 ], areaW );
  }

  for ( LO j = 0; j < ar.getLength( ); ++j ) {
    for ( LO i = 0; i < ar.getLength( ); ++i ) {
      DlapSCVT->add( i, j, ar.get( i ) * ar.get( j ) );
    }
  }
  ProgressMonitor::step( );

  ProgressMonitor::init( "Copying to complex D (dirty, dirty hack)" );
  FullMatrix< LO, SC > * Dlap = new FullMatrix< LO, SC >( 0, 0 );
  DlapSCVT->copyToComplex( *Dlap );
  ProgressMonitor::step( );
  delete DlapSCVT;

  FullMatrix< LO, SC > * V = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtz1Layer< LO, SC > formV( &bespace00, quadOrder,
      kappa, quadType, quadDisjoint );
  formV.assemble( *V );

  FullMatrix< LO, SC > * K = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtz2Layer< LO, SC > formK( &bespace10, quadOrder, kappa,
      quadType, quadDisjoint );
  formK.assemble( *K );

  IdentityOperator< LO, SC > M( &bespace10 );

  SCVT eta = 1.0;
  HelmholtzRegularizedExteriorDirichletOperator< LO, SC > op( &bespace00, V,
      K, Dlap, &M );
  op.setEta( eta );
  op.setEpsHypersingular( cgprec );
  op.setMaxItHypersingular( cgmaxit );

  Vector< LO, SC > rhs( nElems );
  SCVT x1[3], x2[3], x3[3];
  SCVT * xc;
  int rhsOrder = 5;
  int qsize = quadSizes[ rhsOrder ];
  SCVT * qNodes = new SCVT[ 3 * qsize ];
  SC iUnit( 0.0, 1.0 );
  SCVT xs[ 3 ] = { 0.0, 0.0, 0.9 };
  SCVT area, norm;
  SC val;
  for ( LO i = 0; i < nElems; ++i ) {
    val = 0.0;
    mesh.getNodes( i, x1, x2, x3 );
    mesh.getQuadratureNodes( x1, x2, x3, qNodes, rhsOrder );
    area = mesh.getElemArea( i );
    for ( int j = 0; j < qsize; ++j ) {
      xc = qNodes + 3 * j;
      norm = std::sqrt( ( xc[0] - xs[0] )*( xc[0] - xs[0] ) +
          ( xc[1] - xs[1] )*( xc[1] - xs[1] ) +
          ( xc[2] - xs[2] )*( xc[2] - xs[2] ) );
      val += (SCVT) quadWeights[ rhsOrder ][ j ] *
          std::exp( iUnit * kappa * norm ) / norm;
    }
    val *= area;
    rhs.set( i, val );
  }
  delete [] qNodes;

  Vector< LO, SC > w( nElems );
  ProgressMonitor::init( "Solving the system" );
  op.GMRESSolve( rhs, w, gmresprec, gmresmaxit, gmresmaxit );
  ProgressMonitor::step( );

  Vector< LO, SC > v;
  op.getV( w, v );

  delete V;
  delete K;
  delete Dlap;

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
  PotentialsHelmholtz< LO, SC > slp( &bespace00, kappa, &w, repOrder );
  PotentialsHelmholtz< LO, SC > dlp( &bespace00, kappa, &v, repOrder );
  slp.singleLayerPotential( xe, nPoints, res_slp );
  dlp.doubleLayerPotential( xe, nPoints, res_dlp );
  ProgressMonitor::step( );

  Vector< LO, SC > res( res_slp );
  res.add( res_dlp, iUnit * eta );

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

  string nodeNames[] = { "v" };
  string elemNames[] = { "rhs", "w" };
  Vector< LO, SC >* nodalData[] = { &v };
  Vector< LO, SC >* elemData[] = { &rhs, &w };
  mesh.printParaviewVtu( "output/output.vtu", 1, nodeNames, nodalData, 2,
      elemNames, elemData );

}