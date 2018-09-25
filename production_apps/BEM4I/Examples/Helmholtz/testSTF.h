#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../ProgressMonitor.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../BESpace.h"
#include "../../Vector.h"
#include "../../BEIntegratorHelmholtz.h"
#include "../../BEBilinearFormHelmholtzHypersingular.h"
#include "../../BEBilinearFormHelmholtz1Layer.h"
#include "../../BEBilinearFormHelmholtz2Layer.h"
#include "../../IdentityOperator.h"
#include "../../RepresentationFormulaHelmholtz.h"
#include "../../STFOperator.h"
#include "../../STFPreconditioner.h"

using namespace std;
using namespace bem4i;

void testSTF(
    string const &filename,
    int refine4,
    int refine9,
    bool mapToUnitBall
    );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  string filename = "input/cube_12.txt";
  testSTF( filename, 4, 0, false );

  return 0;
}

void testSTF(
    const string& filename,
    int refine4,
    int refine9,
    bool mapToUnitBall
    ) {
  //MPI_Init(nullptr, nullptr);

  typedef complex<double> SC;
  typedef int LO;
  typedef typename SC::value_type SCVT;

  SCVT gmresprec = 1e-8;
  LO gmresmaxit = 3000;

  ProgressMonitor::init( "Initializing mesh" );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );

  //  SurfaceMesh3D< LO, SC > mesh2;
  //  mesh2.load( filename.c_str( ) );
  //  mesh2.move( 0.0, 0.0, -4.0 );
  //  mesh2.scale( 1.0, 1.0, 0.5 );
  //  mesh.append( mesh2 );
  //  mesh.move( 0.0, 0.0, 1.0 );
  mesh.refine( refine4, 2 );
  mesh.refine( refine9, 3 );
  if ( mapToUnitBall ) mesh.mapToUnitBall( );
  ProgressMonitor::step( );
  mesh.printInfo( );

  quadratureType quadType = SauterSchwab;
  int orderNear = 3;
  int orderFar = 4;
  int quadOrder[] = { orderNear, orderNear, orderNear, orderNear };
  int quadDisjoint[] = { orderFar, orderFar };

  SC iUnit( 0.0, 1.0 );
  SC etai = 1.0;
  SC etae( 5.0, -0.1 );
  SCVT lambda = 0.8;
  SCVT kappa = 2.0 * M_PI / lambda;
  SC refrIndex( std::sqrt( ( std::abs( etae ) + etae.real( ) ) / 2.0 ),
      std::sqrt( ( std::abs( etae ) - etae.real( ) ) / 2.0 ) );
  SC ki = kappa * std::sqrt( etai );
  SC ke = kappa * refrIndex;

  ki = 2.0;
  ke = SC( 4.0, 0.1 );

  std::cout << "kappa_i = " << ki << std::endl;
  std::cout << "kappa_e = " << ke << std::endl;

  //  Tree<BECluster<LO, SC>*, LO> tree;
  //  mesh.nestedDissection( tree, 20 );
  //  FastBESpace< LO, SC > bespace11( &mesh, p1, p1, &tree, 1.2, 8, 5, 0 );
  //  bespace11.setScaleACA( 1e-12 );
  //  bespace11.setEpsilonACA( 1e-2 );

  BESpace< LO, SC > bespace11( &mesh, p1, p1 );

  FullMatrix< LO, SC > * Vi = new FullMatrix< LO, SC >( 0, 0 );
  //ACAMatrix< LO, SC > * Vi = new ACAMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtz1Layer< LO, SC > formVi( &bespace11, quadOrder,
      ki, quadType, quadDisjoint );
  formVi.assemble( *Vi );

  FullMatrix< LO, SC > * Ve = new FullMatrix< LO, SC >( 0, 0 );
  //ACAMatrix< LO, SC > * Ve = new ACAMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtz1Layer< LO, SC > formVe( &bespace11, quadOrder,
      ke, quadType, quadDisjoint );
  formVe.assemble( *Ve );

  FullMatrix< LO, SC > * Ki = new FullMatrix< LO, SC >( 0, 0 );
  //ACAMatrix< LO, SC > * Ki = new ACAMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtz2Layer< LO, SC > formKi( &bespace11, quadOrder,
      ki, quadType, quadDisjoint );
  formKi.assemble( *Ki );

  FullMatrix< LO, SC > * Ke = new FullMatrix< LO, SC >( 0, 0 );
  //ACAMatrix< LO, SC > * Ke = new ACAMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtz2Layer< LO, SC > formKe( &bespace11, quadOrder,
      ke, quadType, quadDisjoint );
  formKe.assemble( *Ke );

  FullMatrix< LO, SC > * Di = new FullMatrix< LO, SC >( 0, 0 );
  //ACAMatrix< LO, SC > * Di = new ACAMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtzHypersingular< LO, SC > formDi( &bespace11, quadOrder,
      ki, quadType, quadDisjoint );
  formDi.assemble( *Di );

  FullMatrix< LO, SC > * De = new FullMatrix< LO, SC >( 0, 0 );
  //ACAMatrix< LO, SC > * De = new ACAMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtzHypersingular< LO, SC > formDe( &bespace11, quadOrder,
      ke, quadType, quadDisjoint );
  formDe.assemble( *De );

  SparseMatrix<LO, SC> M;
  IdentityOperator< LO, SC > id( &bespace11 );
  id.assemble( M );

  STFOperator< LO, SC > stf( &mesh, ki, ke, Vi, Ve, Ki, Ke, Di, De, &M );

  LO nNodes = mesh.getNNodes( );

  Vector< LO, SC > f( nNodes );
  Vector< LO, SC > g( nNodes );

  SCVT y[ 3 ] = { 0, 0, 0.3 };
  SCVT x[ 3 ], n[ 3 ];
  SCVT norm, dot;

  for ( LO i = 0; i < nNodes; ++i ) {
    mesh.getNode( i, x );
    norm = DIST3( x, y );
    f.set( i, ( SCVT ) - PI_FACT * std::exp( iUnit * ki * norm ) / norm );
  }

  for ( LO i = 0; i < nNodes; ++i ) {
    mesh.getNormalNodal( i, n );
    mesh.getNode( i, x );
    norm = DIST3( x, y );
    dot = n[ 0 ] * ( x[ 0 ] - y[ 0 ] ) + n[ 1 ] * ( x[ 1 ] - y[ 1 ] ) +
        n[ 2 ] * ( x[ 2 ] - y[ 2 ] );
    g.set( i, -( iUnit * ki * norm - (SCVT) 1.0 ) *
        std::exp( iUnit * ki * norm ) * dot / ( norm * norm * norm ) );
  }

  Vector< LO, SC > rhs;
  stf.getRHS( g, f, rhs );

  //rhs.print();

  ProgressMonitor::init( "Solving the system" );
  Vector< LO, SC > res( nNodes + nNodes );
  STFPreconditioner< LO, SC > stfPrec( &mesh, Vi, Ve, Ki, Ke, Di, De, &M );
  stf.FGMRESSolve( rhs, res, gmresprec, gmresmaxit, gmresmaxit, &stfPrec );
  //stf.FGMRESSolve( rhs, res, gmresprec, gmresmaxit, gmresmaxit, nullptr );
  //Vector< LO, SC > res( rhs );
  //stf.LUSolve( res );
  ProgressMonitor::step( );

  //res.print(); return;

  Vector< LO, SC > ui, ue, ti, te;
  stf.getTraces( res, f, g, ui, ue, ti, te );
  ui.add( f, -1.0 );
  ti.add( g, -1.0 );

  delete Vi;
  delete Ve;
  delete Ki;
  delete Ke;
  delete Di;
  delete De;

  string file_scatter = "input/scatter_grid_y.txt";
  SurfaceMesh3D< LO, SC > ms;
  SurfaceMesh3D< LO, SC > mi;
  SurfaceMesh3D< LO, SC > me;
  ms.load( file_scatter.c_str( ) );
  ms.scale( 2.5 );
  ms.refine( 5 );
  ms.printInfo( );

  mesh.trimEvaluationGrid( ms, mi, me );

  ProgressMonitor::init( "Evaluating inside" );
  LO npi = mi.getNNodes( );
  Vector< LO, SC > resi( npi );
  RepresentationFormulaHelmholtz<LO, SC> fi( &bespace11, &ui, &ti, ki );
  fi.evaluate( mi, true, resi );
  ProgressMonitor::step( );

  ProgressMonitor::init( "Evaluating outside" );
  LO npe = me.getNNodes( );
  Vector< LO, SC > rese( npe );
  RepresentationFormulaHelmholtz<LO, SC> fe( &bespace11, &ue, &te, ke );
  fe.evaluate( me, false, rese );
  ProgressMonitor::step( );



  string nodeNames[] = { "ui", "ue", "ti", "te" };
  //string elemNames[] = { "ti", "te" };
  Vector< LO, SC >* nodalData[] = { &ui, &ue, &ti, &te };
  //Vector< LO, SC >* elemData[] = { &ti, &te };
  mesh.printParaviewVtu( "output/output.vtu", 4, nodeNames, nodalData, 0,
      nullptr, nullptr );

  string namesi[] = { "ui" };
  Vector< LO, SC >* datai[] = { &resi };
  mi.printParaviewVtu( "output/output_i.vtu", 1, namesi, datai, 0, nullptr,
      nullptr );

  string namese[] = { "ue" };
  Vector< LO, SC >* datae[] = { &rese };
  me.printParaviewVtu( "output/output_e.vtu", 1, namese, datae, 0, nullptr,
      nullptr );

  //MPI_Finalize();
}

//void testSTF(
//    const string& filename,
//    int refine4,
//    int refine9,
//    bool mapToUnitBall
//    ) {
//
//  typedef complex<double> SC;
//  typedef int LO;
//  typedef typename SC::value_type SCVT;
//
//  SCVT gmresprec = 0.02;
//  LO gmresmaxit = 100;
//
//  ProgressMonitor::init( "Initializing mesh" );
//  SurfaceMesh3D< LO, SC > mesh;
//  SurfaceMesh3D< LO, SC > mesh2;
//  mesh.load( filename.c_str( ) );
//  mesh.scale(2.0);
//  mesh2.load( filename.c_str( ) );
//  mesh2.move( 0.0, 0.0, -3.8 );
//  mesh2.scale( 2.0, 2.0, 1.0 );
//  mesh.append( mesh2 );
//  mesh.move( 0.0, 0.0, 1.0 );
//  mesh.refine( refine4, 2 );
//  mesh.refine( refine9, 3 );
//  mesh.printParaviewVtu( "output/output.vtu");
//  if ( mapToUnitBall ) mesh.mapToUnitBall( );
//  ProgressMonitor::step( );
//  mesh.printInfo( );
//
//  quadratureType quadType = SauterSchwab;
//  int orderNear = 3;
//  int orderFar = 3;
//  int quadOrder[] = { orderNear, orderNear, orderNear, orderNear };
//  int quadDisjoint[] = { orderFar, orderFar };
//
//  SC iUnit( 0.0, 1.0 );
//  SC etai = 1.0;
//  SC etae( 5.0, -0.5 );
//  SCVT lambda = 0.25;
//  SCVT kappa = 2.0 * M_PI / lambda;
//  SC refrIndex( std::sqrt( ( std::abs( etae ) + etae.real( ) ) / 2.0 ),
//      std::sqrt( ( std::abs( etae ) - etae.real( ) ) / 2.0 ) );
//  SC ki = kappa * std::sqrt( etai );
//  SC ke = kappa * refrIndex;
//
//  //ki = 0.5;
//  //ke = 0.7;
//
//  std::cout << "kappa_i = " << ki << std::endl;
//  std::cout << "kappa_e = " << ke << std::endl;
//
//  Tree<BECluster<LO, SC>*, LO> tree, tree2, tree3;
//  mesh.nestedDissection( tree, 60 );
//  mesh.nestedDissection( tree2, 60 );
//  mesh.nestedDissection( tree3, 60 );
//
//  FastBESpace< LO, SC > bespaceK( &mesh, p1, p0, &tree2, 1.2, 8, 5, 0 );
//
//  FastBESpace< LO, SC > bespace00( &mesh, p0, p0, &tree, 1.2, 8, 5, 0 );
//  FastBESpace< LO, SC > bespace10( &mesh, p1, p0, &tree2, 1.2, 8, 5, 0 );
//  FastBESpace< LO, SC > bespace11( &mesh, p1, p1, &tree3, 1.2, 8, 5, 0 );
//
//  bespace00.setScaleACA( 1e-6 );
//  bespace10.setScaleACA( 1e-6 );
//  bespace11.setScaleACA( 1e-12 );
//  bespace00.setEpsilonACA( 1e-2 );
//  bespace10.setEpsilonACA( 1e-2 );
//  bespace11.setEpsilonACA( 1e-2 );
//
//  //FullMatrix< LO, SC > * Vi = new FullMatrix< LO, SC >( 0, 0 );
//  ACAMatrix< LO, SC > * Vi = new ACAMatrix< LO, SC >( 0, 0 );
//  BEBilinearFormHelmholtz1Layer< LO, SC > formVi( &bespace00, quadOrder,
//      ki, quadType, quadDisjoint );
//  formVi.assemble( *Vi );
//
//  //FullMatrix< LO, SC > * Ve = new FullMatrix< LO, SC >( 0, 0 );
//  ACAMatrix< LO, SC > * Ve = new ACAMatrix< LO, SC >( 0, 0 );
//  BEBilinearFormHelmholtz1Layer< LO, SC > formVe( &bespace00, quadOrder,
//      ke, quadType, quadDisjoint );
//  formVe.assemble( *Ve );
//
//  //FullMatrix< LO, SC > * Ki = new FullMatrix< LO, SC >( 0, 0 );
//  ACAMatrix< LO, SC > * Ki = new ACAMatrix< LO, SC >( 0, 0 );
//  BEBilinearFormHelmholtz2Layer< LO, SC > formKi( &bespace10, quadOrder,
//      ki, quadType, quadDisjoint );
//  formKi.assemble( *Ki );
//
//  //FullMatrix< LO, SC > * Ke = new FullMatrix< LO, SC >( 0, 0 );
//  ACAMatrix< LO, SC > * Ke = new ACAMatrix< LO, SC >( 0, 0 );
//  BEBilinearFormHelmholtz2Layer< LO, SC > formKe( &bespace10, quadOrder,
//      ke, quadType, quadDisjoint );
//  formKe.assemble( *Ke );
//
//  //FullMatrix< LO, SC > * Di = new FullMatrix< LO, SC >( 0, 0 );
//  ACAMatrix< LO, SC > * Di = new ACAMatrix< LO, SC >( 0, 0 );
//  BEBilinearFormHelmholtzHypersingular< LO, SC > formDi( &bespace11, quadOrder,
//      ki, quadType, quadDisjoint );
//  formDi.assemble( *Di );
//
//  //FullMatrix< LO, SC > * De = new FullMatrix< LO, SC >( 0, 0 );
//  ACAMatrix< LO, SC > * De = new ACAMatrix< LO, SC >( 0, 0 );
//  BEBilinearFormHelmholtzHypersingular< LO, SC > formDe( &bespace11, quadOrder,
//      ke, quadType, quadDisjoint );
//  formDe.assemble( *De );
//
//  IdentityOperator< LO, SC > M( &bespace10 );
//
//  STFOperator< LO, SC > stf( &mesh, ki, ke, Vi, Ve, Ki, Ke, Di, De, &M );
//
//  LO nNodes = mesh.getNNodes( );
//  LO nElems = mesh.getNElements( );
//
//  Vector< LO, SC > f( nNodes );
//  Vector< LO, SC > g( nElems );
//
//  SCVT y[ 3 ] = { 0, 0, 2.0 };
//  SCVT x[ 3 ], n[ 3 ];
//  SCVT norm, dot;
//
//  for ( LO i = 0; i < nNodes; ++i ) {
//    mesh.getNode( i, x );
//    norm = DIST3( x, y );
//    f.set( i, -PI_FACT * std::exp( iUnit * ki * norm ) / norm );
//  }
//
//  for ( LO i = 0; i < nElems; ++i ) {
//    mesh.getNormal( i, n );
//    mesh.getCentroid( i, x );
//    norm = DIST3( x, y );
//    dot = n[ 0 ] * ( x[ 0 ] - y[ 0 ] ) + n[ 1 ] * ( x[ 1 ] - y[ 1 ] ) +
//        n[ 2 ] * ( x[ 2 ] - y[ 2 ] );
//    g.set( i, -( iUnit * ki * norm - 1.0 ) *
//        std::exp( iUnit * ki * norm ) * dot / ( norm * norm * norm ) );
//  }
//
//  Vector< LO, SC > rhs;
//  stf.getRHS( g, f, rhs );
//
//  ProgressMonitor::init( "Solving the system" );
//  Vector< LO, SC > res( nNodes + nNodes );
//  STFPreconditioner< LO, SC > stfPrec( &mesh, Vi, Ve, Ki, Ke, Di, De );
//  stf.GMRESSolve( rhs, res, gmresprec, gmresmaxit, gmresmaxit, &stfPrec );
//  //Vector< LO, SC > res( rhs );
//  //stf.LUSolve( res );
//  ProgressMonitor::step( );
//
//  //res.print(); return;
//
//  Vector< LO, SC > ui, ue, ti, te;
//  stf.getTraces( res, f, g, ui, ue, ti, te );
//
//  delete Vi;
//  delete Ve;
//  delete Ki;
//  delete Ke;
//  delete Di;
//  delete De;
//
//  string file_scatter = "input/scatter_grid_y.txt";
//  SurfaceMesh3D< LO, SC > ms;
//  SurfaceMesh3D< LO, SC > mi;
//  SurfaceMesh3D< LO, SC > me;
//  ms.load( file_scatter.c_str( ) );
//  ms.scale( 5.0 );
//  ms.refine( 8 );
//  ms.printInfo( );
//  mesh.trimEvaluationGrid( ms, mi, false );
//  mesh.trimEvaluationGrid( ms, me, true );
//
//  ProgressMonitor::init( "Evaluating inside" );
//  LO npi = mi.getNNodes( );
//  Vector< LO, SC > resi( npi );
//  RepresentationFormulaHelmholtz<LO, SC> fi( &bespace10, &ui, &ti, ki );
//  fi.evaluate( mi, true, resi );
//  ProgressMonitor::step( );
//
//  ProgressMonitor::init( "Evaluating outside" );
//  LO npe = me.getNNodes( );
//  Vector< LO, SC > rese( npe );
//  RepresentationFormulaHelmholtz<LO, SC> fe( &bespace10, &ue, &te, ke );
//  fe.evaluate( me, false, rese );
//  ProgressMonitor::step( );
//
//  ui.add( f, -1.0 );
//  ti.add( g, -1.0 );
//
//  string nodeNames[] = { "ui", "ue" };
//  string elemNames[] = { "ti", "te" };
//  Vector< LO, SC >* nodalData[] = { &ui, &ue };
//  Vector< LO, SC >* elemData[] = { &ti, &te };
//  mesh.printParaviewVtu( "output/output.vtu", 2, nodeNames, nodalData, 2,
//      elemNames, elemData );
//
//  string namesi[] = { "ui" };
//  Vector< LO, SC >* datai[] = { &resi };
//  mi.printParaviewVtu( "output/output_i.vtu", 1, namesi, datai, 0, nullptr,
//      nullptr );
//
//  string namese[] = { "ue" };
//  Vector< LO, SC >* datae[] = { &rese };
//  me.printParaviewVtu( "output/output_e.vtu", 1, namese, datae, 0, nullptr,
//      nullptr );
//}