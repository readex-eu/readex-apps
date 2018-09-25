#include "../../Settings.h"

#include <iostream>
#include <numeric>

#include "../auxiliary.h"
#include "../../ProgressMonitor.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../BESpace.h"
#include "../../Vector.h"
#include "../../IdentityOperator.h"
#include "../../BEBilinearFormLame1Layer.h"
#include "../../BEBilinearFormLaplace2Layer.h"
#include "../../BEBilinearFormLame2Layer.h"
#include "../../BEBilinearFormLameHypersingular.h"
#include "../../CompoundLinearOperator.h"
#include "../../Lame1LayerP0P0MultilvlPrecond.h"
#include "../../BlockLinearOperator.h"
#include "../../SumLinearOperator.h"
#include "../../InverseLinearOperator.h"

using namespace std;
using namespace bem4i;

void testLameMixed(
    string const & filename,
    const string & bc,
    int refine4,
    int refine9,
    bool mapToUnitBall,
    int quadType,
    int orderNear,
    int orderFar
    );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  string filename = "input/spring_17k.txt";
  string bc = "input/spring_17k_bc.txt";
  testLameMixed( filename, bc, 2, 0, false, 1, 4, 4 );

  return 0;
}

void testLameMixed(
    const string & filename,
    const string & bc,
    int refine4,
    int refine9,
    bool mapToUnitBall,
    int quadTypeInt,
    int orderNear,
    int orderFar
    ) {

  typedef double SC;
  typedef int LO;

  //std::cout.precision( 2 );
  //std::cout.setf( std::ios::scientific );

  SC nu = 0.33;
  SC E = 1.1e5;
  SC CGprec = 1e-8;
  LO CGmaxit = 1000;
  SC CGprecInner = 1e-8;
  LO CGmaxitInner = CGmaxit;

  ProgressMonitor::init( "Reading mesh" );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  mesh.refine( refine4, 2 );
  mesh.refine( refine9, 3 );
  if ( mapToUnitBall ) mesh.mapToUnitBall( );
  mesh.printInfo( );
  ProgressMonitor::step( );

  LO nNodes = mesh.getNNodes( );
  LO nElems = mesh.getNElements( );

  std::vector< LO > dirElems;
  std::vector< LO > neuElems;
  std::vector< LO > dirNodes;
  std::vector< LO > neuNodes;
  /*
  Vector< LO, SC > dir;
  Vector< LO, SC > neu;
  readBC< LO, SC >( bc.c_str( ), nElems, nNodes, dirElems,
      neuElems, dirNodes, neuNodes, dir, neu );
   */
  ///*
  Vector< LO, SC > dir( 3 * nNodes, true );
  Vector< LO, SC > neu( 3 * nElems, true );
  std::vector< LO > allElems( 3 * nElems );
  std::iota( allElems.begin( ), allElems.end( ), 0 );
  std::vector< LO > allNodes( 3 * nNodes );
  std::iota( allNodes.begin( ), allNodes.end( ), 0 );

  LO power = pow( 4, refine4 );
  LO elem [ 3 ];
  for ( LO i = 0; i < 2 * power; ++i ) {
    dirElems.push_back( i );
    dirElems.push_back( i + nElems );
    dirElems.push_back( i + 2 * nElems );
    mesh.getElement( i, elem );
    dirNodes.push_back( elem[ 0 ] );
    dirNodes.push_back( elem[ 1 ] );
    dirNodes.push_back( elem[ 2 ] );
    dirNodes.push_back( elem[ 0 ] + nNodes );
    dirNodes.push_back( elem[ 1 ] + nNodes );
    dirNodes.push_back( elem[ 2 ] + nNodes );
    dirNodes.push_back( elem[ 0 ] + 2 * nNodes );
    dirNodes.push_back( elem[ 1 ] + 2 * nNodes );
    dirNodes.push_back( elem[ 2 ] + 2 * nNodes );
    // potahnout v ose x
    dir.set( elem[ 0 ], 0.2 );
    dir.set( elem[ 1 ], 0.2 );
    dir.set( elem[ 2 ], 0.2 );
  }
  for ( LO i = 4 * power; i < 6 * power; ++i ) {
    dirElems.push_back( i );
    dirElems.push_back( i + nElems );
    dirElems.push_back( i + 2 * nElems );
    mesh.getElement( i, elem );
    dirNodes.push_back( elem[ 0 ] );
    dirNodes.push_back( elem[ 1 ] );
    dirNodes.push_back( elem[ 2 ] );
    dirNodes.push_back( elem[ 0 ] + nNodes );
    dirNodes.push_back( elem[ 1 ] + nNodes );
    dirNodes.push_back( elem[ 2 ] + nNodes );
    dirNodes.push_back( elem[ 0 ] + 2 * nNodes );
    dirNodes.push_back( elem[ 1 ] + 2 * nNodes );
    dirNodes.push_back( elem[ 2 ] + 2 * nNodes );
  }
  std::sort( dirNodes.begin( ), dirNodes.end( ) );
  std::sort( dirElems.begin( ), dirElems.end( ) );
  auto last = std::unique( dirNodes.begin( ), dirNodes.end( ) );
  dirNodes.erase( last, dirNodes.end( ) );

  std::set_difference( allElems.begin( ), allElems.end( ), dirElems.begin( ),
      dirElems.end( ), std::inserter( neuElems, neuElems.begin( ) ) );
  std::set_difference( allNodes.begin( ), allNodes.end( ), dirNodes.begin( ),
      dirNodes.end( ), std::inserter( neuNodes, neuNodes.begin( ) ) );
  //*/
  //IOHelper< LO, LO >::print( dirNodes );
  //IOHelper< LO, LO >::print( neuNodes );
  //IOHelper< LO, LO >::print( dirElems );
  //IOHelper< LO, LO >::print( neuElems );
  //dir.print( );
  //neu.print( );

  dir.scale( 0.1 );
  mesh.scale( 0.1 );
  /*
  Vector< LO, SC > bcn( 3 * nNodes );
  bcn.setAll( 0.0 );
  for ( auto it = neuNodes.begin( ); it != neuNodes.end( ); ++it ) {
    bcn.set( 3 * ( *it % nNodes ) + *it / nNodes, 1.0 );
  }
  Vector< LO, SC > bce( 3 * nElems );
  bce.setAll( 0.0 );
  for ( auto it = neuElems.begin( ); it != neuElems.end( ); ++it ) {
    bce.set( 3 * ( *it % nElems ) + *it / nElems, 1.0 );
  }
  std::vector< std::string > vnN{ "bc" };
  std::vector< Vector< LO, SC > * > vnD{ &bcn };
  std::vector< std::string > veN{ "bc" };
  std::vector< Vector< LO, SC > * > veD{ &bce };
  mesh.printParaviewVtu( "output/output.vtu", nullptr, nullptr, nullptr,
      nullptr, &vnN, &vnD, &veN, &veD );
  return;
   */

  std::vector< LO > row, col;
  std::vector< SC > val;

  ProgressMonitor::init( "Assembling P" );
  for ( LO i = 0; i < dirElems.size( ); ++i ) {
    col.push_back( i );
    row.push_back( dirElems[ i ] );
    val.push_back( 1.0 );
  }
  SparseMatrix< LO, SC > P( 3 * nElems, dirElems.size( ), row, col, val );
  ProgressMonitor::step( );

  row.clear( );
  col.clear( );
  val.clear( );

  ProgressMonitor::init( "Assembling R" );
  for ( LO i = 0; i < neuElems.size( ); ++i ) {
    col.push_back( i );
    row.push_back( neuElems[ i ] );
    val.push_back( 1.0 );
  }
  SparseMatrix< LO, SC > R( 3 * nElems, neuElems.size( ), row, col, val );
  ProgressMonitor::step( );

  row.clear( );
  col.clear( );
  val.clear( );

  ProgressMonitor::init( "Assembling S" );
  for ( LO i = 0; i < dirNodes.size( ); ++i ) {
    col.push_back( i );
    row.push_back( dirNodes[ i ] );
    val.push_back( 1.0 );
  }
  SparseMatrix< LO, SC > S( 3 * nNodes, dirNodes.size( ), row, col, val );
  ProgressMonitor::step( );

  row.clear( );
  col.clear( );
  val.clear( );

  ProgressMonitor::init( "Assembling Q" );
  for ( LO i = 0; i < neuNodes.size( ); ++i ) {
    col.push_back( i );
    row.push_back( neuNodes[ i ] );
    val.push_back( 1.0 );
  }
  SparseMatrix< LO, SC > Q( 3 * nNodes, neuNodes.size( ), row, col, val );
  ProgressMonitor::step( );

  int quadNear[ 4 ];
  quadratureType quadType;
  if ( quadTypeInt == 0 ) {
    quadType = Steinbach;
    quadNear[ 0 ] = quadNear[ 1 ] = orderNear;
    std::cout << "Using Steinbach quadrature, order " << orderNear << "."
        << std::endl;
  } else {
    quadType = SauterSchwab;
    quadNear[ 0 ] = quadNear[ 1 ] = quadNear[ 2 ] = quadNear[ 3 ] = orderNear;
    std::cout << "Using Sauter-Schwab quadrature, order " << orderNear << "."
        << std::endl;
  }

  int * quadFar = nullptr;
  if ( orderFar > 0 ) {
    quadFar = new int[ 2 ];
    quadFar[ 0 ] = quadFar[ 1 ] = orderFar;
  }

  BESpace< LO, SC > bespace00( &mesh, p0, p0 );
  BESpace< LO, SC > bespace10( &mesh, p1, p0 );
  BESpace< LO, SC > bespace11( &mesh, p1, p1 );

  SparseMatrix<LO, SC> T12, T13, T23, T;
  mesh.assembleT( T12, T13, T23, T );
  std::vector< SparseMatrix< LO, SC > * > Tv;
  Tv.push_back( &T12 );
  Tv.push_back( &T13 );
  Tv.push_back( &T23 );
  Tv.push_back( &T );

  FullMatrix< LO, SC > * Vlap = new FullMatrix< LO, SC >( 0, 0 );
  FullMatrix< LO, SC > * V = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLame1Layer< LO, SC > formV( &bespace00, quadNear, quadType,
      quadFar, false );
  formV.setE( E );
  formV.setNu( nu );
#if N_MIC < 1
  ProgressMonitor::init( "V, Vlap" );
  formV.assemble( *V, *Vlap );
  ProgressMonitor::step( );
#endif

  FullMatrix< LO, SC > * Klap = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLaplace2Layer< LO, SC > formKlap( &bespace10, quadNear,
      quadType, quadFar );
#if N_MIC < 1
  ProgressMonitor::init( "Klap" );
  formKlap.assemble( *Klap );
#else
  ProgressMonitor::init( "V, VLap, Klap" );
  formV.assembleAllMIC( *V, *Vlap, *Klap, bespace10 );
#endif
  ProgressMonitor::step( );

  FullMatrix< LO, SC > * K = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLame2Layer< LO, SC > formK( &bespace10, quadNear, quadType,
      quadFar );
  formK.setE( E );
  formK.setNu( nu );
  ProgressMonitor::init( "K" );
  formK.assemble( *K, *Vlap, *V, *Klap, Tv );
  ProgressMonitor::step( );

  delete Klap;

  FullMatrix< LO, SC > * D = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLameHypersingular< LO, SC > formD( &bespace11, quadNear,
      quadType, quadFar );
  formD.setE( E );
  formD.setNu( nu );
  ProgressMonitor::init( "D" );
  formD.assemble( *D, *Vlap, *V, Tv );
  ProgressMonitor::step( );

  delete Vlap;

  IdentityOperator< LO, SC > id( &bespace10 );

  Zero< LO, SC > O( nElems, nNodes );

  CompoundLinearOperator< LO, SC > pvp;
  pvp.push_back( &P );
  pvp.push_back( V );
  pvp.push_back( &P, true );
  if ( !pvp.isValid( ) ) std::cout << "pvp invalid!" << std::endl;

  CompoundLinearOperator< LO, SC > pkq;
  pkq.push_back( &Q );
  pkq.push_back( K );
  pkq.push_back( &P, true );
  if ( !pvp.isValid( ) ) std::cout << "pkq invalid!" << std::endl;

  CompoundLinearOperator< LO, SC > qdq;
  qdq.push_back( &Q );
  qdq.push_back( D );
  qdq.push_back( &Q, true );
  if ( !qdq.isValid( ) ) std::cout << "qdq invalid!" << std::endl;

  CompoundLinearOperator< LO, SC > pvr;
  pvr.push_back( &R );
  pvr.push_back( V );
  pvr.push_back( &P, true );
  if ( !pvr.isValid( ) ) std::cout << "pvr invalid!" << std::endl;

  CompoundLinearOperator< LO, SC > qds;
  qds.push_back( &S );
  qds.push_back( D );
  qds.push_back( &Q, true );
  if ( !qds.isValid( ) ) std::cout << "qds invalid!" << std::endl;

  Lame1LayerP0P0MultilvlPrecond<LO, SC> C( &mesh, 15 );
  /*
  Laplace1LayerP0P0MultilvlPrecond<LO, SC> Vprecond( &mesh, 15 );
  Zero< LO, SC > O_V( nElems, nElems );

  BlockLinearOperator< LO, SC > C( 3, 3 );
  C.setBlock( 0, 0, &Vprecond );
  C.setBlock( 0, 1, &O_V );
  C.setBlock( 0, 2, &O_V );
  C.setBlock( 1, 0, &O_V );
  C.setBlock( 1, 1, &Vprecond );
  C.setBlock( 1, 2, &O_V );
  C.setBlock( 2, 0, &O_V );
  C.setBlock( 2, 1, &O_V );
  C.setBlock( 2, 2, &Vprecond );
  if ( !C.isValid( ) ) std::cout << "C invalid!" << std::endl;
   */

  CompoundLinearOperator< LO, SC > pcp;
  pcp.push_back( &P );
  pcp.push_back( &C );
  pcp.push_back( &P, true );
  if ( !pcp.isValid( ) ) std::cout << "pcp invalid!" << std::endl;

  BlockLinearOperator< LO, SC > M( 3, 3 );
  M.setBlock( 0, 0, &id );
  M.setBlock( 0, 1, &O );
  M.setBlock( 0, 2, &O );
  M.setBlock( 1, 0, &O );
  M.setBlock( 1, 1, &id );
  M.setBlock( 1, 2, &O );
  M.setBlock( 2, 0, &O );
  M.setBlock( 2, 1, &O );
  M.setBlock( 2, 2, &id );
  if ( !M.isValid( ) ) std::cout << "M invalid!" << std::endl;

  SumLinearOperator< LO, SC > mplusk;
  mplusk.push_back( &M, false, 0.5 );
  mplusk.push_back( K );
  if ( !mplusk.isValid( ) ) std::cout << "mplusk invalid!" << std::endl;

  SumLinearOperator< LO, SC > mminusk;
  mminusk.push_back( &M, false, 0.5 );
  mminusk.push_back( K, false, -1.0 );
  if ( !mminusk.isValid( ) ) std::cout << "mminusk invalid!" << std::endl;

  CompoundLinearOperator< LO, SC > pmks;
  pmks.push_back( &S );
  pmks.push_back( &mplusk );
  pmks.push_back( &P, true );
  if ( !pmks.isValid( ) ) std::cout << "pmks invalid!" << std::endl;

  CompoundLinearOperator< LO, SC > rmkq;
  rmkq.push_back( &Q );
  rmkq.push_back( &mminusk );
  rmkq.push_back( &R, true );
  if ( !rmkq.isValid( ) ) std::cout << "rmkq invalid!" << std::endl;

  SolverParameters< LO, SC > params;
  params.solver = SolverParameters< LO, SC >::CG;
  params.maxIter = CGmaxitInner;
  params.precision = CGprecInner;
  params.msg = "pvpinv";
  //InverseLinearOperator< LO, SC > pvpinv( &pvp, &params );
  InverseLinearOperator< LO, SC > pvpinv( &pvp, &params, &pcp );

  CompoundLinearOperator< LO, SC > sp1;
  sp1.push_back( &pkq );
  sp1.push_back( &pvpinv );
  sp1.push_back( &pkq, true );
  if ( !sp1.isValid( ) ) std::cout << "sp1 invalid!" << std::endl;

  SumLinearOperator< LO, SC > sp;
  sp.push_back( &sp1 );
  sp.push_back( &qdq );
  if ( !sp.isValid( ) ) std::cout << "sp invalid!" << std::endl;

  BlockLinearOperator< LO, SC > rhsblock( 2, 2 );
  rhsblock.setBlock( 0, 0, &pvr, false, -1.0 );
  rhsblock.setBlock( 0, 1, &pmks );
  rhsblock.setBlock( 1, 0, &rmkq, true );
  rhsblock.setBlock( 1, 1, &qds, false, -1.0 );
  if ( !rhsblock.isValid( ) ) std::cout << "rhsblock invalid!" << std::endl;

  Vector< LO, SC > h( neuElems.size( ) );
  Vector< LO, SC > g( dirNodes.size( ) );
  R.apply( neu, h, true );
  S.apply( dir, g, true );
  Vector< LO, SC > hg( h );
  hg.append( g );

  Vector< LO, SC > rhs( rhsblock.getDimRange( ) );
  rhsblock.apply( hg, rhs );

  Vector< LO, SC > u( dirElems.size( ), rhs.getData( ) );
  Vector< LO, SC > v( neuNodes.size( ), rhs.getData( ) + dirElems.size( ) );

  CompoundLinearOperator< LO, SC > usol;
  usol.push_back( &pvpinv );
  usol.push_back( &pkq, true, -1.0 );
  if ( !usol.isValid( ) ) std::cout << "usol invalid!" << std::endl;

  usol.apply( u, v, false, 1.0, 1.0 );

  FullMatrix< LO, SC > * Vlap11 = new FullMatrix< LO, SC >( 0, 0 );
  FullMatrix< LO, SC > * V11 = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLame1Layer< LO, SC > formV11( &bespace11, quadNear,
      quadType, quadFar, false );
  formV11.setE( E );
  formV11.setNu( nu );
  ProgressMonitor::init( "V11" );
  formV11.assemble( *V11, *Vlap11 );
  ProgressMonitor::step( );

  delete Vlap11;

  IdentityOperator< LO, SC > id11( &bespace11 );
  SolverParameters< LO, SC > params_prec;
  params_prec.solver = SolverParameters< LO, SC >::CG;
  params_prec.maxIter = CGmaxitInner;
  params_prec.precision = CGprecInner;
  params_prec.msg = "qid11qinv";
  CompoundLinearOperator< LO, SC > qv11q;
  qv11q.push_back( &Q );
  qv11q.push_back( V11 );
  qv11q.push_back( &Q, true );
  if ( !qv11q.isValid( ) ) std::cout << "qv11q invalid!" << std::endl;
  Zero< LO, SC > O11( nNodes, nNodes );
  BlockLinearOperator< LO, SC > M11( 3, 3 );
  M11.setBlock( 0, 0, &id11 );
  M11.setBlock( 0, 1, &O11 );
  M11.setBlock( 0, 2, &O11 );
  M11.setBlock( 1, 0, &O11 );
  M11.setBlock( 1, 1, &id11 );
  M11.setBlock( 1, 2, &O11 );
  M11.setBlock( 2, 0, &O11 );
  M11.setBlock( 2, 1, &O11 );
  M11.setBlock( 2, 2, &id11 );
  if ( !M11.isValid( ) ) std::cout << "M11 invalid!" << std::endl;
  CompoundLinearOperator< LO, SC > qid11q;
  qid11q.push_back( &Q );
  qid11q.push_back( &M11 );
  qid11q.push_back( &Q, true );
  if ( !qid11q.isValid( ) ) std::cout << "qid11q invalid!" << std::endl;
  InverseLinearOperator< LO, SC > qid11qinv( &qid11q, &params_prec );
  CompoundLinearOperator< LO, SC > prec;
  prec.push_back( &qid11qinv );
  prec.push_back( &qv11q );
  prec.push_back( &qid11qinv );
  if ( !qid11q.isValid( ) ) std::cout << "qid11qinv invalid!" << std::endl;

  Vector< LO, SC > t( neuNodes.size( ) );
  //sp.CGSolve( v, t, CGprec, CGmaxit, nullptr );
  //t.setAll( 0.0 );
  sp.CGSolve( v, t, CGprec, CGmaxit, &prec, "SP solve" );

  pkq.apply( t, u, false, 1.0, 1.0 );

  Vector< LO, SC > s( dirElems.size( ) );
  pvpinv.apply( u, s );

  P.apply( s, neu, false, 1.0, 1.0 );
  Q.apply( t, dir, false, 1.0, 1.0 );

  Vector< LO, SC > dir2( 3 * nNodes );
  Vector< LO, SC > neu2( 3 * nElems );
  for ( LO i = 0; i < nElems; ++i ) {
    neu2.set( 3 * i, neu.get( i ) );
    neu2.set( 3 * i + 1, neu.get( i + nElems ) );
    neu2.set( 3 * i + 2, neu.get( i + 2 * nElems ) );
  }
  for ( LO i = 0; i < nNodes; ++i ) {
    dir2.set( 3 * i, dir.get( i ) );
    dir2.set( 3 * i + 1, dir.get( i + nNodes ) );
    dir2.set( 3 * i + 2, dir.get( i + 2 * nNodes ) );
  }

  Vector< LO, SC > bc_nodes( 3 * nNodes );
  bc_nodes.setAll( 0.0 );
  for ( auto it = neuNodes.begin( ); it != neuNodes.end( ); ++it ) {
    bc_nodes.set( 3 * ( *it % nNodes ) + *it / nNodes, 1.0 );
  }
  Vector< LO, SC > bc_elems( 3 * nElems );
  bc_elems.setAll( 0.0 );
  for ( auto it = neuElems.begin( ); it != neuElems.end( ); ++it ) {
    bc_elems.set( 3 * ( *it % nElems ) + *it / nElems, 1.0 );
  }

  std::vector< std::string > vnName{ "bc", "displacement" };
  std::vector< Vector< LO, SC > * > vnData{ &bc_nodes, &dir2 };
  std::vector< std::string > veName{ "bc", "stress" };
  std::vector< Vector< LO, SC > * > veData{ &bc_elems, &neu2 };

  mesh.printParaviewVtu( "output/output.vtu", nullptr, nullptr, nullptr,
      nullptr, &vnName, &vnData, &veName, &veData );

  delete V;
  delete K;
  delete D;
  delete V11;
  if ( quadFar ) delete [] quadFar;
}

//void testLameMixed(
//    const string & filename,
//    const string & bc,
//    int refine4,
//    int refine9,
//    bool mapToUnitBall,
//    int quadTypeInt,
//    int orderNear,
//    int orderFar
//    ) {
//
//  typedef double SC;
//  typedef int LO;
//
//  //std::cout.precision( 2 );
//  //std::cout.setf( std::ios::scientific );
//
//  SC nu = 0.33;
//  SC E = 1.1e5;
//
//  ProgressMonitor::init( "Reading mesh" );
//  SurfaceMesh3D< LO, SC > mesh;
//  mesh.load( filename.c_str( ) );
//  mesh.refine( refine4, 2 );
//  mesh.refine( refine9, 3 );
//  if ( mapToUnitBall ) mesh.mapToUnitBall( );
//  mesh.printInfo( );
//  ProgressMonitor::step( );
//
//  LO nNodes = mesh.getNNodes( );
//  LO nElems = mesh.getNElements( );
//
//  std::vector< LO > dirElems;
//  std::vector< LO > neuElems;
//  std::vector< LO > dirNodes;
//  std::vector< LO > neuNodes;
//  /*
//  Vector< LO, SC > dir;
//  Vector< LO, SC > neu;
//  readBC< LO, SC >( bc.c_str( ), nElems, nNodes, dirElems,
//      neuElems, dirNodes, neuNodes, dir, neu );
//   */
//  ///*
//  Vector< LO, SC > dir( 3 * nNodes, true );
//  Vector< LO, SC > neu( 3 * nElems, true );
//  std::vector< LO > allElems( 3 * nElems );
//  std::iota( allElems.begin( ), allElems.end( ), 0 );
//  std::vector< LO > allNodes( 3 * nNodes );
//  std::iota( allNodes.begin( ), allNodes.end( ), 0 );
//
//  LO power = pow( 4, refine4 );
//  LO elem [ 3 ];
//  for ( LO i = 0; i < 2 * power; ++i ) {
//    dirElems.push_back( i );
//    dirElems.push_back( i + nElems );
//    dirElems.push_back( i + 2 * nElems );
//    mesh.getElement( i, elem );
//    dirNodes.push_back( elem[ 0 ] );
//    dirNodes.push_back( elem[ 1 ] );
//    dirNodes.push_back( elem[ 2 ] );
//    dirNodes.push_back( elem[ 0 ] + nNodes );
//    dirNodes.push_back( elem[ 1 ] + nNodes );
//    dirNodes.push_back( elem[ 2 ] + nNodes );
//    dirNodes.push_back( elem[ 0 ] + 2 * nNodes );
//    dirNodes.push_back( elem[ 1 ] + 2 * nNodes );
//    dirNodes.push_back( elem[ 2 ] + 2 * nNodes );
//    // potahnout v ose x
//    dir.set( elem[ 0 ], 0.2 );
//    dir.set( elem[ 1 ], 0.2 );
//    dir.set( elem[ 2 ], 0.2 );
//  }
//  for ( LO i = 4 * power; i < 6 * power; ++i ) {
//    dirElems.push_back( i );
//    dirElems.push_back( i + nElems );
//    dirElems.push_back( i + 2 * nElems );
//    mesh.getElement( i, elem );
//    dirNodes.push_back( elem[ 0 ] );
//    dirNodes.push_back( elem[ 1 ] );
//    dirNodes.push_back( elem[ 2 ] );
//    dirNodes.push_back( elem[ 0 ] + nNodes );
//    dirNodes.push_back( elem[ 1 ] + nNodes );
//    dirNodes.push_back( elem[ 2 ] + nNodes );
//    dirNodes.push_back( elem[ 0 ] + 2 * nNodes );
//    dirNodes.push_back( elem[ 1 ] + 2 * nNodes );
//    dirNodes.push_back( elem[ 2 ] + 2 * nNodes );
//  }
//  std::sort( dirNodes.begin( ), dirNodes.end( ) );
//  std::sort( dirElems.begin( ), dirElems.end( ) );
//  auto last = std::unique( dirNodes.begin( ), dirNodes.end( ) );
//  dirNodes.erase( last, dirNodes.end( ) );
//
//  std::set_difference( allElems.begin( ), allElems.end( ), dirElems.begin( ),
//      dirElems.end( ), std::inserter( neuElems, neuElems.begin( ) ) );
//  std::set_difference( allNodes.begin( ), allNodes.end( ), dirNodes.begin( ),
//      dirNodes.end( ), std::inserter( neuNodes, neuNodes.begin( ) ) );
//  //*/
//  //IOHelper< LO, LO >::print( dirNodes );
//  //IOHelper< LO, LO >::print( neuNodes );
//  //IOHelper< LO, LO >::print( dirElems );
//  //IOHelper< LO, LO >::print( neuElems );
//  //dir.print( );
//  //neu.print( );
//
//  /*
//  Vector< LO, SC > bcn( 3 * nNodes );
//  bcn.setAll( 0.0 );
//  for ( auto it = neuNodes.begin( ); it != neuNodes.end( ); ++it ) {
//    bcn.set( 3 * ( *it % nNodes ) + *it / nNodes, 1.0 );
//  }
//  Vector< LO, SC > bce( 3 * nElems );
//  bce.setAll( 0.0 );
//  for ( auto it = neuElems.begin( ); it != neuElems.end( ); ++it ) {
//    bce.set( 3 * ( *it % nElems ) + *it / nElems, 1.0 );
//  }
//  std::vector< std::string > vnN{ "bc" };
//  std::vector< Vector< LO, SC > * > vnD{ &bcn };
//  std::vector< std::string > veN{ "bc" };
//  std::vector< Vector< LO, SC > * > veD{ &bce };
//  mesh.printParaviewVtu( "output/output.vtu", nullptr, nullptr, nullptr,
//      nullptr, &vnN, &vnD, &veN, &veD );
//  return;
//   */
//
//  std::vector< LO > row, col;
//  std::vector< SC > val;
//
//  ProgressMonitor::init( "Assembling P" );
//  for ( LO i = 0; i < dirElems.size( ); ++i ) {
//    col.push_back( i );
//    row.push_back( dirElems[ i ] );
//    val.push_back( 1.0 );
//  }
//  SparseMatrix< LO, SC > P( 3 * nElems, dirElems.size( ), row, col, val );
//  ProgressMonitor::step( );
//
//  row.clear( );
//  col.clear( );
//  val.clear( );
//
//  ProgressMonitor::init( "Assembling R" );
//  for ( LO i = 0; i < neuElems.size( ); ++i ) {
//    col.push_back( i );
//    row.push_back( neuElems[ i ] );
//    val.push_back( 1.0 );
//  }
//  SparseMatrix< LO, SC > R( 3 * nElems, neuElems.size( ), row, col, val );
//  ProgressMonitor::step( );
//
//  row.clear( );
//  col.clear( );
//  val.clear( );
//
//  ProgressMonitor::init( "Assembling S" );
//  for ( LO i = 0; i < dirNodes.size( ); ++i ) {
//    col.push_back( i );
//    row.push_back( dirNodes[ i ] );
//    val.push_back( 1.0 );
//  }
//  SparseMatrix< LO, SC > S( 3 * nNodes, dirNodes.size( ), row, col, val );
//  ProgressMonitor::step( );
//
//  row.clear( );
//  col.clear( );
//  val.clear( );
//
//  ProgressMonitor::init( "Assembling Q" );
//  for ( LO i = 0; i < neuNodes.size( ); ++i ) {
//    col.push_back( i );
//    row.push_back( neuNodes[ i ] );
//    val.push_back( 1.0 );
//  }
//  SparseMatrix< LO, SC > Q( 3 * nNodes, neuNodes.size( ), row, col, val );
//  ProgressMonitor::step( );
//
//  int quadNear[ 4 ];
//  quadratureType quadType;
//  if ( quadTypeInt == 0 ) {
//    quadType = Steinbach;
//    quadNear[ 0 ] = quadNear[ 1 ] = orderNear;
//    std::cout << "Using Steinbach quadrature, order " << orderNear << "."
//        << std::endl;
//  } else {
//    quadType = SauterSchwab;
//    quadNear[ 0 ] = quadNear[ 1 ] = quadNear[ 2 ] = quadNear[ 3 ] = orderNear;
//    std::cout << "Using Sauter-Schwab quadrature, order " << orderNear << "."
//        << std::endl;
//  }
//
//  int * quadFar = nullptr;
//  if ( orderFar > 0 ) {
//    quadFar = new int[ 2 ];
//    quadFar[ 0 ] = quadFar[ 1 ] = orderFar;
//  }
//
//  BESpace< LO, SC > bespace00( &mesh, p0, p0 );
//  BESpace< LO, SC > bespace10( &mesh, p1, p0 );
//  BESpace< LO, SC > bespace11( &mesh, p1, p1 );
//
//  SparseMatrix<LO, SC> T12, T13, T23, T;
//  mesh.assembleT( T12, T13, T23, T );
//  std::vector< SparseMatrix< LO, SC > * > Tv;
//  Tv.push_back( &T12 );
//  Tv.push_back( &T13 );
//  Tv.push_back( &T23 );
//  Tv.push_back( &T );
//
//  bool symmetricV = false;
//  FullMatrix< LO, SC > * Vlap = new FullMatrix< LO, SC >( 0, 0 );
//  FullMatrix< LO, SC > * V = new FullMatrix< LO, SC >( 0, 0 );
//  BEBilinearFormLame1Layer< LO, SC > formV( &bespace00, quadNear, quadType,
//      quadFar, symmetricV );
//  formV.setE( E );
//  formV.setNu( nu );
//#if N_MIC < 1
//  ProgressMonitor::init( "V, Vlap" );
//  formV.assemble( *V, *Vlap );
//  ProgressMonitor::step( );
//#endif
//
//  FullMatrix< LO, SC > * Klap = new FullMatrix< LO, SC >( 0, 0 );
//  BEBilinearFormLaplace2Layer< LO, SC > formKlap( &bespace10, quadNear,
//      quadType, quadFar );
//#if N_MIC < 1
//  ProgressMonitor::init( "Klap" );
//  formKlap.assemble( *Klap );
//#else
//  ProgressMonitor::init( "V, VLap, Klap" );
//  formV.assembleAllMIC( *V, *Vlap, *Klap, bespace10 );
//#endif
//  ProgressMonitor::step( );
//
//  FullMatrix< LO, SC > * K = new FullMatrix< LO, SC >( 0, 0 );
//  BEBilinearFormLame2Layer< LO, SC > formK( &bespace10, quadNear, quadType,
//      quadFar );
//  formK.setE( E );
//  formK.setNu( nu );
//  ProgressMonitor::init( "K" );
//  formK.assemble( *K, *Vlap, *V, *Klap, Tv );
//  ProgressMonitor::step( );
//
//  delete Klap;
//
//  FullMatrix< LO, SC > * D = new FullMatrix< LO, SC >( 0, 0 );
//  BEBilinearFormLameHypersingular< LO, SC > formD( &bespace11, quadNear,
//      quadType, quadFar );
//  formD.setE( E );
//  formD.setNu( nu );
//  ProgressMonitor::init( "D" );
//  formD.assemble( *D, *Vlap, *V, Tv );
//  ProgressMonitor::step( );
//
//  delete Vlap;
//
//  IdentityOperator< LO, SC > id( &bespace10 );
//  SparseMatrix< LO, SC > Mlap;
//  id.assemble( Mlap );
//
//  Zero< LO, SC > O( nElems, nNodes );
//
//  CompoundLinearOperator< LO, SC > pvr;
//  pvr.push_back( &R );
//  pvr.push_back( V );
//  pvr.push_back( &P, true );
//  if ( !pvr.isValid( ) ) std::cout << "pvr invalid!" << std::endl;
//
//  CompoundLinearOperator< LO, SC > qds;
//  qds.push_back( &S );
//  qds.push_back( D );
//  qds.push_back( &Q, true );
//  if ( !qds.isValid( ) ) std::cout << "qds invalid!" << std::endl;
//
//  BlockLinearOperator< LO, SC > M( 3, 3 );
//  M.setBlock( 0, 0, &Mlap );
//  M.setBlock( 0, 1, &O );
//  M.setBlock( 0, 2, &O );
//  M.setBlock( 1, 0, &O );
//  M.setBlock( 1, 1, &Mlap );
//  M.setBlock( 1, 2, &O );
//  M.setBlock( 2, 0, &O );
//  M.setBlock( 2, 1, &O );
//  M.setBlock( 2, 2, &Mlap );
//  if ( !M.isValid( ) ) std::cout << "M invalid!" << std::endl;
//
//  SumLinearOperator< LO, SC > mplusk;
//  mplusk.push_back( &M, false, 0.5 );
//  mplusk.push_back( K );
//  if ( !mplusk.isValid( ) ) std::cout << "mplusk invalid!" << std::endl;
//
//  SumLinearOperator< LO, SC > mminusk;
//  mminusk.push_back( &M, false, 0.5 );
//  mminusk.push_back( K, false, -1.0 );
//  if ( !mminusk.isValid( ) ) std::cout << "mminusk invalid!" << std::endl;
//
//  CompoundLinearOperator< LO, SC > pmks;
//  pmks.push_back( &S );
//  pmks.push_back( &mplusk );
//  pmks.push_back( &P, true );
//  if ( !pmks.isValid( ) ) std::cout << "pmks invalid!" << std::endl;
//
//  CompoundLinearOperator< LO, SC > rmkq;
//  rmkq.push_back( &Q );
//  rmkq.push_back( &mminusk );
//  rmkq.push_back( &R, true );
//  if ( !rmkq.isValid( ) ) std::cout << "rmkq invalid!" << std::endl;
//
//  BlockLinearOperator< LO, SC > rhsblock( 2, 2 );
//  rhsblock.setBlock( 0, 0, &pvr, false, -1.0 );
//  rhsblock.setBlock( 0, 1, &pmks );
//  rhsblock.setBlock( 1, 0, &rmkq, true );
//  rhsblock.setBlock( 1, 1, &qds, false, -1.0 );
//  if ( !rhsblock.isValid( ) ) std::cout << "rhsblock invalid!" << std::endl;
//
//  Vector< LO, SC > h( neuElems.size( ) );
//  Vector< LO, SC > g( dirNodes.size( ) );
//  R.apply( neu, h, true );
//  S.apply( dir, g, true );
//  Vector< LO, SC > hg( h );
//  hg.append( g );
//
//  Vector< LO, SC > rhs( rhsblock.getDimRange( ) );
//  rhsblock.apply( hg, rhs );
//
//  FullMatrix< LO, SC > block( dirElems.size( ) + neuNodes.size( ),
//      dirElems.size( ) + neuNodes.size( ) );
//  for ( LO i = 0; i < dirElems.size( ); ++i ) {
//    for ( LO j = 0; j < dirElems.size( ); ++j ) {
//      block.set( i, j, V->get( dirElems[ i ], dirElems[ j ] ) );
//    }
//  }
//  for ( LO i = 0; i < dirElems.size( ); ++i ) {
//    for ( LO j = 0; j < neuNodes.size( ); ++j ) {
//      block.set( i, j + dirElems.size( ),
//          -K->get( dirElems[ i ], neuNodes[ j ] ) );
//      block.set( j + dirElems.size( ), i,
//          K->get( dirElems[ i ], neuNodes[ j ] ) );
//    }
//  }
//  for ( LO i = 0; i < neuNodes.size( ); ++i ) {
//    for ( LO j = 0; j < neuNodes.size( ); ++j ) {
//      block.set( i + dirElems.size( ), j + dirElems.size( ), 
//          D->get( neuNodes[ i ], neuNodes[ j ] ) );
//    }
//  }
//
//  delete V;
//  delete K;
//  delete D;
//
//  Vector< LO, SC > sol( rhs );
//  block.LUSolve( sol );
//
//  Vector< LO, SC > s( dirElems.size( ), sol.getData( ) );
//  Vector< LO, SC > t( neuNodes.size( ),
//      sol.getData( ) + dirElems.size( ) );
//  P.apply( s, neu, false, 1.0, 1.0 );
//  Q.apply( t, dir, false, 1.0, 1.0 );
//
//  Vector< LO, SC > dir2( 3 * nNodes );
//  Vector< LO, SC > neu2( 3 * nElems );
//  for ( LO i = 0; i < nElems; ++i ) {
//    neu2.set( 3 * i, neu.get( i ) );
//    neu2.set( 3 * i + 1, neu.get( i + nElems ) );
//    neu2.set( 3 * i + 2, neu.get( i + 2 * nElems ) );
//  }
//  for ( LO i = 0; i < nNodes; ++i ) {
//    dir2.set( 3 * i, dir.get( i ) );
//    dir2.set( 3 * i + 1, dir.get( i + nNodes ) );
//    dir2.set( 3 * i + 2, dir.get( i + 2 * nNodes ) );
//  }
//
//  Vector< LO, SC > bc_nodes( 3 * nNodes );
//  bc_nodes.setAll( 0.0 );
//  for ( auto it = neuNodes.begin( ); it != neuNodes.end( ); ++it ) {
//    bc_nodes.set( 3 * ( *it % nNodes ) + *it / nNodes, 1.0 );
//  }
//  Vector< LO, SC > bc_elems( 3 * nElems );
//  bc_elems.setAll( 0.0 );
//  for ( auto it = neuElems.begin( ); it != neuElems.end( ); ++it ) {
//    bc_elems.set( 3 * ( *it % nElems ) + *it / nElems, 1.0 );
//  }
//
//  std::vector< std::string > vnName{ "bc", "displacement" };
//  std::vector< Vector< LO, SC > * > vnData{ &bc_nodes, &dir2 };
//  std::vector< std::string > veName{ "bc", "stress" };
//  std::vector< Vector< LO, SC > * > veData{ &bc_elems, &neu2 };
//
//  mesh.printParaviewVtu( "output/output.vtu", nullptr, nullptr, nullptr,
//      nullptr, &vnName, &vnData, &veName, &veData );
//
//  if ( quadFar ) delete [] quadFar;
//}

//template< class LO, class SC >
//bool readBC(
//    std::string const & filename,
//    LO nElems,
//    LO nNodes,
//    std::vector< LO > & dirElems,
//    std::vector< LO > & neuElems,
//    std::vector< LO > & dirNodes,
//    std::vector< LO > & neuNodes,
//    Vector< LO, SC > & dir,
//    Vector< LO, SC > & neu
//    ) {
//
//  std::cout << "Reading file " << filename << " ... ";
//  std::ifstream file( filename.c_str( ) );
//
//  if ( !file.is_open( ) ) {
//    std::cout << "File could not be opened!" << std::endl;
//    return false;
//  }
//
//  dirElems.clear( );
//  neuElems.clear( );
//  dirNodes.clear( );
//  neuNodes.clear( );
//  dir.resize( 3 * nNodes, true );
//  neu.resize( 3 * nElems, true );
//
//  LO nDir, nNeu, ind;
//  SC bc;
//
//  file >> nDir;
//  for ( LO i = 0; i < nDir; ++i ) {
//    file >> ind;
//    for ( int j = 0; j < 3; ++j ) {
//      file >> bc;
//      if ( file.fail( ) ) {
//        file.clear( );
//        file.ignore( );
//        continue;
//      }
//      dirNodes.push_back( ind + j * nNodes );
//      dir.set( ind + j * nNodes, bc );
//    }
//  }
//
//  file >> nNeu;
//  for ( LO i = 0; i < nNeu; ++i ) {
//    file >> ind;
//    for ( int j = 0; j < 3; ++j ) {
//      file >> bc;
//      if ( file.fail( ) ) {
//        file.clear( );
//        file.ignore( );
//        continue;
//      }
//      neuElems.push_back( ind + j * nElems );
//      neu.set( ind + j * nElems, bc );
//    }
//  }
//
//  file.close( );
//
//  std::sort( dirNodes.begin( ), dirNodes.end( ) );
//  std::sort( neuElems.begin( ), neuElems.end( ) );
//
//  LO pos = 0;
//  bool check = true;
//  for ( LO i = 0; i < 3 * nElems; ++i ) {
//    if ( check && neuElems[ pos ] == i ) {
//      ++pos;
//      if ( pos == neuElems.size( ) ) check = false;
//    } else {
//      dirElems.push_back( i );
//    }
//  }
//
//  pos = 0;
//  check = true;
//  for ( LO i = 0; i < 3 * nNodes; ++i ) {
//    if ( check && dirNodes[ pos ] == i ) {
//      ++pos;
//      if ( pos == dirNodes.size( ) ) check = false;
//    } else {
//      neuNodes.push_back( i );
//    }
//  }
//
//  std::cout << "done." << std::endl;
//
//  return true;
//}
//
//
//
//
//
//template< class LO, class SC >
//void assembleT(
//    SurfaceMesh3D< LO, SC > & mesh,
//    SparseMatrix< LO, SC > & T1,
//    SparseMatrix< LO, SC > & T2,
//    SparseMatrix< LO, SC > & T3
//    ) {
//
//  LO nElems = mesh.getNElements( );
//  LO nNodes = mesh.getNNodes( );
//
//  std::vector< LO > rowInd;
//  rowInd.reserve( 3 * nElems );
//  std::vector< LO > colInd;
//  colInd.reserve( 3 * nElems );
//  std::vector< SC > values1;
//  values1.reserve( 3 * nElems );
//  std::vector< SC > values2;
//  values2.reserve( 3 * nElems );
//  std::vector< SC > values3;
//  values3.reserve( 3 * nElems );
//  LO elem[ 3 ];
//  Vector< LO, SC > * curls = mesh.getCurls( );
//
//  for ( LO i = 0; i < nElems; ++i ) {
//
//    mesh.getElement( i, elem );
//
//    for ( int node = 0; node < 3; ++node ) {
//
//      rowInd.push_back( i );
//      colInd.push_back( elem[ node ] );
//      values1.push_back( curls->get( 9 * i + 3 * node ) );
//      values2.push_back( curls->get( 9 * i + 3 * node + 1 ) );
//      values3.push_back( curls->get( 9 * i + 3 * node + 2 ) );
//    }
//  }
//
//  T1.setFromTriplets( nElems, nNodes, rowInd, colInd, values1 );
//  T2.setFromTriplets( nElems, nNodes, rowInd, colInd, values2 );
//  T3.setFromTriplets( nElems, nNodes, rowInd, colInd, values3 );
//
//}

//void testLameMixed(
//    const string & filename,
//    const string & bc,
//    int refine4,
//    int refine9,
//    bool mapToUnitBall,
//    int quadTypeInt,
//    int orderNear,
//    int orderFar
//    ) {
//
//  typedef double SC;
//  typedef int LO;
//
//  //std::cout.precision( 2 );
//  //std::cout.setf( std::ios::scientific );
//
//  SC nu = 0.33;
//  SC E = 1.1e5;
//  SC GMRESprec = 1e-8;
//  LO GMRESmaxit = 1000;
//
//  ProgressMonitor::init( "Reading mesh" );
//  SurfaceMesh3D< LO, SC > mesh;
//  mesh.load( filename.c_str( ) );
//  mesh.refine( refine4, 2 );
//  mesh.refine( refine9, 3 );
//  if ( mapToUnitBall ) mesh.mapToUnitBall( );
//  mesh.printInfo( );
//  ProgressMonitor::step( );
//
//  LO nNodes = mesh.getNNodes( );
//  LO nElems = mesh.getNElements( );
//
//  std::vector< LO > dirElems;
//  std::vector< LO > neuElems;
//  std::vector< LO > dirNodes;
//  std::vector< LO > neuNodes;
//  /*
//  Vector< LO, SC > dir;
//  Vector< LO, SC > neu;
//  readBC< LO, SC >( bc.c_str( ), nElems, nNodes, dirElems,
//      neuElems, dirNodes, neuNodes, dir, neu );
//  */
//  ///*
//  Vector< LO, SC > dir( 3 * nNodes, true );
//  Vector< LO, SC > neu( 3 * nElems, true );
//  std::vector< LO > allElems( 3 * nElems );
//  std::iota( allElems.begin( ), allElems.end( ), 0 );
//  std::vector< LO > allNodes( 3 * nNodes );
//  std::iota( allNodes.begin( ), allNodes.end( ), 0 );
//
//  LO power = pow( 4, refine4 );
//  LO elem [ 3 ];
//  for ( LO i = 0; i < 2 * power; ++i ) {
//    dirElems.push_back( i );
//    dirElems.push_back( i + nElems );
//    dirElems.push_back( i + 2 * nElems );
//    mesh.getElement( i, elem );
//    dirNodes.push_back( elem[ 0 ] );
//    dirNodes.push_back( elem[ 1 ] );
//    dirNodes.push_back( elem[ 2 ] );
//    dirNodes.push_back( elem[ 0 ] + nNodes );
//    dirNodes.push_back( elem[ 1 ] + nNodes );
//    dirNodes.push_back( elem[ 2 ] + nNodes );
//    dirNodes.push_back( elem[ 0 ] + 2 * nNodes );
//    dirNodes.push_back( elem[ 1 ] + 2 * nNodes );
//    dirNodes.push_back( elem[ 2 ] + 2 * nNodes );
//    // potahnout v ose x
//    dir.set( elem[ 0 ], 0.2 );
//    dir.set( elem[ 1 ], 0.2 );
//    dir.set( elem[ 2 ], 0.2 );
//  }
//  for ( LO i = 4 * power; i < 6 * power; ++i ) {
//    dirElems.push_back( i );
//    dirElems.push_back( i + nElems );
//    dirElems.push_back( i + 2 * nElems );
//    mesh.getElement( i, elem );
//    dirNodes.push_back( elem[ 0 ] );
//    dirNodes.push_back( elem[ 1 ] );
//    dirNodes.push_back( elem[ 2 ] );
//    dirNodes.push_back( elem[ 0 ] + nNodes );
//    dirNodes.push_back( elem[ 1 ] + nNodes );
//    dirNodes.push_back( elem[ 2 ] + nNodes );
//    dirNodes.push_back( elem[ 0 ] + 2 * nNodes );
//    dirNodes.push_back( elem[ 1 ] + 2 * nNodes );
//    dirNodes.push_back( elem[ 2 ] + 2 * nNodes );
//  }
//  std::sort( dirNodes.begin( ), dirNodes.end( ) );
//  std::sort( dirElems.begin( ), dirElems.end( ) );
//  auto last = std::unique( dirNodes.begin( ), dirNodes.end( ) );
//  dirNodes.erase( last, dirNodes.end( ) );
//
//  std::set_difference( allElems.begin( ), allElems.end( ), dirElems.begin( ),
//      dirElems.end( ), std::inserter( neuElems, neuElems.begin( ) ) );
//  std::set_difference( allNodes.begin( ), allNodes.end( ), dirNodes.begin( ),
//      dirNodes.end( ), std::inserter( neuNodes, neuNodes.begin( ) ) );
//  //*/
//  //IOHelper< LO, LO >::print( dirNodes );
//  //IOHelper< LO, LO >::print( neuNodes );
//  //IOHelper< LO, LO >::print( dirElems );
//  //IOHelper< LO, LO >::print( neuElems );
//  //dir.print( );
//  //neu.print( );
//
//  /*
//  Vector< LO, SC > bcn( 3 * nNodes );
//  bcn.setAll( 0.0 );
//  for ( auto it = neuNodes.begin( ); it != neuNodes.end( ); ++it ) {
//    bcn.set( 3 * ( *it % nNodes ) + *it / nNodes, 1.0 );
//  }
//  Vector< LO, SC > bce( 3 * nElems );
//  bce.setAll( 0.0 );
//  for ( auto it = neuElems.begin( ); it != neuElems.end( ); ++it ) {
//    bce.set( 3 * ( *it % nElems ) + *it / nElems, 1.0 );
//  }
//  std::vector< std::string > vnN{ "bc" };
//  std::vector< Vector< LO, SC > * > vnD{ &bcn };
//  std::vector< std::string > veN{ "bc" };
//  std::vector< Vector< LO, SC > * > veD{ &bce };
//  mesh.printParaviewVtu( "output/output.vtu", nullptr, nullptr, nullptr,
//      nullptr, &vnN, &vnD, &veN, &veD );
//  return;
//  */
//
//  std::vector< LO > row, col;
//  std::vector< SC > val;
//
//  ProgressMonitor::init( "Assembling P" );
//  for ( LO i = 0; i < dirElems.size( ); ++i ) {
//    col.push_back( i );
//    row.push_back( dirElems[ i ] );
//    val.push_back( 1.0 );
//  }
//  SparseMatrix< LO, SC > P( 3 * nElems, dirElems.size( ), row, col, val );
//  ProgressMonitor::step( );
//
//  row.clear( );
//  col.clear( );
//  val.clear( );
//
//  ProgressMonitor::init( "Assembling R" );
//  for ( LO i = 0; i < neuElems.size( ); ++i ) {
//    col.push_back( i );
//    row.push_back( neuElems[ i ] );
//    val.push_back( 1.0 );
//  }
//  SparseMatrix< LO, SC > R( 3 * nElems, neuElems.size( ), row, col, val );
//  ProgressMonitor::step( );
//
//  row.clear( );
//  col.clear( );
//  val.clear( );
//
//  ProgressMonitor::init( "Assembling S" );
//  for ( LO i = 0; i < dirNodes.size( ); ++i ) {
//    col.push_back( i );
//    row.push_back( dirNodes[ i ] );
//    val.push_back( 1.0 );
//  }
//  SparseMatrix< LO, SC > S( 3 * nNodes, dirNodes.size( ), row, col, val );
//  ProgressMonitor::step( );
//
//  row.clear( );
//  col.clear( );
//  val.clear( );
//
//  ProgressMonitor::init( "Assembling Q" );
//  for ( LO i = 0; i < neuNodes.size( ); ++i ) {
//    col.push_back( i );
//    row.push_back( neuNodes[ i ] );
//    val.push_back( 1.0 );
//  }
//  SparseMatrix< LO, SC > Q( 3 * nNodes, neuNodes.size( ), row, col, val );
//  ProgressMonitor::step( );
//
//  int quadNear[ 4 ];
//  quadratureType quadType;
//  if ( quadTypeInt == 0 ) {
//    quadType = Steinbach;
//    quadNear[ 0 ] = quadNear[ 1 ] = orderNear;
//    std::cout << "Using Steinbach quadrature, order " << orderNear << "."
//        << std::endl;
//  } else {
//    quadType = SauterSchwab;
//    quadNear[ 0 ] = quadNear[ 1 ] = quadNear[ 2 ] = quadNear[ 3 ] = orderNear;
//    std::cout << "Using Sauter-Schwab quadrature, order " << orderNear << "."
//        << std::endl;
//  }
//
//  int * quadFar = nullptr;
//  if ( orderFar > 0 ) {
//    quadFar = new int[ 2 ];
//    quadFar[ 0 ] = quadFar[ 1 ] = orderFar;
//  }
//
//  BESpace< LO, SC > bespace00( &mesh, p0, p0 );
//  BESpace< LO, SC > bespace10( &mesh, p1, p0 );
//  BESpace< LO, SC > bespace11( &mesh, p1, p1 );
//
//  SparseMatrix<LO, SC> T12, T13, T23, T;
//  mesh.assembleT( T12, T13, T23, T );
//  std::vector< SparseMatrix< LO, SC > * > Tv;
//  Tv.push_back( &T12 );
//  Tv.push_back( &T13 );
//  Tv.push_back( &T23 );
//  Tv.push_back( &T );
//
//  bool symmetricV = false;
//  FullMatrix< LO, SC > * Vlap = new FullMatrix< LO, SC >( 0, 0 );
//  FullMatrix< LO, SC > * V = new FullMatrix< LO, SC >( 0, 0 );
//  BEBilinearFormLame1Layer< LO, SC > formV( &bespace00, quadNear, quadType,
//      quadFar, symmetricV );
//  formV.setE( E );
//  formV.setNu( nu );
//#if N_MIC < 1
//  ProgressMonitor::init( "V, Vlap" );
//  formV.assemble( *V, *Vlap );
//  ProgressMonitor::step( );
//#endif
//
//  FullMatrix< LO, SC > * Klap = new FullMatrix< LO, SC >( 0, 0 );
//  BEBilinearFormLaplace2Layer< LO, SC > formKlap( &bespace10, quadNear,
//      quadType, quadFar );
//#if N_MIC < 1
//  ProgressMonitor::init( "Klap" );
//  formKlap.assemble( *Klap );
//#else
//  ProgressMonitor::init( "V, VLap, Klap" );
//  formV.assembleAllMIC( *V, *Vlap, *Klap, bespace10 );
//#endif
//  ProgressMonitor::step( );
//
//  FullMatrix< LO, SC > * K = new FullMatrix< LO, SC >( 0, 0 );
//  BEBilinearFormLame2Layer< LO, SC > formK( &bespace10, quadNear, quadType,
//      quadFar );
//  formK.setE( E );
//  formK.setNu( nu );
//  ProgressMonitor::init( "K" );
//  formK.assemble( *K, *Vlap, *V, *Klap, Tv );
//  ProgressMonitor::step( );
//
//  delete Klap;
//
//  FullMatrix< LO, SC > * D = new FullMatrix< LO, SC >( 0, 0 );
//  BEBilinearFormLameHypersingular< LO, SC > formD( &bespace11, quadNear,
//      quadType, quadFar );
//  formD.setE( E );
//  formD.setNu( nu );
//  ProgressMonitor::init( "D" );
//  formD.assemble( *D, *Vlap, *V, Tv );
//  ProgressMonitor::step( );
//
//  delete Vlap;
//
//  IdentityOperator< LO, SC > id( &bespace10 );
//
//  Zero< LO, SC > O( nElems, nNodes );
//
//  CompoundLinearOperator< LO, SC > pvp;
//  pvp.push_back( &P );
//  pvp.push_back( V );
//  pvp.push_back( &P, true );
//  if ( !pvp.isValid( ) ) std::cout << "pvp invalid!" << std::endl;
//
//  CompoundLinearOperator< LO, SC > pkq;
//  pkq.push_back( &Q );
//  pkq.push_back( K );
//  pkq.push_back( &P, true );
//  if ( !pvp.isValid( ) ) std::cout << "pkq invalid!" << std::endl;
//
//  CompoundLinearOperator< LO, SC > qdq;
//  qdq.push_back( &Q );
//  qdq.push_back( D );
//  qdq.push_back( &Q, true );
//  if ( !qdq.isValid( ) ) std::cout << "qdq invalid!" << std::endl;
//
//  CompoundLinearOperator< LO, SC > pvr;
//  pvr.push_back( &R );
//  pvr.push_back( V );
//  pvr.push_back( &P, true );
//  if ( !pvr.isValid( ) ) std::cout << "pvr invalid!" << std::endl;
//
//  CompoundLinearOperator< LO, SC > qds;
//  qds.push_back( &S );
//  qds.push_back( D );
//  qds.push_back( &Q, true );
//  if ( !qds.isValid( ) ) std::cout << "qds invalid!" << std::endl;
//
//  BlockLinearOperator< LO, SC > M( 3, 3 );
//  M.setBlock( 0, 0, &id );
//  M.setBlock( 0, 1, &O );
//  M.setBlock( 0, 2, &O );
//  M.setBlock( 1, 0, &O );
//  M.setBlock( 1, 1, &id );
//  M.setBlock( 1, 2, &O );
//  M.setBlock( 2, 0, &O );
//  M.setBlock( 2, 1, &O );
//  M.setBlock( 2, 2, &id );
//  if ( !M.isValid( ) ) std::cout << "M invalid!" << std::endl;
//
//  SumLinearOperator< LO, SC > mplusk;
//  mplusk.push_back( &M, false, 0.5 );
//  mplusk.push_back( K );
//  if ( !mplusk.isValid( ) ) std::cout << "mplusk invalid!" << std::endl;
//
//  SumLinearOperator< LO, SC > mminusk;
//  mminusk.push_back( &M, false, 0.5 );
//  mminusk.push_back( K, false, -1.0 );
//  if ( !mminusk.isValid( ) ) std::cout << "mminusk invalid!" << std::endl;
//
//  CompoundLinearOperator< LO, SC > pmks;
//  pmks.push_back( &S );
//  pmks.push_back( &mplusk );
//  pmks.push_back( &P, true );
//  if ( !pmks.isValid( ) ) std::cout << "pmks invalid!" << std::endl;
//
//  CompoundLinearOperator< LO, SC > rmkq;
//  rmkq.push_back( &Q );
//  rmkq.push_back( &mminusk );
//  rmkq.push_back( &R, true );
//  if ( !rmkq.isValid( ) ) std::cout << "rmkq invalid!" << std::endl;
//
//  BlockLinearOperator< LO, SC > block( 2, 2 );
//  block.setBlock( 0, 0, &pvp );
//  block.setBlock( 0, 1, &pkq, false, -1.0 );
//  block.setBlock( 1, 0, &pkq, true );
//  block.setBlock( 1, 1, &qdq );
//  if ( !block.isValid( ) ) std::cout << "block invalid!" << std::endl;
//
//  BlockLinearOperator< LO, SC > rhsblock( 2, 2 );
//  rhsblock.setBlock( 0, 0, &pvr, false, -1.0 );
//  rhsblock.setBlock( 0, 1, &pmks );
//  rhsblock.setBlock( 1, 0, &rmkq, true );
//  rhsblock.setBlock( 1, 1, &qds, false, -1.0 );
//  if ( !rhsblock.isValid( ) ) std::cout << "rhsblock invalid!" << std::endl;
//
//  Vector< LO, SC > h( neuElems.size( ) );
//  Vector< LO, SC > g( dirNodes.size( ) );
//  R.apply( neu, h, true );
//  S.apply( dir, g, true );
//  Vector< LO, SC > hg( h );
//  hg.append( g );
//
//  Vector< LO, SC > rhs( rhsblock.getDimRange( ) );
//  rhsblock.apply( hg, rhs );
//
//  Vector< LO, SC > sol( block.getDimRange( ) );
//  block.GMRESSolve( rhs, sol, GMRESprec, GMRESmaxit, GMRESmaxit );
//
//  Vector< LO, SC > s( dirElems.size( ), sol.getData( ) );
//  Vector< LO, SC > t( neuNodes.size( ),
//      sol.getData( ) + dirElems.size( ) );
//  P.apply( s, neu, false, 1.0, 1.0 );
//  Q.apply( t, dir, false, 1.0, 1.0 );
//
//  Vector< LO, SC > dir2( 3 * nNodes );
//  Vector< LO, SC > neu2( 3 * nElems );
//  for ( LO i = 0; i < nElems; ++i ) {
//    neu2.set( 3 * i, neu.get( i ) );
//    neu2.set( 3 * i + 1, neu.get( i + nElems ) );
//    neu2.set( 3 * i + 2, neu.get( i + 2 * nElems ) );
//  }
//  for ( LO i = 0; i < nNodes; ++i ) {
//    dir2.set( 3 * i, dir.get( i ) );
//    dir2.set( 3 * i + 1, dir.get( i + nNodes ) );
//    dir2.set( 3 * i + 2, dir.get( i + 2 * nNodes ) );
//  }
//
//  Vector< LO, SC > bc_nodes( 3 * nNodes );
//  bc_nodes.setAll( 0.0 );
//  for ( auto it = neuNodes.begin( ); it != neuNodes.end( ); ++it ) {
//    bc_nodes.set( 3 * ( *it % nNodes ) + *it / nNodes, 1.0 );
//  }
//  Vector< LO, SC > bc_elems( 3 * nElems );
//  bc_elems.setAll( 0.0 );
//  for ( auto it = neuElems.begin( ); it != neuElems.end( ); ++it ) {
//    bc_elems.set( 3 * ( *it % nElems ) + *it / nElems, 1.0 );
//  }
//
//  std::vector< std::string > vnName{ "bc", "displacement" };
//  std::vector< Vector< LO, SC > * > vnData{ &bc_nodes, &dir2 };
//  std::vector< std::string > veName{ "bc", "stress" };
//  std::vector< Vector< LO, SC > * > veData{ &bc_elems, &neu2 };
//
//  mesh.printParaviewVtu( "output/output.vtu", nullptr, nullptr, nullptr,
//      nullptr, &vnName, &vnData, &veName, &veData );
//
//  delete V;
//  delete K;
//  delete D;
//  if ( quadFar ) delete [] quadFar;
//}