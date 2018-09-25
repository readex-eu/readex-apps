#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../ProgressMonitor.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../FastBESpace.h"
#include "../../Vector.h"
#include "../../BEIntegratorLaplace.h"
#include "../../BEBilinearFormLaplace1Layer.h"
#include "../../BEBilinearFormLaplace2Layer.h"
#include "../../BEBilinearFormLaplaceHypersingular.h"
#include "../../IdentityOperator.h"
#include "../../RepresentationFormulaLaplace.h"

using namespace std;
using namespace bem4i;

void testLaplaceNeumann(
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

  string filename = "input/icosahedron.txt";
  testLaplaceNeumann( filename, 3, 0, true, 0, 5, 1 );

  return 0;
}

void testLaplaceNeumann(
  string const &filename,
  int refine4,
  int refine9,
  bool mapToUnitBall,
  int quadTypeInt,
  int order,
  int nPoints
  ) {

  typedef double SC;
  typedef long LO;
  std::cout.precision( 8 );
  std::cout.setf( std::ios::scientific );

  ProgressMonitor::init( "Loading mesh" );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  mesh.refine( refine4, 2 );
  mesh.refine( refine9, 3 );
  if ( mapToUnitBall ) mesh.mapToUnitBall( );
  ProgressMonitor::step( );
  mesh.printInfo( );

  SC CGeps = 1e-12;

  int quadOrder[ 4 ];
  quadratureType quadType;
  if ( quadTypeInt == 0 ) {
    quadType = Steinbach;
    quadOrder[ 0 ] = quadOrder[ 1 ] = order;
    std::cout << "Using Steinbach quadrature, order " << order << "."
      << std::endl;
  } else {
    quadType = SauterSchwab;
    quadOrder[ 0 ] = quadOrder[ 1 ] = quadOrder[ 2 ] = quadOrder[ 3 ] = order;
    std::cout << "Using Sauter-Schwab quadrature, order " << order << "."
      << std::endl;
  }

  //int qOrdDisjoint[ 2 ] = { 4, 4 };
  int * qOrdDisjoint = nullptr;
  /*

  // create a cluster tree from mesh
  Tree<BECluster<LO, SC>*, LO> tree, tree2;
  ProgressMonitor::init( "Nested dissection" );
  mesh.nestedDissection( tree, 30 );
  mesh.nestedDissection( tree2, 30 );
  ProgressMonitor::step( );

  // create a space for V, assemble V using bilinear form
  FastBESpace< LO, SC > bespaceV( &mesh, p0, p0, &tree, 1.2, 3, 5, 0, false );

  BESpace< LO, SC > bespaceV( &mesh, p0, p0 );
  bespaceV.setEpsilonACA( 1e-4 );
  ACAMatrix< LO, SC > V( 0, 0 );
  BEBilinearFormLaplace1Layer< LO, SC > formV( &bespaceV, quadOrder, quadType,
      qOrdDisjoint );
  formV.assemble( V );
  //V.print();
   */

  BESpace< LO, SC > bespaceD( &mesh, p1, p1 );
  FullMatrix< LO, SC > D( 0, 0 );
  //LaplaceHypersingularOperator< LO, SC > D( &bespaceD, &V );
  //D.setRegularized( );
  BEBilinearFormLaplaceHypersingular< LO, SC > formD( &bespaceD, quadOrder,
    quadType, qOrdDisjoint );
  ProgressMonitor::init( "Assembling D" );
  //formD.assemble( D, V );
  formD.assemble( D );
  ProgressMonitor::step( );
  //D.print( );

  ///*
  //FastBESpace< LO, SC > bespaceK( &mesh, p1, p0, &tree2, 1.2, 8, 5, 0, false );
  //ACAMatrix< LO, SC > K;
  //bespaceK.setEpsilonACA( 1e-4 );
  BESpace< LO, SC > bespaceK( &mesh, p1, p0 );
  FullMatrix< LO, SC > K( 0, 0 );
  BEBilinearFormLaplace2Layer< LO, SC > formK( &bespaceK, quadOrder, quadType,
    qOrdDisjoint );
  ProgressMonitor::init( "Assembling K" );
  formK.assemble( K );
  ProgressMonitor::step( );
  //K.print();

  IdentityOperator< LO, SC > id( &bespaceK );

  LO nNodes = mesh.getNNodes( );
  LO nElems = mesh.getNElements( );
  Vector< LO, SC > dir( nNodes );
  Vector< LO, SC > rhs( nNodes );
  Vector< LO, SC > aux( nNodes );
  Vector< LO, SC > neu( nElems );
  SC x1[ 3 ], x2[ 3 ], x3[ 3 ], y[ 3 ], n[ 3 ];

  for ( LO i = 0; i < nNodes; i++ ) {
    mesh.getNode( i, x1 );
    //dir.set( i, x1[ 0 ] * x1[ 1 ] * x1[ 2 ] );
    dir.set( i, ( 1.0 + x1[0] ) * std::exp( 2.0 * M_PI * x1[1] ) *
      std::cos( 2.0 * M_PI * x1[2] ) );
  }
  //dir.print();

  for ( LO i = 0; i < nElems; i++ ) {
    mesh.getNodes( i, x1, x2, x3 );
    y[0] = 1.0 / 3.0 * ( x1[0] + x2[0] + x3[0] );
    y[1] = 1.0 / 3.0 * ( x1[1] + x2[1] + x3[1] );
    y[2] = 1.0 / 3.0 * ( x1[2] + x2[2] + x3[2] );
    //neu.set( i,  1.0/9.0 * (x1[0]+x2[0]+x3[0]) * (x1[1]+x2[1]+x3[1])
    //  * (x1[2]+x2[2]+x3[2]) );
    mesh.getNormal( i, n );
    neu.set( i, std::exp( 2.0 * M_PI * y[1] )
      * ( n[0] * std::cos( 2.0 * M_PI * y[2] )
      + 2.0 * M_PI * ( 1.0 + y[0] ) * n[1] * std::cos( 2.0 * M_PI * y[2] )
      - 2.0 * M_PI * ( 1.0 + y[0] ) * n[2] * std::sin( 2.0 * M_PI * y[2] ) ) );
  }
  //neu.print();

  ///* // for hypersingular equation
  ProgressMonitor::init( "Setting up RHS" );
  id.apply( neu, rhs, true, 0.5, 0.0 );
  K.apply( neu, aux, true );
  rhs.add( aux, -1.0 );
  ProgressMonitor::step( );
  // */

  /* // for steklov poincare equation
  ProgressMonitor::init( "Setting up RHS" );
  id.apply( neu, rhs, true );
  ProgressMonitor::step( );
   */

  ProgressMonitor::init( "Stabilizing" );
  Vector<LO, SC> a( nNodes, true );
  LO idx[ 3 ];
  SC areaW;
  for ( LO i = 0; i < nElems; i++ ) {
    mesh.getElement( i, idx );
    areaW = mesh.getElemArea( i ) / 3.0;
    a.add( idx[ 0 ], areaW );
    a.add( idx[ 1 ], areaW );
    a.add( idx[ 2 ], areaW );
  }

  SC integral = 0.0;
  SC tmp;
  LO ind[3];
  for ( LO i = 0; i < nElems; i++ ) {
    tmp = 0.0;
    mesh.getElement( i, ind );
    for ( int j = 0; j < 7; j++ ) {
      tmp += quadWeights5[ j ] * (
        ( dir.get( ind[0] ) ) * ( 1.0 - quadPoints5[ 2 * j ]
        - quadPoints5[ 2 * j + 1 ] ) +
        ( dir.get( ind[1] ) ) * quadPoints5[ 2 * j ] +
        ( dir.get( ind[2] ) ) * quadPoints5[ 2 * j + 1 ]
        );
    }
    integral += tmp * mesh.getElemArea( i );
  }
  std::cout << "(Solution integral: " << integral << ".) ";
  rhs.add( a, integral );

  ///*
  for ( LO i = 0; i < a.getLength( ); ++i ) {
    for ( LO j = 0; j < a.getLength( ); ++j ) {
      D.add( i, j, a.get( i ) * a.get( j ) );
    }
  }
  //*/
  ProgressMonitor::step( );


  ProgressMonitor::init( "Solving the system" );
  Vector<LO, SC> rhs2( rhs );
  rhs.setAll( 0.0 );
  D.CGSolve( rhs2, rhs, CGeps );
  //LaplaceSteklovPoincareOperator< LO, SC > DtN( &bespaceD, &V, &K, &D, &id );
  //DtN.CGSolve( rhs2, rhs, CGeps );
  ProgressMonitor::step( );
  //rhs.print();

  std::cout << "L2 relative error: "
    << mesh.l2RelativeErrorLin( rhs ) << "." << std::endl;

  string nodeNames[] = { "Dirichlet_anal", "Dirichlet_comp", "Dirichlet_err" };
  string elemNames[] = { "Neumann_anal" };
  Vector< LO, SC > err( rhs );
  rhs.add( dir, err, -1.0 );
  Vector< LO, SC >* nodalData[] = { &dir, &rhs, &err };
  Vector< LO, SC >* elemData[] = { &neu };
  mesh.printParaviewVtu( "output/output.vtu", 3, nodeNames, nodalData, 1, 
    elemNames, elemData );


  SC * evalPoint = new SC[ 3 * nPoints ];
  for ( LO i = 0; i < nPoints; i++ ) {

    evalPoint[3 * i] = 0.250685;
    evalPoint[3 * i + 1] = 0.417808;
    evalPoint[3 * i + 2] = 0.584932;
  }
  SC exact = ( 1.0 + evalPoint[ 0 ] )
    * std::exp( 2.0 * M_PI * evalPoint[ 1 ] )
    * std::cos( 2.0 * M_PI * evalPoint[ 2 ] );
  Vector< LO, SC > res( nPoints );

  RepresentationFormulaLaplace<LO, SC> formula( &bespaceK, &rhs, &neu );
  ProgressMonitor::init( "Evaluating representation formula" );
  formula.evaluate( evalPoint, nPoints, true, res );
  ProgressMonitor::step( );
  std::cout << "Point evaluation in ["
    << evalPoint[ 0 ] << ", "
    << evalPoint[ 1 ] << ", "
    << evalPoint[ 2 ] << "]: "
    << res.get( 0 ) << "." << std::endl;

  std::cout << "Absolute error: "
    << std::abs( res.get( 0 ) - exact ) << std::endl;

  std::cout << "Relative error: "
    << std::abs( ( res.get( 0 ) - exact ) / exact ) << std::endl;
  //*/
}
