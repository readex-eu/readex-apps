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
#include "../../IdentityOperator.h"

#include "readex.h"

using namespace std;
using namespace bem4i;

void testLaplaceDirichletACA(
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
  testLaplaceDirichletACA( filename, 5, 0, true, 0, 5, 1 );

  return 0;
}

void testLaplaceDirichletACA(
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

  READEX_INIT()
  READEX_REGION_DEFINE(main)
  READEX_PHASE_START(main,"main",SCOREP_USER_REGION_TYPE_COMMON);

  std::cout.precision( 8 );
  std::cout.setf( std::ios::scientific );

  ProgressMonitor::init( "Initializing mesh" );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );

  mesh.refine( refine4, 2 );
  mesh.refine( refine9, 3 );
  if ( mapToUnitBall ) mesh.mapToUnitBall( );
  mesh.printInfo( );
  ProgressMonitor::step( );

  LO nNodes = mesh.getNNodes( );
  LO nElems = mesh.getNElements( );

  const SC CGeps = 1e-12;

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

  //int quadDisj [] = { 4, 4 };
  int * quadDisj = nullptr;

  SC ACAeps = 1e-4;
  SC ACAscale = 1e-3;
  SC ACAeta = 1.2;
  //  bool ACAgroup = true;
  LO ACAgroupSize = 500;
  LO dummy = 0;
  LO maxElems = 50;

  // create a cluster tree from mesh
  Tree<BECluster<LO, SC>*, LO> tree, tree2;
  ProgressMonitor::init( "Nested dissection" );
  mesh.nestedDissection( tree, maxElems );
  mesh.nestedDissection( tree2, maxElems );
  ProgressMonitor::step( );

  READEX_REGION_DEFINE(assemble_v)
  READEX_REGION_START(assemble_v,"assemble_v",SCOREP_USER_REGION_TYPE_COMMON);
  
  §// create a space for V, assemble V using bilinear form
  FastBESpace< LO, SC > bespaceV( &mesh, p0, p0, &tree, ACAeta, dummy, dummy,
      dummy, ACAgroupSize );

  bespaceV.setEpsilonACA( ACAeps );
  bespaceV.setScaleACA( ACAscale );
  ACAMatrix< LO, SC > V;
  BEBilinearFormLaplace1Layer< LO, SC > formV( &bespaceV, quadOrder,
      quadType, quadDisj );
  ProgressMonitor::init( "Assembling V" );
  formV.assemble( V );
  ProgressMonitor::step( );
  //return;

  READEX_REGION_STOP(assemble_v) 
  
  std::cout << "compression: " << V.getCompressionRatio( ) << std::endl;
  std::cout << "row, col: " << formV.getRowCount( ) << ", " <<
      formV.getColCount( ) << std::endl;

  READEX_REGION_DEFINE(assemble_k)
  READEX_REGION_START(assemble_k,"assemble_k",SCOREP_USER_REGION_TYPE_COMMON);
  
  // create a space for K, assemble K using bilinear form
  FastBESpace< LO, SC > bespaceK( &mesh, p1, p0, &tree, ACAeta, dummy, dummy,
      dummy, ACAgroupSize );
  ACAMatrix< LO, SC > K;
  bespaceK.setEpsilonACA( ACAeps );
  bespaceK.setScaleACA( ACAscale );
  BEBilinearFormLaplace2Layer< LO, SC > formK( &bespaceK, quadOrder,
      quadType, quadDisj );
  ProgressMonitor::init( "Assembling K" );
  formK.assemble( K );
  ProgressMonitor::step( );
  //return;

  READEX_REGION_STOP(assemble_k) 
  
  std::cout << "compression: " << K.getCompressionRatio( ) << std::endl;
  std::cout << "row, col: " << formK.getRowCount( ) << ", " <<
      formK.getColCount( ) << std::endl;

  /*
  Vector< LO, SC > vin( nNodes );
  vin.setAll( 1.0 );
  Vector< LO, SC > vout( nElems );
  ProgressMonitor::init( "K apply" );
  for( int i = 0; i < 20; ++i )
    K.apply( vin, vout );
  ProgressMonitor::step( );
  return;
  */

  BESpace< LO, SC > bespace10( &mesh, p1, p0 );
  IdentityOperator< LO, SC > id( &bespace10 );

  Vector< LO, SC > dir( nNodes );
  Vector< LO, SC > rhs( nElems );
  Vector< LO, SC > aux( nElems );
  Vector< LO, SC > neu( nElems );
  SC x1[ 3 ], x2[ 3 ], x3[ 3 ], y[ 3 ], n[ 3 ];

  //SC y1[] = { 0.15, 0.24, 0.7 };
  //SC y2[] = { -0.125, 0.23, 0.74 };
  //SC norm1, norm2, dot1, dot2;
  for ( LO i = 0; i < nNodes; i++ ) {
    mesh.getNode( i, x1 );
    //dir.set( i, x1[ 0 ] * x1[ 1 ] * x1[ 2 ] );
    dir.set( i, ( 1.0 + x1[0] ) * std::exp( 2.0 * M_PI * x1[1] )
        * std::cos( 2.0 * M_PI * x1[2] ) );
    //norm1 = std::sqrt( (x1[0]-y1[0])*(x1[0]-y1[0]) + (x1[1]-y1[1])*(x1[1]-y1[1]) + (x1[2]-y1[2])*(x1[2]-y1[2]) );
    //norm2 = std::sqrt( (x1[0]-y2[0])*(x1[0]-y2[0]) + (x1[1]-y2[1])*(x1[1]-y2[1]) + (x1[2]-y2[2])*(x1[2]-y2[2]) );
    //dir.set( i, 1.0 / (4.0*M_PI) * ( 1.0 / norm1 + 1.0 / norm2 ) ); //sum of two fund. sols with y1, y2
    //TEST PRO OPTIMALIZACI
    //        if ( i >= 642 ) {
    //          dir.set( i, 1.0 );
    //        }
  }
  //dir.print();

  for ( LO i = 0; i < nElems; i++ ) {

    mesh.getNodes( i, x1, x2, x3 );
    y[0] = 1.0 / 3.0 * ( x1[0] + x2[0] + x3[0] );
    y[1] = 1.0 / 3.0 * ( x1[1] + x2[1] + x3[1] );
    y[2] = 1.0 / 3.0 * ( x1[2] + x2[2] + x3[2] );
    //neu.set( i,  1.0/9.0 * (x1[0]+x2[0]+x3[0]) * (x1[1]+x2[1]+x3[1]) * (x1[2]+x2[2]+x3[2]) );
    //neu.set( i, std::exp(2.0*M_PI*y[1]) * ( y[0]*std::cos(2.0*M_PI*y[2]) + 2.0*M_PI*(1.0+y[0])*y[1]*std::cos(2.0*M_PI*y[2]) - 2.0*M_PI*(1.0+y[0])*y[2]*std::sin(2.0*M_PI*y[2]) ) );
    mesh.getNormal( i, n );
    //neu.set( i, n[ 0 ] * y[ 1 ] * y[ 2 ] + y[ 0 ] * n[ 1 ] * y[ 2 ] + y[ 0 ] * y[ 1 ] * n[ 2 ] );
    neu.set( i, std::exp( 2.0 * M_PI * y[1] ) * ( n[0] * std::cos( 2.0 * M_PI * y[2] )
        + 2.0 * M_PI * ( 1.0 + y[0] ) * n[1] * std::cos( 2.0 * M_PI * y[2] )
        - 2.0 * M_PI * ( 1.0 + y[0] ) * n[2] * std::sin( 2.0 * M_PI * y[2] ) ) );
    //norm1 = std::sqrt( (y[0]-y1[0])*(y[0]-y1[0]) + (y[1]-y1[1])*(y[1]-y1[1]) + (y[2]-y1[2])*(y[2]-y1[2]) );
    //norm2 = std::sqrt( (y[0]-y2[0])*(y[0]-y2[0]) + (y[1]-y2[1])*(y[1]-y2[1]) + (y[2]-y2[2])*(y[2]-y2[2]) );
    //dot1 = y[0]*(y[0]-y1[0]) + y[1]*(y[1]-y1[1]) + y[2]*(y[2]-y1[2]);
    //dot2 = y[0]*(y[0]-y2[0]) + y[1]*(y[1]-y2[1]) + y[2]*(y[2]-y2[2]);
    //neu.set( i, - 1.0 / (4.0*M_PI) * ( dot1 / (norm1*norm1*norm1) + dot2 / (norm2*norm2*norm2) ) );
    // TEST PRO OPTIMALIZACI
    //    if ( y[ 2 ] <= 1.0 ) {
    //      neu.set( i, -0.5 * n[ 2 ] );
    //    }
  }

  ProgressMonitor::init( "Setting up rhs" );
  id.apply( dir, aux, false, 0.5, 0.0 );
  //aux.print();
  //Vector<LO, SC> dir2(3 * mesh.getNElements());
  //transformator.apply(dir, dir2);


  K.apply( dir, rhs );
  rhs.add( aux );
  //rhs.print();
  ProgressMonitor::step( );

  READEX_REGION_DEFINE(cg_solve)
  READEX_REGION_START(cg_solve,"cg_solve",SCOREP_USER_REGION_TYPE_COMMON);
  
  ProgressMonitor::init( "Solving the system" );
  Vector<LO, SC> rhs2( rhs );
  rhs.setAll( 0.0 );
  V.CGSolve( rhs2, rhs, CGeps );
  ProgressMonitor::step( );

  READEX_REGION_STOP(cg_solve)
  
  //std::cout << "L2 relative error: " << mesh.l2RelativeErrorConst( neu, rhs ) << "." << std::endl;
  std::cout << "L2 relative error: " << mesh.l2RelativeErrorConst( rhs ) << "."
      << std::endl;

  READEX_REGION_DEFINE(print_vtu)
  READEX_REGION_START(print_vtu,"print_vtu",SCOREP_USER_REGION_TYPE_COMMON);
  
  string nodeNames[] = { "Dirichlet_anal" };
  Vector< LO, SC >* nodalData[] = { &dir };
  string elemNames[] = { "Neumann_anal", "Neumann_comp", "Neumann_err" };
  Vector< LO, SC > err( rhs );
  rhs.add( neu, err, -1.0 );
  Vector< LO, SC >* elemData[] = { &neu, &rhs, &err };
  mesh.printParaviewVtu( "output/output.vtu", 1, nodeNames, nodalData, 3, elemNames, elemData );
  
  READEX_REGION_STOP(print_vtu)

  //printData( "input/target_sphere_5120_2.txt", rhs, 0, 5119 );

  //  std::cout << "Evaluating representation formula ... ";
  //  std::cout.flush( );
  //  LO nPoints = 10000;
  //  SC evalPoint[30000];
  //  for ( int i = 0; i < nPoints; i++ ) {
  //    evalPoint[3 * i] = 0.250685;
  //    evalPoint[3 * i + 1] = 0.417808;
  //    evalPoint[3 * i + 2] = 0.584932;
  //  }
  //  SC exact = ( 1.0 + evalPoint[ 0 ] )
  //      * std::exp( 2.0 * M_PI * evalPoint[ 1 ] )
  //      * std::cos( 2.0 * M_PI * evalPoint[ 2 ] );
  //  Vector< LO, SC > res( nPoints );
  //  RepresentationFormulaLaplace<LO, SC> formula( &bespaceK, &dir, &rhs );
  //  gettimeofday( &start, nullptr );
  //  formula.evaluate( evalPoint, nPoints, true, res );
  //  gettimeofday( &stop, nullptr );
  //  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;
  //  std::cout << "Point evaluation in ["
  //      << evalPoint[ 0 ] << ", "
  //      << evalPoint[ 1 ] << ", "
  //      << evalPoint[ 2 ] << "]: "
  //      << res.get( 0 ) << "." << std::endl;
  //
  //  std::cout << "Absolute error: "
  //      << std::abs( res.get( 0 ) - exact ) << std::endl;
  //  //*/
  
  READEX_PHASE_STOP(main)
  READEX_CLOSE()
}

//void printData( string const &filename, Vector< int, double > const &data, int begin, int end ) {
//  int length = data.getLength( );
//  if ( begin <= 0 || begin >= length ) begin = 0;
//  if ( end <= 0 || end >= length ) end = length - 1;
//
//  std::cout << "Printing  " << filename << " ... ";
//
//  std::ofstream file( filename.c_str( ) );
//
//  file.setf( std::ios::scientific );
//  file.precision( 6 );
//
//  if ( !file.is_open( ) ) {
//    std::cout << "File could not be opened!" << std::endl;
//    return;
//  }
//
//  file << end - begin + 1 << std::endl << std::endl;
//
//  for ( int i = begin; i <= end; ++i ) {
//
//    file << data.get( i ) << std::endl;
//  }
//
//  file.close( );
//
//  std::cout << "done." << std::endl;
//}
