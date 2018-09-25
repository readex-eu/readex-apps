#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../ProgressMonitor.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../BESpace.h"
#include "../../Vector.h"
#include "../../BEIntegratorLaplace.h"
#include "../../BEBilinearFormLaplace1Layer.h"
#include "../../BEBilinearFormLaplace2Layer.h"
#include "../../IdentityOperator.h"
#include "../../SparseMatrix.h"

using namespace std;
using namespace bem4i;

void testLaplaceSteklovPoincare(
    string const & filename,
    int refine4,
    int refine9,
    bool mapToUnitBall,
    int quadType,
    int orderNear,
    int orderFar
    );

template< class LO, class SC >
void assembleT(
    SurfaceMesh3D< LO, SC > & mesh,
    SparseMatrix< LO, SC > & T1,
    SparseMatrix< LO, SC > & T2,
    SparseMatrix< LO, SC > & T3
    );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  string filename = "input/cube_12.txt";
  testLaplaceSteklovPoincare( filename, 3, 0, true, 1, 3, 4 );

  return 0;
}

void testLaplaceSteklovPoincare(
    const string & filename,
    int refine4,
    int refine9,
    bool mapToUnitBall,
    int quadTypeInt,
    int orderNear,
    int orderFar
    ) {

  typedef double SC;
  typedef int LO;

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

  int quadOrder[ 4 ];
  quadratureType quadType;
  if ( quadTypeInt == 0 ) {
    quadType = Steinbach;
    quadOrder[ 0 ] = quadOrder[ 1 ] = orderNear;
    std::cout << "Using Steinbach quadrature, order " << orderNear << "."
        << std::endl;
  } else {
    quadType = SauterSchwab;
    quadOrder[ 0 ] = quadOrder[ 1 ] = quadOrder[ 2 ] = quadOrder[ 3 ] =
        orderNear;
    std::cout << "Using Sauter-Schwab quadrature, order " << orderNear << "."
        << std::endl;
  }

  int quadDisjoint[] = { orderFar, orderFar };

  LO nNodes = mesh.getNNodes( );

  BESpace< LO, SC > bespace00( &mesh, p0, p0 );
  BESpace< LO, SC > bespace10( &mesh, p1, p0 );
  BESpace< LO, SC > bespace11( &mesh, p1, p1 );

  ProgressMonitor::init( "V" );
  FullMatrix< LO, SC > * V = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLaplace1Layer< LO, SC > formV( &bespace00, quadOrder,
      quadType, quadDisjoint );
  formV.assemble( *V );
  ProgressMonitor::step( );

  //V->print();

  ProgressMonitor::init( "K" );
  FullMatrix< LO, SC > * K = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLaplace2Layer< LO, SC > formK( &bespace10, quadOrder,
      quadType, quadDisjoint );
  formK.assemble( *K );
  ProgressMonitor::step( );

  //K->print();

  //  FullMatrix< LO, SC > * D = new FullMatrix< LO, SC >( 0, 0 );
  //  BEBilinearFormLaplaceHypersingular< LO, SC > formD( &bespace11, quadOrder,
  //      quadType, quadDisjoint );
  //  formD.assemble( *D, *V );

  ProgressMonitor::init( "D" );
  SparseMatrix< LO, SC > T1, T2, T3;
  assembleT< LO, SC >( mesh, T1, T2, T3 );
  FullMatrix< LO, SC > * VT1, * VT2, * VT3;
  LO nElems = mesh.getNElements( );

  FullMatrix< LO, SC > * D = new FullMatrix< LO, SC >( nNodes, nNodes );

  VT1 = new FullMatrix< LO, SC >( nElems, nNodes );
  VT1->multiply( *V, T1 );
  D->multiply( T1, *VT1, true );
  delete VT1;

  VT2 = new FullMatrix< LO, SC >( nElems, nNodes );
  VT2->multiply( *V, T2 );
  D->multiply( T2, *VT2, true, false, 1.0, 1.0 );
  delete VT2;

  VT3 = new FullMatrix< LO, SC >( nElems, nNodes );
  VT3->multiply( *V, T3 );
  D->multiply( T3, *VT3, true, false, 1.0, 1.0 );
  delete VT3;

  ProgressMonitor::step( );

  //D->print( );

  IdentityOperator< LO, SC > id( &bespace10 );
  SparseMatrix< LO, SC > M;
  id.assemble( M );

  //M.print();

  K->add( M, 0.5 );
  FullMatrix< LO, SC > * VinvK = new FullMatrix< LO, SC >( *K );
  ProgressMonitor::init( "Applying inverse of V" );
  V->CholeskiSolve( *VinvK );
  ProgressMonitor::step( );

  delete V;

  FullMatrix< LO, SC > S( nNodes, nNodes, false );
  S.multiply( *K, *VinvK, true, false, 1.0, 0.0 );

  delete VinvK;
  delete K;

  S.add( *D, 1.0 );

  //S.print();

  delete D;
}

template< class LO, class SC >
void assembleT(
    SurfaceMesh3D< LO, SC > & mesh,
    SparseMatrix< LO, SC > & T1,
    SparseMatrix< LO, SC > & T2,
    SparseMatrix< LO, SC > & T3
    ) {

  LO nElems = mesh.getNElements( );
  LO nNodes = mesh.getNNodes( );

  std::vector< LO > rowInd;
  rowInd.reserve( 3 * nElems );
  std::vector< LO > colInd;
  colInd.reserve( 3 * nElems );
  std::vector< SC > values1;
  values1.reserve( 3 * nElems );
  std::vector< SC > values2;
  values2.reserve( 3 * nElems );
  std::vector< SC > values3;
  values3.reserve( 3 * nElems );
  LO elem[ 3 ];
  Vector< LO, SC > * curls = mesh.getCurls( );

  for ( LO i = 0; i < nElems; ++i ) {

    mesh.getElement( i, elem );

    for ( int node = 0; node < 3; ++node ) {

      rowInd.push_back( i );
      colInd.push_back( elem[ node ] );
      values1.push_back( curls->get( 9 * i + 3 * node ) );
      values2.push_back( curls->get( 9 * i + 3 * node + 1 ) );
      values3.push_back( curls->get( 9 * i + 3 * node + 2 ) );
    }
  }

  T1.setFromTriplets( nElems, nNodes, rowInd, colInd, values1 );
  T2.setFromTriplets( nElems, nNodes, rowInd, colInd, values2 );
  T3.setFromTriplets( nElems, nNodes, rowInd, colInd, values3 );

}