#include "../../Settings.h"

#include <iostream>

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

using namespace std;
using namespace bem4i;

void testLameSteklovPoincare(
    string const & filename,
    int refine4,
    int refine9,
    bool mapToUnitBall,
    int orderNear,
    int orderFar
    );

namespace bem4i {
template<class LO, class SC>
void getLameSteklovPoincare(
    SC * Sarray,
    LO nNodes,
    const SC * nodes,
    LO nElems,
    const LO * elems,
    SC nu,
    SC E,
    LO orderNear,
    LO orderFar,
    bool verbose = false
    );
}

int main(
    int argc,
    char** argv
    ) {
  intro( );

  string filename = "input/cube_12.txt";
  testLameSteklovPoincare( filename, 1, 0, true, 3, 4 );

  return 0;
}

void testLameSteklovPoincare(
    const string & filename,
    int refine4,
    int refine9,
    bool mapToUnitBall,
    int orderNear,
    int orderFar
    ) {

  typedef double SC;
  typedef int LO;

  std::cout.precision( 8 );
  std::cout.setf( std::ios::scientific );

  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );

  mesh.refine( refine4, 2 );
  mesh.refine( refine9, 3 );
  if ( mapToUnitBall ) mesh.mapToUnitBall( );
  //mesh.printParaviewVtk("mesh_alex.vtu");
  LO nNodes = mesh.getNNodes( );
  LO nElems = mesh.getNElements( );

  LO rows = 3 * nNodes;
  SC * Sarray = new SC[ rows * rows ];
  SC nu = 0.33;
  SC E = 1.1e5;

  bem4i::getLameSteklovPoincare<LO, SC>( Sarray, nNodes,
      mesh.getNodes( )->data( ), nElems, mesh.getElements( )->data( ), nu, E,
      orderNear, orderFar, true );

  ///*
  std::ofstream file( "output/sp_lame.txt" );
  file.precision( 8 );
  file.setf( std::ios::scientific );
  FullMatrix< LO, SC > S( rows, rows, Sarray );
  S.print( file );
  //*/
  //Vector<LO, SC> rhs(456);
  //rhs.load("input/rhs_alex.txt");
  //S.LUSolve(rhs);
  //rhs.print();

  delete [] Sarray;
}

namespace bem4i {

//__attribute__( ( visibility( "default" ) ) )

template<class LO, class SC>
void getLameSteklovPoincare(
    SC * Sarray,
    LO nNodes,
    const SC * nodes,
    LO nElems,
    const LO * elems,
    SC nu,
    SC E,
    LO orderNear,
    LO orderFar,
    bool verbose
    ) {

  typedef typename GetType<LO, SC>::SCVT SCVT;

  if ( orderNear < 2 ) orderNear = 2;
  if ( orderNear > 4 ) orderNear = 4;
  if ( orderFar < 3 ) orderFar = 3;
  if ( orderFar > 5 ) orderFar = 5;

  bool symmetricV = false;
#if N_MIC > 0
  symmetricV = false;
#endif

  std::vector< SCVT > nodesv;
  nodesv.assign( nodes, nodes + 3 * nNodes );
  std::vector< LO > elemsv;
  elemsv.assign( elems, elems + 3 * nElems );
  SurfaceMesh3D< LO, SC > mesh( nodesv, elemsv );
  if ( verbose ) mesh.printInfo( );

  quadratureType quadType = SauterSchwab;
  int quadNear[ 4 ] = { orderNear, orderNear, orderNear, orderNear };
  int quadFar[] = { orderFar, orderFar };
  if ( verbose ) std::cout << "Sauter-Schwab near-field order " << orderNear <<
      ", Gauss far-field order " << orderFar << std::endl;

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
  BEBilinearFormLame1Layer< LO, SC > formV( &bespace00, quadNear,
      quadType, quadFar, symmetricV );
  formV.setE( E );
  formV.setNu( nu );
#if N_MIC < 1
  if ( verbose ) ProgressMonitor::init( "V, Vlap" );
  formV.assemble( *V, *Vlap );
  if ( verbose ) ProgressMonitor::step( );
#endif

  FullMatrix< LO, SC > * Klap = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLaplace2Layer< LO, SC > formKlap( &bespace10, quadNear,
      quadType, quadFar );
#if N_MIC < 1
  if ( verbose ) ProgressMonitor::init( "Klap" );
  formKlap.assemble( *Klap );
#else
  if ( verbose ) ProgressMonitor::init( "V, VLap, Klap" );
  formV.assembleAllMIC( *V, *Vlap, *Klap, bespace10 );
#endif
  if ( verbose ) ProgressMonitor::step( );

  FullMatrix< LO, SC > * K = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLame2Layer< LO, SC > formK( &bespace10, quadNear,
      quadType, quadFar );
  formK.setE( E );
  formK.setNu( nu );
  if ( verbose ) ProgressMonitor::init( "K" );
  formK.assemble( *K, *Vlap, *V, *Klap, Tv );
  if ( verbose ) ProgressMonitor::step( );

  delete Klap;

  FullMatrix< LO, SC > * D = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLameHypersingular< LO, SC > formD( &bespace11, quadNear,
      quadType, quadFar );
  formD.setE( E );
  formD.setNu( nu );
  if ( verbose ) ProgressMonitor::init( "D" );
  formD.assemble( *D, *Vlap, *V, Tv );
  if ( verbose ) ProgressMonitor::step( );

  delete Vlap;

  IdentityOperator< LO, SC > id( &bespace10 );
  SparseMatrix< LO, SC > M;
  id.assemble( M );

  Eigen::SparseMatrix<SC, Eigen::ColMajor, LO> * Me = M.getEigenSparseMatrix( );

  for ( LO j = 0; j < Me->outerSize( ); ++j ) {
    typename Eigen::SparseMatrix<SC, Eigen::ColMajor, LO>::InnerIterator
    it( *Me, j );
    for (; it; ++it ) {
      K->add( it.row( ), it.col( ), 0.5 * it.value( ) );
      K->add( it.row( ) + nElems, it.col( ) + nNodes, 0.5 * it.value( ) );
      K->add( it.row( ) + 2 * nElems, it.col( ) + 2 * nNodes,
          0.5 * it.value( ) );
    }
  }

  FullMatrix< LO, SC > * VinvK = new FullMatrix< LO, SC >( *K );
  if ( verbose ) ProgressMonitor::init( "Applying inverse of V" );
  V->CholeskiSolve( *VinvK );
  if ( verbose ) ProgressMonitor::step( );

  delete V;

  FullMatrix< LO, SC > S( 3 * nNodes, 3 * nNodes, Sarray );
  S.setAll( 0.0 );

  S.multiply( *K, *VinvK, true, false, 1.0, 0.0 );

  delete VinvK;
  delete K;

  S.add( *D, 1.0 );

  delete D;
}

}