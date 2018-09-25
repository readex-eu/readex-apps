#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../SurfaceMesh3D.h"
#include "../../Laplace1LayerP0P0MultilvlPrecond.h"
#include "../../FullMatrix.h"
#include "../../Vector.h"

using namespace std;
using namespace bem4i;

void aux( );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  aux( );

  return 0;
}

void aux( ) {

  typedef int LO;
  typedef double SC;
  //typedef typename GetType<LO, SC>::SCVT SCVT;

  LO n = 1000;
  int nmult = 100;
  FullMatrix< LO, SC > A( n, n );
  Vector< LO, SC > v1( n );
  Vector< LO, SC > v2( n );
  Vector< LO, SC > v3( n );

  for ( LO j = 0; j < n; ++j ) {
    //v1.set( j, SC( ( j * 5 + 7 ) % 9, ( j * 7 + 9 ) % 5 ) );
    v1.set( j, ( j * 5 + 7 ) % 9 );
    for ( LO i = 0; i < n; ++i ) {
      //A.set( i, j, SC( ( i * 5 + 9 ) % 7, ( j * 3 + 11 ) % 7 ) );
      A.set( i, j, ( i * 5 + 9 ) % 7 );
    }
  }
/*
  v1.print( );
  v1.xferToMIC( );
  v1.xferToHost( );
  v1.print( );                                                                                               
  return;
*/
  ProgressMonitor::init( "apply on cpu" );
  for ( int i = 0; i < nmult; ++i )
    A.apply( v1, v2 );
  ProgressMonitor::step( );

  ProgressMonitor::init( "xfer to mic" );
  for ( int i = 0; i < 2; ++i ) {
    A.xferToMIC( i );
    v1.xferToMIC( i );
    v3.xferToMIC( i );
  }
  ProgressMonitor::step( );

  ProgressMonitor::init( "apply on mic" );
  for ( int i = 0; i < nmult; ++i )
    A.applyMIC( v1, v3, false, 1.0, 0.0, 0 );
  ProgressMonitor::step( );

  ProgressMonitor::init( "xfer to host" );
  v3.xferToHost( );
  ProgressMonitor::step( );

  v3.add( v2, -1.0 );
  std::cout << "error: " << v3.norm2( ) << std::endl;

  ProgressMonitor::init( "apply on mic" );
  for ( int i = 0; i < nmult; ++i )
    A.applyMIC( v1, v3, false, 1.0, 0.0, 1 );
  ProgressMonitor::step( );

  ProgressMonitor::init( "xfer to host" );
  v3.xferToHost( );
  ProgressMonitor::step( );

  v3.add( v2, -1.0 );
  std::cout << "error: " << v3.norm2( ) << std::endl;

  ProgressMonitor::init( "delete on mic" );
  A.deleteMIC( );
  v1.deleteMIC( );
  v3.deleteMIC( );
  ProgressMonitor::step( );
}
