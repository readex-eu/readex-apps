#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../BESpace.h"
#include "../../BEBilinearFormHelmholtz1Layer.h"
#include "../../BEBilinearFormHelmholtz2Layer.h"
#include "../../BEBilinearFormHelmholtzHypersingular.h"
#include "../../IdentityOperator.h"
#include "../../BEBilinearFormLaplace1Layer.h"
#include "../../BEBilinearFormLaplace2Layer.h"

using namespace std;
using namespace bem4i;

void testMatrixAssemblyTime(
    string const &filename,
    int refine4,
    int refine9,
    bool mapToUnitBall,
    int quadType,
    int order
    );

int main(
    int argc,
    char** argv
    ) {
  if ( argc == 8 ) {
    omp_set_num_threads( atoi( argv[ 1 ] ) );

    intro( );

    string filename = string( argv[ 2 ] );
    int refine4 = atoi( argv[ 3 ] );
    int refine9 = atoi( argv[ 4 ] );
    bool mapToUnitBall = atoi( argv[ 5 ] );
    int quadType = atoi( argv[ 6 ] );
    int order = atoi( argv[ 7 ] );
    testMatrixAssemblyTime( filename, refine4, refine9, mapToUnitBall, quadType,
        order );
  }

  return 0;
}

void testMatrixAssemblyTime( string const &filename, int refine4, int refine9, bool mapToUnitBall, int quadTypeInt, int order ) {

  typedef complex<double> SC;
  typedef int LO;
  typedef typename SC::value_type SCVT;
  timeval start, stop;

  std::cout.precision( 8 );
  std::cout.setf( std::ios::scientific );

  SurfaceMesh3D< LO, SCVT > mesh;
  mesh.load( filename.c_str( ) );
  mesh.refine( refine4, 2 );
  mesh.refine( refine9, 3 );
  if ( mapToUnitBall ) {
    mesh.mapToUnitBall( );
  }
  mesh.printInfo( );
  int quadOrder[ 4 ];
  quadratureType quadType;
  if ( quadTypeInt == 0 ) {
    quadType = Steinbach;
    quadOrder[ 0 ] = quadOrder[ 1 ] = order;
    std::cout << "Using Steinbach quadrature, order " << order << "." << std::endl;
  } else {
    quadType = SauterSchwab;
    quadOrder[ 0 ] = quadOrder[ 1 ] = quadOrder[ 2 ] = quadOrder[ 3 ] = order;
    std::cout << "Using Sauter-Schwab quadrature, order " << order << "." << std::endl;
  }
  ///*
  BESpace< LO, SCVT > bespaceV( &mesh, p0, p0 );
  FullMatrix< LO, SCVT > V( 0, 0 );
  BEBilinearFormLaplace1Layer< LO, SCVT > formV( &bespaceV, quadOrder, quadType );
  std::cout << "Assembling V ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  formV.assemble( V );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  BESpace< LO, SCVT > bespaceK( &mesh, p1, p0 );
  FullMatrix< LO, SCVT > K( 0, 0 );
  BEBilinearFormLaplace2Layer< LO, SCVT > formK( &bespaceK, quadOrder, quadType );
  std::cout << "Assembling K ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  formK.assemble( K );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  /*
  BESpace< LO, SCVT > bespaceD( &mesh, p1, p1 );
  FullMatrix< LO, SCVT > D( 0, 0 );
  BEBilinearFormLaplaceHypersingular< LO, SCVT > formD( &bespaceD );
  std::cout << "Assembling D ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  formD.assemble( D );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;
   */

  SurfaceMesh3D< LO, SC > meshSCVT;
  meshSCVT.load( filename.c_str( ) );
  meshSCVT.refine( refine4, 2 );
  meshSCVT.refine( refine9, 3 );
  if ( mapToUnitBall ) {

    meshSCVT.mapToUnitBall( );
  }
  meshSCVT.printInfo( );

  SC kappa = 2.0;
  ///*
  BESpace< LO, SC > bespaceVkappa( &meshSCVT, p0, p0 );
  FullMatrix< LO, SC > Vkappa( 0, 0 );
  BEBilinearFormHelmholtz1Layer< LO, SC > formVkappa( &bespaceVkappa, quadOrder, kappa, quadType );
  std::cout << "Assembling V kappa ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  formVkappa.assemble( Vkappa );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  BESpace< LO, SC > bespaceKkappa( &meshSCVT, p1, p0 );
  FullMatrix< LO, SC > Kkappa( 0, 0 );
  BEBilinearFormHelmholtz2Layer< LO, SC > formKkappa( &bespaceKkappa, quadOrder, kappa, quadType );
  std::cout << "Assembling K kappa ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  formKkappa.assemble( Kkappa );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;
  // */
  BESpace< LO, SC > bespaceDkappa( &meshSCVT, p1, p1 );
  FullMatrix< LO, SC > Dkappa( 0, 0 );
  BEBilinearFormHelmholtzHypersingular< LO, SC > formDkappa( &bespaceDkappa, quadOrder, kappa, quadType );
  std::cout << "Assembling D kappa ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  formDkappa.assemble( Dkappa );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;
}