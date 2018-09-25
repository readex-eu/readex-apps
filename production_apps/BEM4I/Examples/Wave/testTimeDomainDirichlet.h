#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../BlockMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../BESpaceTime.h"
#include "../../Vector.h"
#include "../../BEBilinearFormWave1Layer.h"

using namespace std;
using namespace bem4i;

void testTimeDomainDirichlet(
    string const &filename,
    int printPrecision
    );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  string filename = "input/cube_12.txt";
  int printPrecision = 15;
  testTimeDomainDirichlet( filename, printPrecision );

  return 0;
}

void testTimeDomainDirichlet( string const &filename, int printPrecision ) {
  typedef double SC;
  typedef int LO;
  timeval start, stop;
  std::cout.setf( std::ios::showpoint | std::ios::fixed );
  std::cout.precision( printPrecision );

  gettimeofday( &start, nullptr );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  mesh.refine( 3 );
  mesh.mapToUnitBall( );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;
  mesh.printInfo( );

  //mesh.printParaviewVtu( "output/rail_wheel_15k.vtu" );return;
  //  SC node[ 3 ];
  //  mesh.getNode( 2380, node );
  //  std::cout << "[ " << node[0] << ", " << node[1] << ", " << node[2] << " ]" << std::endl;
  //  mesh.getNode( 4299, node );
  //  std::cout << "[ " << node[0] << ", " << node[1] << ", " << node[2] << " ]" << std::endl;
  //  return;

  //quadratureType quadrature = SauterSchwab;
  //int order = 5;
  //int order2 = 5;
  //int quadOrder[] = { order, order };
  //int quadOrder2[] = { order2, order2 };
  // LO nNodes = mesh.getNNodes( );
  LO nElems = mesh.getNElements( );

  ///*

  SC endTime = 2.0;
  int nTimeSteps = 10;
  int legendreOrder = 1;
  SC dt = endTime / ( (double) ( nTimeSteps - 1.0 ) );
  int qOrder[4] = { 3, 3, 3, 3 };
  int qOrderDisj[2] = { 4, 4 };
  BESpaceTime< LO, SC > bespaceVTime( &mesh, p0, p0, legendreOrder,
      endTime, nTimeSteps );
  BEBilinearFormWave1Layer<LO, SC> formV( &bespaceVTime, qOrder, 11,
      qOrderDisj );
  //SparseMatrix<LO, SC> V;
  BlockMatrix<LO, SC> V;
  std::cout << "Assemblin' matrix V, yeah!" << std::endl;

  gettimeofday( &stop, nullptr );
  formV.assemble( V );
  //V.print();
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  BEIntegratorWave<LO, SC> integrator( &bespaceVTime, defaultQuadraturesSauterSchwab, 11 );

  Vector<LO, SC> rhs( 1 );
  integrator.getDirichletRHS( rhs );

  Vector<LO, SC> phi( bespaceVTime.getNTimeSteps( ) );
  phi.setAll( 0.0 );
  SC T, t, res = 0.0;
  for ( int j = 0; j < nTimeSteps; j++ ) {
    T = j*dt;
    for ( int i = 0; i <= std::floor( T / 2.0 ); i++ ) {
      t = T - 2.0 * i;
      res = 2.0 * ( 4.0 * t * t * t - 2.0 * t * t * t * t ) *
          std::exp( -2.0 * t );
      phi.add( j, res );
    }
  }
  phi.print( );
  //rhs.print();
  //V.print();

  Vector<LO, SC> x( ( legendreOrder + 1 ) * nTimeSteps * nElems );
  x.setAll( 0.0 );

  LeftIdentityPreconditioner<LO, SC>* M = new LeftIdentityPreconditioner<LO, SC>;

  SC prec = 1e-12;
  LO maxIt = 3000;
  V.GMRESSolve( rhs, x, prec, maxIt, 1000, M );
  x.print( );

  Vector<LO, SC> *elemDataSol[ 2 ];
  string elemNamesSol[] = { "phi", "solution" };
  std::stringstream file;
  file.fill( '0' );
  file.width( 4 );
  Vector< LO, SC > *vphi;
  Vector< LO, SC > sol( nElems );
  for ( int i = 0; i < nTimeSteps; i++ ) {
    vphi = new Vector< LO, SC >( nElems );
    vphi->setAll( 1.0 );
    vphi->scale( phi.get( i ) );
    for ( int j = 0; j < nElems; j++ ) {
      sol.set( j, x.get( i * nElems + j ) );
    }
    //sol.print();
    elemDataSol[ 0 ] = vphi;
    elemDataSol[ 1 ] = &sol;
    file.str( std::string( ) );
    file.clear( );
    file << "output/timeBEM/timeBEM_" << i << ".vtu";
    mesh.printParaviewVtu( file.str( ).c_str( ), 0, nullptr, nullptr, 2,
        elemNamesSol, elemDataSol );
    delete vphi;
  }

  Vector<LO, SC> x2( nTimeSteps );

  // sum the time basis functions
  SC *bFactors = new SC[legendreOrder + 1];
  for ( int i = 0; i < nTimeSteps; i++ ) {
    for ( int j = 0; j < legendreOrder + 1; j++ ) {
      bFactors[j] = x.get( i * nElems * ( legendreOrder + 1 ) + j * nElems );
      std::cout << bFactors[j] << " ";
    }
    std::cout << std::endl;
    for ( int j = 0; j < legendreOrder + 1; j++ ) {

      x2.add( i, bFactors[j] * integrator.evalB( i*dt, i, j ) );
    }
  }
  delete [] bFactors;
  x2.print( );



  //V.print(std::cout);

  //std::cout << V.nnz() << std::endl;

}
