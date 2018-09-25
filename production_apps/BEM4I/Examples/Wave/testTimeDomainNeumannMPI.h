#include "../../Settings.h"

#include <iostream>
#include <mpi.h>

#include "../auxiliary.h"
#include "../../BlockMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../BESpaceTime.h"
#include "../../Vector.h"
#include "../../BEBilinearFormWaveHypersingular.h"
#include "../../WavePreconditioner.h"
#include "../../PotentialsWave.h"

using namespace std;
using namespace bem4i;

void testTimeDomainNeumannMPI(
    string const &filename,
    int refine4,
    int refine9,
    bool mapToUnitBall,
    int spaceOrder,
    int timeOrder,
    int legendreOrder,
    int nTimeSteps
    );

int main(
    int argc,
    char** argv
    ) {
  if ( argc == 9 ) {
    intro( );

    string filename = string( argv[ 1 ] );
    int refine4 = atoi( argv[ 2 ] );
    int refine9 = atoi( argv[ 3 ] );
    bool mapToUnitBall = atoi( argv[ 4 ] );
    int spaceOrder = atoi( argv[ 5 ] );
    int timeOrder = atoi( argv[ 6 ] );
    int legOrder = atoi( argv[ 7 ] );
    int nTimeSteps = atoi( argv[ 8 ] );
    testTimeDomainNeumannMPI( filename, refine4, refine9, mapToUnitBall,
        spaceOrder, timeOrder, legOrder, nTimeSteps );
  }

  return 0;
}

void testTimeDomainNeumannMPI( string const &filename, int refine4, int refine9,
    bool mapToUnitBall, int spaceOrder, int timeOrder, int legendreOrder,
    int nTimeSteps ) {
  MPI_Init( nullptr, nullptr );
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  typedef double SC;
  typedef int LO;
  timeval start, stop;
  std::cout.setf( std::ios::showpoint | std::ios::fixed );
  std::cout.precision( 10 );

  gettimeofday( &start, nullptr );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  gettimeofday( &stop, nullptr );

  mesh.refine( refine4, 2 );
  mesh.refine( refine9, 3 );
  if ( mapToUnitBall ) mesh.mapToUnitBall( );
  mesh.printInfo( );

  LO nNodes = mesh.getNNodes( );

  SC endTime = 6.0;
  SC dt = endTime / ( nTimeSteps - 1 );
  int quadraturesSauterSchwab[] = { spaceOrder, spaceOrder, spaceOrder,
    spaceOrder };

  BESpaceTime< LO, SC > bespaceDTime( &mesh, p1, p1,
      legendreOrder, endTime, nTimeSteps );

  BEBilinearFormWaveHypersingular<LO, SC> formD( &bespaceDTime,
      quadraturesSauterSchwab, timeOrder );

  if ( rank == 0 ) {
    std::cout << "Using Sauter-Schwab quadrature, order " << spaceOrder
        << "." << std::endl;
    std::cout << "Using temporal quadrature of order " << timeOrder << "."
        << std::endl;
    std::cout << "Using Legendre polynomial of order " << legendreOrder << "."
        << std::endl;
  }

  MPIBlockMatrix<LO, SC> D;

  //BlockMatrix<LO, SC> D;
  if ( rank == 0 ) std::cout << "Assembling matrix D ... ";

  gettimeofday( &stop, nullptr );
  formD.assemble( D );
  MPI_Barrier( MPI_COMM_WORLD );
  gettimeofday( &stop, nullptr );

  if ( rank == 0 ) std::cout << "done in " <<
      (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  //  Vector<LO, SC> testX(D.getNRows());
  //  Vector<LO, SC> testY(D.getNRows());
  //  testX.setAll(1.0);
  //  D.apply(testX, testY);
  //  if ( rank == 0 ) testY.print();
  //  MPI_Finalize();
  //  exit(0);

  //D.print();
  // assemble RHS
  BEIntegratorWave<LO, SC> integrator( &bespaceDTime,
      quadraturesSauterSchwab, timeOrder, 0 );
  Vector<LO, SC> rhs( 1 );
  integrator.getNeumannRHS( rhs );

  // solve the system
  Vector<LO, SC> x( ( legendreOrder + 1 ) *( nTimeSteps ) * nNodes );
  x.setAll( 0.0 );
  
  LeftIdentityPreconditioner<LO, SC>* M =
      new LeftIdentityPreconditioner<LO, SC>;
  //WavePreconditioner<LO, SC>* M = new WavePreconditioner<LO, SC>( &D, 0 );
  gettimeofday( &start, nullptr );
  SC prec = 1e-12;
  LO maxIt = 3000;
  //D.FGMRESSolve( rhs, x, prec, maxIt, 400 , M );
  D.DGMRESSolve( rhs, x, prec, maxIt, 3000, 100, 4, M );


  gettimeofday( &stop, nullptr );
  if ( rank == 0 ) std::cout << "done in " <<
      (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  // evaluate the result in one point on the surface (need to sum the legendre basis)
  if ( rank == 0 ) {
    int intPoints = 4;
    SC locDt = dt / ( (double) ( intPoints + 1 ) );
    Vector<LO, SC> x2( ( nTimeSteps - 1 )*( intPoints + 1 ) + 1 );
    // sum the time basis functions
    SC evalPoint[3];
    mesh.getNode( 0, evalPoint );
    std::cout << evalPoint[2] << std::endl;


    SC *bFactors = new SC[ legendreOrder + 1 ];
    for ( int i = 0; i < nTimeSteps; i++ ) {
      for ( int j = 0; j < legendreOrder + 1; j++ ) {
        bFactors[j] = x.get( i * nNodes * ( legendreOrder + 1 ) + j * nNodes );
      }
      for ( int j = 0; j < legendreOrder + 1; j++ ) {
        //std::cout << i * dt << ": " << integrator.evalB( i*dt, i, j ) << std::endl;
        x2.add( ( intPoints + 1 ) * i, bFactors[j] * integrator.evalB( i*dt, i, j ) );

        if ( i != 0 ) {
          for ( int k = 1; k <= intPoints; k++ ) {
            x2.add( ( intPoints + 1 ) * i - k, bFactors[j] * integrator.evalB( i * dt - k*locDt, i, j ) );
          }
        }
        if ( i != nTimeSteps - 1 ) {
          for ( int k = 1; k <= intPoints; k++ ) {
            x2.add( ( intPoints + 1 ) * i + k, bFactors[j] * integrator.evalB( i * dt + k*locDt, i, j ) );
          }
        }

      }
    }
    x2.print( );
  }

  MPI_Finalize( );
  return;

  if ( rank == 0 ) {
    // compute on the grid
    std::stringstream file;
    file.fill( '0' );
    file.width( 4 );
    string file_grid = "input/scatter_grid_basic.txt";
    SurfaceMesh3D< LO, SC > mesh_scatter;
    mesh_scatter.load( file_grid.c_str( ) );
    mesh_scatter.printInfo( );
    mesh_scatter.scale( 3.0 );
    mesh_scatter.refine( 1 );
    mesh_scatter.refine( 3, 3 );
    mesh_scatter.printInfo( );
    LO nPoints = mesh_scatter.getNNodes( );
    SC *nodes = &( ( *( mesh_scatter.getNodes( ) ) )[0] );

    for ( int i = 0; i < nTimeSteps; i++ ) {
      file.str( std::string( ) );
      file.clear( );
      file << "output/timeBEM/timeBEM_Sol_" << i << ".vtu";

      Vector< LO, SC > scatter( nPoints );
      Vector< LO, SC > incident( nPoints );


      SC A = 1.0;
      SC k[3] = { M_PI, 0.0, 0.0 };
      SC omega = std::sqrt( k[0] * k[0] + k[1] * k[1] + k[2] * k[2] );
      SC mf = 2.0 * M_PI;
      SC mt = 4.0 * M_PI;
      SC phi0 = 3.5 * M_PI;
      SC t = i * dt;

      for ( LO j = 0; j < nPoints; j++ ) {
        SC *x = ( nodes + j * 3 );
        SC kx = DOT3( k, x );
        if ( omega * t - mf >= kx && kx >= omega * t - mt ) {
          incident.set( j, A * std::cos( kx + phi0 - omega * t ) );
        } else {

          incident.set( j, 0.0 );
        }
      }


      std::cout << "Evaluating double layer potential... ";
      std::cout.flush( );
      PotentialsWave<LO, SC> pw( &bespaceDTime, &x );
      pw.doubleLayerPotential( nodes, nPoints, scatter, i * dt );
      Vector< LO, SC > total( incident );
      total.add( scatter );
      string nodeNamesSol[] = { "scatter", "incident", "total" };
      Vector< LO, SC >* nodalDataSol[] = { &scatter, &incident, &total };
      mesh_scatter.printParaviewVtu( file.str( ).c_str( ), 3, nodeNamesSol, nodalDataSol, 0, nullptr, nullptr );

    }
  }

  delete M;
  MPI_Finalize( );

}