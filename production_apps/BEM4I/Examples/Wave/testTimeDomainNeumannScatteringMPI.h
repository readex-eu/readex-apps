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
#include "../../WaveScatteringProblem.h"

using namespace std;
using namespace bem4i;

void testTimeDomainNeumannScatteringMPI(
    string const &filename,
    int printPrecision
    );

double incidentWaveDuDn(
    double t,
    double *x,
    double *n
    );

double incidentWave(
    double t,
    double *x
    );

double incidentWaveDuDn2(
    double t,
    double *x,
    double *n
    );

double incidentWave2(
    double t,
    double *x
    );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  string filename = "input/cube_12.txt";
  int printPrecision = 15;
  testTimeDomainNeumannScatteringMPI( filename, printPrecision );

  return 0;
}

void testTimeDomainNeumannScatteringMPI(
    string const &filename, int printPrecision
    ) {

  MPI_Init( nullptr, nullptr );

  typedef double SC;
  typedef long LO;
  std::cout.setf( std::ios::showpoint | std::ios::fixed );
  std::cout.precision( printPrecision );
  //  int spatialQuad[4] = { 4, 4, 4, 4 };
  /*
    WaveScatteringProblem<LO, SC> wsp;
    wsp.set( INPUT_MESH_FILE, filename );
    wsp.set( OUTPUT_FILE, "output/test/test_file" );
    wsp.set( EVAL_MESH_FILE, "input/scatter_grid_x.txt" );
    wsp.set( PROBLEM_TYPE, NEUMANN );
    wsp.set( END_TIME, 12.0 );
    wsp.set( N_TIME_STEPS, 30 );
    wsp.set( LEGENDRE_ORDER, 1 );
    wsp.set( SPACE_QUAD_ORDER, &spatialQuad[0] );
    wsp.set( TEMP_QUAD_ORDER, 11 );
    wsp.set( INCIDENT_WAVE_DU_DN, &incidentWaveDuDn );
    wsp.set( INCIDENT_WAVE, &incidentWave );
    wsp.set( N_REFINE_INPUT, 1 );
    wsp.set( N_REFINE_EVAL, 6 ); //9
    wsp.set( SCALE_INPUT, 2.0 );
    wsp.set( SCALE_EVAL, 7.5 );
    wsp.set( SOLVER_MAX_IT, 3000 );
    wsp.set( SOLVER_RESTARTS, 1000 );
    wsp.set( ITERATIVE_SOLVER, DGMRES );
    wsp.set( DGMRES_MAX_DIM, 600 );
    wsp.set( DGMRES_MAX_DIM_IT, 8 );
    wsp.set( SOLVER_PRECISION, 1e-4 );
    wsp.solve( );
    MPI_Finalize( );

    return;
   */
  timeval start, stop;
  gettimeofday( &start, NULL );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  gettimeofday( &stop, NULL );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;
  mesh.printInfo( );

  mesh.refine( 2 );
  //mesh.mapToUnitBall( );
  mesh.printInfo( );


  //LO nElems = mesh.getNElements( );
  LO nNodes = mesh.getNNodes( );

  SC endTime = 10.0;
  int nTimeSteps = 5;
  int legendreOrder = 1;
  int nPos = 0;
  SC dt = endTime / ( nTimeSteps - 1 );
  //int nPre = 0;
  int quadraturesSauterSchwab[] = { 3, 3, 3, 3 };
  int quadratureDisjoint[] = { 4, 4 };


  BESpaceTime< LO, SC > bespaceDTime( &mesh, p1, p1, legendreOrder, endTime, nTimeSteps );

  BEBilinearFormWaveHypersingular<LO, SC> formD( &bespaceDTime, quadraturesSauterSchwab, 9, quadratureDisjoint );

  //MPIBlockMatrix<LO, SC> D;
  //SparseMatrix<LO, SC> D;
  BlockMatrix<LO, SC> D;
  std::cout << "Assembling matrix D ... ";

  gettimeofday( &stop, NULL );
  formD.assemble( D );

  //D.print();
  //D.print( );


  gettimeofday( &stop, NULL );

  //return;
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  BEIntegratorWave<LO, SC> integrator( &bespaceDTime, defaultQuadraturesSauterSchwab, defaultTimeQuadrature, 0, nPos );
  Vector<LO, SC> rhs( 1 );
  integrator.getNeumannRHS( rhs );

  Vector<LO, SC> x( ( legendreOrder + 1 ) *( nTimeSteps ) * nNodes );
  x.setAll( 0.0 );


  //LeftIdentityPreconditioner<LO, SC>* M = new LeftIdentityPreconditioner<LO, SC>;
  WavePreconditioner<LO, SC>* M = new WavePreconditioner<LO, SC>( &D, 0 );
  gettimeofday( &start, NULL );
  //LeftIdentityPreconditioner<LO, SC>* M = new LeftIdentityPreconditioner<LO, SC>;
  D.GMRESSolve( rhs, x, 1e-5, 100000, 100, M );

  gettimeofday( &stop, NULL );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;
  //return;
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if ( rank == 0 ) {
    // compute on the grid
    std::stringstream file;
    file.fill( '0' );
    file.width( 4 );
    string file_grid = "input/scatter_grid_x.txt";
    SurfaceMesh3D< LO, SC > mesh_scatter;
    mesh_scatter.load( file_grid.c_str( ) );
    mesh_scatter.printInfo( );
    mesh_scatter.scale( 5.0 );
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


      /*SC A = 1.0;
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
      }*/
      SC s[3] = { 0.0, 3.0, 3.0 };
      SC R = 0.05;
      SC t = i * dt;
      for ( LO j = 0; j < nPoints; j++ ) {
        SC *x = ( nodes + j * 3 );
        SC r = DIST3( s, x );

        if ( ( r - t + 12 * R >= 0.0 && r - t - R <= 0.0 ) && r >= 0.9 * R ) {
          incident.set( j, 1 / r * ( 0.75 - std::cos( M_PI * ( r - t + 3 * R ) / ( 2 * R ) ) + 0.25 * std::cos( M_PI * ( r - t + 3 * R ) / ( R ) ) ) );
        } else {
          incident.set( j, 0.0 );
        }
      }


      std::cout << "Evaluating double layer potential... ";
      std::cout.flush( );
      //x.print();
      PotentialsWave<LO, SC> pw( &bespaceDTime, &x );
      pw.doubleLayerPotential( nodes, nPoints, scatter, i * dt );
      //scatter.print();
      Vector< LO, SC > total( incident );
      total.add( scatter );
      string nodeNamesSol[] = { "scatter", "incident", "total" };
      Vector< LO, SC >* nodalDataSol[] = { &scatter, &incident, &total };
      mesh_scatter.printParaviewVtu( file.str( ).c_str( ), 3, nodeNamesSol, nodalDataSol, 0, NULL, NULL );

    }
  }
  if ( rank == 0 ) {
    mesh.printParaviewVtu( "output/timeBEM/scatterer.vtu" );
  }


  delete M;
  MPI_Finalize( );
}

double incidentWaveDuDn( double t, double *x, double *n ) {
  double A = 1.0;
  double k[3] = { 0.0, -M_PI / std::sqrt( 2.0 ), -M_PI / std::sqrt( 2.0 ) };
  double omega = std::sqrt( k[0] * k[0] + k[1] * k[1] + k[2] * k[2] );
  double mf = 6.0 * M_PI;
  double mt = 8.0 * M_PI;
  double phi0 = 3.5 * M_PI;
  double kx = DOT3( k, x );
  double kn = DOT3( k, n );
  if ( ( omega * t - mf >= kx ) && ( kx >= omega * t - mt ) ) {
    return A * std::sin( kx + phi0 - omega * t ) * kn;
  } else {

    return 0.0;
  }
}

double incidentWave( double t, double *x ) {
  double A = 1.0;
  double k[3] = { 0.0, -M_PI / std::sqrt( 2.0 ), -M_PI / std::sqrt( 2.0 ) };
  double omega = std::sqrt( k[0] * k[0] + k[1] * k[1] + k[2] * k[2] );
  double mf = 6.0 * M_PI;
  double mt = 8.0 * M_PI;
  double phi0 = 3.5 * M_PI;
  double kx = DOT3( k, x );
  if ( ( omega * t - mf >= kx ) && ( kx >= omega * t - mt ) ) {
    return A * std::cos( kx + phi0 - omega * t );
  } else {

    return 0.0;
  }
}

double incidentWave2( double t, double *x ) {

  double s1[3] = { 0.0, 3.0, 3.0 };
  double s2[3] = { 0.0, -3.0, 3.0 };
  double R1 = 0.5;
  double R2 = 0.5;

  double r1 = DIST3( s1, x );
  double r2 = DIST3( s2, x );

  double ret1 = 0.0, ret2 = 0.0;
  if ( ( r1 - t + 12 * R1 >= 0.0 && r1 - t - R1 <= 0.0 ) && r1 >= R1 ) {
    ret1 = 1 / r1 * ( 0.75 - std::cos( M_PI * ( r1 - t + 3 * R1 ) / ( 2 * R1 ) ) + 0.25 * std::cos( M_PI * ( r1 - t + 3 * R1 ) / ( R1 ) ) );
  } else {
    ret1 = 0.0;
  }
  if ( ( r2 - t + 12 * R2 >= 0.0 && r2 - t - R2 <= 0.0 ) && r2 >= R2 ) {
    ret2 = 1 / r2 * ( 0.75 - std::cos( M_PI * ( r2 - t + 3 * R2 ) / ( 2 * R2 ) ) + 0.25 * std::cos( M_PI * ( r2 - t + 3 * R2 ) / ( R2 ) ) );
  } else {

    ret2 = 0.0;
  }
  return ret1 + ret2;
}

double incidentWaveDuDn2( double t, double *x, double *n ) {

  double s1[3] = { 0.0, 3.0, 3.0 };
  double s2[3] = { 0.0, -3.0, 3.0 };
  double R1 = 0.5;
  double R2 = 0.5;

  double r1 = DIST3( s1, x );
  double r2 = DIST3( s2, x );

  double ret1, ret2;
  if ( r1 - t + 12 * R1 >= 0.0 && r1 - t - R1 <= 0.0 ) {
    double dfdx1 = ( ( M_PI * std::sin( ( M_PI * ( 3 * R1 - t + r1 ) ) / R1 )*( 2 * s1[0] - 2 * x[0] ) ) / ( 8 * R1 * r1 ) - ( M_PI * std::sin( ( M_PI * ( 3 * R1 - t + r1 ) ) / ( 2 * R1 ) )*( 2 * s1[0] - 2 * x[0] ) ) / ( 4 * R1 * r1 ) ) / r1 + ( ( 2 * s1[0] - 2 * x[0] )*( std::cos( ( M_PI * ( 3 * R1 - t + r1 ) ) / R1 ) / 4 - std::cos( ( M_PI * ( 3 * R1 - t + r1 ) ) / ( 2 * R1 ) ) + 3 / 4 ) ) / ( 2 * r1 * std::sqrt( r1 ) );
    double dfdx2 = ( ( M_PI * std::sin( ( M_PI * ( 3 * R1 - t + r1 ) ) / R1 )*( 2 * s1[1] - 2 * x[1] ) ) / ( 8 * R1 * r1 ) - ( M_PI * std::sin( ( M_PI * ( 3 * R1 - t + r1 ) ) / ( 2 * R1 ) )*( 2 * s1[1] - 2 * x[1] ) ) / ( 4 * R1 * r1 ) ) / r1 + ( ( 2 * s1[1] - 2 * x[1] )*( std::cos( ( M_PI * ( 3 * R1 - t + r1 ) ) / R1 ) / 4 - std::cos( ( M_PI * ( 3 * R1 - t + r1 ) ) / ( 2 * R1 ) ) + 3 / 4 ) ) / ( 2 * r1 * std::sqrt( r1 ) );
    double dfdx3 = ( ( M_PI * std::sin( ( M_PI * ( 3 * R1 - t + r1 ) ) / R1 )*( 2 * s1[2] - 2 * x[2] ) ) / ( 8 * R1 * r1 ) - ( M_PI * std::sin( ( M_PI * ( 3 * R1 - t + r1 ) ) / ( 2 * R1 ) )*( 2 * s1[2] - 2 * x[2] ) ) / ( 4 * R1 * r1 ) ) / r1 + ( ( 2 * s1[2] - 2 * x[2] )*( std::cos( ( M_PI * ( 3 * R1 - t + r1 ) ) / R1 ) / 4 - std::cos( ( M_PI * ( 3 * R1 - t + r1 ) ) / ( 2 * R1 ) ) + 3 / 4 ) ) / ( 2 * r1 * std::sqrt( r1 ) );

    double dfdn = dfdx1 * n[0] + dfdx2 * n[1] + dfdx3 * n[2];
    ret1 = dfdn;
  } else {
    ret1 = 0.0;
  }
  if ( r2 - t + 12 * R2 >= 0.0 && r2 - t - R2 <= 0.0 ) {
    double dfdx1 = ( ( M_PI * std::sin( ( M_PI * ( 3 * R2 - t + r2 ) ) / R2 )*( 2 * s2[0] - 2 * x[0] ) ) / ( 8 * R2 * r2 ) - ( M_PI * std::sin( ( M_PI * ( 3 * R2 - t + r2 ) ) / ( 2 * R2 ) )*( 2 * s2[0] - 2 * x[0] ) ) / ( 4 * R2 * r2 ) ) / r2 + ( ( 2 * s2[0] - 2 * x[0] )*( std::cos( ( M_PI * ( 3 * R2 - t + r2 ) ) / R2 ) / 4 - std::cos( ( M_PI * ( 3 * R2 - t + r2 ) ) / ( 2 * R2 ) ) + 3 / 4 ) ) / ( 2 * r2 * std::sqrt( r2 ) );
    double dfdx2 = ( ( M_PI * std::sin( ( M_PI * ( 3 * R2 - t + r2 ) ) / R2 )*( 2 * s2[1] - 2 * x[1] ) ) / ( 8 * R2 * r2 ) - ( M_PI * std::sin( ( M_PI * ( 3 * R2 - t + r2 ) ) / ( 2 * R2 ) )*( 2 * s2[1] - 2 * x[1] ) ) / ( 4 * R2 * r2 ) ) / r2 + ( ( 2 * s2[1] - 2 * x[1] )*( std::cos( ( M_PI * ( 3 * R2 - t + r2 ) ) / R2 ) / 4 - std::cos( ( M_PI * ( 3 * R2 - t + r2 ) ) / ( 2 * R2 ) ) + 3 / 4 ) ) / ( 2 * r2 * std::sqrt( r2 ) );
    double dfdx3 = ( ( M_PI * std::sin( ( M_PI * ( 3 * R2 - t + r2 ) ) / R2 )*( 2 * s2[2] - 2 * x[2] ) ) / ( 8 * R2 * r2 ) - ( M_PI * std::sin( ( M_PI * ( 3 * R2 - t + r2 ) ) / ( 2 * R2 ) )*( 2 * s2[2] - 2 * x[2] ) ) / ( 4 * R2 * r2 ) ) / r2 + ( ( 2 * s2[2] - 2 * x[2] )*( std::cos( ( M_PI * ( 3 * R2 - t + r2 ) ) / R2 ) / 4 - std::cos( ( M_PI * ( 3 * R2 - t + r2 ) ) / ( 2 * R2 ) ) + 3 / 4 ) ) / ( 2 * r2 * std::sqrt( r2 ) );

    double dfdn = dfdx1 * n[0] + dfdx2 * n[1] + dfdx3 * n[2];
    ret2 = dfdn;
  } else {

    ret2 = 0.0;
  }
  return ret1 + ret2;
}