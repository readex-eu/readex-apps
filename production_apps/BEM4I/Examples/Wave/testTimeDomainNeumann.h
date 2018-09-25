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

void testTimeDomainNeumann(
    string const &filename,
    int printPrecision
    );

double neumannWavePhi(
    double* x,
    double t,
    int n
    );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  string filename = "input/cube_12.txt";
  int printPrecision = 15;
  testTimeDomainNeumann( filename, printPrecision );

  return 0;
}

void testTimeDomainNeumann( string const &filename, int printPrecision ) {
  MPI_Init( nullptr, nullptr );
  typedef double SC;
  typedef int LO;
  timeval start, stop;
  std::cout.setf( std::ios::showpoint | std::ios::fixed );
  std::cout.precision( printPrecision );

  //      for ( int i = 0; i < 201; i++ ) {
  //        double t = i * 0.02;
  //        double x[3] = { 0.0, 0.0, 1.0 };
  //        std::cout << neumannWavePhi( x, t, 0 ) << " ";
  //      }
  //      std::cout << std::endl;
  //      return;


  gettimeofday( &start, nullptr );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;
  //mesh.printInfo( );


  mesh.refine( 3 );
  mesh.mapToUnitBall( );
  mesh.printInfo( );

  //LO nElems = mesh.getNElements( );
  LO nNodes = mesh.getNNodes( );

  SC endTime = 4.0;
  int nTimeSteps = 10;
  int legendreOrder = 1;
  int nPos = 0;
  int nPre = 0;
  int quadraturesSauterSchwab[] = { 3, 3, 3, 3 };
  int disQOrder[2] = { 4, 4 };

  SC dt = endTime / ( nTimeSteps - 1 );

  BESpaceTime< LO, SC > bespaceDTime( &mesh, p1, p1, legendreOrder, endTime, nTimeSteps );
#ifndef EXPERIMENTAL_WAVE
  BEBilinearFormWaveHypersingular<LO, SC> formD( &bespaceDTime,
      quadraturesSauterSchwab, 11, disQOrder );
#else
  nPos = 10;
  nPre = 10;
  BEBilinearFormWaveHypersingular<LO, SC> formD( &bespaceDTime, nPre, nPos, quadraturesSauterSchwab );

#endif
  //SparseMatrix<LO, SC> D;
  MPIBlockMatrix<LO, SC> D;
  std::cout << "Assemblin' matrix D, yeah!" << std::endl;

  gettimeofday( &stop, nullptr );
  formD.assemble( D );

  //D.saveTripletsBin("D.mat");
  //D.loadTripletsBin("D.mat");

  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  BEIntegratorWave<LO, SC> integrator( &bespaceDTime, quadraturesSauterSchwab,
      11, disQOrder, nPre, nPos );
  Vector<LO, SC> rhs( 1 );
  integrator.getNeumannRHS( rhs );
  //rhs.print( );


#ifndef EXPERIMENTAL_WAVE
  Vector<LO, SC> x( ( legendreOrder + 1 ) *( nTimeSteps ) * nNodes );
#else
  Vector<LO, SC> x( ( nTimeSteps + nPos + nPre ) * nNodes );
#endif
  x.setAll( 0.0 );
  std::cout << x.getLength( ) << " " << rhs.getLength( ) << " " << D.getNRows( ) << " " << D.getNCols( ) << std::endl;


  //LeftIdentityPreconditioner<LO, SC>* M = new LeftIdentityPreconditioner<LO, SC>;
  int iters[2] = { 10, 2 };
  WavePreconditioner<LO, SC>* M = new WavePreconditioner<LO, SC>( &D, 2, &iters[0] );
  gettimeofday( &start, nullptr );
  //LeftIdentityPreconditioner<LO, SC>* M = new LeftIdentityPreconditioner<LO, SC>;
  SC prec = 1e-5;
  LO maxIt = 1000;
  //D.DGMRESSolve( rhs, x, prec, maxIt, 1000, 200, 10, M );
  D.FGMRESSolve( rhs, x, prec, maxIt, 200, M );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;
  //D.LUsolve( rhs, x );
  //D.BiCGStabSolve( rhs, x );
  //D.PARDISOSolve(rhs, x);
  //  x.print( );
  delete M;

  // evaluate the result in one point on the surface (need to sum the legendre basis)
  int intPoints = 8;
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

  // L2 space time error
  int order = 5;
  int timeOrder = 6;
  int qSize = lineQuadSizes[timeOrder];
  double *qPoints = lineQuadPoints[timeOrder];
  double *qWeights = lineQuadWeights[timeOrder];
  double t;
  double startT;
  double aux1 = 0.0;
  double analNorm = 0.0;

  SC *bFactors11 = new SC[ legendreOrder + 1 ];
  SC *bFactors12 = new SC[ legendreOrder + 1 ];
  SC *bFactors13 = new SC[ legendreOrder + 1 ];
  SC *bFactors21 = new SC[ legendreOrder + 1 ];
  SC *bFactors22 = new SC[ legendreOrder + 1 ];
  SC *bFactors23 = new SC[ legendreOrder + 1 ];

  double spaceTimeError = 0.0;
  for ( int i = 0; i < nTimeSteps - 1; i++ ) {
    startT = i * dt;


    for ( int k = 0; k < qSize; k++ ) {
      t = startT + dt * qPoints[k];
      double spaceError = 0.0;
      double spaceAnalNorm = 0.0;
      LO ind[ 3 ];
      double x1[3], x2[3], x3[3], xp[3];
      int size = quadSizes[ order ];
      double* quad = new double[ 3 * size ];

      for ( LO e = 0; e < mesh.getNElements( ); e++ ) {
        mesh.getElement( e, ind );
        mesh.getNodes( e, x1, x2, x3 );
        mesh.getQuadratureNodes( x1, x2, x3, quad, order );

        for ( int j = 0; j < legendreOrder + 1; j++ ) {
          bFactors11[j] = x.get( i * nNodes * ( legendreOrder + 1 ) + j * nNodes + ind[0] );
          bFactors12[j] = x.get( i * nNodes * ( legendreOrder + 1 ) + j * nNodes + ind[1] );
          bFactors13[j] = x.get( i * nNodes * ( legendreOrder + 1 ) + j * nNodes + ind[2] );
          bFactors21[j] = x.get( ( i + 1 ) * nNodes * ( legendreOrder + 1 ) + j * nNodes + ind[0] );
          bFactors22[j] = x.get( ( i + 1 ) * nNodes * ( legendreOrder + 1 ) + j * nNodes + ind[1] );
          bFactors23[j] = x.get( ( i + 1 ) * nNodes * ( legendreOrder + 1 ) + j * nNodes + ind[2] );
        }
        aux1 = 0.0;
        for ( int j = 0; j < size; j++ ) {
          xp[0] = quad[3 * j];
          xp[1] = quad[3 * j + 1];
          xp[2] = quad[3 * j + 2];
          double anal = neumannWavePhi( xp, t, 0 ); // Helmholtz

          double val = 0.0;
          for ( int l = 0; l < legendreOrder + 1; l++ ) {
            val += ( bFactors11[l] * integrator.evalB( t, i, l )*( 1.0 - quadPoints[ order ][ 2 * j ] - quadPoints[ order ][ 2 * j + 1 ] )
                + bFactors12[l] * integrator.evalB( t, i, l ) * quadPoints[ order ][ 2 * j ]
                + bFactors13[l] * integrator.evalB( t, i, l ) * quadPoints[ order ][ 2 * j + 1 ]
                + bFactors21[l] * integrator.evalB( t, i + 1, l )*( 1.0 - quadPoints[ order ][ 2 * j ] - quadPoints[ order ][ 2 * j + 1 ] )
                + bFactors22[l] * integrator.evalB( t, i + 1, l ) * quadPoints[ order ][ 2 * j ]
                + bFactors23[l] * integrator.evalB( t, i + 1, l ) * quadPoints[ order ][ 2 * j + 1 ] );
          }
          //std::cout << anal << " " << val << std::endl;
          aux1 += quadWeights[ order ][ j ] * pow( std::abs( anal - val ), 2 );
          spaceAnalNorm += anal * anal * quadWeights[ order ][ j ] * mesh.getElemArea( e );
        }
        spaceError += aux1 * mesh.getElemArea( e );

      }
      std::cout << spaceError << std::endl;
      delete [] quad;
      spaceTimeError += spaceError * qWeights[k] * dt;
      analNorm += spaceAnalNorm * qWeights[k] * dt;
    }

  }
  std::cout << "aN " << analNorm << std::endl;
  std::cout << std::sqrt( spaceTimeError ) / std::sqrt( analNorm ) << std::endl;
  //  // try evaluating double layer potential
  //  SC dt = endTime / ( nTimeSteps - 1 );
  //  SC testPoint[3] = { 1.2, 1.2, 1.2 };
  //  Vector<LO, SC> values( 1 );
  PotentialsWave<LO, SC> pw( &bespaceDTime, &x );
  //  for ( int i = 0; i < nTimeSteps; i++ ) {
  //    pw.doubleLayerPotential( testPoint, 1, values, i * dt );
  //    std::cout << values.get( 0 ) << " ";
  //  }
  //  std::cout << std::endl;

  // compute on the grid
  std::stringstream file;
  file.fill( '0' );
  file.width( 4 );
  string file_grid = "input/scatter_grid_x.txt";
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
    pw.doubleLayerPotential( nodes, nPoints, scatter, i * dt );
    Vector< LO, SC > total( incident );
    total.add( scatter );
    string nodeNamesSol[] = { "scatter", "incident", "total" };
    Vector< LO, SC >* nodalDataSol[] = { &scatter, &incident, &total };
    mesh_scatter.printParaviewVtu( file.str( ).c_str( ), 3, nodeNamesSol, nodalDataSol, 0, nullptr, nullptr );

  }
  delete [] bFactors;
  delete [] bFactors11;
  delete [] bFactors12;
  delete [] bFactors13;
  delete [] bFactors21;
  delete [] bFactors22;
  delete [] bFactors23;


  //  Vector<LO, SC>** results = new Vector<LO, SC>*[2 * nTimeSteps - 1];
  //  for ( int i = 0; i < 2 * nTimeSteps - 1; i++ ) {
  //    results[i] = new Vector<LO, SC>( nNodes );
  //  }



  // sum the time basis functions

  //  string nodeNames[] = { "neumann" };
  //  Vector< LO, SC >* nodalData[] = { &x };
  //  mesh.printParaviewVtu( "output/sol.vtu", 1, nodeNames, nodalData, 0, nullptr, nullptr );
  //
  //  for ( int i = 0; i < nTimeSteps; i++ ) {
  //    for ( LO k = 0; k < nNodes; k++ ) {
  //
  //      for ( int j = 0; j < legendreOrder + 1; j++ ) {
  //        bFactors[j] = x.get( i * nNodes * ( legendreOrder + 1 ) + j * nNodes + k );
  //      }
  //
  //      for ( int j = 0; j < legendreOrder + 1; j++ ) {
  //        //std::cout << i * dt << ": " << integrator.evalB( i*dt, i, j ) << std::endl;
  //        ( results[2 * i] )->add( k, bFactors[j] * integrator.evalB( i*dt, i, j ) );
  //        //x2.add( 2 * i, bFactors[j] * integrator.evalB( i*dt, i, j ) );
  //
  //        if ( i != 0 ) {
  //          results[2 * i - 1]->add( k, bFactors[j] * integrator.evalB( i * dt - dt / 2.0, i, j ) );
  //          //x2.add( 2 * i - 1, bFactors[j] * integrator.evalB( i * dt - dt / 2.0, i, j ) );
  //        }
  //        if ( i != nTimeSteps - 1 ) {
  //          results[2 * i + 1]->add( k, bFactors[j] * integrator.evalB( i * dt + dt / 2.0, i, j ) );
  //          //x2.add( 2 * i + 1, bFactors[j] * integrator.evalB( i * dt + dt / 2.0, i, j ) );
  //        }
  //      }
  //    }
  //  }
  //
  //
  //
  //  for ( int i = 0; i < 2 * nTimeSteps - 1; i++ ) {
  //
  //    file.str( std::string( ) );
  //    file.clear( );
  //    file << "output/timeBEM/timeBEM_Neu_" << i << ".vtu";
  //    Vector< LO, SC >* nodalData[] = { results[i] };
  //    mesh.printParaviewVtu( file.str( ).c_str( ), 1, nodeNames, nodalData, 0, nullptr, nullptr );
  //  }
  //  for ( int i = 0; i < 2 * nTimeSteps - 1; i++ ) {
  //    delete results[i];
  //  }
}

double neumannWavePhi(
    double* x,
    double T,
    int n
    ) {

  if ( T <= 0.0 ) return 0.0;
  // number of integration intervals per second
  int Ns = 1000;
  // number of integration intervals
  int N = Ns * T;
  // quadrature order
  int order = 5;

  int qSize = lineQuadSizes[order];
  double *qPoints = lineQuadPoints[order];
  double *qWeights = lineQuadWeights[order];

  double t = 0.0;
  double val1 = 0.0, val2 = 0.0;
  double start;
  double length = T / ( (double) N );
  double length2 = 0.0;

  if ( n == 0 ) {
    for ( int i = 0; i < N; i++ ) {
      start = i * length;

      for ( int k = 0; k < qSize; k++ ) {
        t = start + length * qPoints[k];
        val1 += ( -2 * std::sin( 3 * ( T - t ) )*( T - t )*( T - t ) *
            std::exp( -( T - t ) ) * std::cosh( t ) / std::sqrt( M_PI ) / 2.0 ) * qWeights[k];
      }
    }
    int limK = (int) floor( T / 2.0 );

    if ( T > 2.0 ) {
      for ( int k = 1; k <= limK; ++k ) {
        double start2 = 2 * k;
        N = Ns * ( T - start2 );
        if ( N > 0 ) {
          length2 = ( T - start2 ) / ( (double) N );
        } else {
          length2 = 0.0;
        }
        for ( int l = 1; l <= k; ++l ) {

          double locVal = 0.0;

          // precompute factorials
          int km1 = 1, lm1 = 1, kmlp1 = 1, kml = 1;
          for ( int f = 1; f <= k - 1; f++ ) km1 *= f;
          for ( int f = 1; f <= l - 1; f++ ) lm1 *= f;
          for ( int f = 1; f <= k - l + 1; f++ ) kmlp1 *= f;
          for ( int f = 1; f <= k - l; f++ ) kml *= f;
          double c = ( (double) ( km1 * pow( 2.0, k - l ) ) ) /
              ( (double) ( lm1 * kml * kmlp1 ) );

          for ( int i = 0; i < N; i++ ) {
            start = start2 + i * length2;
            for ( int q = 0; q < qSize; q++ ) {

              t = start + length2 * qPoints[q];
              locVal += ( ( pow( t - 2 * k, k - l + 1 )
                  * std::exp( t - 2 * k ) * std::exp( -T + t ) * ( T - t ) * ( 3 * ( T - t ) * std::cos( 3 * ( T - t ) )-
                  ( ( T - t ) - 2 ) * std::sin( 3 * ( T - t ) ) ) ) / std::sqrt( M_PI ) / 2.0 )
                  * qWeights[q];
            }
          }

          locVal *= c * ( pow( -1, k + 1 ) );
          locVal *= length2;
          val2 += locVal;
        }
      }
    }
  } else if ( n == 1 ) {
    for ( int i = 0; i < N; i++ ) {
      start = i * length;

      for ( int k = 0; k < qSize; k++ ) {
        t = start + length * qPoints[k];
        val1 += ( -2 * std::sin( 2 * M_PI * ( T - t ) )*( T - t )*( T - t )*( T - t ) *
            std::exp( -2 * ( T - t ) ) *
            std::cosh( t ) * std::cos( t ) / std::sqrt( 3.0 / M_PI ) / 2.0 * x[2] )
            * qWeights[k];
      }
    }
  }
  val1 *= length;

  return -( val1 + 2.0 * val2 );
}