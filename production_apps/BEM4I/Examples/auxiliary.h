#include <sys/time.h>
#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <string>
#ifdef FENV
#include <fenv.h>
#endif

using namespace std;

void intro( );

void printLogo( );

// Deprecated! Use progress monitor instead! 
double timeDiff(
    timeval start,
    timeval end
    );

double timeDiff( timeval start, timeval end ) {
  return (( end.tv_sec - start.tv_sec ) * 1000000u +
      end.tv_usec - start.tv_usec ) / 1.e6;
}

void intro( ) {
  // enables throwing exception when operation results in NaN or Inf
#ifdef FENV
  feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
#endif
  printLogo( );

  std::cout << "Using " << omp_get_max_threads( ) << " OMP threads."
      << std::endl;

#if N_MIC > 0
  for ( int device = 0; device < N_MIC; ++device ) {
#pragma offload_transfer target(mic:device)
  }
#endif
}

void printLogo( ) {
  string f = "               ";
  std::cout << std::endl;
  std::cout << f << "    _/_/_/    _/_/_/_/  _/      _/  _/  _/  _/_/_/" << f <<
      std::endl;
  std::cout << f << "   _/    _/  _/        _/_/  _/_/  _/  _/    _/   " << f <<
      std::endl;
  std::cout << f << "  _/_/_/    _/_/_/    _/  _/  _/  _/_/_/_/  _/    " << f <<
      std::endl;
  std::cout << f << " _/    _/  _/        _/      _/      _/    _/     " << f <<
      std::endl;
  std::cout << f << "_/_/_/    _/_/_/_/  _/      _/      _/  _/_/_/    " << f <<
      std::endl;
  std::cout << std::endl;
}
