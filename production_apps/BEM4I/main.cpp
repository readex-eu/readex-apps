/*!
 * @file    main.cpp
 * @author  Michal Merta
 * @author  Jan Zapletal
 * @date    July 3, 2013
 * @brief   Main class
 *
 */

/*
 * Set the EXAMPLE variable to one of the following values
 * 
 * LAPLACE SOLVER
 * 100: testLaplaceDirichlet
 * 101: testLaplaceDirichletACA
 * 102: testLaplaceDirichletP0P0
 * 103: testLaplaceDirichletP1P1
 * 104: testLaplaceNeumann
 * 105: testLaplaceSteklovPoincare
 * 106: testFastBEM
 * 
 * HELMHOLTZ SOLVER
 * 200: testHelmholtzDirichlet
 * 201: testHelmholtzDirichletACA
 * 202: testHelmholtzDirichletP1P1
 * 203: testHelmholtzDirichletReg
 * 204: testHelmholtzNeumann
 * 205: testHelmholtzNeumannACA
 * 206: testHelmholtzNeumannP0P0
 * 207: testHelmholtzNeumannReg
 * 208: testHelmholtzSoftScatter
 * 209: testHelmholtzSoftScatterACA
 * 210: testHelmholtzHardScatter
 * 211: testHelmholtzHardScatterReg
 * 212: testSTF
 * 213: testSTFP0P0
 * 
 * LAME SOLVER
 * 300: testLameMixed
 * 301: testLameSteklovPoincare
 * 
 * TIME DOMAIN WAVE SOLVER
 * 400: testTimeDomainDirichlet
 * 401: testTimeDomainNeumann
 * 402: testTimeDomainNeumannMPI
 * 403: testTimeDomainNeumannScatteringMPI
 * 
 * OPTIMIZATION SOLVER
 * 500: testBernoulli
 * 501: testHeatSource
 * 
 * MISCELLANEOUS
 * 900: testHomogenizationProblem
 * 901: testRohan
 * 902: testMatrixAssemblyTime
 * 903: testSparseMatrix
 * 904: testMPIMatrix
 * 905: aux
 */

#define EXAMPLE 200

#if EXAMPLE == 100
#include "Examples/Laplace/testLaplaceDirichlet.h"
#elif EXAMPLE == 101
#include "Examples/Laplace/testLaplaceDirichletACA.h"
#elif EXAMPLE == 102
#include "Examples/Laplace/testLaplaceDirichletP0P0.h"
#elif EXAMPLE == 103
#include "Examples/Laplace/testLaplaceDirichletP1P1.h"
#elif EXAMPLE == 104
#include "Examples/Laplace/testLaplaceNeumann.h"
#elif EXAMPLE == 105
#include "Examples/Laplace/testLaplaceSteklovPoincare.h"
#elif EXAMPLE == 106
#include "Examples/Laplace/testFastBEM.h"
#elif EXAMPLE == 200
#include "Examples/Helmholtz/testHelmholtzDirichlet.h"
#elif EXAMPLE == 201
#include "Examples/Helmholtz/testHelmholtzDirichletACA.h"
#elif EXAMPLE == 202
#include "Examples/Helmholtz/testHelmholtzDirichletP1P1.h"
#elif EXAMPLE == 203
#include "Examples/Helmholtz/testHelmholtzDirichletReg.h"
#elif EXAMPLE == 204
#include "Examples/Helmholtz/testHelmholtzNeumann.h"
#elif EXAMPLE == 205
#include "Examples/Helmholtz/testHelmholtzNeumannACA.h"
#elif EXAMPLE == 206
#include "Examples/Helmholtz/testHelmholtzNeumannP0P0.h"
#elif EXAMPLE == 207
#include "Examples/Helmholtz/testHelmholtzNeumannReg.h"
#elif EXAMPLE == 208
#include "Examples/Helmholtz/testHelmholtzSoftScatter.h"
#elif EXAMPLE == 209
#include "Examples/Helmholtz/testHelmholtzSoftScatterACA.h"
#elif EXAMPLE == 210
#include "Examples/Helmholtz/testHelmholtzHardScatter.h"
#elif EXAMPLE == 211
#include "Examples/Helmholtz/testHelmholtzHardScatterReg.h"
#elif EXAMPLE == 212
#include "Examples/Helmholtz/testSTF.h"
#elif EXAMPLE == 213
#include "Examples/Helmholtz/testSTFP0P0.h"
#elif EXAMPLE == 300
#include "Examples/Lame/testLameMixed.h"
#elif EXAMPLE == 301
#include "Examples/Lame/testLameSteklovPoincare.h"
#elif EXAMPLE == 400
#include "Examples/Wave/testTimeDomainDirichlet.h"
#elif EXAMPLE == 401
#include "Examples/Wave/testTimeDomainNeumann.h"
#elif EXAMPLE == 402
#include "Examples/Wave/testTimeDomainNeumannMPI.h"
#elif EXAMPLE == 403
#include "Examples/Wave/testTimeDomainNeumannScatteringMPI.h"
#elif EXAMPLE == 500
#include "Examples/Optimization/testBernoulli.h"
#elif EXAMPLE == 501
#include "Examples/Optimization/testHeatSource.h"
#elif EXAMPLE == 900
#include "Examples/Miscellaneous/testHomogenizationProblem.h"
#elif EXAMPLE == 901
#include "Examples/Miscellaneous/testRohan.h"
#elif EXAMPLE == 902
#include "Examples/Miscellaneous/testMatrixAssemblyTime.h"
#elif EXAMPLE == 903
#include "Examples/Miscellaneous/testSparseMatrix.h"
#elif EXAMPLE == 904
#include "Examples/Miscellaneous/testMPIMatrix.h"
#elif EXAMPLE == 905
#include "Examples/Miscellaneous/aux.h"
#else 
#include <iostream>

int main(
    int argc,
    char ** argv
    ) {
  std::cout << "The BEM4I Library.\n\n" <<
      "Incorrect problem. See the main.cpp file" << std::endl;
  return 0;
}
#endif
