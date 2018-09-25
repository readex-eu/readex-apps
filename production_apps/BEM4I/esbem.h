#ifndef BEM_H_
#define BEM_H_

#include "Settings.h"
#include <mpi.h>

#ifdef FENV
#include <fenv.h>
#endif

#include "omp.h"

#include <cstdlib>
#include <iostream>
#include "ProgressMonitor.h"
#include "FullMatrix.h"
#include "SurfaceMesh3D.h"
#include "BESpace.h"
#include "Vector.h"
#include "IdentityOperator.h"
#include "SparseMatrix.h"
#include "BEBilinearFormLaplace2Layer.h"
#include "BEBilinearFormLame1Layer.h"
#include "BEBilinearFormLame2Layer.h"
#include "BEBilinearFormLameHypersingular.h"

// intel bug ??? workaround
#include "BEIntegratorWave.h"
#include "Eigen/Eigenvalues"

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

template void bem4i::getLameSteklovPoincare< eslocal, double >(
    double * Sarray,
    eslocal nNodes,
    const double * nodes,
    eslocal nElems,
    const eslocal * elems,
    double nu,
    double E,
    eslocal orderNear,
    eslocal orderFar,
    bool verbose
    );

template void bem4i::getLameSteklovPoincare< eslocal, float >(
    float * Sarray,
    eslocal nNodes,
    const float * nodes,
    eslocal nElems,
    const eslocal * elems,
    float nu,
    float E,
    eslocal orderNear,
    eslocal orderFar,
    bool verbose
    );

}

#endif
