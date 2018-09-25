/*!
 * @file    FMMKernelLaplace1Layer.cpp
 * @author  Michal Merta
 * @date    September 23, 2013
 * 
 */

#ifdef FMM_SCALAR
  #include "Scalar/FMMKernelLaplace1LayerScalar.cpp"
#elif defined(FMM_OMP)
  #include "OMP/FMMKernelLaplace1LayerOMP.cpp"
#endif
