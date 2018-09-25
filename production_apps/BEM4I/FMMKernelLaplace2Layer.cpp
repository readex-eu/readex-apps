/*!
 * @file    FMMKernelLaplace2Layer.cpp
 * @author  Michal Merta 
 * @date    September 3, 2013
 * 
 */

#ifdef FMM_SCALAR
  #include "Scalar/FMMKernelLaplace2LayerScalar.cpp"
#elif defined(FMM_OMP)
  #include "OMP/FMMKernelLaplace2LayerOMP.cpp"
#endif
