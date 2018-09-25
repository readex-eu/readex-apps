#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=icc
CCC=icpc
CXX=icpc
FC=ifort
AS=as

# Macros
CND_PLATFORM=Intel-Linux
CND_DLIB_EXT=so
CND_CONF=debug_intel
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/ACAMatrix.o \
	${OBJECTDIR}/BEBilinearForm.o \
	${OBJECTDIR}/BEBilinearFormHelmholtz1Layer.o \
	${OBJECTDIR}/BEBilinearFormHelmholtz2Layer.o \
	${OBJECTDIR}/BEBilinearFormHelmholtzHypersingular.o \
	${OBJECTDIR}/BEBilinearFormLame1Layer.o \
	${OBJECTDIR}/BEBilinearFormLame2Layer.o \
	${OBJECTDIR}/BEBilinearFormLameHypersingular.o \
	${OBJECTDIR}/BEBilinearFormLaplace1Layer.o \
	${OBJECTDIR}/BEBilinearFormLaplace2Layer.o \
	${OBJECTDIR}/BEBilinearFormLaplaceHypersingular.o \
	${OBJECTDIR}/BEBilinearFormWave1Layer.o \
	${OBJECTDIR}/BEBilinearFormWaveHypersingular.o \
	${OBJECTDIR}/BEIntegrator.o \
	${OBJECTDIR}/BEIntegratorHelmholtz.o \
	${OBJECTDIR}/BEIntegratorLame.o \
	${OBJECTDIR}/BEIntegratorLaplace.o \
	${OBJECTDIR}/BEIntegratorWave.o \
	${OBJECTDIR}/BESpace.o \
	${OBJECTDIR}/BESpaceTime.o \
	${OBJECTDIR}/BernoulliSubproblem.o \
	${OBJECTDIR}/BlockLinearOperator.o \
	${OBJECTDIR}/BlockMatrix.o \
	${OBJECTDIR}/CompoundLinearOperator.o \
	${OBJECTDIR}/DistributedFastBESpace.o \
	${OBJECTDIR}/FMMKernel.o \
	${OBJECTDIR}/FMMKernelLaplace1Layer.o \
	${OBJECTDIR}/FMMKernelLaplace2Layer.o \
	${OBJECTDIR}/FMMMatrix.o \
	${OBJECTDIR}/FastBESpace.o \
	${OBJECTDIR}/FixedTrackingSubproblem.o \
	${OBJECTDIR}/FullMatrix.o \
	${OBJECTDIR}/HeatSourceSubproblem.o \
	${OBJECTDIR}/HelmholtzHypersingularP0P0Operator.o \
	${OBJECTDIR}/HelmholtzNeumannRobinOperator.o \
	${OBJECTDIR}/HelmholtzRegularizedExteriorDirichletOperator.o \
	${OBJECTDIR}/HelmholtzRegularizedExteriorNeumannOperator.o \
	${OBJECTDIR}/HomogenizationProblem.o \
	${OBJECTDIR}/IdentityOperator.o \
	${OBJECTDIR}/InverseLinearOperator.o \
	${OBJECTDIR}/Lame1LayerP0P0MultilvlPrecond.o \
	${OBJECTDIR}/Laplace1LayerP0P0MultilvlPrecond.o \
	${OBJECTDIR}/LaplaceHypersingularOperator.o \
	${OBJECTDIR}/LaplaceSteklovPoincareOperator.o \
	${OBJECTDIR}/MPIACAMatrix.o \
	${OBJECTDIR}/MPIBlockACAMatrix.o \
	${OBJECTDIR}/MPIBlockMatrix.o \
	${OBJECTDIR}/MathFun.o \
	${OBJECTDIR}/Matrix.o \
	${OBJECTDIR}/Mesh.o \
	${OBJECTDIR}/MultiresolutionOptimizer.o \
	${OBJECTDIR}/OMP/FMMKernelLaplace1LayerOMP.o \
	${OBJECTDIR}/OMP/FMMKernelLaplace2LayerOMP.o \
	${OBJECTDIR}/OpenMeshWrapper.o \
	${OBJECTDIR}/PotentialsHelmholtz.o \
	${OBJECTDIR}/PotentialsLaplace.o \
	${OBJECTDIR}/PotentialsWave.o \
	${OBJECTDIR}/RepresentationFormula.o \
	${OBJECTDIR}/RepresentationFormulaHelmholtz.o \
	${OBJECTDIR}/RepresentationFormulaLaplace.o \
	${OBJECTDIR}/STFOperator.o \
	${OBJECTDIR}/STFPreconditioner.o \
	${OBJECTDIR}/Scalar/FMMKernelLaplace1LayerScalar.o \
	${OBJECTDIR}/Scalar/FMMKernelLaplace2LayerScalar.o \
	${OBJECTDIR}/SparseMatrix.o \
	${OBJECTDIR}/SumLinearOperator.o \
	${OBJECTDIR}/SurfaceMesh3D.o \
	${OBJECTDIR}/Tree.o \
	${OBJECTDIR}/Vector.o \
	${OBJECTDIR}/VolumeMesh3D.o \
	${OBJECTDIR}/WavePreconditioner.o \
	${OBJECTDIR}/WaveScatteringProblem.o \
	${OBJECTDIR}/main.o


# C Compiler Flags
CFLAGS=-m64

# CC Compiler Flags
CCFLAGS=-m64 -xHost -qopenmp -O0 -w2
CXXFLAGS=-m64 -xHost -qopenmp -O0 -w2

# Fortran Compiler Flags
FFLAGS=-m64

# Assembler Flags
ASFLAGS=--64

# Link Libraries and Options
LDLIBSOPTIONS=-L/opt/intel/compilers_and_libraries/linux/mpi/intel64/lib -L/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64 -L/opt/intel/compilers_and_libraries/linux/lib/intel64 -L/opt/nlopt/lib -L/opt/OpenMesh/lib -Wl,-rpath,'/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/' -Wl,-rpath,'/opt/intel/compilers_and_libraries/linux/mpi/intel64/lib' -Wl,-rpath,'/opt/intel/compilers_and_libraries/linux/lib/intel64/' -Wl,-rpath,'/opt/nlopt/lib' -Wl,-rpath,'/opt/OpenMesh/lib' -lmpi -lmpicxx -lnlopt_cxx -lOpenMeshCore -lOpenMeshTools -liomp5 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/bem4i

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/bem4i: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/bem4i ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/ACAMatrix.o: ACAMatrix.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ACAMatrix.o ACAMatrix.cpp

${OBJECTDIR}/BEBilinearForm.o: BEBilinearForm.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEBilinearForm.o BEBilinearForm.cpp

${OBJECTDIR}/BEBilinearFormHelmholtz1Layer.o: BEBilinearFormHelmholtz1Layer.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEBilinearFormHelmholtz1Layer.o BEBilinearFormHelmholtz1Layer.cpp

${OBJECTDIR}/BEBilinearFormHelmholtz2Layer.o: BEBilinearFormHelmholtz2Layer.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEBilinearFormHelmholtz2Layer.o BEBilinearFormHelmholtz2Layer.cpp

${OBJECTDIR}/BEBilinearFormHelmholtzHypersingular.o: BEBilinearFormHelmholtzHypersingular.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEBilinearFormHelmholtzHypersingular.o BEBilinearFormHelmholtzHypersingular.cpp

${OBJECTDIR}/BEBilinearFormLame1Layer.o: BEBilinearFormLame1Layer.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEBilinearFormLame1Layer.o BEBilinearFormLame1Layer.cpp

${OBJECTDIR}/BEBilinearFormLame2Layer.o: BEBilinearFormLame2Layer.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEBilinearFormLame2Layer.o BEBilinearFormLame2Layer.cpp

${OBJECTDIR}/BEBilinearFormLameHypersingular.o: BEBilinearFormLameHypersingular.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEBilinearFormLameHypersingular.o BEBilinearFormLameHypersingular.cpp

${OBJECTDIR}/BEBilinearFormLaplace1Layer.o: BEBilinearFormLaplace1Layer.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEBilinearFormLaplace1Layer.o BEBilinearFormLaplace1Layer.cpp

${OBJECTDIR}/BEBilinearFormLaplace2Layer.o: BEBilinearFormLaplace2Layer.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEBilinearFormLaplace2Layer.o BEBilinearFormLaplace2Layer.cpp

${OBJECTDIR}/BEBilinearFormLaplaceHypersingular.o: BEBilinearFormLaplaceHypersingular.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEBilinearFormLaplaceHypersingular.o BEBilinearFormLaplaceHypersingular.cpp

${OBJECTDIR}/BEBilinearFormWave1Layer.o: BEBilinearFormWave1Layer.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEBilinearFormWave1Layer.o BEBilinearFormWave1Layer.cpp

${OBJECTDIR}/BEBilinearFormWaveHypersingular.o: BEBilinearFormWaveHypersingular.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEBilinearFormWaveHypersingular.o BEBilinearFormWaveHypersingular.cpp

${OBJECTDIR}/BEIntegrator.o: BEIntegrator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEIntegrator.o BEIntegrator.cpp

${OBJECTDIR}/BEIntegratorHelmholtz.o: BEIntegratorHelmholtz.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEIntegratorHelmholtz.o BEIntegratorHelmholtz.cpp

${OBJECTDIR}/BEIntegratorLame.o: BEIntegratorLame.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEIntegratorLame.o BEIntegratorLame.cpp

${OBJECTDIR}/BEIntegratorLaplace.o: BEIntegratorLaplace.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEIntegratorLaplace.o BEIntegratorLaplace.cpp

${OBJECTDIR}/BEIntegratorWave.o: BEIntegratorWave.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BEIntegratorWave.o BEIntegratorWave.cpp

${OBJECTDIR}/BESpace.o: BESpace.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BESpace.o BESpace.cpp

${OBJECTDIR}/BESpaceTime.o: BESpaceTime.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BESpaceTime.o BESpaceTime.cpp

${OBJECTDIR}/BernoulliSubproblem.o: BernoulliSubproblem.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BernoulliSubproblem.o BernoulliSubproblem.cpp

${OBJECTDIR}/BlockLinearOperator.o: BlockLinearOperator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BlockLinearOperator.o BlockLinearOperator.cpp

${OBJECTDIR}/BlockMatrix.o: BlockMatrix.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BlockMatrix.o BlockMatrix.cpp

${OBJECTDIR}/CompoundLinearOperator.o: CompoundLinearOperator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CompoundLinearOperator.o CompoundLinearOperator.cpp

${OBJECTDIR}/DistributedFastBESpace.o: DistributedFastBESpace.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/DistributedFastBESpace.o DistributedFastBESpace.cpp

${OBJECTDIR}/FMMKernel.o: FMMKernel.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FMMKernel.o FMMKernel.cpp

${OBJECTDIR}/FMMKernelLaplace1Layer.o: FMMKernelLaplace1Layer.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FMMKernelLaplace1Layer.o FMMKernelLaplace1Layer.cpp

${OBJECTDIR}/FMMKernelLaplace2Layer.o: FMMKernelLaplace2Layer.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FMMKernelLaplace2Layer.o FMMKernelLaplace2Layer.cpp

${OBJECTDIR}/FMMMatrix.o: FMMMatrix.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FMMMatrix.o FMMMatrix.cpp

${OBJECTDIR}/FastBESpace.o: FastBESpace.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FastBESpace.o FastBESpace.cpp

${OBJECTDIR}/FixedTrackingSubproblem.o: FixedTrackingSubproblem.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FixedTrackingSubproblem.o FixedTrackingSubproblem.cpp

${OBJECTDIR}/FullMatrix.o: FullMatrix.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FullMatrix.o FullMatrix.cpp

${OBJECTDIR}/HeatSourceSubproblem.o: HeatSourceSubproblem.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/HeatSourceSubproblem.o HeatSourceSubproblem.cpp

${OBJECTDIR}/HelmholtzHypersingularP0P0Operator.o: HelmholtzHypersingularP0P0Operator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/HelmholtzHypersingularP0P0Operator.o HelmholtzHypersingularP0P0Operator.cpp

${OBJECTDIR}/HelmholtzNeumannRobinOperator.o: HelmholtzNeumannRobinOperator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/HelmholtzNeumannRobinOperator.o HelmholtzNeumannRobinOperator.cpp

${OBJECTDIR}/HelmholtzRegularizedExteriorDirichletOperator.o: HelmholtzRegularizedExteriorDirichletOperator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/HelmholtzRegularizedExteriorDirichletOperator.o HelmholtzRegularizedExteriorDirichletOperator.cpp

${OBJECTDIR}/HelmholtzRegularizedExteriorNeumannOperator.o: HelmholtzRegularizedExteriorNeumannOperator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/HelmholtzRegularizedExteriorNeumannOperator.o HelmholtzRegularizedExteriorNeumannOperator.cpp

${OBJECTDIR}/HomogenizationProblem.o: HomogenizationProblem.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/HomogenizationProblem.o HomogenizationProblem.cpp

${OBJECTDIR}/IdentityOperator.o: IdentityOperator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/IdentityOperator.o IdentityOperator.cpp

${OBJECTDIR}/InverseLinearOperator.o: InverseLinearOperator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/InverseLinearOperator.o InverseLinearOperator.cpp

${OBJECTDIR}/Lame1LayerP0P0MultilvlPrecond.o: Lame1LayerP0P0MultilvlPrecond.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Lame1LayerP0P0MultilvlPrecond.o Lame1LayerP0P0MultilvlPrecond.cpp

${OBJECTDIR}/Laplace1LayerP0P0MultilvlPrecond.o: Laplace1LayerP0P0MultilvlPrecond.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Laplace1LayerP0P0MultilvlPrecond.o Laplace1LayerP0P0MultilvlPrecond.cpp

${OBJECTDIR}/LaplaceHypersingularOperator.o: LaplaceHypersingularOperator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/LaplaceHypersingularOperator.o LaplaceHypersingularOperator.cpp

${OBJECTDIR}/LaplaceSteklovPoincareOperator.o: LaplaceSteklovPoincareOperator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/LaplaceSteklovPoincareOperator.o LaplaceSteklovPoincareOperator.cpp

${OBJECTDIR}/MPIACAMatrix.o: MPIACAMatrix.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MPIACAMatrix.o MPIACAMatrix.cpp

${OBJECTDIR}/MPIBlockACAMatrix.o: MPIBlockACAMatrix.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MPIBlockACAMatrix.o MPIBlockACAMatrix.cpp

${OBJECTDIR}/MPIBlockMatrix.o: MPIBlockMatrix.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MPIBlockMatrix.o MPIBlockMatrix.cpp

${OBJECTDIR}/MathFun.o: MathFun.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MathFun.o MathFun.cpp

${OBJECTDIR}/Matrix.o: Matrix.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Matrix.o Matrix.cpp

${OBJECTDIR}/Mesh.o: Mesh.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Mesh.o Mesh.cpp

${OBJECTDIR}/MultiresolutionOptimizer.o: MultiresolutionOptimizer.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MultiresolutionOptimizer.o MultiresolutionOptimizer.cpp

${OBJECTDIR}/OMP/FMMKernelLaplace1LayerOMP.o: OMP/FMMKernelLaplace1LayerOMP.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/OMP
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/OMP/FMMKernelLaplace1LayerOMP.o OMP/FMMKernelLaplace1LayerOMP.cpp

${OBJECTDIR}/OMP/FMMKernelLaplace2LayerOMP.o: OMP/FMMKernelLaplace2LayerOMP.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/OMP
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/OMP/FMMKernelLaplace2LayerOMP.o OMP/FMMKernelLaplace2LayerOMP.cpp

${OBJECTDIR}/OpenMeshWrapper.o: OpenMeshWrapper.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/OpenMeshWrapper.o OpenMeshWrapper.cpp

${OBJECTDIR}/PotentialsHelmholtz.o: PotentialsHelmholtz.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/PotentialsHelmholtz.o PotentialsHelmholtz.cpp

${OBJECTDIR}/PotentialsLaplace.o: PotentialsLaplace.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/PotentialsLaplace.o PotentialsLaplace.cpp

${OBJECTDIR}/PotentialsWave.o: PotentialsWave.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/PotentialsWave.o PotentialsWave.cpp

${OBJECTDIR}/RepresentationFormula.o: RepresentationFormula.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/RepresentationFormula.o RepresentationFormula.cpp

${OBJECTDIR}/RepresentationFormulaHelmholtz.o: RepresentationFormulaHelmholtz.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/RepresentationFormulaHelmholtz.o RepresentationFormulaHelmholtz.cpp

${OBJECTDIR}/RepresentationFormulaLaplace.o: RepresentationFormulaLaplace.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/RepresentationFormulaLaplace.o RepresentationFormulaLaplace.cpp

${OBJECTDIR}/STFOperator.o: STFOperator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/STFOperator.o STFOperator.cpp

${OBJECTDIR}/STFPreconditioner.o: STFPreconditioner.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/STFPreconditioner.o STFPreconditioner.cpp

${OBJECTDIR}/Scalar/FMMKernelLaplace1LayerScalar.o: Scalar/FMMKernelLaplace1LayerScalar.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/Scalar
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Scalar/FMMKernelLaplace1LayerScalar.o Scalar/FMMKernelLaplace1LayerScalar.cpp

${OBJECTDIR}/Scalar/FMMKernelLaplace2LayerScalar.o: Scalar/FMMKernelLaplace2LayerScalar.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}/Scalar
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Scalar/FMMKernelLaplace2LayerScalar.o Scalar/FMMKernelLaplace2LayerScalar.cpp

${OBJECTDIR}/SparseMatrix.o: SparseMatrix.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SparseMatrix.o SparseMatrix.cpp

${OBJECTDIR}/SumLinearOperator.o: SumLinearOperator.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SumLinearOperator.o SumLinearOperator.cpp

${OBJECTDIR}/SurfaceMesh3D.o: SurfaceMesh3D.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SurfaceMesh3D.o SurfaceMesh3D.cpp

${OBJECTDIR}/Tree.o: Tree.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Tree.o Tree.cpp

${OBJECTDIR}/Vector.o: Vector.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Vector.o Vector.cpp

${OBJECTDIR}/VolumeMesh3D.o: VolumeMesh3D.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/VolumeMesh3D.o VolumeMesh3D.cpp

${OBJECTDIR}/WavePreconditioner.o: WavePreconditioner.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/WavePreconditioner.o WavePreconditioner.cpp

${OBJECTDIR}/WaveScatteringProblem.o: WaveScatteringProblem.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/WaveScatteringProblem.o WaveScatteringProblem.cpp

${OBJECTDIR}/main.o: main.cpp nbproject/Makefile-${CND_CONF}.mk
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/opt/intel/compilers_and_libraries/linux/mpi/intel64/include/ -I/opt/nlopt/include -I/opt/OpenMesh/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
