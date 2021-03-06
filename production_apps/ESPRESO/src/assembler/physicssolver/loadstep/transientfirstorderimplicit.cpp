
#include "transientfirstorderimplicit.h"

#include "../assembler.h"
#include "../timestep/timestepsolver.h"
#include "../../step.h"
#include "../../instance.h"
#include "../../solution.h"
#include "../../../mesh/structures/elementtypes.h"
#include "../../../configuration/physics/transientsolver.h"

#include "../../../solver/generic/FETISolver.h"

#ifdef USE_ATP
#include "atplib.h"
#endif

using namespace espreso;

size_t TransientFirstOrderImplicit::loadStep = 0;
std::vector<Solution*> TransientFirstOrderImplicit::solutions;

TransientFirstOrderImplicit::TransientFirstOrderImplicit(TimeStepSolver &timeStepSolver, const TransientSolver &configuration, double duration)
: LoadStepSolver("TRANSIENT", timeStepSolver, duration), _configuration(configuration), _alpha(0)
{
	if (configuration.time_step < 1e-7) {
		ESINFO(GLOBAL_ERROR) << "Set time step for TRANSIENT solver greater than 1e-7.";
	}
}

Matrices TransientFirstOrderImplicit::updateStructuralMatrices(Step &step, Matrices matrices)
{
	Matrices updatedMatrices = matrices & (Matrices::K | Matrices::M | Matrices::f | Matrices::R | Matrices::B1 | Matrices::B1c | Matrices::B1duplicity);

	if (step.substep && !_timeDependent && !_tempDependent) {
		updatedMatrices &= (Matrices::f | Matrices::B1c);
	}

	return reassembleStructuralMatrices(step, updatedMatrices);
}

Matrices TransientFirstOrderImplicit::reassembleStructuralMatrices(Step &step, Matrices matrices)
{
	_assembler.updateMatrices(step, matrices);
	if (matrices & (Matrices::K | Matrices::M)) {
		_assembler.sum(
				_assembler.instance.K,
				1 / (_alpha * step.timeStep), _assembler.instance.M,
				"K += (1 / alpha * delta T) * M");
	}

	if (matrices & (Matrices::K | Matrices::M | Matrices::f)) {
		_assembler.sum(
				solutions[SolutionIndex::X]->data,
				1 / (_alpha * step.timeStep), solutions[SolutionIndex::U]->data,
				(1 - _alpha) / _alpha, solutions[SolutionIndex::V]->data,
				"x = (1 / alpha * delta T) * U + (1 - alpha) / alpha * V");

		_assembler.multiply(
				solutions[SolutionIndex::Y]->data,
				_assembler.instance.M, solutions[SolutionIndex::X]->data,
				"y = M * x");

		_assembler.sum(_assembler.instance.f,
				1, _assembler.instance.f,
				1, solutions[SolutionIndex::Y]->data,
				"f += y");
	}

	return matrices;
}

void TransientFirstOrderImplicit::initLoadStep(Step &step)
{
	LoadStepSolver::initLoadStep(step);

	_assembler.setEmptyRegularizationCallback();
	_assembler.setB0Callback();

	switch (_configuration.method) {
	case TransientSolver::METHOD::CRANK_NICOLSON:
		_alpha = 0.5;
		break;
	case TransientSolver::METHOD::GALERKIN:
		_alpha = 2 / 3;
		break;
	case TransientSolver::METHOD::BACKWARD_DIFF:
		_alpha = 1;
		break;
	case TransientSolver::METHOD::USER:
		_alpha = _configuration.alpha;
		if (_alpha <= 0 || _alpha > 1) {
			ESINFO(GLOBAL_ERROR) << "Alpha has to be from interval (0, 1>.";
		}
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not supported first order implicit solver method.";
	}

	if (!solutions.size()) {
		solutions.push_back(_assembler.addSolution("trans_U", ElementType::NODES));
		solutions.push_back(_assembler.addSolution("trans_dU", ElementType::NODES));
		solutions.push_back(_assembler.addSolution("trans_V", ElementType::NODES));
		solutions.push_back(_assembler.addSolution("trans_X", ElementType::NODES));
		solutions.push_back(_assembler.addSolution("trans_Y", ElementType::NODES));
	}
	if (loadStep + 1 != step.step) {
		solutions[SolutionIndex::V]->fill(0);
	}
	loadStep = step.step;
	solutions[SolutionIndex::U]->data = _assembler.instance.primalSolution;
}

void TransientFirstOrderImplicit::runNextTimeStep(Step &step)
{
	double last = step.currentTime;
	step.currentTime += _configuration.time_step;
	if (step.currentTime + _precision >= _startTime + _duration) {
		step.currentTime = _startTime + _duration;
	}
	step.timeStep = step.currentTime - last;

/*
// READEX ATP code
#ifdef USE_ATP
	int decomposition = step.substep % 8;
	ATP_PARAM_DECLARE("DECOMP",ATP_PARAM_TYPE_RANGE, 0, "DOM2");
	int decomposition_values[3] = {0,7,1};
	ATP_PARAM_ADD_VALUES("DECOMP", decomposition_values, 3, "DOM2");
	ATP_PARAM_GET("DECOMP", &decomposition, "DOM2");
	std::vector<std::vector<eslocal> > previousDOFMapping, previousDomainMap;
	// domain size in rage <12 .. 1536> == <2^0*12 .. 2^7*12>
	_assembler.changeMeshPartition(12*(1 << decomposition), previousDOFMapping, previousDomainMap);
	_assembler.transformDomainsData(previousDOFMapping, previousDomainMap, solutions[SolutionIndex::U]->data);
	_assembler.transformDomainsData(previousDOFMapping, previousDomainMap, solutions[SolutionIndex::dU]->data);
	_assembler.transformDomainsData(previousDOFMapping, previousDomainMap, solutions[SolutionIndex::V]->data);
	_assembler.transformDomainsData(previousDOFMapping, previousDomainMap, solutions[SolutionIndex::X]->data);
	_assembler.transformDomainsData(previousDOFMapping, previousDomainMap, solutions[SolutionIndex::Y]->data);
#endif
// END OF CODE FOR CHANGE DECOMPOSITION
*/

// READEX ATP code
#ifdef USE_ATP
	int preconditioner = step.substep % 4;
	ATP_PARAM_DECLARE("PRECOND",ATP_PARAM_TYPE_RANGE, 0, "DOM1");
	int precond_values[3] = {0,3,1};
	ATP_PARAM_ADD_VALUES("PRECOND", precond_values, 3, "DOM1");
	ATP_PARAM_GET("PRECOND", &preconditioner, "DOM1");
dynamic_cast<FETISolver&>(_assembler.linearSolver).configuration.preconditioner = static_cast<ESPRESO_PRECONDITIONER>(preconditioner);
#endif
// END OF CODE FOR CHANGE PRECONDITIONER

/*
// READEX ATP code
#ifdef USE_ATP
	int solver = step.substep % 5;
	ATP_PARAM_DECLARE("SOLVER",ATP_PARAM_TYPE_RANGE, 0, "DOM3");
	int solver_values[3] = {0,4,1};
	ATP_PARAM_ADD_VALUES("SOLVER", solver_values, 3, "DOM3");
	ATP_PARAM_GET("SOLVER", &solver, "DOM3");
	dynamic_cast<FETISolver&>(_assembler.linearSolver).configuration.solver = static_cast<ESPRESO_ITERATIVE_SOLVER>(solver);
#endif
// END OF CODE FOR CHANGE SOLVER
*/
	processTimeStep(step);
}

void TransientFirstOrderImplicit::processTimeStep(Step &step)
{
	step.internalForceReduction = 1;
	step.timeIntegrationConstantK = 1;
	step.timeIntegrationConstantM = 1 / (_alpha * step.timeStep);

	_timeStepSolver.solve(step, *this);

	_assembler.sum(
			solutions[SolutionIndex::dU]->data,
			1, _assembler.instance.primalSolution,
			-1, solutions[SolutionIndex::U]->data,
			"delta U = U_i - U_i_1");

	_assembler.sum(
			solutions[SolutionIndex::V]->data,
			1 / (_alpha * step.timeStep), solutions[SolutionIndex::dU]->data,
			- (1 - _alpha) / _alpha, solutions[SolutionIndex::V]->data,
			"V = (1 / alpha * delta T) * delta U - (1 - alpha) / alpha * V");

	solutions[SolutionIndex::U]->data = _assembler.instance.primalSolution;

	_assembler.processSolution(step);
	_assembler.storeSolution(step);
}


