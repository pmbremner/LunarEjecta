#include "solverContext.h"
#include <vector>

using namespace std;

solverContext::initParams(vector<double>& functionParams) {
	currentFunction->initParams(functionParams);
}


solverContext::initSolver(rootSolverStrategy* newSolver){
	delete currentSolver;
	currentSolver = newSolver;
}


MvsV_Context::MvsV_Context(vector<double>& functionParams) {
	currentFunction = new VvsX(functionParams);
	currentSolver   = new zbrent(currentFunction);

	function_MvsX = new MvsX(functionParams);
}

MvsV_Context::~MvsV_Context() {
	delete secondaryFunction;
	delete currentSolver;
	delete currentFunction;
}


MvsV_Context::eval(double v) {
	return function_MvsX->execute( currentSolver->solve(v, , ) );
}