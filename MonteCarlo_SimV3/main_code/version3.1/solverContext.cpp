#include "solverContext.h"
#include <vector>
#include <iostream>

using namespace std;

void solverContext::initParams(vector<double>& functionParams) {
	int N = currentFunction->getNparams();
	currentFunction->initParams(functionParams, N);
}


void solverContext::initSolver(rootSolverStrategy* newSolver){
	delete currentSolver;
	currentSolver = newSolver;
}



MvsV_Context::MvsV_Context(vector<double>& functionParams) {
	cout << "MvsV_Context::MvsV_Context\n";
	currentFunction = new VvsX(adapterFunctionParamsVvsX(functionParams));
	currentSolver   = new zbrent(currentFunction);

	function_MvsX = new MvsX(adapterFunctionParamsMvsX(functionParams));

	crater_vmin = functionParams[vmin];
	crater_vmax = functionParams[vmax];
}

MvsV_Context::~MvsV_Context() {
	cout << "MvsV_Context::~MvsV_Context\n";
	delete function_MvsX;
	delete currentSolver;
	delete currentFunction;
}


double MvsV_Context::eval(double v_over_U) {
	return function_MvsX->execute( currentSolver->solve(v_over_U, crater_vmin, crater_vmax) );
}


vector<double>& MvsV_Context::adapterFunctionParamsVvsX(vector<double>& functionParams){
	if (functionParams.size() != MvsV_Context::paramSize){
		cout << "ERROR: MvsV_Context::adapterFunctionParamsVvsX | parameter vector size does not match expected size\n";
		cout << "       functionParams.size() = " << functionParams.size() << " | Nparams = " << MvsV_Context::paramSize << endl;
	} else{
		functionParams_VvsX.resize(VvsX::paramSize, 0);

		functionParams_VvsX[VvsX::C1]     = functionParams[C1];
		functionParams_VvsX[VvsX::p]      = functionParams[p];
		functionParams_VvsX[VvsX::n1]     = functionParams[n1];
		functionParams_VvsX[VvsX::n2]     = functionParams[n2];
		functionParams_VvsX[VvsX::mu]     = functionParams[mu];
		functionParams_VvsX[VvsX::nu]     = functionParams[nu];
		functionParams_VvsX[VvsX::rho]    = functionParams[rho];
		functionParams_VvsX[VvsX::a]      = functionParams[a];
		functionParams_VvsX[VvsX::Unorm]  = functionParams[Unorm];
		functionParams_VvsX[VvsX::delta]  = functionParams[delta];
		functionParams_VvsX[VvsX::R]      = functionParams[R];
	}


	return functionParams_VvsX;
}


vector<double>& MvsV_Context::adapterFunctionParamsMvsX(vector<double>& functionParams){
	if (functionParams.size() != MvsV_Context::paramSize){
		cout << "ERROR: MvsV_Context::adapterFunctionParamsMvsX | parameter vector size does not match expected size\n";
		cout << "       functionParams.size() = " << functionParams.size() << " | Nparams = " << MvsV_Context::paramSize << endl;

	} else{
		functionParams_MvsX.resize(MvsX::paramSize, 0);

		functionParams_MvsX[MvsX::k]     = functionParams[k];
		functionParams_MvsX[MvsX::n1]    = functionParams[n1];
		functionParams_MvsX[MvsX::n2]    = functionParams[n2];
		functionParams_MvsX[MvsX::rho]   = functionParams[rho];
		functionParams_MvsX[MvsX::a]     = functionParams[a];
		functionParams_MvsX[MvsX::m]     = functionParams[m];
		functionParams_MvsX[MvsX::delta] = functionParams[delta];
		functionParams_MvsX[MvsX::R]     = functionParams[R];
	}

	return functionParams_MvsX;
}
