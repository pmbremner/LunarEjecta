#ifndef SOLVERCONTEXT_H
#define SOLVERCONTEXT_H

#include "functionStrategy.h"
#include "rootSolverStrategy.h"
#include <vector>

using namespace std;

// pure virtual class, cannot intantiate base class, must do so w/ derived classes
class solverContext
{
public:
	//solverContext();
	virtual ~solverContext() = default;

	virtual double eval(double x) = 0;
	void initParams(vector<double>& functionParams);
	void initSolver(rootSolverStrategy* newSolver);

protected:
	functionStrategy* currentFunction;
	rootSolverStrategy* currentSolver;
};


class MvsV_Context : public solverContext
{
public:
	MvsV_Context(vector<double>& functionParams);
	~MvsV_Context();

	double eval(double v_over_U);

	static const int paramSize = 15;
	enum vars_functionParams{C1, k, p, n1, n2, mu, nu, rho, a, m, Unorm, delta, R, vmin, vmax};

private:

	vector<double>& adapterFunctionParamsVvsX(vector<double>& functionParams);
	vector<double>& adapterFunctionParamsMvsX(vector<double>& functionParams);

	functionStrategy* function_MvsX;

	vector<double> functionParams_VvsX;
	vector<double> functionParams_MvsX;

	double crater_vmin;
	double crater_vmax;


};



#endif 