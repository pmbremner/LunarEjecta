#ifndef SOLVERCONTEXT_H
#define SOLVERCONTEXT_H

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


class MvsV_Context
{
public:
	MvsV_Context(vector<double>& functionParams);
	~MvsV_Context();

	double eval(double v);

private:
	functionStrategy* function_MvsX;
	
};



#endif 