#ifndef ROOTSOLVERSTRATEGY_H
#define ROOTSOLVERSTRATEGY_H


// pure virtual class, cannot intantiate base class, must do so w/ derived classes
class rootSolverStrategy
{
public:
	rootSolverStrategy(functionStrategy* specificFunction);
	virtual ~rootSolverStrategy() = default;
	
	void initFunction(functionStrategy* specificFunction);
	virtual double solve(double x_lhs, double xmin, double xmax) = 0;
	void updateMaximumIterations(int n);

protected:
	functionStrategy* function;
	int i_iter_max = 100;
	double eps = 3.E-8;
	double tol = 1.E-10;
};


class zbrent : public rootSolverStrategy
{
public:
	zbrent(functionStrategy* specificFunction) : rootSolverStrategy(specificFunction) {}
	
	double solve(double x_lhs, double xmin, double xmax);
};


#endif 