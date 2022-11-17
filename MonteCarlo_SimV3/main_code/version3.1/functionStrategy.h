#ifndef FUNCTIONSTRATEGY_H
#define FUNCTIONSTRATEGY_H

#include <vector>

using namespace std;

// pure virtual class, cannot intantiate base class, must do so w/ derived classes
class functionStrategy
{
public:
	functionStrategy(vector<double>& functionParams, int Nparams_i);
	virtual ~functionStrategy() = default;
	
	void initParams(vector<double>& functionParams, int Nparams_i);
	virtual double execute(double x) = 0;

protected:
	vector<double> params;
	int Nparams;
};


class myFunc : public functionStrategy
{
public:
	myFunc(vector<double>& functionParams) : functionStrategy(functionParams, 1) {}

	double execute(double x);
	

	enum vars{c0};
};

class myFunc2 : public functionStrategy
{
public:
	myFunc2(vector<double>& functionParams) : functionStrategy(functionParams, 1) {}

	double execute(double x);
	

	enum vars{c0};
};

class VvsX : public functionStrategy
{
public:
	VvsX(vector<double>& functionParams) : functionStrategy(functionParams, 11) {}

	double execute(double v_over_U);

	//enum vars{C1, k, p, n1, n2, mu, nu, rho, a, m, Unorm, delta, R};
	enum vars{C1, p, n1, n2, mu, nu, rho, a, Unorm, delta, R};

};

class MvsX : public functionStrategy
{
public:
	MvsX(vector<double>& functionParams) : functionStrategy(functionParams, 8) {}

	double execute(double v_over_U);


	//enum vars{C1, k, p, n1, n2, mu, nu, rho, a, m, Unorm, delta, R};
	enum vars{k, n1, n2, rho, a, m, delta, R};

};



#endif 