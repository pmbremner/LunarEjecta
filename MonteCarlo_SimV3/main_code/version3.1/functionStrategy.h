#ifndef FUNCTIONSTRATEGY_H
#define FUNCTIONSTRATEGY_H

#include <vector>

using namespace std;

// pure virtual class, cannot intantiate base class, must do so w/ derived classes
class functionStrategy
{
public:
	static const int paramSize = 0;

	functionStrategy(vector<double>& functionParams, int Nparams_i);
	virtual ~functionStrategy() = default;
	
	void initParams(vector<double>& functionParams, int Nparams_i);
	int getNparams();
	virtual double execute(double x) = 0;

protected:
	vector<double> params;
	int Nparams;
};


class myFunc : public functionStrategy
{
public:
	static const int paramSize = 1;

	myFunc(vector<double>& functionParams) : functionStrategy(functionParams, paramSize) {}

	double execute(double x);
	

	enum vars{c0};
};

class myFunc2 : public functionStrategy
{
public:
	static const int paramSize = 1;
	
	myFunc2(vector<double>& functionParams) : functionStrategy(functionParams, paramSize) {}

	double execute(double x);
	

	enum vars{c0};
};

class VvsX : public functionStrategy
{
public:
	static const int paramSize = 11;

	VvsX(vector<double>& functionParams) : functionStrategy(functionParams, paramSize) {}

	double execute(double x);

	enum vars{C1, p, n1, n2, mu, nu, rho, a, Unorm, delta, R};

};

class MvsX : public functionStrategy
{
public:
	static const int paramSize = 8; 

	MvsX(vector<double>& functionParams) : functionStrategy(functionParams, paramSize) {}

	double execute(double x);

	enum vars{k, n1, n2, rho, a, m, delta, R};

};



#endif 