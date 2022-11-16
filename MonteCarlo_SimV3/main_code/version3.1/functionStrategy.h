#ifndef FUNCTIONSTRATEGY_H
#define FUNCTIONSTRATEGY_H

#include <vector>

using namespace std;

// pure virtual class, cannot intantiate base class, must do so w/ derived classes
class functionStrategy
{
public:
	functionStrategy(vector<double>& functionParams);
	virtual ~functionStrategy() = default;
	
	void initParams(vector<double>& functionParams);
	virtual double execute(double x) = 0;

protected:
	vector<double> params;
};


class myFunc : public functionStrategy
{
public:
	myFunc(vector<double>& functionParams) : functionStrategy(functionParams) {}

	double execute(double x);
	

	enum vars{c0, c1, c2};
};

class myFunc2 : public functionStrategy
{
public:
	myFunc2(vector<double>& functionParams) : functionStrategy(functionParams) {}

	double execute(double x);
	

	enum vars{c0, c1, c2};
};

class VvsX : public functionStrategy
{
public:
	VvsX(vector<double>& functionParams) : functionStrategy(functionParams) {}

	double execute(double v_over_U);

};

class MvsX : public functionStrategy
{
public:
	MvsX(vector<double>& functionParams) : functionStrategy(functionParams) {}

	double execute(double v_over_U);

};



#endif 