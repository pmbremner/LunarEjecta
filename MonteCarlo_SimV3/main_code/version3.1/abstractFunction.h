#ifndef ABSTRACTFUNCTION_H
#define ABSTRACTFUNCTION_H

#include "crater.h"

using namespace std;

// pure virtual, need to use derived class
class abstractFunction
{
public:
	abstractFunction();
	~abstractFunction();

	virtual double getYfromX(double x) = 0;
	double getXfromY(double y);

protected:
	double rootSolver(double x);

	double xmin;
	double xmax;
};



class VvsX: public abstractFunction
{
public:
	abstractFunction(crater* ctr);
	~abstractFunction();

	void initializeAbstractFunction(crater* ctr);
	
	double getYfromX(double x);
	//double getXfromY(double y);

	double getMfromV(double v);

private:
	double n2;
	double n1;
	double k;
	double C1;
	double mu;
	double nu;
	double p;
	double rho;
	double a;
	double m;
	double delta;
	double Unorm;

};



#endif 