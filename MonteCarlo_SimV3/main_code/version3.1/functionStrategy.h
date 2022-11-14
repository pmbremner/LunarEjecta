#ifndef FUNCTIONSTRATEGY_H
#define FUNCTIONSTRATEGY_H




class dataSet
{

};


class myData : public dataSet
{
public:
	myData(double ai, double bi);

	double a;
	double b;
	
};

class functionStrategy
{

public:
	dataSet* params;


public:
	virtual double execute(double input);
	



};



class testFunc : public functionStrategy
{

public:

	testFunc(dataSet* testData);

	double execute(double input);
	
};






#endif 