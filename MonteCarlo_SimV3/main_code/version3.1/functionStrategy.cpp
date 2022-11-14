#include "functionStrategy.h"

using namespace std;


double testFunc::execute(double input){
	return params->a * input*input + params->b;
}

testFunc::testFunc(dataSet* testData){
	//cout << "init in testFunc" << endl;
	params = testData;
}


myData::myData(double ai, double bi){
	a = ai;
	b = bi;
}