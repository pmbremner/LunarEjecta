#include "functionStrategy.h"

#include <cmath>
#include <iostream>

using namespace std;

functionStrategy::functionStrategy(vector<double>& functionParams, int Nparams_i){
	initParams(functionParams, Nparams_i);
}

void functionStrategy::initParams(vector<double>& functionParams, int Nparams_i){
	if (functionParams.size() != Nparams_i){
		cout << "ERROR: functionStrategy constructor | parameter vector size does not match expected size\n";
		Nparams = 0;
	} else {
		Nparams = Nparams_i;
		params  = functionParams;
	}
}


double myFunc::execute(double x){
	return sin(x) + params[c0];
}

double myFunc2::execute(double x){
	return cos(x) + params[c0];
}


