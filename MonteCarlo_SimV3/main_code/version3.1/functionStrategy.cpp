#include "functionStrategy.h"

#include <cmath>

using namespace std;

functionStrategy::functionStrategy(vector<double>& functionParams){
	initParams(functionParams);
}

void functionStrategy::initParams(vector<double>& functionParams){
	params = functionParams;
}


double myFunc::execute(double x){
	return sin(x) + params[c0];
}

double myFunc2::execute(double x){
	return cos(x) + params[c0];
}


