#include "functionStrategy.h"

#include <cmath>
#include <iostream>

using namespace std;

functionStrategy::functionStrategy(vector<double>& functionParams, int Nparams_i){
	initParams(functionParams, Nparams_i);
}

void functionStrategy::initParams(vector<double>& functionParams, int Nparams_i){
	cout << "functionStrategy::initParams,";
	cout << " functionParams.size() = " << functionParams.size();
	cout << " | Nparams = " << Nparams_i << endl;
	if (functionParams.size() != Nparams_i){
		cout << "ERROR: functionStrategy constructor | parameter vector size does not match expected size\n";
		Nparams = 0;
	} else {
		Nparams = Nparams_i;
		params  = functionParams;
	}
}

int functionStrategy::getNparams(){
	return Nparams;
}


double myFunc::execute(double x){
	return sin(x) + params[c0];
}

double myFunc2::execute(double x){
	return cos(x) + params[c0];
}



double VvsX::execute(double x){
	return params[C1]
		* pow((x / params[a]) * pow(params[rho] / params[delta], params[nu]), -1./params[mu])
		* pow(1. - x / (params[n2] * params[R]), params[p]);
}


double MvsX::execute(double x){
	return (3. * params[k] / (4.* M_PI))
		* (params[rho] / params[delta])
		* (pow(x / params[a], 3) - pow(params[n1], 3));
}
