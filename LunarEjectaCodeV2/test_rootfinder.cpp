#include "LunarEjecta_params.h"
#include <vector>
#include <cmath>


double testfunc(double x, vector<double>& vars)
{
	return pow(x, -1./vars[0]) * pow(1. - x/vars[1], vars[2]);
}

using namespace std;

int main(int argc, char const *argv[])
{
	int N = atoi(argv[1]);

	vector<double> v, x(N, 0.0), vars;
	vars.push_back(0.4); // mu
	vars.push_back(10.); // n_2*R
	vars.push_back(0.3); // p

	double x_min = 0.1;
	double x_max = vars[1]*0.99;

	double v_U_min = testfunc(x_max, vars);
	double v_U_max = testfunc(x_min, vars);

	
	logspace(v, v_U_min, v_U_max, N, 0, N);

	

	for (int i = 0; i < N; ++i)
	{
		x[i] = findX(v[i], testfunc, x_min, vars[1], vars);
		cout << v[i] << ' ' << x[i] << endl;
	}


	return 0;
}