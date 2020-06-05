#include "lunarEjecta_FractalIntegration.h"
#include "lunarEjecta_AdaptiveMesh.h"
#include "lunarEjecta_GeneralExpressions.h"

using namespace std;

#include <vector>
#include <cmath>
#include <iostream>

int main(int argc, char const *argv[])
{

	vector<double> x, y;
	vector<vector<double>> z;

	int Nx = 200, Ny = 200;

	linspace(x, 0., 1., Nx+1); // bin edges
	linspace(y, 0., 1., Ny+1); // bin edges
	
	// z is bin centers
	z.resize(Nx);
	for (int i = 0; i < Nx; ++i)
		z[i].resize(Ny);

	// lunarEjecta_FractalIntegration integrationTest(0., 1., 0., 1., 0.01, 0.02, 0.0001);
	// integrationTest.printQuarryPointsIfEval();

	lunarEjecta_AdaptiveMesh integrationScheme(x, y, z, 1);

	integrationScheme.evalBin(0.05, 0.3);

	integrationScheme.printDataToFile("Test.txt");

	return 0;
}