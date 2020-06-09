#include "lunarEjecta_FractalIntegration.h"
#include "lunarEjecta_AdaptiveMesh.h"
#include "lunarEjecta_GeneralExpressions.h"

using namespace std;

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

int main(int argc, char const *argv[])
{

	vector<double> x, y;
	vector<vector<double>> z;

	int Nx = 100, Ny = 100;
	int maxLevelMesh = 3;     // 3, 1 seems to be the best combo
	int maxLevelFractal = 1;

	double D0 = 0.001;
	double D1 = 0.0015;

	linspace(x, 0., 1., Nx+1); // bin edges
	linspace(y, 0., 1., Ny+1); // bin edges
	
	// z is bin centers
	z.resize(Nx);
	for (int i = 0; i < Nx; ++i)
		z[i].resize(Ny);

	// lunarEjecta_FractalIntegration integrationTest(0., 1., 0., 1., 0.01, 0.02, 0.0001);
	// integrationTest.printQuarryPointsIfEval();

	lunarEjecta_AdaptiveMesh integrationScheme(Nx, Ny, maxLevelMesh, maxLevelFractal);

	integrationScheme.evalBin(D0, D1);

	integrationScheme.printDataToFile("Test.txt");

	cout << "***************************\n";
	cout << "Reduced Integral = " << setprecision(16) <<
		integrationScheme.getReducedIntegral()
		<< endl;

	cout << " maxLevelMesh = " << maxLevelMesh << endl;
	cout << " maxLevelFractal = " << maxLevelFractal << endl;
	integrationScheme.printEvalCounts();

	// cout << "******************\n";
	// cout << " Total integration points = " << fractIntCount << endl;


	return 0;
}