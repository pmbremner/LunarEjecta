#include "lunarEjecta_FractalIntegration.h"
#include "lunarEjecta_GeneralExpressions.h"

using namespace std;

#include <vector>
#include <cmath>
#include <iostream>

int main(int argc, char const *argv[])
{
	
	lunarEjecta_FractalIntegration integrationTest(0., 1., 0., 1., 0.05);


	integrationTest.printQuarryPointsIfEval();

	return 0;
}