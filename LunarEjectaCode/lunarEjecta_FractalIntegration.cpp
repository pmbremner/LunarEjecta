#include "lunarEjecta_GeneralExpressions.h"
#include "lunarEjecta_FractalIntegration.h"

using namespace std;

#include <vector>
#include <cmath>
#include <iostream>

lunarEjecta_FractalIntegration::lunarEjecta_FractalIntegration
		(double new_xMin,
		 double new_xMax,
		 double new_yMin,
		 double new_yMax,
		 double new_epsError)
{
	xMin = new_xMin;
	xMax = new_xMax;
	yMin = new_yMin;
	yMax = new_yMax;

	epsError = new_epsError;

	levelMin = 2;
	levelMax = 10;
	levelCur = -1;

	// init size of quarry set and reduced sum
	quarrySet.resize(levelMax);
	reducedSum.resize(levelMax, 0.0);
	
	// to level 0
	h_increaseLevel();
	h_evalLevel(levelCur);
}


lunarEjecta_FractalIntegration::~lunarEjecta_FractalIntegration() {}


double lunarEjecta_FractalIntegration::evalIntegral()
{
	double curError = 10.*epsError;



	while(curError > epsError && levelCur <= levelMax)
	{

	}

}
