#ifndef LUNAREJECTA_FRACTALINTEGRATION_H
#define LUNAREJECTA_FRACTALINTEGRATION_H


using namespace std;

#include <vector>
#include <cmath>
#include <iostream>


struct coord
{
	double x;
	double y;
};

struct set
{
	double dist; // units of lunar circumference
	coord loc;   // (x,y) | x = 1 - cos(zenith), y = vel [v_esc]
	bool isEval;   // 0 = no eval, 1 = yes eval
	double evalVal; // evaluation of integrand for integral
};


class lunarEjecta_FractalIntegration
{
public:
	lunarEjecta_FractalIntegration(double new_xMin,
								   double new_xMax,
								   double new_yMin,
								   double new_yMax,
								   double new_epsError);
	~lunarEjecta_FractalIntegration();


	double evalIntegral();

	// the complete integrand is defined in here, to avoid weird function pointer stuff with classes
	double integrand(double x, // location of center
		             double y, //  of eval
		             double dx, // width of eval in x
		             double dy, // width of eval in y
		             double<vector>& vv); // other params needed (should pack beforehand)


private:

	void h_increaseLevel(); // from levelCur -> levelCur++
	void h_evalLevel(int lev);

	void h_initSet(set* s);
	int hh_getNumQuarryPoints(int lev);

	// normalized units
	double xMin;
	double xMax;
	double yMin;
	double yMax;

	double epsError;

	// level of fractal divisions
	// level  | quarry points
	//   0    |  1
	//   1    |  5
	//   2    |  21
	//   3    |  85
	//   4    |  341
	//   5    |  1365
	//   6    |  5461
	//   7    |  21845
	//   8    |  87381
	//   9    |  349525
	//   10   |  1398101
	int levelMin;
	int levelMax;
	int levelCur;

	vector<vector<set>> quarrySet;
	vector<double> reducedSum; // the sum of all evaluated points
};


#endif