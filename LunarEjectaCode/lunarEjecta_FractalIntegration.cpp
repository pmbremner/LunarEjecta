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
		 double new_D0,
		 double new_D1,
		 double new_epsError,
		 double new_levelMax)
{
	xMin = new_xMin;
	xMax = new_xMax;
	yMin = new_yMin;
	yMax = new_yMax;

	D0 = new_D0;
	D1 = new_D1;

	epsError = new_epsError;

	evalCount_skipped = 0;
	evalCount = 0;

	levelMin = 0;
	levelMax = new_levelMax;
	levelCur = -1;

	// init size of quarry set and reduced sum, including the last level
	quarrySet.resize(levelMax+1);
	reducedSum.resize(levelMax+1, 0.0);
	

	///// Note, we won't init level zero here, but when we call evalIntegral in the adaptive mesh class
	// to level 0
	// h_increaseLevel();
	// h_evalLevel_reduce(levelCur);

	// //evalIntegral();
}


lunarEjecta_FractalIntegration::~lunarEjecta_FractalIntegration() {}

// Note, we can probably apply Richardson extrapolation here,
//  but I'm not sure how the error goes... so maybe later
double lunarEjecta_FractalIntegration::evalIntegral(double new_x_azm,
		                							double new_Dbeta,
		                							double new_mu,
									                double new_imp_zenith,
									                double new_excZone)
{
	double curError = 10.*epsError;

	x_azm = new_x_azm;
	Dbeta = new_Dbeta;
	mu = new_mu;
	imp_zenith = new_imp_zenith;
	excZone = new_excZone;

	// moved from init function
	if(levelCur < 0)
	{
		h_increaseLevel();
		h_evalLevel_reduce(levelCur);
	}

	//while(levelCur <= levelMin || (levelCur < levelMax && curError > epsError))
	while(levelCur < levelMax && (levelCur <= levelMin || curError > epsError))
	{
		// increase level and eval the points
		h_increaseLevel();
		h_evalLevel_reduce(levelCur);

		// compute relative error
		curError = fabs((reducedSum[levelCur] - reducedSum[levelCur-1]) / reducedSum[levelCur]);
	}

	// cout << " - - - evalIntegral: level = " << levelCur << " | " << levelMax << endl;
	// cout << " - - - relative error " << curError << " | " << epsError << endl;

	return reducedSum[levelCur];
}


void lunarEjecta_FractalIntegration::printQuarryPoints()
{
	int ilev, jpoint;
	int Npoints;

	for (ilev = 0; ilev <= levelCur; ++ilev)
	{
		Npoints = hh_getNumQuarryPointsAtLevel(ilev);

		for (jpoint = 0; jpoint < Npoints; ++jpoint)
		{
			hh_printSetPoint(&quarrySet[ilev][jpoint]); 
		}
	}
}

void lunarEjecta_FractalIntegration::printQuarryPointsIfEval() {
	int ilev, jpoint;
	int Npoints;

	for (ilev = 0; ilev <= levelCur; ++ilev)
	{
		Npoints = hh_getNumQuarryPointsAtLevel(ilev);

		for (jpoint = 0; jpoint < Npoints; ++jpoint)
		{
			if(quarrySet[ilev][jpoint].isEval)
				hh_printSetPoint(&quarrySet[ilev][jpoint]); 
		}
	}
}


void lunarEjecta_FractalIntegration::printEvalCounts()
{
	cout << "---lunarEjecta_FractalIntegration---\n";
	cout << "Eval counts skipped:\n";
	cout << "---> " << evalCount_skipped << endl;
	cout << "Eval counts total:\n";
	cout << "---> " << evalCount << endl;
}

void lunarEjecta_FractalIntegration::incEvalCounts(int& c, int& sc)
{
	c += evalCount;
	sc += evalCount_skipped;
}



void lunarEjecta_FractalIntegration::h_increaseLevel()
{

	if(levelCur == levelMax) {
		cout << " Increase level failed... reached max level = " << levelMax << endl;
		return;
	}

	levelCur++;
	//cout << "---level increased to: " << levelCur << endl;

	int Npoints = hh_getNumQuarryPointsAtLevel(levelCur);
	int NpointsPrev;
	double dx = (xMax - xMin) / double(2 << levelCur);
	double dy = (yMax - yMin) / double(2 << levelCur);
	double x, y;
	// cout << "   dx | dy = " << dx << ' ' << dy << endl;
	// cout << "   Point in level = " << Npoints << endl;
	// cout << "   Point total = " << hh_getNumQuarryPointsTotal(levelCur) << endl;

//////////////////
	//fractIntCount += Npoints;
//////////////////

	quarrySet[levelCur].resize(Npoints);

	if (levelCur == 0) // base case
	{
		x = xMin + dx;
		y = yMin + dy;

		hh_initSet(&quarrySet[0][0], x, y);
	}
	else // all other cases > level 0
	{
		NpointsPrev = hh_getNumQuarryPointsAtLevel(levelCur-1);

		// loop through locations of all previous points
		//  adding 4 new points arround each
		for (int i = 0; i < NpointsPrev; ++i)
		{
			// get position of previous point
			x = quarrySet[levelCur-1][i].loc.x;
			y = quarrySet[levelCur-1][i].loc.y;

			// add 4 new points around previous point
			hh_initSet(&quarrySet[levelCur][4*i    ], x - dx, y - dy);
			hh_initSet(&quarrySet[levelCur][4*i + 1], x + dx, y - dy);
			hh_initSet(&quarrySet[levelCur][4*i + 2], x - dx, y + dy);
			hh_initSet(&quarrySet[levelCur][4*i + 3], x + dx, y + dy);
		}
	}
}


void lunarEjecta_FractalIntegration::h_evalLevel_reduce(int lev)
{
	int i;
	int Npoints = hh_getNumQuarryPointsAtLevel(levelCur);
	int NpTotal = hh_getNumQuarryPointsTotal(levelCur);

	double x, y;
	double dx = (xMax - xMin) / sqrt(double(NpTotal));
	double dy = (yMax - yMin) / sqrt(double(NpTotal));

	// init sum to zero
	reducedSum[levelCur] = 0.0;

	//cout << " Npoints = " << Npoints << endl;

	// loop through points at current level
	for (i = 0; i < Npoints; ++i)
	{
		x = quarrySet[levelCur][i].loc.x;
		y = quarrySet[levelCur][i].loc.y;

		// Check if we need to evaluate function
		// this in temporary
		//if(y > sqrt(1.-sqr(x)) && y < cos(x)) {
		//if( y < cos(x)) {
		// if(!(quarrySet[levelCur][i].dist >= D0 && quarrySet[levelCur][i].dist <= D1)){
		// 	cout << " D0, D1 " << D0 << ' ' << D1 << ' ' << (quarrySet[levelCur][i].dist >= D0 && quarrySet[levelCur][i].dist <= D1) << endl;
		// 	cout << quarrySet[levelCur][i].dist << endl;
		// }
		// if(1)
		if (quarrySet[levelCur][i].dist >= D0 && quarrySet[levelCur][i].dist <= D1)
		{
			evalCount++;
			quarrySet[levelCur][i].evalVal = integrand(x, y, dx, dy);
			quarrySet[levelCur][i].isEval = 1;
			reducedSum[levelCur] += quarrySet[levelCur][i].evalVal;
		}
		else
		{
			evalCount_skipped++;
			quarrySet[levelCur][i].evalVal = 0.;
			quarrySet[levelCur][i].isEval = 0;
		}
	}

	// renormalize and reduce previous levels
		reducedSum[levelCur] += h_renormReduce(levelCur-1, dx, dy);
	
}

double lunarEjecta_FractalIntegration::h_renormReduce(int lev, double dx, double dy)
{
	int NpTotal = hh_getNumQuarryPointsTotal(lev);
	double dxlev = (xMax - xMin) / sqrt(double(NpTotal));
	double dylev = (yMax - yMin) / sqrt(double(NpTotal));

	return reducedSum[lev] * ((dx * dy) / (dxlev * dylev));
}


void lunarEjecta_FractalIntegration::hh_initSet(set* s, double x, double y)
{
	s->loc.x = x;
	s->loc.y = y;

	s->isEval = 0;
	s->evalVal = 0.0;

	s->dist = HH_calcDist(x, y);

	//s->dist = y; // this is just for a test
	//hh_printSetPoint(s);
	//HH_calcDist(x, v);
}

// evaluation of the integrand in the domain Dx and Dy centered on x and y
double lunarEjecta_FractalIntegration::integrand
                (double x_zenith, // location of center, x_zenith = 1 - cos(zenith)
		         double y, //  of eval, y = v/v_esc
		         double Dx, // width of eval in x_zenith
		         double Dy) // width of eval in y
{
//////////	
	// double y0 = y - Dy/2.;
	// double y1 = y + Dy/2.;

	// double x0 = x_zenith - Dx/2.;
	// double x1 = x_zenith + Dx/2.;
	// //double zenith = acos(1. - x_zenith); // units of rads

	//double a = a_power(imp_zenith, x_azm); // needs impact zenith angle, not secondary zenith angle
	//double a = a_power(0.0, 0.0); // Temp, testing a that doesn't depend on azm

	// return Dbeta * Dx
	// 	* (pow(y0, -3.*mu) - pow(y1, -3.*mu))        // exact integral of speed dist
	// 	* pow(x_zenith, 1./a) * pow(1.-x_zenith, a)  // zenith term, first-order Taylor series approx
	// 	* HH_AzmDist(imp_zenith, x_azm);				 // azimuth term, first-order Taylor series approx
////////////

	///* Case 2  */ return Dbeta * (pow(y0, -3.*mu) - pow(y1, -3.*mu)) * (iBeta(x1, 1.+ 1./a45, 1.+ a45) -iBeta(x0, 1.+ 1./a45, 1.+ a45) );
	
	///*run_equator_A3.txt */ return Dx * Dbeta * (pow(y0, -3.*mu) - pow(y1, -3.*mu)) * pow(x_zenith, 1./a) * pow(1.-x_zenith, a);
	///*Case 1 run_equator_A2.txt */ return Dx * Dbeta * (pow(y0, -3.*mu) - pow(y1, -3.*mu));
	///*run_equator_A1.txt*/ return Dx * Dbeta * Dy;// * (pow(y0, -3.*mu) - pow(y1, -3.*mu));
	// /* Case 0 */ return Dx * Dy;


	return H_betaIntegrate(x_zenith, y, Dx, Dy);

}

// zenith in units of rad, x = \beta - \beta_i
double lunarEjecta_FractalIntegration::HH_AzmDist(double x)
{

	//return 1.;
	if(imp_zenith < PI/3.)
	{
		return (1. + cos(x) * 3.*imp_zenith / (2.*PI - 3.*imp_zenith)) / (2.*PI);

	}
	else
	{ // imp_zenith > PI/3

		return (1. + cos(x) ) / (2.*PI);
	}
}

double lunarEjecta_FractalIntegration::HH_integrand(double x0, double x1, double y0, double y1, double x_beta)
{
	if(x_beta > PI - excZone && x_beta < PI + excZone)
	{ // in the exclusion zone
		return 0.0;
	}
	else
	{
		double a = a_power(imp_zenith, x_beta);
		return HH_AzmDist(x_beta) * (iBeta(x1, 1.+ 1./a, 1.+ a) -iBeta(x0, 1.+ 1./a, 1.+ a) );
	}
}


double lunarEjecta_FractalIntegration::H_betaIntegrate
				(double x_zenith, // location of center, x_zenith = 1 - cos(zenith)
		         double y, //  of eval, y = v/v_esc
		         double Dx, // width of eval in x_zenith
		         double Dy) // width of eval in y
{
	double y0 = y - Dy/2.;
	double y1 = y + Dy/2.;

	double x0 = x_zenith - Dx/2.;
	double x1 = x_zenith + Dx/2.;

	double a = x_azm - Dbeta/2.;
	double b = x_azm + Dbeta/2.;

	double eps = 1.E-2; // epsilon error in numerical integration

	double pow4m, hn, err = 10.*eps, fsum;
	int n = 1, m, k;
	vector<double> R;

	// base case, R(0,0), trapezoid rule
	hn = (b-a) / 2.;
	R.push_back(hn * (HH_integrand(x0, x1, y0, y1, a) + HH_integrand(x0, x1, y0, y1, b)));

	while(n < 4 || err > eps) { // do at least 2^3=8 function evals

		R.push_back(0.); // extend R array
		fsum = 0.;

		// fill in func eval gaps for next level
		for (k = 1; k <= (1 << (n-1)); k++)
			fsum += HH_integrand(x0, x1, y0, y1, a + (2.*k - 1.)*hn);

		R[tri_idx(n,0)] = R[tri_idx(n-1,0)] / 2. + hn * fsum;

		pow4m = 1.;

		// compute n-th row of m's
		for (m = 1; m <= n; m++)
		{
			R.push_back(0.); // extend R array
			pow4m *= 4.;
			R[tri_idx(n,m)] = R[tri_idx(n,m-1)]
			                + (R[tri_idx(n,m-1)] - R[tri_idx(n-1,m-1)])
			                / (pow4m - 1.);

		}

		// compute relative error
		err = fabs((R[tri_idx(n,n)] - R[tri_idx(n-1,n-1)]) / R[tri_idx(n,n)]);
		n++;
		hn /= 2.;
	}
	return R[tri_idx(n-1, n-1)] * (pow(y0, -3.*mu) - pow(y1, -3.*mu));
}


int lunarEjecta_FractalIntegration::hh_getNumQuarryPointsTotal(int lev)
{
	return (double(4 << (2*lev)) - 1.) / 3.;
}

int lunarEjecta_FractalIntegration::hh_getNumQuarryPointsAtLevel(int lev)
{
	return 1 << 2*lev;
}


void lunarEjecta_FractalIntegration::hh_printSetPoint(set* s)
{
	cout << s->loc.x << " " <<  s->loc.y << endl;
}

double lunarEjecta_FractalIntegration::HH_calcDist(double x, double v) {
	return 1./PI * atan2(2.*sqr(v) * (1.-x) * sqrt(x*(2.-x)) , 1. - 2.*sqr(v) * x*(2.-x));
}


//// Same as in lunarEjecta_Assembly
// zenith in units of rad, x = \beta - \beta_i
// double lunarEjecta_FractalIntegration::HH_AzmDist(double impact_zenith, double x) {
// 	if(impact_zenith < PI/3.)
// 	{
// 		return (1. + cos(x) * 3.*impact_zenith / (2.*PI - 3.*impact_zenith)) / (2.*PI);

// 	} else { // impact_zenith > PI/3

// 		double b = (0.05 - 1.) /(PI/2. - PI/3.) * (impact_zenith - PI/3.) + 1.;

// 		return exp(-(x - 2.*(impact_zenith - PI/3.) ) / (PI*b)) + exp(-(x + 2.*(impact_zenith - PI/3.) ) / (PI*b)); // not normalized to itself, it will get normalized with everything later
// 	}
// }

// zenith in units of rad, x = \beta - \beta_i
// x_angle = 1 - cos(zenith)
// double lunarEjecta_FractalIntegration::HH_AltDist(double zenith, double x)
// {
// 	double x_angle = 1. - cos(zenith);
// 	double a = a_power(zenith, x);
	
// 	return pow(x_angle, 1./a) * pow(1. - x_angle, a);
// }

// double lunarEjecta_FractalIntegration::HH_vDist(double v, double mu)
// {
// 	return 3.*mu * pow(v, -(3.*mu+1.));
// }


// beta - beta_i = pi
inline double lunarEjecta_FractalIntegration::HH_zenithDownstream(double impactZenith) { // input and output in degrees
	return 20. + impactZenith * (1.5206 + impactZenith * (-0.036 + impactZenith * 0.0003));
}

// beta - beta_i = 0
inline double lunarEjecta_FractalIntegration::HH_zenithUpstream(double impactZenith) { // input and output in degrees
	return 20. + impactZenith * (0.129 + impactZenith * (0.0236 + impactZenith * -0.00042));
}

// x = beta - beta_i
// Note: beta_i is the impact azimuth (direction seen from impact point), which is defined to be due East (direction or Moon's rotation)
//       
//   beta is the secondary ejecta direction leaving the impact point (not direction seen from secondary impact point, ie ROI)
//
// impactZenith in degrees, x in radians, output in radians
inline double lunarEjecta_FractalIntegration::HH_zenithGeneral(double impactZenith, double x) {
	return (HH_zenithDownstream(impactZenith)
		+ (HH_zenithUpstream(impactZenith) - HH_zenithDownstream(impactZenith)) / PI * fabs(x)) * DtoR;
}

// a = cos(a_max) / (1 - cos(amax)), converted to half sin to avoid subtraction of possible close #'s
// x = \beta - \beta_i
double lunarEjecta_FractalIntegration::a_power(double impactZenith, double x) {
	double alpha_max = HH_zenithGeneral(impactZenith, x); // in units of radians
	return sqrt(cos(alpha_max)/2.) / fabs(sin(alpha_max / 2.));
	//return cos(alpha_max) / (2. * sqr(sin(alpha_max / 2.))); // forgot the sqrt...
}
