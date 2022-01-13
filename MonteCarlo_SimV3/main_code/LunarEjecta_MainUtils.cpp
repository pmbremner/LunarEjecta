#include "LunarEjecta_MainUtils.h"


using namespace std;


double vMax(vector<double>& v){
	double max = v[0];
	for (int i = 1; i < v.size(); ++i)
		if (v[i] > max)
			max = v[i];

	return max;
}

double vMin(vector<double>& v){
	double min = v[0];
	for (int i = 1; i < v.size(); ++i)
		if (v[i] < min)
			min = v[i];

	return min;
}

double vSum(vector<double>& v)
{
	double sum = 0.;
	for (int i = 0; i < v.size(); ++i)
		sum += v[i];
	return sum;
}

void linspace(vector<double>& v, double vmin, double vmax, int Nv)
{
	v.clear();
	for (int i = 0; i < Nv; ++i)
		v.push_back(vmin + (vmax - vmin) * i / double(Nv-1.));
}

void logspace(vector<double>& v, double pmin, double pmax, int Nv)
{
	v.clear();
	for (int i = 0; i < Nv; ++i)
		v.push_back( pow(10., pmin + (pmax - pmin) * i / double(Nv-1.)) );
}

// https://www.learncpp.com/cpp-tutorial/function-pointers/
using RHS_func = double(*)(double, vector<double>&);
//double findX(double LHS, RHS_func RHS, double a, double b, vector<double>& vars);

// Modified False Position, Chapter 4.5 of Numerical Methods for Scientists and Engineers
// Note, if root is complex and goes of to +inf, then findX will return close to b
double findX(double LHS, RHS_func RHS, double a, double b, vector<double>& vars)
{
	double fa = RHS(a, vars) - LHS;
	double fb = RHS(b, vars) - LHS;
	double fx, x;
	//int count = 0;
	
	do{
		x = (a*fb - b*fa) / (fb - fa);
		fx = RHS(x, vars) - LHS;

		//cout << a << ' ' << b << ' ' << x << ' ' << fa << ' ' << fb << ' ' << fx << endl;

		if (fa*fx < 0)
		{
			b = x;
			fb = fx;
			fa /= 2.;
		}
		else if (fa*fx > 0)
		{
			a = x;
			fa = fx;
			fb /= 2.;
		}
		else
		{
			return x;
		}
		//count++;

	} while (fabs(2*(a-b)/(a+b)) > 1.E-5);
	
	//cout << endl;
	//cout << count << ' ';

	return x;
}

// R and D in units of rm
double azm_FOV(double R, double D)
{
	double denom = sqr(sin(D+R)) - sqr(sin(R));

	if (denom >= 0.)
		return 2. * atan2(sin(R), sqrt(denom));
	else
		return PI;
}


double azm_bearing(double lat1, double lon1, double lat2, double lon2, bool short_far_flag)
{
	double phase_shift = (short_far_flag == 0 ? 0. : PI);

	return atan2(sin(lon1 - lon2) * cos(lat1 - lat2) , cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon1 - lon2)) + phase_shift;
}


// return distance in units of rm
double lat_lon_dist(double lat1, double lon1, double lat2, double lon2, bool short_far_flag)
{
	double a = sqr( sin((lat1 - lat2) / 2.) ) + cos(lat1) * cos(lat2) * sqr( sin((lon1 - lon2) / 2.) );
	double pm = (short_far_flag == 0 ? 1. : -1.);

	return 2. * atan2(pm * sqrt(a), sqrt(1. - a));
}


// returns lat2
// d in units of rm
double destination_lat(double lat1, double d, double azm)
{
	return asin(sin(lat1) * cos(d) + cos(lat1) * sin(d) * cos(azm));
}

// returns lon2, needs lat2
// d in units of rm
double destination_lon(double lat1, double lon1, double lat2, double d, double azm)
{
	return lon1 + atan2(sin(azm) * sin(d) * cos(lat1), cos(d) - sin(lat1) * sin(lat2));
}


// final speed at asset
/// a is the asset altitude in units of rm, vp is the ejecta speed at primary impact in units of escape speed
/// returns ejecta speed at asset in units of escape speed
double final_speed(double a, double vp)
{
	return sqrt(1./a + sqr(vp) - 1.);
}



// final zenith at asset (as seen from the asset)
/// d is in units of rm, vp is in units of vesc, g is zenith at ejected point in rads
/// return ejecta zenith at asset in units of radians
double final_zenith(double d, double vp, double g)
{
	return -atan2((sqr(vp) * (cos(d-2.*g) - cos(d)) + cos(d) - 1.), (sqr(vp) * (sin(d-2.*g) - sin(d)) + sin(d)));
}


// the smallest zenith angle to reach asset at ~ escape speed
/// a is the asset altitude in units of rm, d_rm is the projected distance from the impact point to asset in units of rm
/// return ejecta zenith, at ejecta point, in radians
double min_zenith_at_escape(double a, double d_rm)
{
	double x = (1./a - cos(d_rm)) / (1. - cos(d_rm));

	double sqrt_discr = sqrt( 2.*(x-1.) + sqr( 1./cos(d_rm/2.) ) );

	return atan2((1.-x) / tan(d_rm/2.) + fabs(x)*sqrt_discr, x*(x-1.) + fabs(1./tan(d_rm/2.))*sqrt_discr) / 2.;
}