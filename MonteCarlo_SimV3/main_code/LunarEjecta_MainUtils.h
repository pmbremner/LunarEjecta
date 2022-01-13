#ifndef LUNAREJECTA_MAINUTILS_H
#define LUNAREJECTA_MAINUTILS_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <chrono>
#include <random>     // mt19937
#include <algorithm>  // lower_bound
//#include <float.h> // DBL_MAX

// note, -march=native is to allow for vectorization, if possible

//  g++ -O2 -std=c++17 -march=native .\domain_finder_test.cpp

using namespace std;

const double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

// inline functions must be defined in the header
inline double sqr(double x) {return x*x;}
inline double min(double a, double b) {return (a < b ? a : b);}
inline double max(double a, double b) {return (a > b ? a : b);}
inline double mag_s(double x, double y, double z) {return sqrt(sqr(x) + sqr(y) + sqr(z));}
inline double mag2(double x, double y, double z) {return sqr(x) + sqr(y) + sqr(z);}


inline double trapezoid_area(double width, double height0, double height1)
{
	return 0.5 * width * (height0 + height1);
}

inline double uniform(mt19937& rng, double a, double b){
	return a + (b - a) * rng() / double(rng.max());
}

// returns an int from a to b-1
inline int uniformInt(mt19937& rng, int a, int b){
	return int(floor(a + (b - a) * rng() / double(rng.max()+1.)));
}

double vMax(vector<double>& v);
double vMin(vector<double>& v);
double vSum(vector<double>& v);


template <class paramType>
void vCumLow(vector<paramType>& v, vector<paramType>& vCum)
{
	vCum[0] = 0;
	for (int i = 1; i < vCum.size(); ++i)
		vCum[i] = vCum[i-1] + v[i-1];
}

// check this**
// template <class paramType>
// void vCumUp(vector<paramType>& v, vector<paramType>& vCum)
// {
// 	vCum[0] = v[0];
// 	for (int i = 1; i < v.size(); ++i)
// 		vCum[i] = vCum[i-1] + v[i];
// }


void linspace(vector<double>& v, double vmin, double vmax, int Nv);

void logspace(vector<double>& v, double pmin, double pmax, int Nv);

// https://www.learncpp.com/cpp-tutorial/function-pointers/
using RHS_func = double(*)(double, vector<double>&);
//double findX(double LHS, RHS_func RHS, double a, double b, vector<double>& vars);

// Modified False Position, Chapter 4.5 of Numerical Methods for Scientists and Engineers
// Note, if root is complex and goes of to +inf, then findX will return close to b
double findX(double LHS, RHS_func RHS, double a, double b, vector<double>& vars);


double azm_FOV(double R, double D);

// returns azm1, either forward or reverse, depending on short_far_flag
double azm_bearing(double lat1, double lon1, double lat2, double lon2, bool short_far_flag);

// returns distance, d
double lat_lon_dist(double lat1, double lon1, double lat2, double lon2, bool short_far_flag);

// returns lat2
double destination_lat(double lat1, double d, double azm);

// returns lon2, needs lat2
double destination_lon(double lat1, double lon1, double lat2, double d, double azm);


// final speed at asset
/// a is the asset altitude in units of rm, vp is the ejecta speed at primary impact in units of escape speed
/// returns ejecta speed at asset in units of escape speed
double final_speed(double a, double vp);

// final zenith at asset (as seen from the asset)
/// d is in units of rm, vp is in units of vesc, g is zenith at ejected point in rads
/// return ejecta zenith at asset in units of radians
double final_zenith(double d, double vp, double g);


// the smallest zenith angle to reach asset at ~ escape speed
/// a is the asset altitude in units of rm, d_rm is the projected distance from the impact point to asset in units of rm
/// return ejecta zenith, at ejecta point, in radians
double min_zenith_at_escape(double a, double d_rm);



#endif 