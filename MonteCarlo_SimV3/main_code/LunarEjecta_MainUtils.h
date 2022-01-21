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

template <class rng_type>
inline double uniform(rng_type& rng, double a, double b){
	return a + (b - a) * rng() / double(rng.max());
}

// returns an int from a to b
template <class rng_type, class int_type>
inline int_type uniformInt(rng_type& rng, int_type a, int_type b){
	return a + rng() % (b - a);
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
void rlogspace(vector<double>& v, double pmin, double pmax, int Nv);

// https://www.learncpp.com/cpp-tutorial/function-pointers/
using RHS_func = double(*)(double, vector<double>&);
//double findX(double LHS, RHS_func RHS, double a, double b, vector<double>& vars);

// Modified False Position, Chapter 4.5 of Numerical Methods for Scientists and Engineers
// Note, if root is complex and goes of to +inf, then findX will return close to b
double findX(double LHS, RHS_func RHS, double a, double b, vector<double>& vars);


// use the cdf to convert a uniform random number in [0,1] to match the corresponding pdf
int sample_pdf_idx(mt19937& rng, vector<double>& cdf);
int sample_pdf_idx(mt19937& rng, vector<double>& cdf, double& u);

#endif 