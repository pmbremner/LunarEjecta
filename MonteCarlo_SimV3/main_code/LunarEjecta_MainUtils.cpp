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
	for (int i = 0; i < v.size(); ++i){
		sum += v[i];
	}
	return sum;
}

void linspace(vector<double>& v, double vmin, double vmax, int Nv)
{
	v.clear();
	v.resize(Nv);
	//cout << endl <<  vmin << ' ' << vmax << endl;
	for (int i = 0; i < Nv; ++i){
		v[i] = vmin + (vmax - vmin) * i / double(Nv-1.);
		//cout << v[i] << endl;
	}
}

void logspace(vector<double>& v, double pmin, double pmax, int Nv)
{
	v.clear();
	v.resize(Nv);
	//cout << endl <<  pmin << ' ' << pmax << endl;
	for (int i = 0; i < Nv; ++i){
		v[i] = pow(10., pmin + (pmax - pmin) * i / double(Nv-1.)) ;
		//cout << v[i] << endl;
	}
}


void rlogspace(vector<double>& v, double pmin, double pmax, int Nv)
{
	v.clear();
	v.resize(Nv);
	for (int i = Nv-1; i >= 0; i--){
		v[i] =  pow(10., pmin + (pmax - pmin) * i / double(Nv-1.)) ;
	}
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


// use the cdf to convert a uniform random number in [0,1] to match the corresponding pdf
int sample_pdf_idx(mt19937& rng, vector<double>& cdf)
{
	// pull sample from uniform distribution
	double u = uniform(rng, 0., 1.);

	// find index (iterator in this case) of the corresponding location in the cdf

	// Find the index such that cdf(idx-1) <= u <= cdf(idx)
	// If u = 0, use upper_bound, otherwise use lower_bound (both binary search algorithms, logN time)
	//  effectively, this inverts the cdf
	// This guarantees that *(idx_iter-1) <= u <= *(idx_iter) for all values of u in [0,1]
	return (u == 0. ? upper_bound(cdf.begin(), cdf.end(), u) : lower_bound(cdf.begin(), cdf.end(), u)) - cdf.begin();
}

// copy of above, but with returning u as well
int sample_pdf_idx(mt19937& rng, vector<double>& cdf, double& u)
{
	// pull sample from uniform distribution
	u = uniform(rng, 0., 1.);

	// find index (iterator in this case) of the corresponding location in the cdf

	// Find the index such that cdf(idx-1) <= u <= cdf(idx)
	// If u = 0, use upper_bound, otherwise use lower_bound (both binary search algorithms, logN time)
	//  effectively, this inverts the cdf
	// This guarantees that *(idx_iter-1) <= u <= *(idx_iter) for all values of u in [0,1]
	return (u == 0. ? upper_bound(cdf.begin(), cdf.end(), u) : lower_bound(cdf.begin(), cdf.end(), u)) - cdf.begin();
}