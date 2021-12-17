#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
//#include <float.h> // DBL_MAX


using namespace std;

const double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

// inline functions must be defined in the header
inline double sqr(double x) {return x*x;}
inline double min(double a, double b) {return (a < b ? a : b);}
inline double max(double a, double b) {return (a > b ? a : b);}
inline double mag_s(double x, double y, double z) {return sqrt(sqr(x) + sqr(y) + sqr(z));}
inline double mag2(double x, double y, double z) {return sqr(x) + sqr(y) + sqr(z);}

void linspace(vector<double>& v, double vmin, double vmax, int Nv)
{
	v.clear();
	for (int i = 0; i < Nv; ++i)
		v.push_back(vmin + (vmax - vmin) * i / double(Nv-1.));
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
	int count = 0;
	
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
		count++;

	} while (fabs(2*(a-b)/(a+b)) > 1.E-5);
	
	//cout << endl;
	//cout << count << ' ';

	return x;
}


double Fspeed(double g, double h_rm, double d_rm)
{
	double x = 1. - (h_rm / (1. + h_rm)) / (1. - cos(d_rm)); // for h = 0, x = 1

	return 1. / (x * (1. - cos(2.*g)) + sin(2.*g)/tan(d_rm/2.));
}


double Fspeed_v(double v, vector<double>& vars)
{
	return Fspeed(vars[0], vars[1], vars[2]) - sqr(v);

}

double Fspeed_g(double g, vector<double>& vars)
{
	return Fspeed(g, vars[1], vars[2]) - sqr(vars[0]);
}



// if |f(g+dg) - f(g)| > dv, then find dg such that |f(g+dg) - f(g)| = dv, else keep dg
// i.e., sample the zenith angle finer for when the function is steeper
// Note: make sure vars[1] = h and vars[2] = d before this function call, both in units of lunar radii
double find_dg(double g, double dg, double dv, double vmax, vector<double>& vars)
{
	const double alpha = 0.3;
	const double eps   = 1.E-3;

	double f0, f1;
	double dg_1 = dg;

	do
	{
		// compute the function at g and g+dg
		vars[0] = g;
		f0 = findX(0., Fspeed_v, 0., vmax, vars);

		vars[0] = g + dg_1;
		f1 = findX(0., Fspeed_v, 0., vmax, vars);

		// compute a weighted average dg so as not to overshoot intended result
		// (due to taking a linear approximation of the derivative)
		// alpha ranges from 0 to 1, aim for alpha ~ 0.3, possibly higher
		dg_1 *= dv / fabs(f1 - f0) * alpha + (1. - alpha);

		// if the proposed dg is too large, keep the original dg
		dg_1 = (dg_1 < dg ? dg_1 : dg);

		vars[0] = g + dg_1;
		f1 = findX(0., Fspeed_v, 0., vmax, vars);

	// do while the proposed dv is too large
	} while (fabs(f1 - f0) > dv*(1. + eps));

	return dg_1;
}

// For given g, h, d, r, compute the dg which is the minimum dg of the four corners of the wedge
double find_dg_wedge(double g, double dg, double dv, double vlow, double vmax, vector<double>& vars, double h, double d, double r)
{
	int i, j;
	double dg_min = dg;

	for (i = 0; i < 2; ++i)
		for (j = 0; j < 2; ++j)
		{
			vars[0] = g;
			vars[1] = h * j;
			vars[2] = d + r*(2.*i - 1);

			if(findX(0., Fspeed_v, vlow, vmax, vars) < vmax*0.9999)
				dg_min = min(find_dg(g, dg, dv, vmax, vars), dg_min);
		}

	return dg_min;
}

void find_min_max_v(double g, double& vmin, double& vmax, double vlow, double vlim, vector<double>& vars, double h, double d, double r)
{
	int i, j;
	double f;

	vmin = vlim;
	vmax = 0.;

	vars[0] = g;

	for (i = 0; i < 2; ++i)
		for (j = 0; j < 2; ++j)
		{
			vars[1] = h * j;
			vars[2] = d + r*(2.*i - 1);

			f = findX(0., Fspeed_v, vlow, vlim, vars);

			vmin = min(f, vmin);
			vmax = max(f, vmax);
		}
	// force vmin to be at least vlow, and vmax to be at most vlim
	vmin = max(vlow, vmin);
	vmax = max(min(vlim, vmax), vlow);
}


// zenith in units of rads in (0,Pi/2]
// vmin and vmax in units of vesc in (0, vlim]
// h is the height of the wedge
// d is the distance of the wedge, center
// r is the radius of the wedge (front and back)
// the wedge enscribes a sphere-like shape (h can be controlled, so really it's an ellipsoid)
// dg and dv are the maximum grid spacing for the zenith and speed dimensions, respectively
void get_zenith_speed_grid(vector<double>& zenith, vector<double>& vmin, vector<double>& vmax, double vlow, double vlim, double h, double d, double r, double dg, double dv)
{
	zenith.clear();
	vmin.clear();
	vmax.clear();

	vector<double> vars(3, 0.); // size of 3, filled with zeros

	// First, find the smallest zenith angle at the closest point
	vars[0] = vlim;
	vars[1] = 0.;
	vars[2] = d-r;

	double g_min = 1.01*(d-r)/4.;//findX(0., Fspeed_g, 0.000001, PI/2., vars);

	// compute the grid
	double g_cur, v0, v1, dg_new;

	// the first point
	g_cur = g_min;
	find_min_max_v(g_cur, v0, v1, vlow, vlim, vars, h, d, r);
	cout << g_cur << ' ' << v0 << ' ' << v1 << endl;

	// the rest of the points
	while (g_cur < PI/2.)
	{
		// First, find the dg amount, limited by each corner
		vars[1] = h;
		vars[2] = d;
		dg_new = find_dg_wedge(g_cur, dg, dv, vlow, vlim, vars, h, d, r);

		g_cur += dg_new;
		g_cur = (g_cur > PI/2. ? PI/2. : g_cur);

		// Next, find the minimum and maximum speeds of the corners
		find_min_max_v(g_cur, v0, v1, vlow, vlim, vars, h, d, r);

		cout << g_cur << ' ' << v0 << ' ' << v1 << endl;

		zenith.push_back(g_cur);
		vmin.push_back(v0);
		vmax.push_back(v1);
	}

}




int main(int argc, char const *argv[])
{
	
	vector<double> zenith, vminv, vmaxv;

	double vlow = 0.0;
	double v = 5.;
	double h = 0.22;
	double d = 0.56;//0.287;
	double r = 0.21;

	double dg = 0.1;
	double dv = 0.05;
	

	get_zenith_speed_grid(zenith, vminv, vmaxv, vlow, v, h, d, r, dg, dv);
	//get_zenith_speed_grid(zenith, vminv, vmaxv, vlow, v, h, 2.*PI-d, r, dg, dv);


	return 0;
}