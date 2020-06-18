#include "lunarEjecta_GeneralExpressions.h"

using namespace std;

#include <vector>
#include <cmath>
#include <iostream>


void copyVector(vector<double>& xFrom, vector<double>& xTo, int Nx) {
	xTo.resize(Nx);
	for (int i = 0; i < Nx; ++i)
		xTo[i] = xFrom[i];
}


void linspace(vector<double>& x, double xmin, double xmax, int Nx) {
	x.resize(Nx);
	for (int i = 0; i < Nx; ++i)
		x[i] = xmin + (xmax - xmin) * double(i) / double(Nx-1.);
}

void logspace(vector<double>& x, double xmin, double xmax, int Nx, int i0, int im) {
	x.resize(im-i0);
	for (int i = i0; i < im; ++i) {
		x[i-i0] = xmin * pow(xmax/xmin, double(i)/double(Nx-1.));
		//cout << x[i-i0] << endl;
	}
}

void linspaceBinEdges(vector<double>& xl, vector<double>& xr, double xmin, double xmax, int Nbin) {
	double dx = (xmax - xmin) / double(Nbin);
	linspace(xl, xmin, xmax-dx, Nbin);
	linspace(xr, xmin+dx, xmax, Nbin);
}

void logspaceBinEdges(vector<double>& xl, vector<double>& xr, double xmin, double xmax, int Nbin) {
	logspace(xl, xmin, xmax, Nbin+1, 0, Nbin);
	logspace(xr, xmin, xmax, Nbin+1, 1, Nbin+1);
}


double romb_int(double a, double b, double eps, double* vv, double (*f)(double*, double)) {

	double pow4m, hn, err = 10.*eps, fsum;
	int n = 1, m, k;
	vector<double> R;

	// base case, R(0,0), trapezoid rule
	hn = (b-a) / 2.;
	R.push_back(hn * (f(vv, a) + f(vv, b)));

	while(n < 5 || err > eps) { // do at least 2^4=16 function evals

		R.push_back(0.); // extend R array
		fsum = 0.;

		// fill in func eval gaps for next level
		for (k = 1; k <= (1 << (n-1)); k++)
			fsum += f(vv, a + (2.*k - 1.)*hn);

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
	return R[tri_idx(n-1, n-1)];
}


double d_r(int n, double a, double b) {
	return (n % 2 == 0 ? n/2.*(b-n/2.)/((a+n-1.)*(a+n)) : -(a+(n-1.)/2.)*(a+b+(n-1.)/2.)/((a+n-1.)*(a+n)) );
}

double Beta(double a, double b) {
	// helps to avoid overflow errors doing it this way
	return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
}

// doesn't catch negative a and b integers, but doesn't matter, we don't have negative vals
// I wrote my own incomplete beta function to avoid having to use weird libraries
// I used the continued fraction representation of the incomplete beta function
// see-> https://dlmf.nist.gov/8.17
//    along with a recursion relation for subsequent terms in the continued fraction (avoids having to recalculate the whole thing over and over)
// see-> https://en.wikipedia.org/wiki/Continued_fraction#Infinite_continued_fractions_and_convergents
// IMO, it's much better than this implementation, and uses a different approach: https://stackoverflow.com/questions/10927859/incomplete-beta-function-in-raw-c
double iBeta(double x, double a, double b)
{
	if (x > (a+1.) / (a+b+2.))
		return Beta(a,b) - iBeta(1.-x, b, a);

	if(x == 0.) // need to check for zero case (even though it's a float, it happens)
		return 0.;

	double f0 = 1., fn = 1.;
	double a0 = 1.;

	double h0 = 0.;
	double h1 = 1., hn;

	double k0 = 1.;
	double k1 = 1., kn;

	int n = 2;
	int nMax = 100;

	// avoids having to check every loop
	nMax = (fmod(b, 1.) < 0.0001 ? int(2*b + 1) : nMax);
	nMax = (fmod(a, 1.) < 0.0001 ? (nMax < a + 2 ? nMax : int(a + 2) ) : nMax);

	do
	{
		f0 = fn;

		a0 = 1. / (a0 * d_r(n-1, a, b) * x);
		
		hn = a0 * h1 + h0;
		kn = a0 * k1 + k0;

		fn = hn / kn;
		n++;

		// swap
		h0 = h1;
		h1 = hn;
		k0 = k1;
		k1 = kn;

	} while (n < nMax && fabs((f0-fn)/fn) > 0.00001);

	return fn * pow(x, a) * pow(1.-x, b) / a;
}