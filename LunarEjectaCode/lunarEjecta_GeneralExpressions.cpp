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