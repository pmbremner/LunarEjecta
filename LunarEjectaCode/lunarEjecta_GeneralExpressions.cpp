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