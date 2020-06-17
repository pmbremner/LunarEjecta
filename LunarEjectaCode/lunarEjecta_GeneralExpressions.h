#ifndef LUNAREJECTA_GENERALEXPRESSIONS_H
#define LUNAREJECTA_GENERALEXPRESSIONS_H

using namespace std;

#include <vector>
#include <cmath>

const double PI = 3.14159265358979323846;
const double DtoR = PI/180.;
const double EulerGamma = 0.5772156649015328;
inline double sqr(double x) {return x*x;}

void copyVector(vector<double>& xFrom, vector<double>& xTo, int Nx);

void linspace(vector<double>& x, double xmin, double xmax, int Nx);
void logspace(vector<double>& x, double xmin, double xmax, int Nx, int i0, int im);

void linspaceBinEdges(vector<double>& xl, vector<double>& xr, double xmin, double xmax, int Nbin);
void logspaceBinEdges(vector<double>& xl, vector<double>& xr, double xmin, double xmax, int Nbin);

inline int tri_idx(int n, int m) { return ((n * (n+1)) >> 1) + m; }

double romb_int(double a, double b, double eps, double* vv, double (*f)(double*, double));


double d_r(int n, double a, double b);
double Beta(double a, double b);
double iBeta(double x, double a, double b);



#endif 