#ifndef LUNAREJECTA_GENERALEXPRESSIONS_H
#define LUNAREJECTA_GENERALEXPRESSIONS_H

using namespace std;

#include <vector>
#include <cmath>

const double PI = 3.14159265358979323846;
const double DtoR = PI/180.;
inline double sqr(double x) {return x*x;}

void copyVector(vector<double>& xFrom, vector<double>& xTo, int Nx);

void linspace(vector<double>& x, double xmin, double xmax, int Nx);
void logspace(vector<double>& x, double xmin, double xmax, int Nx, int i0, int im);

void linspaceBinEdges(vector<double>& xl, vector<double>& xr, double xmin, double xmax, int Nbin);
void logspaceBinEdges(vector<double>& xl, vector<double>& xr, double xmin, double xmax, int Nbin);
#endif 