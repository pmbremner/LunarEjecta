#ifndef LUNAREJECTA_ADAPTIVEMESH_H
#define LUNAREJECTA_ADAPTIVEMESH_H


using namespace std;

#include <vector>
#include <cmath>
#include <string>

class lunarEjecta_AdaptiveMesh
{
public:
	// new_x and new_y are bin edges
	lunarEjecta_AdaptiveMesh(int new_Nx,
							 int new_Ny,
							 int new_iterMax,
							 int new_levelMax,
							 double new_Vmin,
		 					 double new_Vmax);
	~lunarEjecta_AdaptiveMesh();

	void evalBins(double D0,
	             double D1,
	             double new_x_azm,
	             double new_Dbeta,
	             double new_mu,
	             double new_imp_zenith);

	void restartBins();

	void printDataToFile(string fn);

	double getReducedIntegral();

	void printEvalCounts();

	// yes, we are making these public. It's bad practice, but I need easy access to z
	vector<double> x;
	vector<double> y;
	vector<vector<double>> z; // z[Nx][Ny]

private:

	double H_r_evalBin(double xMin,
		               double xMax,
		               double yMin,
		               double yMax,
		               double D0,
		               double D1,
		               int iter);

	double HH_calcDist(double x, double y);

	// AllSame returns AND(m) OR AND(NOT(m))
	bool AllSame(vector<bool>& m);
	bool AND(vector<bool>& m);
	bool ANDnot(vector<bool>& m);


	// for integrand, initialized when evalIntegral is called
	double x_azm; // = beta - beta_i, units of rads
	double Dbeta; // units of rads
	double mu;
	double imp_zenith;

	int Nx;
	int Ny;

	int iterMax;
	int levelMax;

	double Vmin;
	double Vmax;
	
	int evalCount_skipped;
	int evalCount_easy;
	int evalCount_hard;

	int funcCount_skipped;
	int funcCount;

};

#endif