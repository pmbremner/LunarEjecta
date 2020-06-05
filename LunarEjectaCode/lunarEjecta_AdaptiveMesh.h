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
	lunarEjecta_AdaptiveMesh(vector<double>& new_x,
		                     vector<double>& new_y,
		                     vector<vector<double>>& new_z, // z[Nx][Ny]
							 int new_iterMax,
							 int new_levelMax);
	~lunarEjecta_AdaptiveMesh();

	void evalBin(double D0,
	             double D1);

	void printDataToFile(string fn);

	double getReducedIntegral();

	void printEvalCounts();

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

	vector<double>* x;
	vector<double>* y;
	vector<vector<double>>* z;

	int iterMax;
	int levelMax;
	
	int evalCount_skipped;
	int evalCount_easy;
	int evalCount_hard;

	int funcCount_skipped;
	int funcCount;

};

#endif