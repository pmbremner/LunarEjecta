#include "lunarEjecta_AdaptiveMesh.h"
#include "lunarEjecta_FractalIntegration.h"
#include "lunarEjecta_GeneralExpressions.h"

using namespace std;

#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>

lunarEjecta_AdaptiveMesh::lunarEjecta_AdaptiveMesh
		(int new_Nx,       // number of cell centers
		 int new_Ny,       // number of cell cetners
		 int new_iterMax,
		 int new_levelMax) 
{
	int i,j;
	cout << " Adaptive Mesh Properties: \n";

	// note, we are adding 1 here since Nx and Ny represent edges, not cells
	Nx = new_Nx + 1;
	Ny = new_Ny + 1;

	iterMax = new_iterMax;
	levelMax = new_levelMax;

	cout << "  Nx, Ny = " << Nx << ", " << Ny << endl;
	cout << "  iterMax, levelMax = " << iterMax << ", " << levelMax << endl;

	// init the grid (cell boundary) arrays (x & y) and eval array (cell centers) (z)
	x.resize(Nx);
	y.resize(Ny);
	z.resize(Nx-1);

	for (i = 0; i < Nx-1; ++i)
		z[i].resize(Ny-1);

	// init x (boundaries of x grid)
	// change of variables, linear in zenith angle
	for (i = 0; i < Nx; ++i)
	{
		x[i] = 2. * sqr(sin(PI/4. * double(i)/double(Nx-1.))); // = 1 - cos(alpha)
		cout << i << " | x = " << x[i] << endl;
	}

	// init y (boundaries of y = v/v_esc grid)
	for (i = 0; i < Ny; ++i)
		y[i] = double(i) / double(Ny-1.);

	// init z (cell centered, for integral eval data)
	for (i = 0; i < Nx-1; ++i)
		for (j = 0; j < Ny-1; ++j)
			z[i][j] = 0.;


	evalCount_skipped = 0;
	evalCount_easy = 0;
	evalCount_hard = 0;

	funcCount_skipped = 0;
	funcCount = 0;
}

lunarEjecta_AdaptiveMesh::~lunarEjecta_AdaptiveMesh() {}


void lunarEjecta_AdaptiveMesh::evalBin(double D0, double D1)
{
	int i, j;

	double xMin, xMax, yMin, yMax;

	// break up each bin if needed
	for (i = 0; i < Nx-1; ++i)
	{
		xMin = x[i];
		xMax = x[i+1];

		for (j = 0; j < Ny-1; ++j)
		{
			yMin = y[j];
			yMax = y[j+1];

			z[i][j] = H_r_evalBin(xMin, xMax, yMin, yMax, D0, D1, 0);

			cout << i << ' ' << j;
			cout << "  x-range | y-range | integral = ";
			cout << xMin << ' ' << xMax << " | ";
			cout << yMin << ' ' << yMax << " | ";
			cout << z[i][j] << endl;

		}
	}
}

void lunarEjecta_AdaptiveMesh::printDataToFile(string fn)
{
	ofstream file;
	file.open(fn);

	cout << fn << endl;
	int i, j;


	for (i = 0; i < Nx-1; ++i)
	{
		for (j = 0; j < Ny-1; ++j)
		{
			file << z[i][j] << ' ';
		}
		file << endl;
	}
	file.close();
}


double lunarEjecta_AdaptiveMesh::getReducedIntegral()
{
	int i, j;
	double sum = 0.;

	for (i = 0; i < Nx-1; ++i)
	{
		for (j = 0; j < Ny-1; ++j)
		{
			sum += z[i][j];
		}
	}
	return sum;
}

void lunarEjecta_AdaptiveMesh::printEvalCounts()
{
	cout << "---lunarEjecta_AdaptiveMesh---\n";
	cout << "Fractal calls skipped:\n";
	cout << "---> " << evalCount_skipped << endl;
	cout << "Fractal calls with 0 max depth:\n";
	cout << "---> " << evalCount_easy << endl;
	cout << "Fractal calls with " << levelMax << " max depth:\n";
	cout << "---> " << evalCount_hard << endl;

	cout << "---lunarEjecta_FractalIntegration---\n";
	cout << "Function eval counts skipped:\n";
	cout << "---> " << funcCount_skipped << endl;
	cout << "Function eval counts total:\n";
	cout << "---> " << funcCount << endl;

}



double lunarEjecta_AdaptiveMesh::H_r_evalBin
		(double xMin,
		 double xMax,
		 double yMin,
		 double yMax,
		 double D0,
		 double D1,
		 int iter)
{
	vector<double> D_node(4);
	vector<bool>   maskD0_node(4);
	vector<bool>   maskD1_node(4);

	double binSum = 0.;
	double intEval = 0.;

	iter++;

	// compute distances at each node
	D_node[0] = HH_calcDist(xMin, yMin);
	D_node[1] = HH_calcDist(xMin, yMax);
	D_node[2] = HH_calcDist(xMax, yMin);
	D_node[3] = HH_calcDist(xMax, yMax);
	// D_node[0] = yMin;
	// D_node[1] = yMax;
	// D_node[2] = yMin;
	// D_node[3] = yMax;

	// compute mask for D0 and D1 for each node
	for (int i = 0; i < 4; ++i)
	{
		maskD0_node[i] = D_node[i] < D0 ? 0 : 1;
		maskD1_node[i] = D_node[i] < D1 ? 0 : 1;
	}

	cout << "-------------------------------------------" << endl;
	cout << "---Eval Report -> Iter = " << iter << " | iter max = " << iterMax << endl;
	cout << "AllSame check: " << AllSame(maskD0_node) << ' ' << AllSame(maskD1_node) << ' ' << (AllSame(maskD0_node) && AllSame(maskD1_node)) << endl;
	cout << "AND check: " << AND(maskD0_node) << ' ' << AND(maskD1_node) << ' ' << (AND(maskD0_node) ^ AND(maskD1_node)) << endl;

	cout << "  x-range | y-range ";
	cout << xMin << ' ' << xMax << " | ";
	cout << yMin << ' ' << yMax << endl;

	// if all the nodes are on one side of both of the curves
	if (AllSame(maskD0_node) && AllSame(maskD1_node))
	{
		// if all nodes below D1 and above D0 (completely inside)
		if (AND(maskD0_node) ^ AND(maskD1_node))
		{
			evalCount_easy++;
			lunarEjecta_FractalIntegration scheme(xMin, xMax, yMin, yMax, D0, D1, 0.1, 0);
			intEval = scheme.evalIntegral(); // max steps = 0, so we can eval the domain right away
			scheme.incEvalCounts(funcCount, funcCount_skipped);
			return intEval;
		}
		else // else, completely outside
		{
			evalCount_skipped++;
			return 0.;
		}
	} 
	else // else, the curve crosses inside the domain of the cell
		 //  and we need to further subdivide
	{
		double dx = (xMax - xMin) / 2.;
		double dy = (yMax - yMin) / 2.;

		if (iter < iterMax)
		{
			binSum += H_r_evalBin(xMin, xMin+dx, yMin, yMin+dy, D0, D1, iter); // left bottom quad
			binSum += H_r_evalBin(xMin+dx, xMax, yMin, yMin+dy, D0, D1, iter); // right bottom quad
			binSum += H_r_evalBin(xMin, xMin+dx, yMin+dy, yMax, D0, D1, iter); // left top quad
			binSum += H_r_evalBin(xMin+dx, xMax, yMin+dy, yMax, D0, D1, iter); // right top quad
		}
		else // we still have a crossing, so eval integral with max step size at most 10 (maybe lower)
		{
			evalCount_hard++;
			lunarEjecta_FractalIntegration scheme(xMin, xMax, yMin, yMax, D0, D1, 0.01, levelMax);
			intEval =  scheme.evalIntegral();
			scheme.incEvalCounts(funcCount, funcCount_skipped);
			return intEval;
		}
	}

	return binSum;
}



double lunarEjecta_AdaptiveMesh::HH_calcDist(double x, double y) {
	return 1./PI * atan2(2.*sqr(y) * (1.-x) * sqrt(x*(2.-x)) , 1. - 2.*sqr(y) * x*(2.-x));
}

// AllSame returns AND(m) OR AND(NOT(m))
bool lunarEjecta_AdaptiveMesh::AllSame(vector<bool>& m)
{
	return AND(m) || ANDnot(m);
}

bool lunarEjecta_AdaptiveMesh::AND(vector<bool>& m)
{
	bool m_temp = 1;
	for (int i = 0; i < m.size(); ++i)
	{
		m_temp &= m[i];
	}
	return m_temp;
}

bool lunarEjecta_AdaptiveMesh::ANDnot(vector<bool>& m)
{
	bool m_temp = 1;
	for (int i = 0; i < m.size(); ++i)
	{
		m_temp &= !m[i];
	}
	return m_temp;
}