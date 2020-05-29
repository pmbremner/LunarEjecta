#include "lunarEjecta_SpeedZenithIntegration.h"
#include "lunarEjecta_GeneralExpressions.h"

using namespace std;

#include <vector>
#include <cmath>
#include <iostream>

lunarEjecta_SpeedZenithIntegration::lunarEjecta_SpeedZenithIntegration
		  (double new_D0,
           double new_D1,
           double new_radius,
           double new_Vesc,
           int new_Nx, // this is a goal, won't be exact due to rounding
           int new_Nv)
{
	circumference = 2. * PI * new_radius;
	escapeSpeed   = new_Vesc;

	D0 = new_D0 / circumference;
	D1 = new_D1 / circumference;

	//// compute region case
	if (D0 < 0.5 && D1 < 0.5)
	{ //  Case 1: D0 and D1 < 0.5, all three regions valid
		regionCase = region_I;
	}
	else if (D0 < 0.5 && D1 >= 0.5) 
	{ //  Case 2: D0 < 0.5 and D1 >= 0.5, regions I and II valid
		regionCase = region_II;
	}
	else if (D0 >= 0.5 && D1 >= 0.5)
	{ //  Case 3: D0 and D1 >= 0.5, only region I valid
		regionCase = region_III;
	}
	else
	{
		cout << "ERROR: lunarEjecta_SpeedZenithIntegration, invalid regionCase\n";
		regionCase = -1;
	}

	// domain of x = 1 - cos(zenith)
	xMin = 1. - cos(D0 * PI/2.);
	xMax = 1.;

	// range of v
	if (regionCase == region_I)
	{
		vMinD0 = this->H_calcVMin(D0);
		vMinD1 = this->H_calcVMin(D1);

		x_at_vMinD0 = H_calcX_at_vMin(D0);
		x_at_vMinD1 = H_calcX_at_vMin(D1);

		Nx_R1 = ceil(new_Nx * x_at_vMinD0);
		Nx_R2 = 2;
		Nx_R3 = ceil( new_Nx * (1.-x_at_vMinD0));
	}
	else if (regionCase == region_II)
	{
		vMinD0 = this->H_calcVMin(D0);
		vMinD1 = sqrt(2.)/2.;

		x_at_vMinD0 = H_calcX_at_vMin(D0);
		x_at_vMinD1 = 1.;

		Nx_R1 = ceil(new_Nx * x_at_vMinD0);
		Nx_R2 = 2;
		Nx_R3 = 0;
	}
	else if (regionCase == region_III)
	{
		vMinD0 = sqrt(2.)/2.;
		vMinD1 = sqrt(2.)/2.;

		x_at_vMinD0 = 1.;
		x_at_vMinD1 = 1.;

		Nx_R1 = ceil(new_Nx * x_at_vMinD0);
		Nx_R2 = 0;
		Nx_R3 = 0;
	}

	vMin = vMinD0;
	vMax = 1.;

	Nv = ceil(new_Nv * (1. - vMin));
	
	// init grid (nodes and edges)
	if (regionCase == region_I)
	{
		cout << " Init grid in Region 1: \n";
		cout << "   D0 = " << D0 << endl; 
		cout << "   D1 = " << D1 << endl; 
		initGrid(&region1Grid, Nx_R1, Nv, xMin, x_at_vMinD0, vMinD0, vMax);

		cout << " Init grid in Region 2: \n";
		cout << "   D0 = " << D0 << endl; 
		cout << "   D1 = " << D1 << endl; 
		initGrid(&region2Grid, Nx_R2, Nv, x_at_vMinD0, x_at_vMinD1, vMinD0, vMax);

		cout << " Init grid in Region 3: \n";
		cout << "   D0 = " << D0 << endl; 
		cout << "   D1 = " << D1 << endl; 
		initGrid(&region3Grid, Nx_R3, Nv, x_at_vMinD1, xMax, vMinD1, vMax);
	}
	else if (regionCase == region_II)
	{

	}
	else if (regionCase == region_III)
	{

	}

 	// init cells for each region
	region1Cells.resize((Nx_R1 - 1) * (Nv - 1));
	region2Cells.resize((Nx_R2 - 1) * (Nv - 1));
	region3Cells.resize((Nx_R3 - 1) * (Nv - 1));


}

lunarEjecta_SpeedZenithIntegration::~lunarEjecta_SpeedZenithIntegration() {}



inline int lunarEjecta_SpeedZenithIntegration::cIdx(int ix, int jv) {
	return jv + (Nv-1) * ix;
}

inline int lunarEjecta_SpeedZenithIntegration::gIdx(int ix, int jv) {
	return jv + Nv * ix;
}

double lunarEjecta_SpeedZenithIntegration::H_calcVMin(double D){
	return pow( 1. + fabs(cos(D * PI)) / tan(D * PI) + sin(D * PI) , -0.5);
}

double lunarEjecta_SpeedZenithIntegration::H_calcX_at_vMin(double D) {
	return 1. - sqrt((1. - sin(D * PI)) / 2.);
}

double lunarEjecta_SpeedZenithIntegration::H_calcDist(double x, double v) {
	return 1./PI * atan2(2.*sqr(v) * (1.-x) * sqrt(x*(2.-x)) , 1. - 2.*sqr(v) * x*(2.-x));
}

void lunarEjecta_SpeedZenithIntegration::initGrid
		 (grid* g,
	      int Nx, // number of nodes in x
	      int Nv, // number of nodes in y
	      double xMin,
	      double xMax,
	      double vMin,
	      double vMax) {

	int i, j;
	double x, v;

	int NVE = Nx * (Nv - 1); // Nv like cell
	int NHE = (Nx - 1) * Nv; // Nv like grid
	int NN  = Nx * Nv;       // Nv like grid

	g->node.resize(NN);
	g->vEdge.resize(NVE);
	g->hEdge.resize(NHE);

	// loop through grid nodes
	cout << "...init nodes...\n";
	for (i = 0; i < Nx; ++i)
	{
		x = xMin + (xMax - xMin) * i / double(Nx - 1);

		for (j = 0; j < Nv; ++j)
		{
			v = vMin + (vMax - vMin) * j / double(Nv - 1);

			this->initSet(&(g->node[this->gIdx(i, j)]), x, v);
		}
	}

}


void lunarEjecta_SpeedZenithIntegration::initSet(set* s, double x, double v) {
	s->loc.x = x;
	s->loc.y = v;

	s->dist = this->H_calcDist(x, v);

	cout << " x, v, D | " << x << ' ' << v << ' ' << s->dist << endl;
}

void lunarEjecta_SpeedZenithIntegration::deleteGrid(grid* g) {
}