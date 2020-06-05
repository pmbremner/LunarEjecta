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

		cout << " Init cells in Region 1: \n";
		initCell(&region1Cells, &region1Grid, Nx_R1, Nv);

		cout << " Init grid in Region 2: \n";
		cout << "   D0 = " << D0 << endl; 
		cout << "   D1 = " << D1 << endl; 
		initGrid(&region2Grid, Nx_R2, Nv, x_at_vMinD0, x_at_vMinD1, vMinD0, vMax);

		cout << " Init cells in Region 2: \n";
		initCell(&region2Cells, &region2Grid, Nx_R2, Nv);

		cout << " Init grid in Region 3: \n";
		cout << "   D0 = " << D0 << endl; 
		cout << "   D1 = " << D1 << endl; 
		initGrid(&region3Grid, Nx_R3, Nv, x_at_vMinD1, xMax, vMinD1, vMax);

		cout << " Init cells in Region 3: \n";
		initCell(&region3Cells, &region3Grid, Nx_R3, Nv);
	}
	else if (regionCase == region_II)
	{
		cout << " Init grid in Region 1: \n";
		cout << "   D0 = " << D0 << endl; 
		cout << "   D1 = " << D1 << endl; 
		initGrid(&region1Grid, Nx_R1, Nv, xMin, x_at_vMinD0, vMinD0, vMax);

		cout << " Init cells in Region 1: \n";
		initCell(&region1Cells, &region1Grid, Nx_R1, Nv);

		cout << " Init grid in Region 2: \n";
		cout << "   D0 = " << D0 << endl; 
		cout << "   D1 = " << D1 << endl; 
		initGrid(&region2Grid, Nx_R2, Nv, x_at_vMinD0, x_at_vMinD1, vMinD0, vMax);

		cout << " Init cells in Region 2: \n";
		initCell(&region2Cells, &region2Grid, Nx_R2, Nv);

	}
	else if (regionCase == region_III)
	{
		cout << " Init grid in Region 1: \n";
		cout << "   D0 = " << D0 << endl; 
		cout << "   D1 = " << D1 << endl; 
		initGrid(&region1Grid, Nx_R1, Nv, xMin, x_at_vMinD0, vMinD0, vMax);

		cout << " Init cells in Region 1: \n";
		initCell(&region1Cells, &region1Grid, Nx_R1, Nv);
	}

 	// init cells for each region

	// region1Cells.resize((Nx_R1 - 1) * (Nv - 1));
	// region2Cells.resize((Nx_R2 - 1) * (Nv - 1));
	// region3Cells.resize((Nx_R3 - 1) * (Nv - 1));

	// insert the intersection points into the edges
	if (regionCase == region_I)
	{
		insertD(&region1Grid, 0.1*circumference, 0, region_I);
		insertD(&region2Grid, 0.1*circumference, 0, region_II);
		insertD(&region3Grid, 0.1*circumference, 0, region_III);

		computeIntersectionType(&region1Cells, 0, region_I);
		computeIntersectionType(&region2Cells, 0, region_II);
		computeIntersectionType(&region3Cells, 0, region_III);
	}
	else if (regionCase == region_II)
	{
		insertD(&region1Grid, 0.1*circumference, 0, region_I);
		insertD(&region2Grid, 0.1*circumference, 0, region_II);

		computeIntersectionType(&region1Cells, 0, region_I);
		computeIntersectionType(&region2Cells, 0, region_II);
	}
	else if (regionCase == region_III)
	{
		insertD(&region1Grid, 0.1*circumference, 0, region_I);

		computeIntersectionType(&region1Cells, 0, region_I);
	}
	

}

lunarEjecta_SpeedZenithIntegration::~lunarEjecta_SpeedZenithIntegration() {}



inline int lunarEjecta_SpeedZenithIntegration::cIdx(int ix, int jv) {
	return jv + (Nv-1) * ix;
}

inline int lunarEjecta_SpeedZenithIntegration::gIdx(int ix, int jv) {
	return jv + Nv * ix;
}


// will loop through all edge x and v locations
void lunarEjecta_SpeedZenithIntegration::insertD
            (grid* g,
             double D0km /* in units of km */,
             int iD /* 0 -> D0, 1 - > D1 */,
	         int iR /* i -> region i */ )
{
	int i, j, Nx = 0, idx_min, idx_mid, idx_max;
	double xi, vj;
	double xD, vD;
	coord P;

	double D0 = D0km / circumference; // domain is now D0 -> [0, 1]

	cout << " Inserting intersection points in region ";

	if (iR == region_I)
	{
		Nx = Nx_R1;
		cout << "1...\n";
	}
	else if (iR == region_II)
	{
		Nx = Nx_R2;
		cout << "2...\n";
	}
	else if (iR == region_III)
	{
		Nx = Nx_R3;
		cout << "3...\n";
	}
	else {
		cout << "ERROR: lunarEjecta_SpeedZenithIntegration::insertD, invalid region\n";
	}

	// loop through verticle edge x locations
	// (The same logic regardless the region)
	for (i = 0; i < Nx; ++i)
	{
		xi = g->vEdge[cIdx(i, 0)].loc.x;

		cout << "\ncheck ... " << g->vEdge[cIdx(i, 0)].loc.x << ' ' << g->vEdge[cIdx(i, 1)].loc.x << endl;

		vD = H_calcSpeed(xi, D0); // units of v_esc

		// find index j to insert D position (binary search)
		// see https://en.wikipedia.org/wiki/Binary_search_algorithm
		idx_min = 0;
		idx_max = Nv-1; // # of nodes minus 1 in v-direction for idx

		cout << i << ' ' << idx_min << ' ' << idx_max << "   ";
		cout << vD << ' ' << g->node[gIdx(i, idx_min)].loc.y << ' ' << g->node[gIdx(i, idx_max)].loc.y << endl;
		printSetLoc(&g->node[gIdx(i, idx_min)]);
		printSetLoc(&g->node[gIdx(i, idx_max)]);

		// need to check node vals since vEdge for y is the center of the grid, not edge of grid
		if (vD >= g->node[gIdx(i, idx_min)].loc.y && vD <= g->node[gIdx(i, idx_max)].loc.y)
		{ // Insert only if vD is in the domain of the v's

			while (idx_min <= idx_max) {
				idx_mid = (idx_min + idx_max) / 2;

				if (g->node[gIdx(i, idx_mid)].loc.y < vD)
				{
					idx_min = idx_mid + 1;
				}
				else
				{
					idx_max = idx_mid - 1;
				}
			}
			j = idx_mid-1; // need to set back 1, found node val > vD

			P.x = xi;
			P.y = vD;

			g->vEdge[cIdx(i, j)].P[iD].push_back(P);

			cout << "vEdge " << i << ' ' << j << ' ' << g->vEdge[cIdx(i, j)].P[iD].size() << endl;
			cout << " vD = " << vD << " | ";
			cout << g->node[gIdx(i, j)].loc.y << ' ' << g->node[gIdx(i, j+1)].loc.y << " | ";
			cout << P.x << ' ' << P.y << endl;
			printSetLoc(&g->vEdge[cIdx(i, j)]);

			// // insert crossing point into verticle edge (i,j)
			// if(iD == 0)
			// {
			// 	g->vEdge[cIdx(i, j)].P0.push_back(P);
			// }
			// else if (iD == 1)
			// {
			// 	g->vEdge[cIdx(i, j)].P1.push_back(P);
			// } 
			// else
			// {
			// 	cout << "ERROR: lunarEjecta_SpeedZenithIntegration::insertD, invalid distance index\n";
			// }
		} // END IF, in domain
	} // END FOR, loop through verticle edges


	// loop through horizontal edges (depends on the region we're in)
	for (j = 0; j < Nv; ++j)
	{
		vj = g->hEdge[gIdx(0, j)].loc.y;

		// compute x position of vj and D
		if (iR == region_I)
		{
			xD = H_calcXp(vj, D0);
		}
		else if (iR == region_II)
		{
			if(iD == 0)
			{
				xD = H_calcXm(vj, D0);
			}
			else if (iD == 1)
			{
				xD = H_calcXp(vj, D0);
			} 
			else
			{
				cout << "ERROR: lunarEjecta_SpeedZenithIntegration::insertD, invalid distance index\n";
			}
		}
		else if (iR == region_III)
		{
			xD = H_calcXm(vj, D0);
		}
		else {
			cout << "ERROR: lunarEjecta_SpeedZenithIntegration::insertD, invalid region\n";
		}

		// find index i to insert D position (binary search)
		idx_min = 0;
		idx_max = Nx-1;

		cout << "check hEdge...\n";
		cout << j << ' ' << idx_min << ' ' << idx_max << "   ";
		cout << xD << ' ' << g->node[gIdx(idx_min, j)].loc.x << ' ' << g->node[gIdx(idx_max, j)].loc.x << endl;
		printSetLoc(&g->node[gIdx(idx_min, j)]);
		printSetLoc(&g->node[gIdx(idx_max, j)]);

		if (xD >= g->node[gIdx(idx_min, j)].loc.x && xD <= g->node[gIdx(idx_max, j)].loc.x)
		{ // Insert only if xD is in the domain of the x's

			while (idx_min <= idx_max) {
				idx_mid = (idx_min + idx_max) / 2;

				cout << "  -> " << g->node[gIdx(idx_mid, j)].loc.x  << ' ' << xD;
				cout << "  | idx_mid = " << idx_mid << endl;
				if (g->node[gIdx(idx_mid, j)].loc.x < xD)
				{
					idx_min = idx_mid + 1;
				}
				else
				{
					idx_max = idx_mid - 1;
				}
			}
			i = idx_mid-1;

			P.x = xD;
			P.y = vj;

			// insert crossing point into verticle edge (i,j)
			g->hEdge[gIdx(i, j)].P[iD].push_back(P);

			cout << "hEdge " << i << ' ' << j << ' ' << g->hEdge[gIdx(i, j)].P[iD].size() << endl;
			cout << " xD = " << xD << " | ";
			cout << g->node[gIdx(i, j)].loc.x << ' ' << g->node[gIdx(i+1, j)].loc.x << " | ";
			cout << P.x << ' ' << P.y << endl;
			printSetLoc(&g->hEdge[gIdx(i, j)]);

			// if(iD == 0)
			// {
			// 	g->hEdge[gIdx(i, j)].P0.push_back(P);
			// }
			// else if (iD == 1)
			// {
			// 	g->hEdge[gIdx(i, j)].P1.push_back(P);
			// } 
			// else
			// {
			// 	cout << "ERROR: lunarEjecta_SpeedZenithIntegration::insertD, invalid distance index\n";
			// }

		} // END IF, in domain
	} // END FOR, loop through horizontal edges

}

// will loop through all cells
//  Must do insertD first
void lunarEjecta_SpeedZenithIntegration::computeIntersectionType
             (cells* c,
              int iD /* 0 -> D0, 1 - > D1 */,
	          int iR /* i -> region i */ )
{
	int i, j;
	double Nx, x, v;

	cout << " Computing intersection type in region ";

	if (iR == region_I)
	{
		Nx = Nx_R1;
		cout << "1...\n";
	}
	else if (iR == region_II)
	{
		Nx = Nx_R2;
		cout << "2...\n";
	}
	else if (iR == region_III)
	{
		Nx = Nx_R3;
		cout << "3...\n";
	}
	else {
		cout << "ERROR: lunarEjecta_SpeedZenithIntegration::insertD, invalid region\n";
	}

	// loop through the cells checking the intersection points
	for (i = 0; i < Nx - 1; ++i)
	{
		x = c->cell[cIdx(i, 0)].loc.x;

		for (j = 0; j < Nv - 1; ++j)
		{ 
			v = c->cell[cIdx(0, j)].loc.y;

			printSetLoc(&(c->cell[cIdx(i, j)]));

			H_getIntersectionType(c, i, j, iD, iR);

		}
	}
}

// Note: we didn't store any intersection data in the nodes, even though we inteneded to
// Instead, the algorithm (previoud to this one) will store two edges that really are the same point
// so we will have to check this to get the cases right, otherwise the intergration
// part of the algorithm will probably break down for a zero-sized region (the slope might be infinity)
void lunarEjecta_SpeedZenithIntegration::H_getIntersectionType(cells* c, int i, int j, int iD, int iR)
{
	int idx = cIdx(i, j);

	// reflect left and right depending on which region (iR) and which curve (iD)
	// setting up the edge pointers and node pointers logic
	//
	// Note, for some reason, we can't just use set* .......
	vector<set*> t_EL(1), t_ER(1), t_EU(1), t_ED(1);
	vector<set*> t_NUL(1), t_NUR(1), t_NDL(1), t_NDR(1);

	if (iR == region_I || (iR == region_II && iD == 1))
	{
		t_EL[0] = c->cEdges[idx][EL];
		t_ER[0] = c->cEdges[idx][ER];
		t_EU[0] = c->cEdges[idx][EU];
		t_ED[0] = c->cEdges[idx][ED];

		t_NUL[0] = c->cNodes[idx][NUL];
		t_NUR[0] = c->cNodes[idx][NUR];
		t_NDL[0] = c->cNodes[idx][NDL];
		t_NDR[0] = c->cNodes[idx][NDR]; 
	}
	else if (iR == region_III || (iR == region_II && iD == 0))
	{
		t_EL[0] = c->cEdges[idx][ER]; // reflect
		t_ER[0] = c->cEdges[idx][EL]; // reflect
		t_EU[0] = c->cEdges[idx][EU];
		t_ED[0] = c->cEdges[idx][ED];

		t_NUL[0] = c->cNodes[idx][NUR]; // reflect
		t_NUR[0] = c->cNodes[idx][NUL]; // reflect
		t_NDL[0] = c->cNodes[idx][NDR]; // reflect
		t_NDR[0] = c->cNodes[idx][NDL]; // reflect
	}
	else
	{
		cout << "ERROR: lunarEjecta_SpeedZenithIntegration::H_getIntersectionType, invalid region or distance index\n";
	}



	/////////////////////////////////////////
	// note that we could have made this logic more compact, but it would be harder to read..
	// so we're doing it the "longer" way for clarity 
	//
	// We first need to check the size of the arrays
	//
	//// check NUR_LT_D : 0 ////
	// We will assume if the upper edge equals the node, then the right edge does too
	cout << " EL | " << i << ' ' << j << ' ' << t_EL[0]->P[iD].size() << endl;
	cout << " ER | " << i << ' ' << j << ' ' << t_ER[0]->P[iD].size() << endl;
	cout << " EU | " << i << ' ' << j << ' ' << t_EU[0]->P[iD].size() << endl;
	cout << " ED | " << i << ' ' << j << ' ' << t_ED[0]->P[iD].size() << endl;

	cout << " NUL pos = " << t_NUL[0]->loc.x << ' ' << t_NUL[0]->loc.y << endl;
	cout << " NUR pos = " << t_NUR[0]->loc.x << ' ' << t_NUR[0]->loc.y << endl;
	cout << " NDL pos = " << t_NDL[0]->loc.x << ' ' << t_NDL[0]->loc.y << endl;
	cout << " NDR pos = " << t_NDR[0]->loc.x << ' ' << t_NDR[0]->loc.y << endl;

	if( t_EL[0]->P[iD].size() > 0 ) {
		cout << "  Intersection Point = " << t_EL[0]->P[iD][0].x << ' ' << t_EL[0]->P[iD][0].y << endl;
	}
	if( t_ER[0]->P[iD].size() > 0 ) {
		cout << "  Intersection Point = " << t_ER[0]->P[iD][0].x << ' ' << t_ER[0]->P[iD][0].y << endl;
	}
	if( t_EU[0]->P[iD].size() > 0 ) {
		cout << "  Intersection Point = " << t_EU[0]->P[iD][0].x << ' ' << t_EU[0]->P[iD][0].y << endl;
	}
	if( t_ED[0]->P[iD].size() > 0 ) {
		cout << "  Intersection Point = " << t_ED[0]->P[iD][0].x << ' ' << t_ED[0]->P[iD][0].y << endl;
	}

	// if (t_EU[0]->dist == t_NUR[0]->dist)
	// {
	// 	/* code */
	// }
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

double lunarEjecta_SpeedZenithIntegration::H_calcSpeed(double x, double D){
	return pow(2.*(1.-x)*sqrt(x*(2.-x)) / tan(D * PI) + 2.*x*(2.-x), -0.5);
}


double lunarEjecta_SpeedZenithIntegration::H_calcXp(double v, double D) { // plus solution
	return 1. - sqrt(HH_calcCos2alphap(v, D));
}

double lunarEjecta_SpeedZenithIntegration::H_calcXm(double v, double D){ // minus solution
	return 1. - sqrt(HH_calcCos2alpham(v, D));
}

double lunarEjecta_SpeedZenithIntegration::HH_calcCos2alphap(double v, double D) {
	double tan2_2v_1_term = sqr(tan(D * PI))*(2.*sqr(v)-1.);
	return ( sqr(v) + tan2_2v_1_term + sqrt(sqr(sqr(v)) + tan2_2v_1_term) ) / ( 2.*sqr(v)*(1. + sqr(tan(D * PI))) );
}

double lunarEjecta_SpeedZenithIntegration::HH_calcCos2alpham(double v, double D) {
	double tan2_2v_1_term = sqr(tan(D * PI))*(2.*sqr(v)-1.);
	return ( sqr(v) + tan2_2v_1_term - sqrt(sqr(sqr(v)) + tan2_2v_1_term) ) / ( 2.*sqr(v)*(1. + sqr(tan(D * PI))) );
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
	cout << " Nx by Nv = " << Nx << ' ' << Nv << endl;
	for (i = 0; i < Nx; ++i)
	{
		x = xMin + (xMax - xMin) * i / double(Nx - 1);

		for (j = 0; j < Nv; ++j)
		{
			v = vMin + (vMax - vMin) * j / double(Nv - 1);

			// set the node position (x, v) and corresponding distance
			initSet(&(g->node[gIdx(i, j)]), x, v);
		}
	}

	// loop through grid verticle edges
	cout << "...init verticle edges...\n";
	cout << " Nx by Nv = " << Nx << ' ' << Nv-1 << endl;
	for (i = 0; i < Nx; ++i)
	{
		x = xMin + (xMax - xMin) * i / double(Nx - 1);

		for (j = 0; j < Nv - 1; ++j)
		{   //  the v is the center of the edge
			v = vMin + (vMax - vMin) * (j + 0.5) / double(Nv - 1);

			// set the vert edges position (x, v) and corresponding distance
			initSet(&(g->vEdge[cIdx(i, j)]), x, v);
			cout << i << ' ' << j << ' ';
			printSetLoc(&(g->vEdge[cIdx(i, j)]));
		}
	}

	// loop through grid horizonal edges
	cout << "...init horizontal edges...\n";
	cout << " Nx by Nv = " << Nx-1 << ' ' << Nv << endl;
	for (i = 0; i < Nx - 1; ++i)
	{//  the x is the center of the edge
		x = xMin + (xMax - xMin) * (i + 0.5) / double(Nx - 1);

		for (j = 0; j < Nv; ++j)
		{   
			v = vMin + (vMax - vMin) * j / double(Nv - 1);

			// set the horz edges position (x, v) and corresponding distance
			initSet(&(g->hEdge[gIdx(i, j)]), x, v);
		}
	}
}



void lunarEjecta_SpeedZenithIntegration::initCell
			(cells* c,
		     grid* g,
		     int Nx, // # of nodes
		     int Nv) // # of nodes
{
	double x, v;
	int i, j;
	// cell set data and edge and node pointers
	c->cell.resize(  (Nx - 1) * (Nv - 1));
	c->cEdges.resize((Nx - 1) * (Nv - 1));
	c->cNodes.resize((Nx - 1) * (Nv - 1));

	c->intsecTypeD0.resize((Nx - 1) * (Nv - 1), 0);
	c->intsecTypeD1.resize((Nx - 1) * (Nv - 1), 0);

	// loop through cells, linking the pointers for nodes and edges in each cell
	cout << "...init cells...\n";
	cout << " Nx by Nv = " << Nx-1 << ' ' << Nv-1 << endl;
	for (i = 0; i < Nx - 1; ++i)
	{
		x = g->hEdge[gIdx(i, 0)].loc.x;

		for (j = 0; j < Nv - 1; ++j)
		{   //  the v is the center of the edge
			v = g->vEdge[cIdx(0, j)].loc.y;

			// set the cell position (x, v) and corresponding distance
			initSet(&(c->cell[cIdx(i, j)]), x, v);

			// init the 4 pointers for each the nodes and edges
			c->cNodes[cIdx(i, j)].resize(4);
			c->cEdges[cIdx(i, j)].resize(4);

			// link nodes
			c->cNodes[cIdx(i, j)][NUL] = &(g->node[gIdx(i  , j+1)]);
			c->cNodes[cIdx(i, j)][NUR] = &(g->node[gIdx(i+1, j+1)]);
			c->cNodes[cIdx(i, j)][NDL] = &(g->node[gIdx(i  , j  )]);
			c->cNodes[cIdx(i, j)][NDR] = &(g->node[gIdx(i+1, j  )]);

			// // check locations
			// cout << "...checking linked node locations...\n";
			// printSetLoc(c->cNodes[cIdx(i, j)][NUL]);
			// printSetLoc(c->cNodes[cIdx(i, j)][NUR]);
			// printSetLoc(c->cNodes[cIdx(i, j)][NDL]);
			// printSetLoc(c->cNodes[cIdx(i, j)][NDR]);

			// link edges
			c->cEdges[cIdx(i, j)][EL] = &(g->vEdge[cIdx(i  , j  )]);
			c->cEdges[cIdx(i, j)][ER] = &(g->vEdge[cIdx(i+1, j  )]);
			c->cEdges[cIdx(i, j)][EU] = &(g->hEdge[gIdx(i  , j+1)]);
			c->cEdges[cIdx(i, j)][ED] = &(g->hEdge[gIdx(i  , j  )]);

			// // check locations
			// cout << "...checking linked edge locations...\n";
			// printSetLoc(c->cEdges[cIdx(i, j)][EL]);
			// printSetLoc(c->cEdges[cIdx(i, j)][ER]);
			// printSetLoc(c->cEdges[cIdx(i, j)][EU]);
			// printSetLoc(c->cEdges[cIdx(i, j)][ED]);
		}
	}
}


void lunarEjecta_SpeedZenithIntegration::initSet(set* s, double x, double v) {
	s->loc.x = x;
	s->loc.y = v;

	s->dist = H_calcDist(x, v);

	s->P.resize(2); // for P0 and P1
	//cout << " x, v, D | " << x << ' ' << v << ' ' << s->dist << endl;
}

void lunarEjecta_SpeedZenithIntegration::printSetLoc(set* s) {
	cout << " x = " << s->loc.x << '\t' << s->loc.y << endl;
}

void lunarEjecta_SpeedZenithIntegration::deleteGrid(grid* g) {
}