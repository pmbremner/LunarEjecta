#ifndef LUNAREJECTA_SPEEDZENITHINTEGRATION_H
#define LUNAREJECTA_SPEEDZENITHINTEGRATION_H

using namespace std;

#include <vector>
#include <cmath>
#include <iostream>

struct coord
{
	double x;
	double y;
};


struct set
{
	double dist; // units of lunar circumference
	coord loc;   // (x,y) | x = 1 - cos(zenith), y = vel [v_esc]
	vector<coord> P0;
	vector<coord> P1;
};

// EL = Edge Left
// ER = Edge Right
// EU = Edge Up
// ED = Edge Down
enum cellEdgeType {EL, ER, EU, ED};

// NUL = Node Up Left
// NUR = Node Up Right
// NDL = Node Down Left
// NDR = Node Down Right
enum cellNodeType {NUL, NUR, NDL, NDR};

struct grid
{
	vector<set> node;
	vector<set> vEdge;
	vector<set> hEdge;
};

struct cells
{
	vector<set> cell;
	vector<vector<set*>> cEdges;
	vector<vector<set*>> cNodes;
};


enum regionType {region_I, region_II, region_III};

class lunarEjecta_SpeedZenithIntegration
{
public:
	lunarEjecta_SpeedZenithIntegration(double new_D0, // input units of km
		                               double new_D1, // input units of km
		                               double new_radius,
		                               double new_Vesc,
		                               int new_Nx,
		                               int new_Nv);
	~lunarEjecta_SpeedZenithIntegration();

	inline int cIdx(int ix, int jv); // Nv-1, j inner loop, for cell and VE
	inline int gIdx(int ix, int jv); // Nv, j inner loop, for HE, NN

private:

	double H_calcVMin(double D);
	double H_calcX_at_vMin(double D);
	double H_calcDist(double x, double v);

	void initGrid(grid* g,
		          int Nx,
		          int Nv,
		          double xMin,
		          double xMax,
		          double vMin,
		          double vMax);

	void initCell(cells* c,
		          grid* g,
		          int Nx,
		          int Nv);

	void initSet(set* s, double x, double v);

	void printSetLoc(set* s);

	void deleteGrid(grid* g);

	cells region1Cells;
	cells region2Cells;
	cells region3Cells;

	grid region1Grid;
	grid region2Grid;
	grid region3Grid;

	// number of x divisions in each region
	int Nx_R1;
	int Nx_R2; // will always be 1
	int Nx_R3;

	int Nv;

	// 3 cases:
	//  Case 1: D0 and D1 < 0.5, all three regions valid
	//  Case 2: D0 < 0.5 and D1 >= 0.5, regions I and II valid
	//  Case 3: D0 and D1 >= 0.5, only region I valid
	int regionCase; 

	// dimensions
	double circumference; // km
	double escapeSpeed;   // km/s

	double D0; // distance 0, units of circumference
	double D1; // distance 1, units of circumference

	double xMin; // always the limiting angle at v = v_esc for D0
	double xMax; // always 1

	double vMin; // always the optimal speed for the closer distance D0
	double vMax; // always 1 [units of v_esc]

	// usage depends on regionCase
	double x_at_vMinD0;
	double x_at_vMinD1;

	double vMinD0;
	double vMinD1;

};

#endif 
