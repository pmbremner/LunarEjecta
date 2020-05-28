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

struct nodes
{
	vector<set> node;
};

struct vEdges
{
	vector<set> vEdge;
};

struct hEdges
{
	vector<set> hEdge;
};

enum cellEdgeType {EL, ER, EU, ED};
enum cellNodeType {NUL, NUR, NDL, NDR};

struct cells
{
	set* cEdges[4];
	set* cNodes[4];

};


class lunarEjecta_SpeedZenithIntegration
{
public:
	lunarEjecta_SpeedZenithIntegration();
	~lunarEjecta_SpeedZenithIntegration();


private:
	vector<cells> region1Cells;
	vector<cells> region2Cells; // will always be only 1
	vector<cells> region3Cells;

	// 3 cases:
	//  Case 1: D0 and D1 < 0.5, all three regions valid
	//  Case 2: D0 < 0.5 and D1 >= 0.5, regions I and II valid
	//  Case 3: D0 and D1 >= 0.5, only region I valid
	int regionCase; 

	// dimensions
	double circumference; // km
	double escapeSpeed;   // km/s

	double D0; // distance 0
	double D1; // distance 1

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
