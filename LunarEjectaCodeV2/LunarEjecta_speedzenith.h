#ifndef LUNAREJECTA_SPEED_ZENITH_MAP_H
#define LUNAREJECTA_SPEED_ZENITH_MAP_H

#include <vector>

using namespace std;


struct quad
{
	double x0;
	double x1;
	double y0;
	double y1;
	double w; // weight (0, 1]
};

void getEvalQuads(vector<quad> &q, double D0, double D1, double x0, double x1, double y0, double y1, int level_max, int level);

void printQuad(vector<quad> &q);

double getSumQuad(vector<quad> &q);

#endif 