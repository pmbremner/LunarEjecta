#include "LunarEjecta_speedzenith.h"
#include "LunarEjecta_params.h"

//#include <float.h>
#include <cmath>

using namespace std;


// x is zenith angle in rads, y is speed in units of escape speed, return is distance in units of radii
double getD(double x, double y)
{
	return 2.*atan2(sqr(y)*sin(2.*x), 1.-2.*sqr(y*sin(x)));
}

void push(vector<quad> &q, double x0, double x1, double y0, double y1, double w)
{
	q.resize(q.size()+1);
	q[q.size()-1].x0 = x0;
	q[q.size()-1].x1 = x1;
	q[q.size()-1].y0 = y0;
	q[q.size()-1].y1 = y1;
	q[q.size()-1].w  = w;
}


void getEvalQuads(vector<quad> &q, double D0, double D1, double x0, double x1, double y0, double y1, int level_max, int level)
{

	vector<double> D;
	D.push_back(getD(x0, y0));
	D.push_back(getD(x0, y1));
	D.push_back(getD(x1, y0));
	D.push_back(getD(x1, y1));

	double D_min = vMin(D);
	double D_max = vMax(D);

	//cout << level << ' ' << D_max << ' ' << D_min << endl;

	if (D_max <= D0 || D_min >= D1) // quad is below or above region completely
	{
		return;
	}
	else if(D_min >= D0 && D_max <= D1) // quad is inside region completely
	{
		push(q, x0, x1, y0, y1, 1.);
	}
	else if(level < level_max) // the quad intersects the region, subdivide
	//else if((D_max - D_min) / (D1 - D0) > 0.2) // the quad intersects the region, subdivide
	{
		getEvalQuads(q, D0, D1, x0        , (x0+x1)/2., y0        , (y0+y1)/2., level_max, level+1); // bottom left
		getEvalQuads(q, D0, D1, (x0+x1)/2., x1        , y0        , (y0+y1)/2., level_max, level+1); // bottom right
		getEvalQuads(q, D0, D1, x0        , (x0+x1)/2., (y0+y1)/2., y1        , level_max, level+1); // top left
		getEvalQuads(q, D0, D1, (x0+x1)/2., x1        , (y0+y1)/2., y1        , level_max, level+1); // top right
	}
	else // too many subdivisions, but the quad still intersects. Still save point, but with an approximate weight
	{
		push(q, x0, x1, y0, y1, (min(D1, D_max) - max(D0, D_min)) / (D1 - D0));
	}
}

void printQuad(vector<quad> &q)
{
	for (int i = 0; i < q.size(); ++i)
	{
		cout << (q[i].x0 + q[i].x1) / 2. << ' ' << (q[i].y0 + q[i].y1) / 2. << ' ' << q[i].w << ' ' << (q[i].x1 - q[i].x0) * (q[i].y1 - q[i].y0) * q[i].w << endl;
	}
}

double getSumQuad(vector<quad> &q)
{
	double sum = 0.;
	for (int i = 0; i < q.size(); ++i)
		sum += (q[i].x1 - q[i].x0) * (q[i].y1 - q[i].y0) * q[i].w;
	return sum;
}