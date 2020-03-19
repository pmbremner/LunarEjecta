// Declaration of 2D vector

#ifndef VECTOR2D_H
#define VECTOR2D_H

#include <vector>

using namespace std;


class vector2d
{
public:
	vector2d();
	~vector2d();
	
};



typedef vector<vector<double> > vec2d;


void init_vec2d(vec2d v2, int rows, int cols);


void print_vec2d(vec2d v2, int rows, int cols);



#endif 