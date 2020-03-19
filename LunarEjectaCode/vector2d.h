// Declaration of 2D vector

#ifndef VECTOR2D_H
#define VECTOR2D_H

#include <vector>

using namespace std;

class vector2d
{
public:
	vector2d(int rows, int cols);
	~vector2d();

private:
	template<typename T>
	vector<vector<T>> vec2d;

};


#endif 