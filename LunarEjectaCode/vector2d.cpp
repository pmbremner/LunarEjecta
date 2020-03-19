#include "vector2d.h"
#include <vector>
#include <iostream>

using namespace std;

void init_vec2d(vec2d v2, int rows, int cols)
{
	v2.resize(rows, vector<double>(cols, 0.0));
	cout << v2.size() << endl;
}

void print_vec2d(vec2d v2, int rows, int cols)
{
	int i, j;
	for (i = 0; i < rows; ++i)
	{
		for (j = 0; j < cols; ++j)
		{
			cout << v2[i][j] << ' ';
		}
		cout << endl;
	}
}