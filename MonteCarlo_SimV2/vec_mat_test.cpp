#include <vector>
#include <iostream>

using namespace std;

struct vec3
{
	double x[3];
};

struct mat3x3
{
	vec3 col[3];
};


int main(int argc, char const *argv[])
{
	mat3x3 A;

	A.col[1].x[0] = 1.1;

	cout << A.col[1].x[0] << endl;

	return 0;
}