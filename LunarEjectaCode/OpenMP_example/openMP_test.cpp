#include <iostream>
#include "omp.h"

using namespace std;

int main(int argc, char const *argv[])
{
	#pragma omp parallel
	{
		cout << "Hello World\n";
	}
	return 0;
}