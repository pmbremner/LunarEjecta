#include "functionStrategy.h"
#include "rootSolverStrategy.h"
#include <iostream>

using namespace std;

int main(int argc, char const *argv[])
{
	//  v(number of values, the value)
	vector<double> v(1, 0.);

	functionStrategy* myTestFunc  = new myFunc(v);
	functionStrategy* myTestFunc2 = new myFunc2(v);


	cout << myTestFunc->execute(3.) << endl;
	cout << myTestFunc2->execute(3.) << endl;


	rootSolverStrategy* mySolver = new zbrent(myTestFunc);

	cout << mySolver->solve(0.5, 0., 6.28) << endl;

	mySolver->initFunction(myTestFunc2);

	//mySolver->updateMaximumIterations(3);

	cout << mySolver->solve(0.5, 0., 6.28) << endl;

	return 0;
}