#include "functionStrategy.h"
#include <iostream>

using namespace std;

int main(int argc, char const *argv[])
{
	
	myData* myTestData = new myData(1.1, -3.4);
	functionStrategy* myTestFunc = new testFunc(myTestData);


	cout << myTestFunc->execute(3.) << endl;

	return 0;
}