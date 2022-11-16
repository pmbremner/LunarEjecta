#include <vector>
#include <iostream>

using namespace std;



class A
{
public:
	//A();
	virtual ~A() = default;
	
	vector<double> v;

};

class B : public A
{
public:
	enum var{radius, mass, speed};

	B(double x1, double x2);
	~B();
	
	double a;
	
	double b;

};

B::B(double x1, double x2){
	a = x1;
	b = x2;

	v.push_back(x1);
	v.push_back(x2);
}

B::~B() {}



int main(int argc, char const *argv[])
{
	B* myTest = new B(1.1, -0.9);
	A* myTest2 = myTest;

	cout << myTest->a << ' ' << myTest->b << endl;
	cout << myTest2->v[0] << ' ' << myTest2->v[1] << endl;
	cout << myTest2->v[B::radius] << ' ' << myTest2->v[B::mass] << endl;


	return 0;
}