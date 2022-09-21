#ifndef DIMENSION_H
#define DIMENSION_H

using namespace std;

// virtual class, must create derived classes
class dimension
{
public:
	dimension();
	~dimension();

	// pure virtual function, derived class must implement this function
	virtual void getSample() = 0;

protected:
	vector<attribute*> element;

	vector<dim*> childDim;
	
};



#endif 