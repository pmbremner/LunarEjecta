#ifndef ATTRIBUTE_H
#define ATTRIBUTE_H

using namespace std;


// virtual class, must create derived classes
class attribute
{
public:
	attribute();
	~attribute();
	
	// pure virtual function, derived class must implement this function
	virtual void getSample() = 0;

	// data members
	string name;

};






class massDistribution: public attribute
{
public:
	massDistribution();
	~massDistribution();

	void getSample();

protected:
	cdf* massCDF;
};


class densityDistribution: public attribute
{
public:
	densityDistribution();
	~densityDistribution();

	void getSample();

protected:
	pdf* densityPDF;
}


class speedAngleDistanceDistribution: public attribute
{
public:
	speedAngleDistanceDistribution();
	~speedAngleDistanceDistribution();

	void getSample();

protected:
	
	
};





#endif 