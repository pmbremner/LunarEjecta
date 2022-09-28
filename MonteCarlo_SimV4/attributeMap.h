#ifndef ATTRIBUTE_MAP_H
#define ATTRIBUTE_MAP_H

#include "dimension.h"

#include <vector>

using namespace std;

// virtual class, must create derived classes
class attributeMap
{
public:
	attributeMap();
	~attributeMap();

	// pure virtual function, derived class must implement this function
	virtual void define() = 0;
	void getNextState(); // wrapper to dimenstionCollection's getNextState calls
	void getCurrentState(vector<double>&); // returns stateStack

protected:
	vector<dimension*> dimensionCollection;
	vector<double>     stateStack;
	
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


class att_trajectory : public attributeMap
{
public:
	att_trajectory();
	~att_trajectory();

	void define();

private: 
	dimension* timeStamp;
	dimension* location;
	dimension* orientation;
	dimension* geometry;
};


class att_target : public attributeMap
{
public:
	att_target();
	~att_target();

	void define();

private: 
	dimension* location;
	dimension* density;
	dimension* porosity;
	dimension* strength;
	dimension* scalingLawParams;
};

class att_primaryEnvironmentCollection : public attributeMap
{
public:
	att_primaryEnvironmentCollection();
	~att_primaryEnvironmentCollection();

	void define();
	void getNextState(); // override default to handle multiple components

private: 
	attributeMap* priEnvCollection; // will be an array of att_primaryEnvironment
};



class att_primaryEnvironment : public attributeMap
{
public:
	att_primaryEnvironment();
	~att_primaryEnvironment();

	void define();

private: 
	dimension* fluxSpeedAngle;
	dimension* mass;
	dimension* density;
};

class att_craterEjectaEnvironment : public attributeMap
{
public:
	att_craterEjectaEnvironment();
	~att_craterEjectaEnvironment();

	void define();

private: 
	dimension* fluxSpeedAngleDistance;
	dimension* mass;
	dimension* density;
};


class att_assetEjectaEnvironment : public attributeMap
{
public:
	att_assetEjectaEnvironment();
	~att_assetEjectaEnvironment();

	void define();

private: 
	dimension* fluxSpeedAngleMass;
	dimension* density;
};


#endif 