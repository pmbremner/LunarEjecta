#ifndef STATE_MAP_H
#define STATE_MAP_H

#include "attributeMap.h"

#include <string>

using namespace std;

// virtual class, must create derived classes
class stateMap
{
public:
	stateMap();
	~stateMap();

	// pure virtual function, derived class must implement this function
	virtual void define() = 0;
	void getNextState(); // wrapper to attributeMap's getNextState calls, stores in stateStackCollection
	void getCurrentState(vector<double>&); // returns copy of stateStackCollection

protected:
	string         name;
	attributesMap* attributes;
	vector<double> stateStackCollection;
};


#endif 