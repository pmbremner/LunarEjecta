#ifndef STATE_DIRECTOR_H
#define STATE_DIRECTOR_H

#include "stateMap.h"

#include <vector>

using namespace std;

class stateDirector
{
public:
	stateDirector();
	~stateDirector();

	void buildStateDefinitions(string );

private:
	vector<stateMap*> stateCollection;
	stateMap*         trajectoryState;
	stateMap*         targetState;
	stateMap*         primaryEnvironmentState;
	stateMap*         craterEjectaEnvironmentState;
};


#endif 