#ifndef SAMPLE_DIRECTOR_H
#define SAMPLE_DIRECTOR_H

#include "stateMap.h"

#include <vector>

using namespace std;

class sampleDirector
{
public:
	sampleDirector(string );
	~sampleDirector();

	void sampleAssetEnvironment(string ); 
	void saveSamplesToFile(string );

private:
	void getNextState();
	void getCurrentState(vector<double>&);


	vector<stateMap*> stateMapCollection;
	stateMap*         assetEjectaEnvironmentState;
};


#endif 