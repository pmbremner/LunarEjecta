#include "stateDirector.h"

using namespace std;

stateDirector::stateDirector()
{

}


stateDirector::~stateDirector()
{
	for (vector<stateMap*>::iterator i = stateMapCollection.begin(); i != stateMapCollection.end(); ++i)
		delete (*i);
}




void stateDirector::buildStateDefinitions(string filenameSetup)
{
	for (vector<stateMap*>::iterator i = stateMapCollection.begin(); i != stateMapCollection.end(); ++i)
		{
			
		}
}