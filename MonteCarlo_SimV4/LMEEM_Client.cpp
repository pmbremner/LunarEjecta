#include "LMEEM_headers.h"

using namespace std;


int main(int argc, char const *argv[])
{
	stateDirector* LMEEM_stateDirector   = new stateDirector();
	sampleDirector* LMEEM_sampleDirector = new sampleDirector("defineLMEEM_sampleBins.txt");

	LMEEM_stateDirector.buildStateDefinitions("initLMEEM.txt");

	LMEEM_sampleDirector.sampleAssetEnvironment("defineLMEEM_sampleParams.txt");

	LMEEM_sampleDirector.saveSamplesToFile("run0");

	/* code */
	return 0;
}