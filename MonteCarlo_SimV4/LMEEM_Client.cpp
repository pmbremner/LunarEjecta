#include "LMEEM_headers.h"

using namespace std;


int main(int argc, char const *argv[])
{
	stateDirector* LMEEM_stateDirector = new stateDirector();

	LMEEM_stateDirector.buildStateDefinitions("initLMEEM.txt");

	/* code */
	return 0;
}