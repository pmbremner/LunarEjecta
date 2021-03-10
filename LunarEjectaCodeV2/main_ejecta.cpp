#include <iostream>
#include <fstream>
#include <stdlib.h>     /* atoi */


using namespace std;

int main(int argc, char const *argv[])
{
	int N_proc = atoi(argv[1]); // total number of processes
	int i_proc = atoi(argv[2]); // current process number
	int region_start = atoi(argv[3]); // start index of lat-lon files
	int region_end   = atoi(argv[4]); // end index of lat-lon files



	cout << "Total # of processes = " << N_proc << endl;
	cout << "Current process #    = " << i_proc << endl;
	cout << "  Start and end region index = " << region_start << ", " << region_end << endl;

	for (int i = 0; i < 1000000000000/(i_proc+1); ++i)
	{
		N_proc ++;
	}

	cout << i_proc << " finished\n";


	return 0;
}