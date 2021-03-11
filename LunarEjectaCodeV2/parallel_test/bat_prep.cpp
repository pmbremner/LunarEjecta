#include <iostream>
#include <fstream>
#include <stdlib.h>     /* atoi */
#include <math.h>       /* ceil */


using namespace std;

int main(int argc, char const *argv[])
{
	
	int N_proc = atoi(argv[1]); // total number of processes
	int N_lat  = atoi(argv[2]); // number of latitude files per longitude
	int N_lon  = atoi(argv[3]); // number of longitudes
	int N_tot  = N_lat * N_lon;

	int recompile = atoi(argv[4]); // 0 = no, 1 = yes

	int N_rpp  = ceil(N_tot / double(N_proc)); // number of regions per process
	int region_start, region_end;

	// cout << "Process ID = " << ProcID << endl;
	ofstream batfile;
	batfile.open("setup_parallel_ejecta.bat");

	batfile << "@echo off\n";
	if (recompile)
		batfile << "g++ -O2 main_ejecta.cpp -o MeMoSeE.exe\n";
	for (int i_proc = 0; i_proc < N_proc; ++i_proc)
	{
		region_start = i_proc * N_rpp;
		region_end   = fmin((i_proc + 1) * N_rpp - 1, N_tot - 1);

		batfile << "start MeMoSeE.exe " << N_proc << ' ' << i_proc << ' ' << region_start << ' ' << region_end << endl;
	}
	


	return 0;
}