#include "read_igloo_files.h"
#include "read_parameters.h"

#include <stdlib.h>     /* atoi */
#include <string>
#include <chrono>

using namespace std;

int main(int argc, char const *argv[])
{
	int i, j, idx_latlon;
	double lat_cur, runtime, percentFinished;

	string param_fn = argv[1];
	int N_proc = atoi(argv[2]);
	int i_proc = atoi(argv[3]);

	cout << "Process number = " << i_proc << " / " << N_proc << endl << endl;

	int Nlat, Nlon, NlatlonRegions;

	getParam(param_fn, "Nlat", Nlat, 0);
	getParam(param_fn, "Nlon", Nlon, 0);

	NlatlonRegions = Nlat * Nlon;

	string directory_name[Nlon];
	for (i = 0; i < Nlon; ++i)
		getParam(param_fn, "MEMDirectoryName_" + to_string(i * 360/Nlon), directory_name[i], 0);

	string NEOvel_fn;
	getParam(param_fn, "NEOvel_fn", NEOvel_fn, 0);

	double NEO_mass_min;
	double NEO_mass_max;
	getParam(param_fn, "NEO_mass_min", NEO_mass_min, 0);
	getParam(param_fn, "NEO_mass_max", NEO_mass_max, 0);
	
	// getParam(param_fn, "MEMDirectoryName_90", directory_name, 0);
	// getParam(param_fn, "MEMDirectoryName_180", directory_name, 0);
	// getParam(param_fn, "MEMDirectoryName_270", directory_name, 0);

	// all fluxes in units of #/m^2/yr
	vector<vector<double>> iglooDataHi;
	vector<vector<double>> iglooDataLo;
	vector<vector<double>> iglooDataNEO;

	iglooDataHi.resize(NlatlonRegions);
	iglooDataLo.resize(NlatlonRegions);
	iglooDataNEO.resize(NlatlonRegions);

	vector<double> NEOvel;
	readVelDist(NEOvel_fn, NEOvel);

	// read input igloo files
	auto t_read_latlon1 = std::chrono::high_resolution_clock::now();
	for (i = 0; i < Nlon; ++i)
	{
		for (j = 0; j < Nlat; ++j)
		{
			idx_latlon = j + i*Nlat;
			lat_cur = -90 + (180*j)/(Nlat-1);

			readIgloo(directory_name[i], "/HiDensity", lat_cur, iglooDataHi[idx_latlon]);
			readIgloo(directory_name[i], "/LoDensity", lat_cur, iglooDataLo[idx_latlon]);
			setupNEOIgloo(directory_name[i], lat_cur, NEOvel, iglooDataHi[idx_latlon], iglooDataNEO[idx_latlon], NEO_mass_min, NEO_mass_max);

			if (idx_latlon % 10 == 0 && idx_latlon > 0)
			{
				// timing stuff
				// see: https://stackoverflow.com/questions/12231166/timing-algorithm-clock-vs-time-in-c
				auto t_read_latlon2 = std::chrono::high_resolution_clock::now();
				runtime = std::chrono::duration_cast<std::chrono::milliseconds>(t_read_latlon2-t_read_latlon1).count() / 1000. / 60.;

				percentFinished = idx_latlon / double(NlatlonRegions);

				//cout << "*****************************************************\n";
				cout << " Percent finished = " << 100.*percentFinished  << endl;
				//cout << "*****************************************************\n";
				cout << "  run time = " << runtime << " minutes |  time remaining = " << runtime * (1./percentFinished - 1.) << " minutes\n";
			}
			
		}
	}
	

	//cout << "initError = " << initError << endl;

	return 0;
}