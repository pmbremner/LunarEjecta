#include "read_igloo_files.h"
#include "read_parameters.h"

#include <stdlib.h>     /* atoi */
#include <string>
#include <chrono>
#include <cmath>

using namespace std;

inline double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

const bool skipread = 1;

int main(int argc, char const *argv[])
{
	int i, j, idx_latlon, idx_cur;
	double lat_cur, lon_cur, runtime, percentFinished;

	string param_fn = argv[1];
	int N_proc = atoi(argv[2]); // total number of processes
	int i_proc = atoi(argv[3]); // the current process

	cout << "Process number = " << i_proc << " / " << N_proc << endl << endl;
	if (i_proc >= N_proc)
	{
		cout << "ERROR: i_proc must be between 0 and " << N_proc - 1 << endl;
		return -1;
	}

	int Nlat, Nlon, NlatlonRegions;

	getParam(param_fn, "Nlat", Nlat, 0);
	getParam(param_fn, "Nlon", Nlon, 0);

	NlatlonRegions = Nlat * Nlon;

	int N_rpp = ceil(NlatlonRegions / double(N_proc)); // number of regions per process
	int idx_latlon_min = i_proc * N_rpp; // starts from 0
	int idx_latlon_max = fmin((i_proc + 1)*N_rpp - 1, NlatlonRegions - 1); // inclusive, up to NlatlonRegions - 1
	int N_i_rpp = idx_latlon_max - idx_latlon_min + 1;

	cout << "min/max idx out of total: " << idx_latlon_min << ", " << idx_latlon_max << ", " << NlatlonRegions << endl;
	cout << "lat lon regions responsible for = " << N_i_rpp << endl;

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

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// all fluxes in units of #/m^2/yr
	vector<vector<double>> iglooDataHi;
	vector<vector<double>> iglooDataLo;
	vector<vector<double>> iglooDataNEO;

	iglooDataHi.resize(N_i_rpp);
	iglooDataLo.resize(N_i_rpp);
	iglooDataNEO.resize(N_i_rpp);

	vector<double> NEOvel;
	readVelDist(NEOvel_fn, NEOvel);

	// read input igloo files
	auto t_read_latlon1 = std::chrono::high_resolution_clock::now();
	if(!skipread) {
		for (i = 0; i < Nlon; ++i)
		{
			for (j = 0; j < Nlat; ++j)
			{
				idx_latlon = j + i*Nlat;

				// if this is a region my process is supposed to do
				if (idx_latlon >= idx_latlon_min && idx_latlon <= idx_latlon_max)
				{
					idx_cur = idx_latlon - idx_latlon_min;
					lat_cur = -90 + (180*j)/(Nlat-1);

					readIgloo(directory_name[i], "/HiDensity", lat_cur, iglooDataHi[idx_cur]);
					readIgloo(directory_name[i], "/LoDensity", lat_cur, iglooDataLo[idx_cur]);
					setupNEOIgloo(directory_name[i], lat_cur, NEOvel, iglooDataHi[idx_cur], iglooDataNEO[idx_cur], NEO_mass_min, NEO_mass_max);

					if (idx_latlon % 5 == 0 && idx_latlon > 0)
					{
						// timing stuff
						// see: https://stackoverflow.com/questions/12231166/timing-algorithm-clock-vs-time-in-c
						auto t_read_latlon2 = std::chrono::high_resolution_clock::now();
						runtime = std::chrono::duration_cast<std::chrono::milliseconds>(t_read_latlon2-t_read_latlon1).count() / 1000. / 60.;

						percentFinished = idx_cur / double(N_i_rpp);

						//cout << "*****************************************************\n";
						cout << " Percent finished = " << 100.*percentFinished  << endl;
						//cout << "*****************************************************\n";
						cout << "  run time = " << runtime << " minutes |  time remaining = " << runtime * (1./percentFinished - 1.) << " minutes\n";
					}
				}	
			}
		}
	}
	auto t_read_latlon2 = std::chrono::high_resolution_clock::now();
	runtime = std::chrono::duration_cast<std::chrono::milliseconds>(t_read_latlon2-t_read_latlon1).count() / 1000. / 60.;
	cout << "  reading lat-lon data from file = " << runtime << " minutes\n";
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	auto t_bearingdistance1 = std::chrono::high_resolution_clock::now();


	// for surface area for each lat-lon region
	vector<double> SA_latlon;
	SA_latlon.resize(N_i_rpp);

	double d_lon = 2.*PI/double(Nlon), SAtot = 0.0;
	double d_lat = PI/double(Nlat-1), lat_min, lat_max;

	for (i = 0; i < Nlon; ++i)
	{
		lon_cur = i * d_lon; // center

		for (j = 0; j < Nlat; ++j)
		{
			idx_latlon = j + i*Nlat;

			// if this is a region my process is supposed to do
			if (idx_latlon >= idx_latlon_min && idx_latlon <= idx_latlon_max)
			{
				idx_cur = idx_latlon - idx_latlon_min;
				lat_cur = -PI/2. + j*d_lat; // center

				if (j == 0) // -90 degrees lat
				{
					lat_min = -PI/2.;
					lat_max = lat_min + d_lat/2.;
				}
				else if(j == Nlat-1) // 90 degrees lat
				{
					lat_max = PI/2.;
					lat_min = lat_max - d_lat/2.;
				}
				else
				{
					lat_min = lat_cur - d_lat/2.;
					lat_max = lat_cur + d_lat/2.;
				}

				// check for region that overlaps with equator
				if (fabs(lat_cur) < d_lat/2.)
				{
					SA_latlon[idx_cur] = d_lon * (fabs(sin(lat_max)) + fabs(sin(lat_min)));
				}
				else
				{
					SA_latlon[idx_cur] = d_lon * fabs(sin(lat_max) - sin(lat_min));
				}
				
				SAtot += SA_latlon[idx_cur];

				//cout << lat_cur*180./PI <<' ' << SA_latlon[idx_cur] << endl;
			}
		}
	}
	cout << "Total surface area = " << SAtot / (4.*PI) << " x 4pi r^2\n";

	// For each latlon region, compute the bearing-distance map (2D-histogram)
	int N_D_perRegion, N_azm_perRegion;

	// For number bins for the bearing-distance map
	getParam(param_fn, "N_D_perRegion", N_D_perRegion, 0);
	getParam(param_fn, "N_azm_perRegion", N_azm_perRegion, 0);



	// compute lat and lon points in ROI


	return 0;
}