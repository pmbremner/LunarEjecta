#include "LunarEjecta_MainUtils.h"
#include "LunarEjecta_params.h"
#include "LunarEjecta_igloo.h"
#include "LunarEjecta_SecondaryEjecta.h"
#include "LunarEjecta_PrimaryImpactor.h"
#include "LunarEjecta_Environment.h"



// note, -march=native is to allow for vectorization, if possible
// -Wno-narrowing is to ignore warning about convering a double to a long long int

//  g++ -O2 -std=c++17 -march=native -Wno-narrowing primary_and_secondaries_test.cpp ./main_code/LunarEjecta_MainUtils.cpp ./main_code/LunarEjecta_SecondaryEjecta.cpp ./main_code/LunarEjecta_params.cpp ./main_code/LunarEjecta_igloo.cpp ./main_code/LunarEjecta_PrimaryImpactor.cpp ./main_code/LunarEjecta_Environment.cpp -IC:\Users\AMD-Y500\Documents\GitHub\LunarEjecta\MonteCarlo_SimV3\main_code -o ejecta.exe
//  g++ -O2 -std=c++17 -march=native -Wno-narrowing primary_and_secondaries_test.cpp ./main_code/LunarEjecta_MainUtils.cpp ./main_code/LunarEjecta_SecondaryEjecta.cpp ./main_code/LunarEjecta_params.cpp ./main_code/LunarEjecta_igloo.cpp ./main_code/LunarEjecta_PrimaryImpactor.cpp ./main_code/LunarEjecta_Environment.cpp -IC:\Users\adestefa\Documents\GitHub\LunarEjecta\MonteCarlo_SimV3\main_code -o ejecta.exe

using namespace std;



int main(int argc, char const *argv[])
{
	// timers
	// see https://www.cplusplus.com/reference/chrono/steady_clock/
	chrono::steady_clock::time_point read_igloo_time_0;
	chrono::steady_clock::time_point read_igloo_time_1;

	chrono::steady_clock::time_point gen_secondaries_time_0;
	chrono::steady_clock::time_point gen_secondaries_time_1;

	chrono::steady_clock::time_point gen_primaries_time_0;
	chrono::steady_clock::time_point gen_primaries_time_1;

	chrono::steady_clock::time_point gen_environment_time_0;
	chrono::steady_clock::time_point gen_environment_time_1;

	vector<double> p_sample_azimuth, p_sample_zenith, p_sample_speed, p_sample_flux_weight, p_sample_density, p_sample_mass;
	vector<int> p_sample_type;
	vector<double> sample_latp, sample_lonp, sample_azimuth_0, sample_zenith_0, sample_speed_0;
	vector<double> sample_azimuth_f, sample_zenith_f, sample_speed_f, sample_weight;
	double lat_center, lon_center, dlat, dlon, d;//, N_all_weight;

	vector<double> ejecta_env_speed, ejecta_env_zenith, ejecta_env_azimuth, ejecta_env_size, ejecta_env_flux;

	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	string param_fn = argv[1];
	int N_proc = atoi(argv[2]); // total number of processes
	int i_proc = atoi(argv[3]); // the current process

	cout << "Process number = " << i_proc << " / " << N_proc << endl << endl;
	if (i_proc >= N_proc || i_proc < 0)
	{
		cout << "ERROR: i_proc must be between 0 and " << N_proc - 1 << endl;
		return -1;
	}
	
	// input contains:
	//  parameters read from file, including regolith info, resolution of secondary ejecta igloo file
	input *params = init_input(param_fn, N_proc, i_proc);


	////////////////////////////////////////////
	// HiDensMEM, LoDensMEM, NEO
	// Read the MEM igloo files
	////////////////////////////////////////////

	read_igloo_time_0 = chrono::steady_clock::now();

	iglooSet *MEM_hi_fluxes = read_igloo(params, HiDensMEM); // size p->N_loc
	iglooSet *MEM_lo_fluxes = read_igloo(params, LoDensMEM); // size p->N_loc
	iglooSet *NEO_fluxes;

	cout << "Total Surface Area for process = " << sumSA(params, MEM_hi_fluxes) / (4*PI) << " x 4pi r^2\n";

	//cout << primaryFluxes[HiDensMEM][0].filename << endl;

	////////////////////////////////////////////
	// read or generate the NEO igloo fluxes
	////////////////////////////////////////////
	if (params->readNEO_files)
		NEO_fluxes = read_igloo(params, NEO);
	else {
		NEO_fluxes = generate_NEO_igloo(params, MEM_lo_fluxes); // NESC wants LO, MEM_lo_fluxes
		//if (params->saveNEO_files)
		save_igloo(params, NEO_fluxes, NEO);
		NEO_fluxes = read_igloo(params, NEO);
	}

	//iglooSet *primaryFluxes[3] = {MEM_hi_fluxes, MEM_lo_fluxes, NEO_fluxes};
	vector<iglooSet*> primaryFluxes;
	primaryFluxes.push_back(MEM_hi_fluxes); // pushing in this order is very important!
	primaryFluxes.push_back(MEM_lo_fluxes);
	primaryFluxes.push_back(NEO_fluxes); // need to push AFTER init of the NEO fluxes

	read_igloo_time_1 = chrono::steady_clock::now();
	chrono::duration<double> read_igloo_time_elapse = chrono::duration_cast<chrono::duration<double>>(read_igloo_time_1 - read_igloo_time_0);
	cout << endl << "Total time to read igloo files = " << read_igloo_time_elapse.count() << "s\n\n";


	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	// for each lat-lon location that the process is responsible for
	for (params->latlon_idx_proc = 0; params->latlon_idx_proc < params->N_loc; params->latlon_idx_proc++)
	{
		cout << "\n\n    Process #: " << params->i_proc << " | Location #: " << params->latlon_idx_proc+1 << '/' << params->N_loc << endl;
		params->latlon_idx_cur = params->latlon_idx_proc + params->latlon_idx_min;
		cout << "    Current idx: " << params->latlon_idx_cur << endl;


		lat_center = (primaryFluxes[HiDensMEM][params->latlon_idx_proc].latmin + primaryFluxes[HiDensMEM][params->latlon_idx_proc].latmax) / 2.;
		dlat       = primaryFluxes[HiDensMEM][params->latlon_idx_proc].latmax - primaryFluxes[HiDensMEM][params->latlon_idx_proc].latmin;
		
		lon_center = primaryFluxes[HiDensMEM][params->latlon_idx_proc].lon;
		dlon       = primaryFluxes[HiDensMEM][params->latlon_idx_proc].dlon;


		d = lat_lon_dist(lat_center, lon_center, params->asset_lat, params->asset_lon, 0);

		// get size and arrays of lat-lon locations
		params->N_azm_lat_lon = get_latlon_size(0, params->latlon_idx_cur); // replace 0 with process number

		//N_all_weight = 1.;//pow(0.05 + d/4., -(3./2.*params->HH11_mu + 1));

		gen_secondaries_time_0 = chrono::steady_clock::now();

		get_samples_with_azm_lat_lon( params,
		                              lat_center,   // primary latitude center
		                              lon_center,   // primary longitude center
		                              dlat,  // primary latitude range
		                              dlon,  // primary longitude range
		                              params->asset_lat,   // satellite (asset) latitude center
		                              params->asset_lon,   // satellite (asset) longitude center
		                              params->asset_altitude,      // satellite (asset) altitude [rm]
		                              params->asset_height,      // satellite (asset) height [rm]
		                              params->asset_radius,      // satellite (asset) radius [rm]
		                              params->vel_min,   // minimum ejecta speed [vesc]
		                              params->vel_max,   // maximum ejecta speed [vesc]
		                              params->dg_max,     // maximum zenith grid width
		                              params->dv_max,     // maximum speed grid width
		                              sample_latp,       // primary latitude center
		                              sample_lonp,       // primary longitude center
		                              sample_azimuth_0,  // initial azimuth, zenith, and speed of secondary, at primary impact
		                              sample_zenith_0,
		                              sample_speed_0,    // [vesc]
		                              sample_azimuth_f,  // final azimuth, zenith, and speed of secondary, at asset impact
		                              sample_zenith_f,
		                              sample_speed_f,    // [vesc]
		                              sample_weight,
		                              params->N_azm_lat_lon,   // number of pulls in azimuth-lat-lon sets
		                              params->N_zenith_speed); // number of pulls in zenith-speed sets

		gen_secondaries_time_1 = chrono::steady_clock::now();
		chrono::duration<double> gen_secondaries_time_elapse = chrono::duration_cast<chrono::duration<double>>(gen_secondaries_time_1 - gen_secondaries_time_0);
		cout << endl << "Total time to generate secondaries = " << gen_secondaries_time_elapse.count() << "s\n\n";

		// Next, generate primary impactor samples
		/// We need to construct a CDF for the 3 types, MEM_LO, MEM_HI, and NEO.
		/// First, the net flux for each type is computed and a uniform random number can choose which type to pull from
		/// Then, we take the CDF of that type and pull from using a new unifrom random number
		/// This pull will define the speed, zenith, azimuth, and flux.
		//// Then, we need to sample the density and mass (separately, they don't depend on each other)

		gen_primaries_time_0 = chrono::steady_clock::now();

		get_primary_samples(params,                // need info on density distributions and mass distributions
							primaryFluxes,         // the set of fluxes, MEM_HI, MEM_LO, and NEO
						 	params->latlon_idx_proc,            // the igloo set index
							p_sample_azimuth,      // primary azimuth, at impact [rad]
							p_sample_zenith,       // primary zenith, at impact [rad] (nominally horizon angle in igloo files)
							p_sample_speed,        // primary speed, [km/s]
							p_sample_flux_weight,  // flux weight, [#-impactors/m^2/yr]
							p_sample_density,      // primary density [kg/m^3]
							p_sample_mass,         // primary mass [g]
							p_sample_type,         // (MEM_hi_fluxes, MEM_lo_fluxes, NEO_fluxes)
							params->N_primary_sample );          // number of pulls in igloo-density-mass sets

		gen_primaries_time_1 = chrono::steady_clock::now();				
		chrono::duration<double> gen_primaries_time_elapse = chrono::duration_cast<chrono::duration<double>>(gen_primaries_time_1 - gen_primaries_time_0);
		cout << endl << "Total time to generate primaries = " << gen_primaries_time_elapse.count() << "s\n\n";

		// Finally, we convolute the secondaries and the primaries to get the ejected mass
		gen_environment_time_0 = chrono::steady_clock::now();	

		get_ejecta_environment(params,
							   // secondary ejecta 
							   sample_latp,       // primary latitude center
                               sample_lonp,       // primary longitude center
                               sample_azimuth_0,  // initial azimuth, zenith, and speed of secondary, at primary impact
                               sample_zenith_0,
                               sample_speed_0,    // [vesc]
                               sample_azimuth_f,  // final azimuth, zenith, and speed of secondary, at asset impact
                               sample_zenith_f,
                               sample_speed_f,    // [vesc]
                               sample_weight,    
                               params->N_azm_lat_lon * params->N_zenith_speed, // N_s_sample

                               // primary impactors
                               primaryFluxes,         // the set of fluxes, MEM_HI, MEM_LO, and NEO
							   params->latlon_idx_proc,            // the igloo set index
							   p_sample_azimuth,      // primary azimuth, at impact [rad]
							   p_sample_zenith,       // primary zenith, at impact [rad] (nominally horizon angle in igloo files)
							   p_sample_speed,        // primary speed, [km/s]
							   p_sample_flux_weight,  // flux weight, [#/m^2/yr]
							   p_sample_density,      // primary density [kg/m^3]
							   p_sample_mass,         // primary mass [g]
							   p_sample_type,         // (MEM_hi_fluxes, MEM_lo_fluxes, NEO_fluxes)
							   params->N_primary_sample,

							   // ejecta environment, at asset, sizes of each dimension are in params
							   ejecta_env_speed,    // km/s, all at asset
							   ejecta_env_zenith,   // rad
							   ejecta_env_azimuth,  // rad
							   ejecta_env_size,     // m, diameter
							   ejecta_env_flux      // #-ejecta/yr (> size_i)
			                   );

		gen_environment_time_1 = chrono::steady_clock::now();
		chrono::duration<double> gen_environment_time_elapse = chrono::duration_cast<chrono::duration<double>>(gen_environment_time_1 - gen_environment_time_0);
		cout << endl << "Total time to generate environment = " << gen_environment_time_elapse.count() << "s\n\n";

	}

	return 0;
}