#include "LunarEjecta_MainUtils.h"
#include "LunarEjecta_SecondaryEjecta.h"
#include "LunarEjecta_params.h"
#include "LunarEjecta_igloo.h"
#include "LunarEjecta_PrimaryImpactor.h"



// note, -march=native is to allow for vectorization, if possible

//  g++ -O2 -std=c++17 -march=native primary_and_secondaries_test.cpp ./main_code/LunarEjecta_MainUtils.cpp ./main_code/LunarEjecta_SecondaryEjecta.cpp ./main_code/LunarEjecta_params.cpp ./main_code/LunarEjecta_igloo.cpp ./main_code/LunarEjecta_PrimaryImpactor.cpp -IC:\Users\AMD-Y500\Documents\GitHub\LunarEjecta\MonteCarlo_SimV3\main_code -o ejecta.exe

using namespace std;



int main(int argc, char const *argv[])
{
	// Need to move to param file
	const double Rm = 1737.E3; // m

	double vmin = 0.0;
	double vmax = 3.;
	double a = atof(argv[4])/Rm + 1.; // altitude above lunar center 
	double h = atof(argv[5])/Rm;//120./Rm;
	double r = atof(argv[6])/Rm;//4.5/Rm;

	double dg = 0.05; // [rad]
	double dv = 0.05; // [rm]

	int N_azm_lat_lon  = 100;
	int N_zenith_speed = 100;
	int N_p_sample     = 10000;
	// end of move to param file

	vector<double> p_sample_azimuth, p_sample_zenith, p_sample_speed, p_sample_flux_weight, p_sample_density, p_sample_mass;
	vector<double> sample_latp, sample_lonp, sample_azimuth_0, sample_zenith_0, sample_speed_0;
	vector<double> sample_azimuth_f, sample_zenith_f, sample_speed_f, sample_weight;
	double lat_center, lon_center, dlat, dlon;

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
	iglooSet *MEM_hi_fluxes = read_igloo(params, HiDensMEM); // size p->N_loc
	iglooSet *MEM_lo_fluxes = read_igloo(params, LoDensMEM); // size p->N_loc
	iglooSet *NEO_fluxes;

	cout << "Total Surface Area for process = " << sumSA(params, MEM_hi_fluxes) / (4*PI) << " x 4pi r^2\n";

	//iglooSet *primaryFluxes[3] = {MEM_hi_fluxes, MEM_lo_fluxes, NEO_fluxes};
	vector<iglooSet*> primaryFluxes;
	primaryFluxes.push_back(MEM_hi_fluxes); // pushing in the order is very important!
	primaryFluxes.push_back(MEM_lo_fluxes);
	primaryFluxes.push_back(NEO_fluxes);

	//cout << primaryFluxes[HiDensMEM][0].filename << endl;

	////////////////////////////////////////////
	// read or generate the NEO igloo fluxes
	////////////////////////////////////////////
	if (params->readNEO_files)
		NEO_fluxes = read_igloo(params, NEO);
	else {
		NEO_fluxes = generate_NEO_igloo(params, MEM_lo_fluxes);
		if (params->saveNEO_files)
			save_igloo(params, NEO_fluxes, NEO);
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	// for each lat-lon location that the process is responsible for
	for (params->latlon_idx_proc = 0; params->latlon_idx_proc < params->N_loc; params->latlon_idx_proc++)
	{
		cout << "\n\n    Process #: " << params->i_proc << " | Location #: " << params->latlon_idx_proc+1 << '/' << params->N_loc << endl;
		params->latlon_idx_cur = params->latlon_idx_proc + params->latlon_idx_min;


		lat_center = (primaryFluxes[HiDensMEM][params->latlon_idx_proc].latmin + primaryFluxes[HiDensMEM][params->latlon_idx_proc].latmax) / 2.;
		dlat       = primaryFluxes[HiDensMEM][params->latlon_idx_proc].latmax - primaryFluxes[HiDensMEM][params->latlon_idx_proc].latmin;
		
		lon_center = primaryFluxes[HiDensMEM][params->latlon_idx_proc].lon;
		dlon       = primaryFluxes[HiDensMEM][params->latlon_idx_proc].dlon;


		get_samples_with_azm_lat_lon( lat_center,   // primary latitude center
		                              lon_center,   // primary longitude center
		                              dlat,  // primary latitude range
		                              dlon,  // primary longitude range
		                              15.01 / 180. * PI,   // satellite (asset) latitude center
		                              110.01 / 180. * PI,   // satellite (asset) longitude center
		                              a,      // satellite (asset) altitude [rm]
		                              h,      // satellite (asset) height [rm]
		                              r,      // satellite (asset) radius [rm]
		                              vmin,   // minimum ejecta speed [vesc]
		                              vmax,   // maximum ejecta speed [vesc]
		                              dg,     // maximum zenith grid width
		                              dv,     // maximum speed grid width
		                              sample_latp,       // primary latitude center
		                              sample_lonp,       // primary longitude center
		                              sample_azimuth_0,  // initial azimuth, zenith, and speed of secondary, at primary impact
		                              sample_zenith_0,
		                              sample_speed_0,
		                              sample_azimuth_f,  // final azimuth, zenith, and speed of secondary, at asset impact
		                              sample_zenith_f,
		                              sample_speed_f,
		                              sample_weight,
		                              N_azm_lat_lon,   // number of pulls in azimuth-lat-lon sets
		                              N_zenith_speed); // number of pulls in zenith-speed sets


		// Next, generate primary impactor samples
		/// We need to construct a CDF for the 3 types, MEM_LO, MEM_HI, and NEO.
		/// First, the net flux for each type is computed and a uniform random number can choose which type to pull from
		/// Then, we take the CDF of that type and pull from using a new unifrom random number
		/// This pull will define the speed, zenith, azimuth, and flux.
		//// Then, we need to sample the density and mass (separately, they don't depend on each other)

		get_primary_samples(params,                // need info on density distributions and mass distributions
							primaryFluxes,         // the set of fluxes, MEM_HI, MEM_LO, and NEO
						 	params->latlon_idx_proc,            // the igloo set index
							p_sample_azimuth,      // primary azimuth, at impact [rad]
							p_sample_zenith,       // primary zenith, at impact [rad] (nominally horizon angle in igloo files)
							p_sample_speed,        // primary speed, [km/s]
							p_sample_flux_weight,  // flux weight, [#/m^2/yr]
							p_sample_density,      // primary density [kg/m^3]
							p_sample_mass,         // primary mass [g]
							N_p_sample );          // number of pulls in igloo-density-mass sets
							

	}

	return 0;
}