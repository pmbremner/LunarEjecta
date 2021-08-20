#include "LunarEjecta_params.h"
#include "LunarEjecta_igloo.h"
#include "LunarEjecta_battleshipMonteCarloV2.h"
// #include "LunarEjecta_scalinglaw.h"
// #include "LunarEjecta_bearingdistmap.h"
// #include "LunarEjecta_secondaryejecta.h"

using namespace std;

 // g++ -O2 LunarEjecta_Main.cpp LunarEjecta_params.cpp LunarEjecta_igloo.cpp LunarEjecta_battleshipMonteCarloV2.cpp -o ejecta.exe

int main(int argc, char const *argv[])
{

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
	primaryFluxes.push_back(MEM_hi_fluxes);
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

	vector<int> func_ID;
	func_ID.push_back(0); // speed, const
	func_ID.push_back(1); // zenith, sin
	func_ID.push_back(0); // azimuth, const

	// for each lat-lon location that the process is responsible for
	for (params->latlon_idx_proc = 0; params->latlon_idx_proc < params->N_loc; params->latlon_idx_proc++)
	{
		cout << "\n\n    Process #: " << params->i_proc << " | Location #: " << params->latlon_idx_proc+1 << '/' << params->N_loc << endl;
		params->latlon_idx_cur = params->latlon_idx_proc + params->latlon_idx_min;

		////////////////////////////////////////////
		// generate list of hit shots
		////////////////////////////////////////////
		int hit_count = 0;
		double weight;

		// Define the domain ranges
		vector<double> ph, ph_i; // [speed (m/s), zenith (rad), azimuth (rad)]
		vector<double> dph; // [speed (m/s), zenith (rad), azimuth (rad)]

		ph.push_back(0.5 * params->lunar_escape_speed); // center of speed range
		ph.push_back(0.);  // center of zenith range
		ph.push_back(0.);  // center of azimuth range

		dph.push_back(params->lunar_escape_speed);  // center of speed range
		dph.push_back(PI);  // center of zenith range
		dph.push_back(PI);  // center of azimuth range


		radar_scanner scanner;

		initRadar(scanner, params->alpha_search, params->lifetime_max, params->lifetime_rate, params->dx_rate, ph, dph, func_ID);

		while (hit_count < params->N_hit)
		{

			weight = getSampleScan(scanner, ph_i);

			hit_count++;
		}


	}




	// ////////////////////////////////////////////
	// // compute constants and normalization of the scaling laws
	// ////////////////////////////////////////////
	// scalingLaw *ejectaFactors = compute_constants_and_normalization(params);


	// // for each lat-lon location that the process is responsible for
	// for (params->latlon_idx_proc = 0; params->latlon_idx_proc < params->N_loc; params->latlon_idx_proc++)
	// {
	// 	cout << "\n\n    Process #: " << params->i_proc << " | Location #: " << params->latlon_idx_proc+1 << '/' << params->N_loc << endl;
	// 	params->latlon_idx_cur = params->latlon_idx_proc + params->latlon_idx_min;
	// 	////////////////////////////////////////////
	// 	// compute the bearing-distance map for the current location
	// 	////////////////////////////////////////////
	// 	hist3DSet *bearingDistMap_short = init_bearing_dist_map(params, ejectaShort);
	// 	hist3DSet *bearingDistMap_far   = init_bearing_dist_map(params, ejectaFar);

	// 	////////////////////////////////////////////
	// 	// compute the secondary fluxes for the current location
	// 	////////////////////////////////////////////
	// 	iglooSet *secondary_fluxes_short = compute_ejecta(params, primaryFluxes, ejectaFactors, bearingDistMap_short, ejectaShort);
	// 	iglooSet *secondary_fluxes_far   = compute_ejecta(params, primaryFluxes, ejectaFactors, bearingDistMap_far  , ejectaFar);

	// 	////////////////////////////////////////////
	// 	// save the secondary fluxes for the current location
	// 	////////////////////////////////////////////
	// 	save_igloo(params, secondary_fluxes_short, secEjecta);
	// 	save_igloo(params, secondary_fluxes_far, secEjecta);

	// 	delete bearingDistMap_short, bearingDistMap_far;
	// 	delete secondary_fluxes_short, secondary_fluxes_far;
	// }

	
	/*

	iglooSet *secondary_fluxes = compute_ejecta(params, primaryFluxes, azmDistMap, ejectaFactors);
	save_igloo(params, secondary_fluxes, "secondaryEjecta");

	*/

	return 0;
}