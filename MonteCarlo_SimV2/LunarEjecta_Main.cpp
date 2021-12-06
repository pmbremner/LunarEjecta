#include "LunarEjecta_params.h"
#include "LunarEjecta_igloo.h"
#include "LunarEjecta_battleshipMonteCarloV2.h"

#include "LunarEjecta_track.h"
#include "LunarEjecta_asset.h"

// #include "LunarEjecta_scalinglaw.h"
// #include "LunarEjecta_bearingdistmap.h"
// #include "LunarEjecta_secondaryejecta.h"

using namespace std;

 // g++ -O2 LunarEjecta_Main.cpp LunarEjecta_params.cpp LunarEjecta_igloo.cpp LunarEjecta_battleshipMonteCarloV2.1.cpp -o ejecta.exe
 // g++ -O2 -std=c++17 LunarEjecta_Main.cpp LunarEjecta_params.cpp LunarEjecta_igloo.cpp LunarEjecta_battleshipMonteCarloV2.1.cpp LunarEjecta_track.cpp LunarEjecta_asset.cpp -o ejecta.exe

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

	//// test
	// func_ID.push_back(0);
	// func_ID.push_back(0);

	ofstream phase_hit_file, phase_miss_file;
	phase_hit_file.open("phase_hit_file.txt");
	phase_miss_file.open("phase_miss_file.txt");

	// init asset
	asset lunarLander;
	string asset_fn = "asset_tallCylinder.txt";

	init_asset(lunarLander, asset_fn, params);


	// for each lat-lon location that the process is responsible for
	for (params->latlon_idx_proc = 0; params->latlon_idx_proc < params->N_loc; params->latlon_idx_proc++)
	{
		cout << "\n\n    Process #: " << params->i_proc << " | Location #: " << params->latlon_idx_proc+1 << '/' << params->N_loc << endl;
		params->latlon_idx_cur = params->latlon_idx_proc + params->latlon_idx_min;

		////////////////////////////////////////////
		// generate list of hit shots
		////////////////////////////////////////////
		int hit_count = 0, tot_tries = 0;
		bool hit = 0;
		double weight, sum = 0.;

		// Define the domain ranges
		vector<double> ph, ph_i, ph_f; // [speed (m/s), zenith (rad), azimuth (rad)]
		vector<double> dph; // [speed (m/s), zenith (rad), azimuth (rad)]

		// position of impact point in lat-lon region
		vector<double> loc_latlon; // [lat (rad), lon (rad)]
		vector<double> loc_cart, loc_f;   // [x (m), y (m), z (m)]

		//// need to take into account vmin and vmax***
		ph.push_back(0.5 * (params->vel_min + params->vel_max) ); // center of speed range
		ph.push_back(PI/4.);  // center of zenith range
		ph.push_back(PI);  // center of azimuth range

		dph.push_back(params->vel_max - params->vel_min);  // length of speed range
		dph.push_back(PI/2.);  // length of zenith range
		dph.push_back(2.*PI);  // length of azimuth range

		//// test
		// ph.push_back(0.);
		// ph.push_back(0.);
		// dph.push_back(20.);
		// dph.push_back(20.);


		radar_scanner scanner;

		initRadar(scanner,
			      params->N_max,
			      params->alpha_search,
			      params->lifetime_max,
			      params->lifetime_rate,
			      params->dx_rate,
			      ph,
			      dph,
			      func_ID);

		


		while (hit_count < params->N_hit && tot_tries < params->N_max)
		{

			// randomly pull lat-lon position in lat-lon region, in terms of rad and cartesian meters
			uniformLatLon(scanner.rng,
				          loc_latlon, // rad
				          loc_cart,   // m
				          params->lunar_radius, // m
				          primaryFluxes[HiDensMEM][params->latlon_idx_proc].latmin,
				          primaryFluxes[HiDensMEM][params->latlon_idx_proc].latmax,
				          primaryFluxes[HiDensMEM][params->latlon_idx_proc].lon
				           - primaryFluxes[HiDensMEM][params->latlon_idx_proc].dlon/2.,
				          primaryFluxes[HiDensMEM][params->latlon_idx_proc].lon
				           + primaryFluxes[HiDensMEM][params->latlon_idx_proc].dlon/2.);
			
			//cout << "latlon location = " << loc_latlon[0] * 180./PI << ' ' << loc_latlon[1] * 180./PI << endl;

			//cout << endl << endl;
			weight = getSampleScan(scanner, ph_i);

			// test hit function
			
			//hit = ((ph_i[1] >= sqr(ph_i[0]) && ph_i[1] <= sqrt(ph_i[0])) || (ph_i[1] - 4. >= sqr(ph_i[0]-7.) && ph_i[1] -3.75 <= sqrt(ph_i[0]-7.)) ? 1 : 0);
			//hit = (ph_i[1] >= sqr(ph_i[0]) && ph_i[1] <= sqrt(ph_i[0]) ? 1 : 0);
			//hit =  (sqr(ph_i[0] - 4.) + sqr(ph_i[1] - 2.5) <= 1. ? 1 : 0);
			////hit =  (sqr(ph_i[0] + 3.) + sqr(ph_i[1] - 2.) <= 1. ? 1 : 0);
			//hit = (fabs(ph_i[0] - 0.3) < 0.01 && fabs(ph_i[1] - 0.25) < 0.25 ? 1 : 0); 

			//cout << "ph_i = " << ph_i[0] << " , " << ph_i[1] << " | " << hit << endl;

			// check if the particle trajectory hits the asset or the moon/out of bounds (i.e., a miss)
			hit = runTraj_checkHit(loc_cart,
				             ph_i,
				             loc_f,
				             ph_f,
				             lunarLander,
				             params->lunar_radius,
				             params->lunar_escape_speed,
				             params->lunar_acceleration);


			tallyScan(scanner, hit);

			if (hit)
			{
				cout << "ph_i [(m/s), zen, azm (deg)] = " << ph_i[0] << " , " << ph_i[1]/PI*180. <<  " , " << ph_i[2]/PI*180. << " | ";
				//cout << loc_f[0] << ' ' << loc_f[1] << ' ' << loc_f[2]
				cout << " | lat =  " << (PI/2. - atan2(sqrt(sqr(loc_f[0]) + sqr(loc_f[1])), loc_f[2]))*180./PI <<  endl;
			}


			// print hit information to file to read for next section

			//cout << "weight = " << weight << endl;

			if (hit){
				hit_count++;
				sum += weight;

				for (int i = 0; i < ph.size(); ++i)
					phase_hit_file << scanner.ph_scan[i] << ' ';
				phase_hit_file << (*scanner.idx_scan).generation << endl;
			} else
			{
				for (int i = 0; i < ph.size(); ++i)
					phase_miss_file << scanner.ph_scan[i] << ' ';
				phase_miss_file << (*scanner.idx_scan).generation << endl;
			}


			tot_tries++;
		}

		sum /= double(tot_tries);
		cout << " Sum = " << sum << endl;

		printHitMissReport(scanner);
	}
	phase_hit_file.close();
	phase_miss_file.close();




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