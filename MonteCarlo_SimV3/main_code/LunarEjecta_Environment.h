#ifndef LUNAREJECTA_ENVIRONMENT_H
#define LUNAREJECTA_ENVIRONMENT_H

#include "LunarEjecta_params.h"
#include "LunarEjecta_igloo.h"

using namespace std;



void get_ejecta_environment(input*  params,
					        // secondary ejecta 
					        vector<double>&  sample_latp,       // primary latitude center
	                        vector<double>&  sample_lonp,       // primary longitude center
	                        vector<double>&  sample_azimuth_0,  // initial azimuth, zenith, and speed of secondary, at primary impact
	                        vector<double>&  sample_zenith_0,
	                        vector<double>&  sample_speed_0,
	                        vector<double>&  sample_azimuth_f,  // final azimuth, zenith, and speed of secondary, at asset impact
	                        vector<double>&  sample_zenith_f,
	                        vector<double>&  sample_speed_f,
	                        vector<double>&  sample_altitude_f, // [rm]
		                    vector<double>&  sample_distance_f, // [rm]
	                        vector<double>&  sample_weight,
	                        int              N_s_sample, // N_s_sample

	                        // primary impactors
	                        vector<iglooSet*> primaryFluxes,         // the set of fluxes, MEM_HI, MEM_LO, and NEO
					        int               i_iglooSet,            // the igloo set index
					        vector<double>&   p_sample_azimuth,      // primary azimuth, at impact [rad]
					        vector<double>&   p_sample_zenith,       // primary zenith, at impact [rad] (nominally horizon angle in igloo files)
					        vector<double>&   p_sample_speed,        // primary speed, [km/s]
					        vector<double>&   p_sample_flux_weight,  // flux weight, [#/m^2/yr]
					        vector<double>&   p_sample_density,      // primary density [kg/m^3]
					        vector<double>&   p_sample_mass,         // primary mass [g]
					        vector<int>&      p_sample_type,         // (MEM_hi_fluxes, MEM_lo_fluxes, NEO_fluxes)
					        int               N_p_sample,

					        // ejecta environment, at asset, sizes of each dimension are in params
					        vector<double>& ejecta_env_speed,    // km/s, all at asset
					        vector<double>& ejecta_env_zenith,   // rad
					        vector<double>& ejecta_env_azimuth,  // rad
					        vector<double>& ejecta_env_size,     // m, diameter
					        vector<double>& ejecta_env_flux      // #/m^2/yr (> size_i)
	                        );







#endif 