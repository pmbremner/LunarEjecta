#include "LunarEjecta_Environment.h"
#include "LunarEjecta_MainUtils.h"
#include "LunarEjecta_params.h"


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
	                        )
{
	int i_p i_s;
	double crater_r; // [m]
	double v_crater_max; // maximum speed the crater ejects [m/s]

	// init ejecta vectors (ejecta_env_flux is 4d flattened to 1d, {speed, zenith, azimuth, particle size: #/m^2/yr > particle size i})
	ejecta_env_speed.clear();
    ejecta_env_zenith.clear();
    ejecta_env_azimuth.clear();
    ejecta_env_size.clear();
    ejecta_env_flux.clear();


	ejecta_env_speed.resize(N_env_v, 0.);
    ejecta_env_zenith.resize(N_env_zen, 0.);
    ejecta_env_azimuth.resize(N_env_azm, 0.);
    ejecta_env_size.resize(N_env_size, 0.);
    ejecta_env_flux.resize(N_env_flux, 0.);


    for (i_p = 0; i_p < N_p_sample; ++i_p) // for a given primary impact
    {
    	//compute crater radius [m]
    	crater_r = get_crater_radius(params,
    		                         p_sample_zenith[i_p],
    		                         p_sample_speed[i_p],
    		                         p_sample_density[i_p],
    		                         p_sample_mass[i_p],
    		                         p_sample_type[i_p]);

    	// compute maximum speed ejected by the primary impact [m/s]
    	v_crater_max = get_crater_max_speed(params,
    		                                p_sample_zenith[i_p],
	    		                            p_sample_speed[i_p],
	    		                            p_sample_density[i_p],
	    		                            p_sample_mass[i_p],
	    		                            p_sample_type[i_p],
	    		                            crater_r);
   
    	if (v_crater_max <= 0) // no ejecta created by impact
    	{
    	
    	}
    	else // ejecta is created by impact
    	{
    		// for each secondary created by the i_p'th primary impact, compute fraction of ejected mass and bin it
    		for (i_s = 0; i_s < N_s_sample; ++i_s)
    		{
    			



    		}
    	}


    }

}