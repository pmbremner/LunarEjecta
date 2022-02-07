#include "LunarEjecta_Environment.h"
#include "LunarEjecta_MainUtils.h"
#include "LunarEjecta_params.h"


using namespace std;

//compute crater radius [m]
// note, we will be converting everything to SI (m, kg, s)
double get_crater_radius(input* p,  
	                     double impactor_r,       // [m]
	                     double p_sample_zenith,  // [rad]
	                     double p_sample_speed,   // [km/s]
	                     double p_sample_density, // [kg/m^3]
	                     double p_sample_mass,    // [g]
	                     double p_sample_type,
	                     int&   regime_type)
{
	// regolith bulk density [kg/m^3]
	double reg_bulk_dens = p->regolith_dens * p->regolith_porosity;

	// gravity regime radius, unitless scale of (rho/m)^(1/3)
	double r_grav_norm = p->HH11_H1
	                   * pow(p_sample_density / reg_bulk_dens, (2. + p->HH11_mu - 6.*p->HH11_nu) / (3.*(2. + p->HH11_mu)))
	                   * pow(p->lunar_acceleration * impactor_r / sqr(p_sample_speed * cos(p_sample_zenith) * 1000. /* km/s -> m/s*/), -p->HH11_mu / (2. + p->HH11_mu));

	// crater length scale [m]
	double length_scale = pow(p_sample_density / (p_sample_mass / 1000. /* g -> kg */), -1./3.);

	// Is there a zero tensile strength, i.e., only graviry regime
	if (p->regolith_tensile_strength == 0)
	{
		regime_type = 1;
		return r_grav_norm * length_scale;
	}
	else // need to compare strength and gravity regime crater sizes, return the smaller
	{
		double r_stng_norm = p->HH11_H2
		                   * pow(p_sample_density / reg_bulk_dens, (1. - 3.*p->HH11_nu) / 3.)
		                   * pow(p->regolith_tensile_strength / (p_sample_density * sqr(p_sample_speed * cos(p_sample_zenith) * 1000. /* km/s -> m/s*/)), -p->HH11_mu/2.);
		
		if (r_stng_norm < r_grav_norm)
		{ // in strength regime
			regime_type = 0;
			return r_stng_norm * length_scale;
		}
		else // in gravity regime
		{
			regime_type = 1;
			return r_grav_norm * length_scale;
		}
	}
}



// compute maximum speed ejected by the primary impact [m/s]
double get_crater_max_speed(input*  p,
    		                double  impactor_r,       // [m]
    		                double  p_sample_zenith,  // [rad]
	    		            double  p_sample_speed,   // [km/s]
	    		            double  p_sample_density, // [kg/m^3]
	    		            double  p_sample_mass,    // [g]
	    		            double  p_sample_type,
	    		            double  crater_r,         // [m]
	    		            int     regime_type,
	    		            double& m_crater_max)     // [kg]
{
	// regolith bulk density [kg/m^3]
	double reg_bulk_dens = p->regolith_dens * p->regolith_porosity;

	m_crater_max = (p_sample_mass / 1000. /* g -> kg */) * 3. * p->HH11_k / (4.*PI)
	             * (p_sample_density / reg_bulk_dens)
	             * (pow(((regime_type == 0 ? p->HH11_n2s : p->HH11_n2g) * crater_r) / (impactor_r), 3) - pow(p->HH11_n1, 3));

	return p_sample_speed * cos(p_sample_zenith) * 1000. /* km/s -> m/s*/ * p->HH11_C1
		 * pow(p->HH11_n1 * pow(p_sample_density / reg_bulk_dens, p->HH11_nu), -1./p->HH11_mu)
		 * pow(1. - (p->HH11_n1 * impactor_r) / ((regime_type == 0 ? p->HH11_n2s : p->HH11_n2g) * crater_r), p->HH11_p);
}



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
					        vector<double>& ejecta_env_size,     // m, diameter <--- need to check if it's dia or radius *********
					        vector<double>& ejecta_env_flux      // #/m^2/yr (> size_i)
	                        )
{
	ofstream file;
	file.open("crater_stats.txt");

	int i_p, i_s, j_s_size, crater_regime;
	double crater_r, impactor_r; // [m]
	double v_crater_max; // maximum speed the crater ejects [m/s]
	double m_crater_max, ejected_mass; // maximum mass the crater ejects [kg]

	// init ejecta vectors (ejecta_env_flux is 4d flattened to 1d, {speed, zenith, azimuth, particle size: #/m^2/yr > particle size i})
	ejecta_env_speed.clear();
    ejecta_env_zenith.clear();
    ejecta_env_azimuth.clear();
    ejecta_env_size.clear();
    ejecta_env_flux.clear();


	ejecta_env_speed.resize(params->N_env_v, 0.);
    ejecta_env_zenith.resize(params->N_env_zen, 0.);
    ejecta_env_azimuth.resize(params->N_env_azm, 0.);
    ejecta_env_size.resize(params->N_env_size, 0.);
    ejecta_env_flux.resize(params->N_env_flux, 0.);


    for (i_p = 0; i_p < N_p_sample; ++i_p) // for a given primary impact
    {
    	// compute impactor radius [m]
    	impactor_r = pow(3.*(p_sample_mass[i_p]/1000. /* g- > kg */) / (4.*PI * p_sample_density[i_p]), 1./3.);

    	//compute crater radius [m], and crater_regime (0 = strength, 1 = grav)
    	crater_r = get_crater_radius(params,
    		                         impactor_r,
    		                         p_sample_zenith[i_p],
    		                         p_sample_speed[i_p],
    		                         p_sample_density[i_p],
    		                         p_sample_mass[i_p],
    		                         p_sample_type[i_p],
    		                         crater_regime);

    	// compute maximum speed ejected by the primary impact [m/s], and max mass [kg]
    	v_crater_max = get_crater_max_speed(params,
    		                                impactor_r,
    		                                p_sample_zenith[i_p],
	    		                            p_sample_speed[i_p],
	    		                            p_sample_density[i_p],
	    		                            p_sample_mass[i_p],
	    		                            p_sample_type[i_p],
	    		                            crater_r,
	    		                            crater_regime,
	    		                            m_crater_max);
   
    	//cout << "impactor radius & mass: crater radius, max speed & mass: " << impactor_r << " m | " << p_sample_mass[i_p] / 1000. << " kg : " << crater_r/impactor_r << " r | " << v_crater_max << " m/s | " << m_crater_max/(p_sample_mass[i_p] / 1000.) << " yield\n";
    	if (m_crater_max > 0) // ejecta is created by impact (equivalent to checking if the max speed is >= 0)
    	{
    		file << impactor_r << ' ' << p_sample_mass[i_p] / 1000. << " " << crater_r/impactor_r << " " << v_crater_max << " " << m_crater_max/(p_sample_mass[i_p] / 1000.) << ' ' << p_sample_type[i_p] << endl;
    	


			// for each secondary created by the i_p'th primary impact, compute fraction of ejected mass and bin it
    		for (i_s = 0; i_s < N_s_sample; ++i_s)
    		{

    			if (ejecta_env_speed[i_s] < v_crater_max)
    			{
    				// compute mass ejected (total for all sizes)
	    			ejected_mass = get_ejected_mass(params,
	    				                            impactor_r,
	    				                            crater_r,
	    				                            ejecta_env_speed[i_s]);

	    			// bin ejecta, looping over cumulative sizes
	    			for (j_s_size = 0; j_s_size < params->N_env_size; ++j_s_size)
	    			{
	    				// compute cumulative # of particles > j_s_size-th size
	    				// if ejecta size is greater than impactor size, reject it (no ejecta)
	    				if ( )
	    				{
	    					/* code */
	    				}

	    			}


    			} // no else, if the ejecta speed is greater than the max ejected speed, there is no ejecta
    		
    		}

    	}
    	else // no ejecta created by impact
    	{

    	}


    }
    file.close();
}