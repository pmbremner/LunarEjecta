#include "LunarEjecta_Environment.h"
#include "LunarEjecta_MainUtils.h"
#include "LunarEjecta_params.h"


using namespace std;

//compute crater radius [m]
// note, we will be converting everything to SI (m, kg, s)
double get_crater_radius(input* p,  
	                     double impactor_r,       // [m]
	                     double p_sample_zenith,  // [rad]
	                     double p_norm_speed,     // [m/s]
	                     double p_sample_density, // [kg/m^3]
	                     double p_sample_mass,    // [g]
	                     double p_sample_type,
	                     int&   regime_type,
	                     double reg_bulk_dens)    // [kg/m^3]
{
	// gravity regime radius, unitless scale of (rho/m)^(1/3)
	double r_grav_norm = p->HH11_H1
	                   * pow(p_sample_density / reg_bulk_dens, (2. + p->HH11_mu - 6.*p->HH11_nu) / (3.*(2. + p->HH11_mu)))
	                   * pow(p->lunar_acceleration * impactor_r / sqr(p_norm_speed), -p->HH11_mu / (2. + p->HH11_mu));

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
		                   * pow(p->regolith_tensile_strength / (p_sample_density * sqr(p_norm_speed)), -p->HH11_mu/2.);
		
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
	    		            double  p_norm_speed,     // [m/s]
	    		            double  p_sample_density, // [kg/m^3]
	    		            double  p_sample_mass,    // [g]
	    		            double  p_sample_type,
	    		            double  crater_r,         // [m]
	    		            int     regime_type,
	    		            double  reg_bulk_dens,    // [kg/m^3]
	    		            double& m_crater_max)     // [kg]
{
	m_crater_max = (p_sample_mass / 1000. /* g -> kg */) * 3. * p->HH11_k / (4.*PI)
	             * (p_sample_density / reg_bulk_dens)
	             * (pow(((regime_type == 0 ? p->HH11_n2s : p->HH11_n2g) * crater_r) / (impactor_r), 3) - pow(p->HH11_n1, 3));

	return p_norm_speed * p->HH11_C1
		 * pow(p->HH11_n1 * pow(p_sample_density / reg_bulk_dens, p->HH11_nu), -1./p->HH11_mu)
		 * pow(1. - (p->HH11_n1 * impactor_r) / ((regime_type == 0 ? p->HH11_n2s : p->HH11_n2g) * crater_r), p->HH11_p);
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// Housen Holsapple 2011, eq (14)
double H11_speed(double c1, double a, double rho, double del, double nu, double mu, double n2, double R, double p, double U, double x)
{
	return U * c1 * pow(x/a * pow(rho/del, nu), -1./mu) * pow(1. - x/(n2*R), p);
}

double H11_speed_v(double x, vector<double>& v)
{
	return H11_speed(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], x);
}


double cum_ejected_mass(input* p,
	                    int crater_regime,
	                    double reg_bulk_dens, // kg/m^3
	                    double impactor_r,   // m
	                    double crater_r,     // m
	                    double p_density,    // kg/m^3
	                    double p_norm_v,     // m/s
	                    double p_m,          // kg
	                    double v_i)          // m/s
{
	double n2 = (crater_regime == 0 ? p->HH11_n2s : p->HH11_n2g);

	vector<double> vars{p->HH11_C1,
	                    impactor_r, // m
	                    p_density,  // kg/m^3
	                    reg_bulk_dens, // kg/m^3
	                    p->HH11_nu,
	                    p->HH11_mu,
	                    n2,
	                    crater_r,
	                    p->HH11_p,
	                    p_norm_v}; 

	double x = findX(v_i, H11_speed_v, p->HH11_n1 * impactor_r, n2 * crater_r, vars);

	// Housen Holsapple 2011, eq (19)
	return p_m * 3. * p->HH11_k / (4.*PI)
	     * (p_density/reg_bulk_dens)
	     * (pow(x / impactor_r, 3) - pow(p->HH11_n1, 3));
}


// for a given impactor and ejecta speed, how much mass is ejected?
double get_ejected_mass(input*  p,
						int crater_regime, // 0 = strength, 1 = gravity
						double reg_bulk_dens,
	                    double impactor_r, // m
	                    double crater_r,   // m
	                    double p_sample_density, // kg/m^3
	                    double p_norm_speed,   // m/s (the normal component)
	                    double p_m, // kg
	                    double v_f, // m/s
	                    vector<double>& ejecta_env_speed) // m/s
{
	// find bounding edges around the sample speed v_f
	// Find the index such that ejecta_env_speed(idx-1) <= v_f <= ejecta_env_speed(idx)
	// If v_f = 0, use upper_bound, otherwise use lower_bound (both binary search algorithms, logN time)
	int idx = get_idx(v_f, ejecta_env_speed);
	
	//cout << ejecta_env_speed[idx-1] << ' ' << ejecta_env_speed[idx] << ' ';

	// return the ejected mass in the given speed range
	return ( cum_ejected_mass(p, crater_regime, reg_bulk_dens, impactor_r, crater_r, p_sample_density, p_norm_speed, p_m, ejecta_env_speed[idx-1])
	       - cum_ejected_mass(p, crater_regime, reg_bulk_dens, impactor_r, crater_r, p_sample_density, p_norm_speed, p_m, ejecta_env_speed[idx]) )
	       / (ejecta_env_speed[idx] - ejecta_env_speed[idx-1]);
}

// returns # of particles in 1kg of regolith greater than size 'diameter'
// This is derived from Carrier 2003, assuming a log-normal distribution for the mass-weighted PDF
double number_weighted_regolith_CDF(input* p,
	                                double diameter, // m
	                                double density)  // kg/m^3
{
	return 6. / (PI * density * p->reg_size_sigma * sqrt(2.*PI) * 1.E-9 /* eq, mm^3 to m^3 */)
	     * erfc( (log(diameter * 1000. /* m to mm */) - p->reg_size_mu + 3.*sqr(p->reg_size_sigma)) / (sqrt(2.) * p->reg_size_sigma) )
	     * exp(-3.*p->reg_size_mu + 9.*sqr(p->reg_size_sigma)/2.);
}

// based on Gault & Wedekind 1978, Figure 18 downstream angles (in degrees, but converted to rad at the end)
double H_get_zenith_max(double impactor_zen)
{
	double zen_deg = impactor_zen * 180. / PI;

	return (20. + zen_deg * (1.5206 + zen_deg * (-0.036 + zen_deg * 0.0003))) / 180. * PI;
}


double H_ejecta_zenith_dist(double ejecta_zen, double impactor_zen)
{
	const double s0 = 15./180.*PI; // assuming a full spread of 30 degrees about the peak, rad
	double s, l_bound, r_bound;    // rad

	double mu = H_get_zenith_max(impactor_zen); // peak of distribution, in rad

	l_bound = mu - s0;
	r_bound = mu + s0;

	// find an s such that the normalization we expect still holds
	if (l_bound >= 0. && r_bound <= PI/2.) // can use s0 as is
	{
		s = s0;
	}
	else if (l_bound < 0) // very close to zenith, need to limit s
	{
		s = mu;
	}
	else // if r_bound > PI/2., very close to horizon, need to limit s
	{
		s = PI/2. - mu;
	}

	// update the left and right bounds with the new s
	l_bound = mu - s;
	r_bound = mu + s;

	// check if ejceta is in the distribution
	if (ejecta_zen > l_bound && ejecta_zen < r_bound)
	{
		return (1. + cos((ejecta_zen - mu) * PI / s)) / (2. * s);
	}
	else // outside of distribution
	{
		return 0.;
	}

}

double H_ejecta_azimuth_dist(double ejecta_azm, double impactor_azm, double impactor_zen)
{
	double azm_stream_frame = ejecta_azm - (impactor_azm + PI); // impactor_azm is like wind, it report the direction it's coming from, not going

	if (impactor_zen < PI/3.) // if impactor zenith is less than 60 degrees (PI/3)
	{
		return (1. + (3. * impactor_zen / (2.*PI - 3.*impactor_zen)) * cos(azm_stream_frame)) / (2.*PI);
	}
	else // impactor zenith greater than 60 degrees (PI/3)
	{
		return (1. + cos(azm_stream_frame)) / (2.*PI);
	}
}

// input all in units of rads
double get_fraction_of_ejecta(double ejecta_zen,
	                          double ejecta_azm,
	                          double impactor_zen,
	                          double impactor_azm)
{
	return H_ejecta_zenith_dist(ejecta_zen, impactor_zen) * H_ejecta_azimuth_dist(ejecta_azm, impactor_azm, impactor_zen);
}



////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void get_ejecta_environment(input*  params,
					        // secondary ejecta 
					        vector<double>&  sample_latp,       // primary latitude center
	                        vector<double>&  sample_lonp,       // primary longitude center
	                        vector<double>&  sample_azimuth_0,  // initial azimuth, zenith, and speed of secondary, at primary impact
	                        vector<double>&  sample_zenith_0,   // [rad]
	                        vector<double>&  sample_speed_0,    // [vesc]
	                        vector<double>&  sample_azimuth_f,  // final azimuth, zenith, and speed of secondary, at asset impact
	                        vector<double>&  sample_zenith_f,   // [rad]
	                        vector<double>&  sample_speed_f,    // [vesc]
	                        vector<double>&  sample_weight,
	                        int              N_s_sample, 

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
					        vector<double>& ejecta_env_speed,    // [m/s], all at asset
					        vector<double>& ejecta_env_zenith,   // [rad]
					        vector<double>& ejecta_env_azimuth,  // [rad]
					        vector<double>& ejecta_env_size,     // [m], diameter (comfirmed)
					        vector<double>& ejecta_env_flux      // [#/m^2/yr (> size_i)]
	                        )
{
	chrono::steady_clock::time_point gen_environment_time_0;
	chrono::steady_clock::time_point gen_environment_time_1;

	chrono::duration<double> gen_environment_time_elapse;

	double etime;

	gen_environment_time_0 = chrono::steady_clock::now();

	ofstream file;
	file.open("crater_stats.txt");

	///ofstream idx_file;
    ///idx_file.open("idx.txt");

	double vesc = params->lunar_escape_speed;

	int i_p, i_s, j_s_size, crater_regime;
	double crater_r, impactor_r; // [m]
	double v_crater_max, p_norm_speed; // maximum speed the crater ejects [m/s]
	double m_crater_max, ejected_mass, sum_mass; // maximum mass the crater ejects [kg]
	double ejecta_frac; // unitless, percent of total ejecta (if all ejecta is counted, equals 1, i.e. normalized to 1)

	// init ejecta vectors (ejecta_env_flux is 4d flattened to 1d, {speed, zenith, azimuth, particle size: #/m^2/yr > particle size i})
    ejecta_env_flux.clear();
    ejecta_env_flux.resize(params->N_env_flux, 0.);

    // logspace(ejecta_env_speed,
    // 	     log10(params->vel_min * vesc),
    // 	     log10(params->vel_max * vesc),
    // 	     params->N_env_v);

    linspace(ejecta_env_speed , params->vel_min * vesc, params->vel_max * vesc, params->N_env_v);

    linspace(ejecta_env_zenith , 0., PI   , params->N_env_zen);
    linspace(ejecta_env_azimuth, 0., 2.*PI, params->N_env_azm);

    logspace(ejecta_env_size,
    	     log10(params->reg_min_size),
    	     log10(params->reg_max_size),
    	     params->N_env_size);
    
    // indices
    vector<int> idx_v(N_s_sample, 0.);
    vector<int> idx_g(N_s_sample, 0.);
    vector<int> idx_b(N_s_sample, 0.);
    vector<long long int> idx_f(N_s_sample * params->N_env_size, 0.);

    // We save N_p_sample-1 x recomputing the indices by precomputing them
    // precompute the index into the flux array
    for (i_s = 0; i_s < N_s_sample; ++i_s)
    {
    	// precompute indices idx_v, idx_g, idx_b;
		idx_v[i_s] = get_idx(sample_speed_f[i_s] * vesc, ejecta_env_speed  );
		idx_g[i_s] = get_idx(sample_zenith_f[i_s]      , ejecta_env_zenith );
		idx_b[i_s] = get_idx(sample_azimuth_f[i_s]     , ejecta_env_azimuth);

		for (j_s_size = 0; j_s_size < params->N_env_size; ++j_s_size)
		{
			idx_f[i_s*params->N_env_size + j_s_size] = j_s_size + params->N_env_size * (idx_b[i_s] + params->N_env_azm * (idx_g[i_s] + params->N_env_zen * idx_v[i_s]) );
		}
		
    }

    // regolith bulk density [kg/m^3]
	double reg_bulk_dens = params->regolith_dens * (1. - params->regolith_porosity);

    for (i_p = 0; i_p < N_p_sample; ++i_p) // for a given primary impact
    {
    	// compute impactor radius [m]
    	impactor_r = pow(3.*(p_sample_mass[i_p]/1000. /* g- > kg */) / (4.*PI * p_sample_density[i_p]), 1./3.);

    	p_norm_speed = p_sample_speed[i_p] * 1000. * cos(p_sample_zenith[i_p]);// ;// ; // m/s

    	//compute crater radius [m], and crater_regime (0 = strength, 1 = grav)
    	crater_r = get_crater_radius(params,
    		                         impactor_r,
    		                         p_sample_zenith[i_p],
    		                         p_norm_speed, // m/s
    		                         p_sample_density[i_p],
    		                         p_sample_mass[i_p], // g
    		                         p_sample_type[i_p],
    		                         crater_regime,
    		                         reg_bulk_dens);

    	// compute maximum speed ejected by the primary impact [m/s], and total mass [kg]
    	v_crater_max = get_crater_max_speed(params,
    		                                impactor_r,
    		                                p_sample_zenith[i_p],
	    		                            p_norm_speed, // m/s
	    		                            p_sample_density[i_p],
	    		                            p_sample_mass[i_p], // g
	    		                            p_sample_type[i_p],
	    		                            crater_r,
	    		                            crater_regime,
	    		                            reg_bulk_dens,
	    		                            m_crater_max);

   
    	//cout << "impactor radius & mass: crater radius, max speed & mass: " << impactor_r << " m | " << p_sample_mass[i_p] / 1000. << " kg : " << crater_r/impactor_r << " r | " << v_crater_max << " m/s | " << m_crater_max/(p_sample_mass[i_p] / 1000.) << " yield\n";
    	if (m_crater_max > 0) // ejecta is created by impact (equivalent to checking if the max speed is >= 0)
    	{
			// for each secondary created by the i_p'th primary impact, compute fraction of ejected mass and bin it
			sum_mass = 0.;

    		for (i_s = 0; i_s < N_s_sample; ++i_s)
    		{
    			// the sample ejecta speed must be less than the maximum possible ejected speed and also within the binning bounds
    			if (sample_speed_0[i_s] * vesc < v_crater_max && sample_speed_f[i_s] >= params->vel_min && sample_speed_f[i_s] <= params->vel_max)
    			{
    				


    				// cout << ejecta_env_speed[idx_v-1] << ' ' << sample_speed_f[i_s] * vesc << ' ' << ejecta_env_speed[idx_v] << endl;
    				// cout << ejecta_env_zenith[idx_g-1] << ' ' << sample_zenith_f[i_s]<< ' ' << ejecta_env_zenith[idx_g] << endl;
    				// cout << ejecta_env_azimuth[idx_b-1] << ' ' << sample_azimuth_f[i_s]<< ' ' << ejecta_env_azimuth[idx_b] << endl;

    				// compute mass ejected (total for all sizes) [kg]
	    			ejected_mass = get_ejected_mass(params,
	    				                            crater_regime, // 0 = strength, 1 = gravity
	    				                            reg_bulk_dens,
	    				                            impactor_r, // m
	    				                            crater_r,   // m
	    				                            p_sample_density[i_p], // [kg/m^3]
	    				                            p_norm_speed, // m/s
	    				                            p_sample_mass[i_p] / 1000., // kg
	    				                            sample_speed_0[i_s] * vesc, // m/s
	    				                            ejecta_env_speed); // m/s
	    			//if(ejected_mass > 0)
	    			sum_mass += ejected_mass * sample_weight[i_s];
	    			//cout << ejected_mass << endl;

	    			// compute fraction of ejecta (dependent on both ejecta and impact parameters)
	    			// If the fraction is zero, then skip the loop over sizes
	    			ejecta_frac = get_fraction_of_ejecta(sample_zenith_f[i_s],
	    				                                 sample_azimuth_f[i_s],
	    				                                 p_sample_zenith[i_p],
	    				                                 p_sample_azimuth[i_p]);

	    			if (ejecta_frac > 0.)
	    			{
	    			
		    			// bin ejecta, looping over cumulative sizes
		    			for (j_s_size = 0; j_s_size < params->N_env_size; ++j_s_size)
		    			{
		    				// compute cumulative # of particles > j_s_size-th size
		    				// if ejecta size is greater than impactor size, reject it (no ejecta)
		    				///if (ejecta_env_size[j_s_size]/2. < impactor_r)
		    				///{
		    				
		    					// overall units = #-ejecta/m^2/yr
		    					ejecta_env_flux[idx_f[i_s*params->N_env_size + j_s_size]]
		    					                        += ejected_mass // kg-ejecta / #-impactor
		    					                          * sample_weight[i_s]  // fraction of total ejecta, due to impact
		    					                          * ejecta_frac // weighted fraction of ejecta by angular distribution (zenith * azmimuth)
		    					                          * p_sample_flux_weight[i_p] // #-impactor/m^2/yr
		    					                          * number_weighted_regolith_CDF(params, ejecta_env_size[j_s_size], params->regolith_dens) // #-ejecta/kg-ejecta, we care about the particle density here, not bulk density
		    					                          * primaryFluxes[HiDensMEM][params->latlon_idx_proc].SA //* sqr(params->lunar_radius); // (Rm)^2
		    					                          / (2.*PI*params->asset_radius * (params->asset_radius + params->asset_height)); // surface area of asset (cylinder) 1/(Rm)^2
		    				
		    					// cout << ejecta_env_flux[idx_f] << endl;
		    				///}

		    			} // END loop over ejecta particle size
		    		} // END if ejecta_frac > 0.

    			} // no else, if the ejecta speed is greater than the max ejected speed, there is no ejecta
    		
    		} // END loop over ejecta

    		file << impactor_r << ' ' << p_sample_mass[i_p] / 1000. << " " << crater_r/impactor_r << " " << v_crater_max << " " << m_crater_max/(p_sample_mass[i_p] / 1000.) << ' ' << p_sample_type[i_p] << ' ' << sum_mass * p_sample_flux_weight[i_p] << endl;
    	

    	} // END if there is a crater
    	// else // no ejecta created by impact
    	// {

    	// }
    	if (i_p%5 == 0){
    		gen_environment_time_1 = chrono::steady_clock::now();
			gen_environment_time_elapse = chrono::duration_cast<chrono::duration<double>>(gen_environment_time_1 - gen_environment_time_0);
			etime = gen_environment_time_elapse.count();

			cout << "Generating Environment... " << 100.*(i_p+1.)/double(N_p_sample) << "% finished | time remaining = " << etime * (double(N_p_sample)/(i_p+1.) - 1.)/60. << " min\r";
		}

    } // END loop over primary impactors

    ofstream ej_env_file;
    ej_env_file.open("ejecta_environment_flux.txt");

    cout << "\n\nPrinting ejecta environment...\n";

	ej_env_file << params->N_env_v << ' ' << params->N_env_zen << ' ' << params->N_env_azm << ' ' << params->N_env_size << endl;

    for (int idx_f_i = 0; idx_f_i < params->N_env_flux; ++idx_f_i){
    	ej_env_file << idx_f_i << ' ' << ejecta_env_flux[idx_f_i] << endl;
    }
    ej_env_file.close();
    ///idx_file.close();


    file.close();
}