#include "LunarEjecta_PrimaryImpactor.h"
#include "LunarEjecta_MainUtils.h"
#include "LunarEjecta_params.h"



#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>


using namespace std;


void get_primary_cdf(iglooSet* primaryFlux, vector<double>& cdf_i)
{
	int row_i, col_j;

	// the cdf array will be 1 larger than the pdf, where the first element will always be 0
	cdf_i.clear();
	cdf_i.resize(primaryFlux->N_rows * primaryFlux->N_cols + 1, 0.);

	for (row_i = 0; row_i < primaryFlux->N_rows; ++row_i)
	{
		for (col_j = 0; col_j < primaryFlux->N_cols; ++col_j)
		{
			cdf_i[1 + col_j + primaryFlux->N_cols * row_i]
				= cdf_i[col_j + primaryFlux->N_cols * row_i]
				+ (primaryFlux->iglooData[ig_idx(row_i, 9 + col_j, primaryFlux->N_cols)]) / primaryFlux->netFlux;

			//cout << 1 + col_j + primaryFlux->N_cols * row_i << ' ' << cdf_i[1 + col_j + primaryFlux->N_cols * row_i] << endl;
		}
	}
	cdf_i[primaryFlux->N_rows * primaryFlux->N_cols] = 1.; // force last element to be exactly 1
}

int get_type_sample(mt19937_64& rng_64, vector<long long int>& cdf_type)
{
	const long long int zero = 0;
	long long int ui = uniformInt(rng_64, zero, cdf_type[2]);

	//cout << ui << ' ' << cdf_type[0] << ' ' << cdf_type[1] << ' ' << cdf_type[2] << ' ' << rng_64.max() << endl;

	if (ui <= cdf_type[0])
		return HiDensMEM; // 0
	else if (ui <= cdf_type[1])
		return LoDensMEM; // 1
	else // (ui <= cdf_type[2])
		return NEO; // 2
}


void get_primary_igloo_sample_i(mt19937&        rng,
	                            vector<double>& p_cdf,
	                            iglooSet*       primaryFluxes,
	                            double& p_sample_azimuth,
	                            double& p_sample_zenith,
	                            double& p_sample_speed,
	                            double& p_sample_flux_weight)
{
	do { // do-while the sample is coming from below the horizon
		// pull sample from uniform distribution
		double u = uniform(rng, 0., 1.);

		// find index (iterator in this case) of the corresponding location in the cdf
		vector<double>::iterator idx_iter;
		int idx, row, col;

		// Find the index such that cdf(idx-1) <= u <= cdf(idx)
		// If u = 0, use upper_bound, otherwise use lower_bound (both binary search algorithms, logN time)
		//  effectively, this inverts the cdf
		// This guarantees that *(idx_iter-1) <= u <= *(idx_iter) for all values of u in [0,1]
		idx_iter = (u == 0. ? upper_bound(p_cdf.begin(), p_cdf.end(), u) : lower_bound(p_cdf.begin(), p_cdf.end(), u));

		idx = idx_iter - p_cdf.begin();

		// https://stackoverflow.com/questions/7070346/c-best-way-to-get-integer-division-and-remainder
		row = (int)idx / primaryFluxes->N_cols;
		col = idx % primaryFluxes->N_cols;

		// update index to take into account the info columns (there are 9 of them)
		idx = ig_idx(row, col + 9, primaryFluxes->N_cols);

		// use index to get the igloo data
		double azm_min, azm_max, zenith_min, zenith_max, speed_min, speed_max;

		/// get the azimuth sample, convert to rads
		azm_min = primaryFluxes->iglooData[ig_idx(row, 5, primaryFluxes->N_cols)] * PI / 180.;
		azm_max = primaryFluxes->iglooData[ig_idx(row, 6, primaryFluxes->N_cols)] * PI / 180.;

		p_sample_azimuth = uniform(rng, azm_min, azm_max);

		/// get the zenith sample, convert from horizon angle to zenith angle in rads
		zenith_min = (90. - primaryFluxes->iglooData[ig_idx(row, 3, primaryFluxes->N_cols)]) * PI / 180.;
		zenith_max = (90. - primaryFluxes->iglooData[ig_idx(row, 4, primaryFluxes->N_cols)]) * PI / 180.;

		p_sample_zenith = uniform(rng, zenith_min, zenith_max);
		

		// get the speed sample
		//cout << primaryFluxes->speedEdge[0] << ' ' << primaryFluxes->speedEdge[primaryFluxes->N_cols] << endl;
		speed_min = primaryFluxes->speedEdge[col];
		speed_max = primaryFluxes->speedEdge[col + 1];

		p_sample_speed = uniform(rng, speed_min, speed_max);

		// get the sample weight
		p_sample_flux_weight = primaryFluxes->iglooData[ig_idx(row, col + 9, primaryFluxes->N_cols)];
	} while (p_sample_zenith >= PI/2.); // we don't want primaries hitting from below the horizon x.x

}


// Next, generate primary impactor samples
/// We need to construct a CDF for the 3 types, MEM_LO, MEM_HI, and NEO.
/// First, the net flux for each type is computed and a uniform random number can choose which type to pull from
/// Then, we take the CDF of that type and pull from using a new unifrom random number
/// This pull will define the speed, zenith, azimuth, and flux.
//// Then, we need to sample the density and mass (separately, they don't depend on each other)
// Note:
// iglooSet *primaryFluxes[3] = {HiDensMEM, LoDensMEM, NEO};
void get_primary_samples(input* params,                          // need info on density distributions and mass distributions
						 vector<iglooSet*> primaryFluxes,        // the set of fluxes, MEM_HI, MEM_LO, and NEO
						 int i_iglooSet,                         // the igloo set index
						 vector<double>&  p_sample_azimuth,      // primary azimuth, at impact [rad]
						 vector<double>&  p_sample_zenith,       // primary zenith, at impact [rad] (nominally horizon angle in igloo files)
						 vector<double>&  p_sample_speed,        // primary speed, [km/s]
						 vector<double>&  p_sample_flux_weight,  // flux weight, [#/m^2/yr]
						 vector<double>&  p_sample_density,      // primary density [kg/m^3]
						 vector<double>&  p_sample_mass,         // primary mass [g]
						 vector<int>&  p_sample_type,         // (MEM_hi_fluxes, MEM_lo_fluxes, NEO_fluxes)
						 int N_p_sample )                        // number of pulls in igloo-density-mass sets
{
	ofstream primary_sample_file;
	primary_sample_file.open("primary_samples.txt");

	int i, idx_sample, p_type;
	// init all output arrays
	p_sample_azimuth.clear();
	p_sample_zenith.clear();
	p_sample_speed.clear();
	p_sample_flux_weight.clear();
	p_sample_density.clear();
	p_sample_mass.clear();
	p_sample_type.clear();

	p_sample_azimuth.resize(N_p_sample, 0.);
	p_sample_zenith.resize(N_p_sample, 0.);
	p_sample_speed.resize(N_p_sample, 0.);
	p_sample_flux_weight.resize(N_p_sample, 0.);
	p_sample_density.resize(N_p_sample, 0.);
	p_sample_mass.resize(N_p_sample, 0.);
	p_sample_type.resize(N_p_sample, 0.);

	// generate CDF for the MEM_HI, MEM_LO, and NEO
	/// p_cdf is the cdf of primaries, indexed p_cdf[PRIMARY_TYPE][ID]
	/// where PRIMARY_TYPE = {HiDensMEM, LoDensMEM, NEO, secEjecta};, and ID follows the ID of the igloo column
	vector<vector<double>> p_cdf;
	vector<double> cdf_i;

	for (i = 0; i < 3; ++i)
	{// i = {HiDensMEM, LoDensMEM, NEO};

		//cout << primaryFluxes[i]->filename << endl;

		get_primary_cdf(primaryFluxes[i], cdf_i);

		p_cdf.push_back(cdf_i);
	}


	// pull samples, first figure out which population to pull from {HiDensMEM, LoDensMEM, NEO, secEjecta}
	// then, take the sample and pack the sample parameters
	/// Sampling the cdf will tell us the azm, zen, speed, and flux
	/// separately, we need to sample the density distribution and mass distribution (they are independent of the igloo dist)

	// set up random generators
	mt19937_64 rng_64; // for picking type, need 64 bits for higher resolution
	mt19937 rng; // for picking igloo sample
	unsigned seed_0, seed_1; 

	seed_0 = random_device{}() * chrono::system_clock::now().time_since_epoch().count();
	seed_1 = random_device{}() * chrono::system_clock::now().time_since_epoch().count() * PI;
	rng_64.seed(seed_0);
	rng.seed(seed_1);

	// figure out relative ratios between flux types
	vector<long long int> cdf_type(3, 0); // not normalized, will use uniform random int from 0 to cdf_type[2]

	double flux_min = min(min(primaryFluxes[HiDensMEM]->netFlux, primaryFluxes[LoDensMEM]->netFlux), primaryFluxes[NEO]->netFlux);

	// -Wno-narrowing to get rid of warning, we are assuming these fit into the long long int correctly
	for (i = 0; i < 3; ++i){
		if (i == 0)
			cdf_type[i] = primaryFluxes[i]->netFlux / flux_min;
		else
			cdf_type[i] = cdf_type[i-1] + primaryFluxes[i]->netFlux / flux_min;

		//cout << cdf_type[i] << ' ';
	}
	//cout << endl;

	for (idx_sample = 0; idx_sample < N_p_sample; ++idx_sample)
	{
		// will be either 0, 1, or 2 for HiDensMEM, LoDensMEM, or NEO
		// using the relative weights provided in cdf_type
		p_sample_type[idx_sample] = get_type_sample(rng_64, cdf_type); 
		p_type = p_sample_type[idx_sample]; // temp var w/ shorter name


		get_primary_igloo_sample_i(rng,
			                       p_cdf[p_type],
			                       primaryFluxes[p_type],
			                       p_sample_azimuth[idx_sample],
			                       p_sample_zenith[idx_sample],
			                       p_sample_speed[idx_sample],
			                       p_sample_flux_weight[idx_sample]);

		// normalize the weight by the number of samples
		p_sample_flux_weight[idx_sample] /= double(N_p_sample);

		// get the density sample
		get_primary_density_sample_i(rng, params, p_type, p_sample_density[idx_sample]);

		// get the mass sample
		get_primary_mass_sample_i(rng, params, p_type, p_sample_mass[idx_sample]);



		primary_sample_file << p_type << ' '
		                    << p_sample_azimuth[idx_sample] * 180./PI << ' '
		                    << p_sample_zenith[idx_sample] * 180./PI << ' '
		                    << p_sample_speed[idx_sample] << ' '
		                    << p_sample_flux_weight[idx_sample] << ' '
		                    << p_sample_density[idx_sample] << ' '
		                    << p_sample_mass[idx_sample] << endl;

	}
	primary_sample_file.close();
}