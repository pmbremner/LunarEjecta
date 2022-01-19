#include "LunarEjecta_PrimaryImpactor.h"


#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>


using namespace std;

// Next, generate primary impactor samples
/// We need to construct a CDF for the 3 types, MEM_LO, MEM_HI, and NEO.
/// First, the net flux for each type is computed and a uniform random number can choose which type to pull from
/// Then, we take the CDF of that type and pull from using a new unifrom random number
/// This pull will define the speed, zenith, azimuth, and flux.
//// Then, we need to sample the density and mass (separately, they don't depend on each other)

void get_primary_samples(input* p,                               // need info on density distributions and mass distributions
						 vector<iglooSet*> primaryFluxes,        // the set of fluxes, MEM_HI, MEM_LO, and NEO
						 int i_iglooSet,                         // the igloo set index
						 vector<double>&  p_sample_azimuth,      // primary azimuth, at impact [rad]
						 vector<double>&  p_sample_zenith,       // primary zenith, at impact [rad] (nominally horizon angle in igloo files)
						 vector<double>&  p_sample_speed,        // primary speed, [km/s]
						 vector<double>&  p_sample_flux_weight,  // flux weight, [#/m^2/yr]
						 vector<double>&  p_sample_density,      // primary density [kg/m^3]
						 vector<double>&  p_sample_mass,         // primary mass [g]
						 int N_p_sample )                       // number of pulls in igloo-density-mass sets
{

	


	// generate CDF for the MEM_LO, MEM_HI, and NEO

}