#ifndef LUNAREJECTA_PRIMARYIMPACTOR_H
#define LUNAREJECTA_PRIMARYIMPACTOR_H

#include "LunarEjecta_params.h"
#include "LunarEjecta_igloo.h"



using namespace std;


void get_primary_samples(input* p,                               // need info on density distributions and mass distributions
						 vector<iglooSet*> primaryFluxes,        // the set of fluxes, MEM_HI, MEM_LO, and NEO
						 int i_iglooSet,                         // the igloo set index
						 vector<double>&  p_sample_azimuth,      // primary azimuth, at impact [rad]
						 vector<double>&  p_sample_zenith,       // primary zenith, at impact [rad] (nominally horizon angle in igloo files)
						 vector<double>&  p_sample_speed,        // primary speed, [km/s]
						 vector<double>&  p_sample_flux_weight,  // flux weight, [#/m^2/yr]
						 vector<double>&  p_sample_density,      // primary density [kg/m^3]
						 vector<double>&  p_sample_mass,         // primary mass [g]
						 int N_p_sample );                       // number of pulls in igloo-density-mass sets


#endif 