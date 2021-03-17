#ifndef LUNAREJECTA_SCALINGLAW_H
#define LUNAREJECTA_SCALINGLAW_H

#include "LunarEjecta_params.h"

using namespace std;

struct scalingLaw
{
	double m_factor_MEM; // kg, for HH11
	double m_factor_NEO; // kg, for HH11

	double rho_factor_MEMhi; // unitless, for HH11
	double rho_factor_MEMlo; // unitless, for HH11
	double rho_factor_NEO;   // unitless, for HH11

	// main factor on the number flux and velocity ratio terms
	double main_factor_MEMhi; // units of mass, kg
	double main_factor_MEMlo; // units of mass, kg
	double main_factor_NEO;   // units of mass, kg
};


double HH_Grun(double m);
double HH_MDGrun(double m);
double H_compHH11GrunMassFactor(input* p);
double H_compHH11NEOMassFactor(input* p);
double HH_HH11lognormal(double mu, double sigma, double nu, double regolith_dens);
double HH_HH11dens(double projectile_dense, double target_dense, double nu);
double H_compHH11DensityFactor(input* p, int fluxType);

scalingLaw* compute_constants_and_normalization(input* params);



#endif 