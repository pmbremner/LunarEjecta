#include "LunarEjecta_scalinglaw.h"
#include "LunarEjecta_params.h"

#include <cmath>

using namespace std;


double HH_Grun(double m) { // m = units of g, output in units of 1/(m^2*yr)
	return pow(2.2E3 * pow(m, 0.306) + 15., -4.38)
	     + 1.3E-9 * pow(m + 1.E11 * pow(m, 2.) + 1.E27 * pow(m, 4.), -0.36)
	     + 1.3E-16 * pow(m + 1.E6 * pow(m, 2.), -0.85); 
}

// fit made in SciDAVis, already normalized to 10^-6 g
// double HH_MDGrun(double m) {
// 	return 1. / (1.587E7 * pow(m, 2.3235-1.) + 7.056E4 * pow(m, 1.8189-1.));
// }

// derived in mathematica, includes the mass term (linear) from HH11
double HH_MDGrun(double m) {
	return m * (2948.62 * (pow(15. + 2200.*pow(m, 0.306), -5.38) * pow(m, -0.694))
		 + 1.105E-16 * (1. + 2.E6*m) * pow(m + 1.E6*m*m, -1.85)
		 + 4.68E-10 * (1. + 2.E11*m + 4.E27*pow(m, 3)) * pow(m + 1.E11*m*m + 1.E27*pow(m, 4) , -1.36)) / HH_Grun(1.E-6);
}


double H_compHH11GrunMassFactor(input* p) { // N ~100 is usually sufficient
	int i, N = 200.;
	double massMin = p->MEM_massMin; // g
	double massMax = p->MEM_massMax; // mass limit in MEM is 10 g
	double Ai, Bi, sum = 0.;
	vector<double> m;
	vector<double> f;
	f.resize(N);

	// prepare x and y vectors for power law fits
	logspace(m, massMin, massMax, N, 0, N);
	for (i = 0; i < N; ++i)
		f[i] = HH_MDGrun(m[i]) * (HH_Grun(massMin) / HH_Grun(1.E-6));

	// compute integral by breaking up in log spacings and approx by power law
	for (i = N-2; i >= 0; i--) // start at low end of mass to sum smaller #'s first (less floating-point error)
	{
		Bi = log10(f[i+1] / f[i]) / log10(m[i+1] / m[i]); 
		Ai = (f[i] + f[i+1]) / (pow(m[i], Bi) + pow(m[i+1], Bi));

		sum += Ai/(Bi+1.) * (pow(m[i+1], Bi+1) - pow(m[i], Bi+1));
		//cout << i << ' ' << Bi << endl;
	}
	return sum / 1000.; // units of kg	
}

// mass in grams
double HH_compHH11NEOIntegratedMass(double m)
{
	return 2.601E-10 * pow(m, 0.1);
}

double H_compHH11NEOMassFactor(input* p)
{   // units of kg
	return (HH_compHH11NEOIntegratedMass(p->NEO_massMax) -HH_compHH11NEOIntegratedMass(p->NEO_massMin)) / (NEO_integral_flux(p->NEO_massMin) - NEO_integral_flux(p->NEO_massMax));
}

// density units of kg/m^3
double HH_HH11lognormal(double mu, double sigma, double nu, double regolith_dens)
{
	return exp((3.*nu - 1.)*(mu + (3.*nu - 1.)*sigma*sigma/2)) / pow(regolith_dens, 3.*nu - 1.); // unitless factor
}

double HH_HH11dens(double projectile_dense, double target_dense, double nu)
{
	return pow(projectile_dense/target_dense, 3.*nu-1.); // unitless
}


double H_compHH11DensityFactor(input* p, int fluxType)
{
	double densFactor;

	// HiDensMEM, LoDensMEM, NEO
	switch(fluxType)
	{
		case HiDensMEM:
			densFactor = HH_HH11lognormal(p->MEM_hiDens_mu, p->MEM_hiDens_sigma, p->HH11_nu, p->regolith_dens);
			break;

		case LoDensMEM:
			densFactor = HH_HH11lognormal(p->MEM_loDens_mu, p->MEM_loDens_sigma, p->HH11_nu, p->regolith_dens);
			break;

		case NEO:
			densFactor = HH_HH11dens(p->NEO_dens, p->regolith_dens, p->HH11_nu);
			break;

		default:
			cout << "ERROR: Invalid flux type in H_compH11DensityFactor\n\n";
			densFactor = 0.;
	}


	return densFactor;
}


scalingLaw* compute_constants_and_normalization(input* p)
{
	cout << "\nComputing constants and normalization...\n";
	scalingLaw* factors = new scalingLaw;

	factors->m_factor_MEM = H_compHH11GrunMassFactor(p);
	factors->m_factor_NEO = H_compHH11NEOMassFactor(p);
	cout << "mass factor MEM = " << factors->m_factor_MEM << " kg\n";
	cout << "mass factor NEO = " << factors->m_factor_NEO << " kg\n";

	factors->rho_factor_MEMhi = H_compHH11DensityFactor(p, HiDensMEM);
	factors->rho_factor_MEMlo = H_compHH11DensityFactor(p, LoDensMEM);
	factors->rho_factor_NEO   = H_compHH11DensityFactor(p, NEO);

	cout << "density factor MEM hi dens = " << factors->rho_factor_MEMhi << endl;
	cout << "density factor MEM lo dens = " << factors->rho_factor_MEMlo << endl;
	cout << "density factor NEO dens    = " << factors->rho_factor_NEO << endl;

	factors->main_factor_MEMhi = p->HH11_C4 * factors->m_factor_MEM * factors->rho_factor_MEMhi;
	factors->main_factor_MEMlo = p->HH11_C4 * factors->m_factor_MEM * factors->rho_factor_MEMlo;
	factors->main_factor_NEO   = p->HH11_C4 * factors->m_factor_NEO * factors->rho_factor_NEO; 

	cout << "\nmain factor MEM hi dens = " << factors->main_factor_MEMhi << " kg\n";
	cout << "main factor MEM lo dens = " << factors->main_factor_MEMlo << " kg\n";
	cout << "main factor NEO dens    = " << factors->main_factor_NEO << " kg\n";

	return factors;
}


