/*
Title: LunarEjectaSim.cpp
Project: Lunar Meteoroid Ejecta DSNE Environment
Description: A lunar ejecta model based on MEM 3 + NEAs and scaling laws.
Author: Anthony M. DeStefano
Company: NASA/MSFC/EV44
E-mail: anthony.m.destefano@nasa.gov
Office phone: 256-544-3094
Date last edited: 2/20/2020
*/
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
//#include <time.h>
//#include <algorithm> // max(a,b)
//#include <mpi.h>

using namespace std;

#define ANGLEMODTYPE 0 // 0 = no mod, 1 = sin, 2 = sin w/ offset, used in angleModEq

const double PI = 3.141592653589793238; //https://www.piday.org/million/

const double ALPHA0 = 10.*PI/180.; // 10 degrees to radians, used in angleModEq


// useful functions
inline double sqr(double x) { return x*x; }
inline double mag(double x, double y) { return sqrt(x*x + y*y); }
inline double mag2(double x, double y) { return x*x + y*y; }

void initlogVec(vector<double>& v, double vMin, double vMax, int N)
{
	v.resize(N, 0.0);
	for (int i = 0; i < N; ++i)
		v[i] = vMin * pow(vMax/vMin, double(i) / double(N+1));
}

void initLinearVec(vector<double>& v, double vMin, double vMax, int N)
{
	v.resize(N, 0.0);
	for (int i = 0; i < N; ++i)
		v[i] = vMin + (vMax - vMin) * double(i) / double(N+1)
}


//////////////////////////////////////////////////////////////////////////////////////

struct regolith
{
	double density_rho;         // kg/m^3
	double materialStrength_Y0; // Pa
	double porosity_phi0;       // fraction (0 to 1)
	double 
	// see Housen and Holsapple 2011
	double HH11_k;
	double HH11_C1;
	double HH11_mu;
	double HH11_nu;
	double HH11_C4; // derived from previous HH11 values
	// Fit parameters to regolith size distribution, dervied from Carrier 2003
	//  used in C03_RegolithSizeCDF
	double C03_AA;// = 1.02218; // normalization constant, from 0.001 mm to 10 mm
	double C03_A;//  = 0.05480;
	double C03_B;//  = -1.0147;
	double C03_C;//  = 0.3375;
	double C03_D;//  = -0.25181;
};

struct planetesimal
{
	double gravity_g;     // m/s^2
	double radius_r;      // m
	double escapeSpeed_v; // m/s
};


struct projectile
{
	double density_delta;      // kg/m^3
	double mass_m;             // kg
	double speed_U;            // m/s
	double impactAngle_alpha;  // rad (from horizon)
	double azimuthAngle_betaP; // rad (from East)
	// for MEM 3 size density distributions
	double A_low;
	double mu_low;
	double sigma_low;
	double A_high;
	double mu_high;
	double sigma_high;
};

struct ejecta
{
	double distance_D;         // m
	double speed_v;            // m/s
	double zenithAngle_gamma;  // rad (from zenith)
	double azimuthAngle_betaE; // rad (from East)
};


struct impactParameters
{
	regolith*     regolithP;
	planetesimal* planetesimalP;
	projectile*   projectileP;
	ejecta*       ejectaP;
};





// alpha = projectile impact angle from horizon
double angleModEq(double alpha)
{
	switch(ANGLEMODTYPE) {
		case 0: return 1.0;
		case 1: return sin(alpha);
		case 2: return sin(sqrt(sqr(alpha) + sqr(ALPHA0)));
		default: cerr << "In angleModEq: " << ANGLEMODTYPE << "  Invalid angle mode type!\n";
		return 1.0;
	}

}

// See Eq (1) of Vickery 1986
// v = speed of secondary in m/s
// gamma = angle from zenith of secondary from impact location
double distanceVs_SpeedandZenith(impactParameters* P)
{
	return 2.0 * P.planetesimalP->radius_r
		* atan2(sin(2.*P.ejectaP->zenithAngle_gamma)
			* sqr(P.ejectaP->speed_v / P.planetesimalP->escapeSpeed_v),
			1. - 2.*sqr(P.ejectaP->speed_v * sin(P.ejectaP->zenithAngle_gamma) / P.planetesimalP->escapeSpeed_v));
}


double speedVs_DistanceandZenith(impactParameters* P)
{
	return P.planetesimalP->escapeSpeed_v / sqrt(sin(2.*P.ejectaP->zenithAngle_gamma)
		* (1./tan(P.ejectaP->distance_D / (2.*P.planetesimalP->radius_r)) + tan(P.ejectaP->zenithAngle_gamma)));
} 


void zenithVs_DistanceandSpeed(double& gammaP, double& gammaM, impactParameters* P)
{
	double x = P.ejectaP->speed_v / P.planetesimalP->escapeSpeed_v;
	double t1 = sqr(x) / tan(D/(2.*P.planetesimalP->radius_r));
	double t2 = sqr(t1) + 2.*sqr(x) - 1.;


	if (t2 <= 0.)
	{
		gammaP = gammaM = -PI/2.; // force to horizon
		return;
	}
	t2 = sqrt(t2);

	if (t1 > t2)
	{
		gammaP = atan2(1., t1 + t2);
		gammaM = atan2(1., t1 - t2);
		return;
	} else
	{
		gammaP = atan2(1., t1 + t2);
		gammaM = -PI/2.; // force to horizon
		return;
	}

}

////// These are the raw normal distribution functions
// x = [m/s]
inline double logNormalDistribution_H(double x, double mu, double sigma)
{
	return exp(-sqr((log(x) - mu / sigma) / 2.0));
}

inline double logNormalDistribution(double x, double mu, double sigma, double A)
{
	return A * logNormalDistribution_H(x, mu, sigma) / (sigma * x * sqrt(2.0 * PI));
}

void logNormalDistributionVec(vector<double>& F, vector<double>& x, double mu, double sigma, double A)
{
	F.resize(x.size());
	double coeff = A / (sigma * sqrt(2.0 * PI));

	for (int i = 0; i < x.size(); ++i)
		F[i] = coeff * logNormalDistribution_H(x, mu, sigma) / x[i];
}


inline double momentLogNormalDistribution(double x, double mu, double sigma, double A, double alpha)
{
	return A * exp(alpha * mu + sqr(alpha * sigma) / 2.0);
}

void momentLogNormalDistributionVec(vector<double>& F, vector<double>& x, double mu, double sigma, double A, double alpha)
{
	F.resize(x.size());
	for (int i = 0; i < x.size(); ++i)
		F[i] = momentLogNormalDistribution(x[i], mu, sigma, A, alpha);
}
///////////

// Housen and Holsapple 2011
void HH11_EjectedMassGreaterThanV(impactParameters* P)
{
	return P.regolithP->HH11_C4 * P.projectileP->mass_m
		* pow(P.projectileP->speed_U * angleModEq(P.projectileP->impactAngle_alpha)
			/ P.ejectaP->speed_v, 3.*P.regolithP->HH11_mu)
		* pow(P.projectileP->density_delta / P.regolithP->density_rho, 3.*P.regolithP->HH11_nu - 1.);
}

void HH11_EjectedMassInVBin(impactParameters* P,
							double vMin,  // left boundary of secondary ejecta speed bin
							double vMax)  // right boundary of secondary ejecta speed bin
{
	return HH11_EjectedMassGreaterThanV(P, vMax)
		 - HH11_EjectedMassGreaterThanV(P, vMin);
}

// x = [m], needs to be converted to mm
inline double C03_RegolithSizeCDF(impactParameters* P, double x)
{
	return 1.0 - exp(-1.0 / (P.regolithP->C03_A * pow(x*1000.0, P.regolithP->C03_B)
		+ P.regolithP->C03_C * pow(x*1000.0, P.regolithP->C03_D)));
}

void C03_RegolithSizeCDFVec(vector<double>& F, impactParameters* P, vector<double>& x)
{
	F.resize(x.size());
	for (int i = 0; i < x.size(); ++i)
		F[i] = C03_RegolithSizeCDF(P, x[i]);
}


// x = [m], needs to be converted to mm
//  The normalization constant C03_AA is precomputed for the range 0.001 mm to 10 mm
inline double C03_RegolithSizePDF(impactParameters* P, double x)
{
	double denom = P.regolithP->C03_A * pow(x*1000.0, P.regolithP->C03_B)
		+ P.regolithP->C03_C * pow(x*1000.0, P.regolithP->C03_D);
	return -P.regolithP->C03_AA * (P.regolithP->C03_A*P.regolithP->C03_B * pow(x*1000.0, P.regolithP->C03_B-1.0)
			+ P.regolithP->C03_C*P.regolithP->C03_D * pow(x*1000.0, P.regolithP->C03_D-1.0))
		* exp(-1.0 / denom) / sqr(denom);
}

void C03_RegolithSizePDFVec(vector<double>& F, impactParameters* P, vector<double>& x)
{
	F.resize(x.size());
	for (int i = 0; i < x.size(); ++i)
		F[i] = C03_RegolithSizePDF(P, x[i]);
}


int main(int argc, char const *argv[])
{
	vector<double> v;

	
	
	return 0;
}