/*
Title: LunarEjectaSim.cpp
Project: Lunar Meteoroid Ejecta DSNE Environment
Description: A lunar ejecta model based on MEM 3 + NEAs and scaling laws.
Author: Anthony M. DeStefano
Company: NASA/MSFC/EV44
E-mail: anthony.m.destefano@nasa.gov
Office phone: 256-544-3094
Date last edited: 2/10/2020
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

void initlogVector(vector<double>& v, double vMin, double vMax, int N)
{
	v.resize(N);
	for (int i = 0; i < N; ++i)
		v[i] = vMin * pow(vMax/vMin, double(i)/double(N-1));
}





struct regolith
{
	double density_rho;         // kg/m^3
	double materialStrength_Y0; // Pa
	double porosity_phi0;       // fraction (0 to 1)
	double HH11_k;
	double HH11_C1;
	double HH11_mu;
	double HH11_nu;
	double HH11_C4; // derived from previous HH11 values
};

struct planetesimal
{
	double gravity_g;     // m/s^2
	double radius_r;      // m
	double escapeSpeed_v; // m/s
};

struct target
{
	regolith* regolithT;
	planetesimal* planetesimalT;
};





// alpha = projectile impact angle from horizon
double angleModEq(double alpha)
{
	switch(ANGLEMODTYPE) {
		case 0: return 1.0;
		case 1: return sin(alpha);
		case 2: return sin(sqrt(sqr(alpha) + sqr(ALPHA0)));
		default: cerr << "In angleModEq: " << ANGLEMODTYPE << "  Invalid angle mode type!\n";
	}

}

// See Eq (1) of Vickery 1986
// v = speed of secondary in m/s
// gamma = angle from zenith of secondary from impact location
double distanceVs_SpeedandZenith(planetesimal* p, double v, double gamma)
{
	return 2.0 * p->radius_r
		* atan2(sin(2.*gamma) * sqr(v/p->escapeSpeed_v), 1. - 2.*sqr(v*sin(gamma)/p->escapeSpeed_v));
}


double speedVs_DistanceandZenith(planetesimal* p, double D, double gamma)
{
	return p->escapeSpeed_v / sqrt(sin(2.*gamma) * (1./tan(D/(2.*p->radius_r)) + tan(gamma)));
}


void zenithVs_DistanceandSpeed(double& gammaP, double& gammaM, planetesimal* p, double D, double v)
{
	double x = v / p->escapeSpeed_v;
	double t1 = sqr(x) / tan(D/(2.*p->radius_r));
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

// Housen and Holsapple 2011
void HH11_EjectedMassGreaterThanV(target* T,
								  double v,     // secondary ejecta speed
								  double rho,   // target density
								  double m,     // projectile mass
								  double delta, // projectile density
								  double U,     // projectile speed
								  double alpha) // projectile impact angle (from horizon)
{
	return T.regolithT->HH11_C4 * m * pow(U * angleModEq(alpha) / v, 3.*T.regolithT->HH11_mu)
		* pow(delta / rho, 3.*T.regolithT->HH11_nu - 1.);
}

void HH11_EjectedMassInVBin(target* T,
							double rho,   // target density
							double m,     // projectile mass
							double delta, // projectile density
							double U,     // projectile speed
							double alpha  // projectile impact angle (from horizon)
							double vMin,  // left boundary of secondary ejecta speed bin
							double vMax)  // right boundary of secondary ejecta speed bin
{
	return HH11_EjectedMassGreaterThanV(T, vMax, rho, m, delta, U, alpha)
		 - HH11_EjectedMassGreaterThanV(T, vMin, rho, m, delta, U, alpha);
}


int main(int argc, char const *argv[])
{
	vector<double> v;

	
	
	return 0;
}