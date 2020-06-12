#include "lunarEjecta_NearEarthObjectFlux.h"
#include "lunarEjecta_GeneralExpressions.h"
#include "lunarEjecta_MeteoroidFlux.h"


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

lunarEjecta_NearEarthObjectFlux::lunarEjecta_NearEarthObjectFlux
		(string fn,
		 double new_m_min,
		 double new_m_max,
		 int densType,
		 double userDefDens,
  /* params for MEM_LatData part */
		 string dn, double lMin, double lMax, int NL) : MEM_LatData<MEM_HiDensityIglooAvg>(dn, lMin, lMax, NL)
{
	cout << "--Init lunarEjecta_NearEarthObjectFlux starting--\n";
	filename = fn;

	m_min = new_m_min;
	m_max = new_m_max;

	if (densType == defaultDens)
	{
		density = 3000.; // units of kg/m^3
		cout << " Default NEO density = " << density << " kg/m^3\n";
	}
	else if (densType == userDefined)
	{
		density = userDefDens;
		cout << " User defined NEO density = " << density << " kg/m^3\n";
	}
	else
	{
		density = 0.;
		cout << "ERROR: lunarEjecta_NearEarthObjectFlux invalid density type\n";
	}

	// compute mass factor
	massFactor = H_compNEOmassFactor(100);
	cout << scientific << " m_min, m_max = " << m_min << " kg, " << m_max;
	cout << scientific << " | mass factor = " << massFactor << " kg/m^2/yr\n";


	// read speed distribution from file (similar to )
	int count = 0;
	char C_densLine[64];
	double tdouble1, tdouble2;

	ifstream file;
	file.open(filename);
	cout << filename << endl;

	file.ignore(256, '\n');
	while(file.getline(C_densLine, 64, '\n')) {
		//cout << C_densLine << " | ";
		string tstr(C_densLine);
		stringstream SS_double(tstr);
		SS_double >> tdouble1 >> tdouble2;
		speed.push_back(tdouble1);
		fraction.push_back(tdouble2);
		count++;
		cout << count << ' ' << tdouble1 << '\t' << tdouble2 << endl;
	}
	Nv = count;
	cout << " distribution rows = " << count << endl;

	vMin = speed[0];
	vMax = speed[Nv-1];

	file.close();
}



lunarEjecta_NearEarthObjectFlux::~lunarEjecta_NearEarthObjectFlux() {}


double lunarEjecta_NearEarthObjectFlux::getNEOflux_atSpeed(double v)
{

}

double lunarEjecta_NearEarthObjectFlux::getmMin() {}
double lunarEjecta_NearEarthObjectFlux::getmMax() {}
double lunarEjecta_NearEarthObjectFlux::getDens() {}
double lunarEjecta_NearEarthObjectFlux::getvMin() {}
double lunarEjecta_NearEarthObjectFlux::getvMax() {}
int    lunarEjecta_NearEarthObjectFlux::getNv()   {}




double lunarEjecta_NearEarthObjectFlux::g_NEOflux(double m)
{
	return 2.89E-11 * pow(m, -0.9);
}

// computed in mathematica
double lunarEjecta_NearEarthObjectFlux::Dg_NEOflux(double m)
{
	return 2.601E-11 * pow(m, -0.9);
}

// integrating with HH11 mass term
// similar to H_compH11GrunMassFactor, using a power-law fit on the grid and integrating
double lunarEjecta_NearEarthObjectFlux::H_compNEOmassFactor(int N = 14) // N ~100 is usually sufficient
{
	int i;
	double Ai, Bi, sum = 0.;
	vector<double> m;
	vector<double> f;
	f.resize(N);

	// prepare x and y vectors for power law fits
	logspace(m, m_min, m_max, N, 0, N);
	for (i = 0; i < N; ++i)
		f[i] = Dg_NEOflux(m[i]);

	// compute integral by breaking up in log spacings and approx by power law
	for (i = N-2; i >= 0; i--) // start at low end of mass to sum smaller #'s first (less floating-point error)
	{
		Bi = log10(f[i+1] / f[i]) / log10(m[i+1] / m[i]); 
		Ai = (f[i] + f[i+1]) / (pow(m[i], Bi) + pow(m[i+1], Bi));

		sum += Ai/(Bi+1.) * (pow(m[i+1], Bi+1) - pow(m[i], Bi+1));
	}
	return sum;	
}