#ifndef LUNAREJECTA_NEAREARTHOBJECTFLUX_H
#define LUNAREJECTA_NEAREARTHOBJECTFLUX_H

#include "lunarEjecta_MeteoroidFlux.h"

#include <iostream>
#include <vector>
#include <string>

using namespace std;

enum NEOdensType{defaultDens, userDefined};

// Note: the idea here is to create a new class that inherits
// the MEM_LatData<MEM_HiDensityIglooAvg> class, so that we
// don't have to reinvent the wheel, let alone the whole blasted car.
// This assumes that the data read in is preprocessed to be renormalized
// using a separate script. Therefore, all we need to provide is similar
// init paramters for a MEM_LatData<> type class, which includes the directory
// of latitude dependent, high density igloo flux files
//   See: https://stackoverflow.com/questions/8810224/inheriting-from-a-template-class-in-c
class lunarEjecta_NearEarthObjectFlux : public MEM_LatData<MEM_HiDensityIglooAvg> 
{
public:
	lunarEjecta_NearEarthObjectFlux(string fn,
		                            double new_m_min,
		                            double new_m_max,
		                            int densType,
		                            double userDefDens,
  /* params for MEM_LatData part */ string dn_NEO, double lMin_NEO, double lMax_NEO, int NL_NEO);
	~lunarEjecta_NearEarthObjectFlux();

	double getMassFluxNEO_atAngleVelLat(double alt, double azm, double vel, double lat);

	//double getNEOflux_atSpeed(double v);, not using for now

	double getmMin();
	double getmMax();
	double getDens();
	double getvMin();
	double getvMax();
	int    getNv();

private:
	double g_NEOflux(double m);

	double Dg_NEOflux(double m);

	// integrating with HH11 mass term
	double H_compNEOmassFactor(int N);


	string filename; // file name of speed-fraction distribution

	double massFactor; // units of kg/m^2/yr

	double m_min; // minimum mass range in kg
	double m_max; // maximum mass range in kg

	double density; // kg/m^3, this is the average density of NEO's, should be 3000 kg/m^3 for default (see Brown 2002)

	double vMin; // units of km/s
	double vMax; // units of km/s
	int Nv; // number of vel bins, will be defined once file is read

	vector<double> speed;
	vector<double> fraction;
};


#endif 