#ifndef LUNAREJECTA_NEAREARTHOBJECTFLUX_H
#define LUNAREJECTA_NEAREARTHOBJECTFLUX_H

#include <iostream>
#include <vector>
#include <string>

using namespace std;

class lunarEjecta_NearEarthObjectFlux
{
public:
	lunarEjecta_NearEarthObjectFlux();
	~lunarEjecta_NearEarthObjectFlux();

private:

double m_min; // minimum mass range in kg
double m_max; // maximum mass range in kg

double density; // this is the average density of NEO's, should be 3000 kg/m^3


};


#endif 