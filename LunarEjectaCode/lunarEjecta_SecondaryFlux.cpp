#include "lunarEjecta_GeneralExpressions.h"
#include "lunarEjecta_SecondaryFlux.h"
#include "lunarEjecta_MeteoroidFlux.h"



using namespace std;

#include <vector>
#include <cmath>

latLon::latLon(double new_lat, double new_lon) {
	lat = DtoR*new_lat;
	lon = DtoR*new_lon;
}

latLon::~latLon() {}

double latLon::getLat() {return lat/DtoR;}
double latLon::getLon() {return lon/DtoR;}

// in units of radius
double latLon::getNormDistTo(latLon& toLocation) {
	double a = sqr(sin((toLocation.getLat() - lat)/2.))
		+ cos(lat) * cos(toLocation.getLat())
		* sqr(sin((toLocation.getLon() - lon)/2.));

	return 2. * asin(sqrt(a));
}

// returns tan(D/(2*rm))
double latLon::getTanDistTo(latLon& toLocation) {
	double a = sqr(sin((toLocation.getLat() - lat)/2.))
		+ cos(lon) * cos(toLocation.getLon())
		* sqr(sin((toLocation.getLon() - lon)/2.));

	return sqrt(a / (1.-a));
}



///////////////////////////////

ImpactSites_and_ROI::ImpactSites_and_ROI() {}
ImpactSites_and_ROI::~ImpactSites_and_ROI() {}
