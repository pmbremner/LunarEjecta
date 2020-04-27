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

double latLon::getLatRad() {return lat;}
double latLon::getLonRad() {return lon;}

double latLon::getLatDeg() {return lat/DtoR;}
double latLon::getLonDeg() {return lon/DtoR;}

// in units of radius
double latLon::getNormDistTo(latLon& pos) {
	double a = sqr(sin((pos.getLatRad() - lat)/2.))
		+ cos(lat) * cos(pos.getLatRad())
		* sqr(sin((pos.getLonRad() - lon)/2.));

	return 2. * asin(sqrt(a));
}

// returns tan(D/(2*rm))
double latLon::getTanDistTo(latLon& pos) {
	double a = sqr(sin((pos.getLatRad() - lat)/2.))
		+ cos(lon) * cos(pos.getLonRad())
		* sqr(sin((pos.getLonRad() - lon)/2.));

	return sqrt(a / (1.-a));
}

double latLon::getBearingTo(latLon& pos) { //broken
	return atan2(sin((pos.getLonRad() - lon)/2.)*cos(pos.getLatRad()),
		cos(lat)*sin(pos.getLatRad()) - sin(lat)*cos(pos.getLatRad())
		*cos((pos.getLonRad() - lon)/2.) );
}

double latLon::getBearingFrom(latLon& pos) { //broken
	return atan2(sin((lon - pos.getLonRad())/2.)*cos(lat),
		cos(pos.getLatRad())*sin(lat) - sin(pos.getLatRad())*cos(lat)
		*cos((lon - pos.getLonRad())/2.) );
}

///////////////////////////////

ImpactSites_and_ROI::ImpactSites_and_ROI() {}
ImpactSites_and_ROI::~ImpactSites_and_ROI() {}
