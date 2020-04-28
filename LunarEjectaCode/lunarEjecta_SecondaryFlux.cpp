#include "lunarEjecta_GeneralExpressions.h"
#include "lunarEjecta_SecondaryFlux.h"
#include "lunarEjecta_MeteoroidFlux.h"



using namespace std;

#include <vector>
#include <cmath>

latLon::latLon(double new_lat, double new_lon) {
	this->setLatDeg(new_lat);
	this->setLonDeg(new_lon);
}

latLon::~latLon() {}

inline void latLon::setLatRad(double new_lat) {
	lat = new_lat;
}
inline void latLon::setLonRad(double new_lon) {
	lon = new_lon;
}

inline void latLon::setLatDeg(double new_lat) {
	lat = DtoR*new_lat;
}
inline void latLon::setLonDeg(double new_lon) {
	lon = DtoR*new_lon;
}

inline double latLon::H_a(latLon& pos) {
	return sqr(sin((pos.getLatRad() - lat)/2.))
		+ cos(lat) * cos(pos.getLatRad())
		* sqr(sin((pos.getLonRad() - lon)/2.));
}

double latLon::getLatRad() {return lat;}
double latLon::getLonRad() {return lon;}

double latLon::getLatDeg() {return lat/DtoR;}
double latLon::getLonDeg() {return lon/DtoR;}

void latLon::dispLatLon() {
	cout << " Latitude = " << lat/DtoR << " degrees, ";
	cout << "Longitude = " << lon/DtoR << " degrees\n";
}

void latLon::dispNormDistTo(latLon& pos) {
	cout << " Distance = " << this->getNormDistTo(pos) << " radii\n";
}

void latLon::dispBearingInitial(latLon& pos) {
	cout << " Initial Bearing = " << this->getBearingInitial(pos)/DtoR << " degrees\n"; 
}

void latLon::dispBearingFinal(latLon& pos) {
	cout << " Final Bearing = " << this->getBearingFinal(pos)/DtoR << " degrees\n"; 
}

void latLon::dispAzmInitial(latLon& pos) {
	cout << " Initial Azimuth = " << this->getAzmInitial(pos)/DtoR << " degrees\n"; 
}

void latLon::dispAzmFinal(latLon& pos) {
	cout << " Final Azimuth = " << this->getAzmFinal(pos)/DtoR << " degrees\n"; 
}

// in units of radius
double latLon::getNormDistTo(latLon& pos) {
	return 2. * asin(sqrt(this->H_a(pos)));
}

// returns tan(D/(2*rm))
double latLon::getTanHalfDistTo(latLon& pos) {
	double a = this->H_a(pos);
	return sqrt(a / (1.-a));
}

double latLon::getSinDistTo(latLon& pos) {
	double a = this->H_a(pos);
	return 2.*sqrt(a * (1.-a));
}
double latLon::getCosDistTo(latLon& pos) {
	return 1. - 2.*this->H_a(pos);
}

double latLon::getBearingInitial(latLon& pos) { 
	return fmod(
		atan2(sin((pos.getLonRad() - lon))*cos(pos.getLatRad()),
		cos(lat)*sin(pos.getLatRad()) - sin(lat)*cos(pos.getLatRad())
		*cos((pos.getLonRad() - lon)) ) + 2.*PI
		, 2.*PI);
}

double latLon::getBearingFinal(latLon& pos) { 
	return fmod(
		atan2(sin((lon - pos.getLonRad()))*cos(lat),
		cos(pos.getLatRad())*sin(lat) - sin(pos.getLatRad())*cos(lat)
		*cos((lon - pos.getLonRad())) ) + PI
		, 2.*PI);
}

double latLon::getAzmInitial(latLon& pos) {
	return fmod(2.5*PI - this->getBearingInitial(pos), 2.*PI);
}

double latLon::getAzmFinal(latLon& pos) {
	return fmod(2.5*PI - this->getBearingFinal(pos), 2.*PI);
}



///////////////////////////////

ImpactSites_and_ROI::ImpactSites_and_ROI() {}
ImpactSites_and_ROI::~ImpactSites_and_ROI() {}
