#include "lunarEjecta_GeneralExpressions.h"
#include "lunarEjecta_SecondaryFlux.h"
#include "lunarEjecta_MeteoroidFlux.h"



using namespace std;

#include <vector>
#include <cmath>

latLon::latLon() {
	lat = 0.0;
	lon = 0.0;
}

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

void latLon::copyLatLon(latLon& pos) { // shallow copy only
	lat = pos.getLatRad();
	lon = pos.getLonRad();
}

void latLon::newLatLon(latLon* pos) {
	pos = new latLon(lat, lon);
}


void latLon::deleteLatLon() {
	delete this;
}

double latLon::getLatRad() {return lat;}
double latLon::getLonRad() {return lon;}

double latLon::getLatDeg() {return lat/DtoR;}
double latLon::getLonDeg() {return lon/DtoR;}

inline double latLon::azm2Bearing(double azm){
	return fmod(2.5*PI - azm, 2.*PI);
}

inline double latLon::bearing2Azm(double bearing) {
	return fmod(2.5*PI - bearing, 2.*PI);
}

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
	return bearing2Azm(this->getBearingInitial(pos));
}

double latLon::getAzmFinal(latLon& pos) {
	return bearing2Azm(this->getBearingFinal(pos));
}



///////////////////////////////

ImpactSites_and_ROI::ImpactSites_and_ROI
	                 (double new_ND,
	                  double new_Nazm,
	                  double new_radius,
	                  latLon& new_ROI)
{
	int i_azm, j_dist;
	double temp_bearing, temp_dist;
	double temp_lat, temp_lon;
	latLon temp_pos;
	ND     = new_ND;
	Nazm   = new_Nazm;
	Ntot   = ND * Nazm;
	radius = new_radius; 

	new_ROI.newLatLon(ROI);

	// generate list of impact sites distributed over the globe
	siteLoc.resize(Ntot);
	D.resize(ND);
	for (j_dist = 0; j_dist < ND; ++j_dist) {
		D[j_dist] = new double(5.0);
		cout << D[j_dist] << endl;
		temp_dist = (j_dist + 1.) / double(ND + 1.0) * 2.*PI; // units of radii
		
		for (i_azm = 0; i_azm < Nazm; ++i_azm) {
			if(j_dist == 0){
				siteAzm[i_azm] = new double;
			}

			temp_bearing = double(i_azm) / double(Nazm) * 2.*PI;
			temp_pos.newLatLon(siteLoc[i_azm + Nazm*j_dist]);
		}
	}
}

ImpactSites_and_ROI::~ImpactSites_and_ROI() {
	// int i,j;
	// ROI->deleteLatLon();
	// ROI = NULL;
	// for (j = 0; j < ND; ++j) {
	// 	delete D[j];
	// 	for (i = 0; i < Nazm; ++i) {
	// 		if(j == 0)
	// 			delete siteAzm[i];
	// 		siteLoc[i+Nazm*j]->deleteLatLon();
	// 		siteLoc[i+Nazm*j] = NULL;
	// 	}
	// }
}
