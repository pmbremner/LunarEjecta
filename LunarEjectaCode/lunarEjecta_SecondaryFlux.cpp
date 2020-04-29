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

void latLon::copyLatLon(double new_lat, double new_lon) {
	lat = new_lat;
	lon = new_lon;
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

inline double latLon::init2Final(double angle) {
	return fmod(angle + PI, 2.*PI);
}

inline double latLon::final2Init(double angle) {
	return fmod(angle + PI, 2.*PI);
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

double latLon::getBearingInitial(latLon& pos2) { 
	return fmod(
		atan2(sin((pos2.getLonRad() - lon))*cos(pos2.getLatRad()),
		cos(lat)*sin(pos2.getLatRad()) - sin(lat)*cos(pos2.getLatRad())
		*cos((pos2.getLonRad() - lon)) ) + 2.*PI
		, 2.*PI);
}

double latLon::getBearingFinal(latLon& pos1) { 
	return fmod(
		atan2(sin((lon - pos1.getLonRad()))*cos(lat),
		cos(pos1.getLatRad())*sin(lat) - sin(pos1.getLatRad())*cos(lat)
		*cos((lon - pos1.getLonRad())) ) + PI
		, 2.*PI);
}

double latLon::getAzmInitial(latLon& pos2) {
	return bearing2Azm(this->getBearingInitial(pos2));
}

double latLon::getAzmFinal(latLon& pos1) {
	return bearing2Azm(this->getBearingFinal(pos1));
}

// D is in units of radii
void latLon::getLatLonFromDistanceBearing(double D, double theta, latLon& pos2) {
	double sinLat2 = sin(lat)*cos(D) + cos(lat)*sin(D)*cos(theta);

	pos2.setLatRad(asin(sinLat2));
	pos2.setLonRad(lon + atan2(sin(theta)*sin(D)*cos(lat), cos(D) - sin(lat)*sinLat2) );
}



///////////////////////////////

ImpactSites_and_ROI::ImpactSites_and_ROI
	                 (double new_ND,
	                  double new_Nazm,
	                  double new_radius,
	                  latLon& new_ROI)
{
	int i_azm, j_dist, idx;
	double temp_bearing, temp_dist;
	double temp_lat, temp_lon;

	ND     = new_ND;
	Nazm   = new_Nazm;
	Ntot   = ND * Nazm;
	radius = new_radius; 

	ROI.copyLatLon(new_ROI);

	cout << "ROI latlon:\n";
	ROI.dispLatLon();
	cout << endl;


	// allocate more space in vectors
	siteLoc.resize(Ntot);
	D.resize(ND);
	siteAzm.resize(Ntot);
	ROIAzm.resize(Nazm);
	site_SA.resize(Ntot);

	// generate list of impact sites distributed over the globe
	for (j_dist = 0; j_dist < ND; ++j_dist) {
		D[j_dist] = (j_dist + 1.) / double(ND + 1.0) *PI; // units of radii
		cout << j_dist << " D = " << D[j_dist] << endl;
		for (i_azm = 0; i_azm < Nazm; ++i_azm) {
			idx = i_azm + Nazm*j_dist;
			// will assume intial bearing first, then translate to azm
			//  this bearing is ROI to siteLoc
			temp_bearing = double(i_azm) / double(Nazm) * 2.*PI;

			if(j_dist == 0){
				// takes init bearing point 1-2, converts to azm, then gets final azm point 2-1
				ROIAzm[i_azm] = siteLoc[idx].init2Final(siteLoc[idx].bearing2Azm(temp_bearing));
				//cout << ROIAzm[i_azm]/DtoR << endl;
			}

			// compute lat lon of new site
			ROI.getLatLonFromDistanceBearing(D[j_dist], temp_bearing, siteLoc[idx]);

			// outgoing azimuth from siteLoc to ROI
			siteAzm[idx] = siteLoc[idx].getAzmInitial(ROI);

			siteLoc[idx].dispLatLon();
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
