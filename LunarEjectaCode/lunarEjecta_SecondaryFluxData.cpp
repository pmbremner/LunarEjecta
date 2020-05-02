#include "lunarEjecta_GeneralExpressions.h"
#include "lunarEjecta_SecondaryFluxData.h"
#include "lunarEjecta_MeteoroidFlux.h"

using namespace std;

#include <vector>
#include <cmath>
#include <string>

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

inline int ImpactSites_and_ROI::H_idx(int i_azm, int j_dist) {
	return i_azm + Nazm*j_dist;
}

///////////////////////////////

ImpactSites_and_ROI::ImpactSites_and_ROI
	                 (double new_ND,
	                  double new_Nazm,
	                  double new_radius,
	                  latLon& new_ROI)
{
	int i_azm, j_dist, idx;
	double temp_bearing, temp_dist, d_bearing, d_lat;
	double temp_lat, temp_lon;
	double SA_check = 0.0;

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
	site_SA.resize(ND);

	// generate list of impact sites distributed over the globe
	d_bearing = 2.*PI/ double(Nazm);
	d_lat = PI / double(ND + 1.);

	for (j_dist = 0; j_dist < ND; ++j_dist) {
		D[j_dist] = (j_dist + 1.) * d_lat; // units of radii
		//cout << j_dist << " D = " << D[j_dist] << endl;

		// compute site surface areas
		site_SA[j_dist] = sqr(radius) * d_bearing
			* fabs(cos((0.5+j_dist)*d_lat) - cos((1.5+j_dist)*d_lat));

		SA_check += site_SA[j_dist] * Nazm;

		for (i_azm = 0; i_azm < Nazm; ++i_azm) {
			idx = this->H_idx(i_azm, j_dist);
			// will assume intial bearing first, then translate to azm
			//  this bearing is ROI to siteLoc
			temp_bearing = double(i_azm) * d_bearing;

			if(j_dist == 0){
				// takes init bearing point 1-2, converts to azm, then gets final azm point 2-1
				ROIAzm[i_azm] = siteLoc[idx].init2Final(siteLoc[idx].bearing2Azm(temp_bearing));
			}

			// compute lat lon of new site
			ROI.getLatLonFromDistanceBearing(D[j_dist], temp_bearing, siteLoc[idx]);

			// outgoing azimuth from siteLoc to ROI
			siteAzm[idx] = siteLoc[idx].getAzmInitial(ROI);

			//siteLoc[idx].dispLatLon();
		}
	}
	// compute ROI surface area
	ROI_radius = 2. * radius * sin(0.5*d_lat/2.0);
	ROI_SA = PI * sqr(ROI_radius);

	SA_check += 2.*ROI_SA;

	// cout << " ROI radius = " << ROI_radius << " m" << endl;
	// cout << " ROI SA = " << ROI_SA << " m^2" << endl;
	cout << " Total surface area = " << SA_check << " m^2" << endl;
}

ImpactSites_and_ROI::~ImpactSites_and_ROI() {}

//////////////////////////////////////
//////////////////////////////////////

GeneralIntegralFluxOutput::GeneralIntegralFluxOutput
                             (string oType,
							  double new_xMin,
							  double new_xMax,
							  int new_Nx,
							  int new_xScale,
							  int new_NSetsXY,
							  vector<double> new_setMin,
							  vector<double> new_setMax) {
	outputType = oType;
	xMin       = new_xMin;
	xMax       = new_xMax;
	Nx         = new_Nx;
	xScale     = new_xScale;
	NSetsXY    = new_NSetsXY;
	setMin = new_setMin;
	setMax = new_setMax;

	int NXY = (NSetsXY < 2 ? 1 : NSetsXY);

	xData.resize(NXY);
	for (int i = 0; i < NXY; ++i){
		xData[i].resize(Nx);
		fill(xData[i].begin(), xData[i].end(), 0.0);
	}

	this->dispOutputType();
	this->dispXScaleType();
}

GeneralIntegralFluxOutput::~GeneralIntegralFluxOutput() {}

void GeneralIntegralFluxOutput::dispOutputType() {
	cout << "Output type = " << outputType << endl;
}

void GeneralIntegralFluxOutput::dispXScaleType() {
	cout << "xScale type = ";
	if(xScale == linearScale) {
		cout << "linearScale\n";
	} else if (xScale == log10Scale){
		cout << "log10Scale\n";
	}
}

//////////////////////////////////////
//////////////////////////////////////

MassLimitedIntegralFluxVsMass::MassLimitedIntegralFluxVsMass
 (double new_xMin, double new_xMax, int new_Nx, int new_xScale, int new_NSetsXY, vector<double> new_setMin, vector<double> new_setMax)
 : GeneralIntegralFluxOutput("MassLimitedIntegralFluxVsMass", new_xMin, new_xMax, new_Nx, new_xScale, new_NSetsXY, new_setMin, new_setMax)
{
	// compute normalization constant, fraction greater than m
	fraction_GT_m.resize(Nx);
	fill(fraction_GT_m.begin(), fraction_GT_m.end(), 0.0); // REPLACE!
}

MassLimitedIntegralFluxVsMass::~MassLimitedIntegralFluxVsMass() {}


void MassLimitedIntegralFluxVsMass::saveFluxToFile(string fn) {}


void MassLimitedIntegralFluxVsMass::updateFlux(double flux, double alt, double azm, double speed) {
	// Find which speed range to dump it into, assuming they don't overlap
	int i_vRange, n;
	for (n = 0; n < NSetsXY; ++n)
	{
		if (speed >= setMin[n] && speed <= setMax[n])
		{
			i_vRange = n;
			n = NSetsXY; // to break out of for-loop
		}
	}
	
	// update flux (can probably factor this out to be more efficient...)
	for (n = 0; n < Nx; ++n)
		xData[i_vRange][n] += flux * fraction_GT_m[n];

}


//////////////////////////////////////
//////////////////////////////////////

SizeLimitedIntegralFluxVsSpeed::SizeLimitedIntegralFluxVsSpeed
 (double new_xMin, double new_xMax, int new_Nx, int new_xScale, int new_NSetsXY, vector<double> new_setMin, vector<double> new_setMax)
 : GeneralIntegralFluxOutput("SizeLimitedIntegralFluxVsSpeed", new_xMin, new_xMax, new_Nx, new_xScale, new_NSetsXY, new_setMin, new_setMax)
{}

SizeLimitedIntegralFluxVsSpeed::~SizeLimitedIntegralFluxVsSpeed() {}



void SizeLimitedIntegralFluxVsSpeed::saveFluxToFile(string fn) {}


void SizeLimitedIntegralFluxVsSpeed::updateFlux(double flux, double alt, double azm, double speed) {}


//////////////////////////////////////
//////////////////////////////////////

MassLimitedIglooIntegratedFlux::MassLimitedIglooIntegratedFlux
 (double new_xMin, double new_xMax, int new_Nx, int new_xScale, int new_NSetsXY, vector<double> new_setMin, vector<double> new_setMax)
 : GeneralIntegralFluxOutput("MassLimitedIglooIntegratedFlux", new_xMin, new_xMax, new_Nx, new_xScale, new_NSetsXY, new_setMin, new_setMax)
{}

MassLimitedIglooIntegratedFlux::~MassLimitedIglooIntegratedFlux() {}



void MassLimitedIglooIntegratedFlux::saveFluxToFile(string fn) {}


void MassLimitedIglooIntegratedFlux::updateFlux(double flux, double alt, double azm, double speed) {}

