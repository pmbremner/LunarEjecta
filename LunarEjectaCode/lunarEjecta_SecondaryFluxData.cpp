#include "lunarEjecta_SecondaryFluxData.h"
#include "lunarEjecta_MeteoroidFlux.h"
#include "lunarEjecta_GeneralExpressions.h"

using namespace std;

#include <vector>
#include <cmath>
#include <string>
#include <iomanip>

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

		cout << " D = " << D[j_dist] << endl;

		for (i_azm = 0; i_azm < Nazm; ++i_azm) {
			idx = this->H_idx(i_azm, j_dist);
			// will assume intial bearing first, then translate to azm
			//  this bearing is ROI to siteLoc
			temp_bearing = double(i_azm) * d_bearing;

			if(j_dist == 0){
				// takes init bearing point 1-2, converts to azm, then gets final azm point 2-1
				ROIAzm[i_azm] = siteLoc[idx].init2Final(siteLoc[idx].bearing2Azm(temp_bearing));
				//cout << " ROIAzm = " << ROIAzm[i_azm] << endl;
			}

			// compute lat lon of new site
			ROI.getLatLonFromDistanceBearing(D[j_dist], temp_bearing, siteLoc[idx]);

			// outgoing azimuth from siteLoc to ROI
			siteAzm[idx] = siteLoc[idx].getAzmInitial(ROI);

			//cout << " siteAzm = " << siteAzm[idx] << endl;

			////siteLoc[idx].dispLatLon();
		}
		//cout << endl;
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

int ImpactSites_and_ROI::getND() {
	return ND;
}
int ImpactSites_and_ROI::getNazm() {
	return Nazm;
}
int ImpactSites_and_ROI::getNtot() {
	return Ntot;
}

double ImpactSites_and_ROI::getradius() {
	return radius;
}

double ImpactSites_and_ROI::getsiteLatRad(int i_azm, int j_dist) {
	return siteLoc[H_idx(i_azm, j_dist)].getLatRad();
}

double ImpactSites_and_ROI::getsiteLonRad(int i_azm, int j_dist) {
	return siteLoc[H_idx(i_azm, j_dist)].getLonRad();
}

double ImpactSites_and_ROI::getsiteLatDeg(int i_azm, int j_dist) {
	return siteLoc[H_idx(i_azm, j_dist)].getLatDeg();
}
double ImpactSites_and_ROI::getsiteLonDeg(int i_azm, int j_dist) {
	return siteLoc[H_idx(i_azm, j_dist)].getLonDeg();
}


double ImpactSites_and_ROI::getD(int j_dist) {
	return D[j_dist];
}

double ImpactSites_and_ROI::getsiteAzm(int i_azm, int j_dist) {
	return siteAzm[H_idx(i_azm, j_dist)];
}
double ImpactSites_and_ROI::getROIAzm(int i_azm) {
	return ROIAzm[i_azm];
}

double ImpactSites_and_ROI::getROI_SA() {
	return ROI_SA;
}
double ImpactSites_and_ROI::getROI_radius() {
	return ROI_radius;
}
double ImpactSites_and_ROI::getsite_SA(int j_dist) {
	return site_SA[j_dist];
}


double ImpactSites_and_ROI::getDbeta(double D0, double D1) // D's in units of circumference (2*Pi*rm)
{
	double D = (D0 + D1) / 2.; // units of circ (2pi*rm)
	// A spherical right triangle
	double cos_h_2 = sqr(cos(2.*PI*D + ROI_radius/radius) * cos(ROI_radius/radius));

	return acos((cos(2.*ROI_radius/radius) - cos_h_2) / (1. - cos_h_2));
}


//////////////////////////////////////
//////////////////////////////////////

GeneralIntegralFluxOutput::GeneralIntegralFluxOutput
                             (string oType,
							  double new_xMin,
							  double new_xMax,
							  int new_Nx,
							  int new_xScale,
							  int new_NSetsXY,
							  vector<double>& new_setMin,
							  vector<double>& new_setMax) {
	outputType = oType;
	xMin       = new_xMin;
	xMax       = new_xMax;
	Nx         = new_Nx;
	xScale     = new_xScale;
	NSetsXY    = new_NSetsXY;
	copyVector(new_setMin, setMin, NSetsXY);
	copyVector(new_setMax, setMax, NSetsXY);
	// setMin = new_setMin;
	// setMax = new_setMax;

	int NXY = (NSetsXY < 2 ? 1 : NSetsXY);

	xData.resize(NXY);
	for (int i = 0; i < NXY; ++i){
		xData[i].resize(Nx);
		fill(xData[i].begin(), xData[i].end(), 0.0);
	}

	this->dispOutputType();
	this->dispXScaleType();
	this->dispSetMinMax();
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

void GeneralIntegralFluxOutput::dispSetMinMax() {
	cout << " NSetsXY = " << NSetsXY << endl;
	if(NSetsXY > 0) {
		for (int i = 0; i < NSetsXY; ++i)
			cout << fixed << setprecision(5) << setMin[i] << '\t' << setMax[i] << endl;
	} else {
		cout  << "   No set min-max to display\n\n";
	}
}

void GeneralIntegralFluxOutput::get_xDataPointer(vector<vector<double>>* p_xData)
{
	p_xData = &xData;
}

//////////////////////////////////////
//////////////////////////////////////

MassLimitedIntegralFluxVsMass::MassLimitedIntegralFluxVsMass
 (double new_xMin, double new_xMax, int new_Nx, int new_xScale, int new_NSetsXY, vector<double>& new_setMin, vector<double>& new_setMax)
 : GeneralIntegralFluxOutput("MassLimitedIntegralFluxVsMass", new_xMin, new_xMax, new_Nx, new_xScale, new_NSetsXY, new_setMin, new_setMax)
{
	// compute normalization constant, fraction greater than m
	fraction_GT_m.resize(Nx);
	fill(fraction_GT_m.begin(), fraction_GT_m.end(), 0.0); // REPLACE!
}

MassLimitedIntegralFluxVsMass::~MassLimitedIntegralFluxVsMass() {}


void MassLimitedIntegralFluxVsMass::saveFluxToFile(string fn) {}

// since we integrate out the altitude and azimuth,
//  we just need a reasonable value here
int MassLimitedIntegralFluxVsMass::getNalt()
{
	return 30; // aka, a resolution of 3 degree increments over 0 to 90 degrees
}

// Again, since we're integrating this out,
//  we just need a reasonable value
int MassLimitedIntegralFluxVsMass::getNazm()
{
	return 36; // aka, a resolution of 10 degree increments over 0 to 360 degrees
}

// This is also integrated out
int MassLimitedIntegralFluxVsMass::getNvel()
{
	return 40; 
}


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
 (double new_xMin, double new_xMax, int new_Nx, int new_xScale, int new_NSetsXY, vector<double>& new_setMin, vector<double>& new_setMax)
 : GeneralIntegralFluxOutput("SizeLimitedIntegralFluxVsSpeed", new_xMin, new_xMax, new_Nx, new_xScale, new_NSetsXY, new_setMin, new_setMax)
{// compute normalization constant, fraction greater than m
	fraction_GT_d.resize(new_NSetsXY);
	fill(fraction_GT_d.begin(), fraction_GT_d.end(), 0.0); // REPLACE!
}

SizeLimitedIntegralFluxVsSpeed::~SizeLimitedIntegralFluxVsSpeed() {}



void SizeLimitedIntegralFluxVsSpeed::saveFluxToFile(string fn) {}

// since we integrate out the altitude and azimuth,
//  we just need a reasonable value here
int SizeLimitedIntegralFluxVsSpeed::getNalt()
{
	return 30; // aka, a resolution of 3 degree increments over 0 to 90 degrees
}

// Again, since we're integrating this out,
//  we just need a reasonable value
int SizeLimitedIntegralFluxVsSpeed::getNazm()
{
	return 36; // aka, a resolution of 10 degree increments over 0 to 360 degrees
}

// This is also integrated out
int SizeLimitedIntegralFluxVsSpeed::getNvel()
{
	return Nx; 
}


void SizeLimitedIntegralFluxVsSpeed::updateFlux(double flux, double alt, double azm, double speed) {
	// update flux (can probably factor this out to be more efficient...)
	int n, m;
	for (m = 0; m < NSetsXY; ++m)
		for (n = 0; n < Nx; ++n)
			xData[m][n] += flux * fraction_GT_d[m];
}


//////////////////////////////////////
//////////////////////////////////////

MassLimitedIglooIntegratedFlux::MassLimitedIglooIntegratedFlux
 (double new_xMin, double new_xMax, int new_angleRes, int new_xScale, int new_NSetsXY, vector<double>& new_setMin, vector<double>& new_setMax)
 : GeneralIntegralFluxOutput("MassLimitedIglooIntegratedFlux", new_xMin, new_xMax, H_getIglooNx(new_angleRes), new_xScale, new_NSetsXY, new_setMin, new_setMax)
{
	angleRes = new_angleRes;
	cout << " Angle Resolution = " << angleRes << endl;
	cout << " Nx = " << Nx << endl;

	int cur_ID = 1, i, j, jN;
	igloo_ID.resize(Nx);
	igloo_I.resize(Nx);
	igloo_J.resize(Nx);
	igloo_PHI1.resize(Nx);
	igloo_PHI2.resize(Nx);
	igloo_THETA1.resize(Nx);
	igloo_THETA2.resize(Nx);
	igloo_PHIavg.resize(Nx);
	igloo_THETAavg.resize(Nx);

	for (i = 1; i <= 180/angleRes; ++i)
	{
		jN = H_getIglooJN(angleRes, (i-1)*angleRes,  i*angleRes);
		for (j = 1; j <= jN; ++j)
		{
			igloo_ID[cur_ID-1]       = cur_ID;
			igloo_I[cur_ID-1]        = i;
			igloo_J[cur_ID-1]        = j;
			igloo_PHI1[cur_ID-1]     = -90. + (i-1)*angleRes;
			igloo_PHI2[cur_ID-1]     = -90. + i*angleRes;
			igloo_THETA1[cur_ID-1]   = 360. * double(j-1.) / double(jN);
			igloo_THETA2[cur_ID-1]   = 360. * double(j) / double(jN);
			igloo_PHIavg[cur_ID-1]   = asin(0.5*(sin(igloo_PHI1[cur_ID-1] * DtoR) + sin(igloo_PHI2[cur_ID-1] * DtoR))) / DtoR;
			igloo_THETAavg[cur_ID-1] = 0.5* (igloo_THETA1[cur_ID-1] + igloo_THETA2[cur_ID-1]);
	
			// cout << igloo_ID[cur_ID-1] << ' ' << igloo_I[cur_ID-1] << ' ' << igloo_J[cur_ID-1]
			// 	 << ' ' << igloo_PHI1[cur_ID-1] << ' ' << igloo_PHI2[cur_ID-1] << ' '
			// 	 << igloo_THETA1[cur_ID-1] << ' ' << igloo_THETA2[cur_ID-1] << ' '
			// 	 << igloo_PHIavg[cur_ID-1] << ' ' << igloo_THETAavg[cur_ID-1] << endl;

			cur_ID++;	
		}
	}

	// this->updateFlux(1.0, 23., 12., 0.55);// test
	// this->updateFlux(2.0, 86., 245., 1.2);// test
}

MassLimitedIglooIntegratedFlux::~MassLimitedIglooIntegratedFlux() {}



void MassLimitedIglooIntegratedFlux::saveFluxToFile(string fn)
{
	int i, j;
	ofstream file;
	file.open(fn);

	cout << " Saving MassLimitedIglooIntegratedFlux to file: " << fn << endl;

	file << "#\n#\n#\n#\n#\n#\n#\n#\n"; // placeholder header

	for (i = 0; i < Nx; ++i)
	{
		file << igloo_ID[i] << ' ';
		file << igloo_I[i] << ' ';
		file << igloo_J[i] << ' ';
		file << igloo_PHI1[i] << ' ';
		file << igloo_PHI2[i] << ' ';
		file << igloo_THETA1[i] << ' ';
		file << igloo_THETA2[i] << ' ';
		file << igloo_PHIavg[i] << ' ';
		file << igloo_THETAavg[i] << ' ';

		// output data
		for (j = 0; j < NSetsXY; ++j)
		{
			file << xData[j][i] << ' ';
		}
		file << endl;
	}

	file.close();
}


int MassLimitedIglooIntegratedFlux::getNalt()
{
	return 90/angleRes;
}

// This is inefficient, but it will get it to work without much effort
int MassLimitedIglooIntegratedFlux::getNazm()
{
	return 360/angleRes;
}

// This is also integrated out
int MassLimitedIglooIntegratedFlux::getNvel()
{
	return NSetsXY; 
}

// same as MEM_iglooAvg::getFlux_atAngleVel
// alt, azm units of degrees, vel units of km/s
void MassLimitedIglooIntegratedFlux::updateFlux(double flux, double alt, double azm, double vel) {
	// find index for correct phi (alt)
	int idx_min = 0, idx_max = Nx-1, idx_mid = (idx_max - idx_min)/2;
	int row_idx, col_idx, Nazm;
	double dAzm;

	// make azm from 0 to 360
	azm = fmod(azm + 360.0, 360.0);

	// simple binary search algorithm
	while (igloo_PHI1[idx_min] != igloo_PHI1[idx_mid])
	{

		if (alt < igloo_PHI1[idx_mid])
		{
			idx_max = idx_mid;
		} else
		{
			idx_min = idx_mid;
		}
		idx_mid = (idx_max + idx_min)/2;

		// cout << idx_min << ' ' << idx_mid << ' ' << idx_max << "   ";
		// cout << igloo_PHI1[idx_min] << ' ' << igloo_PHI1[idx_mid] << ' ' << igloo_PHI1[idx_max] << ' ' << alt << endl;
	}

	row_idx = idx_min - igloo_J[idx_min] + 1;
	//cout << "row_idx = " << row_idx << endl;


	// correct index for correct theta (azm)
	dAzm = igloo_THETA2[row_idx];
	Nazm = round(360.0 / dAzm);

	row_idx += int(azm / dAzm);
	// cout << "row_idx = " << row_idx << endl;
	// cout << " NSetsXY = " << NSetsXY << endl;
	// find index for vel and return the flux
	col_idx = 0;
	while (col_idx < NSetsXY && !(vel >= setMin[col_idx] && vel <= setMax[col_idx])) {
		//cout << col_idx << " | " << setMin[col_idx] << ' ' << setMax[col_idx] << endl;
		col_idx++;
	}

	if(col_idx < NSetsXY) {
		// store flux
		xData[col_idx][row_idx] += flux;
	} // else, don't store flux since the speed is out-of-bounds

	// cout << "xData[col_idx][row_idx] = " << xData[col_idx][row_idx] << endl;
	// cout << "col = " << col_idx << " |  row = " << row_idx << endl;
}


int MassLimitedIglooIntegratedFlux::H_getIglooNx(int aRes) {// input in units of degrees
	int i, sum = 0;
	for (i = 1; i <= 180/aRes; ++i) // alt
	{
		//cout << i << " | " << H_getIglooJN(aRes, (i-1)*aRes,  i*aRes) << endl;
		sum += H_getIglooJN(aRes, (i-1)*aRes, i*aRes);
	}

	return sum;
}

inline int MassLimitedIglooIntegratedFlux::H_getIglooJN(int aRes, double azm1, double azm2) {// input in units of degrees
	//cout << "azm1 = " << azm1 << " | azm2 = " << azm2 << endl;
	//cout << "   " << 180. * (azm2 - azm1) * (sin(azm1*DtoR) + sin(azm2*DtoR)) << endl;
	return round(fabs(180. * (azm2 - azm1) * (sin(azm1*DtoR) + sin(azm2*DtoR)) / sqr(double(aRes)) ) );
}