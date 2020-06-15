#ifndef LUNAREJECTA_SECONDARYFLUXDATA_H
#define LUNAREJECTA_SECONDARYFLUXDATA_H

#include <iostream>
#include <vector>
#include <string>

using namespace std;


class latLon
{
public:
	latLon(); // default, 0,0
	latLon(double new_lat, double new_lon); // input in degrees
	~latLon();

	inline void setLatRad(double new_lat);
	inline void setLonRad(double new_lon);
	inline void setLatDeg(double new_lat);
	inline void setLonDeg(double new_lon);

	void copyLatLon(latLon& pos); // from "param" to "this"
	void copyLatLon(double new_lat, double new_lon); // from "param" to "this"
	void newLatLon(latLon* pos); // from "this" to "param"
	void deleteLatLon();

	double getLatRad();
	double getLonRad();

	double getLatDeg();
	double getLonDeg();

	inline double azm2Bearing(double azm);
	inline double bearing2Azm(double bearing);
	inline double init2Final(double angle); // init angle from 1 to 2, to final from 2 to 1
	inline double final2Init(double angle);

	// all angle outputs in degrees for readability, even though lat & lon are in radians
	void dispLatLon();
	void dispNormDistTo(latLon& pos);
	void dispBearingInitial(latLon& pos);
	void dispBearingFinal(latLon& pos);
	void dispAzmInitial(latLon& pos);
	void dispAzmFinal(latLon& pos);

	double getNormDistTo(latLon& pos);
	double getTanHalfDistTo(latLon& pos);
	double getSinDistTo(latLon& pos);
	double getCosDistTo(latLon& pos);

	// outputs in N CW angle
	double getBearingInitial(latLon& pos2);
	double getBearingFinal(latLon& pos1);

	// outputs in E CCW angle
	double getAzmInitial(latLon& pos2);
	double getAzmFinal(latLon& pos1);

	void getLatLonFromDistanceBearing(double D, double theta, latLon& pos2); // returns new lat lon into pos2
	
private:
	inline double H_a(latLon& pos);

	double lat; // radians
	double lon; // radians
};


class ImpactSites_and_ROI
{
public:
	ImpactSites_and_ROI(double new_ND,     // total number of distance increments
		                double new_Nazm,   // total number of azimuth increments
		                double new_radius, // radius of Moon
		                latLon& new_ROI);  // lat-lon location of Region-Of-Interest
	~ImpactSites_and_ROI();
	
	int getND();
	int getNazm();
	int getNtot();

	double getradius();

	double getsiteLatRad(int i_azm, int j_dist);
	double getsiteLonRad(int i_azm, int j_dist);

	double getsiteLatDeg(int i_azm, int j_dist);
	double getsiteLonDeg(int i_azm, int j_dist);

	double getD(int j_dist);
	double getsiteAzm(int i_azm, int j_dist);
	double getROIAzm(int i_azm);

	double getROI_SA();
	double getROI_radius();
	double getsite_SA(int j_dist);

	double getDbeta(double D0, double D1); // D's in units of circumference (2*Pi*rm)

private:
	inline int H_idx(int i_azm, int j_dist);

	int ND;
	int Nazm;
	int Ntot; // ND * Nazm

	double radius; // units of m
	
	latLon ROI;              // region of interest location
	vector<latLon> siteLoc;   // lat-lon of impact sites, size Ntot
	
	vector<double> D;       // distances to impact sites, size ND, units of radii (NOT CIRCUMFERENCE)
	vector<double> siteAzm; // azimuth angle of outgoing secondary debris from impact location, size Ntot (defines the centers)
	vector<double> ROIAzm;  // azimuth angle of incoming direction at ROI from impact location, size Nazm
	
	double ROI_SA;           // units of m^2, surface area of ROI
	double ROI_radius;       // units of m, radius of ROI (assuming a circular region), *will be equal to half the first site distance
	vector<double> site_SA;  // surface area of sites (units of m^2), size ND
};

/////////////////////////////


//////////////////////////////////////
// Here we will define the different output types as base classes
//  and the user class will be a template with the type of the output
//  We will first define an abstract base class which forces
//  a standard set of functions for all output types (don't forget the =0 for virtual functions!)
//////////////////////////////////////

enum xScaleType {linearScale, log10Scale};

// virtual class
class GeneralIntegralFluxOutput
{
public:
	GeneralIntegralFluxOutput(string oType,
							  double new_xMin,
							  double new_xMax,
							  int new_Nx,
							  int new_xScale,
							  int new_NSetsXY,
							  vector<double>& new_setMin,
							  vector<double>& new_setMax);

	~GeneralIntegralFluxOutput();

	// Used for all flux output types, and is the basis for igloo output
	virtual void updateFlux(double flux, double alt, double azm, double speed) = 0;
	
	virtual void saveFluxToFile(string fn) = 0;

	virtual int getNalt() = 0;
	virtual int getNazm() = 0;
	virtual int getNvel() = 0;

	void dispOutputType();
	void dispXScaleType();
	void dispSetMinMax();
	void get_xDataPointer(vector<vector<double>>* p_xData);

protected:
	string outputType; // higher gens will store their type here

	double xMin; // min of x-axis of integral flux
  	double xMax; // max of x-axis of integral flux
  	int Nx;       // number of spacings on x-axis (for igloo, will have to flatten alt-azm into x)
  	int xScale;      // xScaleType = linear or log10
  	int NSetsXY;           // number of sets of x-y data, if 0 will ignore setMin and setMax
  	vector<double> setMin; // minimum of set range i, edges can touch, but do not overlap!
  	vector<double> setMax; // maximum of set range i

  	vector<vector<double>> xData; // xData[j] = j-th setXY for all Nx
};

//////////////////////////////////////

class MassLimitedIntegralFluxVsMass : public GeneralIntegralFluxOutput
{
public:
	MassLimitedIntegralFluxVsMass(double new_xMin,
							 	  double new_xMax,
								  int new_Nx,
								  int new_xScale,
								  int new_NSetsXY,
								  vector<double>& new_setMin,
								  vector<double>& new_setMax);
	~MassLimitedIntegralFluxVsMass();

	void updateFlux(double flux, double alt, double azm, double speed);
	void saveFluxToFile(string fn);

	int getNalt();
	int getNazm();
	int getNvel();

private:
	vector<double> fraction_GT_m; // size Nx (set up during init)
};

//////////////////////////////////////

class SizeLimitedIntegralFluxVsSpeed : public GeneralIntegralFluxOutput
{
public:
	SizeLimitedIntegralFluxVsSpeed(double new_xMin,
							 	   double new_xMax,
								   int new_Nx,
								   int new_xScale,
								   int new_NSetsXY,
								   vector<double>& new_setMin,
								   vector<double>& new_setMax);
	~SizeLimitedIntegralFluxVsSpeed();

	void updateFlux(double flux, double alt, double azm, double speed);
	void saveFluxToFile(string fn);

	int getNalt();
	int getNazm();
	int getNvel();

private:
	vector<double> fraction_GT_d; // size new_NSetsXY (set up during init)
};

// //////////////////////////////////////

class MassLimitedIglooIntegratedFlux : public GeneralIntegralFluxOutput
{
public:
	MassLimitedIglooIntegratedFlux(double new_xMin, // not used here
							 	   double new_xMax, // not used here
								   int new_angleRes, // angular (degrees) resolution of theta and phi
								   int new_xScale,
								   int new_NSetsXY,
								   vector<double>& new_setMin,
								   vector<double>& new_setMax);
	~MassLimitedIglooIntegratedFlux();

	void updateFlux(double flux, double alt, double azm, double speed);
	void saveFluxToFile(string fn);

	int getNalt();
	int getNazm();
	int getNvel();

private:
	int H_getIglooNx(int aRes); // input in units of degrees
	inline int H_getIglooJN(int aRes, double azm1, double azm2);// input in units of degrees

	// new_angleRes = angular resolution of Ntheta and Nphi
	int angleRes; // should be 1, 2, 3, 4, or 5 (also 6, 9, 10, 12, 15, 18, 20, 30, 45, 60, 90)
	// columns for the igloo file, all size Nx
	vector<int> igloo_ID;
	vector<int> igloo_I;
	vector<int> igloo_J;
	vector<double> igloo_PHI1;
	vector<double> igloo_PHI2;
	vector<double> igloo_THETA1;
	vector<double> igloo_THETA2;
	vector<double> igloo_PHIavg;
	vector<double> igloo_THETAavg;
};


// //////////////////////////////////////


// Note: a template class must be definied in the header file
template <class genOutput> 
class SecondaryFluxData
{
public:
	SecondaryFluxData(string fn, // file name
				  double new_xMin, // min of x-axis of integral flux
				  double new_xMax, // max of x-axis of integral flux
				  int new_Nx,          // number of spacings on x-axis
				  int new_xScale,      // xScaleType = linear or log10
				  int new_NSetsXY,           // number of sets of x-y data, if 0 will ignore setMin and setMax
				  vector<double>& new_setMin, // minimum of set range i
				  vector<double>& new_setMax) // maximum of set range i
	{
		filename = fn;
		fluxData.resize(1);
		fluxData[0] = new genOutput(new_xMin, new_xMax, new_Nx, new_xScale, new_NSetsXY, new_setMin, new_setMax);
	}
	~SecondaryFluxData() {
		delete fluxData[0];
	}

	void updateFlux(double flux, double alt, double azm, double speed)
	{
		fluxData[0]->updateFlux(flux, alt, azm, speed);
	}

	void saveFluxToFile()
	{
		fluxData[0]->saveFluxToFile(filename);
	}

	void getOutputFilename() {
		cout << " Output filename = " << filename << endl;
	}

	int getNalt() { return fluxData[0]->getNalt(); }
	int getNazm() { return fluxData[0]->getNazm(); }
	int getNvel() { return fluxData[0]->getNvel(); }


private:
	string filename;
	vector<genOutput*> fluxData;
};

#endif 