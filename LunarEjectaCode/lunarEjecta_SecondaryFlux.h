#ifndef LUNAREJECTA_SECONDARYFLUX_H
#define LUNAREJECTA_SECONDARYFLUX_H

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
	ImpactSites_and_ROI(double new_ND,
		                double new_Nazm,
		                double new_radius,
		                latLon& new_ROI);
	~ImpactSites_and_ROI();
	
private:
	inline int H_idx(int i_azm, int j_dist);

	int ND;
	int Nazm;
	int Ntot;

	double radius; // units of m
	
	latLon ROI;              // region of interest location
	vector<latLon> siteLoc;   // lat-lon of impact sites, size Ntot
	
	vector<double> D;       // distances to impact sites, size ND, units of radii
	vector<double> siteAzm; // azimuth angle of outgoing secondary debris from impact location, size Ntot (defines the centers)
	vector<double> ROIAzm;  // azimuth angle of incoming direction at ROI from impact location, size Nazm
	
	double ROI_SA;           // surface area of ROI
	double ROI_radius;       // radius of ROI (assuming a circular region), *will be equal to half the first site distance
	vector<double> site_SA;  // surface area of sites, size ND
};

/////////////////////////////


//////////////////////////////////////
// Here we will define the different output types as base classes
//  and the user class will be a template with the type of the output
//  We will first define an abstract base class which forces
//  a standard set of functions for all output types (don't forget the =0 for virtual functions!)
//////////////////////////////////////

enum xScaleType {linearScale, log10Scale};


class GeneralIntegralFluxOutput
{
public:
	GeneralIntegralFluxOutput(string oType,
							  double new_xMin,
							  double new_xMax,
							  int new_Nx,
							  int new_xScale,
							  int new_NSetsXY,
							  vector<double> new_setMin,
							  vector<double> new_setMax);

	~GeneralIntegralFluxOutput();

	virtual void computeSecondaryFlux(double D0,   double D1,
		                              double azm0, double azm1) = 0;
	virtual void saveFluxToFile(string fn) = 0;
	void dispOutputType();

protected:
	virtual void H_updateFlux(double flux, double alt, double azm, double speed) = 0;
	string outputType; // higher gens will store their type here

	double xMin; // min of x-axis of integral flux
  	double xMax; // max of x-axis of integral flux
  	int Nx;       // number of spacings on x-axis
  	int xScale;      // xScaleType = linear or log10
  	int NSetsXY;           // number of sets of x-y data, if 0 will ignore setMin and setMax
  	vector<double> setMin; // minimum of set range i
  	vector<double> setMax; // maximum of set range i
};

//////////////////////////////////////

class MassLimitedIntegralFluxVsMass : public GeneralIntegralFluxOutput
{
public:
	MassLimitedIntegralFluxVsMass(string oType,
							  	  double new_xMin,
							 	  double new_xMax,
								  int new_Nx,
								  int new_xScale,
								  int new_NSetsXY,
								  vector<double> new_setMin,
								  vector<double> new_setMax);
	~MassLimitedIntegralFluxVsMass();

	void computeSecondaryFlux(double D0,   double D1,
		                      double azm0, double azm1);
	void saveFluxToFile(string fn);
private:
	void H_updateFlux(double flux, double alt, double azm, double speed);
	
};

// //////////////////////////////////////

// class SizeLimitedIntegralFluxVsSpeed : public GeneralIntegralFluxOutput
// {
// public:
// 	SizeLimitedIntegralFluxVsSpeed();
// 	~SizeLimitedIntegralFluxVsSpeed();

// 	void computeSecondaryFlux(double D0,   double D1,
// 		                           double azm0, double azm1);
// 	void saveFluxToFile(string fn);
// private:
// 	void H_updateFlux(double flux, double alt, double azm, double speed);

// };

// //////////////////////////////////////

// class MassLimitedIglooIntegratedFlux : public GeneralIntegralFluxOutput
// {
// public:
// 	MassLimitedIglooIntegratedFlux();
// 	~MassLimitedIglooIntegratedFlux();

// 	void computeSecondaryFlux(double D0,   double D1,
// 		                           double azm0, double azm1);
// 	void saveFluxToFile(string fn);
// private:
// 	void H_updateFlux(double flux, double alt, double azm, double speed);

// };


// //////////////////////////////////////


// // Note: a template class must be definied in the header file
// template <class genOutput> 
// class SecondaryFlux
// {
// public:
// 	SecondaryFlux(double new_xMin, // min of x-axis of integral flux
// 				  double new_xMax, // max of x-axis of integral flux
// 				  int Nx,          // number of spacings on x-axis
// 				  int xScale,      // xScaleType = linear or log10
// 				  int NSetsXY,           // number of sets of x-y data, if 0 will ignore setMin and setMax
// 				  vector<double> setMin, // minimum of set range i
// 				  vector<double> setMax) // maximum of set range i
// 	~SecondaryFlux();
// private:
// 	string directoryName;
// 	vector<genOutput*> fluxData;
// };

#endif 