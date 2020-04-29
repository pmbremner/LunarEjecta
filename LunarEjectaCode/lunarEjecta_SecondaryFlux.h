#ifndef LUNAREJECTA_SECONDARYFLUX_H
#define LUNAREJECTA_SECONDARYFLUX_H

#include <iostream>
#include <vector>

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
	double getBearingInitial(latLon& pos);
	double getBearingFinal(latLon& pos);

	// outputs in E CCW angle
	double getAzmInitial(latLon& pos);
	double getAzmFinal(latLon& pos);
	
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
	int ND;
	int Nazm;
	int Ntot;

	double radius; // units of km
	
	latLon ROI;              // region of interest location
	vector<latLon> siteLoc;   // lat-lon of impact sites, size Ntot
	
	vector<double> D;       // distances to impact sites, size ND, units of radii
	vector<double> siteAzm; // azimuth angle of outgoing secondary debris from impact location, size Nazm
	vector<double> ROIAzm;  // azimuth angle of incoming direction at ROI from impact location, size Ntot
	
	double ROI_SA;           // surface area of ROI
	double ROI_radius;       // radius of ROI (assuming a circular region), *will be equal to half the first site distance
	vector<double> site_SA;  // surface area of sites, size Ntot
};

#endif 