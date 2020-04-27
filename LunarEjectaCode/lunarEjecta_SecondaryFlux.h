#ifndef LUNAREJECTA_SECONDARYFLUX_H
#define LUNAREJECTA_SECONDARYFLUX_H

#include <iostream>
#include <vector>

using namespace std;


class latLon
{
public:
	latLon(double new_lat, double new_lon); // input in degrees
	~latLon();

	double getLatRad();
	double getLonRad();

	double getLatDeg();
	double getLonDeg();

	void dispLatLon();

	double getNormDistTo(latLon& toLocation);
	double getTanDistTo(latLon& toLocation);
	double getBearingTo(latLon& toLocation);
	double getBearingFrom(latLon& fromLocation);
	
private:
	double lat; // radians
	double lon; // radians
};


class ImpactSites_and_ROI
{
public:
	ImpactSites_and_ROI();
	~ImpactSites_and_ROI();
	
private:
	vector<latLon*> sites;
};

#endif 