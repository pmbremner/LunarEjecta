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

	inline void setLatRad(double new_lat);
	inline void setLonRad(double new_lon);
	inline void setLatDeg(double new_lat);
	inline void setLonDeg(double new_lon);

	double getLatRad();
	double getLonRad();

	double getLatDeg();
	double getLonDeg();

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
	ImpactSites_and_ROI();
	~ImpactSites_and_ROI();
	
private:
	vector<latLon*> sites;
};

#endif 