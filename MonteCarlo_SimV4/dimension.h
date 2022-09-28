#ifndef DIMENSION_H
#define DIMENSION_H

#include "igloo.h"
#include "cdf.h"

#include <string>

using namespace std;

// virtual class, must create derived classes
class dimension
{
public:
	dimension();
	~dimension();

	// pure virtual function, derived class must implement this function
	virtual void define() = 0;
	virtual void getNextState() = 0; 
	void getCurrentState(vector<double>&); // returns dimensionStack

protected:
	string name;
	vector<double> dimensionStack;
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


class dim_timeStamp_mjd : public dimension
{
public:
	dim_timeStamp();
	~dim_timeStamp();
	
	void define();
	void getNextState();

private:
	vector<double> mjd; // modified Julian Date, gives the number of days since midnight on November 17, 1858 GMT
};


class dim_location_rll : public dimension
{
public:
	dim_location_rll();
	~dim_location_rll();

	void define();
	void getNextState();

private:
	vector<double> radius;    // Rm
	vector<double> latitude;  // degrees
	vector<double> longitude; // East, degrees
};


class dim_location_ll : public dimension
{
public:
	dim_location_ll();
	~dim_location_ll();

	void define();
	void getNextState();

private:
	vector<double> latitude;  // degrees
	vector<double> longitude; // East, degrees
};


class dim_orientation_sph : public dimension
{
public:
	dim_orientation();
	~dim_orientation();

	void define();
	void getNextState();

private:
	vector<double> sph_zenith;   // radians, from local zenith
	vector<double> sph_azimuth;  // radians, from east* (or north?)
};


class dim_orientation_vert : public dimension
{
public:
	dim_orientation();
	~dim_orientation();

	void define();
	void getNextState();
};


class dim_geometry_cyl : public dimension
{
public:
	dim_geometry_cyl();
	~dim_geometry_cyl();

	void define();
	void getNextState();
};


class dim_density_const : public dimension
{
public:
	dim_density_const();
	~dim_density_const();

	void define();
	void getNextState();

private:
	double density; // kg/m^3, bulk density
};

class dim_density_logNorm : public dimension
{
public:
	dim_density_const();
	~dim_density_const();

	void define();
	void getNextState();

private:
	double norm;  // normalization factor
	double sigma; // natural log variance
	double mu;    // natural log mean

	double logNormPDF(double); // returns fraction, input density [kg/m^3]
	double logNormCDF(double); // returns cumulative fraction, input density [kg/m^3]
};


class dim_porosity_const : public dimension
{
public:
	dim_porosity_const();
	~dim_porosity_const();

	void define();
	void getNextState();

private:
	double porosity; // fraction
};

class dim_strength_const : public dimension
{
public:
	dim_strength_const();
	~dim_strength_const();

	void define();
	void getNextState();

private:
	double strength; // Pa = kg/m/s^2
};

class dim_scalingLawParams : public dimension
{
public:
	dim_scalingLawParams();
	~dim_scalingLawParams();

	void define();
	void getNextState();

private:
	// all the scaling law parameters go here
};


class dim_mass_func : public dimension
{
public:
	dim_mass_func();
	~dim_mass_func();

	void define();
	void getNextState();

private:
	double (*massCDF_func)(double ); // return #/m^2/s, input mass kg
}



class dim_fluxSpeedAngle_igloo : public dimension
{
public:
	dim_fluxSpeedAngle_igloo();
	~dim_fluxSpeedAngle_igloo();

	void define();
	void getNextState();

private:
	vector<igloo*> fluxPDF; // as a function of location on target
	vector<cdf*>   fluxCDF; // derived from igloo pdf
};

class dim_fluxSpeedAngleDistance : public dimension
{
public:
	dim_fluxSpeedAngleDistance();
	~dim_fluxSpeedAngleDistance();

	void define();
	void getNextState();

private:
	double currentDistance; // rm, regenerate speedAngleCDF if distance is dfferent from currentDistance
	cdf*   speedAngleCDF;   // for a given distance
};



class dim_fluxSpeedAngleMass : public dimension
{
public:
	dim_fluxSpeedAngleMass();
	~dim_fluxSpeedAngleMass();

	void define();
	void getNextState();

private:
	vector<igloo*> fluxPDF; // as a function of mass
	vector<double> mass;
};

#endif 