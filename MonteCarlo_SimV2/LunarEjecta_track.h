#ifndef LUNAREJECTA_TRACK_H
#define LUNAREJECTA_TRACK_H

#include "LunarEjecta_params.h"
#include "LunarEjecta_asset.h"


using namespace std;

const double EPS = 1.E-4;

struct trackVars
{
	double x0;
	double y0;
	double z0;
	double u0;
	double v0;
	double w0;

	double xi;
	double yi;
	double zi;
	double ui;
	double vi;
	double wi;

	// double dpx;
	// double dpy;
	// double Fx;
	// double Fy;

	double h; // timestep
	double dxmin; // minimum distance
};

void init_trackVars(trackVars& tv, double x0, double y0, double z0, double u0, double v0, double w0, double h, double dx0);

void unpackFinalLoc(trackVars& track_i, vector<double> &loc_f);
void unpackFinalPh(trackVars& track_i, vector<double> &ph_f);


// doesn't need to be initialized. Overridden during runtime
struct RK45VarsPosVel
{
	double kx[6];
	double ky[6];
	double kz[6];
	double ku[6];
	double kv[6];
	double kw[6];

	double xi1[2];
	double yi1[2];
	double zi1[2];
	double ui1[2];
	double vi1[2];
	double wi1[2];

	double Rx;
	double Ry;
	double Rz;
	double Ru;
	double Rv;
	double Rw;

	double deltaX;
	double deltaY;
	double deltaZ;
	double deltaU;
	double deltaV;
	double deltaW;
};


const double RK45Coeff[8][6] =
	{ // k1         , k2          , k3          , k4           , k5      , k6
		{0.         , 0.          , 0.          , 0.           ,       0., 0.    },
		{1./4.      , 0.          , 0.          , 0.           ,       0., 0.    },
		{3./32.     , 9./32.      , 0.          , 0.           ,       0., 0.    },
		{1932./2197., -7200./2197., 7296./2197. , 0.           ,       0., 0.    },
		{439./216.  , -8.         , 3680./513.  , -845./4104.  ,       0., 0.    },
		{-8./27.    , 2.          , -3544./2565., 1859./4104.  , -11./40., 0.    },
		{25./216.   , 0.          , 1408./2565. , 2197./4104.  , -1./5.  , 0.    },
		{16./135.   , 0.          , 6656./12825., 28561./56430., -9./50. , 2./55.}
	};

// Gravitational force function
void grav_a(double& ax, double& ay, double& az, double x, double y, double z);

// Used to update the track by one time step
void RK45UpdatePosVel(trackVars& t,
	                  RK45VarsPosVel& r,
	                  const double c[][6], asset& ag,
	                  double Rm,
	                  double gm);

//void printTrack(ofstream& file, trackVars& tv, sums& s);

// After each step update, a collision check is made
// First, we check if the track is within a certain radius,
// if so, then we further check each asset shape against the track location
bool checkCollision(trackVars& t, asset& ag);


// loc = physical location in global cartesian coordinates
// ph  = location of initial ejecta conditions in phase space (i.e., ejecta speed, zenith angle, and azimuth angle)
// 
bool runTraj_checkHit(vector<double> &lat_lon, vector<double> &loc, vector<double> &ph, vector<double> &loc_f, vector<double> &ph_f, asset &ag, double Rm, double vesc, double gm);

#endif 