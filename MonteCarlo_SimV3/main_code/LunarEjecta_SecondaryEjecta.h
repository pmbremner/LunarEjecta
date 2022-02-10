#ifndef LUNAREJECTA_SECONDARYEJECTA_H
#define LUNAREJECTA_SECONDARYEJECTA_H

#include "LunarEjecta_params.h"

#include <vector>
#include <random>     // mt19937


using namespace std;


double Fspeed(double g, double h_rm, double d_rm);
double Fspeed_v(double v, vector<double>& vars);
double Fspeed_g(double g, vector<double>& vars);


double azm_FOV(double R, double D);

// returns azm1, either forward or reverse, depending on short_far_flag
double azm_bearing(double lat1, double lon1, double lat2, double lon2, bool short_far_flag);

// returns distance, d
double lat_lon_dist(double lat1, double lon1, double lat2, double lon2, bool short_far_flag);

// returns lat2
double destination_lat(double lat1, double d, double azm);

// returns lon2, needs lat2
double destination_lon(double lat1, double lon1, double lat2, double d, double azm);


// final speed at asset
/// a is the asset altitude in units of rm, vp is the ejecta speed at primary impact in units of escape speed
/// returns ejecta speed at asset in units of escape speed
double final_speed(double a, double vp);

// final zenith at asset (as seen from the asset)
/// d is in units of rm, vp is in units of vesc, g is zenith at ejected point in rads
/// return ejecta zenith at asset in units of radians
double final_zenith(double d, double vp, double g);


// the smallest zenith angle to reach asset at ~ escape speed
/// a is the asset altitude in units of rm, d_rm is the projected distance from the impact point to asset in units of rm
/// return ejecta zenith, at ejecta point, in radians
double min_zenith_at_escape(double a, double d_rm);



// if |f(g+dg) - f(g)| > dv, then find dg such that |f(g+dg) - f(g)| = dv, else keep dg
// i.e., sample the zenith angle finer for when the function is steeper
// Note: make sure vars[1] = h and vars[2] = d before this function call, both in units of lunar radii
double find_dg(double g, double dg, double dv, double vmax, vector<double>& vars);

// For given g, h, d, r, compute the dg which is the minimum dg of the four corners of the wedge
double find_dg_wedge(double g, double dg, double dv, double vlow, double vmax, vector<double>& vars, double a, double h, double d, double r);

void find_min_max_v(double g, double& vmin, double& vmax, double vlow, double vlim, vector<double>& vars, double a, double h, double d, double r);


// zenith in units of rads in (0,Pi/2]
// vmin and vmax in units of vesc in (0, vlim]
// h is the height of the wedge
// d is the distance of the wedge, center
// r is the radius of the wedge (front and back)
// the wedge enscribes a sphere-like shape (h can be controlled, so really it's an ellipsoid)
// dg and dv are the maximum grid spacing for the zenith and speed dimensions, respectively
void get_zenith_speed_grid(vector<double>& zenith, vector<double>& vmin, vector<double>& vmax, double vlow, double vlim, double a, double h, double d, double r, double dg, double dv);



// Computes the CDF by the trapezoidal area between each zenith grid point
void get_CDF_PDF_from_trapdens(vector<double>& zenith, vector<double>& vminv, vector<double>& vmaxv, vector<double>& cdf, vector<double>& pdf);

// return the left index
int pdf_sample(mt19937& rng, vector<double>& zenith, vector<double>& vminv, vector<double>& vmaxv, vector<double>& cdf, double& g_sample, double& v_sample, double& weight);


// Using Stratified Importance Sampling
// For a given azimuth and lat-lon location
void get_samples(vector<double>& zenith, vector<double>& vminv, vector<double>& vmaxv, double vlow, double vmax, vector<double>& cdf, vector<double>& sample_zenith, vector<double>& sample_speed, vector<double>& sample_weight, int N_sample);

// For a given lat-lon, compute samples both in close and far directions in azimuth, assuming a wedge
//void get_samples_with_azm(double latp, double lonp, double lats, double lons, double a, double h, double r);
void get_samples_with_azm_lat_lon(input* p,
                                  double latp,   // primary latitude center
	                              double lonp,   // primary longitude center
	                              double dlatp,  // primary latitude range
	                              double dlonp,  // primary longitude range
	                              double lats,   // satellite (asset) latitude center
	                              double lons,   // satellite (asset) longitude center
	                              double a,      // satellite (asset) altitude
	                              double h,      // satellite (asset) height
	                              double r,      // satellite (asset) radius
	                              double vmin,   // minimum ejecta speed
	                              double vmax,   // maximum ejecta speed
	                              double dg,     // maximum zenith grid width
	                              double dv,     // maximum speed grid width
	                              vector<double>& sample_latp,       // primary latitude center
	                              vector<double>& sample_lonp,       // primary longitude center
	                              vector<double>& sample_azimuth_0,  // initial azimuth, zenith, and speed of secondary, at primary impact
	                              vector<double>& sample_zenith_0,
	                              vector<double>& sample_speed_0,
	                              vector<double>& sample_azimuth_f,  // final azimuth, zenith, and speed of secondary, at asset impact
	                              vector<double>& sample_zenith_f,
	                              vector<double>& sample_speed_f,
	                              vector<double>& sample_weight,
	                              int N_azm_lat_lon,   // number of pulls in azimuth-lat-lon sets
	                              int N_zenith_speed); // number of pulls in zenith-speed sets




#endif 