#include "LunarEjecta_MainUtils.h"
#include "LunarEjecta_SecondaryEjecta.h"


#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <cmath>
#include <sstream>
#include <chrono>
#include <random>     // mt19937
#include <algorithm>  // lower_bound

using namespace std;


// double Fspeed(double g, double h_rm, double d_rm)
// {
// 	double x = 1. - (h_rm / (1. + h_rm)) / (1. - cos(d_rm)); // for h = 0, x = 1

// 	return 1. / (x * (1. - cos(2.*g)) + sin(2.*g)/tan(d_rm/2.));
// }

double Fspeed(double g, double a_rm, double d_rm)
{
	// double x = (1./a_rm - cos(d_rm)) / (1. - cos(d_rm)); // for a = 1, x = 1

	// return 1. / (x * (1. - cos(2.*g)) + sin(2.*g)/tan(d_rm/2.));

	//double x =  1. - (1. - 1./a_rm) / (2.*sqr(sin(d_rm/2.)));
	double x =  (1./a_rm - cos(d_rm)) / (2.*sqr(sin(d_rm/2.)));

	return 1. / (2.*x*sqr(sin(g)) + sin(2.*g)/tan(d_rm/2.));
}

double Fspeed_direct(double g, double a_rm, double d_rm, double vmax)
{
	double v_sqr = Fspeed(g, a_rm, d_rm);

	return (v_sqr > 0. ? sqrt(v_sqr) : vmax);
}


double Fspeed_v(double v, vector<double>& vars)
{
	return Fspeed(vars[0], vars[1], vars[2]) - sqr(v);

}

double Fspeed_g(double g, vector<double>& vars)
{
	return Fspeed(g, vars[1], vars[2]) - sqr(vars[0]);
}





// R and D in units of rm
double azm_FOV(double R, double D)
{
	double denom = sqr(sin(D+R)) - sqr(sin(R));

	if (denom >= 0.)
		return 2. * atan2(sin(R), sqrt(denom));
	else
		return PI;
}


double azm_bearing(double lat1, double lon1, double lat2, double lon2, bool short_far_flag)
{
	double phase_shift = (short_far_flag == 0 ? 0. : PI);

	return atan2(sin(lon1 - lon2) * cos(lat1 - lat2) , cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon1 - lon2)) + phase_shift;
}


// return distance in units of rm
double lat_lon_dist(double lat1, double lon1, double lat2, double lon2, bool short_far_flag)
{
	double a = sqr( sin((lat1 - lat2) / 2.) ) + cos(lat1) * cos(lat2) * sqr( sin((lon1 - lon2) / 2.) );
	double pm = (short_far_flag == 0 ? 1. : -1.);

	return 2. * atan2(pm * sqrt(a), sqrt(1. - a));
}


// returns lat2
// d in units of rm
double destination_lat(double lat1, double d, double azm)
{
	return asin(sin(lat1) * cos(d) + cos(lat1) * sin(d) * cos(azm));
}

// returns lon2, needs lat2
// d in units of rm
double destination_lon(double lat1, double lon1, double lat2, double d, double azm)
{
	return lon1 + atan2(sin(azm) * sin(d) * cos(lat1), cos(d) - sin(lat1) * sin(lat2));
}


double gp_ap(double vp, double rs)
{
	return acos(sqrt((rs - 1.) * (rs/sqr(vp) - (rs + 1.))));
}


double F(double d, double vp, double g, double rs)
{
	double sign;

	if (d < PI && g > gp_ap(vp, rs)){
		sign = -1.;
	} else {
		sign = 1.;
	}

	return sign * sqrt( 1. + ((rs-1.)/sqr(cos(g))) * (rs + 1. - rs/sqr(vp)) );
}

// final speed at asset
/// rs is the asset altitude in units of rm, vp is the ejecta speed at primary impact in units of escape speed
/// returns ejecta speed at asset in units of escape speed
double final_speed(double rs, double vp)
{
	return sqrt(1./rs + sqr(vp) - 1.);
}



// final zenith at asset (as seen from the asset)
/// d is in units of rm, vp is in units of vesc, g is zenith at ejected point in rads
/// return ejecta zenith at asset in units of radians
double final_zenith(double d, double vp, double g, double rs)
{
	return atan2(tan(g), -F(d, vp, g, rs));
	//return -atan2((sqr(vp) * (cos(d-2.*g) - cos(d)) + cos(d) - 1.), (sqr(vp) * (sin(d-2.*g) - sin(d)) + sin(d))); wrong
}


// final altitude at asset
/// d is in units of rm, vp is in units of vesc, g is zenith at ejected point in rads
// return final altitude at asset in units of rm
double final_altitude(double d, double vp, double g)
{
	return 2.*sqr(vp * sin(g)) / (1. + (sqr(vp) - 1.)*cos(d) - sqr(vp)*cos(d - 2.*g));
}



double distance(double vp, double g, double rs, double d_side)
{
	double Fi = F(d_side, vp, g, rs);
	return fmod(2.*PI + 2. * atan2( sqr(vp) * sin(2.*g) + ((rs - 1.) / (1. - Fi)) * (2.*sqr(vp) - 1.) * tan(g),
		               ((rs - Fi) / (1. - Fi)) - 2.*sqr(vp * sin(g)) ), 2.*PI);
}

// this is half the distance from crater to reimpacting the Moon
// vp is in units of vesc, g is zenith at ejected point in rads
// returns distance in units of rm
double apoapsis_distance(double vp, double g, double rs)
{

	return atan2( sqr(vp) * sin(2.*g), 1. - 2.*sqr(vp * sin(g)) ); 
}



// the smallest zenith angle to reach asset at ~ escape speed
/// rs is the asset altitude in units of rm, d_rm is the projected distance from the impact point to asset in units of rm
/// return ejecta zenith, at ejecta point, in radians
double min_zenith_at_escape(double rs, double d_rm)
{
	return atan2( sin(d_rm/2.), cos(d_rm/2.) + sqrt(1./rs));
}





// if |f(g+dg) - f(g)| > dv, then find dg such that |f(g+dg) - f(g)| = dv, else keep dg
// i.e., sample the zenith angle finer for when the function is steeper
// Note: make sure vars[1] = h and vars[2] = d before this function call, both in units of lunar radii
double find_dg(double g, double dg, double dv, double vmax, vector<double>& vars)
{
	const double alpha = 0.2;
	const double eps   = 1.E-3;

	double f0, f1;
	double dg_1 = dg;

	do
	{
		// compute the function at g and g+dg
		vars[0] = g;

		//f0 = findX(0., Fspeed_v, 0., vmax, vars);
		f0 = Fspeed_direct(vars[0], vars[1], vars[2], vmax);

		vars[0] = g + dg_1;
		//f1 = findX(0., Fspeed_v, 0., vmax, vars);
		f1 = Fspeed_direct(vars[0], vars[1], vars[2], vmax);

		// compute a weighted average dg so as not to overshoot intended result
		// (due to taking a linear approximation of the derivative)
		// alpha ranges from 0 to 1, aim for alpha ~ 0.3, possibly higher
		dg_1 *= dv / fabs(f1 - f0) * alpha + (1. - alpha);

		// if the proposed dg is too large, keep the original dg
		dg_1 = (dg_1 < dg ? dg_1 : dg);

		vars[0] = g + dg_1;
		//f1 = findX(0., Fspeed_v, 0., vmax, vars);
		f1 = Fspeed_direct(vars[0], vars[1], vars[2], vmax);

	// do while the proposed dv is too large
	} while (fabs(f1 - f0) > dv*(1. + eps));

	return dg_1;
}

// For given g, h, d, r, compute the dg which is the minimum dg of the four corners of the wedge
double find_dg_wedge(double g, double dg, double dv, double vlow, double vmax, vector<double>& vars, double a, double h, double d, double r)
{
	int i, j;
	double dg_min = dg;

	for (i = 0; i < 2; ++i)
		for (j = 0; j < 2; ++j)
		{
			vars[0] = g; // zenith angle
			vars[1] = a + h * j; // height
			//vars[2] = d + r*(2.*i - 1); // width
			vars[2] = d + r*2.*i; // width

			// ignore contribution to speed-zenith curves that go above vmax
			//if(findX(0., Fspeed_v, vlow, vmax, vars) < vmax*0.9999)
			if(Fspeed_direct(vars[0], vars[1], vars[2], vmax) < vmax*0.9999)
				dg_min = min(find_dg(g, dg, dv, vmax, vars), dg_min);
		}

	return dg_min;
}

void find_min_max_v(double g, double& vmin, double& vmax, double vlow, double vlim, vector<double>& vars, double a, double h, double d, double r)
{
	int i, j;
	double f;

	vmin = vlim;
	vmax = 0.;

	vars[0] = g;

	for (i = 0; i < 2; ++i)
		for (j = 0; j < 2; ++j)
		{
			vars[1] = a + h * j; // height
			//vars[2] = d + r*(2.*i - 1); // width
			vars[2] = d + r*2.*i; // width

			//f = findX(0., Fspeed_v, vlow, vlim, vars);
			f = Fspeed_direct(vars[0], vars[1], vars[2], vmax);

			vmin = min(f, vmin);
			vmax = max(f, vmax);
		}
	// force vmin to be at least vlow, and vmax to be at most vlim
	vmin = max(vlow, vmin);
	vmax = max(min(vlim, vmax), vlow);

	// if both vmin and vmax are very close to the maximum speed, vlim, then force them to be to avoid sampling noisy unphysical values
	if (fabs(vlim - vmin) < 1E-2 && fabs(vlim - vmax)  < 1E-2)
	{
		vmin = vlim;
		vmax = vlim;
	}
}


// zenith in units of rads in (0,Pi/2]
// vmin and vmax in units of vesc in (0, vlim]
// h is the height of the wedge
// d is the distance of the wedge, front edge
// r is the radius of the wedge (front and back)
// the wedge enscribes a sphere-like shape (h can be controlled, so really it's an ellipsoid)
// dg and dv are the maximum grid spacing for the zenith and speed dimensions, respectively
void get_zenith_speed_grid(vector<double>& zenith, vector<double>& vmin, vector<double>& vmax, double vlow, double vlim, double a, double h, double d, double r, double dg, double dv)
{
	zenith.clear();
	vmin.clear();
	vmax.clear();

	vector<double> vars(3, 0.); // size of 3, filled with zeros

	// First, find the smallest zenith angle at the closest point
	// vars[0] = vlim;
	// vars[1] = 0.;
	// //vars[2] = d-r;
	// vars[2] = d;

	// Note: need to go slightly higher than the bound in order to not get stuck
	//double g_min = 1.01*(d-r)/4.;//findX(0., Fspeed_g, 0.000001, PI/2., vars);
	//double g_min =  d/4.; // works, need to do better ?
	double g_min = min_zenith_at_escape(a, d);


	//cout << "g_min = " << g_min << endl;

	// compute the grid
	double g_cur, v0, v1, dg_new;

	// the first point
	g_cur = g_min;
	find_min_max_v(g_cur, v0, v1, vlow, vlim, vars, a, h, d, r);
	///cout << g_cur << ' ' << v0 << ' ' << v1 << ' ' << v1-v0 << endl;

	zenith.emplace_back(g_cur);
	vmin.emplace_back(v0);
	vmax.emplace_back(v1);

	// the rest of the points
	while (g_cur < PI/2.)
	{
		// First, find the dg amount, limited by each corner
		vars[1] = h;
		vars[2] = d;
		dg_new = find_dg_wedge(g_cur, dg, dv, vlow, vlim, vars, a, h, d, r);

		g_cur += dg_new;
		g_cur = (g_cur > PI/2. ? PI/2. : g_cur); // for last grid point

		// Next, find the minimum and maximum speeds of the corners
		find_min_max_v(g_cur, v0, v1, vlow, vlim, vars, a, h, d, r);

		///cout << g_cur << ' ' << v0 << ' ' << v1 << ' ' << v1-v0 << endl;

		zenith.emplace_back(g_cur);
		vmin.emplace_back(v0);
		vmax.emplace_back(v1);
	}
	///cout << endl << endl;

}



// Computes the CDF by the trapezoidal area between each zenith grid point
void get_CDF_PDF_from_trapdens(vector<double>& zenith, vector<double>& vminv, vector<double>& vmaxv, vector<double>& cdf, vector<double>& pdf)
{
	int N = zenith.size(), i;
	double sum;

	cdf.clear();
	cdf.resize(N);

	pdf.clear();
	pdf.resize(N-1);

	// compute the trapezoidal area of each section
	// So that the density of sampled points is the same no matter the size of the region
	for (i = 0; i < N-1; ++i)
		pdf[i] = trapezoid_area(zenith[i+1] - zenith[i], vmaxv[i] - vminv[i], vmaxv[i+1] - vminv[i+1]);//0.5 * (zenith[i+1] - zenith[i]) * (vmaxv[i] - vminv[i] + vmaxv[i+1] - vminv[i+1]);

	// normalize density distribution (the PDF)
	sum = vSum(pdf);
	for (i = 0; i < N-1; ++i)
		pdf[i] /= sum;

	// cout << "norm sum = " << vSum(pdf) << endl;

	// cout << endl;
	// for (int i = 0; i < pdf.size(); ++i)
	// 	cout << zenith[i] << ' ' << pdf[i] << endl;

	// compute the CFD, starting from 0 to 1
	vCumLow(pdf, cdf);

	// cout << endl;
	// for (int i = 0; i < cdf.size(); ++i)
	// 	cout << zenith[i] << ' ' << cdf[i] << endl;
}

// return the left index
int pdf_sample(mt19937& rng, vector<double>& zenith, vector<double>& vminv, vector<double>& vmaxv, vector<double>& cdf, double& g_sample, double& v_sample, double& weight)
{
	//// Find the zenith sample ///////
	double u;
	int idx = sample_pdf_idx(rng, cdf, u);

	// cout << cdf.size() << " | ";
	// cout << idx << ' ' << (*(idx_iter-1)) << ' ' << u << ' ' << (*(idx_iter)) << endl;

	// linear interpolate to get exact sample points
	g_sample = ((zenith[idx] - zenith[idx-1]) / (cdf[idx] - cdf[idx-1])) * (u - cdf[idx-1]) + zenith[idx-1];
	//cout << idx << ' ' << cdf[idx-1] << ' ' << u << ' ' << cdf[idx] << " | " << zenith[idx-1] << ' ' << g_sample << ' ' << zenith[idx] << endl;
	//cout << u << ' ' << g_sample << endl;

	//// Find the speed sample ///////
	double vlow  = ((vminv[idx] - vminv[idx-1]) / (cdf[idx] - cdf[idx-1])) * (u - cdf[idx-1]) + vminv[idx-1];
	double vhigh = ((vmaxv[idx] - vmaxv[idx-1]) / (cdf[idx] - cdf[idx-1])) * (u - cdf[idx-1]) + vmaxv[idx-1];

	v_sample = uniform(rng, vlow, vhigh);

	// associated with the area of the region (i.e., the normalization factor)
	weight = trapezoid_area(zenith[idx] - zenith[idx-1], vmaxv[idx-1] - vminv[idx-1], vmaxv[idx] - vminv[idx]);

	return idx-1;
}



// Using Stratified Importance Sampling
void get_samples(vector<double>& zenith,
	             vector<double>& vminv,
	             vector<double>& vmaxv,
	             double vlow,
	             double vmax,
	             vector<double>& cdf,
	             double d,
	             double a,
	             double h,
	             double r,
	             vector<double>& sample_zenith,
	             vector<double>& sample_speed,
	             vector<double>& sample_zenith_f,
	             vector<double>& sample_speed_f,
	             vector<double>& sample_altitude_f,
	             vector<double>& sample_distance_f,
	             vector<double>& sample_weight,
	             int N_sample)
{
	//init the random generator
	// https://stackoverflow.com/questions/24334012/best-way-to-seed-mt19937-64-for-monte-carlo-simulations
	mt19937 rng;
	unsigned seed; 

	seed = random_device{}() * chrono::system_clock::now().time_since_epoch().count();
	rng.seed(seed);

	sample_zenith.clear();
	sample_speed.clear();
	sample_zenith_f.clear();
	sample_speed_f.clear();
	sample_altitude_f.clear();
	sample_distance_f.clear();
	sample_weight.clear();

	sample_zenith.resize(N_sample);
	sample_speed.resize(N_sample);
	sample_zenith_f.resize(N_sample);
	sample_speed_f.resize(N_sample);
	sample_altitude_f.resize(N_sample);
	sample_distance_f.resize(N_sample);
	sample_weight.resize(N_sample);

	double r_s;//, d_s; // final altitute and distance


	//init sample_idx, don't need it outside of get_samples function
	vector<int> sample_idx(N_sample, 0);

	// size of the pdf array
	vector<int> stratified_count(cdf.size()-1, 0);

	for (int i = 0; i < N_sample; ++i)
	{
		// pull sample from pdf
		do{
			sample_idx[i] = pdf_sample(rng, zenith, vminv, vmaxv, cdf, sample_zenith[i], sample_speed[i], sample_weight[i]);

			// check altitude of impact at asset at the front face
			// d is the distance from crater to front face of asset (not center of asset)
			r_s = final_altitude(d, sample_speed[i], sample_zenith[i]);


			// Hit the front face (side of "cylinder")
			if (r_s >= a && r_s <= a + h)
			{
				sample_altitude_f[i] = r_s;
				sample_distance_f[i] = d;
			}
			else // top or bottom side
			{

				if(r_s > a + h) // top face
				{
					sample_altitude_f[i] = a + h;
					sample_distance_f[i] = distance(sample_speed[i], sample_zenith[i], a + h, d);

				}
				else // bottom face, r_s < a
				{
					sample_altitude_f[i] = a;
					sample_distance_f[i] = distance(sample_speed[i], sample_zenith[i], a, d);
				}
			}



			// If not at the front face, check top face
			// If not at the top face, check bottom face
			// If no face found, then redo the sample

		} while(!( (sample_distance_f[i] >= d && sample_distance_f[i] <= d + 2.*r) && (sample_altitude_f[i] >= a && sample_altitude_f[i] <= a + h) ));

		stratified_count[sample_idx[i]]++;
	}

	// ofstream file;
	// file.open("samples.txt");

	// after finishing all the samples, we need to divide the weights by the number of counts in each region, NOT the total number of counts overall
	// to normalize the total speed-zenith area to 1, need to divide by the lengths in both dimensions
	for (int i = 0; i < N_sample; ++i){
		// divide by the overall total area, to normalize available area to 1
		// so that later when the ejecta distribution function is taken into account,
		// the area element comes into play there, Or just take into account here, and normalize the ejecta distribution
		sample_weight[i] /= double(stratified_count[sample_idx[i]]);// * PI/2. * (vmax - vlow); 

		//file << sample_zenith[i] << ' ' << sample_speed[i] << ' ' << sample_weight[i] << endl;
	}
	//file.close();
}





void get_samples_with_azm_lat_lon(input* p,
                                  double latp,   // primary latitude center [radians]
	                              double lonp,   // primary longitude center [radians]
	                              double dlatp,  // primary latitude range [radians]
	                              double dlonp,  // primary longitude range [radians]
	                              double lats,   // satellite (asset) latitude center [radians]
	                              double lons,   // satellite (asset) longitude center [radians]
	                              double a,      // satellite (asset) altitude [rm]
	                              double h,      // satellite (asset) height [rm]
	                              double r,      // satellite (asset) radius [rm]
	                              double vmin,   // minimum ejecta speed [vesc]
	                              double vmax,   // maximum ejecta speed [vesc]
	                              double dg,     // maximum zenith grid width [radians]
	                              double dv,     // maximum speed grid width [vesc]
	                              vector<double>& sample_latp,       // primary latitude center
	                              vector<double>& sample_lonp,       // primary longitude center
	                              vector<double>& sample_azimuth_0,  // initial azimuth, zenith, and speed of secondary, at primary impact
	                              vector<double>& sample_zenith_0,
	                              vector<double>& sample_speed_0,
	                              vector<double>& sample_azimuth_f,  // final azimuth, zenith, and speed of secondary, at asset impact
	                              vector<double>& sample_zenith_f,
	                              vector<double>& sample_speed_f,
	                              vector<double>& sample_altitude_f, // final altitude from center of Moon [rm]
	                              vector<double>& sample_distance_f, // final distance from crater to impact at asset [rm]
	                              vector<double>& sample_weight,
	                              int N_azm_lat_lon,   // number of pulls in azimuth-lat-lon sets
	                              int N_zenith_speed) // number of pulls in zenith-speed sets
{
	//init the random generator
	// https://stackoverflow.com/questions/24334012/best-way-to-seed-mt19937-64-for-monte-carlo-simulations
	mt19937 rng;
	unsigned seed; 

	seed = random_device{}() * chrono::system_clock::now().time_since_epoch().count();
	rng.seed(seed);

	ofstream speed_angle_dist_file;
	speed_angle_dist_file.open("speed_vs_angle.txt");

	vector<double> zenith, vminv, vmaxv;
	vector<double> sample_zenith_0_i, sample_speed_0_i, sample_weight_i;
	vector<double> sample_zenith_f_i, sample_speed_f_i, sample_altitude_f_i, sample_distance_f_i;
	double dazm, d, latp_i, lonp_i, azm_i, azm_center, lats_i, lons_i, u_min, u_max;
	int i, j, i_zll, zll_count = 0;//, cur_i;

	// init output arrays
	sample_latp.clear();      
    sample_lonp.clear();     
    sample_azimuth_0.clear(); 
    sample_zenith_0.clear();
    sample_speed_0.clear();
    sample_azimuth_f.clear();  
    sample_zenith_f.clear();
    sample_speed_f.clear();
    sample_altitude_f.clear();
    sample_distance_f.clear();
    sample_weight.clear();

    // open lat-lon location file and get the size
    vector<double> lat_impact, lon_impact;
    get_latlon_arrays(0, p->latlon_idx_cur, N_azm_lat_lon, lat_impact, lon_impact);

  //   for (int i = 0; i < N_azm_lat_lon; ++i)
  //   {
  //   	latp_i = lat_impact[i];
		// lonp_i = lon_impact[i];
		// cout << latp_i << ' ' << lonp_i << endl;
  //   }

    cout << "Azm-lat-lon sample size = " << N_azm_lat_lon << endl;

	for (i_zll = 0; i_zll < N_azm_lat_lon; ++i_zll)
	{
		// do
		// { // make sure primary impact point is not underneath asset
			/////////////////////////////////////////////////////////
			/////// old code ////////////////////////////////////////
			// // randomly sample azimuth, latitude, and longitude
			// /// Note, the full range of latitude should be within +/- pi/2, and should be sampled uniformly on a sphere
			// u_min = (cos(latp + dlatp/2. + PI/2.) + 1.) / 2.;
			// u_max = (cos(latp - dlatp/2. + PI/2.) + 1.) / 2.;

			// //cout << u_min << ' ' << u_max << endl;

			// latp_i = acos( 2.*uniform(rng, u_min, u_max) - 1. ) - PI/2.;
			// //latp_i = uniform(rng, latp - dlatp/2., latp + dlatp/2.);
			// lonp_i = uniform(rng, lonp - dlonp/2., lonp + dlonp/2.);

			/////////////////////////////////////////////////////////
		// generated from a out to the antipode - a
		latp_i = lat_impact[i_zll];
		lonp_i = lon_impact[i_zll];
		//cout << latp_i << ' ' << lonp_i << endl;

			//cout << latp - dlatp/2. << ' ' << latp + dlatp/2. << ' ' << latp_i << ' ' << lonp_i << endl;

			// compute distance between primary and asset (units of rm), to center of asset
		d = lat_lon_dist(latp_i, lonp_i, lats, lons, 0);
		//} while(d <= r);

		d -= r; // shift d to be distance to edge, not center of asset

		// compute azimuth width, distances in units of rm
		dazm = azm_FOV(r, d);

		// compute bearing (i.e., azimuth) direction
		azm_center = azm_bearing(latp_i, lonp_i, lats, lons, 0);

		// compute azimuth sample
		azm_i = uniform(rng, azm_center - dazm/2., azm_center + dazm/2.);


		get_zenith_speed_grid(zenith, vminv, vmaxv, vmin, vmax, a, h, d, r, dg, dv);


		////cout << i_zll << " of " << N_azm_lat_lon << " | d, azm, dazm = " << d << ' ' << azm_center << ' ' << dazm << endl;

		if (zenith.size() == 0)
		{ // need to check to avoid errors with a zero-sized array
			cout << "No grid found for d = " << d << endl;
		
		}
		else
		{
			// CDF is normalized, starting from 0 to 1
			// will be used in conjunction with a uniform number generator that ranges from 0 to 1
			// Note: when pulling from the pdf, the points will have a uniform density
			vector<double> cdf, pdf;
			get_CDF_PDF_from_trapdens(zenith, vminv, vmaxv, cdf, pdf);

			// Next, we need to take samples from the CDF
			get_samples(zenith, // [rad]
				        vminv,  // [vesc]
				        vmaxv,  // [vesc]
				        vmin,   // [vesc]
				        vmax,   // [vesc]
				        cdf,
				        d, // [rm]
				        a, // [rm]
				        h, // [rm]
				        r, // [rm]
				        sample_zenith_0_i, // [rad]
				        sample_speed_0_i,  // [vesc]
				        sample_zenith_f_i, // [rad]
				        sample_speed_f_i,  // [vesc]
				        sample_altitude_f_i, // [rm]
				        sample_distance_f_i, // [rm]
				        sample_weight_i,
				        N_zenith_speed);

			

			// speed_angle_dist_file << d - r << ' ' << vSum(sample_weight) << endl;

			/// Next, the samples need to be checked against the actual asset to see if there is a hit or not

			cout << " d = " << d << " | " << 100.*(i_zll+1.)/double(N_azm_lat_lon) << "% finished | sum = " << vSum(sample_weight_i);
			cout << " | grid size = " << zenith.size() << "                       \r";
		

			// transfer i-th sample set to the main array
			zll_count++; // used to normalize after i_zll loop

			for (j = 0; j < N_zenith_speed; ++j)
			{
				// lat, lon, and azm will be the same for the i_zll-th speed-zenith sample
				sample_latp.emplace_back(latp_i);
				sample_lonp.emplace_back(lonp_i);
				sample_azimuth_0.emplace_back(azm_i);

				// copy from get_samples() output
		    	sample_zenith_0.emplace_back(sample_zenith_0_i[j]);
		    	sample_speed_0.emplace_back(sample_speed_0_i[j]);
		    	sample_altitude_f.emplace_back(sample_altitude_f_i[j]);
		    	sample_distance_f.emplace_back(sample_distance_f_i[j]);
		    	sample_weight.emplace_back(sample_weight_i[j]);

		    	// compute final azimuth (bearing), zenith angle, and speed of secondary as seen from asset
		    	/// For azm, need to compute the exact lat-lon impact given the azm_i sample
		    	lats_i = destination_lat(latp_i, sample_distance_f_i[j], azm_i);
		    	lons_i = destination_lon(latp_i, lonp_i, lats_i, sample_distance_f_i[j], azm_i);

		    	// final bearing (from asset to primary impact location), force to be from 0 to 2PI
				sample_azimuth_f.emplace_back( fmod(azm_bearing(lats_i, lons_i, latp_i, lonp_i, 0) + 2.*PI, 2.*PI) );

				// final zenith (as seen from asset, similar to wind direction), should be from 0 to PI (PI/2 to PI is pointing below the horizon)
				sample_zenith_f.emplace_back( final_zenith(sample_distance_f_i[j], sample_speed_0_i[j], sample_zenith_0_i[j], sample_altitude_f_i[j]) );

				// final speed at asset
				sample_speed_f.emplace_back( final_speed(sample_altitude_f_i[j], sample_speed_0_i[j]) );

				//cur_i = sample_latp.size()-1;

			}
		}
	}

	// need to divide the sample_weight by the number of azm_lat_lon pulls, not counting ones with no grid found
	for (i = 0; i < sample_latp.size(); ++i){
		sample_weight[i] /= double(zll_count); // ~Stratified Sampling
		sample_weight[i] *= p->lunar_escape_speed; // since vmin and vmax were in units of vesc

		speed_angle_dist_file << setprecision(10);

		speed_angle_dist_file << sample_latp[i] << ' ' << sample_lonp[i] << ' ' << sample_azimuth_0[i] << ' '; 
		speed_angle_dist_file << sample_zenith_0[i] << ' ' << sample_speed_0[i] << ' ' << sample_weight[i] << ' '; 
		speed_angle_dist_file << sample_azimuth_f[i] << ' ' << sample_zenith_f[i] << ' ' << sample_speed_f[i] << ' ';
		speed_angle_dist_file << sample_altitude_f[i] - 1. << ' ' << sample_distance_f[i]  << ' ' <<  lat_lon_dist(sample_latp[i], sample_lonp[i], lats, lons, 0) << endl; 
	}

	speed_angle_dist_file.close();

}