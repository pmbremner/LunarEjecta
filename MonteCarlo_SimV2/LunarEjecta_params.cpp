#include "LunarEjecta_params.h"


#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>     /* srand, rand */

using namespace std;


double Beta(double a, double b) {
	// helps to avoid overflow errors doing it this way
	return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
}

double vMax(vector<double>& v){
	double max = v[0];
	for (int i = 1; i < v.size(); ++i)
		if (v[i] > max)
			max = v[i];

	return max;
}

double vMin(vector<double>& v){
	double min = v[0];
	for (int i = 1; i < v.size(); ++i)
		if (v[i] < min)
			min = v[i];

	return min;
}


void shiftAngle(vector<double>& v, double min_ang, double shift_ang)
{
	for (int i = 0; i < v.size(); ++i)
		if (v[i] > min_ang)
			v[i] -= shift_ang;
}

void linspace(vector<double>& x, double xmin, double xmax, int Nx) {
	x.resize(Nx);
	for (int i = 0; i < Nx; ++i)
		x[i] = xmin + (xmax - xmin) * double(i) / double(Nx-1.);
}

void logspace(vector<double>& x, double xmin, double xmax, int Nx, int i0, int im) {
	x.resize(im-i0);
	for (int i = i0; i < im; ++i) {
		x[i-i0] = xmin * pow(xmax/xmin, double(i)/double(Nx-1.));
		//cout << x[i-i0] << endl;
	}
}

double rand_uniform(double min, double max)
{
	return min + (max - min) * rand() / double(RAND_MAX);
}


double calcSA(double cur_lat, double lat_min, double lat_max, double dlat, double dlon)
{
	if(fabs(cur_lat) < dlat/2.)
		return dlon * (fabs(sin(lat_max)) + fabs(sin(lat_min)));
	else
		return dlon * fabs(sin(lat_max) - sin(lat_min));
}

// Nphi is the expected number of azimuth bins at the horizon
void compute_igloo_azm_bin_number(vector<int> &phi, int Nphi)
{
	double dzen = PI/double(2.*phi.size());
	for (int i = 0; i < phi.size(); ++i){
		phi[i] = round(fabs(Nphi * ( sin(i*dzen) + sin((i+1.)*dzen) )/2.));

		//cout << phi[i] << endl;
	}
}

// all in radians
void get_lat_min_max(double cur_lat, double dlat, double& lat_min, double& lat_max)
{
	if(fabs(cur_lat - PI/2.) < 1E-4) // +90 lat
	{
		lat_max = PI/2.;
		lat_min = lat_max - dlat/2.;
	}
	else if(fabs(cur_lat + PI/2.) < 1E-4) // -90 lat
	{
		lat_min = -PI/2.;
		lat_max = lat_min + dlat/2.;
	}
	else
	{
		lat_min = cur_lat - dlat/2.;
		lat_max = cur_lat + dlat/2.;
	}
}


// input mass in units of grams, equation uses mass in kg
double NEO_integral_flux(double m)
{
	return 2.89E-11 * pow(m/1000., -0.9); // units of #/m^2/yr
}

double H_calcH11_C4(input* p)
{
	return 3.*p->HH11_k/(4.*PI) * pow(p->HH11_C1, 3.*p->HH11_mu);
}

input* init_input(string param_fn, int N_proc, int i_proc)
{
	string temp_s;
	input* p = new input;

	p->N_proc = N_proc;
	p->i_proc = i_proc;

	//////////////////////
	// Default parameters
	p->MEM_massMin = 1.E-6; // g
	p->MEM_massMax = 10.; // g
	p->NEO_massMin = 10.; // g
	p->NEO_massMax = 1.5708E15; // g, 1000 m diameter at 3 g/cc

	p->lunar_radius = 1737.4E3; // m (suggested by NESC)
	p->lunar_escape_speed = 2.38E3; // m/s
	p->lunar_acceleration = 1.625;    // m/s^2 at Moon's surface

	p->MEM_hiDens_mu    = 8.241; // log_e (not log_10)
	p->MEM_hiDens_sigma = 0.214;
	p->MEM_loDens_mu    = 6.753;
	p->MEM_loDens_sigma = 0.292;

	p->NEO_dens = 3000.; // kg/m^3
	//////////////////////

	//////////////////////
	// Required parameters

	getParam(param_fn, "Nlat", p->Nlat, 0);
	getParam(param_fn, "Nlon", p->Nlon, 0);

	p->Nlatlon_tot = p->Nlat * p->Nlon;
	p->N_loc       = ceil(p->Nlatlon_tot / double(p->N_proc));

	p->latlon_idx_min = p->i_proc * p->N_loc;
	p->latlon_idx_max = fmin((p->i_proc + 1) * p->N_loc - 1, p->Nlatlon_tot - 1);

	p->N_loc = p->latlon_idx_max - p->latlon_idx_min + 1; // fix number of locations (mostly for the last process)

	cout << "min and max lat-lon index (absolute) = " << p->latlon_idx_min << ", " << p->latlon_idx_max << " | for proc = " << p->N_loc << " | total = " << p->Nlatlon_tot << endl;

	p->lon_directory.resize(4);
	for (int i = 0; i < p->Nlon; ++i)
		getParam(param_fn, "MEMDirectoryName_" + to_string(i * 360/p->Nlon), p->lon_directory[i], 0);

	getParam(param_fn, "NEO_velDist_fn", p->NEO_velDist_fn, 0);
	getParam(param_fn, "readNEO_files", p->readNEO_files, 0);
	getParam(param_fn, "saveNEO_files", p->saveNEO_files, 0);

	getParam(param_fn, "ROI_radius", p->ROI_radius, 0); // km
	getParam(param_fn, "ROI_lat", p->ROI_lat, 0); // read in degrees
	p->ROI_lat *= PI/180.;
	getParam(param_fn, "ROI_lon", p->ROI_lon, 0); // read in degrees
	p->ROI_lon *= PI/180.;

	getParam(param_fn, "regolith_dens", p->regolith_dens, 0);
	getParam(param_fn, "regolith_type", temp_s, 0);


	p->HH11_nu = 0.4; // see footnote 5 of Housen Holsapple 2011
	p->HH11_n1 = 1.2;
	p->HH11_n2s = 1.;

	

	if (temp_s == "rock")
	{
		p->HH11_porosity = 0.0;
		p->HH11_mu       = 0.55;
		p->HH11_C1       = 1.5;
		p->HH11_k        = 0.3;
		p->HH11_p        = 0.5;
		p->HH11_n2g      = 1.5;
		p->HH11_H1       = 0.68; // estimate, using from water
		p->HH11_H2       = 1.1;
	}
	else if (temp_s == "weaklyCementedBasalt")
	{
		p->HH11_porosity = 0.20;
		p->HH11_mu       = 0.46;
		p->HH11_C1       = 0.18;
		p->HH11_k        = 0.3;
		p->HH11_p        = 0.3;
		p->HH11_n2g      = 1.3;
		p->HH11_H1       = 0.5; // estimate, no value given
		p->HH11_H2       = 0.38;
	}
	else if (temp_s == "sand")
	{
		p->HH11_porosity = 0.35;
		p->HH11_mu       = 0.41;
		p->HH11_C1       = 0.55;
		p->HH11_k        = 0.3;
		p->HH11_p        = 0.3;
		p->HH11_n2g      = 1.3;
		p->HH11_H1       = 0.59;
		p->HH11_H2       = 0.4; // estimate, using from SFA
	}
	else if (temp_s == "glassMicroSpheres")
	{
		p->HH11_porosity = 0.36;
		p->HH11_mu       = 0.45;
		p->HH11_C1       = 1.0;
		p->HH11_k        = 0.5;
		p->HH11_p        = 0.3;
		p->HH11_n2g      = 1.3;
		p->HH11_H1       = 0.8;
		p->HH11_H2       = 0.96; // estimate, between rock and PS
	}
	else if (temp_s == "sandFlyAsh")
	{
		p->HH11_porosity = 0.45;
		p->HH11_mu       = 0.4;
		p->HH11_C1       = 0.55;
		p->HH11_k        = 0.3;
		p->HH11_p        = 0.3;
		p->HH11_n2g      = 1.2; // estimate, no value given
		p->HH11_H1       = 0.59; // estimate, using from sand
		p->HH11_H2       = 0.4;
	}
	else if (temp_s == "perliteSandMixture")
	{
		p->HH11_porosity = 0.6;
		p->HH11_mu       = 0.35;
		p->HH11_C1       = 0.6;
		p->HH11_k        = 0.32;
		p->HH11_p        = 0.2;
		p->HH11_n2g      = 1.2; // estimate, no value given
		p->HH11_H1       = 0.59; // estimate, using from sand
		p->HH11_H2       = 0.81;
	}
	else
	{
		cout << "ERROR: invalid reolith type in init_input\n\n";
	}

	p->HH11_C4       = H_calcH11_C4(p);
	cout << "HH11 C4 = " << p->HH11_C4 << endl;

	//////////////////////
	getParam(param_fn, "N_hit", p->N_hit, 0);
	getParam(param_fn, "N_max", p->N_max, 0);
	getParam(param_fn, "alpha_search", p->alpha_search, 0);
	getParam(param_fn, "lifetime_max", p->lifetime_max, 0);
	getParam(param_fn, "lifetime_rate", p->lifetime_rate, 0);
	getParam(param_fn, "dx_rate", p->dx_rate, 0);


	getParam(param_fn, "N_D_perRegion", p->N_D_perRegion, 0);
	getParam(param_fn, "N_bearing_POI", p->N_bearing_POI, 0);
	getParam(param_fn, "N_horizon_ROI", p->N_horizon_ROI, 0);
	getParam(param_fn, "N_bearing_ROI", p->N_bearing_ROI, 0);
	getParam(param_fn, "ROI_sample_points", p->ROI_sample_points, 0);
	getParam(param_fn, "N_vel", p->N_vel, 0);
	getParam(param_fn, "vel_min", p->vel_min, 0);
	getParam(param_fn, "vel_max", p->vel_max, 0);
	//getParam(param_fn, "arc_sample_points", p->arc_sample_points, 0);


	return p;
}

// Modified False Position, Chapter 4.5 of Numerical Methods for Scientists and Engineers
double findX(double LHS, RHS_func RHS, double a, double b, vector<double>& vars)
{
	double fa = RHS(a, vars) - LHS;
	double fb = RHS(b, vars) - LHS;
	double fx, x;
	int count = 0;
	
	do{
		x = (a*fb - b*fa) / (fb - fa);
		fx = RHS(x, vars) - LHS;

		//cout << a << ' ' << b << ' ' << x << ' ' << fa << ' ' << fb << ' ' << fx << endl;

		if (fa*fx < 0)
		{
			b = x;
			fb = fx;
			fa /= 2.;
		}
		else if (fa*fx > 0)
		{
			a = x;
			fa = fx;
			fb /= 2.;
		}
		else
		{
			return x;
		}
		count++;

	} while (fabs(2*(a-b)/(a+b)) > 1.E-8);
	
	//cout << endl;
	//cout << count << ' ';

	return x;
}