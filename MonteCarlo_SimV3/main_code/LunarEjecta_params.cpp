#include "LunarEjecta_params.h"

#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>     /* srand, rand */

using namespace std;


double Beta(double a, double b) {
	// helps to avoid overflow errors doing it this way
	return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
}

void shiftAngle(vector<double>& v, double min_ang, double shift_ang)
{
	for (int i = 0; i < v.size(); ++i)
		if (v[i] > min_ang)
			v[i] -= shift_ang;
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

// mass m is in units of grams
// returns flux [#/m^2/yr]
// See Grun et al 1985, eq A3
double MEM_mass_grun(double m)
{
	const double sec_per_yr = 60.*60.*24.*365.;

	const double c4  = 2.2E3;
	const double c5  = 15.;
	const double c6  = 1.3E-9;
	const double c7  = 1.E11;
	const double c8  = 1.E27;
	const double c9  = 1.3E-16;
	const double c10 = 1.E6;

	const double g4  = 0.306; 
	const double g5  = -4.38;
	const double g6  = 2.;
	const double g7  = 4.;
	const double g8  = -0.36;  
	const double g9  = 2.;
	const double g10 = -0.85;

	return (pow(c4 * pow(m, g4) + c5, g5)
		   + c6 * pow(m + c7 * pow(m, g6) + c8 * pow(m, g7), g8)
		   + c9 * pow(m + c10 * pow(m, g9), g10) ) * sec_per_yr;
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

	cout << "--------------------------------\n";
	cout << "Reading... " << param_fn << endl;

	//////////////////////
	// Default parameters
	p->MEM_massMin = 1.E-6; // g
	p->MEM_massMax = 10.; // g
	p->NEO_massMin = 10.; // g
	p->NEO_massMax = 5524.; // g, mass of a 0.1% chance of hitting the moon in 11 years

	p->lunar_radius = 1737.4E3; // m (suggested by NESC)
	p->lunar_escape_speed = 2.38E3; // m/s
	p->lunar_acceleration = 1.625;    // m/s^2 at Moon's surface

	p->MEM_hiDens_mu    = 8.241; // log_e (not log_10)
	p->MEM_hiDens_sigma = 0.214;
	p->MEM_loDens_mu    = 6.753;
	p->MEM_loDens_sigma = 0.292;

	p->NEO_dens = 3000.; // kg/m^3

	// fit from Carrier 2003
	p->reg_size_mu    = -2.649;
	p->reg_size_sigma = 1.786;
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

	getParam(param_fn, "asset_radius", p->asset_radius, 0); // m to rm
	getParam(param_fn, "asset_altitude", p->asset_altitude, 0); // m to rm
	getParam(param_fn, "asset_height", p->asset_height, 0); // m to rm

	p->asset_radius /= p->lunar_radius;
	p->asset_altitude /= p->lunar_radius;
	p->asset_height /= p->lunar_radius;

	p->asset_altitude += 1.; // input was height above surface, need distance from center of Moon

	getParam(param_fn, "asset_lat", p->asset_lat, 0); // read in degrees
	p->asset_lat *= PI/180.;
	getParam(param_fn, "asset_lon", p->asset_lon, 0); // read in degrees
	p->asset_lon *= PI/180.;

	getParam(param_fn, "regolith_dens", p->regolith_dens, 0);
	getParam(param_fn, "regolith_porosity", p->regolith_porosity, 0);
	getParam(param_fn, "regolith_tensile_strength", p->regolith_tensile_strength, 0);
	getParam(param_fn, "regolith_type", temp_s, 0);
	getParam(param_fn, "reg_min_size", p->reg_min_size, 0);
	getParam(param_fn, "reg_max_size", p->reg_max_size, 0);
	



	getParam(param_fn, "dg_max", p->dg_max, 0);
	getParam(param_fn, "dv_max", p->dv_max, 0);


	getParam(param_fn, "N_azm_lat_lon", p->N_azm_lat_lon, 0);
	getParam(param_fn, "N_zenith_speed", p->N_zenith_speed, 0);
	getParam(param_fn, "N_primary_sample", p->N_primary_sample, 0);


	getParam(param_fn, "N_env_v", p->N_env_v, 0);
	getParam(param_fn, "N_env_zen", p->N_env_zen, 0);
	getParam(param_fn, "N_env_azm", p->N_env_azm, 0);
	getParam(param_fn, "N_env_size", p->N_env_size, 0);

	p->N_env_flux = p->N_env_v * p->N_env_zen * p->N_env_azm * p->N_env_size;




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
		p->HH11_n2g      = 1.3; // estimate, using from sand
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


	getParam(param_fn, "vel_min", p->vel_min, 0);
	getParam(param_fn, "vel_max", p->vel_max, 0);

	p->vel_min /= p->lunar_escape_speed;
	p->vel_max /= p->lunar_escape_speed;



	cout << "--------------------------------\n";
	return p;
}

string get_latlon_fn(int i_proc, int idx)
{
	string s_proc = to_string(i_proc);
	string s_idx = to_string(idx);
	string fn = "latlon_loc_" + string(3 - min(3, s_proc.length()), '0') + s_proc + '_' + string(5 - min(5, s_idx.length()), '0') + s_idx + ".txt";

	return fn;
}

// https://stackoverflow.com/questions/5590381/easiest-way-to-convert-int-to-string-in-c
// https://stackoverflow.com/questions/6143824/add-leading-zeros-to-string-without-sprintf
int get_latlon_size(int i_proc, int idx)
{
	int N = 0;
	string line;
	string fn = get_latlon_fn(i_proc, idx);
	ifstream latlon_file(fn);

	while(getline(latlon_file, line))
		N++;

	latlon_file.close();

	cout << "Number of lat-lon samples: " << N << " in file: " << fn << endl;

	return N;
}

void get_latlon_arrays(int i_proc, int idx, int N, vector<double>& lat, vector<double>& lon)
{
	string fn = get_latlon_fn(i_proc, idx);
	ifstream latlon_file(fn);
	string line, col;

	lat.clear();
	lon.clear();

	lat.resize(N);
	lon.resize(N);

	int i = 0;

	while(getline(latlon_file, line))
	{
		//cout << line << endl;

		istringstream iss_line(line);

		// get lat
		getline(iss_line, col, ' ');
		lat[i] = stof(col);


		// get lon
		getline(iss_line, col, '\n');
		lon[i] = stof(col);

		i++;
		//cout << lat[i] << ' ' << lon[i] << endl;
	}
}