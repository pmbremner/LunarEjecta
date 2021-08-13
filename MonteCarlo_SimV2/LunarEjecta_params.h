#ifndef LUNAREJECTA_PARAMS_H
#define LUNAREJECTA_PARAMS_H


#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
//#include <float.h> // DBL_MAX


using namespace std;

const double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

enum primaryFluxType{HiDensMEM, LoDensMEM, NEO, secEjecta};
enum ejectaDistanceType{ejectaShort, ejectaFar};

// inline functions must be defined in the header
inline double sqr(double x) {return x*x;}
inline double min(double a, double b) {return (a < b ? a : b);}
inline double max(double a, double b) {return (a > b ? a : b);}
inline double mag_s(double x, double y, double z) {return sqrt(sqr(x) + sqr(y) + sqr(z));}
inline double mag2(double x, double y, double z) {return sqr(x) + sqr(y) + sqr(z);}


struct vec3
{
	double x[3];
};

struct mat3x3
{
	vec3 col[3];
};

double Beta(double a, double b);

double vMax(vector<double>& v);
double vMin(vector<double>& v);

template <class paramType>
void vCumLow(vector<paramType>& v, vector<paramType>& vCum)
{
	vCum[0] = 0;
	for (int i = 1; i < v.size(); ++i)
		vCum[i] = vCum[i-1] + v[i-1];
}

template <class paramType>
void vCumUp(vector<paramType>& v, vector<paramType>& vCum)
{
	vCum[0] = v[0];
	for (int i = 1; i < v.size(); ++i)
		vCum[i] = vCum[i-1] + v[i];
}

void shiftAngle(vector<double>& v, double min_ang, double shift_ang);


// https://www.learncpp.com/cpp-tutorial/function-pointers/
using RHS_func = double(*)(double, vector<double>&);
double findX(double LHS, RHS_func RHS, double a, double b, vector<double>& vars);

struct input
{
	int N_proc; // total number of processes
	int i_proc; // current process number
	int N_loc;  // number of locations = latlon_idx_max - latlon_idx_min + 1
	/////////////////////////////
	int latlon_idx_cur;  // current process lat-lon index (absolute)
	int latlon_idx_proc; // current process lat-lon index (relative, from idx_min)
	int latlon_idx_min;  // minimum lat-lon index (absolute)
	int latlon_idx_max;  // maximum lat-lon index (absolute, inclusive)

	int Nlat; // total number of latitudes
	int Nlon; // total number of longitudes
	int Nlatlon_tot; // Nlat*Nlon = total number of locations
	/////////////////////////////
	//bool initError; // 0 = no error, 1 = error
	bool readNEO_files; // 0 = generate the NEOs and save them, 1 = read from file
	bool saveNEO_files; // 0 = no, 1 = yes (where the options are for each lat)

	vector<string> lon_directory; // list of longitude directories

	string NEO_velDist_fn; // NEO velocity distribution filename

	// density params
	double MEM_hiDens_mu;    // log_e (not log_10)
	double MEM_hiDens_sigma; // ''
	double MEM_loDens_mu;    // ''
	double MEM_loDens_sigma; // ''
	double NEO_dens; // kg/m^3

	// regolith params
	double regolith_dens; // kg/m^3
	double HH11_porosity;
	double HH11_nu;
	double HH11_mu;
	double HH11_k;
	double HH11_C1;
	double HH11_C4;

	double HH11_p;
	double HH11_n1;
	double HH11_n2s; // for strength dominated
	double HH11_n2g; // for gravity dominated
	double HH11_H1;  // for strength dominated
	double HH11_H2;  // for gravity dominated


	double MEM_massMin; // minimum mass of MEM primaries, grams (cannot be lower than 1E-6 g)
	double MEM_massMax; // maximum mass of MEM primaries, grams (cannot be higher than 10 g)
	double NEO_massMin; // minimum mass of NEOs, grams
	double NEO_massMax; // maximum mass of NEOs, grams (cannot go to infinity, need a finite cutoff)
	/////////////////////////////
	double lunar_radius; // km
	double lunar_escape_speed; // km/s
	/////////////////////////////
	double ROI_radius;   // km
	double ROI_lat; // rad
	double ROI_lon; // rad
	/////////////////////////////
	int N_D_perRegion;
	int N_bearing_POI;
	int ROI_sample_points;
	//int arc_sample_points;
	/////////////////////////////
	int N_horizon_ROI; // secondary ejecta
	int N_bearing_ROI; // secondary ejecta
	int N_vel;         // secondary ejecta
	double vel_min; // km/s


};

void linspace(vector<double>& x, double xmin, double xmax, int Nx);
void logspace(vector<double>& x, double xmin, double xmax, int Nx, int i0, int im);

double rand_uniform(double min, double max);

double calcSA(double cur_lat, double lat_min, double lat_max, double dlat, double dlon);
void compute_igloo_azm_bin_number(vector<int> &phi, int Nphi);
void get_lat_min_max(double cur_lat, double dlat, double& lat_min, double& lat_max);

double NEO_integral_flux(double m); // m in grams
double H_calcH11_C4(input* p);

// https://www.cplusplus.com/doc/oldtutorial/templates/
template <class paramType>
void getParam(string param_fn, string paramLabel, paramType& param, bool defaultExists) {
	ifstream input_file;

	input_file.open("./param_files/" + param_fn);
	char C_Line[64];
	paramType tparam;
	bool isError = 1;

	int lineNumber = 1;

	while (input_file.getline(C_Line, 64, '#'))
	{
		string paramLabel_file(C_Line);
		
		if(paramLabel_file == paramLabel)
		{
			if(input_file.getline(C_Line, 64, '#'));
				isError = 0;
			string tstr(C_Line);
			stringstream SS_param(tstr);
			SS_param >> param;
			//input_file.ignore(256, '\n'); // temp
			if(defaultExists)
				cout << "  Overriding Default value...\n";
			cout << "Line " << lineNumber << " | " << paramLabel << " = " << param << endl;
			//return;
		}
		else
		{
			input_file.ignore(256, '\n');
		}

		lineNumber++;
	}
	// if there's an error, set the initError flag
	if(!defaultExists && isError){
		//initError = 1;
		cout << "ERROR: Parameter '" << paramLabel << "' not found...\n";
	}

	input_file.close();
}



input* init_input(string param_fn, int N_proc, int i_proc);



#endif 