#include "LunarEjecta_MainUtils.h"
#include "LunarEjecta_SecondaryEjecta.h"


// note, -march=native is to allow for vectorization, if possible

//  g++ -O2 -std=c++17 -march=native secondaries_test.cpp ./main_code/LunarEjecta_MainUtils.cpp ./main_code/LunarEjecta_SecondaryEjecta.cpp -IC:\Users\AMD-Y500\Documents\GitHub\LunarEjecta\MonteCarlo_SimV3\main_code -o ejecta.exe

using namespace std;



int main(int argc, char const *argv[])
{
	//int N_proc = atoi(argv[2])
	//vector<double> zenith, vminv, vmaxv;

	// double vlow = 0.0;
	// double vmax = 5.;
	// double h = 0.22;
	// double d = 0.56;//0.287;
	// double r = 0.21;

	// double dg = 0.1;
	// double dv = 0.05;
	const double Rm = 1737.E3; // m

	double vmin = 0.0;
	double vmax = 3.;
	double a = atof(argv[1])/Rm + 1.; // altitude above lunar surface 
	double h = atof(argv[2])/Rm;//120./Rm;
	//double d; //= //1000./Rm;//0.287;
	double r = atof(argv[3])/Rm;//4.5/Rm;

	double dg = 0.05;
	double dv = 0.05;	

	// int Nd = 50;

	// vector<double> d_vec;
	// logspace(d_vec, log10(r + 10./Rm), log10(2.*PI - (2.*r + 10./Rm)), Nd);

	// for (int i = 0; i < d_vec.size(); ++i)
	// {
	// 	cout << i << ' ' << d_vec[i] << endl;
	// }
	// ofstream ejecta_vs_dist_file;
	// string file_pre("ejecta_vs_dist_");
 //    string file_post("_.txt");
 //    string file_label(to_string(a*Rm) + "_" + to_string(h*Rm) + "_" + to_string(r*Rm));
	// ejecta_vs_dist_file.open(file_pre + file_label + file_post);

	// cout << file_pre + file_label + file_post << endl;


	// ofstream speed_angle_dist_file;
	// string file_pre2("_speed_vs_angle");


	vector<double> sample_latp, sample_lonp, sample_azimuth_0, sample_zenith_0, sample_speed_0;
	vector<double> sample_azimuth_f, sample_zenith_f, sample_speed_f, sample_weight;
	int N_azm_lat_lon = 100, N_zenith_speed = 100;


	get_samples_with_azm_lat_lon( 17 / 180. * PI,   // primary latitude center
	                              110. / 180. * PI,   // primary longitude center
	                              0.005 / 180. * PI,  // primary latitude range
	                              0.005 / 180. * PI,  // primary longitude range
	                              15.01 / 180. * PI,   // satellite (asset) latitude center
	                              110.01 / 180. * PI,   // satellite (asset) longitude center
	                              a,      // satellite (asset) altitude [rm]
	                              h,      // satellite (asset) height [rm]
	                              r,      // satellite (asset) radius [rm]
	                              vmin,   // minimum ejecta speed [vesc]
	                              vmax,   // maximum ejecta speed [vesc]
	                              dg,     // maximum zenith grid width
	                              dv,     // maximum speed grid width
	                              sample_latp,       // primary latitude center
	                              sample_lonp,       // primary longitude center
	                              sample_azimuth_0,  // initial azimuth, zenith, and speed of secondary, at primary impact
	                              sample_zenith_0,
	                              sample_speed_0,
	                              sample_azimuth_f,  // final azimuth, zenith, and speed of secondary, at asset impact
	                              sample_zenith_f,
	                              sample_speed_f,
	                              sample_weight,
	                              N_azm_lat_lon,   // number of pulls in azimuth-lat-lon sets
	                              N_zenith_speed); // number of pulls in zenith-speed sets



	return 0;
}