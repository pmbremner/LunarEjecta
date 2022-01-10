#include "LunarEjecta_MainUtils.h"
#include "LunarEjecta_SecondaryEjecta.h"


// note, -march=native is to allow for vectorization, if possible

//  g++ -O2 -std=c++17 -march=native domain_finder_test.cpp ./main_code/LunarEjecta_MainUtils.cpp ./main_code/LunarEjecta_SecondaryEjecta.cpp -IC:\Users\AMD-Y500\Documents\GitHub\LunarEjecta\MonteCarlo_SimV3\main_code -o ejecta.exe

using namespace std;



int main(int argc, char const *argv[])
{
	//int N_proc = atoi(argv[2])
	vector<double> zenith, vminv, vmaxv;

	// double vlow = 0.0;
	// double vmax = 5.;
	// double h = 0.22;
	// double d = 0.56;//0.287;
	// double r = 0.21;

	// double dg = 0.1;
	// double dv = 0.05;
	const double Rm = 1737.E3; // m

	double vlow = 0.0;
	double vmax = 5.;
	double a = atof(argv[1])/Rm; // altitude above lunar surface 
	double h = atof(argv[2])/Rm;//120./Rm;
	double d; //= //1000./Rm;//0.287;
	double r = atof(argv[3])/Rm;//4.5/Rm;

	double dg = 0.05;
	double dv = 0.05;	

	int Nd = 50;

	vector<double> d_vec;
	logspace(d_vec, log10(r + 10./Rm), log10(2.*PI - (2.*r + 10./Rm)), Nd);

	for (int i = 0; i < d_vec.size(); ++i)
	{
		cout << i << ' ' << d_vec[i] << endl;
	}
	ofstream ejecta_vs_dist_file;
	string file_pre("ejecta_vs_dist_");
    string file_post("_.txt");
    string file_label(to_string(a*Rm) + "_" + to_string(h*Rm) + "_" + to_string(r*Rm));
	ejecta_vs_dist_file.open(file_pre + file_label + file_post);

	cout << file_pre + file_label + file_post << endl;


	ofstream speed_angle_dist_file;
	string file_pre2("_speed_vs_angle");


	vector<double> sample_zenith, sample_speed, sample_weight;
	int N_sample = 5000;

	for (int i = 0; i < Nd; ++i)
	{
		d = d_vec[i];
		//string file_label2(to_string(d*Rm));
		string scount;
		scount = to_string(i);
		auto new_str = string(3 - min(3, scount.length()), '0') + scount;

		speed_angle_dist_file.open(new_str + file_pre2 + file_post);


		// The zenith grid points are at most dg and are made smaller if the speed is changing a lot
		// over a small zenith change. The speed change is limited to at most dv
		// All for corners of a cross-section (r-theta plane) of a wedge are used, but are ignored if the speed is above vmax
		// For each zenith grid point, the minimum and maximum speed extent is computed as well (vminv, and vmaxv)
		// The wedge has a height h, starting at a
		// d is the distance to the center of the wedge, with a thickness 2*r, r in each direction
		get_zenith_speed_grid(zenith, vminv, vmaxv, vlow, vmax, a, h, d, r, dg, dv);
		//get_zenith_speed_grid(zenith, vminv, vmaxv, vlow, v, h, 2.*PI-d, r, dg, dv);

		if (zenith.size() == 0)
		{ // need to check to avoid errors with a zero-sized array
			cout << "No grid found for d = " << d << endl;
		
		}
		else
		{
			// for (int i = 0; i < vminv.size(); ++i)
			// {
			// 	vminv[i] = 0.;
			// 	vmaxv[i] = 10.;
			// }

			// CDF is normalized, starting from 0 to 1
			// will be used in conjunction with a uniform number generator that ranges from 0 to 1
			vector<double> cdf, pdf;
			get_CDF_PDF_from_trapdens(zenith, vminv, vmaxv, cdf, pdf);

			// Next, we need to take samples from the CDF
			get_samples(zenith, vminv, vmaxv, vlow, vmax, cdf, sample_zenith, sample_speed, sample_weight, N_sample);

			for (int j = 0; j < sample_zenith.size(); ++j)
				speed_angle_dist_file << sample_zenith[j] << ' ' << sample_speed[j] << ' ' << sample_weight[j] << endl; 

			ejecta_vs_dist_file << d - r << ' ' << vSum(sample_weight) << endl;

			/// Next, the samples need to be checked against the actual asset to see if there is a hit or not

			cout << "d = " << d -r << " | " << 100.*(i+1.)/double(Nd) << "% finished | sum = " << vSum(sample_weight);
			cout << " | grid size = " << zenith.size() << "                     \r";
		}
		speed_angle_dist_file.close();
		
	}

	ejecta_vs_dist_file.close();

	

	return 0;
}