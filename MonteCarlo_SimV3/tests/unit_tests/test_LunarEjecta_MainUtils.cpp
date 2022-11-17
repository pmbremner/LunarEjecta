/**
 * A simple unit test code meant to test the
 * functionality of utilities contained in
 * LunarEjecta_MainUtils
 * 
 * Vector operations were check via numpy equivalents
 */

#include <string>
#include <limits>
#include <vector>
#include "../../main_code/LunarEjecta_MainUtils.cpp"


using namespace std;

typedef std::numeric_limits< double > dbl;


struct test_container
{
	// Files
	string floatingpoint_in;
	string integer_in;
	string floatingpoint_er;
	string integer_er;
	
	// Create test vectors
	vector<double> fp_vec;
	vector<int> int_vec;

	// Parse the expected results
	double expected_vMax;
	double expected_vMin;
	double expected_vSum;
};



int main()
{
	std::cout.precision(dbl::max_digits10);

	// Setup a tolerance for aritmitic operations
	double tolerance = 1e-10;


	// Informational Message
	string intro_msg = "\nTesting LunarEjecta_MainUtils...";
	cout << intro_msg << "\n\ntest tolerance = " << tolerance << endl;
	


	// Include the test data directory and files
	// ---------------------------------------------------------
	
	// Define the data directory
	string test_data_dir = "../test_data";
	string test_er_dir = "../test_expected_results";

	// Name the data files
	string integer_in_1 = (test_data_dir + "/random_integer_values_1.txt");
	string integer_in_2 = (test_data_dir + "/random_integer_values_2.txt");
	string integer_in_3 = (test_data_dir + "/random_integer_values_3.txt");
	string integer_in_4 = (test_data_dir + "/random_integer_values_4.txt");
	string floatingpoint_in_1 = (test_data_dir + "/random_floating_point_values_1.txt");
	string floatingpoint_in_2 = (test_data_dir + "/random_floating_point_values_2.txt");
	string floatingpoint_in_3 = (test_data_dir + "/random_floating_point_values_3.txt");
	string floatingpoint_in_4 = (test_data_dir + "/random_floating_point_values_4.txt");

	// Name the expected results files
	string integer_er_1 = (test_er_dir + "/random_integer_expectedresults_1.txt");
	string integer_er_2 = (test_er_dir + "/random_integer_expectedresults_2.txt");
	string integer_er_3 = (test_er_dir + "/random_integer_expectedresults_3.txt");
	string integer_er_4 = (test_er_dir + "/random_integer_expectedresults_4.txt");
	string floatingpoint_er_1 = (test_er_dir + "/random_floating_point_expectedresults_1.txt");
	string floatingpoint_er_2 = (test_er_dir + "/random_floating_point_expectedresults_2.txt");
	string floatingpoint_er_3 = (test_er_dir + "/random_floating_point_expectedresults_3.txt");
	string floatingpoint_er_4 = (test_er_dir + "/random_floating_point_expectedresults_4.txt");

	vector<string> test_files, expected_results_files;

	test_files.push_back (floatingpoint_in_1);
	test_files.push_back (floatingpoint_in_2);
	test_files.push_back (floatingpoint_in_3);
	test_files.push_back (floatingpoint_in_4);

	expected_results_files.push_back (floatingpoint_er_1);
	expected_results_files.push_back (floatingpoint_er_2);
	expected_results_files.push_back (floatingpoint_er_3);
	expected_results_files.push_back (floatingpoint_er_4);
	// ---------------------------------------------------------
	

	// Create test result vectors
	//vector<test_container> test_results;
	

	// Read in the test data and the expected values
	// ---------------------------------------------------------
	for (unsigned int i=0; i<test_files.size(); ++i)
	{

		test_container test_bundle;

		ifstream testf_in, testf_er;
		testf_in.open(test_files[i]);
		testf_er.open(expected_results_files[i]);

		
		string temp;
		while (getline(testf_in, temp))
		{
			istringstream line(temp);
			double val;

			line >> val;
			test_bundle.fp_vec.push_back (val);
		}
		while (getline(testf_er, temp))
		{
			istringstream line(temp);
			string descriptor;
			double val;
			
			line >> descriptor >> val;

			if (descriptor.compare("vMax") == 0)
			{
				test_bundle.expected_vMax = val;
			}
			else if (descriptor.compare("vMin") == 0)
			{
				test_bundle.expected_vMin = val;
			}
			else if (descriptor.compare("vSum") == 0)
			{
				test_bundle.expected_vSum = val;
			}
		}
		testf_in.close();
		testf_er.close();
		// ---------------------------------------------------------

		cout << "\n\tTEST " << i+1 << "\n" << endl;

		// Test vMax
		double diff = vMax(test_bundle.fp_vec) - test_bundle.expected_vMax;
		if (abs(diff)<=tolerance)
		{
			cout << "vMax\t\tPASSED" << endl;
		}
		else
		{
			cout << "vMax\t\tFAILED\tDifference = " << diff << endl;
		}


		// Test vMin
		diff = vMin(test_bundle.fp_vec) - test_bundle.expected_vMin;
		if (abs(diff)<=tolerance)
		{
			cout << "vMin\t\tPASSED" << endl;
		}
		else
		{
			cout << "vMin\t\tFAILED\tDifference = " << diff << endl;
		}


		// Test vSum
		diff = vSum(test_bundle.fp_vec) - test_bundle.expected_vSum;
		if (abs(diff)<=tolerance)
		{
			cout << "vSum\t\tPASSED" << endl;
		}
		else
		{
			cout << "vSum\t\tFAILED\tDifference = " << diff << endl;
		}
	}
	
	

	// Test linspace
	vector<double> v;
	linspace(v, 0.2, 20.0, 50);

	vector<double> v_expected = {0.2,  0.60408163,  1.00816327,  1.4122449 ,  1.81632653,
        2.22040816,  2.6244898 ,  3.02857143,  3.43265306,  3.83673469,
        4.24081633,  4.64489796,  5.04897959,  5.45306122,  5.85714286,
        6.26122449,  6.66530612,  7.06938776,  7.47346939,  7.87755102,
        8.28163265,  8.68571429,  9.08979592,  9.49387755,  9.89795918,
       10.30204082, 10.70612245, 11.11020408, 11.51428571, 11.91836735,
       12.32244898, 12.72653061, 13.13061224, 13.53469388, 13.93877551,
       14.34285714, 14.74693878, 15.15102041, 15.55510204, 15.95918367,
       16.36326531, 16.76734694, 17.17142857, 17.5755102 , 17.97959184,
       18.38367347, 18.7877551 , 19.19183673, 19.59591837, 20.};

	//for (const auto i_v : v) cout << i_v << " "; cout << endl;
	double max_diff = 0;
	unsigned int fail_count = 0;
	for (unsigned int i=0; i < v.size(); ++i)
	{
		double diff = abs(v[i] - v_expected[i]);
		if (diff>max_diff) max_diff = diff;
		if (diff>tolerance) ++fail_count;
	}
	if (max_diff<=tolerance)
	{
		cout << "\n\nlinspace\t\tPASSED" << endl;
	}
	else
	{
		cout << "\n\nlinspace\t\tFAILED\t"<< fail_count << "/" << v.size() << " failed elements\tMax Abs Difference = " << max_diff << endl;
	}


	// Test logspace
	vector<double> v1;
	logspace(v1, 0.1, 5.1, 50);
	v_expected = {1.25892541e+00, 1.59235837e+00, 2.01410280e+00, 2.54754843e+00,
       3.22227992e+00, 4.07571757e+00, 5.15519263e+00, 6.52057229e+00,
       8.24757988e+00, 1.04319944e+01, 1.31949626e+01, 1.66897173e+01,
       2.11100760e+01, 2.67011897e+01, 3.37731391e+01, 4.27181312e+01,
       5.40322511e+01, 6.83429746e+01, 8.64439680e+01, 1.09339104e+02,
       1.38298136e+02, 1.74927119e+02, 2.21257479e+02, 2.79858676e+02,
       3.53980707e+02, 4.47734343e+02, 5.66319120e+02, 7.16311693e+02,
       9.06030582e+02, 1.14599751e+03, 1.44952093e+03, 1.83343411e+03,
       2.31902869e+03, 2.93323554e+03, 3.71011828e+03, 4.69276246e+03,
       5.93566508e+03, 7.50775694e+03, 9.49622553e+03, 1.20113504e+04,
       1.51926192e+04, 1.92164637e+04, 2.43060443e+04, 3.07436270e+04,
       3.88862370e+04, 4.91854597e+04, 6.22124853e+04, 7.86897866e+04,
       9.95311870e+04, 1.25892541e+05};
	
	max_diff = 0;
	fail_count = 0;
	for (unsigned int i=0; i < v1.size(); ++i)
	{
		double diff = abs(v1[i] - v_expected[i]);
		if (diff>max_diff) max_diff = diff;
		if (diff>tolerance) ++fail_count;
	}
	if (max_diff<=tolerance)
	{
		cout << "\n\nlogspace\t\tPASSED" << endl;
	}
	else
	{
		cout << "\n\nlogspace\t\tFAILED\t"<< fail_count << "/" << v.size() << " failed elements\tMax Abs Difference = " << max_diff << endl;
	}

	// Test rlogspace
	/*
	void rlogspace(vector<double>& v, double pmin, double pmax, int Nv)
	{
		v.clear();
		v.resize(Nv);
		for (int i = Nv-1; i >= 0; i--){
			v[i] =  pow(10., pmin + (pmax - pmin) * i / double(Nv-1.)) ;
		}
	}
	*/


	// Test RHS_func
	/*
	// https://www.learncpp.com/cpp-tutorial/function-pointers/
	using RHS_func = double(*)(double, vector<double>&);
	//double findX(double LHS, RHS_func RHS, double a, double b, vector<double>& vars);

	// Modified False Position, Chapter 4.5 of Numerical Methods for Scientists and Engineers
	// Note, if root is complex and goes of to +inf, then findX will return close to b
	double findX(double LHS, RHS_func RHS, double a, double b, vector<double>& vars)
	{
		double fa = RHS(a, vars) - LHS;
		double fb = RHS(b, vars) - LHS;
		double fx, x;
		//int count = 0;
		
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
			//count++;

		} while (fabs(2*(a-b)/(a+b)) > 1.E-5);
		
		//cout << endl;
		//cout << count << ' ';

		return x;
	}
	*/


	// Test sample_pdf_idx
	/*
	// use the cdf to convert a uniform random number in [0,1] to match the corresponding pdf
	int sample_pdf_idx(mt19937& rng, vector<double>& cdf)
	{
		// pull sample from uniform distribution
		double u = uniform(rng, 0., 1.);

		// find index (iterator in this case) of the corresponding location in the cdf

		// Find the index such that cdf(idx-1) <= u <= cdf(idx)
		// If u = 0, use upper_bound, otherwise use lower_bound (both binary search algorithms, logN time)
		//  effectively, this inverts the cdf
		// This guarantees that *(idx_iter-1) <= u <= *(idx_iter) for all values of u in [0,1]
		return (u == 0. ? upper_bound(cdf.begin(), cdf.end(), u) : lower_bound(cdf.begin(), cdf.end(), u)) - cdf.begin();
	}
	
	// And Overloaded function

	// copy of above, but with returning u as well
	int sample_pdf_idx(mt19937& rng, vector<double>& cdf, double& u)
	{
		// pull sample from uniform distribution
		u = uniform(rng, 0., 1.);

		// find index (iterator in this case) of the corresponding location in the cdf

		// Find the index such that cdf(idx-1) <= u <= cdf(idx)
		// If u = 0, use upper_bound, otherwise use lower_bound (both binary search algorithms, logN time)
		//  effectively, this inverts the cdf
		// This guarantees that *(idx_iter-1) <= u <= *(idx_iter) for all values of u in [0,1]
		return (u == 0. ? upper_bound(cdf.begin(), cdf.end(), u) : lower_bound(cdf.begin(), cdf.end(), u)) - cdf.begin();
	}
	*/
}