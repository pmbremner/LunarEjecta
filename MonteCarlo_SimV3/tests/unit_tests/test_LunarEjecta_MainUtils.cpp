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

	// Informational Message
	string intro_msg = "Testing LunarEjecta_MainUtils...";
	cout << intro_msg << endl;
	

	// Setup a tolerance for aritmitic operations
	double tolerance = 1e-10;


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
		if (diff<=tolerance)
		{
			cout << "vMax\t\tPASSED" << endl;
		}
		else
		{
			cout << "vMax\t\tFAILED\tDifference = " << diff << endl;
		}


		// Test vMin
		diff = vMin(test_bundle.fp_vec) - test_bundle.expected_vMin;
		if (diff<=tolerance)
		{
			cout << "vMin\t\tPASSED" << endl;
		}
		else
		{
			cout << "vMin\t\tFAILED\tDifference = " << diff << endl;
		}


		// Test vSum
		diff = vSum(test_bundle.fp_vec) - test_bundle.expected_vSum;
		if (diff<=tolerance)
		{
			cout << "vSum\t\tPASSED" << endl;
		}
		else
		{
			cout << "vSum\t\tFAILED\tDifference = " << diff << endl;
		}
	}
	
	

	// Test linspace
	/*
	void linspace(vector<double>& v, double vmin, double vmax, int Nv)
	{
		v.clear();
		v.resize(Nv);
		//cout << endl <<  vmin << ' ' << vmax << endl;
		for (int i = 0; i < Nv; ++i){
			v[i] = vmin + (vmax - vmin) * i / double(Nv-1.);
			//cout << v[i] << endl;
		}
	}
	*/


	// Test logspace
	/*
	void logspace(vector<double>& v, double pmin, double pmax, int Nv)
	{
		v.clear();
		v.resize(Nv);
		//cout << endl <<  pmin << ' ' << pmax << endl;
		for (int i = 0; i < Nv; ++i){
			v[i] = pow(10., pmin + (pmax - pmin) * i / double(Nv-1.)) ;
			//cout << v[i] << endl;
		}
	}
	*/


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