#include "LunarEjecta_params.h"
#include "LunarEjecta_battleshipMonteCarloV2.h"

#include <list>
#include <algorithm> // min and max
#include <chrono>
#include <random>

using namespace std;



// struct scan
// {
// 	vector<double> ph;   // phase space center location
// 	vector<double> dph;  // phase space width
// 	double         prob_dens; // probability density in the phase space domain
// 	int            lifetime;  // current lifetime of scan, in # of iterations
// 	int            lt_dec;    // decrement amount applied to lifetime for each new iteration
// 							  // if zero, the scan will live during the entire simulation
// 	int            generation; // # of generation, base gen is 0 (which will be the search method, lt_dec = 0 usually)
// };


// struct radar_scanner
// {
// 	list<scan>           scan_props; // list of scan properties
// 	list<scan>::iterator idx_scan;   // current scan sampled from
// 	vector<double>       ph_scan;    // current scan phase space location
// 	double               net_prob_dens; // net probability density of ph_scan

// 	vector<long long unsigned> miss_count; // miss count for each generation
// 	vector<long long unsigned> hit_count;  // hit count for each generation

// 	double gen_zero_prob; // the alpha value, all other gens get 1-alpha 
// 	double lifetime_rate; // value between 0 and 1, new loc lifetime based on old * this rate, if from gen zero, just inherit the lt
// 	double dx_rate;       // value between 0 and 1, new loc dx based on old * this rate
// }
double errDens(double a, double b)
{
	return 0.;
}


double constDens(double a, double b)
{
	return b - a;
}

double zenithDens(double a, double b)
{

	if (a < 0. && b < 0. )
	{
		return cos(b) - cos(a);
	}
	else if (a > 0. && b > 0.)
	{
		return cos(a) - cos(b);
	}
	else // a <= 0 and b >= 0
	{
		return 2. - cos(a) - cos(b);
	}
}


// HRO = human readible output
void printScan(scan &s, bool HRO)
{
	if (HRO)
	{
		cout << "Generation " << s.generation << " | Current lifetime = " << s.lifetime << " | Probability density = " << s.prob_dens << endl;
		cout << " Phase space center = ";
		for (int i = 0; i < s.ph.size(); ++i)
			cout << s.ph[i] << " ";
		cout << endl;
		cout << " Phase space width = ";
		for (int i = 0; i < s.dph.size(); ++i)
			cout << s.dph[i] << " ";
		cout << endl;
	}
	// else
	// {

	// }
}


double getLocalProbDens(radar_scanner &rs, vector<double> &ph_min, vector<double> &ph_max)
{
	double dens = 1.;
	for (int i = 0; i < ph_min.size(); ++i){

		cout << ph_min[i] << ' ' << ph_max[i] << " | " << rs.dens_func[i](ph_min[i], ph_max[i]) << endl;

		dens *= 1. / rs.dens_func[i](ph_min[i], ph_max[i]);
	}

	return dens;
}

void insertScan(radar_scanner &rs, vector<double> &ph, vector<double> &dph, int lt, int lt_dec, int gen, double rel_prob)
{
	scan newScan;
	newScan.ph         = ph;
	newScan.dph        = dph;

	for (int i = 0; i < dph.size(); ++i)
	{
		newScan.ph_min.push_back(max(rs.ph_min[i], ph[i] - dph[i]/2.));
		newScan.ph_max.push_back(min(rs.ph_max[i], ph[i] + dph[i]/2.));
	}

	newScan.prob_dens  = getLocalProbDens(rs, newScan.ph_min, newScan.ph_max) * rel_prob;
	newScan.lifetime   = lt;
	newScan.lt_dec     = lt_dec;
	newScan.generation = gen;

	rs.scan_props.push_back(newScan);
}

void setupDensFunc(radar_scanner &rs, vector<int> &dens_func_id)
{
	for (int i = 0; i < dens_func_id.size(); ++i)
	{
		if (dens_func_id[i] == 0)
			rs.dens_func.push_back(constDens);
		else if (dens_func_id[i] == 1)
			rs.dens_func.push_back(zenithDens);
		else
		{
			cerr << "ERROR: invalid density function ID\n";
			rs.dens_func.push_back(errDens);
		}
	}
}

void setupPhaseSpaceBoundary(radar_scanner &rs, vector<double> &ph, vector<double> &dph)
{
	for (int i = 0; i < ph.size(); ++i)
	{
		rs.ph_min.push_back(ph[i] - dph[i]/2.);
		rs.ph_max.push_back(ph[i] + dph[i]/2.);
	}
}


// Init the radar, the zeroth gen is started with no decrement in lifetime, spanning the entire domain
void initRadar(radar_scanner &rs, double alpha, int lt_max, double lt_rate, double dx_rate, vector<double> &ph, vector<double> &dph, vector<int> &dens_func_id)
{
	// sets the density function pointers depending on the IDs supplied
	setupDensFunc(rs, dens_func_id);

	// set up min max phase space
	setupPhaseSpaceBoundary(rs, ph, dph);

	// insert the search mode, no decrement and zeroth gen
	insertScan(rs, ph, dph, lt_max, 0, 0, alpha);

	rs.idx_search = rs.scan_props.begin();

	// init the miss and hit counters for gen 0
	rs.miss_count.push_back(0);
	rs.hit_count.push_back(0);

	rs.gen_zero_prob = alpha;
	rs.lifetime_rate = lt_rate;
	rs.dx_rate       = dx_rate;

	cout << "Initialize radar scanner:\n";
	cout << " percentage = " << rs.gen_zero_prob << " | lifetime rate = " <<  rs.lifetime_rate << " | dx rate = " << rs.dx_rate << endl;

	printScan(*(rs.scan_props.begin()), 1);

	//init the random generator
	// https://stackoverflow.com/questions/24334012/best-way-to-seed-mt19937-64-for-monte-carlo-simulations
	rs.seed = random_device{}() * chrono::system_clock::now().time_since_epoch().count();
	rs.rng.seed(rs.seed);

	cout << " RND seed = " << rs.seed << endl;
}

inline double uniform(mt19937& rng, double a, double b){
	return a + (b - a) * rng() / double(rng.max());
}


// the radar gets a scan sample from all potential scan locations, returns the weight of the scan
double getSampleScan(radar_scanner &rs, vector<double> &ph_scan)
{
	double weight = 0.;

	int Nscan = rs.scan_props.size(), i;
	int Ndim = (*rs.idx_search).ph_min.size();
	//cout << "Nscan = " << Nscan << endl;

	double uni_val = uniform(rs.rng, 0., 1.);

	// clear the phase space scan so we can fill it with a new scan
	rs.ph_scan.clear();

	if (uni_val < rs.gen_zero_prob || Nscan == 1)
	{
		// pull from unifrom distribution for each phase space dimension
		for (i = 0; i < Ndim; ++i)
			rs.ph_scan.push_back( uniform(rs.rng, (*rs.idx_search).ph_min[i], (*rs.idx_search).ph_max[i]) );

		// set the scan index
		rs.idx_scan = rs.idx_search;

		// compute net density, dividing out the alpha if Nscan = 1

	}
	else
	{
		// find a scan to pull from

		// pull from unifrom distribution for each phase space dimension

		// set scan index

		// compute net density

	}

	return weight; // = 1 / net density
}

// The hit or miss is tallied in generation of the current scan
// If hit = 1, a new scan is generated at the scan spot based on the previos scan props
// Also, loop through all scans to apply lt_dec, and remove scans with lifetime = 0
void tallyScan(radar_scanner &rs, bool hit)
{

}