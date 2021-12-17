#include "LunarEjecta_params.h"
#include "LunarEjecta_battleshipMonteCarloV2.h"

#include <list>
#include <algorithm> // min and max
#include <chrono>
#include <random>
#include <math.h>       /* floor */

using namespace std;


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

void printHitMissReport(radar_scanner &rs)
{
	long long unsigned int hit_tot = 0, miss_tot = 0;
	cout << "*************************\n";
	cout << "**** Hit-miss Report ****\n";
	for (int i = 0; i < rs.miss_count.size(); ++i)
	{
		hit_tot += rs.hit_count[i];
		miss_tot += rs.miss_count[i];

		cout << "Generation  = " << i << endl;
		cout << " Hit count  = " << rs.hit_count[i] << endl;
		cout << " Miss count = " << rs.miss_count[i] << endl;
		cout << " Hit % = " << 100. * rs.hit_count[i] / double(rs.hit_count[i] + rs.miss_count[i]) << endl;
	}
	cout << "*************************\n";
	cout << " Total hit count  = " << hit_tot << endl;
	cout << " Total miss count = " << miss_tot << endl;
	cout << " Total hit % = " << 100. * hit_tot / double(hit_tot + miss_tot) << endl;
	cout << " 1 out of every " << double(hit_tot + miss_tot) / double(hit_tot);
}


// HRO = human readible output
void printScan(scan &s, bool HRO)
{	cout << " | " << &s <<  endl;
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
		cout <<  endl;
	}
	// else
	// {

	// }
}


double getLocalProbDens(radar_scanner &rs, vector<double> &ph_min, vector<double> &ph_max)
{
	double dens = 1.;
	for (int i = 0; i < ph_min.size(); ++i){

		///cout << ph_min[i] << ' ' << ph_max[i] << " | " << rs.dens_func[i](ph_min[i], ph_max[i]) << endl;

		dens *= 1. / rs.dens_func[i](ph_min[i], ph_max[i]);
	}

	return dens;
}


double getNetProbDens(radar_scanner &rs)
{
	double dens = 0.;
	bool flag_domain;
	int Nscan = rs.scan_props.size(), j;
	int Ndim = rs.ph_scan.size();

	list<scan>::iterator itr;

	for (itr = rs.scan_props.begin(); itr !=  rs.scan_props.end(); ++itr)
	{
		// check if scan point is in the domain of the i-th listed scan, jump out early if already outside domain in first dims
		flag_domain = 1;
		for (j = 0; j < Ndim && flag_domain == 1; ++j){
			//cout << rs.ph_scan[j] << ' ' << (*itr).ph_min[j] << ' ' << (*itr).ph_max[j] << endl;
			flag_domain = (rs.ph_scan[j] >= (*itr).ph_min[j] && rs.ph_scan[j] <= (*itr).ph_max[j] ? 1 : 0);
		}

		// if in domain, count the probability density
		if (flag_domain)
			dens += (*itr).prob_dens / (itr != rs.scan_props.begin() ? double(Nscan-1.) : 1.);

		//cout << "dens = " << dens << endl;
	}
	return dens;
}


void insertScan(radar_scanner &rs, vector<double> &ph, vector<double> &dph, int lt, int lt_dec, int gen, double rel_prob)
{
	scan newScan;
	newScan.ph         = ph;
	newScan.dph        = dph;

	// if the new gen is further down the line than the deepest-gen-so-far, increase the Ngen
	// also increase the tally size
	if (gen > rs.Ngen-1){
		rs.Ngen++;
		rs.miss_count.push_back(0);
		rs.hit_count.push_back(0);
	}

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
void initRadar(radar_scanner &rs, int N_max, double alpha, int lt_max, double lt_rate, double dx_rate, vector<double> &ph, vector<double> &dph, vector<int> &dens_func_id)
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

	rs.Ngen = 1; // start with gen 0

	rs.N_max         = N_max;
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

// returns an int from a to b-1
inline int uniformInt(mt19937& rng, int a, int b){
	return int(floor(a + (b - a) * rng() / double(rng.max()+1.)));
}

void uniformLatLon(mt19937& rng, vector<double> &lat_lon, vector<double> &lat_lon_cart, double r, double lat_min, double lat_max, double lon_min, double lon_max)
{
	// force an empty vector
	lat_lon.clear();
	lat_lon_cart.clear();

	// generate sampled lat over range
	double U_min = (cos(lat_max + PI/2.) + 1.) / 2.;
	double U_max = (cos(lat_min + PI/2.) + 1.) / 2.;

	lat_lon.push_back( acos( 2.*uniform(rng, U_min, U_max) - 1. ) - PI/2.); // lat point, rad
	lat_lon.push_back( uniform(rng, lon_min, lon_max) ); // lon point, rad

	lat_lon_cart.push_back( r * sin(PI/2. - lat_lon[0]) * cos(lat_lon[1]) ); // x
	lat_lon_cart.push_back( r * sin(PI/2. - lat_lon[0]) * sin(lat_lon[1]) ); // y
	lat_lon_cart.push_back( r * cos(PI/2. - lat_lon[0]) );                   // z
}

// THIS FUNCTION IS NOT WORKING, using advance instead
// N (element #) cannot be larger than Ntot-1
// void getIter(list<scan> &ls, list<scan>::iterator &itr, int N)
// {
// 	int Ntot = ls.size(), i;
// 	if (N < (Ntot-1)/2) // start from beginning
// 	{
// 		itr = ls.begin();
// 		cout << &itr << endl;
// 		for (i = 0; i < N; ++i){ // replace with the function advance***
// 			itr++;
// 			cout << &itr << endl;
// 		}
// 	}	
// 	else // start from end
// 	{
// 		itr = ls.end();
// 		cout << &itr << endl;
// 		for (i = 0; i < Ntot-(N+1); ++i){
// 			itr--;
// 			cout << &itr << endl;
// 		}

// 	}
// }


// the radar gets a scan sample from all potential scan locations, returns the weight of the scan
double getSampleScan(radar_scanner &rs, vector<double> &ph_scan)
{
	double dens = 0.;
	int idx;

	int Nscan = rs.scan_props.size(), i;
	int Ndim = (*rs.idx_search).ph_min.size();
	///cout << "Nscan = " << Nscan << endl;

	double uni_val = uniform(rs.rng, 0., 1.);

	// clear the phase space scan so we can fill it with a new scan
	rs.ph_scan.clear();

	//cout << uni_val << " | " << rs.gen_zero_prob << " | " << Nscan << endl;

	if (uni_val < rs.gen_zero_prob || Nscan == 1)
	{
		//cout << "search:\n";
		// pull from unifrom distribution for each phase space dimension
		for (i = 0; i < Ndim; ++i)
			rs.ph_scan.push_back( uniform(rs.rng, (*rs.idx_search).ph_min[i], (*rs.idx_search).ph_max[i]) );

		// set the scan index
		rs.idx_scan = rs.idx_search;

		///printScan(*rs.idx_scan, 1);

		// compute net density, dividing out the alpha if Nscan = 1
		dens = getNetProbDens(rs) * (Nscan == 1 ? 1./rs.gen_zero_prob : 1.);
	}
	else
	{
		//cout << "destroy:\n";
		// find a scan to pull from
		idx = uniformInt(rs.rng, 1, Nscan);
		///cout << "***Scan index: " << idx << endl;

		// iterate to the idx-th position in the scan list
		rs.idx_scan = rs.idx_search;
		advance(rs.idx_scan, idx);
		//getIter(rs.scan_props, rs.idx_scan, idx);

		///printScan(*rs.idx_scan, 1);

		// pull from unifrom distribution for each phase space dimension
		for (i = 0; i < Ndim; ++i)
			rs.ph_scan.push_back( uniform(rs.rng, (*rs.idx_scan).ph_min[i], (*rs.idx_scan).ph_max[i]) );


		// compute net density
		dens = getNetProbDens(rs);
	}

	ph_scan = rs.ph_scan; // for return

	return 1. / dens; // = 1 / net density
}


void ageScan(list<scan> &ls)
{
	for (list<scan>::iterator i = ls.begin(); i != ls.end(); i++)
	{
		///printScan((*i), 1);

		(*i).lifetime -= (*i).lt_dec;

		if ((*i).lifetime <= 0)
		{
			i = ls.erase(i); // now is at the next, need to step back
			i--;
		}
	}
}

// The hit or miss is tallied in generation of the current scan
// If hit = 1, a new scan is generated at the scan spot based on the previous scan props
// Also, loop through all scans to apply lt_dec, and remove scans with lifetime = 0
void tallyScan(radar_scanner &rs, bool hit)
{
	int cur_gen = (*rs.idx_scan).generation;
	double total_scans_in_gen = rs.hit_count[cur_gen] + rs.miss_count[cur_gen];
	///cout << "current gen = " << cur_gen << endl;
	if (hit)
	{

		//give extra life if the hit rate is near 50%, randomly
		if (rs.idx_scan != rs.idx_search && total_scans_in_gen > 0 && 1. - 2.*fabs(rs.hit_count[cur_gen] / total_scans_in_gen - 0.50) >= uniform(rs.rng, 0., 1.) )
			(*rs.idx_scan).lifetime += 2. * rs.max_lifetime / 200.;

		// if too densly populated, remove life
		if (rs.idx_scan != rs.idx_search && log10((*rs.idx_search).prob_dens / rs.net_prob_dens) <= uniform(rs.rng, -5., 0.) )
			(*rs.idx_scan).lifetime -= 5. * rs.max_lifetime / 200.;
		// // else if(rs.idx_scan != rs.idx_search)
		// //  	(*rs.idx_scan).lifetime += 1;


		// tally the hit in the current gen
		rs.hit_count[cur_gen]++;

		// create new scan point
		// set up dph, dfactor gets a boast going from gen 0 to 1, helps with making dph small enough to be efficient
		vector<double> dph = (*rs.idx_scan).dph;
		//double dfactor = (cur_gen == 0 ? pow(rs.dx_rate, 2) : rs.dx_rate);
		double dfactor =  rs.dx_rate;

		for (int i = 0; i < dph.size(); ++i)
			dph[i] *= dfactor;

		// setup the lifetime, min lifetime is 10 by default
		//int lt = max((*rs.idx_scan).lifetime * rs.lifetime_rate, 10);

		int lt = (*rs.idx_scan).lifetime * rs.lifetime_rate + 1;

		if (lt > 5)
			insertScan(rs, rs.ph_scan, dph, lt, 1, cur_gen + 1, 1. - rs.gen_zero_prob);

		///printScan(*(--rs.scan_props.end()), 1);

	}
	else // miss
	{
		// tally the miss in the current gen
		rs.miss_count[cur_gen]++;
	}

	// apply lifetime decrement and remove scans with lifetime = 0
	ageScan(rs.scan_props);

}