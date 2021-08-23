#ifndef LUNAREJECTA_BATTLESHIPMONTECARLO_H
#define LUNAREJECTA_BATTLESHIPMONTECARLO_H

#include "LunarEjecta_params.h"

#include <list>
#include <vector>
#include <chrono>
#include <random>

using namespace std;

typedef double (*dfunc)(double,double);

//template <typename T> 
struct scan
{
	vector<double> ph;   // phase space center location
	vector<double> dph;  // phase space width
	vector<double> ph_min;  // phase space min, for this scan
	vector<double> ph_max;  // phase space max, for this scan
	double         prob_dens; // probability density in the phase space domain
	int            lifetime;  // current lifetime of scan, in # of iterations
	int            lt_dec;    // decrement amount applied to lifetime for each new iteration
							  // if zero, the scan will live during the entire simulation
	int            generation; // # of generation, base gen is 0 (which will be the search method, lt_dec = 0 usually)
};


struct radar_scanner
{
	list<scan>           scan_props; // list of scan properties
	list<scan>::iterator idx_search;   // search scan, DO NOT MODIFY

	list<scan>::iterator idx_scan;   // current scan sampled from
	vector<double>       ph_scan;    // current scan phase space location
	double               net_prob_dens; // net probability density of ph_scan

	vector<dfunc>  dens_func; // pointers to functions to compute the prob dens for each dimension of phase space  
	vector<double> ph_min; // minimum boundary of phase space
	vector<double> ph_max; // maximum boundary of phase space

	int Ngen; // number of generations
	vector<long long unsigned> miss_count; // miss count for each generation
	vector<long long unsigned> hit_count;  // hit count for each generation

	double gen_zero_prob; // the alpha value, all other gens get 1-alpha 
	double lifetime_rate; // value between 0 and 1, new loc lifetime based on old * this rate, if from gen zero, just inherit the lt
	double dx_rate;       // value between 0 and 1, new loc dx based on old * this rate

	unsigned seed; 
	mt19937 rng;

};

double errDens(double a, double b);
double constDens(double a, double b);
double zenithDens(double a, double b);


void printHitMissReport(radar_scanner &rs);
void printScan(scan &s, bool HRO);

double getLocalProbDens(radar_scanner &rs, vector<double> &ph, vector<double> &dph);

void insertScan(radar_scanner &rs, vector<double> &ph, vector<double> &dph, int lt, int lt_dec, int gen, double rel_prob);

void setupDensFunc(radar_scanner &rs, vector<int> &dens_func_id);

void setupPhaseSpaceBoundary(radar_scanner &rs, vector<double> &ph, vector<double> &dph);


// Init the radar, the zeroth gen is started with no decrement in lifetime, spanning the entire domain
void initRadar(radar_scanner &rs, double alpha, int lt_max, double lt_rate, double dx_rate, vector<double> &ph, vector<double> &dph, vector<int> &dens_func_id);

// the radar gets a scan sample from all potential scan locations, returns the weight of the scan
double getSampleScan(radar_scanner &rs, vector<double> &ph_scan);

// The hit or miss is tallied in generation of the current scan
// If hit = 1, a new scan is generated at the scan spot based on the previos scan props
// Also, loop through all scans to apply lt_dec, and remove scans with lifetime = 0
void tallyScan(radar_scanner &rs, bool hit);


#endif 