#ifndef LUNAREJECTA_BATTLESHIPMONTECARLO_H
#define LUNAREJECTA_BATTLESHIPMONTECARLO_H

#include "LunarEjecta_params.h"

#include <list>

using namespace std;


// Will define a search and destroy scanner separately
// The search scanner will have a location at the center of the domain, with dx to fill the entire domain
// the destroy scanner will generate locations based on previous hits
struct radar_scanner
{
	list<vec3>::iterator   idx_loc;   // if active = 1, this is the iterator that triggered activation
	list<double>::iterator idx_dens;
	list<vec3>::iterator   idx_dx;
	list<uint>::iterator   idx_lifetime;

	list<vec3>   loc;      // location in phase space
	list<double> dens;     // probability phase space density in domain, shot weight = 1/dens
	list<vec3>   dx;       // Domain defined as: loc +/- dx
	list<uint>   lifetime; // lifetime in iterations of the scanner location, decrements each iteration

	bool active;       // 0 no, 1 yes. If active, the miss or hit tally will fall into its counting system
	long long uint miss_count;
	long long uint hit_count;

	double dens_tot;      // alpha or 1- alpha, the overall probablity of the radar scanner. All scanners should add to 1
	double lifetime_max;  // maximum lifetime set for new hit points
	double lifetime_rate; // value between 0 and 1, new loc lifetime based on old * this rate
	double dx_rate;       // value between 0 and 1, new loc dx based on old * this rate
};

void init_radar(radar_scanner &rs, double alpha, double lt_max, double lt_rate, double dx_rate)
{
	rs.active     = 0;
	rs.miss_count = 0;
	rs.hit_count  = 0;

	rs.dens_tot      = alpha;
	rs.lifetime_max  = lt_max;
	rs.lifetime_rate = lt_rate;
	rs.dx_rate       = dx_rate;
}


// get shot weight at location in phase space
// internal memory, store what would be the new lifetime
double get_shot_weight_and_location(vector<radar_scanner> &rs, vec3 &shot_loc);

// make the shot and tally the hit or miss
// if hit, add new hit position to population, returns eval_func at the shot_loc times the weight if a hit
// the hit_func takes the phase space location and decides whether it is a hit or not, returned as the bool
// https://stackoverflow.com/questions/2582161/function-pointer-as-parameter
// https://www.geeksforgeeks.org/passing-a-function-as-a-parameter-in-cpp/
double make_shot_and_tally(vector<radar_scanner> &rs, double weight, vec3 &shot_loc, bool (*hit_func)(vec3&), double (*eval_func)(vec3&));


#endif 