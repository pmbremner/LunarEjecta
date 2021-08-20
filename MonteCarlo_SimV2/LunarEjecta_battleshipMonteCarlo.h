#ifndef LUNAREJECTA_BATTLESHIPMONTECARLO_H
#define LUNAREJECTA_BATTLESHIPMONTECARLO_H

#include "LunarEjecta_params.h"
#include "LunarEjecta_asset.h"

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
	list<int>::iterator    idx_lifetime;

	list<vec3>   ph_loc;   // location in phase space
	list<double> dens;     // probability phase space density in domain, shot weight = 1/dens
	list<vec3>   dx;       // Domain defined as: loc +/- dx
	list<int>    lifetime; // lifetime in iterations of the scanner location, decrements each iteration

	bool active;       // 0 no, 1 yes. If active, the miss or hit tally will fall into its counting system
	long long unsigned miss_count;
	long long unsigned hit_count;

	double dens_tot;      // alpha or 1- alpha, the overall probablity of the radar scanner. All scanners should add to 1
	int lifetime_max;  // maximum lifetime set for new hit points
	double lifetime_rate; // value between 0 and 1, new loc lifetime based on old * this rate
	double dx_rate;       // value between 0 and 1, new loc dx based on old * this rate
};

void initRadar(radar_scanner &rs, double alpha, int lt_max, double lt_rate, double dx_rate);

//void insertScanPoint(radar_scanner &rs, vec3 &ph_loc, )


// get shot weight at location in phase space, weight is returned
double getShotWeightAndLocation(radar_scanner &rs_s, radar_scanner &rs_d, vec3 &shot_ph_loc);

// make the shot and tally the hit or miss, if hit, 1 is returned, if miss, 0 is returned
// if hit, add new hit position to population
// the hit_func takes the phase space and physical location and decides whether it is a hit or not, returned as the bool
// https://stackoverflow.com/questions/2582161/function-pointer-as-parameter
// https://www.geeksforgeeks.org/passing-a-function-as-a-parameter-in-cpp/
bool makeShot(radar_scanner &rs_s,
	          radar_scanner &rs_d,
              double weight,
              vec3 &shot_ph_loc0, // initial phase space location
              vec3 &shot_pos0,    // initial physical space location
              vec3 &shot_ph_loc1, // final phase space location
              vec3 &shot_pos1,    // final physical space location
              bool (*hit_func)(vec3&, vec3&, asset&)); //, double (*eval_func)(vec3&)


#endif 