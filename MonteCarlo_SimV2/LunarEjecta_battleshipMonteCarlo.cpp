#include "LunarEjecta_params.h"
#include "LunarEjecta_asset.h"
#include "LunarEjecta_battleshipMonteCarlo.h"

#include <list>

using namespace std;

void initRadar(radar_scanner &rs, double alpha, int lt_max, double lt_rate, double dx_rate)
{
	rs.active     = 0;
	rs.miss_count = 0;
	rs.hit_count  = 0;

	rs.dens_tot      = alpha;
	rs.lifetime_max  = lt_max;
	rs.lifetime_rate = lt_rate;
	rs.dx_rate       = dx_rate;

	cout << "Initialize radar scanner:\n";
	cout << " percentage = " << rs.dens_tot << " | max lifetime = " << rs.lifetime_max << " | lifetime rate = " <<  rs.lifetime_rate << " | dx rate = " << rs.dx_rate << endl;
}


// get shot weight at location in phase space, weight is returned
// internal memory, store what would be the new lifetime
double getShotWeightAndLocation(radar_scanner &rs_s, radar_scanner &rs_d, vec3 &shot_ph_loc)
{
	return 0.;
}


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
              bool (*hit_func)(vec3&, vec3&, asset&)) //, double (*eval_func)(vec3&)
{
	return 0;
}