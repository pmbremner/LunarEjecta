#ifndef LUNAREJECTA_AZM_DIST_MAP_H
#define LUNAREJECTA_AZM_DIST_MAP_H

#include "LunarEjecta_params.h"

using namespace std;

struct hist2DSet
{
	int N_D_perRegion;
	int N_azm_perRegion;

	vector<vector<double>> D_bin_edge;
	vector<vector<double>> D_bin_center;
	vector<vector<double>> azm_bin_edge;
	vector<vector<double>> azm_bin_center;
	vector<vector<double>> bin_weight;
};

double get_lat_from_pole(double D, double bearing);
double get_lon_from_pole(double D, double bearing);
double get_lat_from_ROI(double D, double bearing, double lat_R, double lon_R);
double get_lon_from_ROI(double D, double bearing, double lat_R, double lon_R);

hist2DSet* init_azm_dist_map(input* params);


#endif 