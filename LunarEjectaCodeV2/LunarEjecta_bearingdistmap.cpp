#include "LunarEjecta_bearingdistmap.h"
#include "LunarEjecta_params.h"

#include <ctime>
#include <stdlib.h>     /* srand, rand */
#include <cmath>
#include <vector>

using namespace std;


double get_lat_from_pole(double D, double bearing){
	return PI/2. - D;
}

double get_lon_from_pole(double D, double bearing){
	return bearing;
}

// D in units of radius, angles in radians
double get_lat_from_ROI(double D, double bearing, double lat_R, double lon_R){
	return asin(sin(lat_R)*cos(D) + cos(lat_R)*sin(D)*cos(bearing));
}

double get_lon_from_ROI(double D, double bearing, double lat_R, double lon_R){
	return lon_R + atan2(sin(bearing)*sin(D), cos(D)*cos(lat_R) - sin(D)*sin(lat_R)*cos(bearing));
}

// 0 deg = N, goes CCW, ranges [0, 2PI)
double initial_bearing(double lat1, double lat2, double lon1, double lon2)
{
	return fmod(atan2(sin(lon1 - lon2)*cos(lat2), cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(lon1-lon2)) + 2.*PI, 2.*PI);
}


// returns distance in units of radii
double distance_latlon(double lat1, double lat2, double lon1, double lon2)
{
	return 2.*asin(sqrt(sqr(sin((lat1-lat2)/2.)) + cos(lat1)*cos(lat2)*sqr(sin((lon1-lon2)/2.))));
}

void generate_lat_lon_points(vector<double> &lat, vector<double> &lon, double lat_center, double lon_center, double D0, double D1, double dlon)
{
	int loc_size = lat.size();

	double z, D, phi;

	for (int i = 0; i < loc_size; ++i)
	{
		z   = rand_uniform(2.*sqr(sin(D0/2.)), 2.*sqr(sin(D1/2.))); // units of radius
		D   = 2.*asin(sqrt(z/2.)); // units of radius
		phi = rand_uniform(-dlon/2. + lon_center, dlon/2. + lon_center); // bearing, lon from pole

		if (fabs(lat_center) == PI/2.)
		{
			lat[i] = get_lat_from_pole(D, phi);
			lon[i] = get_lon_from_pole(D, phi);
		}
		else
		{
			lat[i] = get_lat_from_ROI(D, phi, lat_center, lon_center);
			lon[i] = get_lon_from_ROI(D, phi, lat_center, lon_center);
		}

		//cout << lat[i] * 180./PI << ' ' << lon[i] * 180./PI << endl;
	}
}

// Note: need to make a switch for the longer way around, not just for the shorter way
void compute_mutual_bearings_distances(vector<double> &bearing_ROI, vector<double> &bearing_POI, vector<double> &D, vector<double> &lat_ROI, vector<double> &lon_ROI, vector<double> &lat_POI, vector<double> &lon_POI, int ejectaDistType)
{
	int i, j, idx;
	//ofstream file;
	//file.open("bearing_distance_numbers.txt");


	if (ejectaDistType == ejectaShort)
	{
		for (i = 0; i < lat_ROI.size(); ++i)
		{
			for (j = 0; j < lat_POI.size(); ++j)
			{
				idx = j + i*lat_POI.size();
				D[idx]           = distance_latlon(lat_ROI[i], lat_POI[j], lon_ROI[i], lon_POI[j]);
				bearing_ROI[idx] = initial_bearing(lat_ROI[i], lat_POI[j], lon_ROI[i], lon_POI[j]);
				bearing_POI[idx] = initial_bearing(lat_POI[j], lat_ROI[i], lon_POI[j], lon_ROI[i]);
			}
		}
	}
	else if(ejectaDistType == ejectaFar)
	{
		for (i = 0; i < lat_ROI.size(); ++i)
		{
			for (j = 0; j < lat_POI.size(); ++j)
			{
				idx = j + i*lat_POI.size();
				D[idx]           = PI - distance_latlon(lat_ROI[i], lat_POI[j], lon_ROI[i], lon_POI[j]);
				bearing_ROI[idx] = fmod(initial_bearing(lat_ROI[i], lat_POI[j], lon_ROI[i], lon_POI[j]) + PI, 2.*PI);
				bearing_POI[idx] = fmod(initial_bearing(lat_POI[j], lat_ROI[i], lon_POI[j], lon_ROI[i]) + PI, 2.*PI);
			}
		}
	}

	
	//file.close();
}



void init_bins(hist3DSet *b)
{ // bin[horizon angle][distance][bearing POI][bearing ROI]
	// note: the bearing_ROI_cum is missing the last N, needs to be added
	b->N_bearing_ROI_tot = (b->N_bearing_ROI_cum[b->N_horizon_ROI-1] + b->N_bearing_ROI[b->N_horizon_ROI-1]);
	b->N_bins =  b->N_D * b->N_bearing_POI * b->N_bearing_ROI_tot;
	cout << "Total number of bins = " << b->N_bins  << " " << b->N_bearing_ROI_cum[b->N_horizon_ROI-1] << ' ' << b->N_bearing_ROI[b->N_horizon_ROI-1] << endl;

	b->bin.resize(0);
	b->bin.resize(b->N_bins, 0.0);
}

int idx_bin(hist3DSet *b, int j_zenith, double D, double bearing_POI, double bearing_ROI)
{
	int idx_D     = (b->N_D - 1.) * (D - b->D_min) / (b->D_max - b->D_min);
	int idx_b_POI = (b->N_bearing_POI - 1.) * (bearing_POI - b->bearing_POI_min) / (b->bearing_POI_max - b->bearing_POI_min);
	int idx_b_ROI = (b->N_bearing_ROI[j_zenith] - 1.) * (bearing_ROI) / (2.*PI);

	// if(idx_D >= b->N_D || idx_b_POI >= b->N_bearing_POI || idx_b_ROI >= b->N_bearing_ROI[j_zenith]){
	// 	cout << "ERROR: " << idx_D << "/" << b->N_D << " | " << idx_b_POI << "/" << b->N_bearing_POI << " | " << idx_b_ROI << "/" << b->N_bearing_ROI[j_zenith] << endl;
	// }

	return idx_b_ROI + idx_b_POI * b->N_bearing_ROI[j_zenith] + idx_D * b->N_bearing_ROI[j_zenith] * b->N_bearing_POI + b->N_bearing_ROI_cum[j_zenith] * b->N_bearing_POI * b->N_D;
}

void tally_hist(hist3DSet *b, vector<double> &D, vector<double> &bearing_POI, vector<double> &bearing_ROI)
{
	int i, j, idx;
	cout << " tally hist starting...\n";
	for (i = 0; i < D.size(); ++i){
		for (j = 0; j < b->N_horizon_ROI; ++j){
			idx = idx_bin(b, j, D[i], bearing_POI[i], bearing_ROI[i]);
			// if(idx >= b->N_bins)
			// 	cout << idx << " / " << b->N_bins << endl;
			b->bin[idx] += 1.;
		}
	}

	// normalize bins to 1 over D, and POI/ROI bearings (not over zenith)
	for (i = 0; i < b->N_bins; ++i)
		b->bin[i] /= double(D.size());
}



hist3DSet* init_bearing_dist_map(input* p, int ejectaDistType)
{
	srand(time(0));
	hist3DSet* bearingdist = new hist3DSet;

	bearingdist->N_D              = p->N_D_perRegion;
	bearingdist->N_bearing_POI    = p->N_bearing_POI;
	bearingdist->N_latlon         = p->N_loc;
	bearingdist->N_horizon_ROI    = p->N_horizon_ROI;

	double dlat = PI / double(p->Nlat - 1.); // rad
	double dlon = 2.*PI / double(p->Nlon);   // rad

	double cur_lat = -PI/2. + dlat * (p->latlon_idx_min + p->latlon_idx_proc % p->Nlat);
	double cur_lon = dlon * (int(p->latlon_idx_min + p->latlon_idx_proc) / int(p->Nlat));

	get_lat_min_max(cur_lat, dlat, bearingdist->lat_min, bearingdist->lat_max);

	bearingdist->lon_min = cur_lon - dlon/2.;
	bearingdist->lon_max = cur_lon + dlon/2.;

	double D0 = PI/2. - bearingdist->lat_max; // units of lunar radii
	double D1 = PI/2. - bearingdist->lat_min; // units of lunar radii

	double SA_ratio = calcSA(cur_lat, bearingdist->lat_min, bearingdist->lat_max, dlat, dlon) / calcSA(PI/2., PI/2. - p->ROI_radius/p->lunar_radius, PI/2., p->ROI_radius/p->lunar_radius, 2.*PI);

	cout << "SA_ratio = " << SA_ratio << endl;

	// compute number of azimuth bins as a function of zenith index
	bearingdist->N_bearing_ROI.resize(p->N_horizon_ROI);
	bearingdist->N_bearing_ROI_cum.resize(p->N_horizon_ROI);

	compute_igloo_azm_bin_number(bearingdist->N_bearing_ROI, p->N_bearing_ROI);
	vCumLow(bearingdist->N_bearing_ROI, bearingdist->N_bearing_ROI_cum);

	// generate ROI points
	vector<double> lat_ROI(p->ROI_sample_points, 0.);
	vector<double> lon_ROI(p->ROI_sample_points, 0.);

	generate_lat_lon_points(lat_ROI, lon_ROI, p->ROI_lat, p->ROI_lon, 0.0, p->ROI_radius/p->lunar_radius, 2.*PI);

	// generate POI points
	int POI_sample_points = sqrt(SA_ratio) * p->ROI_sample_points;

	vector<double> lat_POI(POI_sample_points, 0.);
	vector<double> lon_POI(POI_sample_points, 0.);

	generate_lat_lon_points(lat_POI, lon_POI, PI/2., cur_lon, D0, D1, dlon);


	// compute mutual bearings and distances
	int mutual_points = p->ROI_sample_points * POI_sample_points;
	cout << "# of mutual points = " << mutual_points << endl;

	vector<double> D(mutual_points, 0.);
	vector<double> bearing_POI(mutual_points, 0.);
	vector<double> bearing_ROI(mutual_points, 0.);

	compute_mutual_bearings_distances(bearing_ROI, bearing_POI, D, lat_ROI, lon_ROI, lat_POI, lon_POI, ejectaDistType);


	// tally the histrogram bins

	// first, compute the min and max POI bearing. If the center is on 0 and 2PI, phase shift by 2PI so it's only around 0
	bearingdist->bearing_POI_min = vMin(bearing_POI);
	bearingdist->bearing_POI_max = vMax(bearing_POI);

	bearingdist->D_min = vMin(D);
	bearingdist->D_max = vMax(D);

	if (bearingdist->bearing_POI_max - bearingdist->bearing_POI_min > 0.99*2.*PI)
	{
		shiftAngle(bearing_POI, PI, 2.*PI);
		bearingdist->bearing_POI_min = vMin(bearing_POI);
		bearingdist->bearing_POI_max = vMax(bearing_POI);
	}

	cout << "bearing_POI_min/max = " << bearingdist->bearing_POI_min << " | " << bearingdist->bearing_POI_max << endl;


	// next, do the tallying
	init_bins(bearingdist);

	cout << "Average counts per bin = " << mutual_points * bearingdist->N_horizon_ROI / bearingdist->N_bins << endl;

	tally_hist(bearingdist, D, bearing_POI, bearing_ROI);


	return bearingdist;
}


