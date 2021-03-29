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

		cout << lat[i] * 180./PI << ' ' << lon[i] * 180./PI << endl;
	}
}




hist3DSet* init_bearing_dist_map(input* p)
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

	bearingdist->D_min = PI/2. - bearingdist->lat_max; // units of lunar radii
	bearingdist->D_max = PI/2. - bearingdist->lat_min; // units of lunar radii

	double SA_ratio = calcSA(cur_lat, bearingdist->lat_min, bearingdist->lat_max, dlat, dlon) / calcSA(PI/2., PI/2. - p->ROI_radius/p->lunar_radius, PI/2., p->ROI_radius/p->lunar_radius, 2.*PI);

	cout << "SA_ratio = " << SA_ratio << endl;

	// generate ROI points
	vector<double> lat_ROI(p->ROI_sample_points, 0.);
	vector<double> lon_ROI(p->ROI_sample_points, 0.);

	generate_lat_lon_points(lat_ROI, lon_ROI, p->ROI_lat, p->ROI_lon, 0.0, p->ROI_radius/p->lunar_radius, 2.*PI);

	// generate POI points
	int POI_sample_points = sqrt(SA_ratio) * p->ROI_sample_points;

	vector<double> lat_POI(POI_sample_points, 0.);
	vector<double> lon_POI(POI_sample_points, 0.);

	generate_lat_lon_points(lat_POI, lon_POI, PI/2., cur_lon, bearingdist->D_min, bearingdist->D_max, dlon);


	// compute mutual bearings and distances
	int mutual_points = p->ROI_sample_points * POI_sample_points;

	vector<double> D(mutual_points, 0.);
	vector<double> bearing_POI(mutual_points, 0.);
	vector<double> bearing_ROI(mutual_points, 0.);

	compute_mutual_bearings_distances(bearing_ROI, bearing_POI, D, lat_ROI, lon_ROI, lat_POI, lon_POI);


	// talley the histrogram bins

	// first, compute the min and max POI bearing. If the center is on 0 and 2PI, phase shift by 2PI so it's only around 0
	bearingdist->bearing_POI_min = vMin(bearing_POI);
	bearingdist->bearing_POI_max = vMax(bearing_POI);

	if (bearingdist->bearing_POI_max - bearingdist->bearing_POI_min > 0.99*2.*PI)
	{
		shiftAngle(bearing_POI, PI, 2.*PI);
		bearingdist->bearing_POI_min = vMin(bearing_POI);
		bearingdist->bearing_POI_max = vMax(bearing_POI);
	}

	return bearingdist;
}

// hist2DSet* init_azm_dist_map(input* p)
// {
// 	srand(time(0));
// 	hist2DSet* azmdist_map = new hist2DSet;

// 	int i, j, idx, N_sample;
// 	double z, D, phi, D0, D1, a;
// 	azmdist_map->N_D_perRegion    = p->N_D_perRegion;
// 	azmdist_map->N_azm_perRegion  = p->N_azm_perRegion;

// 	azmdist_map->D_bin_edge.resize(p->N_loc);
// 	azmdist_map->D_bin_center.resize(p->N_loc);
// 	azmdist_map->azm_bin_edge.resize(p->N_loc);
// 	azmdist_map->azm_bin_center.resize(p->N_loc);
// 	azmdist_map->bin_weight.resize(p->N_loc);

// 	double cur_lat, cur_lon, dlat, dlon, lat_min, lat_max, lon_min, lon_max;
// 	double d_min, d_max, ang_min, ang_max, d_D, d_ang;

// 	dlat = PI / double(p->Nlat - 1.); // rad
// 	dlon = 2.*PI / double(p->Nlon);   // rad

// 	vector<double> p_ROI_lat;
// 	vector<double> p_ROI_lon;
// 	vector<double> p_arc_lat;
// 	vector<double> p_arc_lon;

// 	p_ROI_lat.resize(p->ROI_sample_points);
// 	p_ROI_lon.resize(p->ROI_sample_points);
// 	p_arc_lat.resize(p->arc_sample_points);
// 	p_arc_lon.resize(p->arc_sample_points);

// 	// loop through each location
// 	for (p->latlon_idx_proc = 0; p->latlon_idx_proc < p->N_loc; p->latlon_idx_proc++)
// 	{
// 		cur_lat = -PI/2. + dlat * (p->latlon_idx_min + p->latlon_idx_proc % p->Nlat);
// 		cur_lon = dlon * (int(p->latlon_idx_min + p->latlon_idx_proc) / int(p->Nlat));

// 		get_lat_min_max(cur_lat, dlat, lat_min, lat_max);

// 		D0 = PI/2. - lat_max; // units of lunar radii
// 		D1 = PI/2. - lat_min; // units of lunar radii

// 		cout << "azm_dist map: lat/lon = " << cur_lat*180./PI << " | " << cur_lon*180./PI << endl;
	
// 		// generate lat & lon points in ROI
// 		for (i = 0; i < p->ROI_sample_points; ++i)
// 		{
// 			z   = rand_uniform(0., 2.*sqr(sin(p->ROI_radius/(2.*p->lunar_radius))) );
// 			D   = 2.*asin(sqrt(z/2.)); // lunar radii
// 			phi = rand_uniform(0., 2.*PI);

// 			p_ROI_lat[i] = get_lat_from_ROI(D, phi, p->ROI_lat, p->ROI_lon);
// 			p_ROI_lon[i] = get_lon_from_ROI(D, phi, p->ROI_lat, p->ROI_lon);

// 			//cout << i << " " << p_ROI_lon[i]*180./PI << " " << p_ROI_lat[i]*180./PI << endl;
// 		}


// 		// compute lat and lon points in arc
// 		for (i = 0; i < p->arc_sample_points; ++i)
// 		{
// 			z   = rand_uniform(2.*sqr(sin(D0/2.)), 2.*sqr(sin(D1/2.))); // units of radius
// 			D   = 2.*asin(sqrt(z/2.)); // units of radius
// 			phi = rand_uniform(-dlon/2. + cur_lon, dlon/2. + cur_lon); // bearing

// 			p_arc_lat[i] = get_lat_from_pole(D, phi);
// 			p_arc_lon[i] = get_lon_from_pole(D, phi);

// 			//cout << i << " " << p_arc_lon[i]*180./PI << " " << p_arc_lat[i]*180./PI << endl;
// 		}

// 		// compute all mutual distances
// 		vector<double> d_array;
// 		vector<double> ang_array, POI_ang_array;

// 		N_sample = p->ROI_sample_points * p->arc_sample_points;
// 		d_array.resize(N_sample);
// 		ang_array.resize(N_sample);
// 		POI_ang_array.resize(N_sample);

// 		d_min = 2.*PI;
// 		d_max = 0.0;
// 		ang_min = 2.*PI;
// 		ang_max = 0.0;

// 		for (i = 0; i < p->ROI_sample_points; ++i)
// 		{
// 			for (j = 0; j < p->arc_sample_points; ++j)
// 			{
// 				idx = j + i*p->arc_sample_points;

// 				dlat = p_ROI_lat[i] - p_arc_lat[j]; // redefining here
// 				dlon = p_ROI_lon[i] - p_arc_lon[j];

// 				a = sqr(sin(dlat/2.)) + cos(p_ROI_lat[i])*cos(p_arc_lat[j])*sqr(sin(dlon/2.));

// 				d_array[idx] = 2.*asin(sqrt(a));

// 				// local azm from ROI to POI
// 				ang_array[idx] = initial_bearing(p_ROI_lat[i], p_arc_lat[j], p_ROI_lon[i], p_arc_lon[j]);
// 				//ang_array[idx] = fmod(atan2(sin(dlon)*cos(p_arc_lat[j]), cos(p_ROI_lat[i])*sin(p_arc_lat[j]) - sin(p_ROI_lat[i])*cos(p_arc_lat[j])*cos(dlon)) + 2.*PI, 2.*PI);
		
// 				POI_ang_array[idx] = initial_bearing(p_arc_lat[j], p_ROI_lat[i], p_arc_lon[j], p_ROI_lon[i]);
// 				//POI_ang_array[idx] = fmod(atan2(-sin(dlon)*cos(p_ROI_lat[i]), cos(p_ROI_lat[i])*sin(p_arc_lat[j]) - sin(p_ROI_lat[i])*cos(p_arc_lat[j])*cos(dlon)) + 2.*PI, 2.*PI);

// 				// update min and max vals
// 				if(d_array[idx] > d_max)
// 					d_max = d_array[idx];
// 				if(d_array[idx] < d_min)
// 					d_min = d_array[idx];
// 				if(ang_array[idx] > ang_max)
// 					ang_max = ang_array[idx];
// 				if(ang_array[idx] < ang_min)
// 					ang_min = ang_array[idx];

// 				cout << d_array[idx] << ' ' << ang_array[idx]*180./PI << ' ' << POI_ang_array[idx]*180./PI << endl;
// 			}
// 		}


// 		// compute 2D-hist
// 		// first, set up arrays
// 		azmdist_map->D_bin_edge[p->latlon_idx_proc].resize(azmdist_map->N_D_perRegion + 1);
// 		azmdist_map->D_bin_center[p->latlon_idx_proc].resize(azmdist_map->N_D_perRegion);

// 		azmdist_map->azm_bin_edge[p->latlon_idx_proc].resize(azmdist_map->N_azm_perRegion + 1);
// 		azmdist_map->azm_bin_center[p->latlon_idx_proc].resize(azmdist_map->N_azm_perRegion);

// 		azmdist_map->bin_weight[p->latlon_idx_proc].resize((azmdist_map->N_D_perRegion)*(azmdist_map->N_azm_perRegion));

// 		// next, define bins
// 		d_D   = (d_max - d_min) / double(azmdist_map->N_D_perRegion);
// 		d_ang = (ang_max - ang_min) / double(azmdist_map->N_azm_perRegion);

// 		azmdist_map->D_bin_edge[p->latlon_idx_proc][0]   = d_min;
// 		azmdist_map->azm_bin_edge[p->latlon_idx_proc][0] = ang_min;

// 		//cout << "d_min, d_max = " << d_min << ' ' << d_max << " | ang_min, ang_max = " << ang_min*180./PI << ' ' << ang_max*180./PI << endl;

// 		for (i = 0; i < azmdist_map->N_D_perRegion; ++i)
// 		{
// 			// bin edges
// 			azmdist_map->D_bin_edge[p->latlon_idx_proc][i+1]   = azmdist_map->D_bin_edge[p->latlon_idx_proc][i]   + d_D;
			
// 			// bin centers
// 			azmdist_map->D_bin_center[p->latlon_idx_proc][i]   = (azmdist_map->D_bin_edge[p->latlon_idx_proc][i] + azmdist_map->D_bin_edge[p->latlon_idx_proc][i+1]) / 2.;
// 		}
// 		for (i = 0; i < azmdist_map->N_azm_perRegion; ++i)
// 		{
// 			// bin edges
// 			azmdist_map->azm_bin_edge[p->latlon_idx_proc][i+1] = azmdist_map->azm_bin_edge[p->latlon_idx_proc][i] + d_ang;

// 			// bin centers
// 			azmdist_map->azm_bin_center[p->latlon_idx_proc][i] = (azmdist_map->azm_bin_edge[p->latlon_idx_proc][i] + azmdist_map->azm_bin_edge[p->latlon_idx_proc][i+1]) / 2.;
// 		}

// 		// for (i = 0; i < N_sample; ++i)
// 		// {
// 		// 	/* code */
// 		// }

// 	}

// 	return azmdist_map;
// }



