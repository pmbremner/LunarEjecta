#include "LunarEjecta_azmdistmap.h"
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

void get_lat_min_max(double cur_lat, double dlat, double& lat_min, double& lat_max)
{
	if(fabs(cur_lat - PI/2.) < 1E-4) // +90 lat
	{
		lat_max = PI/2.;
		lat_min = lat_max - dlat/2.;
	}
	else if(fabs(cur_lat + PI/2.) < 1E-4) // -90 lat
	{
		lat_min = -PI/2.;
		lat_max = lat_min + dlat/2.;
	}
	else
	{
		lat_min = cur_lat - dlat/2.;
		lat_max = cur_lat + dlat/2.;
	}
}

hist3DSet* init_azm_dist_map(input* p)
{
	srand(time(0));
	hist3DSet* azmdist_map = new hist3DSet;

	azmdist_map->N_D_perRegion    = p->N_D_perRegion;
	azmdist_map->N_bearing_POI    = p->N_bearing_POI;
	azmdist_map->N_bearing_ROI    = p->N_bearing_ROI;



}

hist2DSet* init_azm_dist_map(input* p)
{
	srand(time(0));
	hist2DSet* azmdist_map = new hist2DSet;

	int i, j, idx, N_sample;
	double z, D, phi, D0, D1, a;
	azmdist_map->N_D_perRegion    = p->N_D_perRegion;
	azmdist_map->N_azm_perRegion  = p->N_azm_perRegion;

	azmdist_map->D_bin_edge.resize(p->N_loc);
	azmdist_map->D_bin_center.resize(p->N_loc);
	azmdist_map->azm_bin_edge.resize(p->N_loc);
	azmdist_map->azm_bin_center.resize(p->N_loc);
	azmdist_map->bin_weight.resize(p->N_loc);

	double cur_lat, cur_lon, dlat, dlon, lat_min, lat_max, lon_min, lon_max;
	double d_min, d_max, ang_min, ang_max, d_D, d_ang;

	dlat = PI / double(p->Nlat - 1.); // rad
	dlon = 2.*PI / double(p->Nlon);   // rad

	vector<double> p_ROI_lat;
	vector<double> p_ROI_lon;
	vector<double> p_arc_lat;
	vector<double> p_arc_lon;

	p_ROI_lat.resize(p->ROI_sample_points);
	p_ROI_lon.resize(p->ROI_sample_points);
	p_arc_lat.resize(p->arc_sample_points);
	p_arc_lon.resize(p->arc_sample_points);

	// loop through each location
	for (p->latlon_idx_proc = 0; p->latlon_idx_proc < p->N_loc; p->latlon_idx_proc++)
	{
		cur_lat = -PI/2. + dlat * (p->latlon_idx_min + p->latlon_idx_proc % p->Nlat);
		cur_lon = dlon * (int(p->latlon_idx_min + p->latlon_idx_proc) / int(p->Nlat));

		get_lat_min_max(cur_lat, dlat, lat_min, lat_max);

		D0 = PI/2. - lat_max; // units of lunar radii
		D1 = PI/2. - lat_min; // units of lunar radii

		cout << "azm_dist map: lat/lon = " << cur_lat*180./PI << " | " << cur_lon*180./PI << endl;
	
		// generate lat & lon points in ROI
		for (i = 0; i < p->ROI_sample_points; ++i)
		{
			z   = rand_uniform(0., 2.*sqr(sin(p->ROI_radius/(2.*p->lunar_radius))) );
			D   = 2.*asin(sqrt(z/2.)); // lunar radii
			phi = rand_uniform(0., 2.*PI);

			p_ROI_lat[i] = get_lat_from_ROI(D, phi, p->ROI_lat, p->ROI_lon);
			p_ROI_lon[i] = get_lon_from_ROI(D, phi, p->ROI_lat, p->ROI_lon);

			//cout << i << " " << p_ROI_lon[i]*180./PI << " " << p_ROI_lat[i]*180./PI << endl;
		}


		// compute lat and lon points in arc
		for (i = 0; i < p->arc_sample_points; ++i)
		{
			z   = rand_uniform(2.*sqr(sin(D0/2.)), 2.*sqr(sin(D1/2.))); // units of radius
			D   = 2.*asin(sqrt(z/2.)); // units of radius
			phi = rand_uniform(-dlon/2. + cur_lon, dlon/2. + cur_lon); // bearing

			p_arc_lat[i] = get_lat_from_pole(D, phi);
			p_arc_lon[i] = get_lon_from_pole(D, phi);

			//cout << i << " " << p_arc_lon[i]*180./PI << " " << p_arc_lat[i]*180./PI << endl;
		}

		// compute all mutual distances
		vector<double> d_array;
		vector<double> ang_array, POI_ang_array;

		N_sample = p->ROI_sample_points * p->arc_sample_points;
		d_array.resize(N_sample);
		ang_array.resize(N_sample);
		POI_ang_array.resize(N_sample);

		d_min = 2.*PI;
		d_max = 0.0;
		ang_min = 2.*PI;
		ang_max = 0.0;

		for (i = 0; i < p->ROI_sample_points; ++i)
		{
			for (j = 0; j < p->arc_sample_points; ++j)
			{
				idx = j + i*p->arc_sample_points;

				dlat = p_ROI_lat[i] - p_arc_lat[j]; // redefining here
				dlon = p_ROI_lon[i] - p_arc_lon[j];

				a = sqr(sin(dlat/2.)) + cos(p_ROI_lat[i])*cos(p_arc_lat[j])*sqr(sin(dlon/2.));

				d_array[idx] = 2.*asin(sqrt(a));

				// local azm from ROI to POI
				ang_array[idx] = initial_bearing(p_ROI_lat[i], p_arc_lat[j], p_ROI_lon[i], p_arc_lon[j]);
				//ang_array[idx] = fmod(atan2(sin(dlon)*cos(p_arc_lat[j]), cos(p_ROI_lat[i])*sin(p_arc_lat[j]) - sin(p_ROI_lat[i])*cos(p_arc_lat[j])*cos(dlon)) + 2.*PI, 2.*PI);
		
				POI_ang_array[idx] = initial_bearing(p_arc_lat[j], p_ROI_lat[i], p_arc_lon[j], p_ROI_lon[i]);
				//POI_ang_array[idx] = fmod(atan2(-sin(dlon)*cos(p_ROI_lat[i]), cos(p_ROI_lat[i])*sin(p_arc_lat[j]) - sin(p_ROI_lat[i])*cos(p_arc_lat[j])*cos(dlon)) + 2.*PI, 2.*PI);

				// update min and max vals
				if(d_array[idx] > d_max)
					d_max = d_array[idx];
				if(d_array[idx] < d_min)
					d_min = d_array[idx];
				if(ang_array[idx] > ang_max)
					ang_max = ang_array[idx];
				if(ang_array[idx] < ang_min)
					ang_min = ang_array[idx];

				cout << d_array[idx] << ' ' << ang_array[idx]*180./PI << ' ' << POI_ang_array[idx]*180./PI << endl;
			}
		}


		// compute 2D-hist
		// first, set up arrays
		azmdist_map->D_bin_edge[p->latlon_idx_proc].resize(azmdist_map->N_D_perRegion + 1);
		azmdist_map->D_bin_center[p->latlon_idx_proc].resize(azmdist_map->N_D_perRegion);

		azmdist_map->azm_bin_edge[p->latlon_idx_proc].resize(azmdist_map->N_azm_perRegion + 1);
		azmdist_map->azm_bin_center[p->latlon_idx_proc].resize(azmdist_map->N_azm_perRegion);

		azmdist_map->bin_weight[p->latlon_idx_proc].resize((azmdist_map->N_D_perRegion)*(azmdist_map->N_azm_perRegion));

		// next, define bins
		d_D   = (d_max - d_min) / double(azmdist_map->N_D_perRegion);
		d_ang = (ang_max - ang_min) / double(azmdist_map->N_azm_perRegion);

		azmdist_map->D_bin_edge[p->latlon_idx_proc][0]   = d_min;
		azmdist_map->azm_bin_edge[p->latlon_idx_proc][0] = ang_min;

		//cout << "d_min, d_max = " << d_min << ' ' << d_max << " | ang_min, ang_max = " << ang_min*180./PI << ' ' << ang_max*180./PI << endl;

		for (i = 0; i < azmdist_map->N_D_perRegion; ++i)
		{
			// bin edges
			azmdist_map->D_bin_edge[p->latlon_idx_proc][i+1]   = azmdist_map->D_bin_edge[p->latlon_idx_proc][i]   + d_D;
			
			// bin centers
			azmdist_map->D_bin_center[p->latlon_idx_proc][i]   = (azmdist_map->D_bin_edge[p->latlon_idx_proc][i] + azmdist_map->D_bin_edge[p->latlon_idx_proc][i+1]) / 2.;
		}
		for (i = 0; i < azmdist_map->N_azm_perRegion; ++i)
		{
			// bin edges
			azmdist_map->azm_bin_edge[p->latlon_idx_proc][i+1] = azmdist_map->azm_bin_edge[p->latlon_idx_proc][i] + d_ang;

			// bin centers
			azmdist_map->azm_bin_center[p->latlon_idx_proc][i] = (azmdist_map->azm_bin_edge[p->latlon_idx_proc][i] + azmdist_map->azm_bin_edge[p->latlon_idx_proc][i+1]) / 2.;
		}

		// for (i = 0; i < N_sample; ++i)
		// {
		// 	/* code */
		// }

	}

	return azmdist_map;
}