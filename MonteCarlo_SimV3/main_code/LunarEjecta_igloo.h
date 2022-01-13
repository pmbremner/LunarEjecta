#ifndef LUNAREJECTA_IGLOO_H
#define LUNAREJECTA_IGLOO_H

#include <string>

using namespace std;


struct iglooSet // for each location
{
	string filename; // filename of igloo file
	double lat; // radians
	double lon; // radians [0, 2*pi), center
	double dlon; // radians, lonmin = lon - dlon/2, lonmax = lon + dlon/2
	double latmin; // radians
	double latmax; // radians
	double SA;  // surface area of location on Moon

	int N_rows; // number of rows in igloo data
	int N_cols; // number of columns in igloo data (not including the info columns, 9 large)
	vector<double> iglooData;   // includes info columns (extra 9 columns) and data
	vector<double> speedEdge;   // size = N_cols + 1
	vector<double> speedCenter; // size = N_cols
};

double sumSA(input* p, iglooSet* fluxes);
void readVelDist(string vel_fn, vector<double>& velDist);
double read_cube_zenith_flux(string cube_fn);
void H_readInt_FromFile(ifstream& file, int& firstInt);
void H_getRowCol_FromFile(string fileName, iglooSet& ig_data);
void H_read_igloo(string igloo_fn, iglooSet& ig_data);
iglooSet* read_igloo(input* p, int fluxType);
iglooSet* generate_NEO_igloo(input* p, iglooSet* hiDens);
void save_igloo(input* p, iglooSet* fluxes, int fluxType);


#endif 