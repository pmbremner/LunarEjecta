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
	double SA;  // surface area of location on Moon, units of r^2
	double netFlux; // [#/m^2/yr]

	int N_rows; // number of rows in igloo data
	int N_cols; // number of columns in igloo data (not including the info columns, 9 large)
	vector<double> iglooData;   // includes info columns (extra 9 columns) and data
	vector<double> speedEdge;   // size = N_cols + 1
	vector<double> speedCenter; // size = N_cols

	// If MEM, from dens files, if NEO, just a single value
	vector<double> dens_left;
	vector<double> dens_right;
	vector<double> dens_pdf;
	vector<double> dens_cdf; // the first value will always be zero, so the size is one larger than dens_pdf

	// mass cdf, MEM uses Grun and NEO uses Moorhead
	vector<double> mass_edge; // starts high and goes low
	vector<double> mass_cdf;

};

int ig_idx(int row, int col, int N_col);

double sumSA(input* p, iglooSet* fluxes);
void readVelDist(string vel_fn, vector<double>& velDist);
double read_cube_zenith_flux(string cube_fn);
void H_readInt_FromFile(ifstream& file, int& firstInt);
void H_getRowCol_FromFile(string fileName, iglooSet& ig_data);
void H_read_igloo(string igloo_fn, iglooSet& ig_data);
iglooSet* read_igloo(input* p, int fluxType);
iglooSet* generate_NEO_igloo(input* p, iglooSet* hiDens);
void save_igloo(input* p, iglooSet* fluxes, int fluxType);

void H_read_dens(string dens_fn, iglooSet& ig_data);

void net_flux(iglooSet& fluxes);

#endif 