#ifndef LUNAREJECTA_IGLOO_H
#define LUNAREJECTA_IGLOO_H

#include <string>

using namespace std;

enum primaryFluxType
{
	HiDensMEM, LoDensMEM, NEO
};

struct iglooSet // for each location
{
	string filename; // filename of igloo file
	double lat; // radians
	double lon; // radians [0, 2*pi)
	double SA;  // surface area of location on Moon

	int N_rows; // number of rows in igloo data
	int N_cols; // number of columns in igloo data
	vector<double> iglooData;   // includes info columns and data
	vector<double> speedEdge;   // size = N_cols - 8
	vector<double> speedCenter; // size = N_cols - 9
};

void H_readInt_FromFile(ifstream& file, int& firstInt);
void H_getRowCol_FromFile(string fileName, iglooSet& ig_data);
void H_read_igloo(string igloo_fn, iglooSet& ig_data);
iglooSet* read_igloo(input* p, int fluxType);
iglooSet* generate_NEO_igloo(input* p, iglooSet* hiDens);
void save_igloo(input* p, iglooSet* fluxes, int fluxType);


#endif 