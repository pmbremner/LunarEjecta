#include "LunarEjecta_params.h"
#include "LunarEjecta_igloo.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;


// input mass in units of grams, equation uses mass in kg
double NEO_integral_flux(double m)
{
	return 2.89E-11 * pow(m/1000., -0.9); // units of #/m^2/yr
}

void H_readInt_FromFile(ifstream& file, int& firstInt) {
	char C_int[8];
	file.ignore(8, ' ');
	file.get(C_int, 8, ' ');
	string S_int(C_int);
	stringstream SS_int(S_int);
	SS_int >> firstInt;
	file.ignore(32, '\n');
}

void H_getRowCol_FromFile(string fileName, iglooSet& ig_data) {
	ifstream file;
	file.open(fileName);

	cout << fileName << endl;

	// read # of rows
	H_readInt_FromFile(file, ig_data.N_rows);

	cout << " rows = " << ig_data.N_rows << endl;

	// read # of cols
	H_readInt_FromFile(file, ig_data.N_cols);

	cout << " cols = " << ig_data.N_cols << endl;

	file.close();
}

void H_read_igloo(string igloo_fn, iglooSet& ig_data)
{
	int i, j, idx = 0;
	char C_unitName[16];
	stringstream SS_double;
	double D_temp;
	string S_temp;
	int NrowVars = 9;
	//int NcolVars = 1;

	ifstream file;
	//string fileName = lonDirectory + "/lat" + to_string(lat) + dens + "/igloo_avg.txt";
	
	H_getRowCol_FromFile(igloo_fn, ig_data);

	file.open(igloo_fn);

	// skip header
	for (i = 0; i < 7; ++i)
		file.ignore(4096, '\n');

	// read speed bin centers
	file.ignore(57, '\n');

	ig_data.speedCenter.resize(ig_data.N_cols);
	ig_data.speedEdge.resize(ig_data.N_cols + 1);

	ig_data.speedEdge[0] = 0.0;

	for (i = 0; i < ig_data.N_cols; ++i)
	{
		file >> D_temp;
		ig_data.speedCenter[i] = D_temp;
		ig_data.speedEdge[i+1] = 2. * ig_data.speedCenter[i] - ig_data.speedEdge[i];

		//cout << "speed center|edge = " << ig_data.speedCenter[i] << " | " << ig_data.speedEdge[i+1] << " km/s\n";
	}

	
	// start reading data...
	ig_data.iglooData.resize(ig_data.N_rows * (ig_data.N_cols + NrowVars));
	//rowVars.resize(Nrows * NrowVars);
	for (j = 0; j < ig_data.N_rows; ++j)
	{
		//cerr << ig_data.N_rows << ' ' << j << endl;		

		for (i = 0; i < ig_data.N_cols + NrowVars; ++i)
		{
			file >> D_temp;
			ig_data.iglooData[idx++] = D_temp;
			//cout << D_temp << ' ';
		}
		//cout << endl;

	}
	//cout << ig_data.iglooData.size() << endl;
	file.close();

}

iglooSet* read_igloo(input* p, int fluxType)
{
	string cur_filename;
	int idx_lon, cur_lat;
	iglooSet* ig = new iglooSet[p->N_loc];

	for (p->latlon_idx_proc = 0; p->latlon_idx_proc < p->N_loc; p->latlon_idx_proc++)
	{
		p->latlon_idx_cur = p->latlon_idx_min + p->latlon_idx_proc; // absolute index
	
		idx_lon = (p->latlon_idx_cur) / (p->Nlat); // integer division

		//cout << ((p->latlon_idx_cur) % (p->Nlat)) << endl;

		cur_lat = -90 + 180/(p->Nlat - 1)*((p->latlon_idx_cur) % (p->Nlat)); // modulo

		ig[p->latlon_idx_proc].lat = cur_lat * PI / 180.;
		ig[p->latlon_idx_proc].lon = 2. * PI * idx_lon / p->Nlon; 

		cout << "lat-lon = " << ig[p->latlon_idx_proc].lat * 180. / PI << ' ' << ig[p->latlon_idx_proc].lon * 180. / PI << endl;

		if (fluxType == HiDensMEM)
		{
			cur_filename = p->lon_directory[idx_lon] + "/lat" + to_string(cur_lat) + "/HiDensity/igloo_avg.txt";
		}
		else if (fluxType == LoDensMEM)
		{
			cur_filename = p->lon_directory[idx_lon] + "/lat" + to_string(cur_lat) + "/LoDensity/igloo_avg.txt";
		}
		else if (fluxType == NEO)
		{
			cur_filename = p->lon_directory[idx_lon] + "/lat" + to_string(cur_lat) + "/NEO_igloo_avg.txt";
		}

		ig[p->latlon_idx_proc].filename = cur_filename;

		H_read_igloo(cur_filename, ig[p->latlon_idx_proc]);
	}

	return ig;

}



iglooSet* generate_NEO_igloo(input* p, iglooSet* hiDens)
{
	iglooSet* ig = new iglooSet[p->N_loc];
	int NrowVars = 9, i_row, j_col, idx_igloo;


	// compute number flux in mass range
	double NEO_number_flux = NEO_integral_flux(p->NEO_massMin) - NEO_integral_flux(p->NEO_massMax);
	cout << "NEO_number_flux = " << NEO_number_flux << " #/m^2/yr > " << p->NEO_massMin << " g and < " << p->NEO_massMax << " g\n";

	for (p->latlon_idx_proc = 0; p->latlon_idx_proc < p->N_loc; p->latlon_idx_proc++)
	{
		cout << "generate NEO fluxes: " << p->latlon_idx_proc << "/" << p->N_loc << endl;

		ig[p->latlon_idx_proc].lat     = hiDens[p->latlon_idx_proc].lat;
		ig[p->latlon_idx_proc].lon     = hiDens[p->latlon_idx_proc].lon;
		ig[p->latlon_idx_proc].SA      = hiDens[p->latlon_idx_proc].SA;
		ig[p->latlon_idx_proc].N_rows  = hiDens[p->latlon_idx_proc].N_rows;
		ig[p->latlon_idx_proc].N_cols  = hiDens[p->latlon_idx_proc].N_cols;

		ig[p->latlon_idx_proc].speedCenter.resize(ig[p->latlon_idx_proc].N_cols);
		ig[p->latlon_idx_proc].speedEdge.resize(ig[p->latlon_idx_proc].N_cols + 1);

		ig[p->latlon_idx_proc].speedEdge[0] = 0.;
		for (int i = 0; i < ig[p->latlon_idx_proc].N_cols; ++i)
		{
			ig[p->latlon_idx_proc].speedCenter[i] = hiDens[p->latlon_idx_proc].speedCenter[i];
			ig[p->latlon_idx_proc].speedEdge[i+1] = hiDens[p->latlon_idx_proc].speedEdge[i+1];
		}

		ig[p->latlon_idx_proc].iglooData.resize(ig[p->latlon_idx_proc].N_rows * (ig[p->latlon_idx_proc].N_cols + NrowVars));

		for (i_row = 0; i_row < ig[p->latlon_idx_proc].N_rows; ++i_row)
		{
			for (j_col = 0; j_col < ig[p->latlon_idx_proc].N_cols + NrowVars; ++j_col)
			{
				idx_igloo = j_col + i_row * (ig[p->latlon_idx_proc].N_cols + NrowVars);

				if (j_col < NrowVars)
					ig[p->latlon_idx_proc].iglooData[idx_igloo] = hiDens[p->latlon_idx_proc].iglooData[idx_igloo];
				else
					ig[p->latlon_idx_proc].iglooData[idx_igloo] = NEO_number_flux * hiDens[p->latlon_idx_proc].iglooData[idx_igloo];
				
			}
		}
	}

	


	return ig;
}


void save_igloo(input* p, iglooSet* fluxes, int fluxType)
{

}