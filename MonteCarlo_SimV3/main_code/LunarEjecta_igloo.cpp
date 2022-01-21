#include "LunarEjecta_params.h"
#include "LunarEjecta_igloo.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

// N_col doesn't take into account the first 9 columns, need to do that here
int ig_idx(int row, int col, int N_col)
{
	return col + (N_col + 9) * row;
}

double sumSA(input* p, iglooSet* fluxes)
{
	double sum = 0.;
	for (int i = 0; i < p->N_loc; ++i)
		sum += fluxes[i].SA;

	return sum;
}


void readVelDist(string vel_fn, vector<double>& velDist)
{
	cout << vel_fn << endl;
	//read vel weight file
	ifstream NEOvel_file;
	NEOvel_file.open(vel_fn);

	double D_temp;
	int i = 0;
	
	velDist.resize(0);

	NEOvel_file.ignore(256, '\n');
	while(!NEOvel_file.eof())
	{
		NEOvel_file >> D_temp;
		NEOvel_file >> D_temp;

		velDist.push_back(D_temp);
		//cout << D_temp << endl;
	}


	NEOvel_file.close();
}

double read_cube_zenith_flux(string cube_fn)
{
	cout << cube_fn << endl;
	ifstream cube_file;
	double cube_flux;
	cube_file.open(cube_fn);
	for (int i = 0; i < 7; ++i)
		cube_file.ignore(4096, '\n');
	cube_file.ignore(76);
	cube_file >> cube_flux;
	cube_file.close();
	cout << " cube flux = " << cube_flux << " #/m^2/yr\n";
	return cube_flux;
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
void H_read_dens(string dens_fn, iglooSet& ig_data)
{
	int i;
	const int N_rows = 158;

	ifstream file;

	file.open(dens_fn);
	cout << "reading: " << dens_fn << endl;

	// skip header
	for (i = 0; i < 2; ++i)
		file.ignore(4096, '\n');

	ig_data.dens_left.clear();
	ig_data.dens_right.clear();
	ig_data.dens_pdf.clear();
	ig_data.dens_cdf.clear();

	ig_data.dens_left.resize(N_rows, 0.);
	ig_data.dens_right.resize(N_rows, 0.);
	ig_data.dens_pdf.resize(N_rows, 0.);
	ig_data.dens_cdf.resize(N_rows+1, 0.);


	for (i = 0; i < N_rows; ++i)
	{
		file >> ig_data.dens_left[i];
		file >> ig_data.dens_right[i];
		file >> ig_data.dens_pdf[i]; 
		ig_data.dens_cdf[i+1] = ig_data.dens_cdf[i] + ig_data.dens_pdf[i];
	}

	// make sure the cdf is normalized to 1
	for (i = 0; i < N_rows+1; ++i){
		ig_data.dens_cdf[i] /= ig_data.dens_cdf[N_rows];

		// if (i < N_rows)
		// 	cout << ig_data.dens_left[i] << ' '
		// 	     << ig_data.dens_right[i] << ' '
		// 	     << ig_data.dens_cdf[i] << ' '
		// 	     << ig_data.dens_pdf[i] << endl;
	}
	// cout << ig_data.dens_left[N_rows-1] << ' '
	//      << ig_data.dens_right[N_rows-1] << ' '
	//      << ig_data.dens_cdf[N_rows] << ' '
	//      << ig_data.dens_pdf[N_rows-1] << endl;

	file.close();
}

void H_set_NEO_dens(input* p, iglooSet& ig_data)
{
	const int N_rows = 1;

	ig_data.dens_left.clear();
	ig_data.dens_right.clear();
	ig_data.dens_pdf.clear();
	ig_data.dens_cdf.clear();

	ig_data.dens_left.resize(N_rows);
	ig_data.dens_right.resize(N_rows);
	ig_data.dens_pdf.resize(N_rows);
	ig_data.dens_cdf.resize(N_rows+1, 0.);

	ig_data.dens_left[0]   = p->NEO_dens;
	ig_data.dens_right[0]  = p->NEO_dens;
	ig_data.dens_pdf[0] = 1.;
	ig_data.dens_cdf[1] = 1.;

	cout << ig_data.dens_left[0] << ' '
	     << ig_data.dens_right[0] << ' '
	     << ig_data.dens_cdf[1] << ' '
	     << ig_data.dens_pdf[0] << endl;
}



void H_gen_mass_cdf(input* p, iglooSet& ig_data, int fluxType)
{
	const int N = 1000;
	int i;

	ig_data.mass_edge.clear();
	ig_data.mass_cdf.clear();

	ig_data.mass_edge.resize(N);
	ig_data.mass_cdf.resize(N);

	if (fluxType == HiDensMEM || fluxType == LoDensMEM){

		rlogspace(ig_data.mass_edge, log10(p->MEM_massMin), log10(p->MEM_massMax), N);

		for (i = 0; i < N; ++i)
			ig_data.mass_cdf[i] = ( MEM_mass_grun(ig_data.mass_edge[i]) - MEM_mass_grun(ig_data.mass_edge[0]) ) / ( MEM_mass_grun(ig_data.mass_edge[N-1]) - MEM_mass_grun(ig_data.mass_edge[0]) );

	}
	else{ // fluxType == NEO

		rlogspace(ig_data.mass_edge, log10(p->NEO_massMin), log10(p->NEO_massMax), N);

		for (i = 0; i < N; ++i)
			ig_data.mass_cdf[i] = ( NEO_integral_flux(ig_data.mass_edge[i]) - NEO_integral_flux(ig_data.mass_edge[0]) ) / ( NEO_integral_flux(ig_data.mass_edge[N-1]) - NEO_integral_flux(ig_data.mass_edge[0]) );

	}
	// Force the end points to be correct (fixing floating point errors in above calcs)
	ig_data.mass_cdf[0] = 0.;
	ig_data.mass_cdf[N-1] = 1.;
	// for (int i = 0; i < N; ++i)
	// 	cout << ig_data.mass_edge[i] << ' ' << ig_data.mass_cdf[i] << endl;
}


iglooSet* read_igloo(input* p, int fluxType)
{
	string cur_filename, dens_filename;
	int idx_lon, cur_lat;
	//double lat_min, lat_max;
	iglooSet* ig = new iglooSet[p->N_loc];
	double dlat, dlon;

	dlat = PI/double(p->Nlat - 1.); // rad
	dlon = 2.*PI/double(p->Nlon);   // rad

	for (p->latlon_idx_proc = 0; p->latlon_idx_proc < p->N_loc; p->latlon_idx_proc++)
	{
		p->latlon_idx_cur = p->latlon_idx_min + p->latlon_idx_proc; // absolute index
	
		idx_lon = (p->latlon_idx_cur) / (p->Nlat); // integer division

		//cout << ((p->latlon_idx_cur) % (p->Nlat)) << endl;

		cur_lat = -90 + 180/(p->Nlat - 1)*((p->latlon_idx_cur) % (p->Nlat)); // modulo

		ig[p->latlon_idx_proc].lat = cur_lat * PI / 180.;
		ig[p->latlon_idx_proc].lon = 2. * PI * idx_lon / p->Nlon;

		ig[p->latlon_idx_proc].dlon = dlon;

		cout << "lat-lon = " << ig[p->latlon_idx_proc].lat * 180. / PI << ' ' << ig[p->latlon_idx_proc].lon * 180. / PI << endl;

		// compute surface area of location
		// all in rads

		get_lat_min_max(ig[p->latlon_idx_proc].lat, dlat, ig[p->latlon_idx_proc].latmin, ig[p->latlon_idx_proc].latmax);

		// check for region that overlaps with equator
		ig[p->latlon_idx_proc].SA = calcSA(ig[p->latlon_idx_proc].lat, ig[p->latlon_idx_proc].latmin, ig[p->latlon_idx_proc].latmax, dlat, ig[p->latlon_idx_proc].dlon);


		cout << "lat dlat/min/max = " << dlat*180./PI << " | " << ig[p->latlon_idx_proc].latmin*180./PI << " | " << ig[p->latlon_idx_proc].latmax*180./PI << endl;
		cout << "SA = " << ig[p->latlon_idx_proc].SA / (4.*PI) << " x 4pi r^2\n";

		// generate filename for igloo file to read
		if (fluxType == HiDensMEM)
		{
			cur_filename  = p->lon_directory[idx_lon] + "/lat" + to_string(cur_lat) + "/HiDensity/igloo_avg.txt";
			dens_filename = p->lon_directory[idx_lon] + "/lat" + to_string(cur_lat) + "/hidensity.txt";

			H_read_dens(dens_filename, ig[p->latlon_idx_proc]);
		}
		else if (fluxType == LoDensMEM)
		{
			cur_filename  = p->lon_directory[idx_lon] + "/lat" + to_string(cur_lat) + "/LoDensity/igloo_avg.txt";
			dens_filename = p->lon_directory[idx_lon] + "/lat" + to_string(cur_lat) + "/lodensity.txt";

			H_read_dens(dens_filename, ig[p->latlon_idx_proc]);
		}
		else if (fluxType == NEO)
		{
			cur_filename = p->lon_directory[idx_lon] + "/lat" + to_string(cur_lat) + "/NEO_igloo_avg.txt";

			H_set_NEO_dens(p, ig[p->latlon_idx_proc]);
		}

		H_gen_mass_cdf(p, ig[p->latlon_idx_proc], fluxType);

		ig[p->latlon_idx_proc].filename = cur_filename;

		H_read_igloo(cur_filename, ig[p->latlon_idx_proc]);

		net_flux(ig[p->latlon_idx_proc]);

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

	double cube_flux;
	string cube_fn, NEO_fn;

	// read NEO speed distribution
	vector<double> v_NEO;
	//cout << p->NEO_velDist_fn << endl;
	readVelDist(p->NEO_velDist_fn, v_NEO);

	// loop through each location
	for (p->latlon_idx_proc = 0; p->latlon_idx_proc < p->N_loc; p->latlon_idx_proc++)
	{
		cout << "generate NEO fluxes: " << p->latlon_idx_proc << "/" << p->N_loc << endl;

		cube_fn = hiDens[p->latlon_idx_proc].filename;
		cube_fn = hiDens[p->latlon_idx_proc].filename.substr(0, cube_fn.length() - 13) + "cube_avg.txt";

		NEO_fn = hiDens[p->latlon_idx_proc].filename;
		NEO_fn = hiDens[p->latlon_idx_proc].filename.substr(0, cube_fn.length() - 22) + "NEO_igloo_avg.txt";
		ig[p->latlon_idx_proc].filename = NEO_fn;

		//cout << cube_fn << endl;

		cube_flux = read_cube_zenith_flux(cube_fn);

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

		// loop through igloo data
		for (i_row = 0; i_row < ig[p->latlon_idx_proc].N_rows; ++i_row)
		{
			for (j_col = 0; j_col < ig[p->latlon_idx_proc].N_cols + NrowVars; ++j_col)
			{
				idx_igloo = j_col + i_row * (ig[p->latlon_idx_proc].N_cols + NrowVars);

				if (j_col < NrowVars) // info columns
					ig[p->latlon_idx_proc].iglooData[idx_igloo] = hiDens[p->latlon_idx_proc].iglooData[idx_igloo];
				else if ((j_col-NrowVars)/2 < v_NEO.size()-1) // assuming MEM speed bin width is 1 km/s
				{
					// multiply by speed PDF and NEO flux, dividing by the cube flux in the +z direction (NEO fluxes given in planar flux, not igloo style)
					ig[p->latlon_idx_proc].iglooData[idx_igloo] = v_NEO[(j_col-NrowVars)/2] * NEO_number_flux /cube_flux * hiDens[p->latlon_idx_proc].iglooData[idx_igloo];
				}
				else // outside of vel dist range
					ig[p->latlon_idx_proc].iglooData[idx_igloo] = 0.0;
			}
		}
	}

	return ig;
}


void save_igloo(input* p, iglooSet* fluxes, int fluxType)
{
	int i_col, j_row, NrowVars = 9, idx;
	ofstream igloo_file;

	// loop through each location
	for (p->latlon_idx_proc = 0; p->latlon_idx_proc < p->N_loc; p->latlon_idx_proc++)
	{
		igloo_file.open(fluxes[p->latlon_idx_proc].filename);
		cout << "Saving to " << fluxes[p->latlon_idx_proc].filename << endl;

		// header info
		igloo_file << "# " << fluxes[p->latlon_idx_proc].N_rows << " rows of flux data\n";
		igloo_file << "# " << fluxes[p->latlon_idx_proc].N_cols << " columns of flux data\n";
		igloo_file << "# File generated from LMEEM code\n";
		igloo_file << "# flux is distributed by angle and speed:\n";
		igloo_file << "# angle values give boundaries and center of angular bin (degrees)\n";
		igloo_file << "# numeric column labels give midpoint of speed range (km/s)\n";
		igloo_file << "#\n";
		igloo_file << "#  ID   I   J   PHI1   PHI2 THETA1 THETA2 PHIavg THETAavg ";

		// speed bin centers
		for (i_col = 0; i_col < fluxes[p->latlon_idx_proc].N_cols; ++i_col)
			igloo_file << fluxes[p->latlon_idx_proc].speedCenter[i_col] << " ";
		igloo_file << endl;

		// print the data
		for (j_row = 0; j_row < fluxes[p->latlon_idx_proc].N_rows; ++j_row)
		{
			for (i_col = 0; i_col < NrowVars + fluxes[p->latlon_idx_proc].N_cols; ++i_col)
			{
				idx = i_col + j_row * (NrowVars + fluxes[p->latlon_idx_proc].N_cols);

				igloo_file << fluxes[p->latlon_idx_proc].iglooData[idx] << ' ';
			}
			if (j_row < fluxes[p->latlon_idx_proc].N_rows - 1)
				igloo_file << endl;
		}

		igloo_file.close();
	}
	
}


// computes the net flux of the igloo set, and stores it in the igloo set params
void net_flux(iglooSet& fluxes)
{
	int i, j, idx;
	double sum = 0.;
	double partial_sum = 0.;

	for (i = 0; i < fluxes.N_rows; ++i)
	{
		partial_sum = 0;
		for (j = 0; j < fluxes.N_cols; ++j)
		{

			idx = ig_idx(i, 9 + j, fluxes.N_cols);
			//idx = 9 + j + i * (fluxes.N_cols + 9);

			partial_sum += fluxes.iglooData[idx];
		}
		sum += partial_sum;
	}

	fluxes.netFlux = sum; // note, this should correspond to the total cross-sectional flux in the cube_avg file
	cout << "Net flux = " << sum << " #/m^2/yr\n";
}