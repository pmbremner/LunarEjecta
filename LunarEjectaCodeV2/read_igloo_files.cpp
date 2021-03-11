#ifndef READ_IGLOO_FILES_H
#define READ_IGLOO_FILES_H

#include "read_igloo_files.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

inline int Nrows_igloo;
inline int Ncols_igloo;

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

void H_getRowCol_FromFile(string fileName) {
	ifstream file;
	file.open(fileName);

	cout << fileName << endl;

	// read # of rows
	H_readInt_FromFile(file, Nrows_igloo);

	cout << " rows = " << Nrows_igloo << endl;

	// read # of cols
	H_readInt_FromFile(file, Ncols_igloo);

	cout << " cols = " << Ncols_igloo << endl;

	file.close();
}

void readIgloo(string lonDirectory, string dens, int lat, vector<double>& iglooData)
{
	int i, j, idx = 0;
	char C_unitName[16];
	stringstream SS_double;
	double D_temp;
	string S_temp;
	int NrowVars = 9;
	int NcolVars = 1;

	ifstream file;
	string fileName = lonDirectory + "/lat" + to_string(lat) + dens + "/igloo_avg.txt";
	
	H_getRowCol_FromFile(fileName);

	file.open(fileName);

	// skip header
	for (i = 0; i < 8; ++i)
		file.ignore(4096, '\n');

	
	// start reading data...
	iglooData.resize(Nrows_igloo * (Ncols_igloo + NrowVars));
	//rowVars.resize(Nrows * NrowVars);
	for (j = 0; j < Nrows_igloo; ++j)
	{
		//cerr << Nrows_igloo << ' ' << j << endl;		

		for (i = 0; i < Ncols_igloo + NrowVars; ++i)
		{
			file >> D_temp;
			iglooData[idx++] = D_temp;
			//cout << D_temp << ' ';
		}
		//cout << endl;

	}
	//cout << iglooData.size() << endl;
	file.close();

}

void readVelDist(string vel_fn, vector<double>& velDist)
{
	// read vel weight file
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

// mass in units of grams
void setupNEOIgloo(string lonDirectory, int lat, vector<double>& velDist, vector<double>& iglooDataHi, vector<double>& iglooDataNEO, double m_min, double m_max)
{
	
	int i, j, idx = 0;
	int NrowVars = 9;
	int NcolVars = 1;

	// get the z+ (zenith) flux from the cube file
	string cube_fn = lonDirectory + "/lat" + to_string(lat) + "/HiDensity/cube_avg.txt";
	ifstream cube_file;
	double cube_flux;
	cube_file.open(cube_fn);
	for (i = 0; i < 7; ++i)
		cube_file.ignore(4096, '\n');
	cube_file.ignore(76);
	cube_file >> cube_flux;
	cube_file.close();

	cout << cube_fn << " | flux = " << cube_flux << " #/m^2/yr\n";

	// compute number flux in mass range
	double NEO_number_flux = NEO_integral_flux(m_min) - NEO_integral_flux(m_max);
	cout << "NEO_number_flux = " << NEO_number_flux << " #/m^2/yr > " << m_min << " g and < " << m_max << " g\n";


	//copy over the info columns
	iglooDataNEO.resize(Nrows_igloo * (Ncols_igloo + NrowVars));
	cout << "total size = " << Nrows_igloo * (Ncols_igloo + NrowVars) << endl;
	for (j = 0; j < Nrows_igloo; ++j)
	{
		for (i = 0; i < NrowVars; ++i)
		{
			//cout << iglooDataHi[idx] << endl;
			iglooDataNEO[idx] = iglooDataHi[idx];
			idx++;
		}	
		idx += Ncols_igloo;
	}


	// renorm the speed columns and define the fluxes
	// divide by the cube flux, which is the flux to the surface with the zenith as norm, then multiple by the NEO flux and the speed weight
	// the NEO flux equation is also in terms of flux to a surface (same norm), we gain a factor of ~2 vs counting the igloo flux (which is needed)
	// we don't need the factor of 4 according to Althea
	for (i = 0; i < Ncols_igloo; ++i)
	{
		for (j = 0; j < Nrows_igloo; ++j)
		{
			idx = (NrowVars + i) + j * (NrowVars + Ncols_igloo);
			if(i/2 < velDist.size()-1) {
				//cout << idx << ' ' << velDist[i/2] << endl;
				iglooDataNEO[idx] = velDist[i/2] * NEO_number_flux / cube_flux * iglooDataHi[idx];
				//cout << iglooDataNEO[idx] << ' ';
			}
			else
				iglooDataNEO[idx] = 0.0;
		}
		//cout << endl;
	}
}


#endif 