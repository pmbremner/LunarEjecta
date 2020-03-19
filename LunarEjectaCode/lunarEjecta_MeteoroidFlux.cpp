#include "lunarEjecta_MeteoroidFlux.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

MEM_data::MEM_data(string fn) {
	fileName = fn;

	this->H_getRowCol_FromFile();
}

MEM_data::~MEM_data() {}

inline double MEM_data::getFlux(int r, int c)
{
	return fluxData[H_idxFlux(r,c)];
}

inline double MEM_data::getRVar(int r, int c)
{
	return rowVars[H_idxRVar(r,c)];
}

inline double MEM_data::getCVar(int r, int c)
{
	return colVars[H_idxCVar(r,c)];
}

void MEM_data::H_setHeaderLength(int h) { headerLength = h; }


void MEM_data::H_getRowCol_FromFile() {
	ifstream file;
	file.open(fileName);

	cout << fileName << endl;

	// read # of rows
	H_readInt_FromFile(file, Nrows);

	cout << Nrows << endl;

	// read # of cols
	H_readInt_FromFile(file, Ncols);

	cout << Ncols << endl;


}


void MEM_data::H_readInt_FromFile(ifstream& file, int& firstInt) {
	char C_int[8];
	file.ignore(8, ' ');
	file.get(C_int, 8, ' ');
	string S_int(C_int);
	stringstream SS_int(S_int);
	SS_int >> firstInt;
	file.ignore(32, '\n');
}

inline int MEM_data::H_idxFlux(int r, int c)
{
	return r * Ncols + c;
}

inline int MEM_data::H_idxRVar(int r, int c)
{
	return r * NrowVars + c;
}

inline int MEM_data::H_idxCVar(int r, int c)
{
	return r + c * NcolVars;
}

// assumes fluxData has already been initialized
inline void MEM_data::H_storeFlux(int r, int c, double flux)
{
	fluxData[H_idxFlux(r,c)] = flux;
}

inline void MEM_data::H_pushBackFlux(double flux)
{
	fluxData.push_back(flux);
}

inline void MEM_data::H_pushBackRVar(double RVar)
{
	rowVars.push_back(RVar);
}

inline void MEM_data::H_pushBackCVar(double CVar)
{
	colVars.push_back(CVar);
}

/////////////////////////////////////////////////////////////

MEM_cubeAvg::MEM_cubeAvg(string dn)  : MEM_data(dn + "/cube_avg.txt")
{
	this->H_readFile();
}

MEM_cubeAvg::~MEM_cubeAvg() {}

void MEM_cubeAvg::H_readFile(void)
{
	int i, j;
	char C_unitName[16];
	stringstream SS_double;
	double D_temp;
	NrowVars = 1;
	NcolVars = 12;

	ifstream file;
	file.open(fileName);

	// retrieve the column names
	for (i = 0; i < 5; ++i)
		file.ignore(256, '\n');
	file.ignore(23);
	for (i = 0; i < NcolVars; ++i)
	{
		file.get(C_unitName, 14, '\n');
		colVarsUnits.push_back(string(C_unitName));
		cout << colVarsUnits[i] << endl;
	}

	// retrieve the row name

	file.ignore(256, '\n');
	file.ignore(256, '\n');
	file.ignore(256, '\n');
	file.ignore(2, ' ');
	file.get(C_unitName, 16);
	rowVarsUnits.push_back(C_unitName);
	cout << rowVarsUnits[0] << endl;
	file.ignore(256, '\n');

	// start reading data...
	fluxData.resize(Nrows * Ncols);
	for (j = 0; j < Nrows; ++j)
	{
		file >> D_temp;
		this->H_pushBackRVar(D_temp);
		//cout << this->getRVar(j, 0) << ' ';
		for (i = 0; i < Ncols; ++i)
		{
			file >> D_temp;
			this->H_pushBackFlux(D_temp);
			// if(i < 2)
			// 	cout << this->getFlux(j, i) << ' ';
		}
		//cout << endl;
	}

}

/////////////////////////////////////////////////////////////

MEM_fluxAvg::MEM_fluxAvg(string dn)  : MEM_data(dn + "/flux_avg.txt")
{
	this->H_readFile();
}

MEM_fluxAvg::~MEM_fluxAvg() {}

void MEM_fluxAvg::H_readFile(void)
{

}

/////////////////////////////////////////////////////////////

MEM_iglooAvg::MEM_iglooAvg(string dn)  : MEM_data(dn + "/igloo_avg.txt")
{
	this->H_readFile();
}

MEM_iglooAvg::~MEM_iglooAvg() {}

void MEM_iglooAvg::H_readFile(void)
{

}