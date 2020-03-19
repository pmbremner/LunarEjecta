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

/////////////////////////////////////////////////////////////

MEM_cubeAvg::MEM_cubeAvg(string dn)  : MEM_data(dn + "/cube_avg.txt")
{
	this->H_readFile();
}

MEM_cubeAvg::~MEM_cubeAvg() {}

void MEM_cubeAvg::H_readFile(void)
{

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