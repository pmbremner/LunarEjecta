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


void MEM_data::H_getRowCol_FromFile() {
	ifstream file;
	file.open(fileName);

	cout << fileName << endl;

	// read # of rows
	H_readInt_FromFile(file, rows);

	cout << rows << endl;

	// read # of cols
	H_readInt_FromFile(file, cols);

	cout << cols << endl;


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

MEM_CubeAvg::MEM_CubeAvg(string dn)  : MEM_data(dn + "/cube_avg.txt")
{}

MEM_CubeAvg::~MEM_CubeAvg() {}