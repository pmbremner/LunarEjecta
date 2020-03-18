// Declaration of meteoroid flux information

// Have the option to read in cube_avg, flux_avg, and igloo_avg

#ifndef LUNAREJECTA_METEOROIDFLUX_H
#define LUNAREJECTA_METEOROIDFLUX_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;


class MEM_data
{
public:
	MEM_data(string fn);
	~MEM_data();

	void readFile();
	
private:
	void H_getRowCol_FromFile();
	void H_readInt_FromFile(ifstream& file, int& firstInt);

	string fileName;
	int rows;
	int cols;
};


class MEM_CubeAvg : public MEM_data
{
public:
	MEM_CubeAvg(string dn);
	~MEM_CubeAvg();


};




#endif 