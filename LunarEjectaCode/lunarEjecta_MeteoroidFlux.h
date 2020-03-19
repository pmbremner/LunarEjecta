// Declaration of meteoroid flux information

// Have the option to read in cube_avg, flux_avg, and igloo_avg

#ifndef LUNAREJECTA_METEOROIDFLUX_H
#define LUNAREJECTA_METEOROIDFLUX_H

#include <string>
#include <fstream>
#include <vector>
#include "vector2d.h"

using namespace std;

// This is an abstract base class and cannot be used directly
class MEM_data
{
public:
	MEM_data(string fn);
	~MEM_data();
	
protected:
	virtual void H_readFile(void) = 0;
	void H_setHeaderLength(int h);
	void H_getRowCol_FromFile();
	void H_readInt_FromFile(ifstream& file, int& firstInt);

	string fileName;
	int headerLength;
	int Nrows;
	int Ncols;
	int Nvars;
	vector<double> 
	vector<double> fluxData;
};

/////////////////////////////////////////////////////////////
class MEM_cubeAvg : public MEM_data
{
public:
	MEM_cubeAvg(string dn); // dn = directory name
	~MEM_cubeAvg();
private:
	void H_readFile();


};
/////////////////////////////////////////////////////////////
class MEM_fluxAvg : public MEM_data
{
public:
	MEM_fluxAvg(string dn);
	~MEM_fluxAvg();
private:
	void H_readFile();

};
/////////////////////////////////////////////////////////////
class MEM_iglooAvg : public MEM_data
{
public:
	MEM_iglooAvg(string dn);
	~MEM_iglooAvg();
private:
	void H_readFile();

};




#endif 