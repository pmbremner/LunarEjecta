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

	inline double getFlux(int row, int col);
	inline double getRVar(int row, int col);
	inline double getCVar(int row, int col);
	
protected:
	virtual void H_readFile(void) = 0;
	void H_setHeaderLength(int h);
	void H_getRowCol_FromFile();
	void H_readInt_FromFile(ifstream& file, int& firstInt);
	inline int H_idxFlux(int row, int col);
	inline int H_idxRVar(int row, int col);
	inline int H_idxCVar(int row, int col);
	inline void H_storeFlux(int row, int col, double flux); // use only after init
	inline void H_pushBackFlux(double flux); // use only for init
	inline void H_pushBackRVar(double RVar); // use only for init
	inline void H_pushBackCVar(double CVar); // use only for init

	string fileName;
	int headerLength;

	int Nrows; // of data only
	int Ncols; // of data only
	int NrowVars;
	int NcolVars;

	vector<double> fluxData; // total size = Nrows * Ncols
	vector<double> rowVars;  // total size = Nrows * NrowVars
	vector<double> colVars;  // total size = Ncols * NcolVars

	string fluxUnits;            
	vector<string> rowVarsUnits; 
	vector<string> colVarsUnits; 
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