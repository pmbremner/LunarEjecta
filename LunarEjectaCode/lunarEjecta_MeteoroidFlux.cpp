#include "lunarEjecta_MeteoroidFlux.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

MEM_data::MEM_data(string fn) {
	fileName = fn;
	iotype = 1; // input
	vMin = 0.0;  // km/s
	vMax = 80.0; // km/s

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

inline double MEM_data::getvMin() {return vMin;}
inline double MEM_data::getvMax() {return vMax;}


inline double MEM_data::getNrows() {return Nrows;}
inline double MEM_data::getNcols() {return Ncols;}
inline double MEM_data::getNrowVars() {return NrowVars;}
inline double MEM_data::getNcolVars() {return NcolVars;}

void MEM_data::H_setHeaderLength(int h) { headerLength = h; }


void MEM_data::H_getRowCol_FromFile() {
	ifstream file;
	file.open(fileName);

	cout << fileName << endl;

	// read # of rows
	H_readInt_FromFile(file, Nrows);

	cout << " rows = " << Nrows << endl;

	// read # of cols
	H_readInt_FromFile(file, Ncols);

	cout << " cols = " << Ncols << endl;

	file.close();
}

// pulls from the options.txt file
void MEM_data::H_getVelAngleResolution(string fn) {
	ifstream file;
	file.open(fn);
	cout << fn << endl;
	char C_int[64];

	for (int i = 0; i < 39; ++i)
		file.ignore(256, '\n');

	file >> dAngle >> dVel;
	Nvel = (vMax - vMin) / dVel;
	Ntheta = 360/dAngle;
	Nphi = 180/dAngle;

	cout << " angle res: " << dAngle << " | vel res: " << dVel << endl;
	cout << " Nvel = " << Nvel << " | Ntheta = " << Ntheta << " | Nphi = " << Nphi << endl;
	file.close();
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
	return r * NcolVars + c;
}

// assumes fluxData has already been initialized
inline void MEM_data::H_storeFlux(int r, int c, double flux)
{
	fluxData[H_idxFlux(r,c)] = flux;
}

inline void MEM_data::H_storeRVar(int r, int c, double RVar)
{
	rowVars[H_idxRVar(r,c)] = RVar;
}

inline void MEM_data::H_storeCVar(int r, int c, double CVar)
{
	colVars[H_idxCVar(r,c)] = CVar;
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

// assumes alt and azm are in degrees
double MEM_cubeAvg::getFlux_atAngleVel(double alt, double azm, double vel)
{
	double x, y, z;
	int row;
	alt *= PI / 180.0;
	azm *= PI / 180.0;

	x = cos(azm)*cos(alt);
	y = sin(azm)*cos(alt);
	z = sin(alt);

	if (vel < 0.0)
	{
		cerr << "ERROR: MEM_cubeAvg::getFlux_atAngleVel: vel = " << vel << " < 0\n";
		return 0.0;
	}
	else if (vel < dVel/2.0)
	{
		return H_cubeFluxAvg_atRow(x, y, z, 0);
	}
	else if (vel < vMax - dVel/2.0) // linear interp
	{
		row = int(vel/dVel + 0.5);
		return (rowVars[row] - vel)/dVel*H_cubeFluxAvg_atRow(x, y, z, row-1)
			 + (vel - rowVars[row-1]  )/dVel*H_cubeFluxAvg_atRow(x, y, z, row);
	}
	else if (vel <= vMax)
	{
		return H_cubeFluxAvg_atRow(x, y, z, Nrows-1);
	} else if (vel > vMax)
	{
		cerr << "ERROR: MEM_cubeAvg::getFlux_atAngleVel: vel " << vel  << " > " << vMax << endl;
		return 0.0;
	} else
	{
		cerr << "ERROR: MEM_cubeAvg::getFlux_atAngleVel: vel " << vel << " invalid\n";
		return 0.0;
	}
}


double MEM_cubeAvg::H_cubeFluxAvg_atRow(double x, double y, double z, int row)
{ // ignores any nadir fluxes, assumes we are on the surface
	return (x > 0.0 ? this->getFlux(row, p_xRam)       * sqr(x) : 0.0)
		 + (x < 0.0 ? this->getFlux(row, m_xWake)      * sqr(x) : 0.0)
		 + (y > 0.0 ? this->getFlux(row, p_yPort)      * sqr(y) : 0.0)
		 + (y < 0.0 ? this->getFlux(row, m_yStarboard) * sqr(y) : 0.0)
		 + (z > 0.0 ? this->getFlux(row, p_zZenith)    * sqr(z) : 0.0);
}


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

	// set colVars
	colVars.resize(NcolVars);
	for (i = 0; i < NcolVars; ++i)
		this->H_storeCVar(0, i, double(i));


	// retrieve the column names
	for (i = 0; i < 5; ++i)
		file.ignore(256, '\n');
	file.ignore(23);
	for (i = 0; i < NcolVars; ++i)
	{
		file.get(C_unitName, 14, '\n');
		colVarsUnits.push_back(string(C_unitName));
		//cout << colVarsUnits[i] << endl;
	}

	// retrieve the row name
	file.ignore(256, '\n');
	file.ignore(256, '\n');
	file.ignore(256, '\n');
	file.ignore(2, ' ');
	file.get(C_unitName, 16);
	rowVarsUnits.push_back(C_unitName);
	//cout << rowVarsUnits[0] << endl;
	file.ignore(256, '\n');

	// start reading data...
	fluxData.resize(Nrows * Ncols);
	rowVars.resize(Nrows * NrowVars);
	for (j = 0; j < Nrows; ++j)
	{
		file >> D_temp;
		this->H_storeRVar(j, 0, D_temp);
		//cout << this->getRVar(j, 0) << ' ';
		for (i = 0; i < Ncols; ++i)
		{
			file >> D_temp;
			this->H_storeFlux(j, i, D_temp);
			// if(i < 2)
			// 	cout << this->getFlux(j, i) << ' ';
		}
		//cout << endl;
	}
	cout << fluxData.size() << endl;
	file.close();
}

/////////////////////////////////////////////////////////////

MEM_fluxAvg::MEM_fluxAvg(string dn)  : MEM_data(dn + "/flux_avg.txt")
{
	this->H_readFile();
}

MEM_fluxAvg::~MEM_fluxAvg() {}

// no interpolation at this time
double MEM_fluxAvg::getFlux_atAngleVel(double alt, double azm, double vel)
{
	int row_idx = Ntheta * int((alt + 90) / dAngle) + int(azm / dAngle);
	int col_idx = vel / dVel;

	if(row_idx >= 0 && row_idx < Nrows && col_idx >= 0 && col_idx < Ncols)
	{
		return this->getFlux(row_idx, col_idx);
	}
	else // out of bounds
	{
		cerr << "MEM_fluxAvg::getFlux_atAngleVel: index out of bounds " << row_idx << ' ' << col_idx << " invalid\n";
		return 0.0;
	}
}

void MEM_fluxAvg::H_readFile(void)
{
	int i, j;
	char C_unitName[16];
	stringstream SS_double;
	double D_temp;
	NrowVars = 2;
	NcolVars = 1;

	ifstream file;
	file.open(fileName);

	// retrieve the column units
	colVarsUnits.push_back(string("speed (km/s)"));
	//cout << colVarsUnits[0] << endl;

	for (i = 0; i < 7; ++i)
		file.ignore(256, '\n');

	// retrieve the row units
	file.ignore(1);
	file.get(C_unitName, 6);
	rowVarsUnits.push_back(string(C_unitName) + " (degrees)");
	//cout << rowVarsUnits[0] << endl;

	file.get(C_unitName, 8);
	rowVarsUnits.push_back(string(C_unitName) + " (degrees)");
	//cout << rowVarsUnits[1] << endl;

	// retrieve colVars
	colVars.resize(Ncols * NcolVars);
	for (i = 0; i < Ncols; ++i)
	{
		file >> D_temp;
		//cout << "D " << D_temp << endl; 
		this->H_storeCVar(0, i, D_temp);
		//cout << this->getCVar(0, i) << ' ' ;
	}
	//cout << endl;
	file.ignore(256, '\n');
	
	// start reading data...
	fluxData.resize(Nrows * Ncols);
	rowVars.resize(Nrows * NrowVars);
	for (j = 0; j < Nrows; ++j)
	{
		file >> D_temp;
		this->H_storeRVar(j, 0, D_temp);
		// if(j < 10)
		// 	cout << this->getRVar(j, 0) << ' ';

		file >> D_temp;
		this->H_storeRVar(j, 1, D_temp);
		// if(j < 10)
		// 	cout << this->getRVar(j, 1) << ' ';

		for (i = 0; i < Ncols; ++i)
		{
			file >> D_temp;
			this->H_storeFlux(j, i, D_temp);
			// if(i < 2 && j < 10)
			// 	cout << this->getFlux(j, i) << ' ';
		}
		// if(j < 10)
		// 	cout << endl;
	}
	cout << fluxData.size() << endl;
	file.close();
}

/////////////////////////////////////////////////////////////

MEM_iglooAvg::MEM_iglooAvg(string dn)  : MEM_data(dn + "/igloo_avg.txt")
{
	this->H_readFile();
}

MEM_iglooAvg::~MEM_iglooAvg() {}

double MEM_iglooAvg::getFlux_atAngleVel(double alt, double azm, double vel)
{
	return 0.0;
}

void MEM_iglooAvg::H_readFile(void)
{
	int i, j;
	char C_unitName[16];
	stringstream SS_double;
	double D_temp;
	string S_temp;
	NrowVars = 9;
	NcolVars = 1;

	ifstream file;
	file.open(fileName);

	// retrieve the column units
	colVarsUnits.push_back(string("speed (km/s)"));
	//cout << colVarsUnits[0] << endl;

	for (i = 0; i < 7; ++i)
		file.ignore(256, '\n');

	// retrieve the row units
	file.ignore(1);
	for (i = 0; i < NrowVars; ++i)
	{
		file >> S_temp;
		if(i < 3)
			rowVarsUnits.push_back(S_temp);
		else
			rowVarsUnits.push_back(S_temp + " (degrees)");
		//cout << rowVarsUnits[i] << endl;
	}
	
	// retrieve colVars
	colVars.resize(Ncols * NcolVars);
	for (i = 0; i < Ncols; ++i)
	{
		file >> D_temp;
		//cout << "D " << D_temp << endl; 
		this->H_storeCVar(0, i, D_temp);
		//cout << this->getCVar(0, i) << ' ' ;
	}
	//cout << endl;
	file.ignore(256, '\n');
	
	// start reading data...
	fluxData.resize(Nrows * Ncols);
	rowVars.resize(Nrows * NrowVars);
	for (j = 0; j < Nrows; ++j)
	{
		for (i = 0; i < NrowVars; ++i)
		{
			file >> D_temp;
			this->H_storeRVar(j, i, D_temp);
			// if(j < 10)
			// 	cout << this->getRVar(j, i) << ' ';
		}
		

		for (i = 0; i < Ncols; ++i)
		{
			file >> D_temp;
			this->H_storeFlux(j, i, D_temp);
			// if(i < 2 && j < 10)
			// 	cout << this->getFlux(j, i) << ' ';
		}
		// if(j < 10)
		// 	cout << endl;
	}
	cout << fluxData.size() << endl;
	file.close();
}
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MEM_HiDensityCubeAvg::MEM_HiDensityCubeAvg(string dn) : MEM_cubeAvg(dn + "/HiDensity")
{
	this->H_getVelAngleResolution(dn + "/options.txt");
}

MEM_HiDensityCubeAvg::~MEM_HiDensityCubeAvg() {}

/////////////////////////////////////////////////////////////

MEM_HiDensityFluxAvg::MEM_HiDensityFluxAvg(string dn) : MEM_fluxAvg(dn + "/HiDensity")
{
	this->H_getVelAngleResolution(dn + "/options.txt");
}

MEM_HiDensityFluxAvg::~MEM_HiDensityFluxAvg() {}

/////////////////////////////////////////////////////////////

MEM_HiDensityIglooAvg::MEM_HiDensityIglooAvg(string dn) : MEM_iglooAvg(dn + "/HiDensity")
{
	this->H_getVelAngleResolution(dn + "/options.txt");
}

MEM_HiDensityIglooAvg::~MEM_HiDensityIglooAvg() {}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MEM_LoDensityCubeAvg::MEM_LoDensityCubeAvg(string dn) : MEM_cubeAvg(dn + "/LoDensity")
{
	this->H_getVelAngleResolution(dn + "/options.txt");
}

MEM_LoDensityCubeAvg::~MEM_LoDensityCubeAvg() {}

/////////////////////////////////////////////////////////////

MEM_LoDensityFluxAvg::MEM_LoDensityFluxAvg(string dn) : MEM_fluxAvg(dn + "/LoDensity")
{
	this->H_getVelAngleResolution(dn + "/options.txt");
}

MEM_LoDensityFluxAvg::~MEM_LoDensityFluxAvg() {}

/////////////////////////////////////////////////////////////

MEM_LoDensityIglooAvg::MEM_LoDensityIglooAvg(string dn) : MEM_iglooAvg(dn + "/LoDensity")
{
	this->H_getVelAngleResolution(dn + "/options.txt");
}

MEM_LoDensityIglooAvg::~MEM_LoDensityIglooAvg() {}
