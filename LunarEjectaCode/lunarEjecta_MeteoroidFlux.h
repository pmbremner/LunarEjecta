// Declaration of meteoroid flux information

// Have the option to read in cube_avg, flux_avg, and igloo_avg

#ifndef LUNAREJECTA_METEOROIDFLUX_H
#define LUNAREJECTA_METEOROIDFLUX_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


enum cubeDirection {p_xRam, m_xWake, p_yPort, m_yStarboard, p_zZenith, m_zNadir};

enum iglooRowVars {ID, I, J, PHI1, PHI2, THETA1, THETA2, PHIavg, THETAavg};

// This is an abstract base class and cannot be used directly
class MEM_data
{
public:
	MEM_data(string fn);
	~MEM_data();

	virtual double getFlux_atAngleVel(double alt, double azm, double vel) = 0;

	inline double getFlux(int row, int col);
	inline double getRVar(int row, int col);
	inline double getCVar(int row, int col);

	inline double getvMin();
	inline double getvMax();

	inline double getNrows();
	inline double getNcols();
	inline double getNrowVars();
	inline double getNcolVars();

	inline double getNdens();
	inline double getdensLEdge(int idx);
	inline double getdensREdge(int idx);
	inline double getdensFraction(int idx);

	inline double getGrunMinMass();

protected:
	virtual void H_readFile(void) = 0;
	void H_setHeaderLength(int h);
	void H_getRowCol_FromFile();
	void H_getVelAngleResolution(string fn); // used two levels up
	void H_getGrunMinMass(string fn);  // used two levels up
	void H_readDensityFile(string fn); // used two levels up
	void H_readInt_FromFile(ifstream& file, int& firstInt);
	inline int H_idxFlux(int row, int col);
	inline int H_idxRVar(int row, int col);
	inline int H_idxCVar(int row, int col);
	inline void H_storeFlux(int row, int col, double flux);
	inline void H_storeRVar(int row, int col, double RVar);
	inline void H_storeCVar(int row, int col, double RCar);

	// don't use these now
	inline void H_pushBackFlux(double flux); 
	inline void H_pushBackRVar(double RVar);
	inline void H_pushBackCVar(double CVar); 

	bool iotype; // 1 = input, 0 = output
	string fileName;
	int headerLength;

	double vMin;
	double vMax;
	int dVel;   // = 1, 2
	int dAngle; // = 1, 2, 3, 4, 5
	double log10GrunMinMass; // between [-6, 1]

	int Nrows; // of data only
	int Ncols; // of data only
	int NrowVars;
	int NcolVars;
	int Nvel;
	int Ntheta; // azm
	int Nphi;   // alt

	vector<double> fluxData; // total size = Nrows * Ncols
	vector<double> rowVars;  // total size = Nrows * NrowVars
	vector<double> colVars;  // total size = Ncols * NcolVars

	string fluxUnits;            
	vector<string> rowVarsUnits; 
	vector<string> colVarsUnits;

	// density data (not init here, but in 2 gens down)
	//  init in H_readDensityFile
	int Ndens;
	vector<double> densLEdge; // kg/m^3
	vector<double> densREdge; // kg/m^3
	vector<double> densFraction; // units of per 50 kg/m^3
};

/////////////////////////////////////////////////////////////
class MEM_cubeAvg : public MEM_data
{
public:
	MEM_cubeAvg(string dn); // dn = directory name
	~MEM_cubeAvg();
	double getFlux_atAngleVel(double alt, double azm, double vel);

private:
	double H_cubeFluxAvg_atRow(double x, double y, double z, int row);
	void H_readFile();

};
/////////////////////////////////////////////////////////////
class MEM_fluxAvg : public MEM_data
{
public:
	MEM_fluxAvg(string dn);
	~MEM_fluxAvg();
	double getFlux_atAngleVel(double alt, double azm, double vel);

private:
	void H_readFile();
};
/////////////////////////////////////////////////////////////
class MEM_iglooAvg : public MEM_data
{
public:
	MEM_iglooAvg(string dn);
	~MEM_iglooAvg();
	double getFlux_atAngleVel(double alt, double azm, double vel);

private:
	void H_readFile();
};
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
class MEM_HiDensityCubeAvg : public MEM_cubeAvg
{
public:
	MEM_HiDensityCubeAvg(string dn);
	~MEM_HiDensityCubeAvg();
	
};
/////////////////////////////////////////////////////////////
class MEM_HiDensityFluxAvg : public MEM_fluxAvg
{
public:
	MEM_HiDensityFluxAvg(string dn);
	~MEM_HiDensityFluxAvg();
	
};
/////////////////////////////////////////////////////////////
class MEM_HiDensityIglooAvg : public MEM_iglooAvg
{
public:
	MEM_HiDensityIglooAvg(string dn);
	~MEM_HiDensityIglooAvg();
	
};
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
class MEM_LoDensityCubeAvg : public MEM_cubeAvg
{
public:
	MEM_LoDensityCubeAvg(string dn);
	~MEM_LoDensityCubeAvg();
	
};
/////////////////////////////////////////////////////////////
class MEM_LoDensityFluxAvg : public MEM_fluxAvg
{
public:
	MEM_LoDensityFluxAvg(string dn);
	~MEM_LoDensityFluxAvg();
	
};
/////////////////////////////////////////////////////////////
class MEM_LoDensityIglooAvg : public MEM_iglooAvg
{
public:
	MEM_LoDensityIglooAvg(string dn);
	~MEM_LoDensityIglooAvg();
	
};
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// Note: a template class must be definied in the header file
template <class genMEMdata> 
class MEM_LatData 
{
public:
	MEM_LatData(string dn, double lMin, double lMax, int NL) // string dn, double lMin, double lMax, int NL
	{
		directoryName = dn;
		latMin        = lMin;
		latMax        = lMax;
		NLat          = NL;
		cout << "MEM_LatData template class init \n";

		dataSet.resize(NL);
		string latDirectoryName;

		for (int i = 0; i < NLat; ++i)
		{
			latDirectoryName = "/lat" + to_string(int(latMin + (latMax-latMin)*i/(NL-1.0)));
			//cout << latDirectoryName << endl;
			dataSet[i] = new genMEMdata(dn + latDirectoryName);
		}
	}

	~MEM_LatData()
	{
		for (int i = 0; i < NLat; ++i)
		{
			delete dataSet[i];
			//cout << "delete " << i << endl;
		}
	}

	inline int getLatIdx(double cur_lat) {
		return int(round((cur_lat - latMin)/(latMax-latMin) * NLat));
	}

	inline double getGrunMinMass(int i = 0) {
		return dataSet[i].getGrunMinMass();
	}

	//void print();

private:
	string directoryName;
	vector<genMEMdata*> dataSet;

	double latMin;
	double latMax;
	int NLat;
};


#endif 