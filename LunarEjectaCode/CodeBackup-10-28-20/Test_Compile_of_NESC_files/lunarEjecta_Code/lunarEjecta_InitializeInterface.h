#ifndef LUNAREJECTA_INITINTERFACE_H
#define LUNAREJECTA_INITINTERFACE_H

#include "lunarEjecta_Assembly.h"
#include "lunarEjecta_Headers_Assembly.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>


using namespace std;

enum DefaultExists {noDefault, isDefault};

template <class genMEMdataHi, class genMEMdataLo, class genOutput>
class lunarEjecta_InitializeInterface
{
public:
	// Here we will read from file the init parameters, the assembly is not made yet
	lunarEjecta_InitializeInterface(string new_input_fn, string new_run_filename)
	{
		initError = 0; 
		input_fn = new_input_fn;
		run_filename = new_run_filename;
		double ROI_Lat = 0.0; // degrees
		double ROI_Lon = 0.0; // degrees

		// set default values
		// will be used if no keyword exists in the file
		// Other keywords must be defined
		//
		/*  For lunarEjecta_Regolith */
		HH11_targetMaterial = sandFlyAsh;
		regolithDensType = DSNE; // DSNE
		new_lowDensity   = 0.0;
		new_avgDensity   = 0.0;
		new_highDensity  = 0.0;
		new_escapeSpeed  = 2375.89; // m/s
		/*  For ImpactSites_and_ROI */
		new_radius       = 1737.1E3; // m
		/* For lunarEjecta_NearEarthObjectFlux */
		new_m_min        = 1.E-2; // kg
		new_m_max        = 1.57E12; // kg
		densType         = defaultDens;  // defaultDens
		userDefDens      = 0.0;
		/* For SecondaryFluxData */
		new_xMin         = 0.;
		new_xMax         = 1.;
		new_xScale       = linearScale; // linearScale
		new_NSetsXY      = 0; // not a default
		/* For lunarEjecta_AdaptiveMesh */
		new_maxLevelMesh    = 7;
		new_maxLevelFractal = 2;

		cout << "---lunarEjecta_InitializeInterface: Reading " << new_input_fn << " ...\n";

		/*  For lunarEjecta_Regolith */
		getParam("HH11_targetMaterial", HH11_targetMaterial, isDefault);
		getParam("regolithDensType"   , regolithDensType   , isDefault);
		getParam("regLowDensity"      , new_lowDensity     , isDefault);
		getParam("regAvgDensity"      , new_avgDensity     , isDefault);
		getParam("regHighDensity"     , new_highDensity    , isDefault);
		getParam("LunarEscapeSpeed"   , new_escapeSpeed    , isDefault);

		/*  For ImpactSites_and_ROI */
		getParam("NDistance"    , new_ND     , noDefault);
		getParam("Nazm"         , new_Nazm   , noDefault);
		getParam("LunarRadius"  , new_radius , isDefault);
		getParam("ROI_Lat"      , ROI_Lat    , isDefault);
		getParam("ROI_Lon"      , ROI_Lon    , isDefault);

		new_ROI = new latLon(ROI_Lat, ROI_Lon);

		/* For MEM_LatData */
		getParam("MEMDirectoryName"    , dn   , noDefault);
		getParam("MEMLatMin"           , lMin , noDefault);
		getParam("MEMLatMax"           , lMax , noDefault);
		getParam("MEMNLat"             , NL   , noDefault);

		/* For lunarEjecta_NearEarthObjectFlux */
		getParam("NEOSpeedFName"    , NEOfn       , noDefault);
		getParam("NEO_m_min"        , new_m_min   , isDefault);
		getParam("NEO_m_max"        , new_m_max   , isDefault);
		getParam("NEODensType"      , densType    , isDefault);
		getParam("NEOuserDefDens"   , userDefDens , isDefault);
		getParam("NEODirectoryName" , dn_NEO      , noDefault);
		getParam("NEOLatMin"        , lMin_NEO    , noDefault);
		getParam("NEOLatMax"        , lMax_NEO    , noDefault);
		getParam("NEONLat"          , NL_NEO      , noDefault);

		/* For SecondaryFluxData */
		//getParam("SecFluxFName" , SecFluxfn   , noDefault);
		getParam("SecxMin"      , new_xMin    , isDefault);
		getParam("SecxMax"      , new_xMax    , isDefault);
		getParam("SecNx"        , new_Nx      , noDefault);
		getParam("SecxScale"    , new_xScale  , isDefault);
		getParam("SecNSetsXY"   , new_NSetsXY , noDefault);
		getParam("SecvMin"      , new_vMin    , noDefault);
		getParam("SecvMax"      , new_vMax    , noDefault);

		if(new_NSetsXY > 0)
			linspaceBinEdges(new_setMin, new_setMax, new_vMin, new_vMax, new_NSetsXY);

		/* For lunarEjecta_AdaptiveMesh */
		getParam("maxLevelMesh"    , new_maxLevelMesh   , isDefault);
		getParam("maxLevelFractal" , new_maxLevelFractal, isDefault);


	}

	// Commit the intialization and make the assembly
	// Will check for "out of bounds" parameters
	void commitInit()
	{
		//cout << initError << endl;
		if(initError)
		// check error flag
		{
			cout << "ERROR: Commit initialization failed...\n";
		}
		else
		{
			cout << "---Commit initialization successful!---\n";
			ejectaModel = new lunarEjecta_Assembly<genMEMdataHi, genMEMdataLo, genOutput>(
				/*  For lunarEjecta_Regolith */
				HH11_targetMaterial,
				regolithDensType,
				new_lowDensity,
				new_avgDensity,
				new_highDensity,
			    new_escapeSpeed, // m/s

				/*  For ImpactSites_and_ROI */
				new_ND,     // total number of distance increments
			    new_Nazm,   // total number of azimuth increments
			    new_radius, // radius of Moon
			    *new_ROI,   // lat-lon location of Region-Of-Interest

				/* For MEM_LatData */
				dn,   // directory name of lat data
				lMin, // Minimum lat (most likely -90)
				lMax, // Maximum lat (most likely +90)
				NL,      // Assumes lat data equally spaced with NL # of lat files

				/* For lunarEjecta_NearEarthObjectFlux */
				NEOfn,
				new_m_min,
				new_m_max,
				densType,
				userDefDens,
				dn_NEO,
				lMin_NEO,
				lMax_NEO,
				NL_NEO,

				/* For SecondaryFluxData */
				run_filename, // file name
				new_xMin, // min of x-axis of integral flux
				new_xMax, // max of x-axis of integral flux
				new_Nx,          // number of spacings on x-axis
				new_xScale,      // xScaleType = linear or log10
				new_NSetsXY,           // number of sets of x-y data, if 0 will ignore setMin and setMax
				new_setMin, // minimum of set range i
				new_setMax, // maximum of set range i
				new_vMin, // km/s
				new_vMax, // km/s

				/* For lunarEjecta_AdaptiveMesh */
				new_maxLevelMesh,    // the division level of the integration mesh
				new_maxLevelFractal); // the division level of the integrand-domain probing)
		}
	}


	void computeSecondaryFlux()
	{
		if(!initError)
			ejectaModel->computeSecondaryFlux();
		else
			cout << "ERROR: Compute Secondary Flux failed --> Need a successful initialization...\n";
	}

	//////////////
	// Functions to modify individual parameters before commiting

	//////////////


	~lunarEjecta_InitializeInterface()
	{
		if(!initError)
			delete ejectaModel;
		delete new_ROI;
	}
	

private:

	lunarEjecta_Assembly<genMEMdataHi, genMEMdataLo, genOutput>* ejectaModel;

	string input_fn;
	string run_filename;
	ifstream input_file;
	bool initError;
	// assumes input_file is open
	template <class paramType>
	void getParam(string paramLabel, paramType& param, bool defaultExists)
	{
		input_file.open(input_fn);
		char C_Line[64];
		paramType tparam;
		bool isError = 1;

		int lineNumber = 1;

		while (input_file.getline(C_Line, 64, '#'))
		{
			string paramLabel_file(C_Line);
			
			if(paramLabel_file == paramLabel)
			{
				if(input_file.getline(C_Line, 64, '#'));
					isError = 0;
				string tstr(C_Line);
				stringstream SS_param(tstr);
				SS_param >> param;
				//input_file.ignore(256, '\n'); // temp
				if(defaultExists)
					cout << "  Overriding Default value...\n";
				cout << "Line " << lineNumber << " | " << paramLabel << " = " << param << endl;
				//return;
			}
			else
			{
				input_file.ignore(256, '\n');
			}

			lineNumber++;
		}
		// if there's an error, set the initError flag
		if(!defaultExists && isError){
			initError = 1;
			cout << "ERROR: Parameter '" << paramLabel << "' not found...\n";
		}

		input_file.close();
	}

	// These are the inputs to init the assembly class //

	/*  For lunarEjecta_Regolith */
	int HH11_targetMaterial;
	int regolithDensType;
	double new_lowDensity;
	double new_avgDensity;
	double new_highDensity;
    double new_escapeSpeed; // m/s

	/*  For ImpactSites_and_ROI */
	int new_ND;     // total number of distance increments
    int new_Nazm;   // total number of azimuth increments
    double new_radius; // radius of Moon
    latLon* new_ROI;   // lat-lon location of Region-Of-Interest

	/* For MEM_LatData */
	string dn;   // directory name of lat data
	double lMin; // Minimum lat (most likely -90)
	double lMax; // Maximum lat (most likely +90)
	int NL;      // Assumes lat data equally spaced with NL # of lat files

	/* For lunarEjecta_NearEarthObjectFlux */
	string NEOfn;
	double new_m_min;
	double new_m_max;
	int densType;
	double userDefDens;
	string dn_NEO;
	double lMin_NEO;
	double lMax_NEO;
	int NL_NEO;

	/* For SecondaryFluxData */
	string SecFluxfn; // file name
	double new_xMin; // min of x-axis of integral flux
	double new_xMax; // max of x-axis of integral flux
	int new_Nx;          // number of spacings on x-axis
	int new_xScale;      // xScaleType = linear or log10
	int new_NSetsXY;           // number of sets of x-y data; if 0 will ignore setMin and setMax
	vector<double> new_setMin; // minimum of set range i
	vector<double> new_setMax; // maximum of set range i
	double new_vMin; // km/s
	double new_vMax; // km/s

	/* For lunarEjecta_AdaptiveMesh */
	int new_maxLevelMesh;    // the division level of the integration mesh
	int new_maxLevelFractal; // the division level of the integrand-domain probing

};


#endif 