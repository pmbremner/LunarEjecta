#ifndef LUNAREJECTA_ASSEMBLY_H
#define LUNAREJECTA_ASSEMBLY_H


#include "lunarEjecta_Headers_Assembly.h"
#include "lunarEjecta_Assembly.h"

#include <string>

using namespace std;

class lunarEjecta_InitializeInterface
{
public:
	lunarEjecta_InitializeInterface(string input_fn);
	~lunarEjecta_InitializeInterface();
	

private:

	/*  For lunarEjecta_Regolith */
	int HH11_targetMaterial;
	int regolithDensType;
	double new_lowDensity;
	double new_avgDensity;
	double new_highDensity;
    double new_escapeSpeed; // m/s

	/*  For ImpactSites_and_ROI */
	double new_ND;     // total number of distance increments
    double new_Nazm;   // total number of azimuth increments
    double new_radius; // radius of Moon
    latLon& new_ROI;   // lat-lon location of Region-Of-Interest

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