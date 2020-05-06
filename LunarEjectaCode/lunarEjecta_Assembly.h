#ifndef LUNAREJECTA_ASSEMBLY_H
#define LUNAREJECTA_ASSEMBLY_H

#include "lunarEjecta_Regolith.h"
#include "lunarEjecta_SecondaryFluxData.h"
#include "lunarEjecta_MeteoroidFlux.h"

#include <iostream>
#include <vector>
#include <string>

using namespace std;

// Note: a template class must be definied in the header file
//  There is another solution...
//  https://stackoverflow.com/questions/495021/why-can-templates-only-be-implemented-in-the-header-file?rq=1
template <class genMEMdataHi, class genMEMdataLo, class genOutput>
class lunarEjecta_Assembly
{
public:
	lunarEjecta_Assembly(
		/*  For lunarEjecta_Regolith */
		int HH11_targetMaterial,
		int regolithDensType,
		double new_lowDensity,
		double new_highDensity,

		/*  For ImpactSites_and_ROI */
		double new_ND,     // total number of distance increments
        double new_Nazm,   // total number of azimuth increments
        double new_radius, // radius of Moon
        latLon& new_ROI,   // lat-lon location of Region-Of-Interest

		/* For MEM_LatData */
		string dn,   // directory name of lat data
		double lMin, // Minimum lat (most likely -90)
		double lMax, // Maximum lat (most likely +90)
		int NL,      // Assumes lat data equally spaced with NL # of lat files

		/* For SecondaryFluxData */
		string fn, // file name
		double new_xMin, // min of x-axis of integral flux
		double new_xMax, // max of x-axis of integral flux
		int new_Nx,          // number of spacings on x-axis
		int new_xScale,      // xScaleType = linear or log10
		int new_NSetsXY,           // number of sets of x-y data, if 0 will ignore setMin and setMax
		vector<double> new_setMin, // minimum of set range i
		vector<double> new_setMax) // maximum of set range i
	{
		cout << "----------------------------------------\n";
		cout << "lunarEjecta_Assembly template class init \n";
		cout << "----------------------------------------\n";

		// init regolith
		RegolithProperties = new lunarEjecta_Regolith(HH11_targetMaterial, regolithDensType, new_lowDensity, new_highDensity);

		// init site and ROI locations
		ImpactSitesROILoc = new ImpactSites_and_ROI(new_ND, new_Nazm, new_radius, new_ROI);

		// establish high density and low density MEM input data
		MEMLatDataHi = new MEM_LatData<genMEMdataHi>(dn, lMin, lMax, NL);
		MEMLatDataLo = new MEM_LatData<genMEMdataLo>(dn, lMin, lMax, NL);

		// establish secondary flux data
		SecFluxOutputData = new SecondaryFluxData<genOutput>(fn, new_xMin, new_xMax, new_Nx, new_xScale, new_NSetsXY, new_setMin, new_setMax);
	}

	~lunarEjecta_Assembly() {
		delete RegolithProperties;
		delete ImpactSitesROILoc;
		delete MEMLatDataHi;
		delete MEMLatDataLo;
		delete SecFluxOutputData;
	}

	void computeSecondaryFlux() { // All the magic happens here!
		int i_siteDist, j_siteAzm;

		// For each location, outer loop distance, inner loop azm
		for (i_siteDist = 0; i_siteDist < ImpactSitesROILoc->getND(); ++i_siteDist)
		{
			// all distance dependent only terms should be computed here


			for (j_siteAzm = 0; j_siteAzm < ImpactSitesROILoc->getNazm(); ++j_siteAzm)
			{
				
			}
		}
	}

private:

	double H_compH11RegDensFactor(int lowHigh){
		double dens = 0.0;
		switch(lowHigh){
			case 0: // low
				dens = RegolithProperties->getlowDensity();
				break;
			case 1: // high
				dens = RegolithProperties->gethighDensity();
				break;
			default:
				cerr << "ERROR: H_compH11RegDensFactor invalid density type selection\n";
		}
		return pow(dens, -(3.*RegolithProperties->getHH11_nu()-1.));
	}

	double H_compH11HiDensFactor() {
		double mass_sum = 0.0, cur_dens;

		for (int i = 0; i < MEMLatDataHi->getNdens(); ++i)
		{
			cur_dens = (MEMLatDataHi->getdensLEdge(i) + MEMLatDataHi->getdensREdge(i)) / 2.0; // kg/m^3
			mass_sum += MEMLatDataHi->getdensFraction * pow(cur_dens, 3.*RegolithProperties->getHH11_nu()-1.) * 50.; // bins are always 50 per kg/m^3 large
		}
		return mass_sum;
	}

	double H_compH11LiDensFactor() {
		double mass_sum = 0.0, cur_dens;

		for (int i = 0; i < MEMLatDataLo->getNdens(); ++i)
		{
			cur_dens = (MEMLatDataLo->getdensLEdge(i) + MEMLatDataLo->getdensREdge(i)) / 2.0; // kg/m^3
			mass_sum += MEMLatDataLo->getdensFraction * pow(cur_dens, 3.*RegolithProperties->getHH11_nu()-1.) * 50.; // bins are always 50 per kg/m^3 large
		}
		return mass_sum;
	}

	lunarEjecta_Regolith*         RegolithProperties;
	ImpactSites_and_ROI*          ImpactSitesROILoc;
	MEM_LatData<genMEMdataHi>*    MEMLatDataHi;
	MEM_LatData<genMEMdataLo>*    MEMLatDataLo;
	SecondaryFluxData<genOutput>* SecFluxOutputData;
};


#endif 