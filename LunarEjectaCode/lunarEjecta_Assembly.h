#ifndef LUNAREJECTA_ASSEMBLY_H
#define LUNAREJECTA_ASSEMBLY_H

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

		// init site and ROI locations
		ImpactSitesROILoc = new ImpactSites_and_ROI(new_ND, new_Nazm, new_radius, new_ROI);

		// establish high density and low density MEM input data
		MEMLatDataHi = new MEM_LatData<genMEMdataHi>(dn, lMin, lMax, NL);
		MEMLatDataLo = new MEM_LatData<genMEMdataLo>(dn, lMin, lMax, NL);

		// establish secondary flux data
		SecFluxOutputData = new SecondaryFluxData<genOutput>(fn, new_xMin, new_xMax, new_Nx, new_xScale, new_NSetsXY, new_setMin, new_setMax);
	}

	~lunarEjecta_Assembly() {
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


			// Separately, we need to integrate over a finer azm grid in order to 
			//  figure out the normalization for the small azm wedge,
			//  since the exponent 'a' depends on the azm direction
			// Will have to integrate the whole speed-x(angle) between the curves
			//  for each azm wedge

			for (j_siteAzm = 0; j_siteAzm < ImpactSitesROILoc->getNazm(); ++j_siteAzm)
			{
				
			}
		}
	}

private:
	ImpactSites_and_ROI*          ImpactSitesROILoc;
	MEM_LatData<genMEMdataHi>*    MEMLatDataHi;
	MEM_LatData<genMEMdataLo>*    MEMLatDataLo;
	SecondaryFluxData<genOutput>* SecFluxOutputData;
};


#endif 