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
	inline double HH_DGrunM(double m) { // m = units of g, output in units of 1/(m^2*yr)
		return pow(2.2E3 * pow(m, 0.306) + 15., -4.38)
		     + 1.3E-9 * pow(m + 1.E11 * pow(m, 2.) + 1.E27 * pow(m, 4.), -0.36)
		     + 1.3E-16 * pow(m + 1.E6 * pow(m, 2.), -0.85); //## NEED TO MAKE DG, not G
	}

	double H_compH11GrunMassFactor(int lowHigh, int N = 14) {
		int i;
		double GrunNormalization, massMin, massMax = 10.;
		double Ai, Bi, sum = 0.;
		vector<double> m;
		vector<double> f;
		f.resize(N);

		switch(lowHigh){
			case 0: // low density MEM pop
				massMin = MEMLatDataLo->getGrunMinMass();
				GrunNormalization = HH_DGrunM(massMin);
				break;
			case 1: // high density MEM pop
				massMin = MEMLatDataHi->getGrunMinMass();
				GrunNormalization = HH_DGrunM(massMin);
				break;
			default:
				cerr << "ERROR: H_compH11GrunMassFactor invalid density type selection\n";
		}

		// prepare x and y vectors for power law fits
		logspace(m, massMin, massMax, N, 0, N);
		for (i = 0; i < N; ++i)
			f[i] = HH_DGrunM(m[i]) / GrunNormalization;

		// compute integral by breaking up in log spacings and approx by power law
		for (i = N-2; i >= 0; i--) // start at low end of mass to sum smaller #'s first
		{
			Bi = log10(f[i+1] / f[i]) / log10(m[i+1] / m[i]); 
			Ai = (f[i] + f[i+1]) / (pow(m[i], Bi) + pow(m[i+1], Bi));

			sum += Ai/(Bi+1.) * (pow(m[i+1], Bi+1) - pow(m[i], Bi+1));
		}
		return sum;
		
	}

	double H_compH11RegDensFactor(int lowHigh){
		double dens = 0.0;
		switch(lowHigh){
			case 0: // low density regolith
				dens = RegolithProperties->getlowDensity();
				break;
			case 1: // high density regolith
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

	double H_compH11LoDensFactor() {
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