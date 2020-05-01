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
		/* For MEM_LatData*/
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
		cout << "lunarEjecta_Assembly template class init \n";
		// establish high density and low density MEM input data
		MEMLatDataHi = new MEM_LatData<genMEMdataHi>(dn, lMin, lMax, NL);
		MEMLatDataLo = new MEM_LatData<genMEMdataLo>(dn, lMin, lMax, NL);

		// establish secondary flux data
		SecFluxOutputData = new SecondaryFluxData<genOutput>(fn, new_xMin, new_xMax, new_Nx, new_xScale, new_NSetsXY, new_setMin, new_setMax);
	}

	~lunarEjecta_Assembly() {
		delete MEMLatDataHi;
		delete MEMLatDataLo;
		delete SecFluxOutputData;
	}

private:
	MEM_LatData<genMEMdataHi>* MEMLatDataHi;
	MEM_LatData<genMEMdataLo>* MEMLatDataLo;
	SecondaryFluxData<genOutput>* SecFluxOutputData;
};


#endif 