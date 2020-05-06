#include "lunarEjecta_GeneralExpressions.h"
#include "lunarEjecta_Regolith.h"

#include <iostream>
#include <cmath>

using namespace std;


lunarEjecta_Regolith::lunarEjecta_Regolith(int HH11_targetMaterial) {
	switch(HH11_targetMaterial) {
		case rock:
			HH11_porosity = 0.0;
			HH11_mu       = 0.55;
			HH11_C1       = 1.5;
			HH11_k        = 0.3;
			HH11_C4       = this->H_calcHH11_C4();
			break;
		case weaklyCementedBasalt:
			HH11_porosity = 0.20;
			HH11_mu       = 0.46;
			HH11_C1       = 0.18;
			HH11_k        = 0.3;
			HH11_C4       = this->H_calcHH11_C4();
			break;
		case sand:
			HH11_porosity = 0.35;
			HH11_mu       = 0.41;
			HH11_C1       = 0.55;
			HH11_k        = 0.3;
			HH11_C4       = this->H_calcHH11_C4();
			break;
		case glassMicroSpheres:
			HH11_porosity = 0.36;
			HH11_mu       = 0.45;
			HH11_C1       = 1.0;
			HH11_k        = 0.5;
			HH11_C4       = this->H_calcHH11_C4();
			break;
		case sandFlyAsh:
			HH11_porosity = 0.45;
			HH11_mu       = 0.4;
			HH11_C1       = 0.55;
			HH11_k        = 0.3;
			HH11_C4       = this->H_calcHH11_C4();
			break;
		case perliteSandMixture:
			HH11_porosity = 0.6;
			HH11_mu       = 0.35;
			HH11_C1       = 0.6;
			HH11_k        = 0.32;
			HH11_C4       = this->H_calcHH11_C4();
			break;
		default:
			cerr << "ERROR: lunarEjecta_Regolith incorrect HH11_targetMaterial\n";
	}

	// density valid for regolith down to 30 cm
	lowDensity = 1580.; // kg/m^3, DSNE vals

	// density valid for regolith between 30 - 60 cm
	highDensity = 1740.; // kg/m^3, DSNE vals

}

lunarEjecta_Regolith::~lunarEjecta_Regolith() {}


inline double lunarEjecta_Regolith::H_calcHH11_C4() {
	return 3.*HH11_k / (4.*PI) * pow(HH11_C1, 3.*HH11_mu);
}