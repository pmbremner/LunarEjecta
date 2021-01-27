#include "lunarEjecta_GeneralExpressions.h"
#include "lunarEjecta_Regolith.h"

#include <iostream>
#include <cmath>

using namespace std;


lunarEjecta_Regolith::lunarEjecta_Regolith
                      (int HH11_targetMaterial,
					   int regolithDensType,
					   double new_lowDensity,
					   double new_avgDensity,
					   double new_highDensity,
					   double new_radius,
					   double new_escapeSpeed) {
	cout << " Regolith Cratering Properties: \n";

	lunarRadius = new_radius;
	lunarEscapeSpeed = new_escapeSpeed;

	cout << "lunarRadius = " << lunarRadius << " m\n";
	cout << "lunarEscapeSpeed = " << lunarEscapeSpeed << " m/s\n";

	HH11_nu = 0.4; // see footnote 5 of Housen Holsapple 2011

	switch(HH11_targetMaterial) {
		case rock:
			cout << "rock";
			HH11_porosity = 0.0;
			HH11_mu       = 0.55;
			HH11_C1       = 1.5;
			HH11_k        = 0.3;
			HH11_C4       = this->H_calcHH11_C4();
			break;
		case weaklyCementedBasalt:
			cout << "weaklyCementedBasalt";
			HH11_porosity = 0.20;
			HH11_mu       = 0.46;
			HH11_C1       = 0.18;
			HH11_k        = 0.3;
			HH11_C4       = this->H_calcHH11_C4();
			break;
		case sand:
			cout << "sand";
			HH11_porosity = 0.35;
			HH11_mu       = 0.41;
			HH11_C1       = 0.55;
			HH11_k        = 0.3;
			HH11_C4       = this->H_calcHH11_C4();
			break;
		case glassMicroSpheres:
			cout << "glassMicroSpheres";
			HH11_porosity = 0.36;
			HH11_mu       = 0.45;
			HH11_C1       = 1.0;
			HH11_k        = 0.5;
			HH11_C4       = this->H_calcHH11_C4();
			break;
		case sandFlyAsh:  
			cout << "sandFlyAsh";
			HH11_porosity = 0.45;
			HH11_mu       = 0.4;
			HH11_C1       = 0.55;
			HH11_k        = 0.3;
			HH11_C4       = this->H_calcHH11_C4();
			break;
		case perliteSandMixture:
			cout << "perliteSandMixture";
			HH11_porosity = 0.6;
			HH11_mu       = 0.35;
			HH11_C1       = 0.6;
			HH11_k        = 0.32;
			HH11_C4       = this->H_calcHH11_C4();
			break;
		default:
			cerr << "ERROR: lunarEjecta_Regolith incorrect HH11_targetMaterial\n";
	}
	cout << endl;


	switch(regolithDensType){
		case DSNE:
			// Note, low and high dens taken from Table 3.4.2.3-1 with Table 3.4.2.3.4-1
			// density valid for regolith down to 30 cm
			lowDensity = (1580. - 50.) / (1. - (0.49 - 0.02)); // kg/m^3, DSNE vals

			// Average density for regolith due to Carrier et al 1991, DSNE 3.4.2.3.3
			avgDensity = 3100.;

			// density valid for regolith between 30 - 60 cm
			highDensity = (1740. - 50.) / (1. - (0.44 + 0.02)); // kg/m^3, DSNE vals
			cout << " Regolith Density Model: DSNE\n";
			cout << "  lowDensity = " << lowDensity << " kg/m^3" << endl;
			cout << "  avgDensity = " << avgDensity << " kg/m^3" << endl;
			cout << "  highDensity = " << highDensity << " kg/m^3" << endl;
			break;

		case UserDefRegDens:
			lowDensity  = new_lowDensity;
			avgDensity  = new_avgDensity;
			highDensity = new_highDensity;
			cout << " Regolith Density Model: UserDefRegDens\n";
			cout << "  lowDensity = " << lowDensity << " kg/m^3" << endl;
			cout << "  avgDensity = " << avgDensity << " kg/m^3" << endl;
			cout << "  highDensity = " << highDensity << " kg/m^3" << endl;
			break;
		default:
			cerr << "ERROR: lunarEjecta_Regolith incorrect regolithDensType\n";
	}
	

}

lunarEjecta_Regolith::~lunarEjecta_Regolith() {}

inline double lunarEjecta_Regolith::getHH11_porosity()    {return HH11_porosity;}
       double lunarEjecta_Regolith::getHH11_mu()          {return HH11_mu;}
inline double lunarEjecta_Regolith::getHH11_C1()          {return HH11_C1;}
inline double lunarEjecta_Regolith::getHH11_k()           {return HH11_k;}
	   double lunarEjecta_Regolith::getHH11_C4()          {return HH11_C4;}
       double lunarEjecta_Regolith::getHH11_nu()          {return HH11_nu;}
	   double lunarEjecta_Regolith::getlowDensity()       {return lowDensity;}
	   double lunarEjecta_Regolith::getavgDensity()       {return avgDensity;}
	   double lunarEjecta_Regolith::gethighDensity()      {return highDensity;}
	   double lunarEjecta_Regolith::getlunarRadius()      {return lunarRadius;}
	   double lunarEjecta_Regolith::getlunarEscapeSpeed() {return lunarEscapeSpeed;}

inline double lunarEjecta_Regolith::H_calcHH11_C4() {
	return 3.*HH11_k / (4.*PI) * pow(HH11_C1, 3.*HH11_mu);
}