// Declaration of regolith information

#ifndef LUNAREJECTA_REGOLITH_H
#define LUNAREJECTA_REGOLITH_H

using namespace std;

// HH11_targetMaterial
enum HH11_Table3Params {rock, weaklyCementedBasalt, sand, glassMicroSpheres, sandFlyAsh, perliteSandMixture};




class lunarEjecta_Regolith
{
public:
	lunarEjecta_Regolith(int HH11_targetMaterial);
	~lunarEjecta_Regolith();
	

private:
	inline double H_calcHH11_C4();

	double HH11_porosity;
	double HH11_mu;
	double HH11_C1;
	double HH11_k;
	double HH11_C4;

	double lowDensity;  // kg/m^3
	double highDensity; // kg/m^3
};
#endif 