#include "lunarEjecta_Headers.h"

using namespace std;

int main(int argc, char const *argv[])
{
	/* code */

	//MEM_data test_data("../LatRunData/lat0/HiDensity/igloo_avg.txt"); // cube_avg.txt, flux_avg.txt
	
	// MEM_cubeAvg cube_data("../LatRunData/lat0/HiDensity");

	// MEM_fluxAvg flux_data("../LatRunData/lat0/HiDensity");

	// MEM_iglooAvg igloo_data("../LatRunData/lat0/HiDensity");


	// MEM_HiDensityCubeAvg  HiD_cube_data("../LatRunData/lat0");
	// cout << HiD_cube_data.getFlux_atAngleVel(0.0, 0.0, 5.6) << endl;
	// cout << HiD_cube_data.getFlux_atAngleVel(0.0, 180.0, 5.6) << endl;
	// cout << HiD_cube_data.getFlux_atAngleVel(0.0, 90.0, 5.6) << endl;
	// cout << HiD_cube_data.getFlux_atAngleVel(0.0, 270.0, 5.6) << endl;
	// cout << HiD_cube_data.getFlux_atAngleVel(90.0, 0.0, 5.6) << endl;

	// MEM_HiDensityFluxAvg  HiD_flux_data("../LatRunData/lat0");
	// cout << HiD_flux_data.getFlux_atAngleVel(0.0, 0.0, 5.6) << endl;
	// cout << HiD_flux_data.getFlux_atAngleVel(0.0, 180.0, 5.6) << endl;
	// cout << HiD_flux_data.getFlux_atAngleVel(0.0, 90.0, 5.6) << endl;
	// cout << HiD_flux_data.getFlux_atAngleVel(0.0, 270.0, 5.6) << endl;
	// cout << HiD_flux_data.getFlux_atAngleVel(90.0, 0.0, 5.6) << endl;

	// MEM_HiDensityIglooAvg HiD_igloo_data("../LatRunData/lat0");
	// cout << HiD_igloo_data.getFlux_atAngleVel(0.0, 0.0, 5.6) << endl;
	// cout << HiD_igloo_data.getFlux_atAngleVel(0.0, 180.0, 5.6) << endl;
	// cout << HiD_igloo_data.getFlux_atAngleVel(0.0, 90.0, 5.6) << endl;
	// cout << HiD_igloo_data.getFlux_atAngleVel(0.0, 270.0, 5.6) << endl;
	// cout << HiD_igloo_data.getFlux_atAngleVel(90.0, 0.0, 5.6) << endl;

	// MEM_LoDensityCubeAvg  LoD_cube_data("../LatRunData/lat0");
	// MEM_LoDensityFluxAvg  LoD_flux_data("../LatRunData/lat0");
	// MEM_LoDensityIglooAvg LoD_igloo_data("../LatRunData/lat0");

	//MEM_LatData<MEM_HiDensityFluxAvg> lat_data("../LatRunData", -90.0, 90.0, 37); //"../LatRunData", -90.0, 90.0, 37
	//lat_data().print();

	latLon MN_ll(45.5679, -93.593);
	latLon AL_ll(34.937286,  -86.828090);

	MN_ll.dispLatLon();
	AL_ll.dispLatLon();

	MN_ll.dispNormDistTo(AL_ll);

	MN_ll.dispBearingInitial(AL_ll);
	MN_ll.dispBearingFinal(AL_ll);
	MN_ll.dispAzmInitial(AL_ll);
	MN_ll.dispAzmFinal(AL_ll);

	AL_ll.dispBearingInitial(MN_ll);
	AL_ll.dispBearingFinal(MN_ll);
	AL_ll.dispAzmInitial(MN_ll);
	AL_ll.dispAzmFinal(MN_ll);

	cout << endl << endl;
	ImpactSites_and_ROI siteList(3, 3, 1737.1E3, MN_ll);

	vector<double> emptyVector;
	MassLimitedIntegralFluxVsMass fluxM(8.E-13, 0.8, 36, log10Scale, 0, emptyVector, emptyVector);
	SizeLimitedIntegralFluxVsSpeed fluxD(1.E-6, 1.E-2, 36, log10Scale, 0, emptyVector, emptyVector);
	//MassLimitedIglooIntegratedFlux fluxIgloo(1.E-2, 2.5, 36, linearScale, 0, emptyVector, emptyVector);
	
	//SecondaryFluxData<MassLimitedIntegralFluxVsMass> fluxSecM("run0", 8.E-13, 0.8, 36, log10Scale, 0, emptyVector, emptyVector);

	cout << endl << endl;

	latLon Equator_ll(0.0, 0.0);
	latLon Polar_ll(89.0, 0.0);
	latLon Lat45_ll(45.0, 0.0);

	int Nv = 13;
	double vMin = 0.1; // km/s
	double vMax = 2.4; // km/s
	vector<double> vLeft, vRight;
	//logspaceBinEdges(vLeft, vRight, 0.1, 2.4, Nv);
	linspaceBinEdges(vLeft, vRight, vMin, vMax, Nv);

	lunarEjecta_Assembly<MEM_HiDensityIglooAvg,
						 MEM_LoDensityIglooAvg,
						 MassLimitedIglooIntegratedFlux>
		lunarEjecta(
			// sandFlyAsh, DSNE, 0.0, 0.0, 0.0, 2375.89,
			//         21, 12, 1737.1E3, MN_ll,
		 //            "../LatRunData", -90.0, 90.0, 37, // 37
		 //            "../NEA_Brown/vmass.txt", 1.E-2, 1.57E12, defaultDens, 0.,
		 //            "../LatNEOData", -90.0, 90.0, 37,
			// 		"run4.txt", 0., 1., 5, linearScale, Nv, vLeft, vRight,
			// 		/*51, 51,*/ 3, 1);


		// /*  For lunarEjecta_Regolith */
		sandFlyAsh, // int HH11_targetMaterial,
		DSNE, // int regolithDensType,
		0.0, // double new_lowDensity,
		0.0, // double new_avgDensity,
		0.0, // double new_highDensity,
	    2375.89, // double new_escapeSpeed, units of m/s

		// /*  For ImpactSites_and_ROI */
		21,//21, // double new_ND,     // total number of distance increments
        15,//15, // double new_Nazm,   // total number of azimuth increments
        1737.1E3, // double new_radius, // radius of Moon, units of m
        Equator_ll,//Polar_ll,//Equator_ll,Lat45_ll    // latLon& new_ROI,   // lat-lon location of Region-Of-Interest

		// /* For MEM_LatData */
		"../LatRunData", // string dn,   // directory name of lat data
		-90.0, // double lMin, // Minimum lat (most likely -90)
		90.0, // double lMax, // Maximum lat (most likely +90)
		37, // int NL,      // Assumes lat data equally spaced with NL # of lat files

		// /* For lunarEjecta_NearEarthObjectFlux */
		"../NEA_Brown/vmass.txt", // string NEOfn,
		1.E-2, // double new_m_min, units of kg
		1.57E12, // double new_m_max, units of kg
		defaultDens, // int densType,
		0., // double userDefDens,
		"../LatNEOData", // string dn_NEO,
		-90.0, // double lMin_NEO,
		90.0, // double lMax_NEO,
		37, // int NL_NEO,

		// /* For SecondaryFluxData */
		"run_equator_A3.txt", // string fn, // file name
		0., // double new_xMin, // min of x-axis of integral flux
		1., // double new_xMax, // max of x-axis of integral flux
		5, // int new_Nx,          // number of spacings on x-axis (will be resolution for igloo)
		linearScale, // int new_xScale,      // xScaleType = linear or log10
		Nv, // int new_NSetsXY,           // number of sets of x-y data, if 0 will ignore setMin and setMax
		vLeft, // vector<double> new_setMin, // minimum of set range i
		vRight, // vector<double> new_setMax, // maximum of set range i
		vMin, // double new_vMin,
		vMax, // double new_vMax,


		// /* For lunarEjecta_AdaptiveMesh */
		3, // int new_maxLevelMesh,    // the division level of the integration mesh
		1);// int new_maxLevelFractal // the division level of the integrand-domain probing



	lunarEjecta.computeSecondaryFlux();

	return 0;
} 