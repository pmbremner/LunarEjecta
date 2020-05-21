#include "lunarEjecta_Headers.h"

using namespace std;

int main(int argc, char const *argv[])
{
	/* code */

	//MEM_data test_data("../LatRunData/lat0/HiDensity/igloo_avg.txt"); // cube_avg.txt, flux_avg.txt
	
	// MEM_cubeAvg cube_data("../LatRunData/lat0/HiDensity");

	// MEM_fluxAvg flux_data("../LatRunData/lat0/HiDensity");

	// MEM_iglooAvg igloo_data("../LatRunData/lat0/HiDensity");


	MEM_HiDensityCubeAvg  HiD_cube_data("../LatRunData/lat0");
	cout << HiD_cube_data.getFlux_atAngleVel(0.0, 0.0, 5.6) << endl;
	cout << HiD_cube_data.getFlux_atAngleVel(0.0, 180.0, 5.6) << endl;
	cout << HiD_cube_data.getFlux_atAngleVel(0.0, 90.0, 5.6) << endl;
	cout << HiD_cube_data.getFlux_atAngleVel(0.0, 270.0, 5.6) << endl;
	cout << HiD_cube_data.getFlux_atAngleVel(90.0, 0.0, 5.6) << endl;

	MEM_HiDensityFluxAvg  HiD_flux_data("../LatRunData/lat0");
	cout << HiD_flux_data.getFlux_atAngleVel(0.0, 0.0, 5.6) << endl;
	cout << HiD_flux_data.getFlux_atAngleVel(0.0, 180.0, 5.6) << endl;
	cout << HiD_flux_data.getFlux_atAngleVel(0.0, 90.0, 5.6) << endl;
	cout << HiD_flux_data.getFlux_atAngleVel(0.0, 270.0, 5.6) << endl;
	cout << HiD_flux_data.getFlux_atAngleVel(90.0, 0.0, 5.6) << endl;

	MEM_HiDensityIglooAvg HiD_igloo_data("../LatRunData/lat0");
	cout << HiD_igloo_data.getFlux_atAngleVel(0.0, 0.0, 5.6) << endl;
	cout << HiD_igloo_data.getFlux_atAngleVel(0.0, 180.0, 5.6) << endl;
	cout << HiD_igloo_data.getFlux_atAngleVel(0.0, 90.0, 5.6) << endl;
	cout << HiD_igloo_data.getFlux_atAngleVel(0.0, 270.0, 5.6) << endl;
	cout << HiD_igloo_data.getFlux_atAngleVel(90.0, 0.0, 5.6) << endl;

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
	
	SecondaryFluxData<MassLimitedIntegralFluxVsMass> fluxSecM("run0", 8.E-13, 0.8, 36, log10Scale, 0, emptyVector, emptyVector);

	cout << endl << endl;

	int Nv = 11;
	vector<double> vMin, vMax;
	logspaceBinEdges(vMin, vMax, 0.1, 2.4, Nv);

	lunarEjecta_Assembly<MEM_HiDensityIglooAvg,
						 MEM_LoDensityIglooAvg,
						 MassLimitedIglooIntegratedFlux>
		lunarEjecta(sandFlyAsh, DSNE, 0.0, 0.0, 0.0, 
			        50, 72, 1737.1E3, MN_ll,
		            "../LatRunData", -90.0, 90.0, 2, // 37
					"run0", 0., 1., 5, log10Scale, Nv, vMin, vMax);

	return 0;
} 