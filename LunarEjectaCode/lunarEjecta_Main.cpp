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


	return 0;
}