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
	cout << HiD_cube_data.getFlux_atAngleVel(0.0, 0.0, 3.0) << endl;
	cout << HiD_cube_data.getFlux_atAngleVel(0.0, 180.0, 3.0) << endl;
	cout << HiD_cube_data.getFlux_atAngleVel(0.0, 90.0, 3.0) << endl;
	cout << HiD_cube_data.getFlux_atAngleVel(0.0, 270.0, 3.0) << endl;
	cout << HiD_cube_data.getFlux_atAngleVel(90.0, 0.0, 3.0) << endl;

	MEM_HiDensityFluxAvg  HiD_flux_data("../LatRunData/lat0");
	cout << HiD_flux_data.getFlux_atAngleVel(5.0, 345, 3.0) << endl;
	cout << HiD_flux_data.getFlux_atAngleVel(5.0, 350, 3.0) << endl;
	cout << HiD_flux_data.getFlux_atAngleVel(5.0, 355, 3.0) << endl;

	// MEM_HiDensityIglooAvg HiD_igloo_data("../LatRunData/lat0");

	// MEM_LoDensityCubeAvg  LoD_cube_data("../LatRunData/lat0");
	// MEM_LoDensityFluxAvg  LoD_flux_data("../LatRunData/lat0");
	// MEM_LoDensityIglooAvg LoD_igloo_data("../LatRunData/lat0");

	//MEM_LatData<MEM_HiDensityFluxAvg> lat_data("../LatRunData", -90.0, 90.0, 37); //"../LatRunData", -90.0, 90.0, 37
	//lat_data().print();


	return 0;
}