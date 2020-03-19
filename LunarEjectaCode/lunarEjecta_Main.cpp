#include "lunarEjecta_Headers.h"

using namespace std;

int main(int argc, char const *argv[])
{
	/* code */

	//MEM_data test_data("../LatRunData/lat0/HiDensity/igloo_avg.txt"); // cube_avg.txt, flux_avg.txt
	
	MEM_cubeAvg cube_data("../LatRunData/lat0/HiDensity");

	MEM_fluxAvg flux_data("../LatRunData/lat0/HiDensity");

	MEM_iglooAvg igloo_data("../LatRunData/lat0/HiDensity");

	return 0;
}