#include "lunarEjecta_Headers.h"

using namespace std;

#include <string>

int main(int argc, char const *argv[])
{
	string parameter_filename(argv[1]);
	string run_filename(argv[2]);

	lunarEjecta_InitializeInterface<MEM_HiDensityIglooAvg,
						 MEM_LoDensityIglooAvg,
						 MassLimitedIglooIntegratedFlux> setupModel(parameter_filename, run_filename);

	setupModel.commitInit();
	setupModel.computeSecondaryFlux();
	return 0;
}