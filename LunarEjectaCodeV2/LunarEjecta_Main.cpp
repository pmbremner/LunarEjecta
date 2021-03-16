#include "LunarEjecta_params.h"
#include "LunarEjecta_igloo.h"

using namespace std;

int main(int argc, char const *argv[])
{

	string param_fn = argv[1];
	int N_proc = atoi(argv[2]); // total number of processes
	int i_proc = atoi(argv[3]); // the current process

	cout << "Process number = " << i_proc << " / " << N_proc << endl << endl;
	if (i_proc >= N_proc || i_proc < 0)
	{
		cout << "ERROR: i_proc must be between 0 and " << N_proc - 1 << endl;
		return -1;
	}
	
	// input contains:
	//  parameters read from file, including regolith info, resolution of secondary ejecta igloo file
	input *params = init_input(param_fn, N_proc, i_proc);


	
	// HiDensMEM, LoDensMEM, NEO

	iglooSet *MEM_hi_fluxes = read_igloo(params, HiDensMEM);
	iglooSet *MEM_lo_fluxes = read_igloo(params, LoDensMEM);
	iglooSet *NEO_fluxes;

	iglooSet *primaryFluxes[3] = {MEM_hi_fluxes, MEM_lo_fluxes, NEO_fluxes};



	if (params->readNEO_files)
		NEO_fluxes = read_igloo(params, NEO);
	else {
		NEO_fluxes = generate_NEO_igloo(params, MEM_hi_fluxes);
		if (params->saveNEO_files)
			save_igloo(params, NEO_fluxes, NEO);
	}

	/*

	// note: for each lat-lon location, there will be a map for each zenith angle bin
	hist2DSet *azmDistMap = init_azm_dist_map(params);
	
	// includes surface area of lat-lon locations, as well as scaling law stuff
	scalingLaw *ejectaFactors = compute_constants_and_normalization(params);

	iglooSet *secondary_fluxes = compute_ejecta(params, primaryFluxes, azmDistMap, ejectaFactors);
	save_igloo(params, secondary_fluxes, "secondaryEjecta");

	*/

	return 0;
}