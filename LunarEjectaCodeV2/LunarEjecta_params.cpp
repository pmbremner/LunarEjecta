#include "LunarEjecta_params.h"


#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;


input* init_input(string param_fn, int N_proc, int i_proc)
{
	input* p = new input;

	p->N_proc = N_proc;
	p->i_proc = i_proc;

	//////////////////////
	// Default parameters
	p->MEM_massMin = 1.E-6; // g
	p->MEM_massMax = 10.; // g
	p->NEO_massMin = 10.; // g
	p->NEO_massMax = 1.5708E15; // g, 1000 m diameter at 3 g/cc

	p->lunar_radius = 1737.4; // km
	p->lunar_escape_speed = 2.38; // km/s
	//////////////////////

	//////////////////////
	// Required parameters

	getParam(param_fn, "Nlat", p->Nlat, 0);
	getParam(param_fn, "Nlon", p->Nlon, 0);

	p->Nlatlon_tot = p->Nlat * p->Nlon;
	p->N_loc       = ceil(p->Nlatlon_tot / double(p->N_proc));

	p->latlon_idx_min = p->i_proc * p->N_loc;
	p->latlon_idx_max = fmin((p->i_proc + 1) * p->N_loc - 1, p->Nlatlon_tot - 1);

	p->N_loc = p->latlon_idx_max - p->latlon_idx_min + 1; // fix number of locations (mostly for the last process)

	cout << "min and max lat-lon index (absolute) = " << p->latlon_idx_min << ", " << p->latlon_idx_max << " | for proc = " << p->N_loc << " | total = " << p->Nlatlon_tot << endl;

	p->lon_directory.resize(4);
	for (int i = 0; i < p->Nlon; ++i)
		getParam(param_fn, "MEMDirectoryName_" + to_string(i * 360/p->Nlon), p->lon_directory[i], 0);

	getParam(param_fn, "readNEO_files", p->readNEO_files, 0);
	getParam(param_fn, "saveNEO_files", p->saveNEO_files, 0);

	getParam(param_fn, "ROI_radius", p->ROI_radius, 0);

	//////////////////////
	return p;
}