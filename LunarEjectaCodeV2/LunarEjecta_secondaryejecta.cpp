#include "LunarEjecta_secondaryejecta.h"
#include "LunarEjecta_params.h"
#include "LunarEjecta_igloo.h"
#include "LunarEjecta_scalinglaw.h"
#include "LunarEjecta_bearingdistmap.h"

using namespace std;



// struct iglooSet // for each location
// {
// 	string filename; // filename of igloo file
// 	double lat; // radians
// 	double lon; // radians [0, 2*pi)
// 	double SA;  // surface area of location on Moon

// 	int N_rows; // number of rows in igloo data
// 	int N_cols; // number of columns in igloo data
// 	vector<double> iglooData;   // includes info columns (extra 9 columns) and data
// 	vector<double> speedEdge;   // size = N_cols - 8
// 	vector<double> speedCenter; // size = N_cols - 9
// };


// called per lat-lon location
iglooSet* compute_ejecta(input *p, vector<iglooSet*> &primaryFluxes, scalingLaw *ejectaFactors, hist3DSet* bearingDistMap, int ejectaDistType)
{
	iglooSet* secEjectaFluxes = new iglooSet;

	cout << "compute_ejecta idx = " << p->latlon_idx_proc << endl;
	//cout << primaryFluxes[HiDensMEM][p->latlon_idx_proc].filename << endl;
	string ejectaDistType_name;
	if (ejectaDistType == ejectaShort)
		ejectaDistType_name = "Short";
	else if(ejectaDistType == ejectaFar)
		ejectaDistType_name = "Far";

	secEjectaFluxes->filename = primaryFluxes[HiDensMEM][p->latlon_idx_proc].filename.substr(0, primaryFluxes[HiDensMEM][p->latlon_idx_proc].filename.length() - 23) + "SecEjecta" + ejectaDistType_name + "_igloo_avg.txt";
	cout << secEjectaFluxes->filename << endl;

	double dlat = PI / double(p->Nlat - 1.); // rad
	double dlon = 2.*PI / double(p->Nlon);   // rad

	secEjectaFluxes->lat = -PI/2. + dlat * (p->latlon_idx_min + p->latlon_idx_proc % p->Nlat);
	secEjectaFluxes->lon = dlon * (int(p->latlon_idx_min + p->latlon_idx_proc) / int(p->Nlat));

	secEjectaFluxes->N_rows = bearingDistMap->N_bearing_ROI_tot;
	secEjectaFluxes->N_cols = p->N_vel + 9;
	cout << "igloo N_rows/N_cols: " << secEjectaFluxes->N_rows << ' ' << secEjectaFluxes->N_cols << endl;


	return secEjectaFluxes;
}