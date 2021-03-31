#include "LunarEjecta_secondaryejecta.h"
#include "LunarEjecta_params.h"
#include "LunarEjecta_igloo.h"
#include "LunarEjecta_scalinglaw.h"
#include "LunarEjecta_bearingdistmap.h"

#include <cmath>

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
	int NrowVars = 9, i, j, idx, count = 0;

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
	secEjectaFluxes->N_cols = p->N_vel; // NrowVars not counted here
	cout << "igloo N_rows/N_cols: " << secEjectaFluxes->N_rows << ' ' << secEjectaFluxes->N_cols << endl;

	// setup speed bin edges and centers
	linspace(secEjectaFluxes->speedEdge, p->vel_min, p->lunar_escape_speed, p->N_vel+1);
	double dv = secEjectaFluxes->speedEdge[1] - secEjectaFluxes->speedEdge[0];
	linspace(secEjectaFluxes->speedCenter, p->vel_min + dv/2., p->lunar_escape_speed - dv/2., p->N_vel);

	// setup info columns
	double dhor_ang;
	secEjectaFluxes->iglooData.resize(0);
	secEjectaFluxes->iglooData.resize(secEjectaFluxes->N_rows * (secEjectaFluxes->N_cols + NrowVars), 0.);

	for (i = 0; i < bearingDistMap->N_horizon_ROI; ++i)
	{
		for (j = 0; j < bearingDistMap->N_bearing_ROI[bearingDistMap->N_horizon_ROI-i-1]; ++j)
		{
			idx = count * (secEjectaFluxes->N_cols + NrowVars); 

			secEjectaFluxes->iglooData[idx    ] = count + 1;
			secEjectaFluxes->iglooData[idx + 1] = i + 1;
			secEjectaFluxes->iglooData[idx + 2] = j + 1;
			secEjectaFluxes->iglooData[idx + 3] = i * 90. / double(bearingDistMap->N_horizon_ROI); // left bound, horizon angle, degrees
			secEjectaFluxes->iglooData[idx + 4] = (i + 1.) * 90. / double(bearingDistMap->N_horizon_ROI); // right bound, horizon angle, degrees
			secEjectaFluxes->iglooData[idx + 5] = j * 360. / double(bearingDistMap->N_bearing_ROI[bearingDistMap->N_horizon_ROI-i-1]); // left bound, azimuth angle, degrees
			secEjectaFluxes->iglooData[idx + 6] = (j + 1.) * 360. / double(bearingDistMap->N_bearing_ROI[bearingDistMap->N_horizon_ROI-i-1]); // right bound, azimuth angle, degrees
			secEjectaFluxes->iglooData[idx + 7] = (180./PI) * asin((sin(secEjectaFluxes->iglooData[idx + 3]*PI/180.) + sin(secEjectaFluxes->iglooData[idx + 4]*PI/180.))/2.); // average horizon angle in bin
			secEjectaFluxes->iglooData[idx + 8] = (secEjectaFluxes->iglooData[idx + 5] + secEjectaFluxes->iglooData[idx + 6]) / 2.; // average azimuth angle in bin

			count++;
		}

	}


	// for each ejecta zenith angle
	// for each ejecta speed

	// for each distance -> then compute eval points, if none, skip all following loops
	// for each POI bearing (outgoing ejecta azimuth angle)
	// for each ROI bearing (incoming ejecta azimuth angle)

	// for each impact horizon angle
	// for each impact azimuth angle
	// for each impact speed


	return secEjectaFluxes;
}