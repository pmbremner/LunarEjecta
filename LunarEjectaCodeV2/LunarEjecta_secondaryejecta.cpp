#include "LunarEjecta_secondaryejecta.h"
#include "LunarEjecta_params.h"
#include "LunarEjecta_igloo.h"
#include "LunarEjecta_scalinglaw.h"
#include "LunarEjecta_bearingdistmap.h"
#include "LunarEjecta_speedzenith.h"

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
			secEjectaFluxes->iglooData[idx + 5] = fmod(j * 360. / double(bearingDistMap->N_bearing_ROI[bearingDistMap->N_horizon_ROI-i-1]) + 90., 360.); // left bound, azimuth angle, degrees
			secEjectaFluxes->iglooData[idx + 6] = fmod((j + 1.) * 360. / double(bearingDistMap->N_bearing_ROI[bearingDistMap->N_horizon_ROI-i-1]) + 90., 360.); // right bound, azimuth angle, degrees
			secEjectaFluxes->iglooData[idx + 7] = (180./PI) * asin((sin(secEjectaFluxes->iglooData[idx + 3]*PI/180.) + sin(secEjectaFluxes->iglooData[idx + 4]*PI/180.))/2.); // average horizon angle in bin
			secEjectaFluxes->iglooData[idx + 8] = (secEjectaFluxes->iglooData[idx + 5] + secEjectaFluxes->iglooData[idx + 6]) / 2.; // average azimuth angle in bin

			count++;
		}

	}

	
	int i0_ej_ang, j0_ej_v, i1_d, j1_POI, k1_ROI, bin_idx, i2_imp_row, j2_imp_v, k2_imp_azm, idx_row, i2_imp_row_start, idx_offset, h_ang_idx, N_imp_azm;
	int i_flux;
	double D_min_ej, D_max_ej, D_min_map, D_max_map;
	double z_ang_min_ej, z_ang_max_ej, v_min_ej, v_max_ej;
	double azm_POI_min, azm_POI_max, azm_ROI_min, azm_ROI_max; 
	double h_ang_min_imp, h_ang_max_imp, h_ang_avg_imp, v_avg_imp, azm_min_imp, azm_max_imp, d_h_ang;
	vector<double> cur_flux(3, 0.0);

	// loop through ejecta zenith angle and ejecta speed
	for (i0_ej_ang = 0; i0_ej_ang < bearingDistMap->N_horizon_ROI; ++i0_ej_ang)
	{
		z_ang_min_ej = PI/2. * i0_ej_ang        / double(bearingDistMap->N_horizon_ROI); // rad
		z_ang_max_ej = PI/2. * (i0_ej_ang + 1.) / double(bearingDistMap->N_horizon_ROI); // rad

		for (j0_ej_v = 0; j0_ej_v < p->N_vel; ++j0_ej_v)
		{
			v_min_ej = secEjectaFluxes->speedEdge[j0_ej_v] / p->lunar_escape_speed; // units of escape speed
			v_max_ej = secEjectaFluxes->speedEdge[j0_ej_v + 1] / p->lunar_escape_speed; // units of escape speed


			// compute what range of distance the ejecta can go to
			vector<double> D;
			D.push_back(getD(z_ang_min_ej, v_min_ej)); // needs units of rads, radii
			D.push_back(getD(z_ang_min_ej, v_max_ej));
			D.push_back(getD(z_ang_max_ej, v_min_ej));
			D.push_back(getD(z_ang_max_ej, v_max_ej));

			D_min_ej = vMin(D); // units of radii
			D_max_ej = vMax(D);


			// Loop through bearing distance map
			for (i1_d = 0; i1_d < bearingDistMap->N_D; ++i1_d)
			{
				// units of radii
				D_min_map = bearingDistMap->D_min + (bearingDistMap->D_max - bearingDistMap->D_min) * i1_d        / double(bearingDistMap->N_D);
				D_max_map = bearingDistMap->D_min + (bearingDistMap->D_max - bearingDistMap->D_min) * (i1_d + 1.) / double(bearingDistMap->N_D);

				// if the valid ejecta distance range is within/overlapping the specific range in the distance map
				// if(!(D_max_ej <= D_min_map || D_min_ej >= D_max_map))
				// {
					// compute points to evaluate
					vector<quad> q;
					getEvalQuads(q, D_min_map, D_max_map, z_ang_min_ej, z_ang_max_ej, v_min_ej, v_max_ej, 10, 0);

					cout << q.size() << ' ';

					// if the # of points is greater than zero
				 	// evaluate the points here if no directional azimuth
					if (q.size() > 0)
					{
						for (j1_POI = 0; j1_POI < bearingDistMap->N_bearing_POI; ++j1_POI)
						{
							azm_POI_min = bearingDistMap->bearing_POI_min + (bearingDistMap->bearing_POI_max - bearingDistMap->bearing_POI_min) * j1_POI        / double(bearingDistMap->N_bearing_POI);
							azm_POI_max = bearingDistMap->bearing_POI_min + (bearingDistMap->bearing_POI_max - bearingDistMap->bearing_POI_min) * (j1_POI + 1.) / double(bearingDistMap->N_bearing_POI);

							// convert from bearing to azimuth
							azm_POI_min = fmod(azm_POI_min + PI/2., 2.*PI);
							azm_POI_max = fmod(azm_POI_max + PI/2., 2.*PI);

							for (k1_ROI = 0; k1_ROI < bearingDistMap->N_bearing_ROI[i0_ej_ang]; ++k1_ROI)
							{
								azm_ROI_min = 2.*PI * k1_ROI        / double(bearingDistMap->N_bearing_ROI[i0_ej_ang]);
								azm_ROI_max = 2.*PI * (k1_ROI + 1.) / double(bearingDistMap->N_bearing_ROI[i0_ej_ang]);

								// convert from bearing to azimuth
								azm_ROI_min = fmod(azm_ROI_min + PI/2., 2.*PI);
								azm_ROI_max = fmod(azm_ROI_max + PI/2., 2.*PI);


								bin_idx = idx_bin(bearingDistMap, i0_ej_ang, i1_d, j1_POI, k1_ROI);


								// for (i_flux = 0; i_flux < primaryFluxes.size(); ++i_flux)
								// {
								// 	cur_flux[i_flux] = compute_ejecta_at(z_ang_min_ej, z_ang_max_ej, v_min_ej, v_max_ej)
								// }

								// loop through impactor vars
								// for each row in the igloo data
								// i2_imp_row = 0;
								// // error down here
								// while (i2_imp_row < primaryFluxes[NEO][p->latlon_idx_proc].N_rows)
								// {
								// 	idx_row = i2_imp_row * (primaryFluxes[NEO][p->latlon_idx_proc].N_cols + NrowVars);

								// 	h_ang_min_imp = primaryFluxes[NEO][p->latlon_idx_proc].iglooData[idx_row + 3] * PI/180.; // rad
								// 	h_ang_max_imp = primaryFluxes[NEO][p->latlon_idx_proc].iglooData[idx_row + 4] * PI/180.; // rad
								// 	h_ang_avg_imp = primaryFluxes[NEO][p->latlon_idx_proc].iglooData[idx_row + 7] * PI/180.; // rad

								// 	h_ang_idx = primaryFluxes[NEO][p->latlon_idx_proc].iglooData[idx_row + 1] - 1;
								// 	d_h_ang   = (h_ang_max_imp-h_ang_min_imp);
								// 	N_imp_azm = round(fabs(d_h_ang * ( cos(h_ang_idx*d_h_ang) + cos((h_ang_idx+1.)*d_h_ang) )/2.));

								// 	cout << i2_imp_row << ' ' << h_ang_min_imp << ' ' << h_ang_max_imp << ' ' << h_ang_idx << endl;

								// 	i2_imp_row += N_imp_azm;



								// 	//i2_imp_row_start = i2_imp_row;

								// 	// for each impact speed
								// 	// for (j2_imp_v = 0; j2_imp_v < primaryFluxes[NEO][p->latlon_idx_proc].N_cols; ++j2_imp_v)
								// 	// {
								// 	// 	v_avg_imp = primaryFluxes[NEO][p->latlon_idx_proc].speedCenter[j2_imp_v];

								// 	// 	// for each impact azimuth, just sum the fluxes
								// 	// 	cur_flux[0] = cur_flux[1] = cur_flux[2] = 0.0;
								// 	// 	i2_imp_row = i2_imp_row_start;

								// 	// 	for (k2_imp_azm = 0; k2_imp_azm < N_imp_azm; ++k2_imp_azm)
								// 	// 	{//HiDensMEM, LoDensMEM, NEO

								// 	// 		idx_row = i2_imp_row * (primaryFluxes[NEO][p->latlon_idx_proc].N_cols + NrowVars);
								// 	// 		idx_offset = primaryFluxes[NEO][p->latlon_idx_proc].N_rows * (primaryFluxes[NEO][p->latlon_idx_proc].N_cols + NrowVars); // used to skip below horizon

								// 	// 		cur_flux[HiDensMEM] += primaryFluxes[HiDensMEM][p->latlon_idx_proc].iglooData[idx_row + NrowVars + j2_imp_v + idx_offset];
								// 	// 		cur_flux[LoDensMEM] += primaryFluxes[LoDensMEM][p->latlon_idx_proc].iglooData[idx_row + NrowVars + j2_imp_v + idx_offset];
								// 	// 		cur_flux[NEO      ] += primaryFluxes[NEO      ][p->latlon_idx_proc].iglooData[idx_row + NrowVars + j2_imp_v];

								// 	// 		i2_imp_row++;
								// 	// 	}							

								// 	// 	// for each impact mass
								// 	// 	// for each density

								// 	// 	// compute ejecta
								// 	// }
								// 	//cout << i2_imp_row << " / " << primaryFluxes[NEO][p->latlon_idx_proc].N_rows << endl;
								// }

								// sec_ej_flux[???] = bin[bin_idx](fraction of flux) * (normalized integral within the bounds)


							}
						}

					//} // END IF in distance range
				}
			}
			cout << endl;


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