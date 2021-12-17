#include "LunarEjecta_params.h"
#include "LunarEjecta_track.h"
#include "LunarEjecta_asset.h"

using namespace std;

void init_trackVars(trackVars& tv, double x0, double y0, double z0, double u0, double v0, double w0, double h, double dx0)
{
	tv.x0 = x0;
	tv.y0 = y0;
	tv.z0 = z0;
	tv.u0 = u0;
	tv.v0 = v0;
	tv.w0 = w0;
	tv.xi = x0;
	tv.yi = y0;
	tv.zi = z0;
	tv.ui = u0;
	tv.vi = v0;
	tv.wi = w0;

	tv.h = h;
	tv.dxmin = dx0;

	// if(VERBOSE_INIT > 1){
		// cout << "Initializing trackVars paramerters:\n";
		// cout << "  x0 = " << x0 << " m\n";
		// cout << "  y0 = " << y0 << " m\n";
		// cout << "  z0 = " << z0 << " m\n";
		// cout << "  u0 = " << u0 << " m/s\n";
		// cout << "  v0 = " << v0 << " m/s\n";
		// cout << "  w0 = " << w0 << " m/s\n";
		// cout << "  h  = " << h << " s\n\n"; 
	// }
}

void grav_a(double& ax, double& ay, double& az, double x, double y, double z, double Rm, double gm)
{
	double grav_coeff = -gm * sqr(Rm) * pow(mag2(x, y, z), -1.5);

	ax = grav_coeff * x;
	ay = grav_coeff * y;
	az = grav_coeff * z;
}

// see https://math.okstate.edu/people/yqwang/teaching/math4513_fall11/Notes/rungekutta.pdf
void RK45UpdatePosVel(trackVars& t,
	                  RK45VarsPosVel& r,
	                  const double c[][6], asset& ag,
	                  double Rm,
	                  double gm)
{
	int i, j;
	double ax, ay, az, rcur, xcur, ycur, zcur, vcur, ucur, wcur, maxR;
	double x_a, y_a, z_a;

	double t_dot_c, c_x0, c_z0, ht, st, dt_temp, rt;

	do {
		for (i = 0; i < 6; ++i) // RK45 steps 1-6
		{
			r.kx[i] = t.ui;
			r.ky[i] = t.vi;
			r.kz[i] = t.wi;
			xcur = t.xi;
			ycur = t.yi;
			zcur = t.zi;
			for (j = 0; j < i; ++j)
			{
				r.kx[i] += c[i][j] * r.ku[j];
				r.ky[i] += c[i][j] * r.kv[j];
				r.kz[i] += c[i][j] * r.kw[j];
				xcur    += c[i][j] * r.kx[j];
				ycur    += c[i][j] * r.ky[j];
				zcur    += c[i][j] * r.kz[j];
			}
			r.kx[i] *= t.h;
			r.ky[i] *= t.h;
			r.kz[i] *= t.h;

			grav_a(ax, ay, az, xcur, ycur, zcur, Rm, gm);

			r.ku[i] = t.h * ax;
			r.kv[i] = t.h * ay;
			r.kw[i] = t.h * az;
		}

		// check error and update timestep, t.h
		for (i = 0; i < 2; ++i)
		{
			r.xi1[i] = t.xi;
			r.yi1[i] = t.yi;
			r.zi1[i] = t.zi;
			r.ui1[i] = t.ui;
			r.vi1[i] = t.vi;
			r.wi1[i] = t.wi;
			for (j = 0; j < 6; ++j)
			{
				r.xi1[i] += c[i+6][j] * r.kx[j];
				r.yi1[i] += c[i+6][j] * r.ky[j];
				r.zi1[i] += c[i+6][j] * r.kz[j];
				r.ui1[i] += c[i+6][j] * r.ku[j];
				r.vi1[i] += c[i+6][j] * r.kv[j];
				r.wi1[i] += c[i+6][j] * r.kw[j];
			}
		}

		r.Rx = fabs(r.xi1[1] - r.xi1[0]) / t.h;
		r.Ry = fabs(r.yi1[1] - r.yi1[0]) / t.h;
		r.Rz = fabs(r.zi1[1] - r.zi1[0]) / t.h;
		// r.Ru = fabs(r.ui1[1] - r.ui1[0]) / t.h;
		// r.Rv = fabs(r.vi1[1] - r.vi1[0]) / t.h;
		// r.Rw = fabs(r.wi1[1] - r.wi1[0]) / t.h;

		maxR = max(max(r.Rx, r.Ry), r.Rz);
		t.h *= (maxR > 0.0 ? 0.84 * pow(EPS / maxR, 0.25) : 1.2);


		// shift position to asset frame (no rotation)
		x_a = r.xi1[0] - ag.origin.x[0];
		y_a = r.yi1[0] - ag.origin.x[1];
		z_a = r.zi1[0] - ag.origin.x[2];

		// distance from asset origin divided by speed
		dt_temp = mag_s(x_a, y_a, z_a) / mag_s(r.ui1[0], r.vi1[0], r.wi1[0]); 

		// if position is within collision boundary, reduce timestep to 2% of the radius
		// otherize, reduce timestep by 50%
		dt_temp /= (mag_s(x_a, y_a, z_a) <= ag.collision_radius_boundary ? 50. : 2.);

		// replace timestep to not overstep the asset
		t.h = (dt_temp < t.h ? dt_temp : t.h);


		////////// OLD CODE ////////////////
		// Check if within boundary sphere

		// If outside sphere limit dt to slighlty half of distance to sphere


		// If inside sphere, limit dt to a few percent of the boundary sphere size



		// c_x0 = sin(ll.distance/Rm);
		// c_z0 = cos(ll.distance/Rm);

		// //t_dot_c = t.xi * c_x0 + t.zi * c_z0;

		// // height of traj wrt cylinder
		// //ht = t_dot_c - Rm;

		// // cylindrical radius wrt cylinder
		// //st = mag_s(t.xi - t_dot_c * c_x0, 0., t.zi - t_dot_c * c_z0);

		// rt = (mag_s(t.xi - c_x0*Rm, 0.0, t.zi - c_z0*Rm) - 0.9*ll.radius)/10.;

		// dt_temp = rt / mag_s(t.ui, t.vi, t.wi);

		// t.h = (dt_temp < t.h ? dt_temp : t.h);

		// if(DEBUG > 1){
		// 	cout << "Rx = " << r.Rx << endl;
		// 	cout << "Ry = " << r.Ry << endl;
		// 	cout << "Rz = " << r.Rz << endl;
		// 	cout << "Ru = " << r.Ru << endl;
		// 	cout << "Rv = " << r.Rv << endl;
		// 	cout << "Rw = " << r.Rw << endl;
		// 	cout << "timestep = " << t.h << endl;
		// 	cout << "maxR = " << maxR << ' ' << EPS << endl << endl;
		// }
		///// END OF OLD CODE /////////

	} while(maxR > EPS);

	// update position and velocity
	t.xi = r.xi1[0];
	t.yi = r.yi1[0];
	t.zi = r.zi1[0];
	t.ui = r.ui1[0];
	t.vi = r.vi1[0];
	t.wi = r.wi1[0];

	// if(DEBUG > 1) {
	// 	cout << "xi = " << t.xi << endl;
	// 	cout << "yi = " << t.yi << endl;
	// 	cout << "zi = " << t.zi << endl;
	// 	cout << "ui = " << t.ui << endl;
	// 	cout << "vi = " << t.vi << endl;
	// 	cout << "wi = " << t.wi << endl << endl;
	// }
}



// void printTrack(ofstream& file, trackVars& tv, sums& s)
// {
// 	double speed = mag_s(tv.u0, 0.0, tv.w0) / vesc; // units of escape speed
// 	double zang = atan2(tv.u0, tv.w0);
// 	double dist = (zang > 0. ? atan2(tv.xi, tv.zi) : 2.*PI - atan2(tv.xi, tv.zi)); // units of lunar radius
// 	file << scientific << setprecision(14) << tv.xi  << ' ' << tv.zi << ' ' << dist << ' ' << speed << ' ' << zang
// 		 << ' ' << s.sum_miss << ' ' << s.sum_hit << ' ' << s.N_miss << ' ' << s.N_hit << endl;
// }

bool checkCollision(trackVars& t, asset& ag)
{
	// double t_dot_c;

	// double c_x0 = sin(c.distance/Rm);
	// double c_z0 = cos(c.distance/Rm);

	// t_dot_c = t.xi * c_x0 + t.zi * c_z0;

	// // height of traj wrt cylinder
	// double ht = t_dot_c - Rm;

	// // cylindrical radius wrt cylinder
	// double st = mag_s(t.xi - t_dot_c * c_x0, 0., t.zi - t_dot_c * c_z0);

	// //cout << "ht = " << ht << " | st = " << st << endl;

	// if (ht > 0. && ht <= c.height && st <= c.radius)
	// 	return 1;
	return 0;
}

void unpackFinalLoc(trackVars& track_i, vector<double> &loc_f)
{
	loc_f.clear();
	loc_f.push_back(track_i.xi);
	loc_f.push_back(track_i.yi);
	loc_f.push_back(track_i.zi);
}


void unpackFinalPh(trackVars& track_i, vector<double> &ph_f)
{
	ph_f.clear();
	ph_f.push_back(track_i.ui);
	ph_f.push_back(track_i.vi);
	ph_f.push_back(track_i.wi);
}


// loc = physical location in global cartesian coordinates
// ph  = location of initial ejecta conditions in phase space (i.e., ejecta speed, zenith angle, and bearing angle)
// 
bool runTraj_checkHit(vector<double> &lat_lon, vector<double> &loc, vector<double> &ph, vector<double> &loc_f, vector<double> &ph_f, asset &ag, double Rm, double vesc, double gm)
{
	// init track and RK45 vars
	RK45VarsPosVel RK45Vars;
	trackVars track_i;
	bool moon_hit_flag = 0, asset_hit_flag = 0, escape_flag = 0;

	// define velocity vector and init
	// vel_before_rot = velocity as if at the north pole
	// rot_m = rotation matrix to move velocity vectory from pole to lat-lon location
	// vel_after_rot = velocity in the local frame at the lat-lon location
	mat3x3 rot_m;
	vec3 vel_before_rot, vel_after_rot;
	vel_before_rot.x[0] = -ph[0] * cos(ph[2]) * sin(ph[1]); // bearing angle introduces a minus sign here
	vel_before_rot.x[1] = ph[0] * sin(ph[2]) * sin(ph[1]);
	vel_before_rot.x[2] = ph[0] * cos(ph[1]);

	//h_rot_m_from_angs( rot_m, atan2(sqrt(sqr(loc[0]) + sqr(loc[1])), loc[2]), atan2(loc[1], loc[0]) );
	h_rot_m_from_angs( rot_m, PI/2. - lat_lon[0], lat_lon[1] );

	h_matrix_vector_multiply(vel_after_rot, rot_m, vel_before_rot);


	//cout << " vel mag before and after rot = " << mag_s(vel_before_rot) << " | " << mag_s(vel_after_rot) << endl;

	init_trackVars(track_i,
	/* x0 */       loc[0], // m
	/* y0 */       loc[1], // m
	/* z0 */       loc[2], // m
	/* u0 */       vel_after_rot.x[0], // m/s
	/* v0 */       vel_after_rot.x[1], // m/s
	/* w0 */       vel_after_rot.x[2], // m/s
	/* h  */       1e-4, // initial dt, s
	/* dxmin */    -1.); // not used at the moment

	//cout << endl << endl;
	/// step through trajectory
	while (!(moon_hit_flag || escape_flag || asset_hit_flag))
	{
		// update position and velocity
		RK45UpdatePosVel(track_i, RK45Vars, RK45Coeff, ag, Rm, gm);
		// cout.precision(12);
		// cout <<  track_i.xi/Rm << ' ' << track_i.yi/Rm << ' ' << track_i.zi/Rm << ' ' << mag_s(track_i.xi/Rm, track_i.yi/Rm, track_i.zi/Rm) << ' ' << track_i.h << endl;

		// check for collision with moon
		moon_hit_flag  = check_collision_moon(Rm, track_i.xi, track_i.yi, track_i.zi);

		if (moon_hit_flag){
			unpackFinalLoc(track_i, loc_f);
			unpackFinalPh(track_i, ph_f);
			return 0; // missed asset
		}

		// if speed greater than escape speed and position is beyond asset origin + sphere boundary
		escape_flag    = check_escape(ag, vesc, track_i.xi, track_i.yi, track_i.zi, track_i.ui, track_i.vi, track_i.wi);

		if (escape_flag){
			unpackFinalLoc(track_i, loc_f);
			unpackFinalPh(track_i, ph_f);
			return 0; // missed asset, won't return to hit moon
		}

		// if pos is outside boundary radius, don't check asset, otherwise check asset
		asset_hit_flag = check_collision_asset(ag, track_i.xi, track_i.yi, track_i.zi);

		if (asset_hit_flag){
			unpackFinalLoc(track_i, loc_f);
			unpackFinalPh(track_i, ph_f);
			return 1; // hit the asset
		}
	}
	

	return 0;
}