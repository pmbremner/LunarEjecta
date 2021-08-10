#include "LunarEjecta_params.h"
#include "LunarEjecta_track.h"

using namespace std;

void grav_a(double& ax, double& ay, double& az, double x, double y, double z)
{
	double grav_coeff = -gm * sqr(Rm) * pow(mag2(x, y, z), -1.5);

	ax = grav_coeff * x;
	ay = grav_coeff * y;
	az = grav_coeff * z;
}

// see https://math.okstate.edu/people/yqwang/teaching/math4513_fall11/Notes/rungekutta.pdf
void RK45UpdatePosVel(trackVars& t,
	                  RK45VarsPosVel& r,
	                  const double c[][6], asset_geometry& ag)
{
	int i, j;
	double ax, ay, az, rcur, xcur, ycur, zcur, vcur, ucur, wcur, maxR;

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

			grav_a(ax, ay, az, xcur, ycur, zcur);

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
		r.Ru = fabs(r.ui1[1] - r.ui1[0]) / t.h;
		r.Rv = fabs(r.vi1[1] - r.vi1[0]) / t.h;
		r.Rw = fabs(r.wi1[1] - r.wi1[0]) / t.h;

		maxR = max(max(r.Rx, r.Ry), r.Rz);
		t.h *= (maxR > 0.0 ? 0.84 * pow(EPS / maxR, 0.25) : 1.2);


		c_x0 = sin(ll.distance/Rm);
		c_z0 = cos(ll.distance/Rm);

		//t_dot_c = t.xi * c_x0 + t.zi * c_z0;

		// height of traj wrt cylinder
		//ht = t_dot_c - Rm;

		// cylindrical radius wrt cylinder
		//st = mag_s(t.xi - t_dot_c * c_x0, 0., t.zi - t_dot_c * c_z0);

		rt = (mag_s(t.xi - c_x0*Rm, 0.0, t.zi - c_z0*Rm) - 0.9*ll.radius)/10.;

		dt_temp = rt / mag_s(t.ui, t.vi, t.wi);

		t.h = (dt_temp < t.h ? dt_temp : t.h);

		if(DEBUG > 1){
			cout << "Rx = " << r.Rx << endl;
			cout << "Ry = " << r.Ry << endl;
			cout << "Rz = " << r.Rz << endl;
			cout << "Ru = " << r.Ru << endl;
			cout << "Rv = " << r.Rv << endl;
			cout << "Rw = " << r.Rw << endl;
			cout << "timestep = " << t.h << endl;
			cout << "maxR = " << maxR << ' ' << EPS << endl << endl;
		}

	} while(maxR > EPS);

	// update position and velocity
	t.xi = r.xi1[0];
	t.yi = r.yi1[0];
	t.zi = r.zi1[0];
	t.ui = r.ui1[0];
	t.vi = r.vi1[0];
	t.wi = r.wi1[0];

	if(DEBUG > 1) {
		cout << "xi = " << t.xi << endl;
		cout << "yi = " << t.yi << endl;
		cout << "zi = " << t.zi << endl;
		cout << "ui = " << t.ui << endl;
		cout << "vi = " << t.vi << endl;
		cout << "wi = " << t.wi << endl << endl;
	}
}



void printTrack(ofstream& file, trackVars& tv, sums& s)
{
	double speed = mag_s(tv.u0, 0.0, tv.w0) / vesc; // units of escape speed
	double zang = atan2(tv.u0, tv.w0);
	double dist = (zang > 0. ? atan2(tv.xi, tv.zi) : 2.*PI - atan2(tv.xi, tv.zi)); // units of lunar radius
	file << scientific << setprecision(14) << tv.xi  << ' ' << tv.zi << ' ' << dist << ' ' << speed << ' ' << zang
		 << ' ' << s.sum_miss << ' ' << s.sum_hit << ' ' << s.N_miss << ' ' << s.N_hit << endl;
}

bool checkCollision(trackVars& t, asset_geometry& ag)
{
	double t_dot_c;

	double c_x0 = sin(c.distance/Rm);
	double c_z0 = cos(c.distance/Rm);

	t_dot_c = t.xi * c_x0 + t.zi * c_z0;

	// height of traj wrt cylinder
	double ht = t_dot_c - Rm;

	// cylindrical radius wrt cylinder
	double st = mag_s(t.xi - t_dot_c * c_x0, 0., t.zi - t_dot_c * c_z0);

	//cout << "ht = " << ht << " | st = " << st << endl;

	if (ht > 0. && ht <= c.height && st <= c.radius)
		return 1;
	return 0;
}