#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <time.h>
#include <algorithm> // max(a,b)
//#include <mpi.h>
#include <stdlib.h> // rand(), atof()

using namespace std;

#define DEBUG 0 // 0 = no output, 1 = per track, 2 = per time step
#define VERBOSE_INIT 0 // 0 = no output, 1 = all but tracks, 2 = all
#define HELP_INIT 0

const double PI = 3.14159;
const double gm = 1.625;    // m/s^2 at Moon's surface
const double Rm = 1737.1E3; // m, lunar radius
const double vesc = 2.38E3; // m/s, Moon escape speed
const double EPS = 1E-4;
const int NMAX = 1E4;
const double mu_HH11 = 0.4;

inline double sqr(double a) {return a*a;}

inline double mag_s(double x, double y, double z) {return sqrt(sqr(x) + sqr(y) + sqr(z));}
inline double mag2(double x, double y, double z) {return sqr(x) + sqr(y) + sqr(z);}

struct trackVars
{
	double x0;
	double y0;
	double z0;
	double u0;
	double v0;
	double w0;

	double xi;
	double yi;
	double zi;
	double ui;
	double vi;
	double wi;

	// double dpx;
	// double dpy;
	// double Fx;
	// double Fy;

	double h; // timestep
	double dxmin; // minimum distance
};

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
	// tv.dpx = 0.0;
	// tv.dpy = 0.0;
	// tv.Fx = 0.0;
	// tv.Fy = 0.0;
	tv.h = h;
	tv.dxmin = dx0;

	if(VERBOSE_INIT > 1){
		cout << "Initializing trackVars paramerters:\n";
		cout << "  x0 = " << x0 << " m\n";
		cout << "  y0 = " << y0 << " m\n";
		cout << "  z0 = " << z0 << " m\n";
		cout << "  u0 = " << u0 << " m/s\n";
		cout << "  v0 = " << v0 << " m/s\n";
		cout << "  w0 = " << w0 << " m/s\n";
		cout << "  h  = " << h << " s\n\n"; 
	}
}

// void print_track(trackVars& tv)
// {
// 	cout << "Track info:\n";
// 	cout << "  x0 = " << tv.x0 << " m\n";
// 	cout << "  y0 = " << tv.y0 << " m\n";
// 	cout << "  u0 = " << tv.u0 << " m/s\n";
// 	cout << "  v0 = " << tv.v0 << " m/s\n";
// 	cout << "  xi = " << tv.xi << " m\n";
// 	cout << "  yi = " << tv.yi << " m\n";
// 	cout << "  ui = " << tv.ui << " m/s\n";
// 	cout << "  vi = " << tv.vi << " m/s\n";
// 	// cout << "  dpx = " << tv.dpx << " kg*m\n";
// 	// cout << "  dpy = " << tv.dpy << " kg*m\n";
// 	// cout << "  Fx = " << tv.Fx << " N/m\n";
// 	// cout << "  Fy = " << tv.Fy << " N/m\n";
// 	cout << "  h  = " << tv.h << " s\n\n";
// }

// doesn't need to be initialized. Overridden during runtime
struct RK45VarsPosVel
{
	double kx[6];
	double ky[6];
	double kz[6];
	double ku[6];
	double kv[6];
	double kw[6];

	double xi1[2];
	double yi1[2];
	double zi1[2];
	double ui1[2];
	double vi1[2];
	double wi1[2];

	double Rx;
	double Ry;
	double Rz;
	double Ru;
	double Rv;
	double Rw;

	double deltaX;
	double deltaY;
	double deltaZ;
	double deltaU;
	double deltaV;
	double deltaW;
};


const double RK45Coeff[8][6] =
	{ // k1         , k2          , k3          , k4           , k5      , k6
		{0.         , 0.          , 0.          , 0.           ,       0., 0.    },
		{1./4.      , 0.          , 0.          , 0.           ,       0., 0.    },
		{3./32.     , 9./32.      , 0.          , 0.           ,       0., 0.    },
		{1932./2197., -7200./2197., 7296./2197. , 0.           ,       0., 0.    },
		{439./216.  , -8.         , 3680./513.  , -845./4104.  ,       0., 0.    },
		{-8./27.    , 2.          , -3544./2565., 1859./4104.  , -11./40., 0.    },
		{25./216.   , 0.          , 1408./2565. , 2197./4104.  , -1./5.  , 0.    },
		{16./135.   , 0.          , 6656./12825., 28561./56430., -9./50. , 2./55.}
	};


void grav_a(double& ax, double& ay, double& az, double x, double y, double z)
{
	double grav_coeff = -gm * sqr(Rm) * pow(mag2(x, y, z), -1.5);

	ax = grav_coeff * x;
	ay = grav_coeff * y;
	az = grav_coeff * z;
}

struct cylinder
{
	double distance; // from impact point to center of lander, m
	double radius; // m
	double height; // starting from surface of Moon, m
};

void init_cylinder(cylinder& c, double d, double r, double h)
{
	c.distance = d;
	c.radius   = r;
	c.height   = h;
}

// functions ///////////////////////////////////////
// see https://math.okstate.edu/people/yqwang/teaching/math4513_fall11/Notes/rungekutta.pdf
void RK45UpdatePosVel(trackVars& t,
	                  RK45VarsPosVel& r,
	                  const double c[][6], cylinder& ll)
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


struct sums
{
	long double sum_miss;
	long double sum_hit;
	long long int N_miss;
	long long int N_hit;
};

void init_sums(sums& s)
{
	s.sum_miss = 0.0;
	s.sum_hit  = 0.0;
	s.N_miss = 0.;
	s.N_hit  = 0.;
}

double ejecta_distribution(trackVars& t)
{
	return pow(mag_s(t.u0, t.v0, t.w0), -(3.*mu_HH11 + 1.)) * fabs(sin(atan2(t.u0, t.w0)));
}

void inc_miss(sums& s, trackVars& t)
{
	s.N_miss++;

	// using Housen Holsapple 2011 speed distribution with isotropic term
	s.sum_miss += ejecta_distribution(t);
}

void inc_hit(sums& s, trackVars& t)
{
	s.N_hit++;

	// using Housen Holsapple 2011 speed distribution with isotropic term
	s.sum_hit += ejecta_distribution(t);
}


void printTrack(ofstream& file, trackVars& tv, sums& s)
{
	double dist = atan2(tv.xi, tv.zi); // units of lunar radius
	double speed = mag_s(tv.u0, 0.0, tv.w0) / vesc; // units of escape speed
	double zang = atan2(tv.u0, tv.w0);
	file << scientific << setprecision(14) << tv.xi  << ' ' << tv.zi << ' ' << dist << ' ' << speed << ' ' << zang
		 << ' ' << s.sum_miss << ' ' << s.sum_hit << ' ' << s.N_miss << ' ' << s.N_hit << endl;
}

bool checkCollision(trackVars& t, cylinder& c)
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





int main(int argc, char const *argv[])
{
	srand(time(0));

	double vmag, zang, dist;

	int N    = atoi(argv[1]);
	dist     = atof(argv[2]);
	int proc = atoi(argv[3]);
	double lander_radius = atof(argv[4]);
	double lander_height = atof(argv[5]);
	double vmin = atof(argv[6]); // m/s

	cout << "lander_radius = " << lander_radius << " m | lander_height = " << lander_height << " m | vmin = " << vmin << " m/s\n";

	RK45VarsPosVel RK45Vars;
	trackVars track_i;

	sums n;
	init_sums(n);

	cylinder lunarLander;
	bool hit = 0;

	init_cylinder(lunarLander,
	/* distance*/ dist,//Rm/252.3,
	/* radius  */ lander_radius, // 4.5, // m
	/* height  */ lander_height); // 50.);// m

	string sdist(to_string(proc));

	ofstream trackFile, distFile;
	//trackFile.open("tracks_" + sdist + ".txt");
	distFile.open("dists_" + sdist + ".txt");

	while(n.N_hit < N)
	{
		do{
			vmag = vmin/vesc + (1. - vmin/vesc) * rand()/double(RAND_MAX);
		} while(vmag == 0. || vmag == 1.);

		do {
			zang = 2*(rand()/double(RAND_MAX) - 0.5) * PI / 2.;
		} while(zang == -PI/2. || zang == PI/2. || zang == 0.);

		//cout << "vmag = " << vmag << " | zang = " << zang << endl;

		init_trackVars(track_i,
		/* x0 */       0.0,
		/* y0 */       0.0,
		/* z0 */       Rm,
		/* u0 */       vmag * sin(zang)*vesc,//0.5*vesc,
		/* v0 */       0.0,
		/* w0 */       vmag * cos(zang)*vesc, //0.8*vesc,,
		/* h  */       1e-4,
		/* dxmin */    lunarLander.radius/1000.);


		//printTrack(trackFile, track_i);

		for (int i = 0; i < NMAX; ++i)
		{
			RK45UpdatePosVel(track_i, RK45Vars, RK45Coeff, lunarLander);

			hit = checkCollision(track_i, lunarLander);

			//printTrack(trackFile, track_i);
			if(mag_s(track_i.xi, track_i.yi, track_i.zi) < Rm || hit)
				i = NMAX+1;

		}

		if (hit)
		{
			inc_hit(n, track_i);
			init_trackVars(track_i,
			/* x0 */       0.0,
			/* y0 */       0.0,
			/* z0 */       Rm,
			/* u0 */       vmag * sin(zang)*vesc,//0.5*vesc,
			/* v0 */       0.0,
			/* w0 */       vmag * cos(zang)*vesc, //0.8*vesc,,
			/* h  */       1.e-4,
			/* dxmin */    lunarLander.radius/1000.);

			//printTrack(trackFile, track_i);

			for (int i = 0; i < NMAX; ++i)
			{
				RK45UpdatePosVel(track_i, RK45Vars, RK45Coeff, lunarLander);

				hit = checkCollision(track_i, lunarLander);

				//printTrack(trackFile, track_i);
				if(mag_s(track_i.xi, track_i.yi, track_i.zi) < Rm || hit)
					i = NMAX+1;

			}
			//cout << "hit #: " << n.N_hit << " out of " << n.N_hit + n.N_miss << endl;
			printTrack(distFile, track_i, n);
		}
		else
		{
			inc_miss(n, track_i);
		}
		
	}

	
	
	//trackFile.close();
	distFile.close();
	return 0;
}