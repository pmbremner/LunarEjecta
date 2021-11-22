#include "LunarEjecta_asset.h"

using namespace std;

// struct vec3
// {
// 	double x[3];
// };

// struct mat3x3
// {
// 	vec3 col[3];
// };


void init_asset(asset &a, string fn)
{
	int N_sphere, N_cylinder, N_rect_prism, i, n;
	double asset_lon, asset_lat, asset_height;
	mat3x3 temp_rot_m_asset, temp_rot_m_asset_latlon;

	cout << "--------------------------------\n";
	cout << "Reading... " << fn << endl;

	//getParam(fn, "N_shapes", a.N_shapes, 0);
	
	getParam(fn, "asset_lat", asset_lat, 0); // deg
	getParam(fn, "asset_lon", asset_lon, 0); // deg
	getParam(fn, "asset_height", asset_height, 0); // m

	getParam(fn, "collision_radius_boundary", a.collision_radius_boundary, 0); // m
	
 
	//// get asset orientation 
	getParam(fn, "orientation", a.orientation, 0);

	if (a.orientation == 1)
	{ // z axis is user defined wrt global frame
		getParam(fn, "y_axis_rot_theta", a.y_axis_rot_theta, 0);
		getParam(fn, "z_axis_rot_phi", a.z_axis_rot_phi, 0);

		// compute asset rot matrix directly
		h_rot_m_from_angs(a.rot_m_asset, a.y_axis_rot_theta, a.z_axis_rot_phi);

	}
	else if (a.orientation == 2)
	{ // z axis is user defined wrt local surface plane
		getParam(fn, "y_axis_rot_theta", a.y_axis_rot_theta, 0);
		getParam(fn, "z_axis_rot_phi", a.z_axis_rot_phi, 0);

		// compute asset temp rot matrix of asset, then lat-lon contribution, and then put together
		h_rot_m_from_angs(temp_rot_m_asset, a.y_axis_rot_theta, a.z_axis_rot_phi);
		h_rot_m_from_angs(temp_rot_m_asset_latlon, PI/2. - asset_lat, asset_lon);
		h_matrix_matrix_multiply(a.rot_m_asset, temp_rot_m_asset_latlon, temp_rot_m_asset);
	}
	else // assuming a.orientation == 0
	{ // z axis is normal to the local surface
		if (a.orientation != 0)
			cout << "WARNING: Assuming asset orientation to be normal to local surface!\n";

		a.y_axis_rot_theta = 0.;
		a.z_axis_rot_phi = 0.;

		h_rot_m_from_angs(a.rot_m_asset, PI/2. - asset_lat, asset_lon);
	}

	//// compute asset origin
	a.origin.x[0] = asset_height * sin(PI/2. - asset_lat) * cos(asset_lon);
	a.origin.x[1] = asset_height * sin(PI/2. - asset_lat) * sin(asset_lon);
	a.origin.x[2] = asset_height * cos(PI/2. - asset_lat);

	// temp get number of each shape and sum to N_shapes
	getParam(fn, "N_sphere", N_sphere, 0);
	getParam(fn, "N_cylinder", N_cylinder, 0);
	getParam(fn, "N_rect_prism", N_rect_prism, 0);

	a.N_shapes = N_sphere + N_cylinder + N_rect_prism;

	// start with fresh vectors
	a.shapes.clear();
	a.shape_idx.clear();
	a.sp.clear();
	a.cy.clear();
	a.rp.clear();

	// allocate space for vectors, still need to init values in them
	a.shapes.resize(a.N_shapes);
	a.shape_idx.resize(a.N_shapes);
	a.sp.resize(N_sphere);
	a.cy.resize(N_cylinder);
	a.rp.resize(N_rect_prism);

	n = 0;
	for (i = 0; i < N_sphere; ++i)
	{
		a.shapes[n] = 0;

		a.shape_idx[n] = i;
		n++;
	}

	for (i = 0; i < N_cylinder; ++i)
	{
		a.shapes[n] = 1;

		a.shape_idx[n] = i;
		n++;
	}

	for (i = 0; i < N_rect_prism; ++i)
	{
		a.shapes[n] = 2;

		a.shape_idx[n] = i;
		n++;
	}



	cout << "--------------------------------\n";
}


//// Helper functions
// prod_m = left_m * right_m
void h_matrix_matrix_multiply(mat3x3 &prod_m, mat3x3 &left_m, mat3x3 &right_m)
{
	for (int i = 0; i < 3; ++i)
	{
		h_matrix_vector_multiply(prod_m.col[i], left_m, right_m.col[i]);
	}

}


// prod_v = left_m * right_v
void h_matrix_vector_multiply(vec3 &prod_v, mat3x3 &left_m, vec3 &right_v)
{
	int j, k;
	for (j = 0; j < 3; ++j)
	{
		prod_v.x[j] = 0.;
		for (k = 0; k < 3; ++k)
		{
			prod_v.x[j] += left_m.col[k].x[j] * right_v.x[k];
		}
	}
}

// https://en.wikipedia.org/wiki/Rotation_matrix
// theta rotation about RH y-axis, phi rotation about RH z-axis 
void h_rot_m_from_angs(mat3x3 &rot_m, double theta, double phi)
{
	double cost = cos(theta);
	double sint = sin(theta);
	double cosp = cos(phi);
	double sinp = sin(phi);

	rot_m.col[0].x[0] = cosp * cost;
	rot_m.col[0].x[1] = sinp * cost;
	rot_m.col[0].x[2] = -sint;
	
	rot_m.col[1].x[0] = -sinp;
	rot_m.col[1].x[1] = cosp;
	rot_m.col[1].x[2] = 0.;
	
	rot_m.col[2].x[0] = cosp * sint;
	rot_m.col[2].x[1] = sinp * sint;
	rot_m.col[2].x[2] = cost;
}

void h_rot_m_from_ang_axis(mat3x3 &rot_m, double alpha, vec3 &u)
{
	double cosa = cos(alpha);
	double cosa_one_m = 1. - cos(alpha);
	double sina = sin(alpha);

	rot_m.col[0].x[0] = cosa + sqr(u.x[2]) * cosa_one_m;
	rot_m.col[0].x[1] = u.x[1]*u.x[0] * cosa_one_m + u.x[2] * sina;
	rot_m.col[0].x[2] = u.x[2]*u.x[0] * cosa_one_m - u.x[1] * sina;
	
	rot_m.col[1].x[0] = u.x[0]*u.x[1] * cosa_one_m - u.x[2] * sina;
	rot_m.col[1].x[1] = cosa + sqr(u.x[1]) * cosa_one_m;
	rot_m.col[1].x[2] = u.x[2]*u.x[1] * cosa_one_m + u.x[0] * sina;
	
	rot_m.col[2].x[0] = u.x[0]*u.x[2] * cosa_one_m + u.x[1] * sina;
	rot_m.col[2].x[1] = u.x[1]*u.x[2] * cosa_one_m - u.x[0] * sina;
	rot_m.col[2].x[2] = cosa + sqr(u.x[2]) * cosa_one_m;
}


void h_rot_m_from_angs_axes(mat3x3 &rot_m, double alpha, vec3 &u, double beta, vec3 &v)
{
	mat3x3 Ru, Rv;

	h_rot_m_from_ang_axis(Ru, alpha, u);
	h_rot_m_from_ang_axis(Rv, beta , v);

	h_matrix_matrix_multiply(rot_m, Ru, Rv);
}



bool check_collision_asset(asset &a, double x, double y, double z)
{
	return 0;
}

bool check_collision_moon(double r, double x, double y, double z)
{
	return (mag_s(x, y, z) <= r ? 1 : 0);
}

bool check_escape(asset &a, double vesc, double x, double y, double z, double u, double v, double w)
{
	if (mag_s(u, v, w) >= vesc)
		return (mag_s(x, y, z) > mag_s(a.origin) + a.collision_radius_boundary ? 1 : 0);
	else
		return 0;
}