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



//// Helper functions
// prod_m = left_m * right_m
void h_matrix_matrix_multiply(mat3x3 &prod_m, mat3x3 &left_m, mat3x3 &right_m)
{
	for (int i = 0; i < 3; ++i)
	{
		h_matrix_vector_multiply(prod_m.col[i], left_m, right_m.col[i])
	}

}


// prod_v = left_m * right_v
void h_matrix_vector_multiply(vec3 &prod_v, mat3x3 &left_m, vec3 &right_v);
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
	rot_m.col[2].x[1] = sinp * sint
	rot_m.col[2].x[2] = cost;
}

void h_rot_m_from_ang_axis(mat3x3 &rot_m, double alpha, vec3 &u)
{
	double cosa = cos(alpha);
	double cosa_one_m = 1. - cos(alpha);
	double sina = sin(alpha);

	rot_m.col[0].x[0] = cosa + sqr(u.x[2]) * cosa_one_m;
	rot_m.col[0].x[1] = u[1]*u[0] * cosa_one_m + u[2] * sina;
	rot_m.col[0].x[2] = u[2]*u[0] * cosa_one_m - u[1] * sina;
	
	rot_m.col[1].x[0] = u[0]*u[1] * cosa_one_m - u[2] * sina;
	rot_m.col[1].x[1] = cosa + sqr(u[1]) * cosa_one_m;
	rot_m.col[1].x[2] = u[2]*u[1] * cosa_one_m + u[0] * sina;
	
	rot_m.col[2].x[0] = u[0]*u[2] * cosa_one_m + u[1] * sina;
	rot_m.col[2].x[1] = u[1]*u[2] * cosa_one_m - u[0] * sina;
	rot_m.col[2].x[2] = cosa + sqr(u[2]) * cosa_one_m;
}


void h_rot_m_from_angs_axes(mat3x3 &rot_m, double alpha, vec3 &u, double beta, vec3 &v)
{
	mat3x3 Ru, Rv;

	h_rot_m_from_ang_axis(Ru, alpha, u);
	h_rot_m_from_ang_axis(Rv, beta , v);

	h_matrix_matrix_multiply(rot_m, Ru, Rv);
}