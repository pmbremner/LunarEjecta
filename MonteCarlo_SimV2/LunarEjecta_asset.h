#ifndef LUNAREJECTA_ASSET_H
#define LUNAREJECTA_ASSET_H

#include "LunarEjecta_params.h"
//#include "LunarEjecta_trajectory.h"

#include <vector>
#include <string>

using namespace std;


struct shape
{
	vec3 origin_offset_nom; // m, offset from asset main origin, with no rotations
	// origin_offset_rot = rot_m_asset * origin_offset_nom
	vec3 origin_offset_rot; // m, offset from asset main origin, with asset rotations

	double z_axis_tilt_theta; // rad, angle from asset z-axis to local z-axis
	double z_axis_tilt_phi;   // rad, angle from asset x-axis to local x-axis

	// rot_m_tot.col[0] = x'' axis
	// rot_m_tot.col[1] = y'' axis
	// rot_m_tot.col[2] = z'' axis
	mat3x3 rot_m_tot; // rotation matrix 3x3, columns first, both asset and shape together
};

struct rect_prism : shape
{
	double x_width; // m
	double y_width; // m
	double z_width; // m
};

struct cylinder : shape
{
	double z_height; // m
	double s_radius; // m
};

struct sphere : shape
{
	double r_radius; // m
};


struct asset
{
	unsigned int N_shapes; // number of shapes that define the asset
	
	vector<unsigned int> shapes; // id of shape, 0 = sphere, 1 = cylinder, 2 = rect_prism
	vector<unsigned int> shape_idx; // index of particular shape in the sphere, cylinder, or rect_prism vectors

	vector<sphere> sp;
	vector<cylinder> cy;
	vector<rect_prism> rp;

	//// asset origin and orientation are defined by current state of the trajectory
	// asset main origin
	vec3 origin; // m, from global origin
	double collision_radius_boundary; // m, only check collisions if within this radius
	
	double z_axis_tilt_theta; // rad, angle from global z-axis to asset z-axis
	double z_axis_tilt_phi;   // rad, angle from glabal z-axis to asset x-axis

	// rot_m_asset.col[0] = x' axis
	// rot_m_asset.col[1] = y' axis
	// rot_m_asset.col[2] = z' axis
	mat3x3 rot_m_asset; // rotation matrix 3x3, columns first, only asset
};

// prod_m = left_m * right_m
void h_matrix_matrix_multiply(mat3x3 &prod_m, mat3x3 &left_m, mat3x3 &right_m);
// prod_v = left_m * right_v
void h_matrix_vector_multiply(vec3 &prod_v, mat3x3 &left_m, vec3 &right_v);
// wrt global axis
void h_rot_m_from_angs(mat3x3 &rot_m, double theta, double phi);
// wrt to asset axis, u
void h_rot_m_from_ang_axis(mat3x3 &rot_m, double alpha, vec3 &u);
// wrt to asset axes, u and v
void h_rot_m_from_angs_axes(mat3x3 &rot_m, double alpha, vec3 &u, double beta, vec3 &v);

// read asset definition from file and initialize the asset data, compute shape rot_m
void init_asset(asset &a, string asset_fn);

// set asset origin and orientation, compute rot_m for asset and tot in each shape
// Note: Not now, but may be possible to read new asset as a function of trajectory
////void set_asset_trajectory_state(asset &a, trajectory &t, int i);


bool check_collision_rect_prism(rect_prism &rp, vec3 &pos);
bool check_collision_cylinder(cylinder &cy,     vec3 &pos);
bool check_collision_sphere(sphere &sp,         vec3 &pos);

// loops through each shape of the asset to check collisions, will stop early if collision detected
bool check_collision_asset(asset &a, vec3 &pos);




#endif 