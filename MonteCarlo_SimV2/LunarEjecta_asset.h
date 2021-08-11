#ifndef LUNAREJECTA_ASSET_H
#define LUNAREJECTA_ASSET_H

#include <vector>

using namespace std;

struct shape
{
	double x_origin_offset; // m, offset from asset main origin
	double y_origin_offset; // m, offset from asset main origin
	double z_origin_offset; // m, offset from asset main origin

	double z_axis_tilt_theta; // rad, angle from global z-axis to local z-axis, wrt asset main origin
	double z_axis_tilt_phi;   // rad, angle from glabal z-axis to local x-axis, wrt asset main origin
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
	vector<rect_prism> re;
}



bool check_collision_asset(asset &a, double x, double y, double z);

#endif 