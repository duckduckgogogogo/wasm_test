// math_lib.h: math libraries

//////////////////////////////////////////////////////////////////////////////////
// Shin Dae-Seock (dsshin@midasIT.com)
// Don't believe in magic;
// understand what your libraries do, how they do it, and at what cost they do it.
//////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <cmath>
#include <float.h>
#include "MathFunc.h"

/////////////////////////////////////////////////////////////////////
namespace mesher_math_lib
{    
    const static double PI   = CDgnMathFunc::Get_pi();
    const static double ZERO = 1.0e-8; // DBL_EPSILON
    /////////////////////////////////////////////
    template<typename T> inline T mean_2(T t_value_1, T t_value_2) 
    {
        return (t_value_1 + t_value_2) / 2.0;
    } // end: mean_2()
    /////////////////////////////////////////////
    template<typename T> inline void mean_2(const T* a_array_1, const T* a_array_2, T* a_mean_array, int n_dimension)
    {
        for (int i=0; i<n_dimension; ++i) a_mean_array[i] = (a_array_1[i] + a_array_2[i]) / 2.0;
    } // end: mean_2()
    /////////////////////////////////////////////
    inline double radian(double d_degree) { return d_degree * PI / 180.0; }
    /////////////////////////////////////////////
    inline double degree(double d_radian) { return d_radian * 180.0 / PI; }
    /////////////////////////////////////////////
    inline int round_off(double d_value)
    {
        int n_value = (int)d_value;
        return ((d_value - n_value > 0.5 - ZERO) ?  n_value+1 : n_value);
    } // end: round_off()
    /////////////////////////////////////////////////////////////////
    inline int round_up(const double d_value)
    {
        int n_value = (int)d_value;
        return ((d_value - n_value > ZERO) ?  n_value+1 : n_value);
    } // end: round_up()
    /////////////////////////////////////////////
    inline bool in_range(double d_value, double d_limit_1, double d_limit_2)
    {
        return (d_value > d_limit_1 && d_value < d_limit_2);
    } // end: in_range()
    /////////////////////////////////////////////
    inline double norm_2d(const double* a_uv)
    {
        return sqrt(a_uv[0] * a_uv[0] + a_uv[1] * a_uv[1]);
    } // end: norm_2d()
    /////////////////////////////////////////////
    inline double squared_distance_2d(const double* a_uv_1, const double* a_uv_2)
    {
        return (a_uv_2[0]-a_uv_1[0]) * (a_uv_2[0]-a_uv_1[0]) + (a_uv_2[1]-a_uv_1[1]) * (a_uv_2[1]-a_uv_1[1]);
    } // end: squared_distance_2d()
    /////////////////////////////////////////////
    inline double squared_distance_3d(const double* a_xyz_1, const double* a_xyz_2)
    {
        return (a_xyz_2[0]-a_xyz_1[0]) * (a_xyz_2[0]-a_xyz_1[0]) + (a_xyz_2[1]-a_xyz_1[1]) * (a_xyz_2[1]-a_xyz_1[1]) + (a_xyz_2[2]-a_xyz_1[2]) * (a_xyz_2[2]-a_xyz_1[2]);
    } // end: squared_distance_2d()
    /////////////////////////////////////////////
    inline double distance_2d(const double* a_uv_1, const double* a_uv_2)
    {
        return sqrt(squared_distance_2d(a_uv_1, a_uv_2));
    } // end: distance_2d()
    /////////////////////////////////////////////
    inline double distance_3d(const double* a_xyz_1, const double* a_xyz_2)
    {
        return sqrt(squared_distance_3d(a_xyz_1, a_xyz_2));
    } // end: distance_2d()
    /////////////////////////////////////////////////////////////////////
    inline double triangle_area_2d(const double* a_uv_1, const double* a_uv_2, const double* a_uv_3)
    {
        return fabs(0.5 * (a_uv_1[0] * (a_uv_2[1] - a_uv_3[1]) + a_uv_2[0] * (a_uv_3[1] - a_uv_1[1]) + a_uv_3[0] * (a_uv_1[1] - a_uv_2[1])));
    }
    /////////////////////////////////////////////////////////////////////
    inline bool is_ccw_rotation_2d(const double* a_vector_1, const double* a_vector_2)
    {
        return (a_vector_1[0] * a_vector_2[1] - a_vector_2[0] * a_vector_1[1]) >= 0.0;
    } // end: is_ccw_rotation_2d()
    /////////////////////////////////////////////////////////////////////
    inline bool is_ccw_rotation_2d(const double* a_uv_1, const double* a_uv_2, const double* a_uv_3)
    {
        return (a_uv_1[0] * (a_uv_2[1] - a_uv_3[1]) + a_uv_2[0] * (a_uv_3[1] - a_uv_1[1]) + a_uv_3[0] * (a_uv_1[1] - a_uv_2[1])) >= 0.0;
    } // end: is_ccw_rotation_2d()
    /////////////////////////////////////////////////////////////////////
    double angle_2d(const double* a_vector_1, const double* a_vector_2);
    double angle_2d(const double* a_uv_0, const double* a_uv_1, const double* a_uv_2);
    double positive_angle_2d(const double* a_vector_1, const double* a_vector_2);
    double positive_angle_2d(const double* a_uv_0, const double* a_uv_1, const double* a_uv_2);
    void   translation_2d(double* a_uv, double d_delta_u, double d_delta_v);
    void   rotation_2d(double* a_uv, double d_angle);
    void   rotation_2d(double* a_uv, const double* a_origin_uv, double d_angle);
    void   mirror_2d(double* a_uv, const double* a_mirror_uv_1, const double* a_mirror_uv_2);
    bool   is_on_line_point_2d(const double* a_uv_1, const double* a_uv_2, const double* a_uv);
    bool   line_cross_2d(const double* a_uv_11, const double* a_uv_12, const double* a_uv_21, const double* a_uv_22);
    bool   line_intersection_2d(const double* a_uv_11, const double* a_uv_12, const double* a_uv_21, const double* a_uv_22, double& d_xi_1, double& d_xi_2);
    bool   line_intersection_2d(const double* a_uv_11, const double* a_uv_12, const double* a_uv_21, const double* a_uv_22, double* a_uv);
    bool   line_intersection_2d(const double* a_uv_11, const double* a_uv_12, const double* a_uv_21, const double* a_uv_22);
    /////////////////////////////////////////////////////////////////////
    void   circumscribed_circle_2d(const double* a_uv_1, const double* a_uv_2, const double* a_uv_3, double* a_center_uv, double& d_radius);
    double quadrilateral_angle_deviation_2d(const double* a_uv_1, const double* a_uv_2, const double* a_uv_3, const double* a_uv_4, double& d_min_angle, double& d_max_angle);
    /////////////////////////////////////////////////////////////////////
    double distortion_metric_2d(const double* a_uv_1, const double* a_uv_2, const double* a_uv_3);
    double distortion_metric_2d(const double* a_uv_1, const double* a_uv_2, const double* a_uv_3, const double* a_uv_4);
    /////////////////////////////////////////////////////////////////////
} // end: namespace mesher_math_lib
/////////////////////////////////////////////////////////////////////

