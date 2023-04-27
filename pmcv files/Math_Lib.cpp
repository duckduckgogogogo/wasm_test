// math_lib.cpp: math libraries

//////////////////////////////////////////////////////////////////////////////////
// Shin Dae-Seock (dsshin@midasIT.com)
// Don't believe in magic;
// understand what your libraries do, how they do it, and at what cost they do it.
//////////////////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Math_Lib.h"

#include <functional>
#include <algorithm>

using namespace std;

/////////////////////////////////////////////////////////////////////
namespace mesher_math_lib
{
    /////////////////////////////////////////////////////////////////
    double angle_2d(const double* a_vector_1, const double* a_vector_2)
    {
        double a_vector_3[2];
        transform(a_vector_2, a_vector_2+2, a_vector_1, a_vector_3, minus<double>());
        double d_norm_1      = norm_2d(a_vector_1);
        double d_norm_2      = norm_2d(a_vector_2);
        double d_norm_3      = norm_2d(a_vector_3);
        double d_numerator   = pow(d_norm_1, 2.0) + pow(d_norm_2, 2.0) - pow(d_norm_3, 2.0);
        double d_denominator = 2.0 * d_norm_1 * d_norm_2;
        /////////////////////////////////////////
        if (d_denominator < ZERO) return 0.0;
        /////////////////////////////////////////
        double d_cosine = d_numerator / d_denominator;
        if (d_cosine > ( 1.0 - ZERO)) return   0.0;
        if (d_cosine < (-1.0 + ZERO)) return 180.0;
        /////////////////////////////////////////
        double d_angle = fabs(degree(acos(d_cosine)));
        return is_ccw_rotation_2d(a_vector_1, a_vector_2) ? d_angle : -d_angle;
    } // end: angle_2d()
    /////////////////////////////////////////////////////////////////
    double angle_2d(const double* a_uv_0, const double* a_uv_1, const double* a_uv_2)
    {
        double a_vector_01[2], a_vector_02[2];
        transform(a_uv_1, a_uv_1+2, a_uv_0, a_vector_01, minus<double>());
        transform(a_uv_2, a_uv_2+2, a_uv_0, a_vector_02, minus<double>());
        return angle_2d(a_vector_01, a_vector_02);
    } // end: angle_2d()
    /////////////////////////////////////////////////////////////////
    double positive_angle_2d(const double* a_vector_1, const double* a_vector_2)
    {
        double d_angle = angle_2d(a_vector_1, a_vector_2);
        if (d_angle < 0.0) d_angle += 360.0;
        return d_angle;
    } // end: positive_angle_2d()
    /////////////////////////////////////////////////////////////////
    double positive_angle_2d(const double* a_uv_0, const double* a_uv_1, const double* a_uv_2)
    {
        double a_vector_01[2], a_vector_02[2];
        transform(a_uv_1, a_uv_1+2, a_uv_0, a_vector_01, minus<double>());
        transform(a_uv_2, a_uv_2+2, a_uv_0, a_vector_02, minus<double>());
        double d_angle = angle_2d(a_vector_01, a_vector_02);
        return d_angle < 0.0 ? d_angle + 360.0 : d_angle;
    } // end: positive_angle_2d()
    /////////////////////////////////////////////////////////////////
    void translation_2d(double* a_uv, double d_delta_u, double d_delta_v)
    {
        a_uv[0] += d_delta_u;
        a_uv[1] += d_delta_v;
    } // end: translation_2d()
    /////////////////////////////////////////////////////////////////
    void rotation_2d(double* a_uv, double d_angle)
    {
        double       a_uv_0[2] = { a_uv[0], a_uv[1] };
        const double d_sine    = sin(radian(d_angle));
        const double d_cosine  = cos(radian(d_angle));
        a_uv[0]                = d_cosine * a_uv_0[0] - d_sine   * a_uv_0[1];
        a_uv[1]                = d_sine   * a_uv_0[0] + d_cosine * a_uv_0[1];
    } // end: rotation_2d()
    /////////////////////////////////////////////////////////////////
    void rotation_2d(double* a_uv, const double* a_origin_uv, double d_angle)
    {
        double d_cosine  = cos(radian(d_angle));
        double d_sine    = sin(radian(d_angle));
        double a_delta_u = a_uv[0] - a_origin_uv[0];
        double a_delta_v = a_uv[1] - a_origin_uv[1];
        a_uv[0]          = d_cosine * a_delta_u - d_sine   * a_delta_v + a_origin_uv[0];
        a_uv[1]          = d_sine   * a_delta_u + d_cosine * a_delta_v + a_origin_uv[1];
    } // end: rotation_2d()
    /////////////////////////////////////////////////////////////////
    void mirror_2d(double* a_uv, const double* a_mirror_uv_1, const double* a_mirror_uv_2)
    {
        if (squared_distance_2d(a_mirror_uv_1, a_mirror_uv_2) < ZERO) return;
        /////////////////////////////////////////////
        double a_origin_uv[2];
        if      (fabs(a_mirror_uv_2[0]-a_mirror_uv_1[0]) < ZERO)
        {
            a_origin_uv[0] = a_mirror_uv_1[0];
            a_origin_uv[1] = a_uv[1];
        }
        else if (fabs(a_mirror_uv_2[1]-a_mirror_uv_1[1]) < ZERO)
        {
            a_origin_uv[0] = a_uv[0];
            a_origin_uv[1] = a_mirror_uv_1[1];
        }
        else
        {
            double a       = (a_mirror_uv_2[1] - a_mirror_uv_1[1]) / (a_mirror_uv_2[0] - a_mirror_uv_1[0]);
            double b       = a_mirror_uv_1[1] - a * a_mirror_uv_1[0];
            a_origin_uv[0] = (a_uv[1] + a_uv[0] / a - b) / (a + 1.0 / a);
            a_origin_uv[1] = a * a_origin_uv[0] + b;
        }
        rotation_2d(a_uv, a_origin_uv, 180.0);
    } // end: mirror_2d()
    /////////////////////////////////////////////////////////////////
    bool is_on_line_point_2d(const double* a_uv_1, const double* a_uv_2, const double* a_uv)
    {
        if (squared_distance_2d(a_uv_1, a_uv) < ZERO) return true;
        if (squared_distance_2d(a_uv_2, a_uv) < ZERO) return true;
        double d_angle = positive_angle_2d(a_uv, a_uv_1, a_uv_2);
        return (d_angle > 179.0 && d_angle < 181.0);
    } // end: is_on_line_point_2d()
    /////////////////////////////////////////////////////////////////
    bool line_cross_2d(const double* a_uv_11, const double* a_uv_12, const double* a_uv_21, const double* a_uv_22)
    {
        double d_cross_product_1 = (a_uv_11[0]-a_uv_21[0]) * (a_uv_12[1]-a_uv_21[1]) - (a_uv_11[1]-a_uv_21[1]) * (a_uv_12[0]-a_uv_21[0]);
        double d_cross_product_2 = (a_uv_11[0]-a_uv_22[0]) * (a_uv_12[1]-a_uv_22[1]) - (a_uv_11[1]-a_uv_22[1]) * (a_uv_12[0]-a_uv_22[0]);
        double d_product         = d_cross_product_1 * d_cross_product_2;
        ////////////////////////////////////
        if (fabs(d_product) < ZERO || d_product >= 0.0) return false;
        ////////////////////////////////////
        d_cross_product_1 = (a_uv_21[0]-a_uv_11[0]) * (a_uv_22[1]-a_uv_11[1]) - (a_uv_21[1]-a_uv_11[1]) * (a_uv_22[0]-a_uv_11[0]);
        d_cross_product_2 = (a_uv_21[0]-a_uv_12[0]) * (a_uv_22[1]-a_uv_12[1]) - (a_uv_21[1]-a_uv_12[1]) * (a_uv_22[0]-a_uv_12[0]);
        d_product         = d_cross_product_1 * d_cross_product_2;
        return (fabs(d_product) > ZERO && d_product < 0.0);
    } // end: line_cross_2d()
    /////////////////////////////////////////////////////////////////
    bool line_intersection_2d(const double* a_uv_11, const double* a_uv_12, const double* a_uv_21, const double* a_uv_22, double& d_xi_1, double& d_xi_2)
    {
        double d_determinant = (a_uv_21[0] - a_uv_22[0]) * (a_uv_11[1] - a_uv_12[1]) - (a_uv_11[0] - a_uv_12[0]) * (a_uv_21[1] - a_uv_22[1]);
        if (fabs(d_determinant) < ZERO) return false;
        d_xi_1 = ((a_uv_22[0] - a_uv_21[0]) * (a_uv_21[1] - a_uv_11[1]) + (a_uv_21[0] - a_uv_11[0]) * (a_uv_21[1] - a_uv_22[1])) / d_determinant;
        /////////////////////////////////////////
        if (d_xi_1 <= -ZERO || d_xi_1 >= 1.0 + ZERO) return false;
        /////////////////////////////////////////
        d_xi_2 = ((a_uv_21[0] - a_uv_11[0]) * (a_uv_11[1] - a_uv_12[1]) + (a_uv_12[0] - a_uv_11[0]) * (a_uv_21[1] - a_uv_11[1])) / d_determinant;
        return (d_xi_2 > -ZERO && d_xi_2 < 1.0 + ZERO);
    } // end: line_intersection_2d()
    /////////////////////////////////////////////////////////////////
    bool line_intersection_2d(const double* a_uv_11, const double* a_uv_12, const double* a_uv_21, const double* a_uv_22, double* a_uv)
    {
        double d_xi_1, d_xi_2;
        if (!line_intersection_2d(a_uv_11, a_uv_12, a_uv_21, a_uv_22, d_xi_1, d_xi_2)) return false;
        a_uv[0] = a_uv_11[0] + d_xi_1 * (a_uv_12[0] - a_uv_11[0]);
        a_uv[1] = a_uv_11[1] + d_xi_1 * (a_uv_12[1] - a_uv_11[1]);
        return true;
    } // end: line_intersection_2d()
    /////////////////////////////////////////////////////////////////
    bool line_intersection_2d(const double* a_uv_11, const double* a_uv_12, const double* a_uv_21, const double* a_uv_22)
    {
        double d_xi_1, d_xi_2;
        return line_intersection_2d(a_uv_11, a_uv_12, a_uv_21, a_uv_22, d_xi_1, d_xi_2);
    } // end: line_intersection_2d()
    /////////////////////////////////////////////////////////////////
    void circumscribed_circle_2d(const double* a_uv_1, const double* a_uv_2, const double* a_uv_3, double* a_center_uv, double& d_radius)
    {
        double a_middle_12[2], a_middle_23[2];
        mean_2<double>(a_uv_1, a_uv_2, a_middle_12, 2);
        mean_2<double>(a_uv_2, a_uv_3, a_middle_23, 2);
        double d_numerator   = (a_middle_12[0] - a_middle_23[0]) * (a_uv_2[0] - a_uv_1[0]) + (a_middle_12[1] - a_middle_23[1]) * (a_uv_2[1] - a_uv_1[1]);
        double d_denominator = (a_uv_2[0] - a_uv_1[0]) * (a_uv_3[1] - a_uv_2[1]) - (a_uv_3[0] - a_uv_2[0]) * (a_uv_2[1] - a_uv_1[1]);
        /////////////////////////////////////////
        if (fabs(d_denominator) < ZERO)
        {
            fill(a_center_uv, a_center_uv+2, 0.0);
            d_radius = 0.0;
            return;
        }
        /////////////////////////////////////////
        a_center_uv[0]  = a_middle_23[0] + d_numerator / d_denominator * (a_uv_3[1] - a_uv_2[1]);
        a_center_uv[1]  = a_middle_23[1] - d_numerator / d_denominator * (a_uv_3[0] - a_uv_2[0]);
        d_radius        = distance_2d(a_uv_1, a_center_uv);
    } // end: circumscribed_circle_2d()
    /////////////////////////////////////////////////////////////////
    double quadrilateral_angle_deviation_2d(const double* a_uv_1, const double* a_uv_2, const double* a_uv_3, const double* a_uv_4, double& d_min_angle, double& d_max_angle)
    {
        double a_angle[4];
        a_angle[3]  = a_angle[0] = angle_2d(a_uv_1, a_uv_2, a_uv_4);
        a_angle[3] += a_angle[1] = angle_2d(a_uv_2, a_uv_3, a_uv_1);
        a_angle[3] += a_angle[2] = angle_2d(a_uv_3, a_uv_4, a_uv_2);
        a_angle[3]               = 360.0 - a_angle[3];
        d_min_angle              = *min_element(a_angle, a_angle+4);
        d_max_angle              = *max_element(a_angle, a_angle+4);
        double d_angle_deviation = 0.0;
        for (int i=0; i<4; ++i) d_angle_deviation += fabs(90.0 - a_angle[i]);
        return d_angle_deviation / 4.0;
    } // end: quadrilateral_angle_deviation_2d()
    /////////////////////////////////////////////////////////////////
    double distortion_metric_2d(const double* a_uv_1, const double* a_uv_2, const double* a_uv_3)
    {
        const static double d_normal_factor = 2.0 * sqrt(3.0);
        /////////////////////////////////////////
        double a_vector_31[2], a_vector_32[2];
        transform(a_uv_1, a_uv_1+2, a_uv_3, a_vector_31, minus<double>());
        transform(a_uv_2, a_uv_2+2, a_uv_3, a_vector_32, minus<double>());
        double d_cross_product     = a_vector_31[0] * a_vector_32[1] - a_vector_32[0] * a_vector_31[1];
        double d_squared_length_12 = squared_distance_2d(a_uv_1, a_uv_2);
        double d_squared_length_23 = squared_distance_2d(a_uv_2, a_uv_3);
        double d_squared_length_31 = squared_distance_2d(a_uv_3, a_uv_1);
        return d_normal_factor * d_cross_product / (d_squared_length_12 + d_squared_length_23 + d_squared_length_31);
    } // end: distortion_metric_2d()
    /////////////////////////////////////////////////////////////////
    double distortion_metric_2d(const double* a_uv_1, const double* a_uv_2, const double* a_uv_3, const double* a_uv_4)
    {
        double a_triagle_metric[4];
        a_triagle_metric[0] = distortion_metric_2d(a_uv_1, a_uv_2, a_uv_3);
        a_triagle_metric[1] = distortion_metric_2d(a_uv_1, a_uv_3, a_uv_4);
        a_triagle_metric[2] = distortion_metric_2d(a_uv_1, a_uv_2, a_uv_4);
        a_triagle_metric[3] = distortion_metric_2d(a_uv_2, a_uv_3, a_uv_4);
        sort(a_triagle_metric, a_triagle_metric+4, greater<double>());
        double d_sign = (a_triagle_metric[3] >= 0.0) ? 1.0 : -1.0;
        return d_sign * (a_triagle_metric[2] * a_triagle_metric[3]) / (a_triagle_metric[0] * a_triagle_metric[1]);
    } // end: distortion_metric_2d()
} // end: namespace mesher_math_lib
/////////////////////////////////////////////////////////////////////
