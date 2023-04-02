#include <iostream>
#include <iomanip>
#include "quartic.h"
#include "outliers.h"
#include <math.h>
#include <chrono>
#include <thread>

//**********************************
// We assume that the input terrain data spread on x and y coordinates, that is,
// max_x - min_x >> max_z - min_z and max_y - min_y >> max_z - min_z,
// and the normal vector for the whole terrain data is (0,0,1), which is pointing
// to z coordinate. We need this because we need to define the left and right
// on the terrain data
//**********************************

double round_if_close(double a)
{
    double b = round(a);
    if (((a - b) >= -1e-6 && a < b) ||
        ((a - b) <= 1e-6 && a > b))
    {
        a = b;
    }
    return a;
}

double max_1(double a, double b)
{
    if (a >= b)
    {
        return a;
    }
    else
    {
        return b;
    }
}
double min_1(double a, double b)
{
    if (a <= b)
    {
        return a;
    }
    else
    {
        return b;
    }
}

// check whether a vertex is in the edge
bool vertex_in_edge(double edge_1_x, double edge_2_x,
                    double edge_1_y, double edge_2_y,
                    double edge_1_z, double edge_2_z,
                    double vertex_x, double vertex_y,
                    double vertex_z)
{
    return (vertex_x >= min_1(edge_1_x, edge_2_x) &&
            vertex_x <= max_1(edge_1_x, edge_2_x) &&
            vertex_y >= min_1(edge_1_y, edge_2_y) &&
            vertex_y <= max_1(edge_1_y, edge_2_y) &&
            vertex_z >= min_1(edge_1_z, edge_2_z) &&
            vertex_z <= max_1(edge_1_z, edge_2_z));
}

// check whether a vertex is in the edge (not include two endpoint vertices)
bool vertex_in_edge_not_endpoint(double edge_1_x, double edge_2_x,
                                 double edge_1_y, double edge_2_y,
                                 double edge_1_z, double edge_2_z,
                                 double vertex_x, double vertex_y,
                                 double vertex_z)
{
    return (vertex_x >= min_1(edge_1_x, edge_2_x) &&
            vertex_x <= max_1(edge_1_x, edge_2_x) &&
            vertex_y >= min_1(edge_1_y, edge_2_y) &&
            vertex_y <= max_1(edge_1_y, edge_2_y) &&
            vertex_z >= min_1(edge_1_z, edge_2_z) &&
            vertex_z <= max_1(edge_1_z, edge_2_z) &&
            (!(vertex_x == edge_1_x && vertex_y == edge_1_y && vertex_z == edge_1_z)) &&
            (!(vertex_x == edge_2_x && vertex_y == edge_2_y && vertex_z == edge_2_z)));
}

double distance_xyz(double x_1, double y_1, double z_1,
                    double x_2, double y_2, double z_2)
{
    double d_x = x_2 - x_1;
    double d_y = y_2 - y_1;
    double d_z = z_2 - z_1;

    return sqrt(d_x * d_x + d_y * d_y + d_z * d_z);
}

// check whether a vertex is on the edge
bool vertex_on_edge(double edge_1_x, double edge_2_x,
                    double edge_1_y, double edge_2_y,
                    double edge_1_z, double edge_2_z,
                    double vertex_x, double vertex_y,
                    double vertex_z)
{
    if (!vertex_in_edge(edge_1_x, edge_2_x, edge_1_y, edge_2_y,
                        edge_1_z, edge_2_z, vertex_x, vertex_y, vertex_z))
    {
        return false;
    }
    if (distance_xyz(edge_1_x, edge_1_y, edge_1_z, vertex_x, vertex_y, vertex_z) +
            distance_xyz(edge_2_x, edge_2_y, edge_2_z, vertex_x, vertex_y, vertex_z) -
            distance_xyz(edge_1_x, edge_1_y, edge_1_z, edge_2_x, edge_2_y, edge_2_z) >
        0)
    {
        return false;
    }
    return true;
}

// calcualte the cross product for two vectors
void vector_cross_product(double a_x, double a_y, double a_z,
                          double b_x, double b_y, double b_z,
                          double &result_x, double &result_y, double &result_z)
{
    result_x = a_y * b_z - a_z * b_y;
    result_y = a_z * b_x - a_x * b_z;
    result_z = a_x * b_y - a_y * b_x;
}

// calcualte the dot product for two vectors
void vector_dot_product(double a_x, double a_y, double a_z,
                        double b_x, double b_y, double b_z,
                        double &result)
{
    result = a_x * b_x + a_y * b_y + a_z * b_z;
}

// calculate the normal vector for a given face
void calculate_normal_vector_of_face(double face_x_1, double face_y_1, double face_z_1,
                                     double face_x_2, double face_y_2, double face_z_2,
                                     double face_x_3, double face_y_3, double face_z_3,
                                     double &normal_vector_x, double &normal_vector_y,
                                     double &normal_vector_z)
{
    double vector_a_x = face_x_2 - face_x_1;
    double vector_a_y = face_y_2 - face_y_1;
    double vector_a_z = face_z_2 - face_z_1;
    double vector_b_x = face_x_3 - face_x_1;
    double vector_b_y = face_y_3 - face_y_1;
    double vector_b_z = face_z_3 - face_z_1;

    vector_cross_product(vector_a_x, vector_a_y, vector_a_z,
                         vector_b_x, vector_b_y, vector_b_z,
                         normal_vector_x, normal_vector_y, normal_vector_z);
}

// this calculate whether the checking point is on the left of the ray, we first need
// to calculate the cross product from the ray_from_point_to_ray_to_vector to the
// ray_from_point_to_checking_point_vector, if this cross product is the same as the
// direction of the whole terrain data normal vector (0,0,1), then use the right hand
// rule, the checking_point is on the left of the ray_from_point_to_ray_to_vector
bool point_on_left_of_ray(double ray_to_point_x, double ray_to_point_y, double ray_to_point_z,
                          double ray_from_point_x, double ray_from_point_y,
                          double ray_from_point_z, double checking_point_x,
                          double checking_point_y, double checking_point_z)
{

    // we first need to make sure that the direction for last_face_normal_vector and first_face_normal_vector
    // is the same as the whole terrain data normal vector (0,0,1)

    // ray_from_point_to_ray_to_vector is denoted as vector A
    double ray_from_point_to_ray_to_vector_x = ray_to_point_x - ray_from_point_x;
    double ray_from_point_to_ray_to_vector_y = ray_to_point_y - ray_from_point_y;
    double ray_from_point_to_ray_to_vector_z = ray_to_point_z - ray_from_point_z;

    // ray_from_point_to_checking_point_vector is denoted as vector B
    double ray_from_point_to_checking_point_vector_x = checking_point_x - ray_from_point_x;
    double ray_from_point_to_checking_point_vector_y = checking_point_y - ray_from_point_y;
    double ray_from_point_to_checking_point_vector_z = checking_point_z - ray_from_point_z;

    double A_and_B_vector_cross_product_x;
    double A_and_B_vector_cross_product_y;
    double A_and_B_vector_cross_product_z;

    vector_cross_product(ray_from_point_to_ray_to_vector_x, ray_from_point_to_ray_to_vector_y,
                         ray_from_point_to_ray_to_vector_z, ray_from_point_to_checking_point_vector_x,
                         ray_from_point_to_checking_point_vector_y, ray_from_point_to_checking_point_vector_z,
                         A_and_B_vector_cross_product_x, A_and_B_vector_cross_product_y, A_and_B_vector_cross_product_z);

    // if the A_and_B_vector_cross_product is pointing out, then use the right hand rule,
    // the checking_point is on the left of the ray_from_point_to_ray_to_vector
    if (A_and_B_vector_cross_product_z > 0)
    {
        return true;
    }
    return false;
}

// this calculate whether the checking point is on the left, right, or vertical of the ray,
// we first need to calculate the cross product from the ray_from_point_to_ray_to_vector to the
// ray_from_point_to_checking_point_vector, if this cross product is the same as the
// direction of the whole terrain data normal vector (0,0,1), i.e. > 0, then use the right hand
// rule, the checking_point is on the left of the ray_from_point_to_ray_to_vector, if < 0, then
// on the right, if = 0, then vertical
double point_on_left_right_vertical_of_ray(double ray_to_point_x, double ray_to_point_y, double ray_to_point_z,
                                           double ray_from_point_x, double ray_from_point_y,
                                           double ray_from_point_z, double checking_point_x,
                                           double checking_point_y, double checking_point_z)
{
    // we first need to make sure that the direction for last_face_normal_vector and first_face_normal_vector
    // is the same as the whole terrain data normal vector (0,0,1)

    // ray_from_point_to_ray_to_vector is denoted as vector A
    double ray_from_point_to_ray_to_vector_x = ray_to_point_x - ray_from_point_x;
    double ray_from_point_to_ray_to_vector_y = ray_to_point_y - ray_from_point_y;
    double ray_from_point_to_ray_to_vector_z = ray_to_point_z - ray_from_point_z;

    // ray_from_point_to_checking_point_vector is denoted as vector B
    double ray_from_point_to_checking_point_vector_x = checking_point_x - ray_from_point_x;
    double ray_from_point_to_checking_point_vector_y = checking_point_y - ray_from_point_y;
    double ray_from_point_to_checking_point_vector_z = checking_point_z - ray_from_point_z;

    double A_and_B_vector_cross_product_x;
    double A_and_B_vector_cross_product_y;
    double A_and_B_vector_cross_product_z;

    vector_cross_product(ray_from_point_to_ray_to_vector_x, ray_from_point_to_ray_to_vector_y,
                         ray_from_point_to_ray_to_vector_z, ray_from_point_to_checking_point_vector_x,
                         ray_from_point_to_checking_point_vector_y, ray_from_point_to_checking_point_vector_z,
                         A_and_B_vector_cross_product_x, A_and_B_vector_cross_product_y, A_and_B_vector_cross_product_z);

    // if the A_and_B_vector_cross_product is pointing out, i.e. > 0, then use the right hand rule,
    // the checking_point is on the left of the ray_from_point_to_ray_to_vector, if < 0, then
    // on the right, if = 0, then vertical
    return A_and_B_vector_cross_product_z;
}

// https://physics.stackexchange.com/questions/435512/snells-law-in-vector-form
// https://www.quora.com/How-can-you-make-a-projection-of-a-vector-onto-a-plane
// calculate the out ray vector that passes face 2
// this calculated vector is just a direction, it doesn't passes the incident_point on edge 1
// if not out at critical angle (common case), out_at_critical_angle = 1
// if out at critical angle (noncommon case), out_at_critical_angle = 2
void calculate_out_ray(double face_1_weight, double face_2_weight,
                       double source_x, double source_y, double source_z,
                       double edge_1_x_1, double edge_1_y_1, double edge_1_z_1,
                       double edge_1_x_2, double edge_1_y_2, double edge_1_z_2,
                       double face_1_x_1, double face_1_y_1, double face_1_z_1,
                       double face_1_x_2, double face_1_y_2, double face_1_z_2,
                       double face_1_x_3, double face_1_y_3, double face_1_z_3,
                       double incident_point_x, double incident_point_y, double incident_point_z,
                       double face_2_x_1, double face_2_y_1, double face_2_z_1,
                       double face_2_x_2, double face_2_y_2, double face_2_z_2,
                       double face_2_x_3, double face_2_y_3, double face_2_z_3,
                       double &out_ray_on_face_2_vector_x, double &out_ray_on_face_2_vector_y,
                       double &out_ray_on_face_2_vector_z, int &out_at_critical_angle)
{
    // we first need to calculate the normal of the edge that the ray will pass on,
    // we need to make sure that this normal is on the face that the in ray comes from

    double face_1_normal_vector_x;
    double face_1_normal_vector_y;
    double face_1_normal_vector_z;

    calculate_normal_vector_of_face(face_1_x_1, face_1_y_1, face_1_z_1,
                                    face_1_x_2, face_1_y_2, face_1_z_2,
                                    face_1_x_3, face_1_y_3, face_1_z_3,
                                    face_1_normal_vector_x, face_1_normal_vector_y,
                                    face_1_normal_vector_z);

    // face_1_edge_1_normal_vector is the normal vector vertical to the edge 1 and on face 1
    double face_1_edge_1_normal_vector_x;
    double face_1_edge_1_normal_vector_y;
    double face_1_edge_1_normal_vector_z;

    vector_cross_product(edge_1_x_2 - edge_1_x_1, edge_1_y_2 - edge_1_y_1, edge_1_z_2 - edge_1_z_1,
                         face_1_normal_vector_x, face_1_normal_vector_y, face_1_normal_vector_z,
                         face_1_edge_1_normal_vector_x, face_1_edge_1_normal_vector_y,
                         face_1_edge_1_normal_vector_z);

    // we need to check the angle between the face_1_edge_1_normal_vector and the ray
    // vector is larger than 90 degree or less than 90 degree
    // if their angle is larger than 90 degree, we need to let face_1_edge_1_normal_vector
    // be the opposite direction
    double angle_between_face_1_edge_1_normal_vector_and_ray_vector;

    double in_ray_vector_x = incident_point_x - source_x;
    double in_ray_vector_y = incident_point_y - source_y;
    double in_ray_vector_z = incident_point_z - source_z;

    vector_dot_product(face_1_edge_1_normal_vector_x, face_1_edge_1_normal_vector_y, face_1_edge_1_normal_vector_z,
                       in_ray_vector_x, in_ray_vector_y, in_ray_vector_z,
                       angle_between_face_1_edge_1_normal_vector_and_ray_vector);
    if (angle_between_face_1_edge_1_normal_vector_and_ray_vector < 0)
    {
        face_1_edge_1_normal_vector_x *= -1;
        face_1_edge_1_normal_vector_y *= -1;
        face_1_edge_1_normal_vector_z *= -1;
    }

    // using the result from https://physics.stackexchange.com/questions/435512/snells-law-in-vector-form
    // we then need to calculate the out ray that lies on the face 1

    double in_ray_vector_length = sqrt(pow(in_ray_vector_x, 2) + pow(in_ray_vector_y, 2) + pow(in_ray_vector_z, 2));

    in_ray_vector_x /= in_ray_vector_length;
    in_ray_vector_y /= in_ray_vector_length;
    in_ray_vector_z /= in_ray_vector_length;

    double face_1_edge_1_normal_vector_length = sqrt(pow(face_1_edge_1_normal_vector_x, 2) + pow(face_1_edge_1_normal_vector_y, 2) + pow(face_1_edge_1_normal_vector_z, 2));

    face_1_edge_1_normal_vector_x /= face_1_edge_1_normal_vector_length;
    face_1_edge_1_normal_vector_y /= face_1_edge_1_normal_vector_length;
    face_1_edge_1_normal_vector_z /= face_1_edge_1_normal_vector_length;

    double mu = face_1_weight / face_2_weight;
    double dot_product_in_ray_vector_and_face_1_edge_1_normal_vector;

    vector_dot_product(in_ray_vector_x, in_ray_vector_y, in_ray_vector_z,
                       face_1_edge_1_normal_vector_x, face_1_edge_1_normal_vector_y,
                       face_1_edge_1_normal_vector_z,
                       dot_product_in_ray_vector_and_face_1_edge_1_normal_vector);

    double first_term_square = 1 - pow(mu, 2) * (1 - pow(dot_product_in_ray_vector_and_face_1_edge_1_normal_vector, 2));
    double first_term;
    // the angle is not larger than the critical angle
    if (first_term_square >= 0)
    {
        first_term = sqrt(first_term_square);
        out_at_critical_angle = 1;
    }
    // the angle is larger than the critical angle
    else
    {
        first_term = 0;
        out_at_critical_angle = 2;
    }

    double out_ray_on_face_1_vector_x = first_term * face_1_edge_1_normal_vector_x +
                                        mu * (in_ray_vector_x - dot_product_in_ray_vector_and_face_1_edge_1_normal_vector * face_1_edge_1_normal_vector_x);
    double out_ray_on_face_1_vector_y = first_term * face_1_edge_1_normal_vector_y +
                                        mu * (in_ray_vector_y - dot_product_in_ray_vector_and_face_1_edge_1_normal_vector * face_1_edge_1_normal_vector_y);
    double out_ray_on_face_1_vector_z = first_term * face_1_edge_1_normal_vector_z +
                                        mu * (in_ray_vector_z - dot_product_in_ray_vector_and_face_1_edge_1_normal_vector * face_1_edge_1_normal_vector_z);

    // we next need to project this out_ray_vector onto the face_2 using
    // https://www.quora.com/How-can-you-make-a-projection-of-a-vector-onto-a-plane

    double face_2_normal_vector_x;
    double face_2_normal_vector_y;
    double face_2_normal_vector_z;

    calculate_normal_vector_of_face(face_2_x_1, face_2_y_1, face_2_z_1,
                                    face_2_x_2, face_2_y_2, face_2_z_2,
                                    face_2_x_3, face_2_y_3, face_2_z_3,
                                    face_2_normal_vector_x, face_2_normal_vector_y,
                                    face_2_normal_vector_z);

    double cross_product_out_ray_on_face_1_vector_and_face_2_normal_vector_x;
    double cross_product_out_ray_on_face_1_vector_and_face_2_normal_vector_y;
    double cross_product_out_ray_on_face_1_vector_and_face_2_normal_vector_z;

    vector_cross_product(out_ray_on_face_1_vector_x, out_ray_on_face_1_vector_y, out_ray_on_face_1_vector_z,
                         face_2_normal_vector_x, face_2_normal_vector_y, face_2_normal_vector_z,
                         cross_product_out_ray_on_face_1_vector_and_face_2_normal_vector_x,
                         cross_product_out_ray_on_face_1_vector_and_face_2_normal_vector_y,
                         cross_product_out_ray_on_face_1_vector_and_face_2_normal_vector_z);

    vector_cross_product(face_2_normal_vector_x, face_2_normal_vector_y, face_2_normal_vector_z,
                         cross_product_out_ray_on_face_1_vector_and_face_2_normal_vector_x,
                         cross_product_out_ray_on_face_1_vector_and_face_2_normal_vector_y,
                         cross_product_out_ray_on_face_1_vector_and_face_2_normal_vector_z,
                         out_ray_on_face_2_vector_x, out_ray_on_face_2_vector_y, out_ray_on_face_2_vector_z);
}

// calculate the intersection point between the out ray and the edge 2
// we seperate it from the calculate_out_ray function is because in the last edge sequence,
// there is no need to use this function, but we need calculate_out_ray
bool calculate_out_ray_intersection_point(
    double incident_point_x, double incident_point_y, double incident_point_z,
    double edge_2_x_1, double edge_2_y_1, double edge_2_z_1,
    double edge_2_x_2, double edge_2_y_2, double edge_2_z_2,
    double out_ray_on_face_2_vector_x, double out_ray_on_face_2_vector_y,
    double out_ray_on_face_2_vector_z,
    double &intersection_point_between_out_ray_and_edge_2_x,
    double &intersection_point_between_out_ray_and_edge_2_y,
    double &intersection_point_between_out_ray_and_edge_2_z)
{
    double edge_2_vector_x = edge_2_x_2 - edge_2_x_1;
    double edge_2_vector_y = edge_2_y_2 - edge_2_y_1;
    double edge_2_vector_z = edge_2_z_2 - edge_2_z_1;

    // based on these three equations, we calculate a and b
    // incident_point_x + a * out_ray_on_face_2_vector_x = edge_2_x_1 + b * edge_2_vector_x
    // incident_point_y + a * out_ray_on_face_2_vector_y = edge_2_y_1 + b * edge_2_vector_y
    // incident_point_z + a * out_ray_on_face_2_vector_z = edge_2_z_1 + b * edge_2_vector_z

    double a = (incident_point_x * edge_2_vector_y - incident_point_y * edge_2_vector_x +
                edge_2_y_1 * edge_2_vector_x - edge_2_x_1 * edge_2_vector_y) /
               (edge_2_vector_x * out_ray_on_face_2_vector_y - edge_2_vector_y * out_ray_on_face_2_vector_x);

    double b = (edge_2_y_1 * out_ray_on_face_2_vector_x - edge_2_x_1 * out_ray_on_face_2_vector_y +
                incident_point_x * out_ray_on_face_2_vector_y - incident_point_y * out_ray_on_face_2_vector_x) /
               (edge_2_vector_x * out_ray_on_face_2_vector_y - edge_2_vector_y * out_ray_on_face_2_vector_x);

    // because in the factor that we use later, we have
    // edge_2_x_2 + intersection_point_between_out_ray_and_edge_2_factor * edge_2_vector_x
    // edge_2_y_2 + intersection_point_between_out_ray_and_edge_2_factor * edge_2_vector_y
    // edge_2_z_2 + intersection_point_between_out_ray_and_edge_2_factor * edge_2_vector_z
    // intersection_point_between_out_ray_and_edge_2_factor = b - 1;
    intersection_point_between_out_ray_and_edge_2_x = round_if_close(edge_2_x_1 + b * edge_2_vector_x);
    intersection_point_between_out_ray_and_edge_2_y = round_if_close(edge_2_y_1 + b * edge_2_vector_y);
    intersection_point_between_out_ray_and_edge_2_z = round_if_close(edge_2_z_1 + b * edge_2_vector_z);

    if (!vertex_in_edge_not_endpoint(edge_2_x_1, edge_2_x_2, edge_2_y_1, edge_2_y_2, edge_2_z_1, edge_2_z_2,
                                     intersection_point_between_out_ray_and_edge_2_x,
                                     intersection_point_between_out_ray_and_edge_2_y,
                                     intersection_point_between_out_ray_and_edge_2_z))
    {
        return false;
    }
    return true;
}

// calculate the intersection point between the out ray and the edge just next to the destination point
void calculate_out_ray_intersection_point_on_edge_next_to_destination(
    double incident_point_x, double incident_point_y, double incident_point_z,
    double edge_x_1, double edge_y_1, double edge_z_1, // this is the first vertex of the edge just next to the destination point
    double edge_x_2, double edge_y_2, double edge_z_2, // this is the second vertex of the edge just next to the destination point
    double out_ray_vector_x, double out_ray_vector_y, double out_ray_vector_z,
    double &intersection_point_between_out_ray_and_edge_x,
    double &intersection_point_between_out_ray_and_edge_y,
    double &intersection_point_between_out_ray_and_edge_z)
{
    double edge_vector_x = edge_x_2 - edge_x_1;
    double edge_vector_y = edge_y_2 - edge_y_1;
    double edge_vector_z = edge_z_2 - edge_z_1;

    // based on these three equations, we calculate a and b
    // incident_point_x + a * out_ray_vector_x = edge_x_1 + b * edge_vector_x
    // incident_point_y + a * out_ray_vector_y = edge_y_1 + b * edge_vector_y
    // incident_point_z + a * out_ray_vector_z = edge_z_1 + b * edge_vector_z

    double a = (incident_point_x * edge_vector_y - incident_point_y * edge_vector_x +
                edge_y_1 * edge_vector_x - edge_x_1 * edge_vector_y) /
               (edge_vector_x * out_ray_vector_y - edge_vector_y * out_ray_vector_x);

    double b = (edge_y_1 * out_ray_vector_x - edge_x_1 * out_ray_vector_y +
                incident_point_x * out_ray_vector_y - incident_point_y * out_ray_vector_x) /
               (edge_vector_x * out_ray_vector_y - edge_vector_y * out_ray_vector_x);

    intersection_point_between_out_ray_and_edge_x = edge_x_1 + b * edge_vector_x;
    intersection_point_between_out_ray_and_edge_y = edge_y_1 + b * edge_vector_y;
    intersection_point_between_out_ray_and_edge_z = edge_z_1 + b * edge_vector_z;
}

// calculate the effective weight based on the first face and the remaining face using the given points
void calculate_effective_weight_with_given_points(
    double source_x, double source_y, double source_z,
    double first_edge_intersection_point_x,
    double first_edge_intersection_point_y,
    double first_edge_intersection_point_z,
    double edge_next_to_destination_and_ray_intersection_point_x,
    double edge_next_to_destination_and_ray_intersection_point_y,
    double edge_next_to_destination_and_ray_intersection_point_z,
    double first_face_x_1, double first_face_y_1, double first_face_z_1,
    double first_face_x_2, double first_face_y_2, double first_face_z_2,
    double first_face_x_3, double first_face_y_3, double first_face_z_3,
    double first_edge_x_1, double first_edge_y_1, double first_edge_z_1,
    double first_edge_x_2, double first_edge_y_2, double first_edge_z_2,
    double &effective_weight_ratio_of_remaining_face_over_first_face)
{
    double effective_out_ray_vector_x = edge_next_to_destination_and_ray_intersection_point_x - first_edge_intersection_point_x;
    double effective_out_ray_vector_y = edge_next_to_destination_and_ray_intersection_point_y - first_edge_intersection_point_y;
    double effective_out_ray_vector_z = edge_next_to_destination_and_ray_intersection_point_z - first_edge_intersection_point_z;

    // we then need to project the effective_out_ray_vector to the first face using
    // https://www.quora.com/How-can-you-make-a-projection-of-a-vector-onto-a-plane
    double first_face_normal_vector_x;
    double first_face_normal_vector_y;
    double first_face_normal_vector_z;

    calculate_normal_vector_of_face(first_face_x_1, first_face_y_1, first_face_z_1,
                                    first_face_x_2, first_face_y_2, first_face_z_2,
                                    first_face_x_3, first_face_y_3, first_face_z_3,
                                    first_face_normal_vector_x, first_face_normal_vector_y,
                                    first_face_normal_vector_z);

    double cross_product_effective_out_ray_vector_and_first_face_normal_vector_x;
    double cross_product_effective_out_ray_vector_and_first_face_normal_vector_y;
    double cross_product_effective_out_ray_vector_and_first_face_normal_vector_z;

    vector_cross_product(effective_out_ray_vector_x, effective_out_ray_vector_y, effective_out_ray_vector_z,
                         first_face_normal_vector_x, first_face_normal_vector_y, first_face_normal_vector_z,
                         cross_product_effective_out_ray_vector_and_first_face_normal_vector_x,
                         cross_product_effective_out_ray_vector_and_first_face_normal_vector_y,
                         cross_product_effective_out_ray_vector_and_first_face_normal_vector_z);

    double effective_out_ray_on_first_face_vector_x;
    double effective_out_ray_on_first_face_vector_y;
    double effective_out_ray_on_first_face_vector_z;

    vector_cross_product(first_face_normal_vector_x, first_face_normal_vector_y, first_face_normal_vector_z,
                         cross_product_effective_out_ray_vector_and_first_face_normal_vector_x,
                         cross_product_effective_out_ray_vector_and_first_face_normal_vector_y,
                         cross_product_effective_out_ray_vector_and_first_face_normal_vector_z,
                         effective_out_ray_on_first_face_vector_x, effective_out_ray_on_first_face_vector_y,
                         effective_out_ray_on_first_face_vector_z);

    // we then need to calculate the vector that is vertical to the first edge, and this vector should
    // lie on the first face

    // first_face_first_edge_normal_vector is the normal vector vertical to the first edge and on the first face
    double first_face_first_edge_normal_vector_x;
    double first_face_first_edge_normal_vector_y;
    double first_face_first_edge_normal_vector_z;

    vector_cross_product(first_edge_x_2 - first_edge_x_1, first_edge_y_2 - first_edge_y_1, first_edge_z_2 - first_edge_z_1,
                         first_face_normal_vector_x, first_face_normal_vector_y, first_face_normal_vector_z,
                         first_face_first_edge_normal_vector_x, first_face_first_edge_normal_vector_y,
                         first_face_first_edge_normal_vector_z);

    // we need to check the angle between the first_face_first_edge_normal_vector and the ray
    // vector is larger than 90 degree or less than 90 degree
    // if their angle is larger than 90 degree, we need to let first_face_first_edge_normal_vector
    // be the opposite direction
    double angle_between_first_face_first_edge_normal_vector_and_ray_vector;

    double in_ray_vector_x = first_edge_intersection_point_x - source_x;
    double in_ray_vector_y = first_edge_intersection_point_y - source_y;
    double in_ray_vector_z = first_edge_intersection_point_z - source_z;

    vector_dot_product(first_face_first_edge_normal_vector_x, first_face_first_edge_normal_vector_y, first_face_first_edge_normal_vector_z,
                       in_ray_vector_x, in_ray_vector_y, in_ray_vector_z,
                       angle_between_first_face_first_edge_normal_vector_and_ray_vector);
    if (angle_between_first_face_first_edge_normal_vector_and_ray_vector < 0)
    {
        first_face_first_edge_normal_vector_x *= -1;
        first_face_first_edge_normal_vector_y *= -1;
        first_face_first_edge_normal_vector_z *= -1;
    }

    // using the result from https://physics.stackexchange.com/questions/435512/snells-law-in-vector-form
    // we calculate the effective weight ratio

    // normalize in_ray_vector and effective_out_ray_on_first_face_vector and first_face_first_edge_normal_vector
    double in_ray_vector_length = sqrt(pow(in_ray_vector_x, 2) + pow(in_ray_vector_y, 2) + pow(in_ray_vector_z, 2));
    in_ray_vector_x /= in_ray_vector_length;
    in_ray_vector_y /= in_ray_vector_length;
    in_ray_vector_z /= in_ray_vector_length;

    double effective_out_ray_on_first_face_vector_length = sqrt(pow(effective_out_ray_on_first_face_vector_x, 2) + pow(effective_out_ray_on_first_face_vector_y, 2) + pow(effective_out_ray_on_first_face_vector_z, 2));
    effective_out_ray_on_first_face_vector_x /= effective_out_ray_on_first_face_vector_length;
    effective_out_ray_on_first_face_vector_y /= effective_out_ray_on_first_face_vector_length;
    effective_out_ray_on_first_face_vector_z /= effective_out_ray_on_first_face_vector_length;

    double first_face_first_edge_normal_vector_length = sqrt(pow(first_face_first_edge_normal_vector_x, 2) + pow(first_face_first_edge_normal_vector_y, 2) + pow(first_face_first_edge_normal_vector_z, 2));
    first_face_first_edge_normal_vector_x /= first_face_first_edge_normal_vector_length;
    first_face_first_edge_normal_vector_y /= first_face_first_edge_normal_vector_length;
    first_face_first_edge_normal_vector_z /= first_face_first_edge_normal_vector_length;

    double cross_product_between_first_face_first_edge_normal_vector_and_in_ray_vector_x;
    double cross_product_between_first_face_first_edge_normal_vector_and_in_ray_vector_y;
    double cross_product_between_first_face_first_edge_normal_vector_and_in_ray_vector_z;

    double cross_product_between_first_face_first_edge_normal_vector_and_effective_out_ray_on_first_face_vector_x;
    double cross_product_between_first_face_first_edge_normal_vector_and_effective_out_ray_on_first_face_vector_y;
    double cross_product_between_first_face_first_edge_normal_vector_and_effective_out_ray_on_first_face_vector_z;

    vector_cross_product(first_face_first_edge_normal_vector_x, first_face_first_edge_normal_vector_y, first_face_first_edge_normal_vector_z,
                         in_ray_vector_x, in_ray_vector_y, in_ray_vector_z,
                         cross_product_between_first_face_first_edge_normal_vector_and_in_ray_vector_x,
                         cross_product_between_first_face_first_edge_normal_vector_and_in_ray_vector_y,
                         cross_product_between_first_face_first_edge_normal_vector_and_in_ray_vector_z);

    vector_cross_product(first_face_first_edge_normal_vector_x, first_face_first_edge_normal_vector_y, first_face_first_edge_normal_vector_z,
                         effective_out_ray_on_first_face_vector_x, effective_out_ray_on_first_face_vector_y, effective_out_ray_on_first_face_vector_z,
                         cross_product_between_first_face_first_edge_normal_vector_and_effective_out_ray_on_first_face_vector_x,
                         cross_product_between_first_face_first_edge_normal_vector_and_effective_out_ray_on_first_face_vector_y,
                         cross_product_between_first_face_first_edge_normal_vector_and_effective_out_ray_on_first_face_vector_z);

    effective_weight_ratio_of_remaining_face_over_first_face = cross_product_between_first_face_first_edge_normal_vector_and_in_ray_vector_x / cross_product_between_first_face_first_edge_normal_vector_and_effective_out_ray_on_first_face_vector_x;
}

// calculate the intersection point on the first edge using effective weight
void calculate_intersection_point_on_first_edge_using_effective_weight(
    double source_x, double source_y, double source_z,
    double destination_x, double destination_y, double destination_z,
    double effective_weight_ratio_of_remaining_face_over_first_face,
    double first_face_x_1, double first_face_y_1, double first_face_z_1,
    double first_face_x_2, double first_face_y_2, double first_face_z_2,
    double first_face_x_3, double first_face_y_3, double first_face_z_3,
    double first_edge_x_1, double first_edge_y_1, double first_edge_z_1,
    double first_edge_x_2, double first_edge_y_2, double first_edge_z_2,
    double curr_edge_one_side_vertex_factor,
    double curr_edge_another_side_vertex_factor,
    double factor_error, // we move left and right of the intersection poiont factor by this factor_error and also store the corresponding factor in the result list
    std::vector<double> &first_edge_intersection_point_factor_list,
    std::vector<double> &first_edge_intersection_point_x_list,
    std::vector<double> &first_edge_intersection_point_y_list,
    std::vector<double> &first_edge_intersection_point_z_list)
{

    double edge_vector_x = first_edge_x_2 - first_edge_x_1;
    double edge_vector_y = first_edge_y_2 - first_edge_y_1;
    double edge_vector_z = first_edge_z_2 - first_edge_z_1;

    // we need to rotate the destination point on the first face,
    // so we first calculate the normal vector of the first face
    double first_face_normal_vector_x;
    double first_face_normal_vector_y;
    double first_face_normal_vector_z;

    calculate_normal_vector_of_face(first_face_x_1, first_face_y_1, first_face_z_1,
                                    first_face_x_2, first_face_y_2, first_face_z_2,
                                    first_face_x_3, first_face_y_3, first_face_z_3,
                                    first_face_normal_vector_x, first_face_normal_vector_y,
                                    first_face_normal_vector_z);

    // we then need to project the destination point on the first face
    double cross_product_destination_to_first_face_1_vector_and_first_face_normal_vector_x;
    double cross_product_destination_to_first_face_1_vector_and_first_face_normal_vector_y;
    double cross_product_destination_to_first_face_1_vector_and_first_face_normal_vector_z;

    vector_cross_product(destination_x - first_face_x_1, destination_y - first_face_y_1, destination_z - first_face_z_1,
                         first_face_normal_vector_x, first_face_normal_vector_y, first_face_normal_vector_z,
                         cross_product_destination_to_first_face_1_vector_and_first_face_normal_vector_x,
                         cross_product_destination_to_first_face_1_vector_and_first_face_normal_vector_y,
                         cross_product_destination_to_first_face_1_vector_and_first_face_normal_vector_z);

    double destination_to_first_face_1_on_first_face_vector_x;
    double destination_to_first_face_1_on_first_face_vector_y;
    double destination_to_first_face_1_on_first_face_vector_z;

    vector_cross_product(first_face_normal_vector_x, first_face_normal_vector_y, first_face_normal_vector_z,
                         cross_product_destination_to_first_face_1_vector_and_first_face_normal_vector_x,
                         cross_product_destination_to_first_face_1_vector_and_first_face_normal_vector_y,
                         cross_product_destination_to_first_face_1_vector_and_first_face_normal_vector_z,
                         destination_to_first_face_1_on_first_face_vector_x, destination_to_first_face_1_on_first_face_vector_y,
                         destination_to_first_face_1_on_first_face_vector_z);

    double first_face_normal_vector_length_square = pow(first_face_normal_vector_x, 2) + pow(first_face_normal_vector_y, 2) + pow(first_face_normal_vector_z, 2);
    destination_to_first_face_1_on_first_face_vector_x /= first_face_normal_vector_length_square;
    destination_to_first_face_1_on_first_face_vector_y /= first_face_normal_vector_length_square;
    destination_to_first_face_1_on_first_face_vector_z /= first_face_normal_vector_length_square;

    double destination_on_first_face_x = first_face_x_1 + destination_to_first_face_1_on_first_face_vector_x;
    double destination_on_first_face_y = first_face_y_1 + destination_to_first_face_1_on_first_face_vector_y;
    double destination_on_first_face_z = first_face_z_1 + destination_to_first_face_1_on_first_face_vector_z;

    // we regard the first edge's two vertices and the destination point as the second (fake) face
    double second_face_normal_vector_x;
    double second_face_normal_vector_y;
    double second_face_normal_vector_z;

    calculate_normal_vector_of_face(first_edge_x_1, first_edge_y_1, first_edge_z_1,
                                    first_edge_x_2, first_edge_y_2, first_edge_z_2,
                                    destination_x, destination_y, destination_z,
                                    second_face_normal_vector_x, second_face_normal_vector_y,
                                    second_face_normal_vector_z);

    // second_face_first_edge_normal_vector is the normal vector vertical to the first edge and on the second face
    double second_face_first_edge_normal_vector_x;
    double second_face_first_edge_normal_vector_y;
    double second_face_first_edge_normal_vector_z;

    vector_cross_product(edge_vector_x, edge_vector_y, edge_vector_z,
                         second_face_normal_vector_x, second_face_normal_vector_y, second_face_normal_vector_z,
                         second_face_first_edge_normal_vector_x, second_face_first_edge_normal_vector_y,
                         second_face_first_edge_normal_vector_z);

    // we then calculate the intersection point between the second_face_first_edge_normal_vector and the edge_vector
    // the line of this intersection point and destination is on second face and vertical to the edge_vector

    // based on these three equations, we calculate a_1 and b_1
    // destination_x + a_1 * second_face_first_edge_normal_vector_x = first_edge_x_2 + b_1 * edge_vector_x
    // destination_y + a_1 * second_face_first_edge_normal_vector_y = first_edge_y_2 + b_1 * edge_vector_y
    // destination_z + a_1 * second_face_first_edge_normal_vector_z = first_edge_z_2 + b_1 * edge_vector_z

    double a_1 = (destination_x * edge_vector_y - destination_y * edge_vector_x +
                  first_edge_y_2 * edge_vector_x - first_edge_x_2 * edge_vector_y) /
                 (edge_vector_x * second_face_first_edge_normal_vector_y - edge_vector_y * second_face_first_edge_normal_vector_x);

    double b_1 = (first_edge_y_2 * second_face_first_edge_normal_vector_x - first_edge_x_2 * second_face_first_edge_normal_vector_y +
                  destination_x * second_face_first_edge_normal_vector_y - destination_y * second_face_first_edge_normal_vector_x) /
                 (edge_vector_x * second_face_first_edge_normal_vector_y - edge_vector_y * second_face_first_edge_normal_vector_x);

    double foot_point_of_destination_on_first_edge_x = first_edge_x_2 + b_1 * edge_vector_x;
    double foot_point_of_destination_on_first_edge_y = first_edge_y_2 + b_1 * edge_vector_y;
    double foot_point_of_destination_on_first_edge_z = first_edge_z_2 + b_1 * edge_vector_z;

    double foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector_x = destination_on_first_face_x - foot_point_of_destination_on_first_edge_x;
    double foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector_y = destination_on_first_face_y - foot_point_of_destination_on_first_edge_y;
    double foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector_z = destination_on_first_face_z - foot_point_of_destination_on_first_edge_z;

    // the rotation endpoint of foor_point_of_destination_on_first_edge_to_destination_vector to foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector is
    // foot_point_of_destination_on_first_edge_x + a_2 * foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector_x
    // foot_point_of_destination_on_first_edge_x + a_2 * foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector_y
    // foot_point_of_destination_on_first_edge_x + a_2 * foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector_z

    // we need to let the length of foor_point_of_destination_on_first_edge_to_destination_on_first_face to be the same as foor_point_of_destination_on_first_edge_to_destination

    // pow(a_2, 2) * (pow(foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector_x, 2) + pow(foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector_y, 2) + pow(foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector_z, 2))
    // = pow((destination_x - foot_point_of_destination_on_first_edge_x), 2) + pow((destination_y - foot_point_of_destination_on_first_edge_y), 2) + pow((destination_z - foot_point_of_destination_on_first_edge_z), 2)

    double a_2 = sqrt((pow((destination_x - foot_point_of_destination_on_first_edge_x), 2) + pow((destination_y - foot_point_of_destination_on_first_edge_y), 2) + pow((destination_z - foot_point_of_destination_on_first_edge_z), 2)) /
                      (pow(foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector_x, 2) + pow(foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector_y, 2) + pow(foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector_z, 2)));

    destination_on_first_face_x = foot_point_of_destination_on_first_edge_x + a_2 * foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector_x;
    destination_on_first_face_y = foot_point_of_destination_on_first_edge_y + a_2 * foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector_y;
    destination_on_first_face_z = foot_point_of_destination_on_first_edge_z + a_2 * foor_point_of_destination_on_first_edge_to_destination_on_first_face_vector_z;

    // the calculated effective weight intersection point could be expressed as
    // first_edge_x_2 + a * edge_vector_x;
    // first_edge_y_2 + a * edge_vector_y;
    // first_edge_z_2 + a * edge_vector_z;

    // we then need to calculate the vector that is vertical to the first edge, and this vector should
    // lie on the first face

    // first_face_first_edge_normal_vector is the normal vector vertical to the first edge and on the first face
    double first_face_first_edge_normal_vector_x;
    double first_face_first_edge_normal_vector_y;
    double first_face_first_edge_normal_vector_z;

    vector_cross_product(edge_vector_x, edge_vector_y, edge_vector_z,
                         first_face_normal_vector_x, first_face_normal_vector_y, first_face_normal_vector_z,
                         first_face_first_edge_normal_vector_x, first_face_first_edge_normal_vector_y,
                         first_face_first_edge_normal_vector_z);

    // normalize first_face_first_edge_normal_vector
    double first_face_first_edge_normal_vector_length = sqrt(pow(first_face_first_edge_normal_vector_x, 2) + pow(first_face_first_edge_normal_vector_y, 2) + pow(first_face_first_edge_normal_vector_z, 2));
    first_face_first_edge_normal_vector_x /= first_face_first_edge_normal_vector_length;
    first_face_first_edge_normal_vector_y /= first_face_first_edge_normal_vector_length;
    first_face_first_edge_normal_vector_z /= first_face_first_edge_normal_vector_length;

    // using the result from https://physics.stackexchange.com/questions/435512/snells-law-in-vector-form
    // we calculate the intersection point on the current edge using effective weight

    // in order to save the space, we rename the variable as follows
    double mu = effective_weight_ratio_of_remaining_face_over_first_face;
    double sx = source_x;
    double sy = source_y;
    double sz = source_z;
    double dx = destination_on_first_face_x;
    double dy = destination_on_first_face_y;
    double dz = destination_on_first_face_z;
    double x1 = first_edge_x_1;
    double y1 = first_edge_y_1;
    double z1 = first_edge_z_1;
    double x2 = first_edge_x_2;
    double y2 = first_edge_y_2;
    double z2 = first_edge_z_2;
    double nx = first_face_first_edge_normal_vector_x;
    double ny = first_face_first_edge_normal_vector_y;
    double nz = first_face_first_edge_normal_vector_z;
    double ex = edge_vector_x;
    double ey = edge_vector_y;
    double ez = edge_vector_z;
    double a1 = 0;
    double a2 = 0;
    double a3 = 0;
    double a4 = 0;

    double x3 = dx - x2;
    double y3 = dy - y2;
    double z3 = dz - z2;
    double x4 = x2 - sx;
    double y4 = y2 - sy;
    double z4 = z2 - sz;

    double p1 = ny * z3 - nz * y3;
    double q1 = nz * ey - ny * ez;
    double p2 = ny * z4 - nz * y4;
    double q2 = ny * ez - nz * ey;

    double k0 = pow(ex, 2) + pow(ey, 2) + pow(ez, 2);
    double k1 = -2 * x3 * ex - 2 * y3 * ey - 2 * z3 * ez;
    double k2 = pow(x3, 2) + pow(y3, 2) + pow(z3, 2);
    double k3 = 2 * x4 * ex + 2 * y4 * ey + 2 * z4 * ez;
    double k4 = pow(x4, 2) + pow(y4, 2) + pow(z4, 2);

    double m1 = pow(q1, 2);
    double m2 = 2 * p1 * q1;
    double m3 = pow(p1, 2);
    double n1 = pow(q2, 2);
    double n2 = 2 * p2 * q2;
    double n3 = pow(p2, 2);

    double v1 = pow(mu, 2) * m1 * k0 - n1 * k0;
    double v2 = pow(mu, 2) * m1 * k3 + pow(mu, 2) * m2 * k0 - n1 * k1 - n2 * k0;
    double v3 = pow(mu, 2) * m1 * k4 + pow(mu, 2) * m2 * k3 + pow(mu, 2) * m3 * k0 - n1 * k2 - n2 * k1 - n3 * k0;
    double v4 = pow(mu, 2) * m2 * k4 + pow(mu, 2) * m3 * k3 - n2 * k2 - n3 * k1;
    double v5 = pow(mu, 2) * m3 * k4 - n3 * k2;

    // we then use the following function to calculate the solution for quartic equation
    // https://github.com/sasamil/Quartic
    std::complex<double> *solutions = solve_quartic(v2 / v1, v3 / v1, v4 / v1, v5 / v1);

    if (solutions[0].imag() == 0.0)
    {
        a1 = solutions[0].real();
    }
    if (solutions[1].imag() == 0.0)
    {
        a2 = solutions[1].real();
    }
    if (solutions[2].imag() == 0.0)
    {
        a3 = solutions[2].real();
    }
    if (solutions[3].imag() == 0.0)
    {
        a4 = solutions[3].real();
    }

    bool at_least_one_value_in_factor_list = false;

    // store the calculated factor
    if ((curr_edge_one_side_vertex_factor > curr_edge_another_side_vertex_factor && a1 < curr_edge_one_side_vertex_factor && a1 > curr_edge_another_side_vertex_factor) ||
        (curr_edge_one_side_vertex_factor < curr_edge_another_side_vertex_factor && a1 > curr_edge_one_side_vertex_factor && a1 < curr_edge_another_side_vertex_factor))
    {
        first_edge_intersection_point_factor_list.push_back(a1);
        at_least_one_value_in_factor_list = true;
    }
    if ((curr_edge_one_side_vertex_factor > curr_edge_another_side_vertex_factor && a2 < curr_edge_one_side_vertex_factor && a2 > curr_edge_another_side_vertex_factor) ||
        (curr_edge_one_side_vertex_factor < curr_edge_another_side_vertex_factor && a2 > curr_edge_one_side_vertex_factor && a2 < curr_edge_another_side_vertex_factor))
    {
        first_edge_intersection_point_factor_list.push_back(a2);
        at_least_one_value_in_factor_list = true;
    }
    if ((curr_edge_one_side_vertex_factor > curr_edge_another_side_vertex_factor && a3 < curr_edge_one_side_vertex_factor && a3 > curr_edge_another_side_vertex_factor) ||
        (curr_edge_one_side_vertex_factor < curr_edge_another_side_vertex_factor && a3 > curr_edge_one_side_vertex_factor && a3 < curr_edge_another_side_vertex_factor))
    {
        first_edge_intersection_point_factor_list.push_back(a3);
        at_least_one_value_in_factor_list = true;
    }
    if ((curr_edge_one_side_vertex_factor > curr_edge_another_side_vertex_factor && a4 < curr_edge_one_side_vertex_factor && a4 > curr_edge_another_side_vertex_factor) ||
        (curr_edge_one_side_vertex_factor < curr_edge_another_side_vertex_factor && a4 > curr_edge_one_side_vertex_factor && a4 < curr_edge_another_side_vertex_factor))
    {
        first_edge_intersection_point_factor_list.push_back(a4);
        at_least_one_value_in_factor_list = true;
    }

    // store the calculated factor plus the factor_error
    if ((a1 != 0 && curr_edge_one_side_vertex_factor > curr_edge_another_side_vertex_factor && a1 + factor_error < curr_edge_one_side_vertex_factor && a1 + factor_error > curr_edge_another_side_vertex_factor) ||
        (a1 != 0 && curr_edge_one_side_vertex_factor < curr_edge_another_side_vertex_factor && a1 + factor_error > curr_edge_one_side_vertex_factor && a1 + factor_error < curr_edge_another_side_vertex_factor))
    {
        first_edge_intersection_point_factor_list.push_back(a1 + factor_error);
        at_least_one_value_in_factor_list = true;
    }
    if ((a2 != 0 && curr_edge_one_side_vertex_factor > curr_edge_another_side_vertex_factor && a2 + factor_error < curr_edge_one_side_vertex_factor && a2 + factor_error > curr_edge_another_side_vertex_factor) ||
        (a2 != 0 && curr_edge_one_side_vertex_factor < curr_edge_another_side_vertex_factor && a2 + factor_error > curr_edge_one_side_vertex_factor && a2 + factor_error < curr_edge_another_side_vertex_factor))
    {
        first_edge_intersection_point_factor_list.push_back(a2 + factor_error);
        at_least_one_value_in_factor_list = true;
    }
    if ((a3 != 0 && curr_edge_one_side_vertex_factor > curr_edge_another_side_vertex_factor && a3 + factor_error < curr_edge_one_side_vertex_factor && a3 + factor_error > curr_edge_another_side_vertex_factor) ||
        (a3 != 0 && curr_edge_one_side_vertex_factor < curr_edge_another_side_vertex_factor && a3 + factor_error > curr_edge_one_side_vertex_factor && a3 + factor_error < curr_edge_another_side_vertex_factor))
    {
        first_edge_intersection_point_factor_list.push_back(a3 + factor_error);
        at_least_one_value_in_factor_list = true;
    }
    if ((a4 != 0 && curr_edge_one_side_vertex_factor > curr_edge_another_side_vertex_factor && a4 + factor_error < curr_edge_one_side_vertex_factor && a4 + factor_error > curr_edge_another_side_vertex_factor) ||
        (a4 != 0 && curr_edge_one_side_vertex_factor < curr_edge_another_side_vertex_factor && a4 + factor_error > curr_edge_one_side_vertex_factor && a4 + factor_error < curr_edge_another_side_vertex_factor))
    {
        first_edge_intersection_point_factor_list.push_back(a4 + factor_error);
        at_least_one_value_in_factor_list = true;
    }

    // store the calculated factor minus the factor_error
    if ((a1 != 0 && curr_edge_one_side_vertex_factor > curr_edge_another_side_vertex_factor && a1 - factor_error < curr_edge_one_side_vertex_factor && a1 - factor_error > curr_edge_another_side_vertex_factor) ||
        (a1 != 0 && curr_edge_one_side_vertex_factor < curr_edge_another_side_vertex_factor && a1 - factor_error > curr_edge_one_side_vertex_factor && a1 - factor_error < curr_edge_another_side_vertex_factor))
    {
        first_edge_intersection_point_factor_list.push_back(a1 - factor_error);
        at_least_one_value_in_factor_list = true;
    }
    if ((a2 != 0 && curr_edge_one_side_vertex_factor > curr_edge_another_side_vertex_factor && a2 - factor_error < curr_edge_one_side_vertex_factor && a2 - factor_error > curr_edge_another_side_vertex_factor) ||
        (a2 != 0 && curr_edge_one_side_vertex_factor < curr_edge_another_side_vertex_factor && a2 - factor_error > curr_edge_one_side_vertex_factor && a2 - factor_error < curr_edge_another_side_vertex_factor))
    {
        first_edge_intersection_point_factor_list.push_back(a2 - factor_error);
        at_least_one_value_in_factor_list = true;
    }
    if ((a3 != 0 && curr_edge_one_side_vertex_factor > curr_edge_another_side_vertex_factor && a3 - factor_error < curr_edge_one_side_vertex_factor && a3 - factor_error > curr_edge_another_side_vertex_factor) ||
        (a3 != 0 && curr_edge_one_side_vertex_factor < curr_edge_another_side_vertex_factor && a3 - factor_error > curr_edge_one_side_vertex_factor && a3 - factor_error < curr_edge_another_side_vertex_factor))
    {
        first_edge_intersection_point_factor_list.push_back(a3 - factor_error);
        at_least_one_value_in_factor_list = true;
    }
    if ((a4 != 0 && curr_edge_one_side_vertex_factor > curr_edge_another_side_vertex_factor && a4 - factor_error < curr_edge_one_side_vertex_factor && a4 - factor_error > curr_edge_another_side_vertex_factor) ||
        (a4 != 0 && curr_edge_one_side_vertex_factor < curr_edge_another_side_vertex_factor && a4 - factor_error > curr_edge_one_side_vertex_factor && a4 - factor_error < curr_edge_another_side_vertex_factor))
    {
        first_edge_intersection_point_factor_list.push_back(a4 - factor_error);
        at_least_one_value_in_factor_list = true;
    }

    // if no value in factor list
    if (!at_least_one_value_in_factor_list)
    {
        first_edge_intersection_point_factor_list.push_back((curr_edge_one_side_vertex_factor + curr_edge_another_side_vertex_factor) / 2);
    }

    assert(curr_edge_one_side_vertex_factor != curr_edge_another_side_vertex_factor);
    if (curr_edge_one_side_vertex_factor > curr_edge_another_side_vertex_factor)
    {
        std::sort(std::begin(first_edge_intersection_point_factor_list), std::end(first_edge_intersection_point_factor_list), std::greater<double>{});
    }
    else if (curr_edge_one_side_vertex_factor < curr_edge_another_side_vertex_factor)
    {
        std::sort(std::begin(first_edge_intersection_point_factor_list), std::end(first_edge_intersection_point_factor_list), std::less<double>{});
    }

    for (int i = 0; i < first_edge_intersection_point_factor_list.size(); ++i)
    {
        first_edge_intersection_point_x_list.push_back(first_edge_x_2 + first_edge_intersection_point_factor_list[i] * edge_vector_x);
        first_edge_intersection_point_y_list.push_back(first_edge_y_2 + first_edge_intersection_point_factor_list[i] * edge_vector_y);
        first_edge_intersection_point_z_list.push_back(first_edge_z_2 + first_edge_intersection_point_factor_list[i] * edge_vector_z);
    }
}

// this function is different from the snell_law() function
// in this function, we will give the position of the edge intersection point on the edge
bool snell_law_given_edge_intersection_point(
    geodesic::Mesh *mesh,
    std::vector<geodesic::Edge> &edge_sequence,
    std::vector<geodesic::Face> &face_sequence,
    double source_x, double source_y, double source_z,
    geodesic::SurfacePoint &destination,
    double first_intersection_point_x,
    double first_intersection_point_y,
    double first_intersection_point_z,
    std::vector<geodesic::SurfacePoint> &snell_law_path,
    int &exit_edge_index, // the snell law path pass edge_sequence[exit_edge_index], but not pass edge_sequence[exit_edge_index - 1], i.e. the next edge
    int &return_out_at_critical_angle,
    int segment_path_length)
{
    geodesic::Edge *first_edge;
    first_edge = &edge_sequence[edge_sequence.size() - 1]; // calculate the first edge that is just opposite of the source

    assert(face_sequence.size() - 1 == edge_sequence.size());
    assert(edge_sequence.size() >= 1);

    snell_law_path.clear();

    geodesic::Vertex *snell_law_path_point = new geodesic::Vertex();

    snell_law_path_point->set(source_x, source_y, source_z);
    snell_law_path.push_back(geodesic::SurfacePoint(snell_law_path_point));

    double incident_point_x = first_intersection_point_x;
    double incident_point_y = first_intersection_point_y;
    double incident_point_z = first_intersection_point_z;

    snell_law_path_point->set(incident_point_x, incident_point_y, incident_point_z);
    snell_law_path.push_back(geodesic::SurfacePoint(snell_law_path_point));

    double out_ray_on_face_2_vector_x;
    double out_ray_on_face_2_vector_y;
    double out_ray_on_face_2_vector_z;

    double intersection_point_between_out_ray_and_edge_2_factor;
    double intersection_point_between_out_ray_and_edge_2_x;
    double intersection_point_between_out_ray_and_edge_2_y;
    double intersection_point_between_out_ray_and_edge_2_z;

    bool snell_law_path_inside_edge_region = true;
    int out_at_critical_angle;

    exit_edge_index = 0;

    for (int i = edge_sequence.size() - 1; i >= 0; --i)
    {
        int iteration = 0;
        if (mesh->vertices().size() < 4000)
        {
            iteration = 100;
        }
        else if (mesh->vertices().size() < 100000)
        {
            iteration = 500;
        }
        else
        {
            iteration = 10000;
        }
        for (int j = 0; j < iteration * segment_path_length; ++j)
        {
            double face_1_weight = face_sequence[i + 1].weight();
            double face_2_weight = face_sequence[i].weight();

            calculate_out_ray(face_1_weight, face_2_weight, source_x, source_y, source_z,
                              edge_sequence[i].adjacent_vertices()[0]->x(),
                              edge_sequence[i].adjacent_vertices()[0]->y(),
                              edge_sequence[i].adjacent_vertices()[0]->z(),
                              edge_sequence[i].adjacent_vertices()[1]->x(),
                              edge_sequence[i].adjacent_vertices()[1]->y(),
                              edge_sequence[i].adjacent_vertices()[1]->z(),
                              face_sequence[i + 1].adjacent_vertices()[0]->x(),
                              face_sequence[i + 1].adjacent_vertices()[0]->y(),
                              face_sequence[i + 1].adjacent_vertices()[0]->z(),
                              face_sequence[i + 1].adjacent_vertices()[1]->x(),
                              face_sequence[i + 1].adjacent_vertices()[1]->y(),
                              face_sequence[i + 1].adjacent_vertices()[1]->z(),
                              face_sequence[i + 1].adjacent_vertices()[2]->x(),
                              face_sequence[i + 1].adjacent_vertices()[2]->y(),
                              face_sequence[i + 1].adjacent_vertices()[2]->z(),
                              incident_point_x, incident_point_y, incident_point_z,
                              face_sequence[i].adjacent_vertices()[0]->x(),
                              face_sequence[i].adjacent_vertices()[0]->y(),
                              face_sequence[i].adjacent_vertices()[0]->z(),
                              face_sequence[i].adjacent_vertices()[1]->x(),
                              face_sequence[i].adjacent_vertices()[1]->y(),
                              face_sequence[i].adjacent_vertices()[1]->z(),
                              face_sequence[i].adjacent_vertices()[2]->x(),
                              face_sequence[i].adjacent_vertices()[2]->y(),
                              face_sequence[i].adjacent_vertices()[2]->z(),
                              out_ray_on_face_2_vector_x, out_ray_on_face_2_vector_y,
                              out_ray_on_face_2_vector_z, out_at_critical_angle);
        }
        return_out_at_critical_angle = out_at_critical_angle;

        if (i == 0)
        {
            double out_ray_on_final_face_x = incident_point_x + out_ray_on_face_2_vector_x;
            double out_ray_on_final_face_y = incident_point_y + out_ray_on_face_2_vector_y;
            double out_ray_on_final_face_z = incident_point_z + out_ray_on_face_2_vector_z;

            snell_law_path_point->set(out_ray_on_final_face_x, out_ray_on_final_face_y, out_ray_on_final_face_z);
            snell_law_path.push_back(geodesic::SurfacePoint(snell_law_path_point));
            break;
        }

        snell_law_path_inside_edge_region =
            calculate_out_ray_intersection_point(incident_point_x, incident_point_y, incident_point_z,
                                                 edge_sequence[i - 1].adjacent_vertices()[0]->x(),
                                                 edge_sequence[i - 1].adjacent_vertices()[0]->y(),
                                                 edge_sequence[i - 1].adjacent_vertices()[0]->z(),
                                                 edge_sequence[i - 1].adjacent_vertices()[1]->x(),
                                                 edge_sequence[i - 1].adjacent_vertices()[1]->y(),
                                                 edge_sequence[i - 1].adjacent_vertices()[1]->z(),
                                                 out_ray_on_face_2_vector_x, out_ray_on_face_2_vector_y,
                                                 out_ray_on_face_2_vector_z,
                                                 intersection_point_between_out_ray_and_edge_2_x,
                                                 intersection_point_between_out_ray_and_edge_2_y,
                                                 intersection_point_between_out_ray_and_edge_2_z);

        if (!snell_law_path_inside_edge_region)
        {
            exit_edge_index = i;

            double out_ray_on_final_face_x = incident_point_x + out_ray_on_face_2_vector_x;
            double out_ray_on_final_face_y = incident_point_y + out_ray_on_face_2_vector_y;
            double out_ray_on_final_face_z = incident_point_z + out_ray_on_face_2_vector_z;

            snell_law_path_point->set(out_ray_on_final_face_x, out_ray_on_final_face_y, out_ray_on_final_face_z);
            snell_law_path.push_back(geodesic::SurfacePoint(snell_law_path_point));
            return false;
        }

        snell_law_path_point->set(intersection_point_between_out_ray_and_edge_2_x,
                                  intersection_point_between_out_ray_and_edge_2_y,
                                  intersection_point_between_out_ray_and_edge_2_z);
        snell_law_path.push_back(geodesic::SurfacePoint(snell_law_path_point));

        source_x = incident_point_x;
        source_y = incident_point_y;
        source_z = incident_point_z;

        incident_point_x = intersection_point_between_out_ray_and_edge_2_x;
        incident_point_y = intersection_point_between_out_ray_and_edge_2_y;
        incident_point_z = intersection_point_between_out_ray_and_edge_2_z;
    }
    return true;
}

void baseline_binary_search_one_time_of_one_edge(geodesic::Mesh *mesh,
                                                 std::vector<geodesic::Edge> &edge_sequence,
                                                 std::vector<geodesic::Face> &face_sequence,
                                                 //  geodesic::SurfacePoint &source,
                                                 double source_x, double source_y, double source_z,
                                                 geodesic::SurfacePoint &destination,
                                                 std::vector<geodesic::SurfacePoint> &snell_law_path,
                                                 double delta,
                                                 double &next_source_x,
                                                 double &next_source_y,
                                                 double &next_source_z,
                                                 int &binary_search_of_snell_law_path_count_return,
                                                 double &snell_law_memory_size,
                                                 int segment_path_length)
{
    geodesic::Edge *first_edge;
    first_edge = &edge_sequence[edge_sequence.size() - 1]; // calculate the first edge that is just opposite of the source

    double curr_edge_one_side_vertex_x;
    double curr_edge_one_side_vertex_y;
    double curr_edge_one_side_vertex_z;
    double curr_edge_another_side_vertex_x;
    double curr_edge_another_side_vertex_y;
    double curr_edge_another_side_vertex_z;

    double curr_edge_one_side_vertex_factor = -1;
    double curr_edge_another_side_vertex_factor = 0;

    bool snell_law_path_inside_edge_region = true;

    curr_edge_one_side_vertex_x = first_edge->adjacent_vertices()[0]->x();
    curr_edge_one_side_vertex_y = first_edge->adjacent_vertices()[0]->y();
    curr_edge_one_side_vertex_z = first_edge->adjacent_vertices()[0]->z();
    curr_edge_another_side_vertex_x = first_edge->adjacent_vertices()[1]->x();
    curr_edge_another_side_vertex_y = first_edge->adjacent_vertices()[1]->y();
    curr_edge_another_side_vertex_z = first_edge->adjacent_vertices()[1]->z();

    double curr_edge_vector_x = curr_edge_another_side_vertex_x - curr_edge_one_side_vertex_x;
    double curr_edge_vector_y = curr_edge_another_side_vertex_y - curr_edge_one_side_vertex_y;
    double curr_edge_vector_z = curr_edge_another_side_vertex_z - curr_edge_one_side_vertex_z;

    double curr_edge_start_vertex_x = first_edge->adjacent_vertices()[1]->x();
    double curr_edge_start_vertex_y = first_edge->adjacent_vertices()[1]->y();
    double curr_edge_start_vertex_z = first_edge->adjacent_vertices()[1]->z();

    double curr_edge_intersection_point_x;
    double curr_edge_intersection_point_y;
    double curr_edge_intersection_point_z;

    double curr_edge_intersection_point_factor = 0;

    int exit_edge_index;
    int out_at_critical_angle;

    int binary_search_of_snell_law_path_count = 0;

    // calcualte longest edge
    double max_edge_length = -1e100;
    std::vector<double> edge_length_list;
    std::vector<double> edge_length_without_outliers_list;
    for (unsigned i = 0; i < mesh->edges().size(); ++i)
    {
        geodesic::Edge &e = mesh->edges()[i];
        double edge_length = e.length();
        edge_length_list.push_back(edge_length);
    }
    remove_outliers(edge_length_list, edge_length_without_outliers_list);
    for (int i = 0; i < edge_length_without_outliers_list.size(); i++)
    {
        max_edge_length = std::max(max_edge_length, edge_length_without_outliers_list[i]);
    }

    int max_iteration = floor(log(max_edge_length / delta) / log(2));

    while (distance_xyz(curr_edge_one_side_vertex_x, curr_edge_one_side_vertex_y,
                        curr_edge_one_side_vertex_z, curr_edge_another_side_vertex_x,
                        curr_edge_another_side_vertex_y, curr_edge_another_side_vertex_z) >= delta)
    {
        if (binary_search_of_snell_law_path_count > max_iteration)
        {
            break;
        }

        binary_search_of_snell_law_path_count++;

        snell_law_path.clear();

        curr_edge_intersection_point_factor = (curr_edge_one_side_vertex_factor + curr_edge_another_side_vertex_factor) / 2;

        curr_edge_intersection_point_x = curr_edge_intersection_point_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
        curr_edge_intersection_point_y = curr_edge_intersection_point_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
        curr_edge_intersection_point_z = curr_edge_intersection_point_factor * curr_edge_vector_z + curr_edge_start_vertex_z;

        auto start = std::chrono::high_resolution_clock::now();

        snell_law_path_inside_edge_region = snell_law_given_edge_intersection_point(
            mesh, edge_sequence, face_sequence, source_x, source_y, source_z, destination, curr_edge_intersection_point_x, curr_edge_intersection_point_y, curr_edge_intersection_point_z, snell_law_path, exit_edge_index, out_at_critical_angle, segment_path_length);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

        snell_law_memory_size += snell_law_path.size() * sizeof(geodesic::SurfacePoint);

        // check the left or right side of the current edge
        // if the v0 of the current edge is on the left of the ray
        bool v0_of_curr_edge_is_on_left_of_ray = false;
        if (point_on_left_of_ray(
                curr_edge_intersection_point_x, curr_edge_intersection_point_y, curr_edge_intersection_point_z,
                source_x, source_y, source_z,
                curr_edge_one_side_vertex_x, curr_edge_one_side_vertex_y, curr_edge_one_side_vertex_z))
        {
            v0_of_curr_edge_is_on_left_of_ray = true;
        }

        // if the ray pass the whole edge sequence
        if (snell_law_path_inside_edge_region)
        {
            double last_ray_point_x = snell_law_path[snell_law_path.size() - 1].x();
            double last_ray_point_y = snell_law_path[snell_law_path.size() - 1].y();
            double last_ray_point_z = snell_law_path[snell_law_path.size() - 1].z();
            double last_edge_intersection_point_x = snell_law_path[snell_law_path.size() - 2].x();
            double last_edge_intersection_point_y = snell_law_path[snell_law_path.size() - 2].y();
            double last_edge_intersection_point_z = snell_law_path[snell_law_path.size() - 2].z();

            // if the ray pass the destination point
            if (vertex_on_edge(last_ray_point_x, last_ray_point_y, last_ray_point_z,
                               last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                               destination.x(), destination.y(), destination.z()))
            {
                return;
            }
            // if the destination point is on the left of the ray
            else if (point_on_left_of_ray(
                         last_ray_point_x, last_ray_point_y, last_ray_point_z,
                         last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                         destination.x(), destination.y(), destination.z()))
            {
                if (v0_of_curr_edge_is_on_left_of_ray)
                {
                    curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                    curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                    curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                    curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                }
                else
                {
                    curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                    curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                    curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                    curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                }
            }
            // if the destination point is on the right of the ray
            else if (!point_on_left_of_ray(
                         last_ray_point_x, last_ray_point_y, last_ray_point_z,
                         last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                         destination.x(), destination.y(), destination.z()))
            {
                if (!v0_of_curr_edge_is_on_left_of_ray)
                {
                    curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                    curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                    curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                    curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                }
                else
                {
                    curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                    curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                    curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                    curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                }
            }
        }
        // if the ray doesn't pass the whole edge sequence and leave at edge_sequence[exit_edge_index - 1]
        // note that this ray does pass edge_sequence[exit_edge_index]
        else
        {
            // not out at critical angle (common case)
            if (out_at_critical_angle == 1)
            {
                double last_ray_point_x = snell_law_path[snell_law_path.size() - 1].x();
                double last_ray_point_y = snell_law_path[snell_law_path.size() - 1].y();
                double last_ray_point_z = snell_law_path[snell_law_path.size() - 1].z();
                double last_edge_intersection_point_x = snell_law_path[snell_law_path.size() - 2].x();
                double last_edge_intersection_point_y = snell_law_path[snell_law_path.size() - 2].y();
                double last_edge_intersection_point_z = snell_law_path[snell_law_path.size() - 2].z();

                double e_v0_ray_position = point_on_left_right_vertical_of_ray(
                    last_ray_point_x, last_ray_point_y, last_ray_point_z,
                    last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                    edge_sequence[exit_edge_index - 1].adjacent_vertices()[0]->x(),
                    edge_sequence[exit_edge_index - 1].adjacent_vertices()[0]->y(),
                    edge_sequence[exit_edge_index - 1].adjacent_vertices()[0]->z());

                double e_v1_ray_position = point_on_left_right_vertical_of_ray(
                    last_ray_point_x, last_ray_point_y, last_ray_point_z,
                    last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                    edge_sequence[exit_edge_index - 1].adjacent_vertices()[1]->x(),
                    edge_sequence[exit_edge_index - 1].adjacent_vertices()[1]->y(),
                    edge_sequence[exit_edge_index - 1].adjacent_vertices()[1]->z());

                // if the two points on the edge_sequence[exit_edge_index - 1] is on the left of the ray
                if ((e_v0_ray_position > 0 && e_v1_ray_position > 0) ||
                    (e_v0_ray_position == 0 && e_v1_ray_position > 0) ||
                    (e_v0_ray_position > 0 && e_v1_ray_position == 0) ||
                    (e_v0_ray_position < 0 && e_v0_ray_position > -1e-3 && e_v1_ray_position > 0) ||
                    (e_v0_ray_position > 0 && e_v1_ray_position < 0 && e_v1_ray_position > -1e-3))
                {
                    if (v0_of_curr_edge_is_on_left_of_ray)
                    {
                        curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                        curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                        curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                        curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                    }
                    else
                    {
                        curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                        curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                        curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                        curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                    }
                }
                // if the two points on the edge_sequence[exit_edge_index - 1] is on the right of the ray
                else if ((e_v0_ray_position < 0 && e_v1_ray_position < 0) ||
                         (e_v0_ray_position == 0 && e_v1_ray_position < 0) ||
                         (e_v0_ray_position < 0 && e_v1_ray_position == 0) ||
                         (e_v0_ray_position > 0 && e_v0_ray_position < 1e-3 && e_v1_ray_position < 0) ||
                         (e_v0_ray_position < 0 && e_v1_ray_position > 0 && e_v1_ray_position < 1e-3))
                {
                    if (!v0_of_curr_edge_is_on_left_of_ray)
                    {
                        curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                        curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                        curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                        curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                    }
                    else
                    {
                        curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                        curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                        curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                        curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                    }
                }
            }
            // out at critical angle
            else if (out_at_critical_angle == 2)
            {
                double last_ray_point_x = snell_law_path[snell_law_path.size() - 1].x();
                double last_ray_point_y = snell_law_path[snell_law_path.size() - 1].y();
                double last_ray_point_z = snell_law_path[snell_law_path.size() - 1].z();
                double last_edge_intersection_point_x = snell_law_path[snell_law_path.size() - 2].x();
                double last_edge_intersection_point_y = snell_law_path[snell_law_path.size() - 2].y();
                double last_edge_intersection_point_z = snell_law_path[snell_law_path.size() - 2].z();
                double last_second_edge_intersection_point_x = snell_law_path[snell_law_path.size() - 3].x();
                double last_second_edge_intersection_point_y = snell_law_path[snell_law_path.size() - 3].y();
                double last_second_edge_intersection_point_z = snell_law_path[snell_law_path.size() - 3].z();

                // if the ray at critical angle pass the exit edge v0
                if (vertex_in_edge_not_endpoint(last_ray_point_x, last_edge_intersection_point_x,
                                                last_ray_point_y, last_edge_intersection_point_y,
                                                last_ray_point_z, last_edge_intersection_point_z,
                                                edge_sequence[exit_edge_index].adjacent_vertices()[0]->x(),
                                                edge_sequence[exit_edge_index].adjacent_vertices()[0]->y(),
                                                edge_sequence[exit_edge_index].adjacent_vertices()[0]->z()))
                {
                    // if the exit edge v0 is on the left of the ray
                    if (point_on_left_of_ray(
                            last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                            last_second_edge_intersection_point_x, last_second_edge_intersection_point_y, last_second_edge_intersection_point_z,
                            edge_sequence[exit_edge_index].adjacent_vertices()[0]->x(),
                            edge_sequence[exit_edge_index].adjacent_vertices()[0]->y(),
                            edge_sequence[exit_edge_index].adjacent_vertices()[0]->z()))
                    {
                        if (v0_of_curr_edge_is_on_left_of_ray)
                        {
                            curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                        else
                        {
                            curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                    }
                    // if the exit edge v0 is on the right of the ray
                    else if (point_on_left_right_vertical_of_ray(
                                 last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                                 last_second_edge_intersection_point_x, last_second_edge_intersection_point_y, last_second_edge_intersection_point_z,
                                 edge_sequence[exit_edge_index].adjacent_vertices()[0]->x(),
                                 edge_sequence[exit_edge_index].adjacent_vertices()[0]->y(),
                                 edge_sequence[exit_edge_index].adjacent_vertices()[0]->z()) < 0)
                    {
                        if (v0_of_curr_edge_is_on_left_of_ray)
                        {
                            curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                        else
                        {
                            curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                    }
                }

                // if the ray at critical angle pass the exit edge v1
                else if (vertex_in_edge_not_endpoint(last_ray_point_x, last_edge_intersection_point_x,
                                                     last_ray_point_y, last_edge_intersection_point_y,
                                                     last_ray_point_z, last_edge_intersection_point_z,
                                                     edge_sequence[exit_edge_index].adjacent_vertices()[1]->x(),
                                                     edge_sequence[exit_edge_index].adjacent_vertices()[1]->y(),
                                                     edge_sequence[exit_edge_index].adjacent_vertices()[1]->z()))
                {
                    // if the exit edge v1 is on the left of the ray
                    if (point_on_left_of_ray(
                            last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                            last_second_edge_intersection_point_x, last_second_edge_intersection_point_y, last_second_edge_intersection_point_z,
                            edge_sequence[exit_edge_index].adjacent_vertices()[1]->x(),
                            edge_sequence[exit_edge_index].adjacent_vertices()[1]->y(),
                            edge_sequence[exit_edge_index].adjacent_vertices()[1]->z()))
                    {
                        if (v0_of_curr_edge_is_on_left_of_ray)
                        {
                            curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                        else
                        {
                            curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                    }
                    // if the exit edge v1 is on the right of the ray
                    else if (point_on_left_right_vertical_of_ray(
                                 last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                                 last_second_edge_intersection_point_x, last_second_edge_intersection_point_y, last_second_edge_intersection_point_z,
                                 edge_sequence[exit_edge_index].adjacent_vertices()[1]->x(),
                                 edge_sequence[exit_edge_index].adjacent_vertices()[1]->y(),
                                 edge_sequence[exit_edge_index].adjacent_vertices()[1]->z()) < 0)
                    {
                        if (v0_of_curr_edge_is_on_left_of_ray)
                        {
                            curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                        else
                        {
                            curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                    }
                }
            }
        }
    }

    curr_edge_intersection_point_factor = (curr_edge_one_side_vertex_factor + curr_edge_another_side_vertex_factor) / 2;

    next_source_x = curr_edge_intersection_point_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
    next_source_y = curr_edge_intersection_point_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
    next_source_z = curr_edge_intersection_point_factor * curr_edge_vector_z + curr_edge_start_vertex_z;

    binary_search_of_snell_law_path_count_return = binary_search_of_snell_law_path_count;
}

void baseline_binary_search_multiple_times_of_each_edge(geodesic::Mesh *mesh,
                                                        std::vector<geodesic::Edge> &edge_sequence,
                                                        std::vector<geodesic::Face> &face_sequence,
                                                        geodesic::SurfacePoint &source,
                                                        geodesic::SurfacePoint &destination,
                                                        std::vector<geodesic::SurfacePoint> &result_path,
                                                        double delta,
                                                        int &total_binary_search_of_snell_law_path_count,
                                                        double &total_snell_law_memory_size,
                                                        int segment_path_length)
{
    double source_x = source.x();
    double source_y = source.y();
    double source_z = source.z();

    double next_source_x;
    double next_source_y;
    double next_source_z;

    std::vector<geodesic::SurfacePoint> snell_law_path;
    geodesic::Vertex *snell_law_path_point = new geodesic::Vertex();

    while (edge_sequence.size() > 0)
    {
        int binary_search_of_snell_law_path_count = 0;
        double snell_law_memory_size = 0;
        snell_law_path.clear();
        snell_law_path_point->set(source_x, source_y, source_z);
        result_path.push_back(geodesic::SurfacePoint(snell_law_path_point));
        baseline_binary_search_one_time_of_one_edge(mesh, edge_sequence, face_sequence, source_x, source_y, source_z, destination, snell_law_path, delta, next_source_x, next_source_y, next_source_z, binary_search_of_snell_law_path_count, snell_law_memory_size, segment_path_length);
        total_binary_search_of_snell_law_path_count += binary_search_of_snell_law_path_count;
        total_snell_law_memory_size += snell_law_memory_size;
        edge_sequence.pop_back();
        face_sequence.pop_back();
        source_x = next_source_x;
        source_y = next_source_y;
        source_z = next_source_z;
    }

    for (int i = 1; i < snell_law_path.size() - 1; i++)
    {
        result_path.push_back(snell_law_path[i]);
    }
    result_path.push_back(destination);
}

// in this effective weight algorithm, when the ray reach the destination point, we just use this ray to calculate the effective
// weight, and will not consider two rays
void effective_weight_binary_search_one_time_one_ray_of_one_edge(geodesic::Mesh *mesh,
                                                                 std::vector<geodesic::Edge> &edge_sequence,
                                                                 std::vector<geodesic::Face> &face_sequence,
                                                                 //  geodesic::SurfacePoint &source,
                                                                 double source_x, double source_y, double source_z,
                                                                 geodesic::SurfacePoint &destination,
                                                                 std::vector<geodesic::SurfacePoint> &snell_law_path,
                                                                 double delta,
                                                                 double &next_source_x,
                                                                 double &next_source_y,
                                                                 double &next_source_z,
                                                                 int &binary_search_and_effective_weight_of_snell_law_path_count_return,
                                                                 double &snell_law_memory_size,
                                                                 int segment_path_length)
{

    geodesic::Edge *first_edge;
    first_edge = &edge_sequence[edge_sequence.size() - 1]; // calculate the first edge that is just opposite of the source

    geodesic::Edge *last_edge;
    last_edge = &edge_sequence[0];

    geodesic::Face *first_face;
    first_face = &face_sequence[face_sequence.size() - 1];

    geodesic::Face *last_face;
    last_face = &face_sequence[0];

    geodesic::vertex_pointer last_edge_v0 = last_edge->adjacent_vertices()[0];
    geodesic::vertex_pointer last_edge_v1 = last_edge->adjacent_vertices()[1];

    double curr_edge_one_side_vertex_x;
    double curr_edge_one_side_vertex_y;
    double curr_edge_one_side_vertex_z;
    double curr_edge_another_side_vertex_x;
    double curr_edge_another_side_vertex_y;
    double curr_edge_another_side_vertex_z;

    double curr_edge_one_side_vertex_factor = -1;
    double curr_edge_another_side_vertex_factor = 0;

    bool snell_law_path_inside_edge_region = false;

    curr_edge_one_side_vertex_x = first_edge->adjacent_vertices()[0]->x();
    curr_edge_one_side_vertex_y = first_edge->adjacent_vertices()[0]->y();
    curr_edge_one_side_vertex_z = first_edge->adjacent_vertices()[0]->z();
    curr_edge_another_side_vertex_x = first_edge->adjacent_vertices()[1]->x();
    curr_edge_another_side_vertex_y = first_edge->adjacent_vertices()[1]->y();
    curr_edge_another_side_vertex_z = first_edge->adjacent_vertices()[1]->z();

    double curr_edge_vector_x = curr_edge_another_side_vertex_x - curr_edge_one_side_vertex_x;
    double curr_edge_vector_y = curr_edge_another_side_vertex_y - curr_edge_one_side_vertex_y;
    double curr_edge_vector_z = curr_edge_another_side_vertex_z - curr_edge_one_side_vertex_z;

    double curr_edge_start_vertex_x = first_edge->adjacent_vertices()[1]->x();
    double curr_edge_start_vertex_y = first_edge->adjacent_vertices()[1]->y();
    double curr_edge_start_vertex_z = first_edge->adjacent_vertices()[1]->z();

    double curr_edge_intersection_point_x;
    double curr_edge_intersection_point_y;
    double curr_edge_intersection_point_z;

    double curr_edge_intersection_point_factor = 0;

    int exit_edge_index;
    int out_at_critical_angle;

    double intersection_point_between_out_ray_and_edge_next_to_destination_x;
    double intersection_point_between_out_ray_and_edge_next_to_destination_y;
    double intersection_point_between_out_ray_and_edge_next_to_destination_z;

    double effective_weight_ratio = 0;

    // this count the number of times that we use the effective weight algorithm
    int use_effective_weight_times = 0;

    // we don't want to run the effective weight algorithm for too many times
    int max_use_effective_weight_times = 1;

    int binary_search_of_snell_law_path_count = 0;
    int effective_weight_of_snell_law_path_count = 0;

    // calcualte longest edge
    double max_edge_length = -1e100;
    std::vector<double> edge_length_list;
    std::vector<double> edge_length_without_outliers_list;
    for (unsigned i = 0; i < mesh->edges().size(); ++i)
    {
        geodesic::Edge &e = mesh->edges()[i];
        double edge_length = e.length();
        edge_length_list.push_back(edge_length);
    }
    remove_outliers(edge_length_list, edge_length_without_outliers_list);
    for (int i = 0; i < edge_length_without_outliers_list.size(); i++)
    {
        max_edge_length = std::max(max_edge_length, edge_length_without_outliers_list[i]);
    }

    int max_iteration = floor(log(max_edge_length / delta) / log(2));

    while (distance_xyz(curr_edge_one_side_vertex_x, curr_edge_one_side_vertex_y,
                        curr_edge_one_side_vertex_z, curr_edge_another_side_vertex_x,
                        curr_edge_another_side_vertex_y, curr_edge_another_side_vertex_z) >= delta)
    {
        if (binary_search_of_snell_law_path_count + effective_weight_of_snell_law_path_count > max_iteration)
        {
            break;
        }

        // if the ray doesn't pass the whole edge sequence or the effective weight ratio is unknown,
        // we use the midpoint to calculate the intersection point
        if (!snell_law_path_inside_edge_region || effective_weight_ratio == 0)
        {
            binary_search_of_snell_law_path_count++;

            snell_law_path.clear();

            curr_edge_intersection_point_factor = (curr_edge_one_side_vertex_factor + curr_edge_another_side_vertex_factor) / 2;

            curr_edge_intersection_point_x = curr_edge_intersection_point_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
            curr_edge_intersection_point_y = curr_edge_intersection_point_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
            curr_edge_intersection_point_z = curr_edge_intersection_point_factor * curr_edge_vector_z + curr_edge_start_vertex_z;

            auto start = std::chrono::high_resolution_clock::now();

            snell_law_path_inside_edge_region = snell_law_given_edge_intersection_point(
                mesh, edge_sequence, face_sequence, source_x, source_y, source_z, destination, curr_edge_intersection_point_x, curr_edge_intersection_point_y, curr_edge_intersection_point_z, snell_law_path, exit_edge_index, out_at_critical_angle, segment_path_length);

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

            snell_law_memory_size += snell_law_path.size() * sizeof(geodesic::SurfacePoint);

            // check the left or right side of the current edge
            // if the v0 of the current edge is on the left of the ray
            bool v0_of_curr_edge_is_on_left_of_ray = false;
            if (point_on_left_of_ray(
                    curr_edge_intersection_point_x, curr_edge_intersection_point_y, curr_edge_intersection_point_z,
                    source_x, source_y, source_z,
                    curr_edge_one_side_vertex_x, curr_edge_one_side_vertex_y, curr_edge_one_side_vertex_z))
            {
                v0_of_curr_edge_is_on_left_of_ray = true;
            }

            // if the ray pass the whole edge sequence
            if (snell_law_path_inside_edge_region)
            {
                double last_ray_point_x = snell_law_path[snell_law_path.size() - 1].x();
                double last_ray_point_y = snell_law_path[snell_law_path.size() - 1].y();
                double last_ray_point_z = snell_law_path[snell_law_path.size() - 1].z();

                double last_edge_intersection_point_x = snell_law_path[snell_law_path.size() - 2].x();
                double last_edge_intersection_point_y = snell_law_path[snell_law_path.size() - 2].y();
                double last_edge_intersection_point_z = snell_law_path[snell_law_path.size() - 2].z();

                double out_ray_vector_x = last_ray_point_x - last_edge_intersection_point_x;
                double out_ray_vector_y = last_ray_point_y - last_edge_intersection_point_y;
                double out_ray_vector_z = last_ray_point_z - last_edge_intersection_point_z;

                geodesic::edge_pointer last_ray_intersection_right_edge_next_to_destination;
                geodesic::edge_pointer last_ray_intersection_left_edge_next_to_destination;

                // if the ray pass the destination point
                if (vertex_on_edge(last_ray_point_x, last_ray_point_y, last_ray_point_z,
                                   last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                                   destination.x(), destination.y(), destination.z()))
                {
                    return;
                }
                // if the destination point is on the left of the ray
                else if (point_on_left_of_ray(
                             last_ray_point_x, last_ray_point_y, last_ray_point_z,
                             last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                             destination.x(), destination.y(), destination.z()))
                {
                    double last_e_v0_ray_position = point_on_left_right_vertical_of_ray(
                        last_ray_point_x, last_ray_point_y, last_ray_point_z,
                        last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                        last_edge_v0->x(), last_edge_v0->y(), last_edge_v0->z());

                    double last_e_v1_ray_position = point_on_left_right_vertical_of_ray(
                        last_ray_point_x, last_ray_point_y, last_ray_point_z,
                        last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                        last_edge_v1->x(), last_edge_v1->y(), last_edge_v1->z());

                    // if the v0 vertex of the last edge is on the left of the ray
                    // the opposite edge of v0 on the last face is the intersection edge (next to destination) of the ray
                    if (last_e_v0_ray_position >= 0)
                    {
                        last_ray_intersection_right_edge_next_to_destination = last_face->opposite_edge(last_edge_v0);
                        assert(last_ray_intersection_right_edge_next_to_destination->belongs(static_cast<geodesic::vertex_pointer>(destination.base_element())));
                    }
                    // if the v1 vertex of the last edge is on the left of the ray
                    // the opposite edge of v1 on the last face is the intersection edge (next to destination) of the ray
                    else if (last_e_v1_ray_position >= 0)
                    {
                        last_ray_intersection_right_edge_next_to_destination = last_face->opposite_edge(last_edge_v1);
                        assert(last_ray_intersection_right_edge_next_to_destination->belongs(static_cast<geodesic::vertex_pointer>(destination.base_element())));
                    }
                    // we then calculate the ray intersection point on the edge (next to destination)
                    calculate_out_ray_intersection_point_on_edge_next_to_destination(
                        last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                        last_ray_intersection_right_edge_next_to_destination->v0()->x(),
                        last_ray_intersection_right_edge_next_to_destination->v0()->y(),
                        last_ray_intersection_right_edge_next_to_destination->v0()->z(),
                        last_ray_intersection_right_edge_next_to_destination->v1()->x(),
                        last_ray_intersection_right_edge_next_to_destination->v1()->y(),
                        last_ray_intersection_right_edge_next_to_destination->v1()->z(),
                        out_ray_vector_x, out_ray_vector_y, out_ray_vector_z,
                        intersection_point_between_out_ray_and_edge_next_to_destination_x,
                        intersection_point_between_out_ray_and_edge_next_to_destination_y,
                        intersection_point_between_out_ray_and_edge_next_to_destination_z);

                    calculate_effective_weight_with_given_points(source_x, source_y, source_z,
                                                                 curr_edge_intersection_point_x,
                                                                 curr_edge_intersection_point_y,
                                                                 curr_edge_intersection_point_z,
                                                                 intersection_point_between_out_ray_and_edge_next_to_destination_x,
                                                                 intersection_point_between_out_ray_and_edge_next_to_destination_y,
                                                                 intersection_point_between_out_ray_and_edge_next_to_destination_z,
                                                                 first_face->adjacent_vertices()[0]->x(),
                                                                 first_face->adjacent_vertices()[0]->y(),
                                                                 first_face->adjacent_vertices()[0]->z(),
                                                                 first_face->adjacent_vertices()[1]->x(),
                                                                 first_face->adjacent_vertices()[1]->y(),
                                                                 first_face->adjacent_vertices()[1]->z(),
                                                                 first_face->adjacent_vertices()[2]->x(),
                                                                 first_face->adjacent_vertices()[2]->y(),
                                                                 first_face->adjacent_vertices()[2]->z(),
                                                                 first_edge->adjacent_vertices()[0]->x(),
                                                                 first_edge->adjacent_vertices()[0]->y(),
                                                                 first_edge->adjacent_vertices()[0]->z(),
                                                                 first_edge->adjacent_vertices()[1]->x(),
                                                                 first_edge->adjacent_vertices()[1]->y(),
                                                                 first_edge->adjacent_vertices()[1]->z(),
                                                                 effective_weight_ratio);

                    if (v0_of_curr_edge_is_on_left_of_ray)
                    {
                        curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                        curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                        curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                        curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                    }
                    else
                    {
                        curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                        curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                        curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                        curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                    }
                }
                // if the destination point is on the right of the ray
                else if (!point_on_left_of_ray(
                             last_ray_point_x, last_ray_point_y, last_ray_point_z,
                             last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                             destination.x(), destination.y(), destination.z()))
                {
                    double last_e_v0_ray_position = point_on_left_right_vertical_of_ray(
                        last_ray_point_x, last_ray_point_y, last_ray_point_z,
                        last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                        last_edge_v0->x(), last_edge_v0->y(), last_edge_v0->z());

                    double last_e_v1_ray_position = point_on_left_right_vertical_of_ray(
                        last_ray_point_x, last_ray_point_y, last_ray_point_z,
                        last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                        last_edge_v1->x(), last_edge_v1->y(), last_edge_v1->z());

                    // if the v0 vertex of the last edge is on the left of the ray
                    // the opposite edge of v1 on the last face is the intersection edge (next to destination) of the ray
                    if (last_e_v0_ray_position >= 0)
                    {
                        last_ray_intersection_left_edge_next_to_destination = last_face->opposite_edge(last_edge_v1);
                        assert(last_ray_intersection_left_edge_next_to_destination->belongs(static_cast<geodesic::vertex_pointer>(destination.base_element())));
                    }
                    // if the v1 vertex of the last edge is on the left of the ray
                    // the opposite edge of v0 on the last face is the intersection edge (next to destination) of the ray
                    else if (last_e_v1_ray_position >= 0)
                    {
                        last_ray_intersection_left_edge_next_to_destination = last_face->opposite_edge(last_edge_v0);
                        assert(last_ray_intersection_left_edge_next_to_destination->belongs(static_cast<geodesic::vertex_pointer>(destination.base_element())));
                    }

                    // we then calculate the ray intersection point on the edge (next to destination)
                    calculate_out_ray_intersection_point_on_edge_next_to_destination(
                        last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                        last_ray_intersection_left_edge_next_to_destination->v0()->x(),
                        last_ray_intersection_left_edge_next_to_destination->v0()->y(),
                        last_ray_intersection_left_edge_next_to_destination->v0()->z(),
                        last_ray_intersection_left_edge_next_to_destination->v1()->x(),
                        last_ray_intersection_left_edge_next_to_destination->v1()->y(),
                        last_ray_intersection_left_edge_next_to_destination->v1()->z(),
                        out_ray_vector_x, out_ray_vector_y, out_ray_vector_z,
                        intersection_point_between_out_ray_and_edge_next_to_destination_x,
                        intersection_point_between_out_ray_and_edge_next_to_destination_y,
                        intersection_point_between_out_ray_and_edge_next_to_destination_z);

                    calculate_effective_weight_with_given_points(source_x, source_y, source_z,
                                                                 curr_edge_intersection_point_x,
                                                                 curr_edge_intersection_point_y,
                                                                 curr_edge_intersection_point_z,
                                                                 intersection_point_between_out_ray_and_edge_next_to_destination_x,
                                                                 intersection_point_between_out_ray_and_edge_next_to_destination_y,
                                                                 intersection_point_between_out_ray_and_edge_next_to_destination_z,
                                                                 first_face->adjacent_vertices()[0]->x(),
                                                                 first_face->adjacent_vertices()[0]->y(),
                                                                 first_face->adjacent_vertices()[0]->z(),
                                                                 first_face->adjacent_vertices()[1]->x(),
                                                                 first_face->adjacent_vertices()[1]->y(),
                                                                 first_face->adjacent_vertices()[1]->z(),
                                                                 first_face->adjacent_vertices()[2]->x(),
                                                                 first_face->adjacent_vertices()[2]->y(),
                                                                 first_face->adjacent_vertices()[2]->z(),
                                                                 first_edge->adjacent_vertices()[0]->x(),
                                                                 first_edge->adjacent_vertices()[0]->y(),
                                                                 first_edge->adjacent_vertices()[0]->z(),
                                                                 first_edge->adjacent_vertices()[1]->x(),
                                                                 first_edge->adjacent_vertices()[1]->y(),
                                                                 first_edge->adjacent_vertices()[1]->z(),
                                                                 effective_weight_ratio);

                    if (!v0_of_curr_edge_is_on_left_of_ray)
                    {
                        curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                        curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                        curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                        curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                    }
                    else
                    {
                        curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                        curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                        curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                        curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                    }
                }
            }
            // if the ray doesn't pass the whole edge sequence and leave at edge_sequence[exit_edge_index - 1]
            // note that this ray does pass edge_sequence[exit_edge_index]
            else
            {
                // not out at critical angle (common case)
                if (out_at_critical_angle == 1)
                {
                    double last_ray_point_x = snell_law_path[snell_law_path.size() - 1].x();
                    double last_ray_point_y = snell_law_path[snell_law_path.size() - 1].y();
                    double last_ray_point_z = snell_law_path[snell_law_path.size() - 1].z();
                    double last_edge_intersection_point_x = snell_law_path[snell_law_path.size() - 2].x();
                    double last_edge_intersection_point_y = snell_law_path[snell_law_path.size() - 2].y();
                    double last_edge_intersection_point_z = snell_law_path[snell_law_path.size() - 2].z();

                    double e_v0_ray_position = point_on_left_right_vertical_of_ray(
                        last_ray_point_x, last_ray_point_y, last_ray_point_z,
                        last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                        edge_sequence[exit_edge_index - 1].adjacent_vertices()[0]->x(),
                        edge_sequence[exit_edge_index - 1].adjacent_vertices()[0]->y(),
                        edge_sequence[exit_edge_index - 1].adjacent_vertices()[0]->z());

                    double e_v1_ray_position = point_on_left_right_vertical_of_ray(
                        last_ray_point_x, last_ray_point_y, last_ray_point_z,
                        last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                        edge_sequence[exit_edge_index - 1].adjacent_vertices()[1]->x(),
                        edge_sequence[exit_edge_index - 1].adjacent_vertices()[1]->y(),
                        edge_sequence[exit_edge_index - 1].adjacent_vertices()[1]->z());

                    // if the two points on the edge_sequence[exit_edge_index - 1] is on the left of the ray
                    if ((e_v0_ray_position > 0 && e_v1_ray_position > 0) ||
                        (e_v0_ray_position == 0 && e_v1_ray_position > 0) ||
                        (e_v0_ray_position > 0 && e_v1_ray_position == 0) ||
                        (e_v0_ray_position < 0 && e_v0_ray_position > -1e-3 && e_v1_ray_position > 0) ||
                        (e_v0_ray_position > 0 && e_v1_ray_position < 0 && e_v1_ray_position > -1e-3))
                    {
                        if (v0_of_curr_edge_is_on_left_of_ray)
                        {
                            curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                        else
                        {
                            curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                    }
                    // if the two points on the edge_sequence[exit_edge_index - 1] is on the right of the ray
                    else if ((e_v0_ray_position < 0 && e_v1_ray_position < 0) ||
                             (e_v0_ray_position == 0 && e_v1_ray_position < 0) ||
                             (e_v0_ray_position < 0 && e_v1_ray_position == 0) ||
                             (e_v0_ray_position > 0 && e_v0_ray_position < 1e-3 && e_v1_ray_position < 0) ||
                             (e_v0_ray_position < 0 && e_v1_ray_position > 0 && e_v1_ray_position < 1e-3))
                    {
                        if (!v0_of_curr_edge_is_on_left_of_ray)
                        {
                            curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                        else
                        {
                            curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                    }
                }
                // out at critical angle
                else if (out_at_critical_angle == 2)
                {
                    double last_ray_point_x = snell_law_path[snell_law_path.size() - 1].x();
                    double last_ray_point_y = snell_law_path[snell_law_path.size() - 1].y();
                    double last_ray_point_z = snell_law_path[snell_law_path.size() - 1].z();
                    double last_edge_intersection_point_x = snell_law_path[snell_law_path.size() - 2].x();
                    double last_edge_intersection_point_y = snell_law_path[snell_law_path.size() - 2].y();
                    double last_edge_intersection_point_z = snell_law_path[snell_law_path.size() - 2].z();
                    double last_second_edge_intersection_point_x = snell_law_path[snell_law_path.size() - 3].x();
                    double last_second_edge_intersection_point_y = snell_law_path[snell_law_path.size() - 3].y();
                    double last_second_edge_intersection_point_z = snell_law_path[snell_law_path.size() - 3].z();
                    // if the ray at critical angle pass the exit edge v0
                    if (vertex_in_edge_not_endpoint(last_ray_point_x, last_edge_intersection_point_x,
                                                    last_ray_point_y, last_edge_intersection_point_y,
                                                    last_ray_point_z, last_edge_intersection_point_z,
                                                    edge_sequence[exit_edge_index].adjacent_vertices()[0]->x(),
                                                    edge_sequence[exit_edge_index].adjacent_vertices()[0]->y(),
                                                    edge_sequence[exit_edge_index].adjacent_vertices()[0]->z()))
                    {
                        // if the exit edge v0 is on the left of the ray
                        if (point_on_left_of_ray(
                                last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                                last_second_edge_intersection_point_x, last_second_edge_intersection_point_y, last_second_edge_intersection_point_z,
                                edge_sequence[exit_edge_index].adjacent_vertices()[0]->x(),
                                edge_sequence[exit_edge_index].adjacent_vertices()[0]->y(),
                                edge_sequence[exit_edge_index].adjacent_vertices()[0]->z()))
                        {
                            if (v0_of_curr_edge_is_on_left_of_ray)
                            {
                                curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                            else
                            {
                                curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                        }
                        // if the exit edge v0 is on the right of the ray
                        else if (point_on_left_right_vertical_of_ray(
                                     last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                                     last_second_edge_intersection_point_x, last_second_edge_intersection_point_y, last_second_edge_intersection_point_z,
                                     edge_sequence[exit_edge_index].adjacent_vertices()[0]->x(),
                                     edge_sequence[exit_edge_index].adjacent_vertices()[0]->y(),
                                     edge_sequence[exit_edge_index].adjacent_vertices()[0]->z()) < 0)
                        {
                            if (v0_of_curr_edge_is_on_left_of_ray)
                            {
                                curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                            else
                            {
                                curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                        }
                    }

                    // if the ray at critical angle pass the exit edge v1
                    else if (vertex_in_edge_not_endpoint(last_ray_point_x, last_edge_intersection_point_x,
                                                         last_ray_point_y, last_edge_intersection_point_y,
                                                         last_ray_point_z, last_edge_intersection_point_z,
                                                         edge_sequence[exit_edge_index].adjacent_vertices()[1]->x(),
                                                         edge_sequence[exit_edge_index].adjacent_vertices()[1]->y(),
                                                         edge_sequence[exit_edge_index].adjacent_vertices()[1]->z()))
                    {
                        // if the exit edge v1 is on the left of the ray
                        if (point_on_left_of_ray(
                                last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                                last_second_edge_intersection_point_x, last_second_edge_intersection_point_y, last_second_edge_intersection_point_z,
                                edge_sequence[exit_edge_index].adjacent_vertices()[1]->x(),
                                edge_sequence[exit_edge_index].adjacent_vertices()[1]->y(),
                                edge_sequence[exit_edge_index].adjacent_vertices()[1]->z()))
                        {
                            if (v0_of_curr_edge_is_on_left_of_ray)
                            {
                                curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                            else
                            {
                                curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                        }
                        // if the exit edge v1 is on the right of the ray
                        else if (point_on_left_right_vertical_of_ray(
                                     last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                                     last_second_edge_intersection_point_x, last_second_edge_intersection_point_y, last_second_edge_intersection_point_z,
                                     edge_sequence[exit_edge_index].adjacent_vertices()[1]->x(),
                                     edge_sequence[exit_edge_index].adjacent_vertices()[1]->y(),
                                     edge_sequence[exit_edge_index].adjacent_vertices()[1]->z()) < 0)
                        {
                            if (v0_of_curr_edge_is_on_left_of_ray)
                            {
                                curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                            else
                            {
                                curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                        }
                    }
                }
            }
        }

        // if we have run the effecive weight algorithm for too many times,
        // we use the midpoint to calculate the intersection point, and don't calculate the effective weight for the first face and
        // remaining face anymore
        else if ((use_effective_weight_times >= max_use_effective_weight_times) &&
                 (snell_law_path_inside_edge_region && effective_weight_ratio != 0))
        {
            binary_search_of_snell_law_path_count++;

            snell_law_path.clear();

            curr_edge_intersection_point_factor = (curr_edge_one_side_vertex_factor + curr_edge_another_side_vertex_factor) / 2;

            curr_edge_intersection_point_x = curr_edge_intersection_point_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
            curr_edge_intersection_point_y = curr_edge_intersection_point_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
            curr_edge_intersection_point_z = curr_edge_intersection_point_factor * curr_edge_vector_z + curr_edge_start_vertex_z;

            auto start = std::chrono::high_resolution_clock::now();

            snell_law_path_inside_edge_region = snell_law_given_edge_intersection_point(
                mesh, edge_sequence, face_sequence, source_x, source_y, source_z, destination, curr_edge_intersection_point_x, curr_edge_intersection_point_y, curr_edge_intersection_point_z, snell_law_path, exit_edge_index, out_at_critical_angle, segment_path_length);

            snell_law_memory_size += snell_law_path.size() * sizeof(geodesic::SurfacePoint);

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

            // check the left or right side of the current edge
            // if the v0 of the current edge is on the left of the ray
            bool v0_of_curr_edge_is_on_left_of_ray = false;
            if (point_on_left_of_ray(
                    curr_edge_intersection_point_x, curr_edge_intersection_point_y, curr_edge_intersection_point_z,
                    source_x, source_y, source_z,
                    curr_edge_one_side_vertex_x, curr_edge_one_side_vertex_y, curr_edge_one_side_vertex_z))
            {
                v0_of_curr_edge_is_on_left_of_ray = true;
            }

            // if the ray pass the whole edge sequence
            if (snell_law_path_inside_edge_region)
            {
                double last_ray_point_x = snell_law_path[snell_law_path.size() - 1].x();
                double last_ray_point_y = snell_law_path[snell_law_path.size() - 1].y();
                double last_ray_point_z = snell_law_path[snell_law_path.size() - 1].z();

                double last_edge_intersection_point_x = snell_law_path[snell_law_path.size() - 2].x();
                double last_edge_intersection_point_y = snell_law_path[snell_law_path.size() - 2].y();
                double last_edge_intersection_point_z = snell_law_path[snell_law_path.size() - 2].z();

                // if the ray pass the destination point
                if (vertex_on_edge(last_ray_point_x, last_ray_point_y, last_ray_point_z,
                                   last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                                   destination.x(), destination.y(), destination.z()))
                {
                    return;
                }
                // if the destination point is on the left of the ray
                else if (point_on_left_of_ray(
                             last_ray_point_x, last_ray_point_y, last_ray_point_z,
                             last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                             destination.x(), destination.y(), destination.z()))
                {
                    if (v0_of_curr_edge_is_on_left_of_ray)
                    {
                        curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                        curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                        curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                        curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                    }
                    else
                    {
                        curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                        curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                        curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                        curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                    }
                }
                // if the destination point is on the right of the ray
                else if (!point_on_left_of_ray(
                             last_ray_point_x, last_ray_point_y, last_ray_point_z,
                             last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                             destination.x(), destination.y(), destination.z()))
                {
                    if (!v0_of_curr_edge_is_on_left_of_ray)
                    {
                        curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                        curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                        curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                        curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                    }
                    else
                    {
                        curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                        curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                        curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                        curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                    }
                }
            }
            // if the ray doesn't pass the whole edge sequence and leave at edge_sequence[exit_edge_index - 1]
            // note that this ray does pass edge_sequence[exit_edge_index]
            else
            {
                // not out at critical angle (common case)
                if (out_at_critical_angle == 1)
                {
                    double last_ray_point_x = snell_law_path[snell_law_path.size() - 1].x();
                    double last_ray_point_y = snell_law_path[snell_law_path.size() - 1].y();
                    double last_ray_point_z = snell_law_path[snell_law_path.size() - 1].z();
                    double last_edge_intersection_point_x = snell_law_path[snell_law_path.size() - 2].x();
                    double last_edge_intersection_point_y = snell_law_path[snell_law_path.size() - 2].y();
                    double last_edge_intersection_point_z = snell_law_path[snell_law_path.size() - 2].z();

                    double e_v0_ray_position = point_on_left_right_vertical_of_ray(
                        last_ray_point_x, last_ray_point_y, last_ray_point_z,
                        last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                        edge_sequence[exit_edge_index - 1].adjacent_vertices()[0]->x(),
                        edge_sequence[exit_edge_index - 1].adjacent_vertices()[0]->y(),
                        edge_sequence[exit_edge_index - 1].adjacent_vertices()[0]->z());

                    double e_v1_ray_position = point_on_left_right_vertical_of_ray(
                        last_ray_point_x, last_ray_point_y, last_ray_point_z,
                        last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                        edge_sequence[exit_edge_index - 1].adjacent_vertices()[1]->x(),
                        edge_sequence[exit_edge_index - 1].adjacent_vertices()[1]->y(),
                        edge_sequence[exit_edge_index - 1].adjacent_vertices()[1]->z());

                    // if the two points on the edge_sequence[exit_edge_index - 1] is on the left of the ray
                    if ((e_v0_ray_position > 0 && e_v1_ray_position > 0) ||
                        (e_v0_ray_position == 0 && e_v1_ray_position > 0) ||
                        (e_v0_ray_position > 0 && e_v1_ray_position == 0) ||
                        (e_v0_ray_position < 0 && e_v0_ray_position > -1e-3 && e_v1_ray_position > 0) ||
                        (e_v0_ray_position > 0 && e_v1_ray_position < 0 && e_v1_ray_position > -1e-3))
                    {
                        if (v0_of_curr_edge_is_on_left_of_ray)
                        {
                            curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                        else
                        {
                            curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                    }
                    // if the two points on the edge_sequence[exit_edge_index - 1] is on the right of the ray
                    else if ((e_v0_ray_position < 0 && e_v1_ray_position < 0) ||
                             (e_v0_ray_position == 0 && e_v1_ray_position < 0) ||
                             (e_v0_ray_position < 0 && e_v1_ray_position == 0) ||
                             (e_v0_ray_position > 0 && e_v0_ray_position < 1e-3 && e_v1_ray_position < 0) ||
                             (e_v0_ray_position < 0 && e_v1_ray_position > 0 && e_v1_ray_position < 1e-3))
                    {
                        if (!v0_of_curr_edge_is_on_left_of_ray)
                        {
                            curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                        else
                        {
                            curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                            curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                            curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                            curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                        }
                    }
                }
                // out at critical angle
                else if (out_at_critical_angle == 2)
                {
                    double last_ray_point_x = snell_law_path[snell_law_path.size() - 1].x();
                    double last_ray_point_y = snell_law_path[snell_law_path.size() - 1].y();
                    double last_ray_point_z = snell_law_path[snell_law_path.size() - 1].z();
                    double last_edge_intersection_point_x = snell_law_path[snell_law_path.size() - 2].x();
                    double last_edge_intersection_point_y = snell_law_path[snell_law_path.size() - 2].y();
                    double last_edge_intersection_point_z = snell_law_path[snell_law_path.size() - 2].z();
                    double last_second_edge_intersection_point_x = snell_law_path[snell_law_path.size() - 3].x();
                    double last_second_edge_intersection_point_y = snell_law_path[snell_law_path.size() - 3].y();
                    double last_second_edge_intersection_point_z = snell_law_path[snell_law_path.size() - 3].z();

                    // if the ray at critical angle pass the exit edge v0
                    if (vertex_in_edge_not_endpoint(last_ray_point_x, last_edge_intersection_point_x,
                                                    last_ray_point_y, last_edge_intersection_point_y,
                                                    last_ray_point_z, last_edge_intersection_point_z,
                                                    edge_sequence[exit_edge_index].adjacent_vertices()[0]->x(),
                                                    edge_sequence[exit_edge_index].adjacent_vertices()[0]->y(),
                                                    edge_sequence[exit_edge_index].adjacent_vertices()[0]->z()))
                    {
                        // if the exit edge v0 is on the left of the ray
                        if (point_on_left_of_ray(
                                last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                                last_second_edge_intersection_point_x, last_second_edge_intersection_point_y, last_second_edge_intersection_point_z,
                                edge_sequence[exit_edge_index].adjacent_vertices()[0]->x(),
                                edge_sequence[exit_edge_index].adjacent_vertices()[0]->y(),
                                edge_sequence[exit_edge_index].adjacent_vertices()[0]->z()))
                        {
                            if (v0_of_curr_edge_is_on_left_of_ray)
                            {
                                curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                            else
                            {
                                curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                        }
                        // if the exit edge v0 is on the right of the ray
                        else if (point_on_left_right_vertical_of_ray(
                                     last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                                     last_second_edge_intersection_point_x, last_second_edge_intersection_point_y, last_second_edge_intersection_point_z,
                                     edge_sequence[exit_edge_index].adjacent_vertices()[0]->x(),
                                     edge_sequence[exit_edge_index].adjacent_vertices()[0]->y(),
                                     edge_sequence[exit_edge_index].adjacent_vertices()[0]->z()) < 0)
                        {
                            if (v0_of_curr_edge_is_on_left_of_ray)
                            {
                                curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                            else
                            {
                                curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                        }
                    }

                    // if the ray at critical angle pass the exit edge v1
                    else if (vertex_in_edge_not_endpoint(last_ray_point_x, last_edge_intersection_point_x,
                                                         last_ray_point_y, last_edge_intersection_point_y,
                                                         last_ray_point_z, last_edge_intersection_point_z,
                                                         edge_sequence[exit_edge_index].adjacent_vertices()[1]->x(),
                                                         edge_sequence[exit_edge_index].adjacent_vertices()[1]->y(),
                                                         edge_sequence[exit_edge_index].adjacent_vertices()[1]->z()))
                    {
                        // if the exit edge v1 is on the left of the ray
                        if (point_on_left_of_ray(
                                last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                                last_second_edge_intersection_point_x, last_second_edge_intersection_point_y, last_second_edge_intersection_point_z,
                                edge_sequence[exit_edge_index].adjacent_vertices()[1]->x(),
                                edge_sequence[exit_edge_index].adjacent_vertices()[1]->y(),
                                edge_sequence[exit_edge_index].adjacent_vertices()[1]->z()))
                        {
                            if (v0_of_curr_edge_is_on_left_of_ray)
                            {
                                curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                            else
                            {
                                curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                        }
                        // if the exit edge v1 is on the right of the ray
                        else if (point_on_left_right_vertical_of_ray(
                                     last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                                     last_second_edge_intersection_point_x, last_second_edge_intersection_point_y, last_second_edge_intersection_point_z,
                                     edge_sequence[exit_edge_index].adjacent_vertices()[1]->x(),
                                     edge_sequence[exit_edge_index].adjacent_vertices()[1]->y(),
                                     edge_sequence[exit_edge_index].adjacent_vertices()[1]->z()) < 0)
                        {
                            if (v0_of_curr_edge_is_on_left_of_ray)
                            {
                                curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                            else
                            {
                                curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                            }
                        }
                    }
                }
            }
        }
        // if the ray passes the whole edge sequence and the effective weight ratio is known,
        // and we haven't run the effective weight algorithm for too many times,
        // we use the effective weight to calculate the intersection point
        else if ((snell_law_path_inside_edge_region && effective_weight_ratio != 0) &&
                 (use_effective_weight_times < max_use_effective_weight_times))
        {
            use_effective_weight_times++;

            double last_ray_point_x = snell_law_path[snell_law_path.size() - 1].x();
            double last_ray_point_y = snell_law_path[snell_law_path.size() - 1].y();
            double last_ray_point_z = snell_law_path[snell_law_path.size() - 1].z();

            double last_edge_intersection_point_x = snell_law_path[snell_law_path.size() - 2].x();
            double last_edge_intersection_point_y = snell_law_path[snell_law_path.size() - 2].y();
            double last_edge_intersection_point_z = snell_law_path[snell_law_path.size() - 2].z();

            // if the destination point is on the left of the ray
            bool destination_is_on_left_of_prev_ray = false;
            if (point_on_left_of_ray(
                    last_ray_point_x, last_ray_point_y, last_ray_point_z,
                    last_edge_intersection_point_x, last_edge_intersection_point_y, last_edge_intersection_point_z,
                    destination.x(), destination.y(), destination.z()))
            {
                destination_is_on_left_of_prev_ray = true;
            }

            snell_law_path.clear();

            // calculate the intersection point on the first edge using the effective weight

            std::vector<double> curr_edge_intersection_point_factor_list;
            std::vector<double> curr_edge_intersection_point_x_list;
            std::vector<double> curr_edge_intersection_point_y_list;
            std::vector<double> curr_edge_intersection_point_z_list;
            std::vector<bool> snell_law_path_inside_edge_region_list;
            std::vector<double> distance_of_intersection_point_between_out_ray_and_edge_next_to_destination_and_destination_list;

            calculate_intersection_point_on_first_edge_using_effective_weight(
                source_x, source_y, source_z,
                destination.x(), destination.y(), destination.z(),
                effective_weight_ratio,
                first_face->adjacent_vertices()[0]->x(),
                first_face->adjacent_vertices()[0]->y(),
                first_face->adjacent_vertices()[0]->z(),
                first_face->adjacent_vertices()[1]->x(),
                first_face->adjacent_vertices()[1]->y(),
                first_face->adjacent_vertices()[1]->z(),
                first_face->adjacent_vertices()[2]->x(),
                first_face->adjacent_vertices()[2]->y(),
                first_face->adjacent_vertices()[2]->z(),
                first_edge->adjacent_vertices()[0]->x(),
                first_edge->adjacent_vertices()[0]->y(),
                first_edge->adjacent_vertices()[0]->z(),
                first_edge->adjacent_vertices()[1]->x(),
                first_edge->adjacent_vertices()[1]->y(),
                first_edge->adjacent_vertices()[1]->z(),
                curr_edge_one_side_vertex_factor,
                curr_edge_another_side_vertex_factor,
                0.02,
                curr_edge_intersection_point_factor_list,
                curr_edge_intersection_point_x_list,
                curr_edge_intersection_point_y_list,
                curr_edge_intersection_point_z_list);

            for (int i = 0; i < curr_edge_intersection_point_factor_list.size(); ++i)
            {
                if (binary_search_of_snell_law_path_count + effective_weight_of_snell_law_path_count > max_iteration)
                {
                    break;
                }
                effective_weight_of_snell_law_path_count++;

                snell_law_path_inside_edge_region_list.push_back(snell_law_given_edge_intersection_point(
                    mesh, edge_sequence, face_sequence, source_x, source_y, source_z, destination, curr_edge_intersection_point_x_list[i], curr_edge_intersection_point_y_list[i], curr_edge_intersection_point_z_list[i], snell_law_path, exit_edge_index, out_at_critical_angle, segment_path_length));

                snell_law_memory_size += snell_law_path.size() * sizeof(geodesic::SurfacePoint);

                double ew_last_ray_point_x = snell_law_path[snell_law_path.size() - 1].x();
                double ew_last_ray_point_y = snell_law_path[snell_law_path.size() - 1].y();
                double ew_last_ray_point_z = snell_law_path[snell_law_path.size() - 1].z();

                double ew_last_edge_intersection_point_x = snell_law_path[snell_law_path.size() - 2].x();
                double ew_last_edge_intersection_point_y = snell_law_path[snell_law_path.size() - 2].y();
                double ew_last_edge_intersection_point_z = snell_law_path[snell_law_path.size() - 2].z();

                snell_law_path_inside_edge_region = snell_law_path_inside_edge_region_list[i];
                curr_edge_intersection_point_factor = curr_edge_intersection_point_factor_list[i];
                curr_edge_intersection_point_x = curr_edge_intersection_point_x_list[i];
                curr_edge_intersection_point_y = curr_edge_intersection_point_y_list[i];
                curr_edge_intersection_point_z = curr_edge_intersection_point_z_list[i];

                // check the left or right side of the current edge
                // if the v0 of the current edge is on the left of the ray
                bool v0_of_curr_edge_is_on_left_of_ray = false;
                if (point_on_left_of_ray(
                        curr_edge_intersection_point_x, curr_edge_intersection_point_y, curr_edge_intersection_point_z,
                        source_x, source_y, source_z,
                        first_edge->adjacent_vertices()[0]->x(), first_edge->adjacent_vertices()[0]->y(), first_edge->adjacent_vertices()[0]->z()))
                {
                    v0_of_curr_edge_is_on_left_of_ray = true;
                }

                // if the snell law path calculated using the effective weight passes the whole edge sequence
                if (snell_law_path_inside_edge_region)
                {

                    bool destination_is_on_left_of_curr_ray = false;
                    if (point_on_left_of_ray(ew_last_ray_point_x, ew_last_ray_point_y, ew_last_ray_point_z,
                                             ew_last_edge_intersection_point_x, ew_last_edge_intersection_point_y, ew_last_edge_intersection_point_z,
                                             destination.x(), destination.y(), destination.z()))
                    {
                        destination_is_on_left_of_curr_ray = true;
                    }

                    if (destination_is_on_left_of_prev_ray)
                    {
                        if (v0_of_curr_edge_is_on_left_of_ray)
                        {
                            // the last step is moving curr_edge_another_side_vertex to curr_edge_intersection_point

                            // if the ray pass the destination point
                            if (vertex_on_edge(ew_last_ray_point_x, ew_last_ray_point_y, ew_last_ray_point_z,
                                               ew_last_edge_intersection_point_x, ew_last_edge_intersection_point_y, ew_last_edge_intersection_point_z,
                                               destination.x(), destination.y(), destination.z()))
                            {
                                return;
                            }

                            else if (destination_is_on_left_of_curr_ray)
                            {
                                curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                break;
                            }

                            else if (!destination_is_on_left_of_curr_ray)
                            {
                                curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                continue;
                            }
                        }
                        else
                        {
                            // the last step is moving curr_edge_one_side_vertex to curr_edge_intersection_point

                            // if the ray pass the destination point
                            if (vertex_on_edge(ew_last_ray_point_x, ew_last_ray_point_y, ew_last_ray_point_z,
                                               ew_last_edge_intersection_point_x, ew_last_edge_intersection_point_y, ew_last_edge_intersection_point_z,
                                               destination.x(), destination.y(), destination.z()))
                            {
                                return;
                            }

                            else if (destination_is_on_left_of_curr_ray)
                            {
                                curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                continue;
                            }

                            else if (!destination_is_on_left_of_curr_ray)
                            {
                                curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                break;
                            }
                        }
                    }
                    else
                    {
                        if (v0_of_curr_edge_is_on_left_of_ray)
                        {
                            // the last step is moving curr_edge_one_side_vertex to curr_edge_intersection_point

                            // if the ray pass the destination point
                            if (vertex_on_edge(ew_last_ray_point_x, ew_last_ray_point_y, ew_last_ray_point_z,
                                               ew_last_edge_intersection_point_x, ew_last_edge_intersection_point_y, ew_last_edge_intersection_point_z,
                                               destination.x(), destination.y(), destination.z()))
                            {
                                return;
                            }

                            else if (destination_is_on_left_of_curr_ray)
                            {
                                curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                break;
                            }

                            else if (!destination_is_on_left_of_curr_ray)
                            {
                                curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                continue;
                            }
                        }
                        else
                        {
                            // the last step is moving curr_edge_another_side_vertex to curr_edge_intersection_point

                            // if the ray pass the destination point
                            if (vertex_on_edge(ew_last_ray_point_x, ew_last_ray_point_y, ew_last_ray_point_z,
                                               ew_last_edge_intersection_point_x, ew_last_edge_intersection_point_y, ew_last_edge_intersection_point_z,
                                               destination.x(), destination.y(), destination.z()))
                            {
                                return;
                            }

                            else if (destination_is_on_left_of_curr_ray)
                            {
                                curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                continue;
                            }

                            else if (!destination_is_on_left_of_curr_ray)
                            {
                                curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                break;
                            }
                        }
                    }
                }

                // if the snell law path calculated using the effective weight doesn't pass the whole edge sequence
                else
                {
                    if (out_at_critical_angle == 1)
                    {
                        bool exit_edge_is_on_left_of_ray = false;
                        if ((point_on_left_of_ray(
                                ew_last_ray_point_x, ew_last_ray_point_y, ew_last_ray_point_z,
                                ew_last_edge_intersection_point_x, ew_last_edge_intersection_point_y, ew_last_edge_intersection_point_z,
                                edge_sequence[exit_edge_index - 1].adjacent_vertices()[0]->x(),
                                edge_sequence[exit_edge_index - 1].adjacent_vertices()[0]->y(),
                                edge_sequence[exit_edge_index - 1].adjacent_vertices()[0]->z())) &&
                            (point_on_left_of_ray(
                                ew_last_ray_point_x, ew_last_ray_point_y, ew_last_ray_point_z,
                                ew_last_edge_intersection_point_x, ew_last_edge_intersection_point_y, ew_last_edge_intersection_point_z,
                                edge_sequence[exit_edge_index - 1].adjacent_vertices()[1]->x(),
                                edge_sequence[exit_edge_index - 1].adjacent_vertices()[1]->y(),
                                edge_sequence[exit_edge_index - 1].adjacent_vertices()[1]->z())))
                        {
                            exit_edge_is_on_left_of_ray = true;
                        }

                        if (destination_is_on_left_of_prev_ray)
                        {
                            if (v0_of_curr_edge_is_on_left_of_ray)
                            {
                                // the last step is moving curr_edge_another_side_vertex to curr_edge_intersection_point
                                if (exit_edge_is_on_left_of_ray)
                                {
                                    curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                    curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                    curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                    curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                    break;
                                }

                                else if (!exit_edge_is_on_left_of_ray)
                                {
                                    curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                    curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                    curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                    curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                    continue;
                                }
                            }
                            else
                            {
                                // the last step is moving curr_edge_one_side_vertex to curr_edge_intersection_point
                                if (exit_edge_is_on_left_of_ray)
                                {
                                    curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                    curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                    curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                    curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                    continue;
                                }

                                else if (!exit_edge_is_on_left_of_ray)
                                {
                                    curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                    curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                    curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                    curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                    break;
                                }
                            }
                        }
                        else
                        {
                            if (v0_of_curr_edge_is_on_left_of_ray)
                            {
                                // the last step is moving curr_edge_one_side_vertex to curr_edge_intersection_point
                                if (exit_edge_is_on_left_of_ray)
                                {
                                    curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                    curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                    curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                    curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                    break;
                                }

                                if (!exit_edge_is_on_left_of_ray)
                                {
                                    curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                    curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                    curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                    curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                    continue;
                                }
                            }
                            else
                            {
                                // the last step is moving curr_edge_another_side_vertex to curr_edge_intersection_point
                                if (exit_edge_is_on_left_of_ray)
                                {
                                    curr_edge_one_side_vertex_factor = curr_edge_intersection_point_factor;

                                    curr_edge_one_side_vertex_x = curr_edge_one_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                    curr_edge_one_side_vertex_y = curr_edge_one_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                    curr_edge_one_side_vertex_z = curr_edge_one_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                    continue;
                                }

                                if (!exit_edge_is_on_left_of_ray)
                                {
                                    curr_edge_another_side_vertex_factor = curr_edge_intersection_point_factor;

                                    curr_edge_another_side_vertex_x = curr_edge_another_side_vertex_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
                                    curr_edge_another_side_vertex_y = curr_edge_another_side_vertex_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
                                    curr_edge_another_side_vertex_z = curr_edge_another_side_vertex_factor * curr_edge_vector_z + curr_edge_start_vertex_z;
                                    break;
                                }
                            }
                        }
                    }
                    else if (out_at_critical_angle == 2)
                    {
                        break;
                    }
                }
            }
        }
    }

    curr_edge_intersection_point_factor = (curr_edge_one_side_vertex_factor + curr_edge_another_side_vertex_factor) / 2;

    next_source_x = curr_edge_intersection_point_factor * curr_edge_vector_x + curr_edge_start_vertex_x;
    next_source_y = curr_edge_intersection_point_factor * curr_edge_vector_y + curr_edge_start_vertex_y;
    next_source_z = curr_edge_intersection_point_factor * curr_edge_vector_z + curr_edge_start_vertex_z;

    binary_search_and_effective_weight_of_snell_law_path_count_return = binary_search_of_snell_law_path_count + effective_weight_of_snell_law_path_count;
}

void effective_weight_binary_search_multiple_times_of_each_edge(geodesic::Mesh *mesh,
                                                                std::vector<geodesic::Edge> &edge_sequence,
                                                                std::vector<geodesic::Face> &face_sequence,
                                                                geodesic::SurfacePoint &source,
                                                                geodesic::SurfacePoint &destination,
                                                                std::vector<geodesic::SurfacePoint> &result_path,
                                                                double delta,
                                                                int &total_binary_search_and_effective_weight_of_snell_law_path_count,
                                                                double &total_snell_law_memory_size,
                                                                int segment_path_length)
{
    // the first num_of_edges_using_effective_weight edges use effective weight (and also binary search)
    // the remaining edges only use binary search
    int num_of_edges_using_effective_weight = edge_sequence.size();

    double source_x = source.x();
    double source_y = source.y();
    double source_z = source.z();

    std::vector<geodesic::SurfacePoint> snell_law_path;
    geodesic::Vertex *snell_law_path_point = new geodesic::Vertex();

    double next_source_x;
    double next_source_y;
    double next_source_z;

    if (num_of_edges_using_effective_weight <= edge_sequence.size())
    {
        for (int i = 0; i < num_of_edges_using_effective_weight; ++i)
        {
            int binary_search_and_effective_weight_of_snell_law_path_count = 0;
            double snell_law_memory_size = 0;
            snell_law_path.clear();
            snell_law_path_point->set(source_x, source_y, source_z);
            result_path.push_back(geodesic::SurfacePoint(snell_law_path_point));
            effective_weight_binary_search_one_time_one_ray_of_one_edge(mesh, edge_sequence, face_sequence, source_x, source_y, source_z, destination, snell_law_path, delta, next_source_x, next_source_y, next_source_z, binary_search_and_effective_weight_of_snell_law_path_count, snell_law_memory_size, segment_path_length);
            total_binary_search_and_effective_weight_of_snell_law_path_count += binary_search_and_effective_weight_of_snell_law_path_count;
            total_snell_law_memory_size += snell_law_memory_size;
            edge_sequence.pop_back();
            face_sequence.pop_back();
            source_x = next_source_x;
            source_y = next_source_y;
            source_z = next_source_z;
        }
        while (edge_sequence.size() > 0)
        {
            int binary_search_and_effective_weight_of_snell_law_path_count = 0;
            double snell_law_memory_size = 0;
            snell_law_path.clear();
            snell_law_path_point->set(source_x, source_y, source_z);
            result_path.push_back(geodesic::SurfacePoint(snell_law_path_point));
            baseline_binary_search_one_time_of_one_edge(mesh, edge_sequence, face_sequence, source_x, source_y, source_z, destination, snell_law_path, delta, next_source_x, next_source_y, next_source_z, binary_search_and_effective_weight_of_snell_law_path_count, snell_law_memory_size, segment_path_length);
            total_binary_search_and_effective_weight_of_snell_law_path_count += binary_search_and_effective_weight_of_snell_law_path_count;
            total_snell_law_memory_size += snell_law_memory_size;
            edge_sequence.pop_back();
            face_sequence.pop_back();
            source_x = next_source_x;
            source_y = next_source_y;
            source_z = next_source_z;
        }
    }
    else
    {
        while (edge_sequence.size() > 0)
        {
            int binary_search_and_effective_weight_of_snell_law_path_count = 0;
            double snell_law_memory_size = 0;
            snell_law_path.clear();
            snell_law_path_point->set(source_x, source_y, source_z);
            result_path.push_back(geodesic::SurfacePoint(snell_law_path_point));
            effective_weight_binary_search_one_time_one_ray_of_one_edge(mesh, edge_sequence, face_sequence, source_x, source_y, source_z, destination, snell_law_path, delta, next_source_x, next_source_y, next_source_z, binary_search_and_effective_weight_of_snell_law_path_count, snell_law_memory_size, segment_path_length);
            total_binary_search_and_effective_weight_of_snell_law_path_count += binary_search_and_effective_weight_of_snell_law_path_count;
            total_snell_law_memory_size += snell_law_memory_size;
            edge_sequence.pop_back();
            face_sequence.pop_back();
            source_x = next_source_x;
            source_y = next_source_y;
            source_z = next_source_z;
        }
    }

    for (int i = 1; i < snell_law_path.size() - 1; i++)
    {
        result_path.push_back(snell_law_path[i]);
    }
    result_path.push_back(destination);
}
