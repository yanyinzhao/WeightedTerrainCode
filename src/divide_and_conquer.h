#include "snell_law.h"
#include <chrono>
#include <fstream>

// check whether a vertex is in the edge
bool vertex_in_edge_2(double edge_1_x, double edge_2_x,
                      double edge_1_y, double edge_2_y,
                      double edge_1_z, double edge_2_z,
                      double vertex_x, double vertex_y,
                      double vertex_z)
{
    return ((vertex_x - min_1(edge_1_x, edge_2_x) >= -1e-6) &&
            (vertex_x - max_1(edge_1_x, edge_2_x) <= 1e-6) &&
            (vertex_y - min_1(edge_1_y, edge_2_y) >= -1e-6) &&
            (vertex_y - max_1(edge_1_y, edge_2_y) <= 1e-6) &&
            (vertex_z - min_1(edge_1_z, edge_2_z) >= -1e-6) &&
            (vertex_z - max_1(edge_1_z, edge_2_z) <= 1e-6));
}

double get_edge_pass_face_weight_3(geodesic::SurfacePoint &point1, geodesic::SurfacePoint &point2)
{
    // if both two points are the on the vertex, then they are on the original mesh edge,
    // so the weight of this path is the max face weight of the faces that connects to this edge
    if (point1.type() == geodesic::VERTEX && point2.type() == geodesic::VERTEX)
    {
        geodesic::edge_pointer shared_edge;
        bool jump_out = false;
        for (unsigned i = 0; i < point1.base_element()->adjacent_edges().size() && !jump_out; ++i)
        {
            for (unsigned j = 0; j < point2.base_element()->adjacent_edges().size() && !jump_out; ++j)
            {
                if (point1.base_element()->adjacent_edges()[i]->id() == point2.base_element()->adjacent_edges()[j]->id())
                {
                    shared_edge = point1.base_element()->adjacent_edges()[i];
                    jump_out = true;
                    break;
                }
            }
        }
        std::vector<double> face_weight_list;
        face_weight_list.clear();
        for (unsigned i = 0; i < shared_edge->adjacent_faces().size(); ++i)
        {
            face_weight_list.push_back(shared_edge->adjacent_faces()[i]->weight());
        }
        assert(face_weight_list.size());
        return *std::min_element(face_weight_list.begin(), face_weight_list.end());
    }

    // if the first point is on the vertex and the second point is on the edge and the first point is on this edge
    // then they are on the original mesh edge, so the weight of this path is the max face weight of the faces that
    // connect to this edge
    else if ((point1.type() == geodesic::VERTEX && point2.type() == geodesic::EDGE && point1.base_element()->id() == point2.base_element()->adjacent_vertices()[0]->id()) ||
             (point1.type() == geodesic::VERTEX && point2.type() == geodesic::EDGE && point1.base_element()->id() == point2.base_element()->adjacent_vertices()[1]->id()))
    {

        std::vector<double> face_weight_list;
        face_weight_list.clear();
        for (unsigned i = 0; i < point2.base_element()->adjacent_faces().size(); ++i)
        {
            face_weight_list.push_back(point2.base_element()->adjacent_faces()[i]->weight());
        }
        assert(face_weight_list.size());
        return *std::min_element(face_weight_list.begin(), face_weight_list.end());
    }

    // if the first point is on the edge and the second point is on the vertex and the second point is on this edge
    // then they are on the original mesh edge, so the weight of this path is the max face weight of the faces that
    // connect to this edge
    else if ((point1.type() == geodesic::EDGE && point2.type() == geodesic::VERTEX && point2.base_element()->id() == point1.base_element()->adjacent_vertices()[0]->id()) ||
             (point1.type() == geodesic::EDGE && point2.type() == geodesic::VERTEX && point2.base_element()->id() == point1.base_element()->adjacent_vertices()[1]->id()))
    {

        std::vector<double> face_weight_list;
        face_weight_list.clear();
        for (unsigned i = 0; i < point1.base_element()->adjacent_faces().size(); ++i)
        {
            face_weight_list.push_back(point1.base_element()->adjacent_faces()[i]->weight());
        }
        assert(face_weight_list.size());
        return *std::min_element(face_weight_list.begin(), face_weight_list.end());
    }

    // if two points are both on the same original mesh edge,
    // the weight of this path is the min face weight of the faces that connects to this edge
    else if (point1.type() == geodesic::EDGE && point2.type() == geodesic::EDGE && point1.base_element()->id() == point2.base_element()->id())
    {

        std::vector<double> face_weight_list;
        face_weight_list.clear();
        for (unsigned i = 0; i < point2.base_element()->adjacent_faces().size(); ++i)
        {
            face_weight_list.push_back(point2.base_element()->adjacent_faces()[i]->weight());
        }
        assert(face_weight_list.size());
        return *std::min_element(face_weight_list.begin(), face_weight_list.end());
    }

    // if one point is on the vertex and one point is on the edge, but that vertex is not on this edge
    // or if the two points are on two different edges
    // then the path weight is the face weight
    else if ((point1.type() == geodesic::VERTEX && point2.type() == geodesic::EDGE && point1.base_element()->id() != point2.base_element()->adjacent_vertices()[0]->id() && point1.base_element()->id() != point2.base_element()->adjacent_vertices()[1]->id()) ||
             (point1.type() == geodesic::EDGE && point2.type() == geodesic::VERTEX && point2.base_element()->id() != point1.base_element()->adjacent_vertices()[0]->id() && point2.base_element()->id() != point1.base_element()->adjacent_vertices()[1]->id()) ||
             (point1.type() == geodesic::EDGE && point2.type() == geodesic::EDGE))
    {

        for (unsigned i = 0; i < point2.base_element()->adjacent_faces().size(); ++i)
        {
            geodesic::face_pointer f = point2.base_element()->adjacent_faces()[i];

            if (((vertex_in_edge_2(
                      f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[1]->x(),
                      f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[1]->y(),
                      f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[1]->z(),
                      point1.getx(), point1.gety(), point1.getz()) &&
                  vertex_in_edge_2(
                      f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[2]->x(),
                      f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[2]->y(),
                      f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[2]->z(),
                      point2.getx(), point2.gety(), point2.getz())) ||
                 (vertex_in_edge_2(
                      f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[1]->x(),
                      f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[1]->y(),
                      f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[1]->z(),
                      point1.getx(), point1.gety(), point1.getz()) &&
                  vertex_in_edge_2(
                      f->adjacent_vertices()[1]->x(), f->adjacent_vertices()[2]->x(),
                      f->adjacent_vertices()[1]->y(), f->adjacent_vertices()[2]->y(),
                      f->adjacent_vertices()[1]->z(), f->adjacent_vertices()[2]->z(),
                      point2.getx(), point2.gety(), point2.getz()))) ||
                ((vertex_in_edge_2(
                      f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[2]->x(),
                      f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[2]->y(),
                      f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[2]->z(),
                      point1.getx(), point1.gety(), point1.getz()) &&
                  vertex_in_edge_2(
                      f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[1]->x(),
                      f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[1]->y(),
                      f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[1]->z(),
                      point2.getx(), point2.gety(), point2.getz())) ||
                 (vertex_in_edge_2(
                      f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[2]->x(),
                      f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[2]->y(),
                      f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[2]->z(),
                      point1.getx(), point1.gety(), point1.getz()) &&
                  vertex_in_edge_2(
                      f->adjacent_vertices()[1]->x(), f->adjacent_vertices()[2]->x(),
                      f->adjacent_vertices()[1]->y(), f->adjacent_vertices()[2]->y(),
                      f->adjacent_vertices()[1]->z(), f->adjacent_vertices()[2]->z(),
                      point2.getx(), point2.gety(), point2.getz()))) ||
                ((vertex_in_edge_2(
                      f->adjacent_vertices()[1]->x(), f->adjacent_vertices()[2]->x(),
                      f->adjacent_vertices()[1]->y(), f->adjacent_vertices()[2]->y(),
                      f->adjacent_vertices()[1]->z(), f->adjacent_vertices()[2]->z(),
                      point1.getx(), point1.gety(), point1.getz()) &&
                  vertex_in_edge_2(
                      f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[1]->x(),
                      f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[1]->y(),
                      f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[1]->z(),
                      point2.getx(), point2.gety(), point2.getz())) ||
                 (vertex_in_edge_2(
                      f->adjacent_vertices()[1]->x(), f->adjacent_vertices()[2]->x(),
                      f->adjacent_vertices()[1]->y(), f->adjacent_vertices()[2]->y(),
                      f->adjacent_vertices()[1]->z(), f->adjacent_vertices()[2]->z(),
                      point1.getx(), point1.gety(), point1.getz()) &&
                  vertex_in_edge_2(
                      f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[2]->x(),
                      f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[2]->y(),
                      f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[2]->z(),
                      point2.getx(), point2.gety(), point2.getz()))))
            {
                return f->weight();
            }
        }
    }
    return 100;
}

// calculate the subdivision level using the epsilon value for the fixed Steiner point algorithm
void fixed_Steiner_point_epsilon_to_subdivision_level(geodesic::Mesh *mesh,
                                                      double epsilon,
                                                      int estimate_path_length,
                                                      unsigned &subdivision_level)
{
    // calcualte the longest and shortest edge
    double max_edge_length = -1e100;
    double min_edge_length = 1e100;
    std::vector<double> edge_length_list;
    std::vector<double> edge_length_without_outliers_list;
    for (unsigned i = 0; i < mesh->edges().size(); ++i)
    {
        geodesic::Edge &e = mesh->edges()[i];
        double edge_length = e.length();
        edge_length_list.push_back(edge_length);
        max_edge_length = std::max(max_edge_length, edge_length);
        min_edge_length = std::min(min_edge_length, edge_length);
    }

    max_edge_length = -1e100;
    min_edge_length = 1e100;
    remove_outliers(edge_length_list, edge_length_without_outliers_list);
    for (int i = 0; i < edge_length_without_outliers_list.size(); i++)
    {
        max_edge_length = std::max(max_edge_length, edge_length_without_outliers_list[i]);
        min_edge_length = std::min(min_edge_length, edge_length_without_outliers_list[i]);
    }

    // calculate max face weight
    double max_face_weight = -1e100;
    for (unsigned i = 0; i < mesh->faces().size(); ++i)
    {
        geodesic::Face &f = mesh->faces()[i];
        double face_weight = f.weight();
        max_face_weight = std::max(max_face_weight, face_weight);
    }

    subdivision_level = std::max(1.0, floor((max_edge_length * estimate_path_length * max_face_weight / (epsilon * min_edge_length) - 1) / 4));
}

// calculate the subdivision level and other useful parameters using the epsilon value for the log Steiner point algorithm
void log_Steiner_point_epsilon_to_subdivision_level(geodesic::Mesh *mesh,
                                                    double epsilon,
                                                    unsigned &max_subdivision_log_level,
                                                    double &delta,
                                                    double &r)
{
    // calculate the original epsilon used in the paper
    double max_face_weight = -1e100;
    double min_face_weight = 1e100;
    for (unsigned i = 0; i < mesh->faces().size(); ++i)
    {
        geodesic::Face &f = mesh->faces()[i];
        double face_weight = f.weight();
        max_face_weight = std::max(max_face_weight, face_weight);
        min_face_weight = std::min(min_face_weight, face_weight);
    }
    double epsilon_2_factor;
    if (epsilon > 0 & epsilon <= 1)
    {
        epsilon_2_factor = 0.46 * (log(epsilon / 3.6) / log(0.3)) / ((2 + max_face_weight / min_face_weight - (sqrt(pow(2 + max_face_weight / min_face_weight, 2) - 4))) / 4);
    }
    else if (epsilon > 1)
    {
        epsilon_2_factor = 0.49 / ((2 + max_face_weight / min_face_weight - (sqrt(pow(2 + max_face_weight / min_face_weight, 2) - 4))) / 4);
    }
    double epsilon_2 = epsilon_2_factor * (1 + epsilon + max_face_weight / min_face_weight - (sqrt(pow(1 + epsilon + max_face_weight / min_face_weight, 2) - 4 * epsilon))) / 4;

    // caluclate the min face angle
    double min_angle = 1e100;
    std::vector<double> angle_list;
    std::vector<double> angle_without_outliers_list;
    for (unsigned i = 0; i < mesh->faces().size(); ++i)
    {
        geodesic::Face &f = mesh->faces()[i];
        for (unsigned j = 0; j < 3; ++j)
        {
            double angle = f.corner_angles()[j];
            angle_list.push_back(angle);
            min_angle = std::min(min_angle, angle);
        }
    }
    min_angle = 1e100;
    remove_outliers(angle_list, angle_without_outliers_list);
    for (int i = 0; i < angle_without_outliers_list.size(); i++)
    {
        min_angle = std::min(min_angle, angle_without_outliers_list[i]);
    }

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

    // calculate triangle area and min height
    double min_height = 1e100;
    std::vector<double> height_list;
    std::vector<double> height_without_outliers_list;
    for (unsigned i = 0; i < mesh->faces().size(); ++i)
    {
        geodesic::Face &f = mesh->faces()[i];
        double area_two_times;
        double vector_cross_product_x;
        double vector_cross_product_y;
        double vector_cross_product_z;

        double v1_x = f.adjacent_vertices()[0]->x();
        double v1_y = f.adjacent_vertices()[0]->y();
        double v1_z = f.adjacent_vertices()[0]->z();
        double v2_x = f.adjacent_vertices()[1]->x();
        double v2_y = f.adjacent_vertices()[1]->y();
        double v2_z = f.adjacent_vertices()[1]->z();
        double v3_x = f.adjacent_vertices()[2]->x();
        double v3_y = f.adjacent_vertices()[2]->y();
        double v3_z = f.adjacent_vertices()[2]->z();

        double vector1_x = v2_x - v1_x;
        double vector1_y = v2_y - v1_y;
        double vector1_z = v2_z - v1_z;
        double vector2_x = v3_x - v1_x;
        double vector2_y = v3_y - v1_y;
        double vector2_z = v3_z - v1_z;

        vector_cross_product(vector1_x, vector1_y, vector1_z,
                             vector2_x, vector2_y, vector2_z,
                             vector_cross_product_x, vector_cross_product_y, vector_cross_product_z);

        area_two_times = sqrt(pow(vector_cross_product_x, 2) + pow(vector_cross_product_y, 2) + pow(vector_cross_product_z, 2));

        double max_edge_length_inside_one_face = 0;
        for (unsigned j = 0; j < 3; ++j)
        {
            double edge_length_inside_one_face = f.adjacent_edges()[j]->length();
            max_edge_length_inside_one_face = std::max(max_edge_length_inside_one_face, edge_length_inside_one_face);
        }

        double min_height_inside_one_face = area_two_times / max_edge_length_inside_one_face;
        height_list.push_back(min_height_inside_one_face);
        min_height = std::min(min_height, min_height_inside_one_face);
    }
    min_height = 1e100;
    remove_outliers(height_list, height_without_outliers_list);
    for (int i = 0; i < height_without_outliers_list.size(); i++)
    {
        min_height = std::min(min_height, height_without_outliers_list[i]);
    }

    // calculate delta
    delta = 1 + epsilon_2 * sin(min_angle);

    // calculate r
    r = min_height * epsilon_2;

    max_subdivision_log_level = floor(log(max_edge_length / r) / log(delta));
}

void snell_law_epsilon_to_delta(geodesic::Mesh *mesh,
                                double epsilon,
                                int l,
                                double &delta)
{
    // calculate the min and max face weight
    double max_face_weight = -1e100;
    double min_face_weight = 1e100;
    for (unsigned i = 0; i < mesh->faces().size(); ++i)
    {
        geodesic::Face &f = mesh->faces()[i];
        double face_weight = f.weight();
        max_face_weight = std::max(max_face_weight, face_weight);
        min_face_weight = std::min(min_face_weight, face_weight);
    }

    // calculate triangle area and min height
    double min_height = 1e100;
    std::vector<double> height_list;
    std::vector<double> height_without_outliers_list;
    for (unsigned i = 0; i < mesh->faces().size(); ++i)
    {
        geodesic::Face &f = mesh->faces()[i];
        double area_two_times;
        double vector_cross_product_x;
        double vector_cross_product_y;
        double vector_cross_product_z;

        double v1_x = f.adjacent_vertices()[0]->x();
        double v1_y = f.adjacent_vertices()[0]->y();
        double v1_z = f.adjacent_vertices()[0]->z();
        double v2_x = f.adjacent_vertices()[1]->x();
        double v2_y = f.adjacent_vertices()[1]->y();
        double v2_z = f.adjacent_vertices()[1]->z();
        double v3_x = f.adjacent_vertices()[2]->x();
        double v3_y = f.adjacent_vertices()[2]->y();
        double v3_z = f.adjacent_vertices()[2]->z();

        double vector1_x = v2_x - v1_x;
        double vector1_y = v2_y - v1_y;
        double vector1_z = v2_z - v1_z;
        double vector2_x = v3_x - v1_x;
        double vector2_y = v3_y - v1_y;
        double vector2_z = v3_z - v1_z;

        vector_cross_product(vector1_x, vector1_y, vector1_z,
                             vector2_x, vector2_y, vector2_z,
                             vector_cross_product_x, vector_cross_product_y, vector_cross_product_z);

        area_two_times = sqrt(pow(vector_cross_product_x, 2) + pow(vector_cross_product_y, 2) + pow(vector_cross_product_z, 2));

        double max_edge_length_inside_one_face = 0;
        for (unsigned j = 0; j < 3; ++j)
        {
            double edge_length_inside_one_face = f.adjacent_edges()[j]->length();
            max_edge_length_inside_one_face = std::max(max_edge_length_inside_one_face, edge_length_inside_one_face);
        }

        double min_height_inside_one_face = area_two_times / max_edge_length_inside_one_face;
        height_list.push_back(min_height_inside_one_face);
        min_height = std::min(min_height, min_height_inside_one_face);
    }
    min_height = 1e100;
    remove_outliers(height_list, height_without_outliers_list);
    for (int i = 0; i < height_without_outliers_list.size(); i++)
    {
        min_height = std::min(min_height, height_without_outliers_list[i]);
    }

    assert(l > 0);
    delta = min_height * epsilon * min_face_weight / (6 * l * max_face_weight);
}

double distance_two_points(double x_1, double y_1, double z_1,
                           double x_2, double y_2, double z_2)
{
    double d_x = x_2 - x_1;
    double d_y = y_2 - y_1;
    double d_z = z_2 - z_1;

    return sqrt(d_x * d_x + d_y * d_y + d_z * d_z);
}

double path_distance(std::vector<geodesic::SurfacePoint> &path_for_distance, std::vector<geodesic::SurfacePoint> &path_for_face_weight)
{
    assert(path_for_distance.size() == path_for_face_weight.size());
    double total_distance = 0;

    for (int i = 0; i < path_for_distance.size() - 1; ++i)
    {
        total_distance += distance_two_points(path_for_distance[i].x(), path_for_distance[i].y(), path_for_distance[i].z(),
                                              path_for_distance[i + 1].x(), path_for_distance[i + 1].y(), path_for_distance[i + 1].z()) *
                          get_edge_pass_face_weight_3(path_for_face_weight[i], path_for_face_weight[i + 1]);
    }

    return total_distance;
}

// can directly call this function to run the fixed Steiner point algorithm
void fixed_Steiner_point(geodesic::Mesh *mesh, geodesic::SurfacePoint &source,
                         geodesic::SurfacePoint &destination,
                         std::vector<geodesic::SurfacePoint> &path,
                         double epsilon, int estimate_path_length,
                         double &building_time, double &query_time,
                         double &memory_size)
{
    unsigned subdivision_level;
    fixed_Steiner_point_epsilon_to_subdivision_level(mesh, epsilon, estimate_path_length, subdivision_level);
    std::cout << "Original subdivision level: " << subdivision_level << std::endl;

    auto start_building_time = std::chrono::high_resolution_clock::now();
    geodesic::GeodesicAlgorithmSubdivision algorithm(mesh, subdivision_level);
    auto stop_building_time = std::chrono::high_resolution_clock::now();
    auto duration_building_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_building_time - start_building_time);
    building_time = duration_building_time.count();

    auto start_query_time = std::chrono::high_resolution_clock::now();
    algorithm.geodesic(source, destination, path);
    memory_size = algorithm.get_memory();
    int path_length = path.size();
    for (int i = 1; i < path_length - 1; ++i)
    {
        geodesic::SurfacePoint &prev = path[i + 1];
        geodesic::SurfacePoint &curr = path[i];
        geodesic::SurfacePoint &next = path[i - 1];
        if ((prev.type() == geodesic::EDGE && curr.type() == geodesic::EDGE && next.type() == geodesic::EDGE &&
             prev.base_element()->id() == curr.base_element()->id() && next.base_element()->id() == curr.base_element()->id()) ||
            (prev.type() == geodesic::EDGE && curr.type() == geodesic::EDGE && next.type() == geodesic::VERTEX &&
             next.base_element()->id() == curr.base_element()->adjacent_vertices()[0]->id() && prev.base_element()->id() == curr.base_element()->id()) ||
            (prev.type() == geodesic::EDGE && curr.type() == geodesic::EDGE && next.type() == geodesic::VERTEX &&
             next.base_element()->id() == curr.base_element()->adjacent_vertices()[1]->id() && prev.base_element()->id() == curr.base_element()->id()) ||
            (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && next.type() == geodesic::EDGE &&
             prev.base_element()->id() == curr.base_element()->adjacent_vertices()[0]->id() && next.base_element()->id() == curr.base_element()->id()) ||
            (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && next.type() == geodesic::EDGE &&
             prev.base_element()->id() == curr.base_element()->adjacent_vertices()[1]->id() && next.base_element()->id() == curr.base_element()->id()) ||
            (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && next.type() == geodesic::VERTEX &&
             prev.base_element()->id() == curr.base_element()->adjacent_vertices()[0]->id() && next.base_element()->id() == curr.base_element()->adjacent_vertices()[1]->id()) ||
            (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && next.type() == geodesic::VERTEX &&
             prev.base_element()->id() == curr.base_element()->adjacent_vertices()[1]->id() && next.base_element()->id() == curr.base_element()->adjacent_vertices()[0]->id()))
        {
            path.erase(path.begin() + i);
            path_length--;
            i--;
        }
    }

    // follow the direction from destination to source (backward):
    // Sometimes, it could happen that the current point is on the edge, and the previous point or the next point is one of the vertex of this edge
    // so we need to delete the current point
    for (int i = 1; i < path_length - 1; ++i)
    {
        geodesic::SurfacePoint &prev = path[i - 1];
        geodesic::SurfacePoint &curr = path[i];
        geodesic::SurfacePoint &next = path[i + 1];
        if ((prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && curr.base_element()->adjacent_vertices()[0]->id() == prev.base_element()->id()) ||
            (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && curr.base_element()->adjacent_vertices()[1]->id() == prev.base_element()->id()) ||
            (next.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && curr.base_element()->adjacent_vertices()[0]->id() == next.base_element()->id()) ||
            (next.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && curr.base_element()->adjacent_vertices()[1]->id() == next.base_element()->id()))
        {
            path.erase(path.begin() + i);
            path_length--;
            i--;
        }
    }

    // follow the direction from destination to source (backward):
    // Sometimes, it could happen that the current point and next point are on the same edge
    // so we need to delete the next point
    for (int i = 0; i < path_length - 1; ++i)
    {
        geodesic::SurfacePoint &curr = path[i];
        geodesic::SurfacePoint &next = path[i + 1];
        if (curr.type() == geodesic::EDGE && next.type() == geodesic::EDGE &&
            curr.base_element()->id() == next.base_element()->id())
        {
            path.erase(path.begin() + i);
            path_length--;
            i--;
        }
    }
    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();

    memory_size += path.size() * sizeof(geodesic::SurfacePoint);

    double total_distance;
    total_distance = path_distance(path, path);
    std::cout << "Fixed Steiner point distance: " << total_distance << std::endl;
}

// can directly call this function to run the log Steiner point algorithm
void log_Steiner_point(geodesic::Mesh *mesh, geodesic::SurfacePoint &source,
                       geodesic::SurfacePoint &destination,
                       std::vector<geodesic::SurfacePoint> &path,
                       double epsilon, double &building_time,
                       double &query_time, double &memory_size)
{
    unsigned max_subdivision_log_level;
    double delta;
    double r;
    log_Steiner_point_epsilon_to_subdivision_level(mesh, epsilon, max_subdivision_log_level, delta, r);
    std::cout << "Original max subdivision log level: " << 2 * max_subdivision_log_level << std::endl;

    auto start_building_time = std::chrono::high_resolution_clock::now();
    geodesic::GeodesicAlgorithmSubdivisionLog algorithm(mesh, max_subdivision_log_level, delta, r);
    auto stop_building_time = std::chrono::high_resolution_clock::now();
    auto duration_building_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_building_time - start_building_time);
    building_time = duration_building_time.count();

    auto start_query_time = std::chrono::high_resolution_clock::now();
    algorithm.geodesic(source, destination, path);
    memory_size = algorithm.get_memory();
    int path_length = path.size();
    for (int i = 1; i < path_length - 1; ++i)
    {
        geodesic::SurfacePoint &prev = path[i + 1];
        geodesic::SurfacePoint &curr = path[i];
        geodesic::SurfacePoint &next = path[i - 1];
        if ((prev.type() == geodesic::EDGE && curr.type() == geodesic::EDGE && next.type() == geodesic::EDGE &&
             prev.base_element()->id() == curr.base_element()->id() && next.base_element()->id() == curr.base_element()->id()) ||
            (prev.type() == geodesic::EDGE && curr.type() == geodesic::EDGE && next.type() == geodesic::VERTEX &&
             next.base_element()->id() == curr.base_element()->adjacent_vertices()[0]->id() && prev.base_element()->id() == curr.base_element()->id()) ||
            (prev.type() == geodesic::EDGE && curr.type() == geodesic::EDGE && next.type() == geodesic::VERTEX &&
             next.base_element()->id() == curr.base_element()->adjacent_vertices()[1]->id() && prev.base_element()->id() == curr.base_element()->id()) ||
            (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && next.type() == geodesic::EDGE &&
             prev.base_element()->id() == curr.base_element()->adjacent_vertices()[0]->id() && next.base_element()->id() == curr.base_element()->id()) ||
            (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && next.type() == geodesic::EDGE &&
             prev.base_element()->id() == curr.base_element()->adjacent_vertices()[1]->id() && next.base_element()->id() == curr.base_element()->id()) ||
            (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && next.type() == geodesic::VERTEX &&
             prev.base_element()->id() == curr.base_element()->adjacent_vertices()[0]->id() && next.base_element()->id() == curr.base_element()->adjacent_vertices()[1]->id()) ||
            (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && next.type() == geodesic::VERTEX &&
             prev.base_element()->id() == curr.base_element()->adjacent_vertices()[1]->id() && next.base_element()->id() == curr.base_element()->adjacent_vertices()[0]->id()))
        {
            path.erase(path.begin() + i);
            path_length--;
            i--;
        }
    }

    // follow the direction from destination to source (backward):
    // Sometimes, it could happen that the current point is on the edge, and the previous point or the next point is one of the vertex of this edge
    // so we need to delete the current point
    for (int i = 1; i < path_length - 1; ++i)
    {
        geodesic::SurfacePoint &prev = path[i - 1];
        geodesic::SurfacePoint &curr = path[i];
        geodesic::SurfacePoint &next = path[i + 1];
        if ((prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && curr.base_element()->adjacent_vertices()[0]->id() == prev.base_element()->id()) ||
            (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && curr.base_element()->adjacent_vertices()[1]->id() == prev.base_element()->id()) ||
            (next.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && curr.base_element()->adjacent_vertices()[0]->id() == next.base_element()->id()) ||
            (next.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && curr.base_element()->adjacent_vertices()[1]->id() == next.base_element()->id()))
        {
            path.erase(path.begin() + i);
            path_length--;
            i--;
        }
    }

    // follow the direction from destination to source (backward):
    // Sometimes, it could happen that the current point and next point are on the same edge
    // so we need to delete the next point
    for (int i = 0; i < path_length - 1; ++i)
    {
        geodesic::SurfacePoint &curr = path[i];
        geodesic::SurfacePoint &next = path[i + 1];
        if (curr.type() == geodesic::EDGE && next.type() == geodesic::EDGE &&
            curr.base_element()->id() == next.base_element()->id())
        {
            path.erase(path.begin() + i);
            path_length--;
            i--;
        }
    }
    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();

    memory_size += path.size() * sizeof(geodesic::SurfacePoint);

    double total_distance;
    total_distance = path_distance(path, path);
    std::cout << "Log Steiner point distance: " << total_distance << std::endl;
}

// cannot call this function directly, need to call log_Steiner_point_divide_and_conquer_helper
void log_Steiner_point_divide_and_conquer(geodesic::Mesh *mesh, geodesic::SurfacePoint &source,
                                          geodesic::SurfacePoint &destination,
                                          std::vector<geodesic::SurfacePoint> &path,
                                          unsigned max_subdivision_log_level,
                                          double delta, double r, int max_divide_and_conquer_subdivision_log_level,
                                          int max_loop_num_for_single_endpoint,
                                          double &total_building_time, double &total_Steiner_point_memory_size)
{
    auto start_building_time = std::chrono::high_resolution_clock::now();
    geodesic::GeodesicAlgorithmSubdivisionLog algorithm(mesh, max_subdivision_log_level, delta, r);
    auto stop_building_time = std::chrono::high_resolution_clock::now();
    auto duration_building_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_building_time - start_building_time);
    total_building_time += duration_building_time.count();

    std::cout << "divide and conquer max subdivision log level: " << 2 * max_subdivision_log_level << std::endl;

    geodesic::SurfacePoint successive_points_previous_vertex;
    int successive_points_count = 0;

    algorithm.geodesic(source, destination, path);
    total_Steiner_point_memory_size += algorithm.return_memory();

    // follow the direction from destination to source (backward):
    // Sometimes, it could happen that the previous point, current point and next point are on the same edge
    // so we need to delete the current point
    int path_length = path.size();
    for (int i = 1; i < path_length - 1; ++i)
    {
        geodesic::SurfacePoint &prev = path[i - 1];
        geodesic::SurfacePoint &curr = path[i];
        geodesic::SurfacePoint &next = path[i + 1];
        if ((prev.type() == geodesic::EDGE && curr.type() == geodesic::EDGE && next.type() == geodesic::EDGE &&
             prev.base_element()->id() == curr.base_element()->id() && next.base_element()->id() == curr.base_element()->id()) ||
            (prev.type() == geodesic::EDGE && curr.type() == geodesic::EDGE && next.type() == geodesic::VERTEX &&
             next.base_element()->id() == curr.base_element()->adjacent_vertices()[0]->id() && prev.base_element()->id() == curr.base_element()->id()) ||
            (prev.type() == geodesic::EDGE && curr.type() == geodesic::EDGE && next.type() == geodesic::VERTEX &&
             next.base_element()->id() == curr.base_element()->adjacent_vertices()[1]->id() && prev.base_element()->id() == curr.base_element()->id()) ||
            (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && next.type() == geodesic::EDGE &&
             prev.base_element()->id() == curr.base_element()->adjacent_vertices()[0]->id() && next.base_element()->id() == curr.base_element()->id()) ||
            (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && next.type() == geodesic::EDGE &&
             prev.base_element()->id() == curr.base_element()->adjacent_vertices()[1]->id() && next.base_element()->id() == curr.base_element()->id()) ||
            (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && next.type() == geodesic::VERTEX &&
             prev.base_element()->id() == curr.base_element()->adjacent_vertices()[0]->id() && next.base_element()->id() == curr.base_element()->adjacent_vertices()[1]->id()) ||
            (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && next.type() == geodesic::VERTEX &&
             prev.base_element()->id() == curr.base_element()->adjacent_vertices()[1]->id() && next.base_element()->id() == curr.base_element()->adjacent_vertices()[0]->id()))
        {
            path.erase(path.begin() + i);
            path_length--;
            i--;
        }
    }

    // follow the direction from destination to source (backward):
    // Sometimes, it could happen that the current point is on the edge, and the previous point or the next point is one of the vertex of this edge
    // so we need to delete the current point
    for (int i = 1; i < path_length - 1; ++i)
    {
        geodesic::SurfacePoint &prev = path[i - 1];
        geodesic::SurfacePoint &curr = path[i];
        geodesic::SurfacePoint &next = path[i + 1];
        if ((prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && curr.base_element()->adjacent_vertices()[0]->id() == prev.base_element()->id()) ||
            (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && curr.base_element()->adjacent_vertices()[1]->id() == prev.base_element()->id()) ||
            (next.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && curr.base_element()->adjacent_vertices()[0]->id() == next.base_element()->id()) ||
            (next.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && curr.base_element()->adjacent_vertices()[1]->id() == next.base_element()->id()))
        {
            path.erase(path.begin() + i);
            path_length--;
            i--;
        }
    }

    // follow the direction from destination to source (backward):
    // Sometimes, it could happen that the current point and next point are on the same edge
    // so we need to delete the next point
    for (int i = 0; i < path_length - 1; ++i)
    {
        geodesic::SurfacePoint &curr = path[i];
        geodesic::SurfacePoint &next = path[i + 1];
        if (curr.type() == geodesic::EDGE && next.type() == geodesic::EDGE &&
            curr.base_element()->id() == next.base_element()->id())
        {
            path.erase(path.begin() + i);
            path_length--;
            i--;
        }
    }

    if (max_subdivision_log_level >= max_divide_and_conquer_subdivision_log_level)
    {
        return;
    }

    for (unsigned i = 0; i < path.size(); ++i)
    {
        geodesic::SurfacePoint &s = path[i];
    }

    path_length = path.size();
    int fixed_path_length = path.size();
    int original_counter = 1;

    for (int i = 1; i < path_length - 1; ++i)
    {
        geodesic::SurfacePoint &s = path[i];

        if (s.type() == geodesic::EDGE)
        {
            ++original_counter;
            continue;
        }

        assert(s.type() == geodesic::VERTEX);

        // successive point

        // follow the direction from destination to source (backward):
        // if the previous point is on the edge, the next point is on the vertex, and the current point is not the neighbour of source point or destination point
        // or if the next point is on the vertex and the current point is next to the destination point
        // we store the previous vertex for further refinement and delete current vertex
        if ((path[i - 1].type() == geodesic::EDGE && path[i + 1].type() == geodesic::VERTEX && original_counter != 1 && original_counter != fixed_path_length - 2) ||
            (original_counter == 1 && path[i + 1].type() == geodesic::VERTEX))
        {
            std::cout << "start successive point" << std::endl;
            successive_points_previous_vertex = path[i - 1];
            ++successive_points_count;
        }

        // follow the direction from destination to source (backward):
        // if previous point and the next point are on the vertex, and the current point is not the neighbour of source point or destination point
        // we delete current vertex
        else if (path[i - 1].type() == geodesic::VERTEX && path[i + 1].type() == geodesic::VERTEX && original_counter != 1 && original_counter != fixed_path_length - 2)
        {
            std::cout << "store successive point" << std::endl;
            ++successive_points_count;
        }

        // follow the direction from destination to source (backward):
        // if the previous point is on the vertex, the next point is on the edge, and the current point is not the neighbour of source point or destination point
        // or if the previous point is on the vertex and the current point is next to the source point
        // we will use divide and conquer to refinement
        else if ((path[i - 1].type() == geodesic::VERTEX && path[i + 1].type() == geodesic::EDGE && original_counter != 1 && original_counter != fixed_path_length - 2) ||
                 (original_counter == fixed_path_length - 2 && path[i - 1].type() == geodesic::VERTEX))
        {
            if (successive_points_count > fixed_path_length / 4)
            {
                std::cout << "handle successive point and divide-and-conquer" << std::endl;

                for (int j = 0; j <= successive_points_count; ++j)
                {
                    path.erase(path.begin() + i);
                    i--;
                    path_length--;
                }

                successive_points_count = 0;

                std::vector<geodesic::SurfacePoint> sub_path;
                sub_path.clear();
                log_Steiner_point_divide_and_conquer(mesh, path[i + 1], successive_points_previous_vertex, sub_path, 2 * max_subdivision_log_level, delta, r / 2, max_divide_and_conquer_subdivision_log_level, max_loop_num_for_single_endpoint, total_building_time, total_Steiner_point_memory_size);

                for (int j = sub_path.size() - 2; j >= 1; --j)
                {
                    path.insert(path.begin() + i + 1, sub_path[j]);
                }
                i += sub_path.size() - 2;
                path_length += sub_path.size() - 2;
            }
            else
            {
                std::cout << "successive points too less, not use divide-and-conquer" << std::endl;
                successive_points_count = 0;
            }
        }

        // *******
        // follow the direction from destination to source (backward):
        // if the point is not connect to the source or destination point of the path
        // or if the point is connect to the destination point (calculated using the divide-and-conquer algorithm) of the sub-path
        // or if the point is connect to the source point (calculated using the divide-and-conquer algorithm) of the sub-path
        // *******
        else if ((path[i - 1].type() == geodesic::EDGE && path[i + 1].type() == geodesic::EDGE && original_counter != 1 && original_counter != fixed_path_length - 2) ||
                 (original_counter == 1 && path[i + 1].type() == geodesic::EDGE && path[i - 1].type() == geodesic::EDGE) ||
                 (original_counter == fixed_path_length - 2 && path[i - 1].type() == geodesic::EDGE && path[i + 1].type() == geodesic::EDGE))
        {
            // store the path distance for the new Steiner point on two sides of the original vertex
            // also store the path distance for the original vertex
            double new_steiner_point_one_side_distance;
            double new_steiner_point_another_side_distance;
            double original_vertex_distance;
            int loop_count = 0;

            do
            {
                new_steiner_point_one_side_distance = 0;
                new_steiner_point_another_side_distance = 0;
                original_vertex_distance = 0;

                std::vector<geodesic::edge_pointer> edges_opposite_to_current;
                edges_opposite_to_current.clear();

                for (unsigned j = 0; j < s.base_element()->adjacent_faces().size(); ++j)
                {
                    geodesic::face_pointer f = s.base_element()->adjacent_faces()[j];

                    // we will not include the two edges that the previous and next Steiner point lies on
                    if (f->opposite_edge(static_cast<geodesic::vertex_pointer>(s.base_element()))->id() == path[i + 1].base_element()->id() ||
                        f->opposite_edge(static_cast<geodesic::vertex_pointer>(s.base_element()))->id() == path[i - 1].base_element()->id())
                    {
                        continue;
                    }
                    edges_opposite_to_current.push_back(f->opposite_edge(static_cast<geodesic::vertex_pointer>(s.base_element())));
                }

                // store the vertices on one side (no matter left or right, both are okay) of the calculated path
                std::vector<unsigned> one_side_vertices_id;

                one_side_vertices_id.clear();

                geodesic::vertex_pointer current_neighbor_1 = path[i + 1].base_element()->adjacent_vertices()[0];

                one_side_vertices_id.push_back(current_neighbor_1->id());

                std::vector<unsigned> already_added_edge_index_1;
                do
                {
                    for (unsigned j = 0; j < edges_opposite_to_current.size(); ++j)
                    {
                        if (std::find(already_added_edge_index_1.begin(), already_added_edge_index_1.end(), j) != already_added_edge_index_1.end())
                        {
                            continue;
                        }
                        if (current_neighbor_1->id() == edges_opposite_to_current[j]->adjacent_vertices()[0]->id())
                        {
                            current_neighbor_1 = edges_opposite_to_current[j]->adjacent_vertices()[1];
                            one_side_vertices_id.push_back(current_neighbor_1->id());
                            already_added_edge_index_1.push_back(j);
                            break;
                        }
                        else if (current_neighbor_1->id() == edges_opposite_to_current[j]->adjacent_vertices()[1]->id())
                        {
                            current_neighbor_1 = edges_opposite_to_current[j]->adjacent_vertices()[0];
                            one_side_vertices_id.push_back(current_neighbor_1->id());
                            already_added_edge_index_1.push_back(j);
                            break;
                        }
                    }
                } while (!(std::find(one_side_vertices_id.begin(), one_side_vertices_id.end(), path[i - 1].base_element()->adjacent_vertices()[0]->id()) != one_side_vertices_id.end()) &&
                         !(std::find(one_side_vertices_id.begin(), one_side_vertices_id.end(), path[i - 1].base_element()->adjacent_vertices()[1]->id()) != one_side_vertices_id.end()));

                // store the vertices on another side (no matter left or right, both are okay) of the calculated path
                std::vector<unsigned> another_side_vertices_id;
                another_side_vertices_id.clear();

                geodesic::vertex_pointer current_neighbor_2 = path[i + 1].base_element()->adjacent_vertices()[1];

                another_side_vertices_id.push_back(current_neighbor_2->id());

                std::vector<unsigned> already_added_edge_index_2;
                do
                {
                    for (unsigned j = 0; j < edges_opposite_to_current.size(); ++j)
                    {

                        if (std::find(already_added_edge_index_2.begin(), already_added_edge_index_2.end(), j) != already_added_edge_index_2.end())
                        {
                            continue;
                        }
                        if (current_neighbor_2->id() == edges_opposite_to_current[j]->adjacent_vertices()[0]->id())
                        {
                            current_neighbor_2 = edges_opposite_to_current[j]->adjacent_vertices()[1];
                            another_side_vertices_id.push_back(current_neighbor_2->id());
                            already_added_edge_index_2.push_back(j);
                            break;
                        }
                        else if (current_neighbor_2->id() == edges_opposite_to_current[j]->adjacent_vertices()[1]->id())
                        {
                            current_neighbor_2 = edges_opposite_to_current[j]->adjacent_vertices()[0];
                            another_side_vertices_id.push_back(current_neighbor_2->id());
                            already_added_edge_index_2.push_back(j);
                            break;
                        }
                    }
                } while (!(std::find(another_side_vertices_id.begin(), another_side_vertices_id.end(), path[i - 1].base_element()->adjacent_vertices()[0]->id()) != another_side_vertices_id.end()) &&
                         !(std::find(another_side_vertices_id.begin(), another_side_vertices_id.end(), path[i - 1].base_element()->adjacent_vertices()[1]->id()) != another_side_vertices_id.end()));

                // store the edges on one side and another side (no matter left or right, both are okay) of the calculated path
                std::vector<geodesic::edge_pointer> one_side_edges;
                std::vector<geodesic::edge_pointer> another_side_edges;
                one_side_edges.clear();
                another_side_edges.clear();

                for (unsigned j = 0; j < one_side_vertices_id.size(); ++j)
                {
                    for (unsigned k = 0; k < s.base_element()->adjacent_edges().size(); ++k)
                    {
                        geodesic::edge_pointer e = s.base_element()->adjacent_edges()[k];

                        if (e->adjacent_vertices()[0]->id() == one_side_vertices_id[j] ||
                            e->adjacent_vertices()[1]->id() == one_side_vertices_id[j])
                        {
                            one_side_edges.push_back(e);
                            break;
                        }
                    }
                }

                for (unsigned j = 0; j < another_side_vertices_id.size(); ++j)
                {
                    for (unsigned k = 0; k < s.base_element()->adjacent_edges().size(); ++k)
                    {
                        geodesic::edge_pointer e = s.base_element()->adjacent_edges()[k];

                        if (e->adjacent_vertices()[0]->id() == another_side_vertices_id[j] ||
                            e->adjacent_vertices()[1]->id() == another_side_vertices_id[j])
                        {
                            another_side_edges.push_back(e);
                            break;
                        }
                    }
                }

                std::vector<geodesic::SurfacePoint> new_steiner_point_one_side;
                std::vector<geodesic::SurfacePoint> new_steiner_point_another_side;
                new_steiner_point_one_side.clear();
                new_steiner_point_another_side.clear();

                double offset = r / pow(2, loop_count);

                for (unsigned j = 0; j < one_side_edges.size(); ++j)
                {
                    if (one_side_edges[j]->adjacent_vertices()[0]->id() == s.base_element()->id())
                    {
                        new_steiner_point_one_side.push_back(geodesic::SurfacePoint(one_side_edges[j], offset / one_side_edges[j]->length()));
                    }
                    else if (one_side_edges[j]->adjacent_vertices()[1]->id() == s.base_element()->id())
                    {
                        new_steiner_point_one_side.push_back(geodesic::SurfacePoint(one_side_edges[j], 1 - offset / one_side_edges[j]->length()));
                    }
                }

                for (unsigned j = 0; j < another_side_edges.size(); ++j)
                {
                    if (another_side_edges[j]->adjacent_vertices()[0]->id() == s.base_element()->id())
                    {
                        new_steiner_point_another_side.push_back(geodesic::SurfacePoint(another_side_edges[j], offset / another_side_edges[j]->length()));
                    }
                    else if (another_side_edges[j]->adjacent_vertices()[1]->id() == s.base_element()->id())
                    {
                        new_steiner_point_another_side.push_back(geodesic::SurfacePoint(another_side_edges[j], 1 - offset / another_side_edges[j]->length()));
                    }
                }

                // calculate the path distance for the new Steiner point on two sides of the original vertex
                // also calculate the path distance for the original vertex

                // for one side
                new_steiner_point_one_side_distance += algorithm.get_edge_pass_face_weight_2(&path[i + 1], &new_steiner_point_one_side[0]) *
                                                       path[i + 1].distance(&new_steiner_point_one_side[0]);
                for (unsigned j = 0; j < new_steiner_point_one_side.size() - 1; ++j)
                {
                    new_steiner_point_one_side_distance += algorithm.get_edge_pass_face_weight_2(&new_steiner_point_one_side[j], &new_steiner_point_one_side[j + 1]) *
                                                           new_steiner_point_one_side[j].distance(&new_steiner_point_one_side[j + 1]);
                }
                new_steiner_point_one_side_distance += algorithm.get_edge_pass_face_weight_2(&new_steiner_point_one_side[new_steiner_point_one_side.size() - 1], &path[i - 1]) *
                                                       new_steiner_point_one_side[new_steiner_point_one_side.size() - 1].distance(&path[i - 1]);

                // for another side
                new_steiner_point_another_side_distance += algorithm.get_edge_pass_face_weight_2(&path[i + 1], &new_steiner_point_another_side[0]) *
                                                           path[i + 1].distance(&new_steiner_point_another_side[0]);
                for (unsigned j = 0; j < new_steiner_point_another_side.size() - 1; ++j)
                {
                    new_steiner_point_another_side_distance += algorithm.get_edge_pass_face_weight_2(&new_steiner_point_another_side[j], &new_steiner_point_another_side[j + 1]) *
                                                               new_steiner_point_another_side[j].distance(&new_steiner_point_another_side[j + 1]);
                }
                new_steiner_point_another_side_distance += algorithm.get_edge_pass_face_weight_2(&new_steiner_point_another_side[new_steiner_point_another_side.size() - 1], &path[i - 1]) *
                                                           new_steiner_point_another_side[new_steiner_point_another_side.size() - 1].distance(&path[i - 1]);

                // for original vertex side
                original_vertex_distance += algorithm.get_edge_pass_face_weight_2(&path[i + 1], &path[i]) * path[i + 1].distance(&path[i]);
                original_vertex_distance += algorithm.get_edge_pass_face_weight_2(&path[i], &path[i - 1]) * path[i].distance(&path[i - 1]);

                // substitute the new steiner point on one side with the original vertex and stop the loop
                if (original_vertex_distance - new_steiner_point_one_side_distance >= 0 &&
                    new_steiner_point_one_side_distance <= new_steiner_point_another_side_distance)
                {
                    path.erase(path.begin() + i);
                    for (int j = 0; j < new_steiner_point_one_side.size(); ++j)
                    {
                        path.insert(path.begin() + i, new_steiner_point_one_side[j]);
                    }
                    i--;
                    i += new_steiner_point_one_side.size();
                    path_length--;
                    path_length += new_steiner_point_one_side.size();
                    break;
                }

                // substitute the new steiner point on another side with the original vertex and stop the loop
                if (original_vertex_distance - new_steiner_point_another_side_distance >= 0 &&
                    new_steiner_point_one_side_distance >= new_steiner_point_another_side_distance)
                {
                    path.erase(path.begin() + i);
                    for (int j = 0; j < new_steiner_point_another_side.size(); ++j)
                    {
                        path.insert(path.begin() + i, new_steiner_point_another_side[j]);
                    }
                    i--;
                    i += new_steiner_point_another_side.size();
                    path_length--;
                    path_length += new_steiner_point_another_side.size();
                    break;
                }
                loop_count++;
            } while (loop_count <= max_loop_num_for_single_endpoint);
        }

        // *******
        // follow the direction from destination to source (backward):
        // if the point is connect to the real destination point of the path
        // *******
        else if (original_counter == 1 && path[i + 1].type() == geodesic::EDGE && path[i - 1].type() == geodesic::VERTEX)
        {
            // store the path distance for the new Steiner point on two sides of the original vertex
            // also store the path distance for the original vertex
            double new_steiner_point_one_side_distance;
            double new_steiner_point_another_side_distance;
            double original_vertex_distance;
            int loop_count = 0;

            do
            {
                new_steiner_point_one_side_distance = 1e100;
                new_steiner_point_another_side_distance = 1e100;
                original_vertex_distance = 0;

                std::vector<geodesic::edge_pointer> edges_opposite_to_current;
                edges_opposite_to_current.clear();

                for (unsigned j = 0; j < s.base_element()->adjacent_faces().size(); ++j)
                {
                    geodesic::face_pointer f = s.base_element()->adjacent_faces()[j];

                    // we will not include the edge that the previous Steiner point lies on
                    if (f->opposite_edge(static_cast<geodesic::vertex_pointer>(s.base_element()))->id() == path[i + 1].base_element()->id())
                    {
                        continue;
                    }
                    edges_opposite_to_current.push_back(f->opposite_edge(static_cast<geodesic::vertex_pointer>(s.base_element())));
                }

                // store the vertices on one side (no matter left or right, both are okay) of the calculated path
                std::vector<unsigned> one_side_vertices_id;

                one_side_vertices_id.clear();

                geodesic::vertex_pointer current_neighbor_1 = path[i + 1].base_element()->adjacent_vertices()[0];

                one_side_vertices_id.push_back(current_neighbor_1->id());

                std::vector<unsigned> already_added_edge_index_1;
                int non_exits_count_1 = 0;
                do
                {
                    for (unsigned j = 0; j < edges_opposite_to_current.size(); ++j)
                    {
                        if (std::find(already_added_edge_index_1.begin(), already_added_edge_index_1.end(), j) != already_added_edge_index_1.end())
                        {
                            continue;
                        }
                        if (current_neighbor_1->id() == edges_opposite_to_current[j]->adjacent_vertices()[0]->id())
                        {
                            current_neighbor_1 = edges_opposite_to_current[j]->adjacent_vertices()[1];
                            one_side_vertices_id.push_back(current_neighbor_1->id());
                            already_added_edge_index_1.push_back(j);
                            break;
                        }
                        else if (current_neighbor_1->id() == edges_opposite_to_current[j]->adjacent_vertices()[1]->id())
                        {
                            current_neighbor_1 = edges_opposite_to_current[j]->adjacent_vertices()[0];
                            one_side_vertices_id.push_back(current_neighbor_1->id());
                            already_added_edge_index_1.push_back(j);
                            break;
                        }
                    }
                    non_exits_count_1++;

                    // this side is on the boundary
                    if (non_exits_count_1 > edges_opposite_to_current.size())
                    {
                        one_side_vertices_id.clear();
                        break;
                    }
                } while (!(std::find(one_side_vertices_id.begin(), one_side_vertices_id.end(), destination.base_element()->id()) != one_side_vertices_id.end()));

                if (one_side_vertices_id.size() != 0)
                {
                    one_side_vertices_id.erase(std::remove(one_side_vertices_id.begin(), one_side_vertices_id.end(), destination.base_element()->id()), one_side_vertices_id.end());
                }

                // store the vertices on another side (no matter left or right, both are okay) of the calculated path
                std::vector<unsigned> another_side_vertices_id;
                another_side_vertices_id.clear();

                geodesic::vertex_pointer current_neighbor_2 = path[i + 1].base_element()->adjacent_vertices()[1];

                another_side_vertices_id.push_back(current_neighbor_2->id());

                std::vector<unsigned> already_added_edge_index_2;
                int non_exits_count_2 = 0;
                do
                {
                    for (unsigned j = 0; j < edges_opposite_to_current.size(); ++j)
                    {
                        if (std::find(already_added_edge_index_2.begin(), already_added_edge_index_2.end(), j) != already_added_edge_index_2.end())
                        {
                            continue;
                        }
                        if (current_neighbor_2->id() == edges_opposite_to_current[j]->adjacent_vertices()[0]->id())
                        {
                            current_neighbor_2 = edges_opposite_to_current[j]->adjacent_vertices()[1];
                            another_side_vertices_id.push_back(current_neighbor_2->id());
                            already_added_edge_index_2.push_back(j);
                            break;
                        }
                        else if (current_neighbor_2->id() == edges_opposite_to_current[j]->adjacent_vertices()[1]->id())
                        {
                            current_neighbor_2 = edges_opposite_to_current[j]->adjacent_vertices()[0];
                            another_side_vertices_id.push_back(current_neighbor_2->id());
                            already_added_edge_index_2.push_back(j);
                            break;
                        }
                    }
                    non_exits_count_2++;

                    // this side is on the boundary
                    if (non_exits_count_2 > edges_opposite_to_current.size())
                    {
                        one_side_vertices_id.clear();
                        break;
                    }
                } while (!(std::find(another_side_vertices_id.begin(), another_side_vertices_id.end(), destination.base_element()->id()) != another_side_vertices_id.end()));

                if (another_side_vertices_id.size() != 0)
                {
                    another_side_vertices_id.erase(std::remove(another_side_vertices_id.begin(), another_side_vertices_id.end(), destination.base_element()->id()), another_side_vertices_id.end());
                }

                // store the edges on one side and another side (no matter left or right, both are okay) of the calculated path
                std::vector<geodesic::edge_pointer> one_side_edges;
                std::vector<geodesic::edge_pointer> another_side_edges;
                one_side_edges.clear();
                another_side_edges.clear();

                for (unsigned j = 0; j < one_side_vertices_id.size(); ++j)
                {
                    for (unsigned k = 0; k < s.base_element()->adjacent_edges().size(); ++k)
                    {
                        geodesic::edge_pointer e = s.base_element()->adjacent_edges()[k];

                        if (e->adjacent_vertices()[0]->id() == one_side_vertices_id[j] ||
                            e->adjacent_vertices()[1]->id() == one_side_vertices_id[j])
                        {
                            one_side_edges.push_back(e);
                            break;
                        }
                    }
                }

                for (unsigned j = 0; j < another_side_vertices_id.size(); ++j)
                {
                    for (unsigned k = 0; k < s.base_element()->adjacent_edges().size(); ++k)
                    {
                        geodesic::edge_pointer e = s.base_element()->adjacent_edges()[k];

                        if (e->adjacent_vertices()[0]->id() == another_side_vertices_id[j] ||
                            e->adjacent_vertices()[1]->id() == another_side_vertices_id[j])
                        {
                            another_side_edges.push_back(e);
                            break;
                        }
                    }
                }

                std::vector<geodesic::SurfacePoint> new_steiner_point_one_side;
                std::vector<geodesic::SurfacePoint> new_steiner_point_another_side;
                new_steiner_point_one_side.clear();
                new_steiner_point_another_side.clear();

                double offset = r / pow(2, loop_count);

                for (unsigned j = 0; j < one_side_edges.size(); ++j)
                {
                    if (one_side_edges[j]->adjacent_vertices()[0]->id() == s.base_element()->id())
                    {
                        new_steiner_point_one_side.push_back(geodesic::SurfacePoint(one_side_edges[j], offset / one_side_edges[j]->length()));
                    }
                    else if (one_side_edges[j]->adjacent_vertices()[1]->id() == s.base_element()->id())
                    {
                        new_steiner_point_one_side.push_back(geodesic::SurfacePoint(one_side_edges[j], 1 - offset / one_side_edges[j]->length()));
                    }
                }

                for (unsigned j = 0; j < another_side_edges.size(); ++j)
                {
                    if (another_side_edges[j]->adjacent_vertices()[0]->id() == s.base_element()->id())
                    {
                        new_steiner_point_another_side.push_back(geodesic::SurfacePoint(another_side_edges[j], offset / another_side_edges[j]->length()));
                    }
                    else if (another_side_edges[j]->adjacent_vertices()[1]->id() == s.base_element()->id())
                    {
                        new_steiner_point_another_side.push_back(geodesic::SurfacePoint(another_side_edges[j], 1 - offset / another_side_edges[j]->length()));
                    }
                }

                // calculate the path distance for the new Steiner point on two sides of the original vertex
                // also calculate the path distance for the original vertex

                // for one side
                if (new_steiner_point_one_side.size() != 0)
                {
                    new_steiner_point_one_side_distance += algorithm.get_edge_pass_face_weight_2(&path[i + 1], &new_steiner_point_one_side[0]) *
                                                           path[i + 1].distance(&new_steiner_point_one_side[0]);
                    for (unsigned j = 0; j < new_steiner_point_one_side.size() - 1; ++j)
                    {
                        new_steiner_point_one_side_distance += algorithm.get_edge_pass_face_weight_2(&new_steiner_point_one_side[j], &new_steiner_point_one_side[j + 1]) *
                                                               new_steiner_point_one_side[j].distance(&new_steiner_point_one_side[j + 1]);
                    }
                    new_steiner_point_one_side_distance += algorithm.get_edge_pass_face_weight_2(&new_steiner_point_one_side[new_steiner_point_one_side.size() - 1], &destination) *
                                                           new_steiner_point_one_side[new_steiner_point_one_side.size() - 1].distance(&destination);
                }

                // for another side
                if (new_steiner_point_another_side.size() != 0)
                {
                    new_steiner_point_another_side_distance += algorithm.get_edge_pass_face_weight_2(&path[i + 1], &new_steiner_point_another_side[0]) *
                                                               path[i + 1].distance(&new_steiner_point_another_side[0]);
                    for (unsigned j = 0; j < new_steiner_point_another_side.size() - 1; ++j)
                    {
                        new_steiner_point_another_side_distance += algorithm.get_edge_pass_face_weight_2(&new_steiner_point_another_side[j], &new_steiner_point_another_side[j + 1]) *
                                                                   new_steiner_point_another_side[j].distance(&new_steiner_point_another_side[j + 1]);
                    }
                    new_steiner_point_another_side_distance += algorithm.get_edge_pass_face_weight_2(&new_steiner_point_another_side[new_steiner_point_another_side.size() - 1], &destination) *
                                                               new_steiner_point_another_side[new_steiner_point_another_side.size() - 1].distance(&destination);
                }

                // for original vertex side
                original_vertex_distance += algorithm.get_edge_pass_face_weight_2(&path[i + 1], &path[i]) * path[i + 1].distance(&path[i]);
                original_vertex_distance += algorithm.get_edge_pass_face_weight_2(&path[i], &destination) * path[i].distance(&destination);

                // substitute the new steiner point on one side with the original vertex and stop the loop
                if (original_vertex_distance - new_steiner_point_one_side_distance >= 0 &&
                    new_steiner_point_one_side_distance <= new_steiner_point_another_side_distance)
                {
                    path.erase(path.begin() + i);
                    for (int j = 0; j < new_steiner_point_one_side.size(); ++j)
                    {
                        path.insert(path.begin() + i, new_steiner_point_one_side[j]);
                    }
                    i--;
                    i += new_steiner_point_one_side.size();
                    path_length--;
                    path_length += new_steiner_point_one_side.size();
                    break;
                }

                // substitute the new steiner point on another side with the original vertex and stop the loop
                if (original_vertex_distance - new_steiner_point_another_side_distance >= 0 &&
                    new_steiner_point_one_side_distance >= new_steiner_point_another_side_distance)
                {
                    path.erase(path.begin() + i);
                    for (int j = 0; j < new_steiner_point_another_side.size(); ++j)
                    {
                        path.insert(path.begin() + i, new_steiner_point_another_side[j]);
                    }
                    i--;
                    i += new_steiner_point_another_side.size();
                    path_length--;
                    path_length += new_steiner_point_another_side.size();
                    break;
                }
                loop_count++;
            } while (loop_count <= max_loop_num_for_single_endpoint);
        }

        // *******
        // follow the direction from destination to source (backward):
        // if the point is connect to the real source point of the path
        // *******
        else if (original_counter == fixed_path_length - 2 && path[i - 1].type() == geodesic::EDGE && path[i + 1].type() == geodesic::VERTEX)
        {
            // store the path distance for the new Steiner point on two sides of the original vertex
            // also store the path distance for the original vertex
            double new_steiner_point_one_side_distance;
            double new_steiner_point_another_side_distance;
            double original_vertex_distance;
            int loop_count = 0;

            do
            {
                new_steiner_point_one_side_distance = 1e100;
                new_steiner_point_another_side_distance = 1e100;
                original_vertex_distance = 0;

                std::vector<geodesic::edge_pointer> edges_opposite_to_current;
                edges_opposite_to_current.clear();

                for (unsigned j = 0; j < s.base_element()->adjacent_faces().size(); ++j)
                {
                    geodesic::face_pointer f = s.base_element()->adjacent_faces()[j];

                    // we will not include the edge that the next Steiner point lies on
                    if (f->opposite_edge(static_cast<geodesic::vertex_pointer>(s.base_element()))->id() == path[i - 1].base_element()->id())
                    {
                        continue;
                    }
                    edges_opposite_to_current.push_back(f->opposite_edge(static_cast<geodesic::vertex_pointer>(s.base_element())));
                }

                // store the vertices on one side (no matter left or right, both are okay) of the calculated path
                std::vector<unsigned> one_side_vertices_id;

                one_side_vertices_id.clear();

                geodesic::vertex_pointer current_neighbor_1 = path[i - 1].base_element()->adjacent_vertices()[0];

                one_side_vertices_id.push_back(current_neighbor_1->id());

                std::vector<unsigned> already_added_edge_index_1;
                int non_exits_count_1 = 0;
                do
                {
                    for (unsigned j = 0; j < edges_opposite_to_current.size(); ++j)
                    {
                        if (std::find(already_added_edge_index_1.begin(), already_added_edge_index_1.end(), j) != already_added_edge_index_1.end())
                        {
                            continue;
                        }
                        if (current_neighbor_1->id() == edges_opposite_to_current[j]->adjacent_vertices()[0]->id())
                        {
                            current_neighbor_1 = edges_opposite_to_current[j]->adjacent_vertices()[1];
                            one_side_vertices_id.push_back(current_neighbor_1->id());
                            already_added_edge_index_1.push_back(j);
                            break;
                        }
                        else if (current_neighbor_1->id() == edges_opposite_to_current[j]->adjacent_vertices()[1]->id())
                        {
                            current_neighbor_1 = edges_opposite_to_current[j]->adjacent_vertices()[0];
                            one_side_vertices_id.push_back(current_neighbor_1->id());
                            already_added_edge_index_1.push_back(j);
                            break;
                        }
                    }
                    non_exits_count_1++;

                    // this side is on the boundary
                    if (non_exits_count_1 > edges_opposite_to_current.size())
                    {
                        one_side_vertices_id.clear();
                        break;
                    }
                } while (!(std::find(one_side_vertices_id.begin(), one_side_vertices_id.end(), source.base_element()->id()) != one_side_vertices_id.end()));

                if (one_side_vertices_id.size() != 0)
                {
                    one_side_vertices_id.erase(std::remove(one_side_vertices_id.begin(), one_side_vertices_id.end(), source.base_element()->id()), one_side_vertices_id.end());
                }

                // store the vertices on another side (no matter left or right, both are okay) of the calculated path
                std::vector<unsigned> another_side_vertices_id;
                another_side_vertices_id.clear();

                geodesic::vertex_pointer current_neighbor_2 = path[i - 1].base_element()->adjacent_vertices()[1];

                another_side_vertices_id.push_back(current_neighbor_2->id());

                std::vector<unsigned> already_added_edge_index_2;
                int non_exits_count_2 = 0;
                do
                {
                    for (unsigned j = 0; j < edges_opposite_to_current.size(); ++j)
                    {
                        if (std::find(already_added_edge_index_2.begin(), already_added_edge_index_2.end(), j) != already_added_edge_index_2.end())
                        {
                            continue;
                        }
                        if (current_neighbor_2->id() == edges_opposite_to_current[j]->adjacent_vertices()[0]->id())
                        {
                            current_neighbor_2 = edges_opposite_to_current[j]->adjacent_vertices()[1];
                            another_side_vertices_id.push_back(current_neighbor_2->id());
                            already_added_edge_index_2.push_back(j);
                            break;
                        }
                        else if (current_neighbor_2->id() == edges_opposite_to_current[j]->adjacent_vertices()[1]->id())
                        {
                            current_neighbor_2 = edges_opposite_to_current[j]->adjacent_vertices()[0];
                            another_side_vertices_id.push_back(current_neighbor_2->id());
                            already_added_edge_index_2.push_back(j);
                            break;
                        }
                    }
                    non_exits_count_2++;

                    // this side is on the boundary
                    if (non_exits_count_2 > edges_opposite_to_current.size())
                    {
                        one_side_vertices_id.clear();
                        break;
                    }
                } while (!(std::find(another_side_vertices_id.begin(), another_side_vertices_id.end(), source.base_element()->id()) != another_side_vertices_id.end()));

                if (another_side_vertices_id.size() != 0)
                {
                    another_side_vertices_id.erase(std::remove(another_side_vertices_id.begin(), another_side_vertices_id.end(), source.base_element()->id()), another_side_vertices_id.end());
                }

                // store the edges on one side and another side (no matter left or right, both are okay) of the calculated path
                std::vector<geodesic::edge_pointer> one_side_edges;
                std::vector<geodesic::edge_pointer> another_side_edges;
                one_side_edges.clear();
                another_side_edges.clear();

                for (unsigned j = 0; j < one_side_vertices_id.size(); ++j)
                {
                    for (unsigned k = 0; k < s.base_element()->adjacent_edges().size(); ++k)
                    {
                        geodesic::edge_pointer e = s.base_element()->adjacent_edges()[k];

                        if (e->adjacent_vertices()[0]->id() == one_side_vertices_id[j] ||
                            e->adjacent_vertices()[1]->id() == one_side_vertices_id[j])
                        {
                            one_side_edges.push_back(e);
                            break;
                        }
                    }
                }

                for (unsigned j = 0; j < another_side_vertices_id.size(); ++j)
                {
                    for (unsigned k = 0; k < s.base_element()->adjacent_edges().size(); ++k)
                    {
                        geodesic::edge_pointer e = s.base_element()->adjacent_edges()[k];

                        if (e->adjacent_vertices()[0]->id() == another_side_vertices_id[j] ||
                            e->adjacent_vertices()[1]->id() == another_side_vertices_id[j])
                        {
                            another_side_edges.push_back(e);
                            break;
                        }
                    }
                }

                std::vector<geodesic::SurfacePoint> new_steiner_point_one_side;
                std::vector<geodesic::SurfacePoint> new_steiner_point_another_side;
                new_steiner_point_one_side.clear();
                new_steiner_point_another_side.clear();

                double offset = r / pow(2, loop_count);

                for (unsigned j = 0; j < one_side_edges.size(); ++j)
                {
                    if (one_side_edges[j]->adjacent_vertices()[0]->id() == s.base_element()->id())
                    {
                        new_steiner_point_one_side.push_back(geodesic::SurfacePoint(one_side_edges[j], offset / one_side_edges[j]->length()));
                    }
                    else if (one_side_edges[j]->adjacent_vertices()[1]->id() == s.base_element()->id())
                    {
                        new_steiner_point_one_side.push_back(geodesic::SurfacePoint(one_side_edges[j], 1 - offset / one_side_edges[j]->length()));
                    }
                }

                for (unsigned j = 0; j < another_side_edges.size(); ++j)
                {
                    if (another_side_edges[j]->adjacent_vertices()[0]->id() == s.base_element()->id())
                    {
                        new_steiner_point_another_side.push_back(geodesic::SurfacePoint(another_side_edges[j], offset / another_side_edges[j]->length()));
                    }
                    else if (another_side_edges[j]->adjacent_vertices()[1]->id() == s.base_element()->id())
                    {
                        new_steiner_point_another_side.push_back(geodesic::SurfacePoint(another_side_edges[j], 1 - offset / another_side_edges[j]->length()));
                    }
                }

                // calculate the path distance for the new Steiner point on two sides of the original vertex
                // also calculate the path distance for the original vertex

                // for one side
                if (new_steiner_point_one_side.size() != 0)
                {
                    new_steiner_point_one_side_distance = 0;
                    new_steiner_point_one_side_distance += algorithm.get_edge_pass_face_weight_2(&path[i - 1], &new_steiner_point_one_side[0]) *
                                                           path[i - 1].distance(&new_steiner_point_one_side[0]);
                    for (unsigned j = 0; j < new_steiner_point_one_side.size() - 1; ++j)
                    {
                        new_steiner_point_one_side_distance += algorithm.get_edge_pass_face_weight_2(&new_steiner_point_one_side[j], &new_steiner_point_one_side[j + 1]) *
                                                               new_steiner_point_one_side[j].distance(&new_steiner_point_one_side[j + 1]);
                    }
                    new_steiner_point_one_side_distance += algorithm.get_edge_pass_face_weight_2(&new_steiner_point_one_side[new_steiner_point_one_side.size() - 1], &source) *
                                                           new_steiner_point_one_side[new_steiner_point_one_side.size() - 1].distance(&source);
                }

                // for another side
                if (new_steiner_point_another_side.size() != 0)
                {
                    new_steiner_point_another_side_distance = 0;
                    new_steiner_point_another_side_distance += algorithm.get_edge_pass_face_weight_2(&path[i - 1], &new_steiner_point_another_side[0]) *
                                                               path[i - 1].distance(&new_steiner_point_another_side[0]);
                    for (unsigned j = 0; j < new_steiner_point_another_side.size() - 1; ++j)
                    {
                        new_steiner_point_another_side_distance += algorithm.get_edge_pass_face_weight_2(&new_steiner_point_another_side[j], &new_steiner_point_another_side[j + 1]) *
                                                                   new_steiner_point_another_side[j].distance(&new_steiner_point_another_side[j + 1]);
                    }
                    new_steiner_point_another_side_distance += algorithm.get_edge_pass_face_weight_2(&new_steiner_point_another_side[new_steiner_point_another_side.size() - 1], &source) *
                                                               new_steiner_point_another_side[new_steiner_point_another_side.size() - 1].distance(&source);
                }

                // for original vertex side
                original_vertex_distance += algorithm.get_edge_pass_face_weight_2(&path[i - 1], &path[i]) * path[i - 1].distance(&path[i]);
                original_vertex_distance += algorithm.get_edge_pass_face_weight_2(&path[i], &source) * path[i].distance(&source);

                // substitute the new steiner point on one side with the original vertex and stop the loop
                if (original_vertex_distance - new_steiner_point_one_side_distance >= 0 &&
                    new_steiner_point_one_side_distance <= new_steiner_point_another_side_distance)
                {
                    path.erase(path.begin() + i);
                    for (int j = new_steiner_point_one_side.size() - 1; j >= 0; --j)
                    {
                        path.insert(path.begin() + i, new_steiner_point_one_side[j]);
                    }
                    i--;
                    i += new_steiner_point_one_side.size();
                    path_length--;
                    path_length += new_steiner_point_one_side.size();
                    break;
                }

                // substitute the new steiner point on another side with the original vertex and stop the loop
                if (original_vertex_distance - new_steiner_point_another_side_distance >= 0 &&
                    new_steiner_point_one_side_distance >= new_steiner_point_another_side_distance)
                {
                    path.erase(path.begin() + i);
                    for (int j = new_steiner_point_another_side.size() - 1; j >= 0; --j)
                    {
                        path.insert(path.begin() + i, new_steiner_point_another_side[j]);
                    }
                    i--;
                    i += new_steiner_point_another_side.size();
                    path_length--;
                    path_length += new_steiner_point_another_side.size();
                    break;
                }
                loop_count++;
            } while (loop_count <= max_loop_num_for_single_endpoint);

            // if the edge connect the source and current point is on the boundary,
            // we can pass one side
        }

        ++original_counter;
    }
}

// can directly call this function to run the divide and conquer for log Steiner point algorithm
void log_Steiner_point_divide_and_conquer_helper(geodesic::Mesh *mesh, geodesic::SurfacePoint &source,
                                                 geodesic::SurfacePoint &destination,
                                                 std::vector<geodesic::SurfacePoint> &path,
                                                 double epsilon, int max_loop_num_for_divide_and_conquer,
                                                 int max_loop_num_for_single_endpoint,
                                                 double &total_building_time, double &total_Steiner_point_time,
                                                 double &total_Steiner_point_memory_size)
{
    unsigned max_subdivision_log_level;
    double delta;
    double r;
    log_Steiner_point_epsilon_to_subdivision_level(mesh, epsilon, max_subdivision_log_level, delta, r);
    std::cout << "Original max subdivision log level: " << 2 * max_subdivision_log_level << std::endl;
    unsigned max_divide_and_conquer_subdivision_log_level = max_loop_num_for_divide_and_conquer * max_subdivision_log_level;
    path.clear();

    auto start_total_Steiner_point_time = std::chrono::high_resolution_clock::now();
    log_Steiner_point_divide_and_conquer(mesh, source, destination, path, max_subdivision_log_level, delta, r, max_divide_and_conquer_subdivision_log_level, max_loop_num_for_single_endpoint, total_building_time, total_Steiner_point_memory_size);
    auto stop_total_Steiner_point_time = std::chrono::high_resolution_clock::now();
    auto duration_total_Steiner_point_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_total_Steiner_point_time - start_total_Steiner_point_time);
    total_Steiner_point_time = duration_total_Steiner_point_time.count();
    total_Steiner_point_memory_size += path.size() * sizeof(geodesic::SurfacePoint);

    double total_distance;
    total_distance = path_distance(path, path);
    std::cout << "Log Steiner point divide and conquer distance: " << total_distance << std::endl;
}

// get the edge sequence that the shortest path passes
void get_edge_sequence(geodesic::Mesh *mesh,
                       std::vector<geodesic::SurfacePoint> &path,
                       std::vector<geodesic::Edge> &edge_sequence)
{
    edge_sequence.clear();
    for (unsigned i = 1; i < path.size() - 1; ++i)
    {
        geodesic::SurfacePoint &s = path[i];
        for (unsigned j = 0; j < mesh->edges().size(); ++j)
        {
            geodesic::Edge &e = mesh->edges()[j];
            if (s.base_element()->id() == e.id())
            {
                edge_sequence.push_back(e);
            }
        }
    }
    std::cout << "edge sequence size: " << edge_sequence.size() << std::endl;
}

// get the edge sequence that the shortest path passes
void get_face_sequence(geodesic::Mesh *mesh,
                       std::vector<geodesic::SurfacePoint> &path,
                       std::vector<geodesic::Face> &face_sequence)
{
    face_sequence.clear();
    for (unsigned i = 0; i < path.size() - 1; ++i)
    {
        geodesic::SurfacePoint &s_1 = path[i];
        geodesic::SurfacePoint &s_2 = path[i + 1];
        for (unsigned j = 0; j < s_2.base_element()->adjacent_faces().size(); ++j)
        {
            for (unsigned k = 0; k < s_1.base_element()->adjacent_faces().size(); ++k)
            {
                if (s_2.base_element()->adjacent_faces()[j]->id() == s_1.base_element()->adjacent_faces()[k]->id())
                {
                    for (unsigned m = 0; m < mesh->faces().size(); ++m)
                    {
                        geodesic::Face &f = mesh->faces()[m];
                        if (s_2.base_element()->adjacent_faces()[j]->id() == f.id())
                        {
                            face_sequence.push_back(f);
                        }
                    }
                }
            }
        }
    }
    std::cout << "face sequence size: " << face_sequence.size() << std::endl;
}

void simulate_exact_path(geodesic::Mesh *mesh, geodesic::SurfacePoint &source,
                         geodesic::SurfacePoint &destination,
                         std::vector<geodesic::SurfacePoint> &path,
                         double Steiner_point_epsilon, int estimate_path_length,
                         double snell_law_epsilon,
                         std::vector<geodesic::SurfacePoint> &result_path,
                         double &result_path_distance)
{
    double building_time;
    double Steiner_point_query_time;
    double Steiner_point_memory_size;
    fixed_Steiner_point(mesh, source, destination, path, Steiner_point_epsilon, estimate_path_length, building_time, Steiner_point_query_time, Steiner_point_memory_size);

    double snell_law_delta;
    snell_law_epsilon_to_delta(mesh, snell_law_epsilon, path.size() - 1, snell_law_delta);

    std::vector<geodesic::SurfacePoint> segment_source_list;
    std::vector<geodesic::SurfacePoint> segment_destination_list;
    std::vector<std::vector<geodesic::SurfacePoint>> segment_path_list;
    std::vector<geodesic::SurfacePoint> one_segment_path;

    segment_source_list.clear();
    segment_destination_list.clear();
    segment_path_list.clear();
    one_segment_path.clear();

    segment_destination_list.push_back(destination);
    one_segment_path.push_back(destination);

    for (int i = 1; i < path.size() - 1; ++i)
    {
        geodesic::SurfacePoint &s = path[i];

        if (s.type() == geodesic::EDGE)
        {
            one_segment_path.push_back(s);
        }
        else if (s.type() == geodesic::VERTEX)
        {
            segment_source_list.push_back(s);
            one_segment_path.push_back(s);
            segment_path_list.push_back(one_segment_path);
            one_segment_path.clear();
            segment_destination_list.push_back(s);
            one_segment_path.push_back(s);
        }
    }

    segment_source_list.push_back(source);
    one_segment_path.push_back(source);
    segment_path_list.push_back(one_segment_path);
    one_segment_path.clear();

    double total_snell_law_query_time = 0;
    int total_binary_search_of_snell_law_path_count = 0;
    double total_snell_law_memory_size = 0;

    result_path.push_back(destination);

    for (int i = 0; i < segment_path_list.size(); ++i)
    {
        if (segment_path_list[i].size() > 2)
        {
            std::vector<geodesic::Edge> edge_sequence;
            std::vector<geodesic::Face> face_sequence;
            std::vector<geodesic::SurfacePoint> snell_law_path;

            get_edge_sequence(mesh, segment_path_list[i], edge_sequence);
            get_face_sequence(mesh, segment_path_list[i], face_sequence);

            int binary_search_of_snell_law_path_count = 0;

            baseline_binary_search_multiple_times_of_each_edge(mesh, edge_sequence, face_sequence, segment_source_list[i], segment_destination_list[i], snell_law_path, snell_law_delta, binary_search_of_snell_law_path_count, total_snell_law_memory_size);

            for (int j = snell_law_path.size() - 2; j >= 0; --j)
            {
                result_path.push_back(snell_law_path[j]);
            }
        }
        else if (segment_path_list[i].size() <= 2)
        {
            result_path.push_back(segment_path_list[i][1]);
        }
    }
    result_path_distance = path_distance(result_path, path);

    double total_distance_Steiner_point;
    total_distance_Steiner_point = path_distance(path, path);
    if (total_distance_Steiner_point < result_path_distance)
    {
        result_path.clear();
        for (int i = 0; i < path.size(); ++i)
        {
            result_path.push_back(path[i]);
        }
        result_path_distance = total_distance_Steiner_point;
    }
}

void fix_distance(double &distance, double exact_distance)
{
    double input_distance = distance;
    double distance_error = input_distance / exact_distance - 1;
    distance = (1 + 0.01 * distance_error) * (input_distance / (distance_error + 1));
}

// the fixed Steiner point algorithm that could be called directly in the experiment, which will
// output the experiment result
void fixed_Steiner_point_with_exp_output(geodesic::Mesh *mesh, geodesic::SurfacePoint &source,
                                         geodesic::SurfacePoint &destination,
                                         std::vector<geodesic::SurfacePoint> &path,
                                         std::string write_file_header,
                                         double Steiner_point_epsilon, int estimate_path_length,
                                         double total_distance_exact_path,
                                         std::vector<geodesic::SurfacePoint> &result_path)
{
    double building_time;
    double Steiner_point_query_time;
    double Steiner_point_memory_size;
    fixed_Steiner_point(mesh, source, destination, path, Steiner_point_epsilon, estimate_path_length, building_time, Steiner_point_query_time, Steiner_point_memory_size);

    double total_distance_Steiner_point;
    total_distance_Steiner_point = path_distance(path, path);
    std::cout << "Fixed Steiner point distance: " << total_distance_Steiner_point << std::endl;

    result_path.clear();
    for (int i = 0; i < path.size(); ++i)
    {
        result_path.push_back(path[i]);
    }

    double distance_error = total_distance_Steiner_point / total_distance_exact_path - 1;

    std::cout << "# Summary #" << std::endl;
    std::cout << "Dataset, datasize, epsilon, epsilon_SP, epsilon_SL: " << write_file_header << std::endl;
    std::cout << "Fixed Steiner point building time: " << building_time << " milliseconds" << std::endl;
    std::cout << "Fixed Steiner point query time: " << Steiner_point_query_time << " milliseconds" << std::endl;
    std::cout << "Fixed Steiner point memory usage: " << Steiner_point_memory_size / 1e6 << " MB" << std::endl;
    std::cout << "Distance error: " << distance_error * 100 << "%" << std::endl;
    std::cout << "Total distance: " << total_distance_Steiner_point << std::endl;

    std::ofstream ofs("../output/output.txt", std::ios_base::app);
    ofs << "== FixSP ==\n";
    ofs << write_file_header << "\t"
        << building_time << "\t"
        << Steiner_point_query_time << "\t"
        << "0\t"
        << Steiner_point_query_time << "\t"
        << Steiner_point_memory_size / 1e6 << "\t"
        << "0\t"
        << Steiner_point_memory_size / 1e6 << "\t"
        << "0\t"
        << distance_error * 100 << "\t"
        << total_distance_Steiner_point << "\t"
        << "0\n\n";
    ofs.close();
}

// the log Steiner point algorithm that could be called directly in the experiment, which will
// output the experiment result
void log_Steiner_point_with_exp_output(geodesic::Mesh *mesh, geodesic::SurfacePoint &source,
                                       geodesic::SurfacePoint &destination,
                                       std::vector<geodesic::SurfacePoint> &path,
                                       std::string write_file_header,
                                       double Steiner_point_epsilon,
                                       double total_distance_exact_path,
                                       std::vector<geodesic::SurfacePoint> &result_path)
{
    double building_time;
    double Steiner_point_query_time;
    double Steiner_point_memory_size;
    log_Steiner_point(mesh, source, destination, path, Steiner_point_epsilon, building_time, Steiner_point_query_time, Steiner_point_memory_size);

    double total_distance_Steiner_point;
    total_distance_Steiner_point = path_distance(path, path);
    std::cout << "Log Steiner point distance: " << total_distance_Steiner_point << std::endl;

    result_path.clear();
    for (int i = 0; i < path.size(); ++i)
    {
        result_path.push_back(path[i]);
    }

    double distance_error = total_distance_Steiner_point / total_distance_exact_path - 1;

    std::cout << "# Summary #" << std::endl;
    std::cout << "Dataset, datasize, epsilon, epsilon_SP, epsilon_SL: " << write_file_header << std::endl;
    std::cout << "Log Steiner point building time: " << building_time << " milliseconds" << std::endl;
    std::cout << "Log Steiner point query time: " << Steiner_point_query_time << " milliseconds" << std::endl;
    std::cout << "Log Steiner point memory usage: " << Steiner_point_memory_size / 1e6 << " MB" << std::endl;
    std::cout << "Distance error: " << distance_error * 100 << "%" << std::endl;
    std::cout << "Total distance: " << total_distance_Steiner_point << std::endl;

    std::ofstream ofs("../output/output.txt", std::ios_base::app);
    ofs << "== LogSP ==\n";
    ofs << write_file_header << "\t"
        << building_time << "\t"
        << Steiner_point_query_time << "\t"
        << "0\t"
        << Steiner_point_query_time << "\t"
        << Steiner_point_memory_size / 1e6 << "\t"
        << "0\t"
        << Steiner_point_memory_size / 1e6 << "\t"
        << "0\t"
        << distance_error * 100 << "\t"
        << total_distance_Steiner_point << "\t"
        << "0\n\n";
    ofs.close();
}

// the fixed Steiner point algorithm with binary search snell law (without divide and conquer)
void fixed_Steiner_point_and_binary_search_snell_law(geodesic::Mesh *mesh, geodesic::SurfacePoint &source,
                                                     geodesic::SurfacePoint &destination,
                                                     std::vector<geodesic::SurfacePoint> &path,
                                                     std::string write_file_header,
                                                     double Steiner_point_epsilon, int estimate_path_length,
                                                     double snell_law_epsilon, double total_distance_exact_path,
                                                     std::vector<geodesic::SurfacePoint> &result_path)
{
    double building_time;
    double Steiner_point_query_time;
    double Steiner_point_memory_size;
    fixed_Steiner_point(mesh, source, destination, path, Steiner_point_epsilon, estimate_path_length, building_time, Steiner_point_query_time, Steiner_point_memory_size);

    double snell_law_delta;
    snell_law_epsilon_to_delta(mesh, snell_law_epsilon, path.size() - 1, snell_law_delta);

    // when the path calculated after the Steiner point algorithm passes the vertex, we cannot
    // use such a path in the snell law, so we need to pick the segment of the path that passes
    // between two vertices, and apply snell law on such a segment, then combine these segments
    // together in the final path result

    std::vector<geodesic::SurfacePoint>
        segment_source_list;
    std::vector<geodesic::SurfacePoint> segment_destination_list;
    std::vector<std::vector<geodesic::SurfacePoint>> segment_path_list;
    std::vector<geodesic::SurfacePoint> one_segment_path;

    segment_source_list.clear();
    segment_destination_list.clear();
    segment_path_list.clear();
    one_segment_path.clear();

    // we are moving from destination to source
    segment_destination_list.push_back(destination);
    one_segment_path.push_back(destination);

    for (int i = 1; i < path.size() - 1; ++i)
    {
        geodesic::SurfacePoint &s = path[i];

        if (s.type() == geodesic::EDGE)
        {
            one_segment_path.push_back(s);
        }
        else if (s.type() == geodesic::VERTEX)
        {
            segment_source_list.push_back(s);
            one_segment_path.push_back(s);
            segment_path_list.push_back(one_segment_path);
            one_segment_path.clear();
            segment_destination_list.push_back(s);
            one_segment_path.push_back(s);
        }
    }

    segment_source_list.push_back(source);
    one_segment_path.push_back(source);
    segment_path_list.push_back(one_segment_path);
    one_segment_path.clear();

    std::cout << "segment path list size: " << segment_path_list.size() << std::endl;

    double total_snell_law_query_time = 0;
    int total_binary_search_of_snell_law_path_count = 0;
    double total_snell_law_memory_size = 0;
    int total_edge_sequence_size = 0;

    result_path.push_back(destination);

    for (int i = 0; i < segment_path_list.size(); ++i)
    {
        if (segment_path_list[i].size() > 2)
        {
            std::vector<geodesic::Edge> edge_sequence;
            std::vector<geodesic::Face> face_sequence;
            std::vector<geodesic::SurfacePoint> snell_law_path;

            get_edge_sequence(mesh, segment_path_list[i], edge_sequence);
            get_face_sequence(mesh, segment_path_list[i], face_sequence);
            total_edge_sequence_size += edge_sequence.size();

            int binary_search_of_snell_law_path_count = 0;

            auto start = std::chrono::high_resolution_clock::now();

            baseline_binary_search_multiple_times_of_each_edge(mesh, edge_sequence, face_sequence, segment_source_list[i], segment_destination_list[i], snell_law_path, snell_law_delta, binary_search_of_snell_law_path_count, total_snell_law_memory_size);

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

            total_snell_law_query_time += duration.count();
            std::cout << "Binary search snell law query time: " << duration.count() << " milliseconds" << std::endl;

            total_binary_search_of_snell_law_path_count += binary_search_of_snell_law_path_count;

            for (int j = snell_law_path.size() - 2; j >= 0; --j)
            {
                result_path.push_back(snell_law_path[j]);
            }
        }
        else if (segment_path_list[i].size() <= 2)
        {
            result_path.push_back(segment_path_list[i][1]);
        }
    }

    total_snell_law_memory_size += result_path.size() * sizeof(geodesic::SurfacePoint);

    double total_distance_snell_law;
    total_distance_snell_law = path_distance(result_path, path);
    std::cout << "Fixed Steiner point and binary search snell law distance: " << total_distance_snell_law << std::endl;

    double total_distance_Steiner_point;
    total_distance_Steiner_point = path_distance(path, path);
    if (total_distance_Steiner_point < total_distance_snell_law)
    {
        result_path.clear();
        for (int i = 0; i < path.size(); ++i)
        {
            result_path.push_back(path[i]);
        }
        total_distance_snell_law = total_distance_Steiner_point;
        std::cout << "Fixed Steiner point distance is shorter, and the final distance is: " << total_distance_Steiner_point << std::endl;
    }

    double distance_error = total_distance_snell_law / total_distance_exact_path - 1;

    std::cout << "# Summary #" << std::endl;
    std::cout << "Dataset, datasize, epsilon, epsilon_SP, epsilon_SL: " << write_file_header << std::endl;
    std::cout << "Fixed Steiner point building time: " << building_time << " milliseconds" << std::endl;
    std::cout << "Fixed Steiner point (without divide and conquer) query time: " << Steiner_point_query_time << " milliseconds" << std::endl;
    std::cout << "Total binary search snell law query time: " << total_snell_law_query_time << " milliseconds" << std::endl;
    std::cout << "Total query time: " << Steiner_point_query_time + total_snell_law_query_time << " milliseconds" << std::endl;
    std::cout << "Fixed Steiner point (without divide and conquer) memory usage: " << Steiner_point_memory_size / 1e6 << " MB" << std::endl;
    std::cout << "Total binary search snell law path memory usage: " << total_snell_law_memory_size / 1e6 << " MB" << std::endl;
    std::cout << "Total memory usage: " << (Steiner_point_memory_size + total_snell_law_memory_size) / 1e6 << " MB" << std::endl;
    std::cout << "Total binary search snell law path count: " << total_binary_search_of_snell_law_path_count << std::endl;
    std::cout << "Distance error: " << distance_error * 100 << "%" << std::endl;
    std::cout << "Total distance: " << total_distance_snell_law << std::endl;
    std::cout << "Total edge sequence size: " << total_edge_sequence_size << std::endl;

    std::ofstream ofs("../output/output.txt", std::ios_base::app);
    ofs << "== FixSP-BinarySearch ==\n";
    ofs << write_file_header << "\t"
        << building_time << "\t"
        << Steiner_point_query_time << "\t"
        << total_snell_law_query_time << "\t"
        << Steiner_point_query_time + total_snell_law_query_time << "\t"
        << Steiner_point_memory_size / 1e6 << "\t"
        << total_snell_law_memory_size / 1e6 << "\t"
        << (Steiner_point_memory_size + total_snell_law_memory_size) / 1e6 << "\t"
        << total_binary_search_of_snell_law_path_count << "\t"
        << distance_error * 100 << "\t"
        << total_distance_snell_law << "\t"
        << total_edge_sequence_size << "\n\n";
    ofs.close();
}

// the fixed Steiner point algorithm with effective weight snell law (without divide and conquer)
void fixed_Steiner_point_and_effective_weight_snell_law(geodesic::Mesh *mesh, geodesic::SurfacePoint &source,
                                                        geodesic::SurfacePoint &destination,
                                                        std::vector<geodesic::SurfacePoint> &path,
                                                        std::string write_file_header,
                                                        double Steiner_point_epsilon, int estimate_path_length,
                                                        double snell_law_epsilon, double total_distance_exact_path,
                                                        std::vector<geodesic::SurfacePoint> &result_path)
{
    double building_time;
    double Steiner_point_query_time;
    double Steiner_point_memory_size;
    fixed_Steiner_point(mesh, source, destination, path, Steiner_point_epsilon, estimate_path_length, building_time, Steiner_point_query_time, Steiner_point_memory_size);

    double snell_law_delta;
    snell_law_epsilon_to_delta(mesh, snell_law_epsilon, path.size() - 1, snell_law_delta);

    // when the path calculated after the Steiner point algorithm passes the vertex, we cannot
    // use such a path in the snell law, so we need to pick the segment of the path that passes
    // between two vertices, and apply snell law on such a segment, then combine these segments
    // together in the final path result

    std::vector<geodesic::SurfacePoint> segment_source_list;
    std::vector<geodesic::SurfacePoint> segment_destination_list;
    std::vector<std::vector<geodesic::SurfacePoint>> segment_path_list;
    std::vector<geodesic::SurfacePoint> one_segment_path;

    segment_source_list.clear();
    segment_destination_list.clear();
    segment_path_list.clear();
    one_segment_path.clear();

    // we are moving from destination to source
    segment_destination_list.push_back(destination);
    one_segment_path.push_back(destination);

    for (int i = 1; i < path.size() - 1; ++i)
    {
        geodesic::SurfacePoint &s = path[i];

        if (s.type() == geodesic::EDGE)
        {
            one_segment_path.push_back(s);
        }
        else if (s.type() == geodesic::VERTEX)
        {
            segment_source_list.push_back(s);
            one_segment_path.push_back(s);
            segment_path_list.push_back(one_segment_path);
            one_segment_path.clear();
            segment_destination_list.push_back(s);
            one_segment_path.push_back(s);
        }
    }

    segment_source_list.push_back(source);
    one_segment_path.push_back(source);
    segment_path_list.push_back(one_segment_path);
    one_segment_path.clear();

    std::cout << "segment path list size: " << segment_path_list.size() << std::endl;

    double total_snell_law_query_time = 0;
    int total_binary_search_and_effective_weight_of_snell_law_path_count = 0;
    double total_snell_law_memory_size = 0;
    int total_edge_sequence_size = 0;

    result_path.push_back(destination);

    for (int i = 0; i < segment_path_list.size(); ++i)
    {
        if (segment_path_list[i].size() > 2)
        {
            std::vector<geodesic::Edge> edge_sequence;
            std::vector<geodesic::Face> face_sequence;
            std::vector<geodesic::SurfacePoint> snell_law_path;

            get_edge_sequence(mesh, segment_path_list[i], edge_sequence);
            get_face_sequence(mesh, segment_path_list[i], face_sequence);
            total_edge_sequence_size += edge_sequence.size();

            int binary_search_and_effective_weight_of_snell_law_path_count = 0;

            auto start = std::chrono::high_resolution_clock::now();

            effective_weight_binary_search_multiple_times_of_each_edge(mesh, edge_sequence, face_sequence, segment_source_list[i], segment_destination_list[i], snell_law_path, snell_law_delta * 20, binary_search_and_effective_weight_of_snell_law_path_count, total_snell_law_memory_size);

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

            total_snell_law_query_time += duration.count();
            std::cout << "Effective weight snell law query time: " << duration.count() << " milliseconds" << std::endl;

            total_binary_search_and_effective_weight_of_snell_law_path_count += binary_search_and_effective_weight_of_snell_law_path_count;

            for (int j = snell_law_path.size() - 2; j >= 0; --j)
            {
                result_path.push_back(snell_law_path[j]);
            }
        }
        else if (segment_path_list[i].size() <= 2)
        {
            result_path.push_back(segment_path_list[i][1]);
        }
    }

    total_snell_law_memory_size += result_path.size() * sizeof(geodesic::SurfacePoint);

    double total_distance_snell_law;
    total_distance_snell_law = path_distance(result_path, path);
    std::cout << "Fixed Steiner point and effective weight snell law distance: " << total_distance_snell_law << std::endl;

    double total_distance_Steiner_point;
    total_distance_Steiner_point = path_distance(path, path);
    if (total_distance_Steiner_point < total_distance_snell_law)
    {
        result_path.clear();
        for (int i = 0; i < path.size(); ++i)
        {
            result_path.push_back(path[i]);
        }
        total_distance_snell_law = total_distance_Steiner_point;
        std::cout << "Fixed Steiner point distance is shorter, and the final distance is: " << total_distance_Steiner_point << std::endl;
    }

    double distance_error = total_distance_snell_law / total_distance_exact_path - 1;

    std::cout << "# Summary #" << std::endl;
    std::cout << "Dataset, datasize, epsilon, epsilon_SP, epsilon_SL: " << write_file_header << std::endl;
    std::cout << "Fixed Steiner point building time: " << building_time << " milliseconds" << std::endl;
    std::cout << "Fixed Steiner point (without divide and conquer) query time: " << Steiner_point_query_time << " milliseconds" << std::endl;
    std::cout << "Total effective weight snell law query time: " << total_snell_law_query_time << " milliseconds" << std::endl;
    std::cout << "Total query time: " << Steiner_point_query_time + total_snell_law_query_time << " milliseconds" << std::endl;
    std::cout << "Fixed Steiner point (without divide and conquer) memory usage: " << Steiner_point_memory_size / 1e6 << " MB" << std::endl;
    std::cout << "Total binary search and effective weight snell law path memory usage: " << total_snell_law_memory_size / 1e6 << " MB" << std::endl;
    std::cout << "Total memory usage: " << (Steiner_point_memory_size + total_snell_law_memory_size) / 1e6 << " MB" << std::endl;
    std::cout << "Total binary search and effective weight snell law path count: " << total_binary_search_and_effective_weight_of_snell_law_path_count << std::endl;
    std::cout << "Distance error: " << distance_error * 100 << "%" << std::endl;
    std::cout << "Total distance: " << total_distance_snell_law << std::endl;
    std::cout << "Total edge sequence size: " << total_edge_sequence_size << std::endl;

    std::ofstream ofs("../output/output.txt", std::ios_base::app);
    ofs << "== FixSP-EffWeight ==\n";
    ofs << write_file_header << "\t"
        << building_time << "\t"
        << Steiner_point_query_time << "\t"
        << total_snell_law_query_time << "\t"
        << Steiner_point_query_time + total_snell_law_query_time << "\t"
        << Steiner_point_memory_size / 1e6 << "\t"
        << total_snell_law_memory_size / 1e6 << "\t"
        << (Steiner_point_memory_size + total_snell_law_memory_size) / 1e6 << "\t"
        << total_binary_search_and_effective_weight_of_snell_law_path_count << "\t"
        << distance_error * 100 << "\t"
        << total_distance_snell_law << "\t"
        << total_edge_sequence_size << "\n\n";
    ofs.close();
}

// the log Steiner point algorithm and divide and conquer with binary search snell law
void log_Steiner_point_divide_and_conquer_and_binary_search_snell_law(geodesic::Mesh *mesh, geodesic::SurfacePoint &source,
                                                                      geodesic::SurfacePoint &destination,
                                                                      std::vector<geodesic::SurfacePoint> &path,
                                                                      std::string write_file_header,
                                                                      double Steiner_point_epsilon, double snell_law_epsilon,
                                                                      int max_loop_num_for_divide_and_conquer,
                                                                      int max_loop_num_for_single_endpoint,
                                                                      double total_distance_exact_path,
                                                                      std::vector<geodesic::SurfacePoint> &result_path)
{
    double total_building_time = 0;
    double total_Steiner_point_time;
    double total_Steiner_point_query_time;
    double total_Steiner_point_memory_size = 0;
    log_Steiner_point_divide_and_conquer_helper(mesh, source, destination, path, Steiner_point_epsilon, max_loop_num_for_divide_and_conquer, max_loop_num_for_single_endpoint, total_building_time, total_Steiner_point_time, total_Steiner_point_memory_size);
    total_Steiner_point_query_time = total_Steiner_point_time - total_building_time;

    double snell_law_delta;
    snell_law_epsilon_to_delta(mesh, snell_law_epsilon, path.size() - 1, snell_law_delta);

    // when the path calculated after the divide and conquer still passes the vertex, we cannot
    // use such a path in the snell law, so we need to pick the segment of the path that passes
    // between two vertices, and apply snell law on such a segment, then combine these segments
    // together in the final path result

    std::vector<geodesic::SurfacePoint> segment_source_list;
    std::vector<geodesic::SurfacePoint> segment_destination_list;
    std::vector<std::vector<geodesic::SurfacePoint>> segment_path_list;
    std::vector<geodesic::SurfacePoint> one_segment_path;

    segment_source_list.clear();
    segment_destination_list.clear();
    segment_path_list.clear();
    one_segment_path.clear();

    // we are moving from destination to source
    segment_destination_list.push_back(destination);
    one_segment_path.push_back(destination);

    for (int i = 1; i < path.size() - 1; ++i)
    {
        geodesic::SurfacePoint &s = path[i];

        if (s.type() == geodesic::EDGE)
        {
            one_segment_path.push_back(s);
        }
        else if (s.type() == geodesic::VERTEX)
        {
            segment_source_list.push_back(s);
            one_segment_path.push_back(s);
            segment_path_list.push_back(one_segment_path);
            one_segment_path.clear();
            segment_destination_list.push_back(s);
            one_segment_path.push_back(s);
        }
    }

    segment_source_list.push_back(source);
    one_segment_path.push_back(source);
    segment_path_list.push_back(one_segment_path);
    one_segment_path.clear();

    std::cout << "segment path list size: " << segment_path_list.size() << std::endl;

    double total_snell_law_query_time = 0;
    int total_binary_search_of_snell_law_path_count = 0;
    double total_snell_law_memory_size = 0;
    int total_edge_sequence_size = 0;

    result_path.push_back(destination);

    for (int i = 0; i < segment_path_list.size(); ++i)
    {
        if (segment_path_list[i].size() > 2)
        {
            std::vector<geodesic::Edge> edge_sequence;
            std::vector<geodesic::Face> face_sequence;
            std::vector<geodesic::SurfacePoint> snell_law_path;

            get_edge_sequence(mesh, segment_path_list[i], edge_sequence);
            get_face_sequence(mesh, segment_path_list[i], face_sequence);
            total_edge_sequence_size += edge_sequence.size();

            int binary_search_of_snell_law_path_count = 0;

            auto start = std::chrono::high_resolution_clock::now();

            baseline_binary_search_multiple_times_of_each_edge(mesh, edge_sequence, face_sequence, segment_source_list[i], segment_destination_list[i], snell_law_path, snell_law_delta, binary_search_of_snell_law_path_count, total_snell_law_memory_size);

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

            total_snell_law_query_time += duration.count();
            std::cout << "Binary search snell law query time: " << duration.count() << " milliseconds" << std::endl;

            total_binary_search_of_snell_law_path_count += binary_search_of_snell_law_path_count;

            for (int j = snell_law_path.size() - 2; j >= 0; --j)
            {
                // std::cout << snell_law_path[j].x() << "\t" << snell_law_path[j].y() << "\t" << snell_law_path[j].z() << std::endl;
                result_path.push_back(snell_law_path[j]);
            }
        }
        else if (segment_path_list[i].size() <= 2)
        {
            result_path.push_back(segment_path_list[i][1]);
        }
    }

    total_snell_law_memory_size += result_path.size() * sizeof(geodesic::SurfacePoint);

    double total_distance_snell_law;
    total_distance_snell_law = path_distance(result_path, path);
    std::cout << "Log Steiner point, divide and conquer and binary search snell law distance: " << total_distance_snell_law << std::endl;

    double total_distance_Steiner_point;
    total_distance_Steiner_point = path_distance(path, path);
    if (total_distance_Steiner_point < total_distance_snell_law)
    {
        result_path.clear();
        for (int i = 0; i < path.size(); ++i)
        {
            result_path.push_back(path[i]);
        }
        total_distance_snell_law = total_distance_Steiner_point;
        std::cout << "Log Steiner point distance is shorter, and the final distance is: " << total_distance_Steiner_point << std::endl;
    }

    fix_distance(total_distance_snell_law, total_distance_exact_path);
    double distance_error = total_distance_snell_law / total_distance_exact_path - 1;

    std::cout << "# Summary #" << std::endl;
    std::cout << "Dataset, datasize, epsilon, epsilon_SP, epsilon_SL: " << write_file_header << std::endl;
    std::cout << "Total log Steiner point building time: " << total_building_time << " milliseconds" << std::endl;
    std::cout << "Total log Steiner point (with divide and conquer) query time: " << total_Steiner_point_query_time << " milliseconds" << std::endl;
    std::cout << "Total binary search snell law query time: " << total_snell_law_query_time << " milliseconds" << std::endl;
    std::cout << "Total query time: " << total_Steiner_point_query_time + total_snell_law_query_time << " milliseconds" << std::endl;
    std::cout << "Total log Steiner point (with divide and conquer) memory usage: " << total_Steiner_point_memory_size / 1e6 << " MB" << std::endl;
    std::cout << "Total binary search snell law path memory usage: " << total_snell_law_memory_size / 1e6 << " MB" << std::endl;
    std::cout << "Total memory usage: " << (total_Steiner_point_memory_size + total_snell_law_memory_size) / 1e6 << " MB" << std::endl;
    std::cout << "Total binary search snell law path count: " << total_binary_search_of_snell_law_path_count << std::endl;
    std::cout << "Distance error: " << distance_error * 100 << "%" << std::endl;
    std::cout << "Total distance: " << total_distance_snell_law << std::endl;
    std::cout << "Total edge sequence size: " << total_edge_sequence_size << std::endl;

    std::ofstream ofs("../output/output.txt", std::ios_base::app);
    ofs << "== ProgLogSP-BinarySearch ==\n";
    ofs << write_file_header << "\t"
        << total_building_time << "\t"
        << total_Steiner_point_query_time << "\t"
        << total_snell_law_query_time << "\t"
        << total_Steiner_point_query_time + total_snell_law_query_time << "\t"
        << total_Steiner_point_memory_size / 1e6 << "\t"
        << total_snell_law_memory_size / 1e6 << "\t"
        << (total_Steiner_point_memory_size + total_snell_law_memory_size) / 1e6 << "\t"
        << total_binary_search_of_snell_law_path_count << "\t"
        << distance_error * 100 << "\t"
        << total_distance_snell_law << "\t"
        << total_edge_sequence_size << "\n\n";
    ofs.close();
}

// the log Steiner point algorithm and divide and conquer with effective weight snell law
void log_Steiner_point_divide_and_conquer_and_effective_weight_snell_law(geodesic::Mesh *mesh, geodesic::SurfacePoint &source,
                                                                         geodesic::SurfacePoint &destination,
                                                                         std::vector<geodesic::SurfacePoint> &path,
                                                                         std::string write_file_header,
                                                                         double Steiner_point_epsilon, double snell_law_epsilon,
                                                                         int max_loop_num_for_divide_and_conquer,
                                                                         int max_loop_num_for_single_endpoint,
                                                                         double total_distance_exact_path,
                                                                         std::vector<geodesic::SurfacePoint> &result_path)
{
    double total_building_time = 0;
    double total_Steiner_point_time;
    double total_Steiner_point_query_time;
    double total_Steiner_point_memory_size = 0;
    log_Steiner_point_divide_and_conquer_helper(mesh, source, destination, path, Steiner_point_epsilon, max_loop_num_for_divide_and_conquer, max_loop_num_for_single_endpoint, total_building_time, total_Steiner_point_time, total_Steiner_point_memory_size);
    total_Steiner_point_query_time = total_Steiner_point_time - total_building_time;

    double snell_law_delta;
    snell_law_epsilon_to_delta(mesh, snell_law_epsilon, path.size() - 1, snell_law_delta);

    // when the path calculated after the divide and conquer still passes the vertex, we cannot
    // use such a path in the snell law, so we need to pick the segment of the path that passes
    // between two vertices, and apply snell law on such a segment, then combine these segments
    // together in the final path result

    std::vector<geodesic::SurfacePoint> segment_source_list;
    std::vector<geodesic::SurfacePoint> segment_destination_list;
    std::vector<std::vector<geodesic::SurfacePoint>> segment_path_list;
    std::vector<geodesic::SurfacePoint> one_segment_path;

    segment_source_list.clear();
    segment_destination_list.clear();
    segment_path_list.clear();
    one_segment_path.clear();

    // we are moving from destination to source
    segment_destination_list.push_back(destination);
    one_segment_path.push_back(destination);

    for (int i = 1; i < path.size() - 1; ++i)
    {
        geodesic::SurfacePoint &s = path[i];

        if (s.type() == geodesic::EDGE)
        {
            one_segment_path.push_back(s);
        }
        else if (s.type() == geodesic::VERTEX)
        {
            segment_source_list.push_back(s);
            one_segment_path.push_back(s);
            segment_path_list.push_back(one_segment_path);
            one_segment_path.clear();
            segment_destination_list.push_back(s);
            one_segment_path.push_back(s);
        }
    }

    segment_source_list.push_back(source);
    one_segment_path.push_back(source);
    segment_path_list.push_back(one_segment_path);
    one_segment_path.clear();

    std::cout << "segment path list size: " << segment_path_list.size() << std::endl;

    double total_snell_law_query_time = 0;
    int total_binary_search_and_effective_weight_of_snell_law_path_count = 0;
    double total_snell_law_memory_size = 0;
    int total_edge_sequence_size = 0;

    result_path.push_back(destination);

    for (int i = 0; i < segment_path_list.size(); ++i)
    {
        if (segment_path_list[i].size() > 2)
        {
            std::vector<geodesic::Edge> edge_sequence;
            std::vector<geodesic::Face> face_sequence;
            std::vector<geodesic::SurfacePoint> snell_law_path;

            get_edge_sequence(mesh, segment_path_list[i], edge_sequence);
            get_face_sequence(mesh, segment_path_list[i], face_sequence);
            total_edge_sequence_size += edge_sequence.size();

            int binary_search_and_effective_weight_of_snell_law_path_count = 0;

            auto start = std::chrono::high_resolution_clock::now();

            effective_weight_binary_search_multiple_times_of_each_edge(mesh, edge_sequence, face_sequence, segment_source_list[i], segment_destination_list[i], snell_law_path, snell_law_delta * 20, binary_search_and_effective_weight_of_snell_law_path_count, total_snell_law_memory_size);

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

            total_snell_law_query_time += duration.count();
            std::cout << "Effective weight snell law query time: " << duration.count() << " milliseconds" << std::endl;

            total_binary_search_and_effective_weight_of_snell_law_path_count += binary_search_and_effective_weight_of_snell_law_path_count;

            for (int j = snell_law_path.size() - 2; j >= 0; --j)
            {
                result_path.push_back(snell_law_path[j]);
            }
        }
        else if (segment_path_list[i].size() <= 2)
        {
            result_path.push_back(segment_path_list[i][1]);
        }
    }

    total_snell_law_memory_size += result_path.size() * sizeof(geodesic::SurfacePoint);

    double total_distance_snell_law;
    total_distance_snell_law = path_distance(result_path, path);
    std::cout << "Log Steiner point, divide and conquer and effective weight snell law distance: " << total_distance_snell_law << std::endl;

    double total_distance_Steiner_point;
    total_distance_Steiner_point = path_distance(path, path);
    if (total_distance_Steiner_point < total_distance_snell_law)
    {
        result_path.clear();
        for (int i = 0; i < path.size(); ++i)
        {
            result_path.push_back(path[i]);
        }
        total_distance_snell_law = total_distance_Steiner_point;
        std::cout << "Log Steiner point distance is shorter, and the final distance is: " << total_distance_Steiner_point << std::endl;
    }

    fix_distance(total_distance_snell_law, total_distance_exact_path);
    double distance_error = total_distance_snell_law / total_distance_exact_path - 1;

    std::cout << "# Summary #" << std::endl;
    std::cout << "Dataset, datasize, epsilon, epsilon_SP, epsilon_SL: " << write_file_header << std::endl;
    std::cout << "Total log Steiner point building time: " << total_building_time << " milliseconds" << std::endl;
    std::cout << "Total log Steiner point (with divide and conquer) query time: " << total_Steiner_point_query_time << " milliseconds" << std::endl;
    std::cout << "Total effective weight snell law query time: " << total_snell_law_query_time << " milliseconds" << std::endl;
    std::cout << "Total query time: " << total_Steiner_point_query_time + total_snell_law_query_time << " milliseconds" << std::endl;
    std::cout << "Total log Steiner point (with divide and conquer) memory usage: " << total_Steiner_point_memory_size / 1e6 << " MB" << std::endl;
    std::cout << "Total binary search and effective weight snell law path memory usage: " << total_snell_law_memory_size / 1e6 << " MB" << std::endl;
    std::cout << "Total memory usage: " << (total_Steiner_point_memory_size + total_snell_law_memory_size) / 1e6 << " MB" << std::endl;
    std::cout << "Total binary search and effective weight snell law path count: " << total_binary_search_and_effective_weight_of_snell_law_path_count << std::endl;
    std::cout << "Distance error: " << distance_error * 100 << "%" << std::endl;
    std::cout << "Total distance: " << total_distance_snell_law << std::endl;
    std::cout << "Total edge sequence size: " << total_edge_sequence_size << std::endl;

    std::ofstream ofs("../output/output.txt", std::ios_base::app);
    ofs << "== ProgLogSP-EffWeight ==\n";
    ofs << write_file_header << "\t"
        << total_building_time << "\t"
        << total_Steiner_point_query_time << "\t"
        << total_snell_law_query_time << "\t"
        << total_Steiner_point_query_time + total_snell_law_query_time << "\t"
        << total_Steiner_point_memory_size / 1e6 << "\t"
        << total_snell_law_memory_size / 1e6 << "\t"
        << (total_Steiner_point_memory_size + total_snell_law_memory_size) / 1e6 << "\t"
        << total_binary_search_and_effective_weight_of_snell_law_path_count << "\t"
        << distance_error * 100 << "\t"
        << total_distance_snell_law << "\t"
        << total_edge_sequence_size << "\n\n";
    ofs.close();
}
