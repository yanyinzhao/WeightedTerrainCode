// Copyright (C) 2008 Danil Kirsanov, MIT License

#ifndef GEODESIC_ALGORITHM_BASE_122806
#define GEODESIC_ALGORITHM_BASE_122806

#include "geodesic_mesh.h"
#include "geodesic_constants_and_simple_functions.h"
#include <iostream>
#include <ctime>
#include <unordered_map>

namespace geodesic
{

	class GeodesicAlgorithmBase
	{
	public:
		enum AlgorithmType
		{
			SUBDIVISION,
			SUBDIVISIONLOG,
			UNDEFINED_ALGORITHM
		};

		GeodesicAlgorithmBase(geodesic::Mesh *mesh) : m_type(UNDEFINED_ALGORITHM),
													  m_max_propagation_distance(1e100),
													  m_mesh(mesh){};
		virtual ~GeodesicAlgorithmBase(){};

		virtual void propagate(std::vector<SurfacePoint> &sources,
							   double max_propagation_distance,
							   std::vector<SurfacePoint> *stop_points,
							   std::unordered_map<int, double> input_dist,
							   std::unordered_map<int, int> input_prev_node,
							   std::unordered_map<int, int> input_src_index,
							   std::unordered_map<int, double> &output_dist,
							   std::unordered_map<int, int> &output_prev_node,
							   std::unordered_map<int, int> &output_src_index) = 0;

		virtual void trace_back(SurfacePoint &destination, // trace back piecewise-linear path
								std::vector<SurfacePoint> &path) = 0;

		virtual double get_memory() = 0;

		void geodesic(SurfacePoint &source,
					  SurfacePoint &destination,
					  std::vector<SurfacePoint> &path,
					  std::unordered_map<int, double> input_dist,
					  std::unordered_map<int, int> input_prev_node,
					  std::unordered_map<int, int> input_src_index,
					  std::unordered_map<int, double> &output_dist,
					  std::unordered_map<int, int> &output_prev_node,
					  std::unordered_map<int, int> &output_src_index); // lazy people can find geodesic path with one function call

		double return_memory();

		virtual unsigned best_source(SurfacePoint &point, // after propagation step is done, quickly find what source this point belongs to and what is the distance to this source
									 double &best_source_distance) = 0;

		virtual void print_statistics() // print info about timing and memory usage in the propagation step of the algorithm
		{
			std::cout << "propagation step took " << m_time_consumed << " seconds " << std::endl;
		};

		AlgorithmType type() { return m_type; };

		virtual std::string name();

		geodesic::Mesh *mesh() { return m_mesh; };

	protected:
		void set_stop_conditions(std::vector<SurfacePoint> *stop_points,
								 double stop_distance);
		double stop_distance()
		{
			return m_max_propagation_distance;
		}

		AlgorithmType m_type; // type of the algorithm

		typedef std::pair<vertex_pointer, double> stop_vertex_with_distace_type;
		std::vector<stop_vertex_with_distace_type> m_stop_vertices; // algorithm stops propagation after covering certain vertices
		double m_max_propagation_distance;							// or reaching the certain distance

		geodesic::Mesh *m_mesh;

		double m_time_consumed;				   // how much time does the propagation step takes
		double m_propagation_distance_stopped; // at what distance (if any) the propagation algorithm stopped
	};

	inline double length(std::vector<SurfacePoint> &path)
	{
		double length = 0;
		if (!path.empty())
		{
			for (unsigned i = 0; i < path.size() - 1; ++i)
			{
				length += path[i].distance(&path[i + 1]);
			}
		}
		return length;
	}

	inline void print_info_about_path(std::vector<SurfacePoint> &path)
	{
		std::cout << "number of the points in the path = " << path.size()
				  << ", length of the path = " << length(path)
				  << std::endl;
	}

	inline std::string GeodesicAlgorithmBase::name()
	{
		switch (m_type)
		{
		case SUBDIVISION:
			return "subdivision";
		case SUBDIVISIONLOG:
			return "subdivision_log";
		default:
		case UNDEFINED_ALGORITHM:
			return "undefined";
		}
	}

	// now using this geodesic function
	inline void GeodesicAlgorithmBase::geodesic(SurfacePoint &source,
												SurfacePoint &destination,
												std::vector<SurfacePoint> &path,
												std::unordered_map<int, double> input_dist,
												std::unordered_map<int, int> input_prev_node,
												std::unordered_map<int, int> input_src_index,
												std::unordered_map<int, double> &output_dist,
												std::unordered_map<int, int> &output_prev_node,
												std::unordered_map<int, int> &output_src_index) // lazy people can find geodesic path with one function call
	{
		std::vector<SurfacePoint> sources(1, source);
		std::vector<SurfacePoint> stop_points(1, destination);
		double const max_propagation_distance = GEODESIC_INF;

		propagate(sources,
				  max_propagation_distance,
				  &stop_points,
				  input_dist, input_prev_node,
				  input_src_index, output_dist,
				  output_prev_node, output_src_index);

		trace_back(destination, path);
	}

	inline double GeodesicAlgorithmBase::return_memory()
	{
		return get_memory();
	}

	inline void GeodesicAlgorithmBase::set_stop_conditions(std::vector<SurfacePoint> *stop_points,
														   double stop_distance)
	{
		m_max_propagation_distance = stop_distance;

		if (!stop_points)
		{
			m_stop_vertices.clear();
			return;
		}

		m_stop_vertices.resize(stop_points->size());

		std::vector<vertex_pointer> possible_vertices;
		for (unsigned i = 0; i < stop_points->size(); ++i)
		{
			SurfacePoint *point = &(*stop_points)[i];

			possible_vertices.clear();
			m_mesh->closest_vertices(point, &possible_vertices);

			vertex_pointer closest_vertex = NULL;
			double min_distance = 1e100;
			for (unsigned j = 0; j < possible_vertices.size(); ++j)
			{
				double distance = point->distance(possible_vertices[j]);
				if (distance < min_distance)
				{
					min_distance = distance;
					closest_vertex = possible_vertices[j];
				}
			}
			assert(closest_vertex);

			m_stop_vertices[i].first = closest_vertex;
			m_stop_vertices[i].second = min_distance;
		}
	}

} // namespace geodesic

#endif // GEODESIC_ALGORITHM_BASE_122806
