#ifndef GEODESIC_ALGORITHM_GRAPH_BASE_010907
#define GEODESIC_ALGORITHM_GRAPH_BASE_010907

#include "geodesic_algorithm_base.h"
#include "geodesic_mesh_elements.h"
#include <vector>
#include <set>
#include <assert.h>

namespace geodesic
{

	template <class Node>
	class GeodesicAlgorithmGraphBase : public GeodesicAlgorithmBase
	{
	public:
		typedef Node *node_pointer;

		GeodesicAlgorithmGraphBase(geodesic::Mesh *mesh) : GeodesicAlgorithmBase(mesh){};

		~GeodesicAlgorithmGraphBase(){};

		void propagate(std::vector<SurfacePoint> &sources,
					   double max_propagation_distance,		   // propagation algorithm stops after reaching the certain distance from the source
					   std::vector<SurfacePoint> *stop_points, // or after ensuring that all the stop_points are covered
					   std::unordered_map<int, double> input_dist,
					   std::unordered_map<int, int> input_prev_node,
					   std::unordered_map<int, int> input_src_index,
					   std::unordered_map<int, double> &output_dist,
					   std::unordered_map<int, int> &output_prev_node,
					   std::unordered_map<int, int> &output_src_index);

		void trace_back(SurfacePoint &destination, // trace back piecewise-linear path
						std::vector<SurfacePoint> &path);

		unsigned best_source(SurfacePoint &point, // quickly find what source this point belongs to and what is the distance to this source
							 double &best_source_distance);

		void print_statistics()
		{
			GeodesicAlgorithmBase::print_statistics();

			double memory = m_nodes.size() * sizeof(Node) / (m_removing_value + 1);
			std::cout << "uses about " << memory / 1e6 << "MB of memory" << std::endl;
		}

		double get_memory();

		double get_edge_pass_face_weight(node_pointer node, SurfacePoint *point)
		{
			// if both two points are the on the vertex, then they are on the original mesh edge,
			// so the weight of this path is the min face weight of the faces that connects to this edge
			if (node->type() == geodesic::VERTEX && point->type() == geodesic::VERTEX)
			{
				geodesic::edge_pointer shared_edge;
				bool jump_out = false;
				for (unsigned i = 0; i < node->base_element()->adjacent_edges().size() && !jump_out; ++i)
				{
					for (unsigned j = 0; j < point->base_element()->adjacent_edges().size() && !jump_out; ++j)
					{
						if (node->base_element()->adjacent_edges()[i]->id() == point->base_element()->adjacent_edges()[j]->id())
						{
							shared_edge = node->base_element()->adjacent_edges()[i];
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
			// then they are on the original mesh edge, so the weight of this path is the min face weight of the faces that
			// connect to this edge
			else if ((node->type() == geodesic::VERTEX && point->type() == geodesic::EDGE && node->base_element()->id() == point->base_element()->adjacent_vertices()[0]->id()) ||
					 (node->type() == geodesic::VERTEX && point->type() == geodesic::EDGE && node->base_element()->id() == point->base_element()->adjacent_vertices()[1]->id()))
			{
				std::vector<double> face_weight_list;
				face_weight_list.clear();
				for (unsigned i = 0; i < point->base_element()->adjacent_faces().size(); ++i)
				{
					face_weight_list.push_back(point->base_element()->adjacent_faces()[i]->weight());
				}
				assert(face_weight_list.size());
				return *std::min_element(face_weight_list.begin(), face_weight_list.end());
			}

			// if the first point is on the edge and the second point is on the vertex and the second point is on this edge
			// then they are on the original mesh edge, so the weight of this path is the min face weight of the faces that
			// connect to this edge
			else if ((node->type() == geodesic::EDGE && point->type() == geodesic::VERTEX && point->base_element()->id() == node->base_element()->adjacent_vertices()[0]->id()) ||
					 (node->type() == geodesic::EDGE && point->type() == geodesic::VERTEX && point->base_element()->id() == node->base_element()->adjacent_vertices()[1]->id()))
			{
				std::vector<double> face_weight_list;
				face_weight_list.clear();
				for (unsigned i = 0; i < node->base_element()->adjacent_faces().size(); ++i)
				{
					face_weight_list.push_back(node->base_element()->adjacent_faces()[i]->weight());
				}
				assert(face_weight_list.size());
				return *std::min_element(face_weight_list.begin(), face_weight_list.end());
			}

			// if two points are both on the same original mesh edge,
			// the weight of this path is the min face weight of the faces that connects to this edge
			else if (node->type() == geodesic::EDGE && point->type() == geodesic::EDGE && node->base_element()->id() == point->base_element()->id())
			{
				std::vector<double> face_weight_list;
				face_weight_list.clear();
				for (unsigned i = 0; i < point->base_element()->adjacent_faces().size(); ++i)
				{
					face_weight_list.push_back(point->base_element()->adjacent_faces()[i]->weight());
				}
				assert(face_weight_list.size());
				return *std::min_element(face_weight_list.begin(), face_weight_list.end());
			}

			// if one point is on the vertex and one point is on the edge, but that vertex is not on this edge
			// or if the two points are on two different edges
			// then the path weight is the face weight
			else if ((node->type() == geodesic::VERTEX && point->type() == geodesic::EDGE && node->base_element()->id() != point->base_element()->adjacent_vertices()[0]->id() && node->base_element()->id() != point->base_element()->adjacent_vertices()[1]->id()) ||
					 (node->type() == geodesic::EDGE && point->type() == geodesic::VERTEX && point->base_element()->id() != node->base_element()->adjacent_vertices()[0]->id() && point->base_element()->id() != node->base_element()->adjacent_vertices()[1]->id()) ||
					 (node->type() == geodesic::EDGE && point->type() == geodesic::EDGE && node->base_element()->id() != point->base_element()->id()))
			{
				bool found_edge_in_face = false;
				for (unsigned i = 0; i < point->base_element()->adjacent_faces().size(); ++i)
				{
					geodesic::face_pointer f = point->base_element()->adjacent_faces()[i];

					if (((vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[1]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[1]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[1]->z(),
							  node->x(), node->y(), node->z()) &&
						  vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[2]->z(),
							  point->getx(), point->gety(), point->getz())) ||
						 (vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[1]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[1]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[1]->z(),
							  node->x(), node->y(), node->z()) &&
						  vertex_in_edge(
							  f->adjacent_vertices()[1]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[1]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[1]->z(), f->adjacent_vertices()[2]->z(),
							  point->getx(), point->gety(), point->getz()))) ||
						((vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[2]->z(),
							  node->x(), node->y(), node->z()) &&
						  vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[1]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[1]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[1]->z(),
							  point->getx(), point->gety(), point->getz())) ||
						 (vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[2]->z(),
							  node->x(), node->y(), node->z()) &&
						  vertex_in_edge(
							  f->adjacent_vertices()[1]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[1]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[1]->z(), f->adjacent_vertices()[2]->z(),
							  point->getx(), point->gety(), point->getz()))) ||
						((vertex_in_edge(
							  f->adjacent_vertices()[1]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[1]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[1]->z(), f->adjacent_vertices()[2]->z(),
							  node->x(), node->y(), node->z()) &&
						  vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[1]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[1]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[1]->z(),
							  point->getx(), point->gety(), point->getz())) ||
						 (vertex_in_edge(
							  f->adjacent_vertices()[1]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[1]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[1]->z(), f->adjacent_vertices()[2]->z(),
							  node->x(), node->y(), node->z()) &&
						  vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[2]->z(),
							  point->getx(), point->gety(), point->getz()))))
					{
						found_edge_in_face = true;
						return f->weight();
					}
				}
			}

			return 100;
		}

		double get_edge_pass_face_weight_2(SurfacePoint *point1, SurfacePoint *point2)
		{
			// if both two points are the on the vertex, then they are on the original mesh edge,
			// so the weight of this path is the max face weight of the faces that connects to this edge
			if (point1->type() == geodesic::VERTEX && point2->type() == geodesic::VERTEX)
			{
				geodesic::edge_pointer shared_edge;
				bool jump_out = false;
				for (unsigned i = 0; i < point1->base_element()->adjacent_edges().size() && !jump_out; ++i)
				{
					for (unsigned j = 0; j < point2->base_element()->adjacent_edges().size() && !jump_out; ++j)
					{
						if (point1->base_element()->adjacent_edges()[i]->id() == point2->base_element()->adjacent_edges()[j]->id())
						{
							shared_edge = point1->base_element()->adjacent_edges()[i];
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
			else if ((point1->type() == geodesic::VERTEX && point2->type() == geodesic::EDGE && point1->base_element()->id() == point2->base_element()->adjacent_vertices()[0]->id()) ||
					 (point1->type() == geodesic::VERTEX && point2->type() == geodesic::EDGE && point1->base_element()->id() == point2->base_element()->adjacent_vertices()[1]->id()))
			{
				std::vector<double> face_weight_list;
				face_weight_list.clear();
				for (unsigned i = 0; i < point2->base_element()->adjacent_faces().size(); ++i)
				{
					face_weight_list.push_back(point2->base_element()->adjacent_faces()[i]->weight());
				}
				assert(face_weight_list.size());
				return *std::min_element(face_weight_list.begin(), face_weight_list.end());
			}

			// if the first point is on the edge and the second point is on the vertex and the second point is on this edge
			// then they are on the original mesh edge, so the weight of this path is the max face weight of the faces that
			// connect to this edge
			else if ((point1->type() == geodesic::EDGE && point2->type() == geodesic::VERTEX && point2->base_element()->id() == point1->base_element()->adjacent_vertices()[0]->id()) ||
					 (point1->type() == geodesic::EDGE && point2->type() == geodesic::VERTEX && point2->base_element()->id() == point1->base_element()->adjacent_vertices()[1]->id()))
			{
				std::vector<double> face_weight_list;
				face_weight_list.clear();
				for (unsigned i = 0; i < point1->base_element()->adjacent_faces().size(); ++i)
				{
					face_weight_list.push_back(point1->base_element()->adjacent_faces()[i]->weight());
				}
				assert(face_weight_list.size());
				return *std::min_element(face_weight_list.begin(), face_weight_list.end());
			}

			// if two points are both on the same original mesh edge,
			// the weight of this path is the min face weight of the faces that connects to this edge
			else if (point1->type() == geodesic::EDGE && point2->type() == geodesic::EDGE && point1->base_element()->id() == point2->base_element()->id())
			{
				std::vector<double> face_weight_list;
				face_weight_list.clear();
				for (unsigned i = 0; i < point2->base_element()->adjacent_faces().size(); ++i)
				{
					face_weight_list.push_back(point2->base_element()->adjacent_faces()[i]->weight());
				}
				assert(face_weight_list.size());
				return *std::min_element(face_weight_list.begin(), face_weight_list.end());
			}

			// if one point is on the vertex and one point is on the edge, but that vertex is not on this edge
			// or if the two points are on two different edges
			// then the path weight is the face weight
			else if ((point1->type() == geodesic::VERTEX && point2->type() == geodesic::EDGE && point1->base_element()->id() != point2->base_element()->adjacent_vertices()[0]->id() && point1->base_element()->id() != point2->base_element()->adjacent_vertices()[1]->id()) ||
					 (point1->type() == geodesic::EDGE && point2->type() == geodesic::VERTEX && point2->base_element()->id() != point1->base_element()->adjacent_vertices()[0]->id() && point2->base_element()->id() != point1->base_element()->adjacent_vertices()[1]->id()) ||
					 (point1->type() == geodesic::EDGE && point2->type() == geodesic::EDGE))
			{
				for (unsigned i = 0; i < point2->base_element()->adjacent_faces().size(); ++i)
				{
					geodesic::face_pointer f = point2->base_element()->adjacent_faces()[i];

					if (((vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[1]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[1]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[1]->z(),
							  point1->getx(), point1->gety(), point1->getz()) &&
						  vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[2]->z(),
							  point2->getx(), point2->gety(), point2->getz())) ||
						 (vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[1]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[1]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[1]->z(),
							  point1->getx(), point1->gety(), point1->getz()) &&
						  vertex_in_edge(
							  f->adjacent_vertices()[1]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[1]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[1]->z(), f->adjacent_vertices()[2]->z(),
							  point2->getx(), point2->gety(), point2->getz()))) ||
						((vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[2]->z(),
							  point1->getx(), point1->gety(), point1->getz()) &&
						  vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[1]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[1]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[1]->z(),
							  point2->getx(), point2->gety(), point2->getz())) ||
						 (vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[2]->z(),
							  point1->getx(), point1->gety(), point1->getz()) &&
						  vertex_in_edge(
							  f->adjacent_vertices()[1]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[1]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[1]->z(), f->adjacent_vertices()[2]->z(),
							  point2->getx(), point2->gety(), point2->getz()))) ||
						((vertex_in_edge(
							  f->adjacent_vertices()[1]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[1]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[1]->z(), f->adjacent_vertices()[2]->z(),
							  point1->getx(), point1->gety(), point1->getz()) &&
						  vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[1]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[1]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[1]->z(),
							  point2->getx(), point2->gety(), point2->getz())) ||
						 (vertex_in_edge(
							  f->adjacent_vertices()[1]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[1]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[1]->z(), f->adjacent_vertices()[2]->z(),
							  point1->getx(), point1->gety(), point1->getz()) &&
						  vertex_in_edge(
							  f->adjacent_vertices()[0]->x(), f->adjacent_vertices()[2]->x(),
							  f->adjacent_vertices()[0]->y(), f->adjacent_vertices()[2]->y(),
							  f->adjacent_vertices()[0]->z(), f->adjacent_vertices()[2]->z(),
							  point2->getx(), point2->gety(), point2->getz()))))
					{
						return f->weight();
					}
				}
			}
			return 100;
		}

		double max(double a, double b)
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
		double min(double a, double b)
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
			return ((vertex_x - min(edge_1_x, edge_2_x) >= -1e-6) &&
					(vertex_x - max(edge_1_x, edge_2_x) <= 1e-6) &&
					(vertex_y - min(edge_1_y, edge_2_y) >= -1e-6) &&
					(vertex_y - max(edge_1_y, edge_2_y) <= 1e-6) &&
					(vertex_z - min(edge_1_z, edge_2_z) >= -1e-6) &&
					(vertex_z - max(edge_1_z, edge_2_z) <= 1e-6));
		}

	protected:
		unsigned node_index(vertex_pointer v) // gives index of the node that corresponds to this vertex
		{
			return v->id();
		};

		void set_sources(std::vector<SurfacePoint> &sources)
		{
			m_sources = sources;
		}

		node_pointer best_first_node(SurfacePoint &point, double &best_total_distance)
		{
			node_pointer best_node = NULL;
			if (point.type() == VERTEX)
			{
				vertex_pointer v = (vertex_pointer)point.base_element();
				best_node = &m_nodes[node_index(v)];
				best_total_distance = best_node->distance_from_source();
			}
			else
			{
				std::vector<node_pointer> possible_nodes;
				list_nodes_visible_from_source(point.base_element(), possible_nodes);

				best_total_distance = GEODESIC_INF;
				for (unsigned i = 0; i < possible_nodes.size(); ++i)
				{
					node_pointer node = possible_nodes[i];

					double distance_from_dest = get_edge_pass_face_weight(node, &point) * node->distance(&point);
					if (node->distance_from_source() + distance_from_dest < best_total_distance)
					{
						best_total_distance = node->distance_from_source() + distance_from_dest;
						best_node = node;
					}
				}
			}

			// assert(best_node);
			// assert(best_total_distance<GEODESIC_INF);
			if (best_total_distance > m_propagation_distance_stopped) // result is unreliable
			{
				best_total_distance = GEODESIC_INF;
				return NULL;
			}
			else
			{
				return best_node;
			}
		}; // quickly find what node will be the next one in geodesic path

		bool check_stop_conditions(unsigned &index); // check when propagation should stop

		virtual void list_nodes_visible_from_source(MeshElementBase *p,
													std::vector<node_pointer> &storage) = 0; // list all nodes that belong to this mesh element

		virtual void list_nodes_visible_from_node(node_pointer node, // list all nodes that belong to this mesh element
												  std::vector<node_pointer> &storage,
												  std::vector<double> &distances,
												  double threshold_distance) = 0; // list only the nodes whose current distance is larger than the threshold

		std::vector<Node> m_nodes; // nodes of the graph
		std::vector<Node> m_face;  // face of the graph
		int m_removing_value;

		typedef std::set<node_pointer, Node> queue_type;
		queue_type m_queue;

		std::vector<SurfacePoint> m_sources; // for simplicity, we keep sources as they are
	};

	// propagation algorithm stops after reaching the certain distance from the source or after ensuring that all the stop_points are covered
	template <class Node>
	void GeodesicAlgorithmGraphBase<Node>::propagate(std::vector<SurfacePoint> &sources,
													 double max_propagation_distance,
													 std::vector<SurfacePoint> *stop_points,
													 std::unordered_map<int, double> input_dist,
													 std::unordered_map<int, int> input_prev_node,
													 std::unordered_map<int, int> input_src_index,
													 std::unordered_map<int, double> &output_dist,
													 std::unordered_map<int, int> &output_prev_node,
													 std::unordered_map<int, int> &output_src_index)
	{
		set_stop_conditions(stop_points, max_propagation_distance);
		set_sources(sources);

		m_queue.clear();
		m_propagation_distance_stopped = GEODESIC_INF;
		for (unsigned i = 0; i < m_nodes.size(); ++i)
		{
			m_nodes[i].clear();
		}

		for (unsigned i = 0; i < m_nodes.size(); ++i)
		{
			if (i != m_nodes[i].node_id())
			{
				std::cout << "i: " << i << ", node id: " << m_nodes[i].node_id() << std::endl;
			}
			assert(i == m_nodes[i].node_id());

			if (input_dist.count(i) != 0)
			{
				m_nodes[i].distance_from_source() = input_dist[i];
			}
			if (input_prev_node.count(i) != 0)
			{
				m_nodes[i].previous() = &m_nodes[input_prev_node[i]];
			}
			if (input_src_index.count(i) != 0)
			{
				m_nodes[i].source_index() = input_src_index[i];
			}
		}

		clock_t start = clock();

		std::vector<node_pointer> visible_nodes; // initialize vertices directly visible from sources
		for (unsigned i = 0; i < m_sources.size(); ++i)
		{
			SurfacePoint *source = &m_sources[i];
			list_nodes_visible_from_source(source->base_element(),
										   visible_nodes);

			for (unsigned j = 0; j < visible_nodes.size(); ++j)
			{
				node_pointer node = visible_nodes[j];
				double distance = get_edge_pass_face_weight(node, source) * node->distance(source);
				if (distance < node->distance_from_source())
				{
					node->distance_from_source() = distance;
					node->source_index() = i;
					node->previous() = NULL;
					m_queue.insert(node);
				}
			}
			visible_nodes.clear();
		}

		unsigned counter = 0;
		unsigned satisfied_index = 0;

		std::vector<double> distances_between_nodes;
		while (!m_queue.empty()) // main cycle
		{
			if (counter++ % 10 == 0) // check if we covered all required vertices
			{
				if (check_stop_conditions(satisfied_index))
				{
					break;
				}
			}

			node_pointer min_node = *m_queue.begin();
			m_queue.erase(m_queue.begin());
			assert(min_node->distance_from_source() < GEODESIC_INF);

			visible_nodes.clear();
			distances_between_nodes.clear();
			list_nodes_visible_from_node(min_node,
										 visible_nodes,
										 distances_between_nodes,
										 min_node->distance_from_source());
			for (unsigned i = 0; i < visible_nodes.size(); ++i) // update all the adgecent vertices
			{
				node_pointer next_node = visible_nodes[i];

				if (next_node->distance_from_source() > min_node->distance_from_source() +
															distances_between_nodes[i])
				{
					next_node->distance_from_source() = min_node->distance_from_source() +
														distances_between_nodes[i];
					next_node->source_index() = min_node->source_index();
					next_node->previous() = min_node;
					m_queue.insert(next_node);
				}
			}
		}

		m_propagation_distance_stopped = m_queue.empty() ? GEODESIC_INF : (*m_queue.begin())->distance_from_source();
		clock_t finish = clock();
		m_time_consumed = (static_cast<double>(finish) - static_cast<double>(start)) / CLOCKS_PER_SEC;

		for (unsigned i = 0; i < m_nodes.size(); ++i)
		{
			if (m_nodes[i].previous() != NULL)
			{
				output_dist[i] = m_nodes[i].distance_from_source();
				output_prev_node[i] = m_nodes[i].previous()->node_id();
				output_src_index[i] = m_nodes[i].source_index();
			}
		}
	}

	template <class Node>
	inline bool GeodesicAlgorithmGraphBase<Node>::check_stop_conditions(unsigned &index)
	{
		double queue_min_distance = (*m_queue.begin())->distance_from_source();

		if (queue_min_distance < m_max_propagation_distance)
		{
			return false;
		}

		while (index < m_stop_vertices.size())
		{
			vertex_pointer v = m_stop_vertices[index].first;
			Node &node = m_nodes[node_index(v)];
			if (queue_min_distance < node.distance_from_source() + m_stop_vertices[index].second)
			{
				return false;
			}
			++index;
		}
		return true;
	}

	template <class Node>
	inline void GeodesicAlgorithmGraphBase<Node>::trace_back(SurfacePoint &destination, // trace back piecewise-linear path
															 std::vector<SurfacePoint> &path)
	{
		path.clear();

		double total_path_length;
		node_pointer node = best_first_node(destination, total_path_length);

		if (total_path_length > GEODESIC_INF / 2.0) // unable to find the path
		{
			return;
		}

		path.push_back(destination);

		if (node->distance(&destination) > 1e-50)
		{
			path.push_back(node->surface_point());
		}

		while (node->previous()) // follow the path
		{
			node = node->previous();
			path.push_back(node->surface_point());
		}

		SurfacePoint &source = m_sources[node->source_index()]; // add source to the path if it is not already there
		if (node->distance(&source) > 1e-50)
		{
			path.push_back(source);
		}
	}

	template <class Node>
	inline double GeodesicAlgorithmGraphBase<Node>::get_memory()
	{
		double memory_size = m_nodes.size() * sizeof(Node) / (m_removing_value + 1);
		return memory_size;
	}

	template <class Node>
	inline unsigned GeodesicAlgorithmGraphBase<Node>::best_source(SurfacePoint &point, // quickly find what source this point belongs to and what is the distance to this source
																  double &best_source_distance)
	{
		node_pointer node = best_first_node(point, best_source_distance);
		return node ? node->source_index() : 0;
	};

} // namespace geodesic

#endif // GEODESIC_ALGORITHM_GRAPH_BASE_010907
