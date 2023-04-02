// Copyright (C) 2008 Danil Kirsanov, MIT License
#ifndef GEODESIC_ALGORITHM_SUBDIVISION_LOG_122806
#define GEODESIC_ALGORITHM_SUBDIVISION_LOG_122806

#include "geodesic_algorithm_graph_base.h"
#include "geodesic_mesh_elements.h"
#include <vector>
#include <set>
#include <assert.h>

namespace geodesic
{

	class SubdivisionLogNode : public SurfacePoint
	{
		typedef SubdivisionLogNode *node_pointer;

	public:
		SubdivisionLogNode(){};

		template <class Pointer>
		SubdivisionLogNode(Pointer p) : SurfacePoint(p),
										m_previous(NULL),
										m_distance(0.0){};

		template <class Pointer, class Parameter>
		SubdivisionLogNode(Pointer p, Parameter param) : SurfacePoint(p, param),
														 m_previous(NULL),
														 m_distance(0.0){};

		~SubdivisionLogNode(){};

		double &distance_from_source() { return m_distance; };
		node_pointer &previous() { return m_previous; };
		unsigned &source_index() { return m_source_index; };
		int &node_id() { return m_id; };

		void clear()
		{
			m_distance = GEODESIC_INF;
			m_previous = NULL;
		}

		bool operator()(node_pointer const s1, node_pointer const s2) const
		{
			if (s1 == s2)
			{
				return false;
			}
			if (s1->distance_from_source() != s2->distance_from_source())
			{
				return s1->distance_from_source() < s2->distance_from_source();
			}

			if (s1->x() != s2->x()) // two nodes cannot be located in the same space
			{
				return s1->x() < s2->x();
			}
			if (s1->y() != s2->y())
			{
				return s1->y() < s2->y();
			}
			if (s1->z() != s2->z())
			{
				return s1->z() < s2->z();
			}

			assert(0);
			return true;
		};

		SurfacePoint &surface_point() { return static_cast<SurfacePoint &>(*this); };

	private:
		double m_distance;		 // distance to the closest source
		unsigned m_source_index; // closest source index
		node_pointer m_previous; // previous node in the geodesic path
		int m_id;
	};

	class GeodesicAlgorithmSubdivisionLog : public GeodesicAlgorithmGraphBase<SubdivisionLogNode>
	{
		typedef SubdivisionLogNode Node;

	public:
		GeodesicAlgorithmSubdivisionLog(geodesic::Mesh *mesh = NULL,
										unsigned max_subdivision_log_level = 100,
										int removing_value = 2,
										double delta = 1.1,
										double r = 1) : GeodesicAlgorithmGraphBase<Node>(mesh)
		{
			m_type = SUBDIVISIONLOG;

			m_nodes.reserve(mesh->vertices().size());
			for (unsigned i = 0; i < mesh->vertices().size(); ++i)
			{
				vertex_pointer v = &mesh->vertices()[i];

				m_nodes.push_back(Node(v)); //!!
				m_nodes[m_nodes.size() - 1].node_id() = i;
			}

			set_subdivision_log_level(max_subdivision_log_level, delta, r);
			m_removing_value = removing_value;
		};

		~GeodesicAlgorithmSubdivisionLog(){};

		unsigned get_total_Steiner_point_on_edge() { return total_Steiner_point_on_edge; };

		void vector_cross_product(double a_x, double a_y, double a_z,
								  double b_x, double b_y, double b_z,
								  double &result_x, double &result_y, double &result_z)
		{
			result_x = a_y * b_z - a_z * b_y;
			result_y = a_z * b_x - a_x * b_z;
			result_z = a_x * b_y - a_y * b_x;
		}

		void set_subdivision_log_level(unsigned max_subdivision_log_level,
									   double delta,
									   double r)
		{
			// note that in order to test the worst case, we use the worst case to calculate the subdivision log level

			subdivision_log_level_each_edge.reserve(m_mesh->edges().size());
			Steiner_point_on_edge_indexx_map.reserve(m_mesh->edges().size());
			total_Steiner_point_on_edge = 0;

			for (unsigned i = 0; i < m_mesh->edges().size(); ++i)
			{
				Steiner_point_on_edge_indexx_map.push_back(total_Steiner_point_on_edge);

				edge_pointer e = &m_mesh->edges()[i];
				unsigned subdivision_log_level_curr_edge = 0;

				for (unsigned j = 0; j < max_subdivision_log_level; ++j)
				{
					double offset = r * pow(delta, j) / e->length();
					if (offset < 0.5 && offset > 0)
					{
						subdivision_log_level_curr_edge++;
					}
					else
					{
						break;
					}
				}
				// for each edge, we count the Steiner point from two sides
				subdivision_log_level_each_edge.push_back(2 * subdivision_log_level_curr_edge);
				total_Steiner_point_on_edge += 2 * subdivision_log_level_curr_edge;
			}

			m_nodes.resize(m_mesh->vertices().size());
			m_nodes.reserve(m_mesh->vertices().size() + total_Steiner_point_on_edge);

			int node_index = m_mesh->vertices().size();
			for (unsigned i = 0; i < m_mesh->edges().size(); ++i)
			{
				edge_pointer e = &m_mesh->edges()[i];

				for (unsigned j = 0; j < subdivision_log_level_each_edge[i] / 2; ++j)
				{
					double offset = r * pow(delta, j) / e->length();
					assert(offset < 0.5 && offset > 0);
					m_nodes.push_back(Node(e, offset));
					m_nodes[m_nodes.size() - 1].node_id() = node_index;
					node_index++;
				}
				for (unsigned j = 0; j < subdivision_log_level_each_edge[i] / 2; ++j)
				{
					double offset = 1 - (r * pow(delta, j) / e->length());
					assert(offset > 0.5 && offset < 1);
					m_nodes.push_back(Node(e, offset));
					m_nodes[m_nodes.size() - 1].node_id() = node_index;
					node_index++;
				}
			}
			std::cout << "total nodes: " << m_nodes.size() << "\n";
		};

	protected:
		void list_nodes_visible_from_source(MeshElementBase *p,
											std::vector<node_pointer> &storage); // list all nodes that belong to this mesh element

		void list_nodes_visible_from_node(node_pointer node, // list all nodes that belong to this mesh element
										  std::vector<node_pointer> &storage,
										  std::vector<double> &distances,
										  double threshold_distance); // list only the nodes whose current distance is larger than the threshold

		unsigned node_indexx(edge_pointer e)
		{
			return Steiner_point_on_edge_indexx_map[e->id()] + m_mesh->vertices().size();
		};

	private:
		void list_nodes(MeshElementBase *p, // list nodes that belong to this mesh element
						std::vector<node_pointer> &storage,
						double threshold_distance = -1.0); // list only the nodes whose current distance is larger than the threshold

		std::vector<unsigned> subdivision_log_level_each_edge;
		std::vector<unsigned> Steiner_point_on_edge_indexx_map;
		unsigned total_Steiner_point_on_edge;
	};

	inline void GeodesicAlgorithmSubdivisionLog::list_nodes(MeshElementBase *p,
															std::vector<node_pointer> &storage,
															double threshold_distance)
	{

		assert(p->type() != UNDEFINED_POINT);

		if (p->type() == VERTEX)
		{
			vertex_pointer v = static_cast<vertex_pointer>(p);
			node_pointer node = &m_nodes[node_index(v)];
			if (node->distance_from_source() > threshold_distance)
			{
				storage.push_back(node);
			}
		}
		else if (p->type() == EDGE)
		{
			edge_pointer e = static_cast<edge_pointer>(p);
			unsigned node_index = node_indexx(e);
			for (int i = 0; i < subdivision_log_level_each_edge[e->id()] / 2; i = i + 1 + m_removing_value)
			{
				node_pointer node = &m_nodes[node_index];
				node_index = node_index + 1 + m_removing_value;
				if (node->distance_from_source() > threshold_distance)
				{
					storage.push_back(node);
				}
			}
			node_index = node_indexx(e) + subdivision_log_level_each_edge[e->id()] / 2;
			for (int i = 0; i < subdivision_log_level_each_edge[e->id()] / 2; i = i + 1 + m_removing_value)
			{
				node_pointer node = &m_nodes[node_index];
				node_index = node_index + 1 + m_removing_value;
				if (node->distance_from_source() > threshold_distance)
				{
					storage.push_back(node);
				}
			}
		}

		// FACE has no nodes
	}

	void GeodesicAlgorithmSubdivisionLog::list_nodes_visible_from_source(MeshElementBase *p,
																		 std::vector<node_pointer> &storage)
	{
		assert(p->type() != UNDEFINED_POINT);

		if (p->type() == FACE)
		{
			face_pointer f = static_cast<face_pointer>(p);
			for (unsigned i = 0; i < 3; ++i)
			{
				list_nodes(f->adjacent_vertices()[i], storage);
				list_nodes(f->adjacent_edges()[i], storage);
			}
		}
		else if (p->type() == EDGE)
		{
			list_nodes(p, storage);
			list_nodes(p->adjacent_vertices()[0], storage);
			list_nodes(p->adjacent_vertices()[1], storage);
		}
		else // VERTEX
		{
			list_nodes(p, storage);
		}
	}

	void GeodesicAlgorithmSubdivisionLog::list_nodes_visible_from_node(node_pointer node, // list all nodes that belong to this mesh element
																	   std::vector<node_pointer> &storage,
																	   std::vector<double> &distances,
																	   double threshold_distance)
	{
		MeshElementBase *p = node->base_element();
		assert(p->type() != UNDEFINED_POINT);
		assert(storage.size() == distances.size());

		if (p->type() == VERTEX)
		{
			vertex_pointer v = static_cast<vertex_pointer>(p);

			for (unsigned i = 0; i < v->adjacent_edges().size(); ++i)
			{
				edge_pointer e = v->adjacent_edges()[i];
				vertex_pointer v_opposite = e->opposite_vertex(v);
				list_nodes(e, storage, threshold_distance);
				list_nodes(v_opposite, storage, threshold_distance);
			}
			for (unsigned i = 0; i < v->adjacent_faces().size(); ++i)
			{
				face_pointer f = v->adjacent_faces()[i];
				edge_pointer e = f->opposite_edge(v);
				list_nodes(e, storage, threshold_distance);
			}
		}
		else if (p->type() == EDGE)
		{
			edge_pointer e = static_cast<edge_pointer>(p);

			vertex_pointer v0 = e->adjacent_vertices()[0];
			vertex_pointer v1 = e->adjacent_vertices()[1];
			list_nodes(v0, storage, threshold_distance);
			list_nodes(v1, storage, threshold_distance);

			for (unsigned i = 0; i < e->adjacent_faces().size(); ++i)
			{
				face_pointer f = e->adjacent_faces()[i];

				list_nodes(f->next_edge(e, v0), storage, threshold_distance);
				list_nodes(f->next_edge(e, v1), storage, threshold_distance);
				list_nodes(f->opposite_vertex(e), storage, threshold_distance);
			}
		}
		else
		{
			assert(0);
		}

		unsigned index = distances.size();
		distances.resize(storage.size());
		for (; index < storage.size(); ++index)
		{
			distances[index] = get_edge_pass_face_weight(node, &storage[index]->surface_point()) * node->distance(&storage[index]->surface_point());
		}
	}

} // namespace geodesic

#endif // GEODESIC_ALGORITHM_SUBDIVISION_LOG_122806
