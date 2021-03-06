#ifndef GRAPHITE_VG_PATHTRACEVISITOR_H
#define GRAPHITE_VG_PATHTRACEVISITOR_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/range/irange.hpp>
#include <boost/pending/indirect_cmp.hpp>

#include <vector>

namespace graphite
{
namespace vg
{
	class PathTraceVisitor : public boost:: default_dfs_visitor
	{
	public:
	    PathTraceVisitor()
		{
		}

		~PathTraceVisitor()
		{
		}

		template< typename V, typename Graph >
			void discover_vertex(V u, Graph& g)
		{
			this->m_current_path.emplace_back(u);
			if (boost::out_degree(u, g) == 0)
			{
				this->m_paths.emplace_back(this->m_current_path);
				this->m_current_path.pop_back();
				std::cout << "end" << std::endl;
			}
			else
			{
				std::cout << "discover" << std::endl;
			}
			std::cout << "paths: " << m_paths.size() << std::endl;
		}

		template< typename V, typename Graph >
			void finish_vertex(V u, Graph& g)
		{
			/*
			std::cout << "finish1" << std::endl;
			// if we are at the end of a path
			if (boost::out_degree(u, g) == 0)
			{
				std::cout << "finish2a: " << this->m_current_path.size() << std::endl;
				this->m_paths.emplace_back(this->m_current_path);
				std::cout << "finish2b: " << this->m_paths.size() << std::endl;
			}
			std::cout << "finish3" << std::endl;
			this->m_current_path.pop_back();
			*/
		}

		std::vector< std::vector< VariantGraph::VariantVertexDescriptor > > getPaths()
		{
			std::cout << "getPaths: " << this->m_paths.size() << std::endl;
			return this->m_paths;
		}

	private:
		std::vector< VariantGraph::VariantVertexDescriptor > m_current_path;
		std::vector< std::vector< VariantGraph::VariantVertexDescriptor > > m_paths;
	};
}
}

#endif //GRAPHITE_VG_PATHTRACEVISITOR_H
