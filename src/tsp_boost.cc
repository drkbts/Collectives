#include "tsp_boost.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/metric_tsp_approx.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>

namespace TSPBoost {

BoostGraph convertToBoostGraph(const UnidirectionalGraph& graph, 
                              const std::vector<std::vector<double>>& weights) {
    std::vector<int> vertices = graph.getVertices();
    
    // Create vertex mapping
    std::map<int, BoostVertex> vertex_map;
    BoostGraph boost_graph(vertices.size());
    
    for (size_t i = 0; i < vertices.size(); ++i) {
        vertex_map[vertices[i]] = boost::vertex(i, boost_graph);
    }
    
    // Add edges
    WeightMap weight_map = boost::get(boost::edge_weight, boost_graph);
    
    for (int from : vertices) {
        std::vector<int> neighbors = graph.getNeighbors(from);
        for (int to : neighbors) {
            if (vertex_map.find(to) != vertex_map.end()) {
                BoostVertex v1 = vertex_map[from];
                BoostVertex v2 = vertex_map[to];
                
                // Check if edge already exists (for undirected graph)
                auto edge_pair = boost::edge(v1, v2, boost_graph);
                if (!edge_pair.second) {
                    auto [edge, inserted] = boost::add_edge(v1, v2, boost_graph);
                    
                    // Set weight
                    double weight = 1.0;  // Default weight
                    if (!weights.empty() && 
                        static_cast<size_t>(from) < weights.size() && 
                        static_cast<size_t>(to) < weights[from].size()) {
                        weight = weights[from][to];
                    }
                    weight_map[edge] = weight;
                }
            }
        }
    }
    
    return boost_graph;
}

UnidirectionalGraph convertFromBoostGraph(const BoostGraph& boost_graph) {
    UnidirectionalGraph graph;
    
    // Add vertices
    auto vertex_range = boost::vertices(boost_graph);
    for (auto it = vertex_range.first; it != vertex_range.second; ++it) {
        graph.addVertex(static_cast<int>(*it));
    }
    
    // Add edges
    auto edge_range = boost::edges(boost_graph);
    for (auto it = edge_range.first; it != edge_range.second; ++it) {
        int from = static_cast<int>(boost::source(*it, boost_graph));
        int to = static_cast<int>(boost::target(*it, boost_graph));
        graph.addEdge(from, to);
    }
    
    return graph;
}

// TSP Visitor class that properly implements the TSP visitor concept
class TSPVisitor {
public:
    TSPVisitor(std::vector<BoostVertex>& tour_) : tour(tour_) {}
    
    template<class Vertex, class Graph>
    void visit_vertex(Vertex v, const Graph& g) {
        // Only add vertex if not already visited
        if (visited.find(v) == visited.end()) {
            tour.push_back(v);
            visited.insert(v);
        }
    }
    
private:
    std::vector<BoostVertex>& tour;
    std::set<BoostVertex> visited;
};

std::vector<int> findApproximateHamiltonianCycle(const UnidirectionalGraph& graph,
                                                const std::vector<std::vector<double>>& weights) {
    std::vector<int> vertices = graph.getVertices();
    
    if (vertices.size() < 3) {
        return {};  // Need at least 3 vertices for a cycle
    }
    
    // Convert to complete graph for TSP
    BoostGraph complete_graph(vertices.size());
    WeightMap weight_map = boost::get(boost::edge_weight, complete_graph);
    
    // Create mapping from original vertex IDs to boost vertex descriptors
    std::map<int, BoostVertex> vertex_map;
    for (size_t i = 0; i < vertices.size(); ++i) {
        vertex_map[vertices[i]] = boost::vertex(i, complete_graph);
    }
    
    // Add all possible edges to make it complete
    for (size_t i = 0; i < vertices.size(); ++i) {
        for (size_t j = i + 1; j < vertices.size(); ++j) {
            BoostVertex v1 = boost::vertex(i, complete_graph);
            BoostVertex v2 = boost::vertex(j, complete_graph);
            auto [edge, inserted] = boost::add_edge(v1, v2, complete_graph);
            
            // Determine weight
            double weight = 1.0;  // Default weight
            int orig_i = vertices[i];
            int orig_j = vertices[j];
            
            if (!weights.empty() && 
                static_cast<size_t>(orig_i) < weights.size() && 
                static_cast<size_t>(orig_j) < weights[orig_i].size()) {
                weight = weights[orig_i][orig_j];
            } else {
                // Use shortest path distance if no direct weight
                // For simplicity, use 1.0 if edge exists in original graph, 2.0 otherwise
                if (graph.hasEdge(orig_i, orig_j) || graph.hasEdge(orig_j, orig_i)) {
                    weight = 1.0;
                } else {
                    weight = 2.0;  // Penalty for non-existing edges
                }
            }
            
            weight_map[edge] = weight;
        }
    }
    
    // Run TSP approximation
    std::vector<BoostVertex> tour;
    try {
        TSPVisitor visitor(tour);
        boost::metric_tsp_approx(complete_graph, visitor);
    } catch (const std::exception& e) {
        return {};  // Failed to find tour
    }
    
    // Convert back to original vertex IDs
    std::vector<int> result;
    for (BoostVertex v : tour) {
        result.push_back(vertices[static_cast<size_t>(v)]);
    }
    
    // Add the starting vertex at the end to complete the cycle
    if (!result.empty()) {
        result.push_back(result[0]);
    }
    
    return result;
}

UnidirectionalGraph computeMSTKruskal(const UnidirectionalGraph& graph,
                                     const std::vector<std::vector<double>>& weights) {
    BoostGraph boost_graph = convertToBoostGraph(graph, weights);
    std::vector<int> vertices = graph.getVertices();
    
    // Run Kruskal's algorithm
    std::vector<BoostEdge> mst_edges;
    boost::kruskal_minimum_spanning_tree(boost_graph, std::back_inserter(mst_edges));
    
    // Convert result back to UnidirectionalGraph
    UnidirectionalGraph mst;
    
    // Add all vertices
    for (int v : vertices) {
        mst.addVertex(v);
    }
    
    // Add MST edges
    for (const BoostEdge& edge : mst_edges) {
        int from = vertices[static_cast<size_t>(boost::source(edge, boost_graph))];
        int to = vertices[static_cast<size_t>(boost::target(edge, boost_graph))];
        mst.addEdge(from, to);
    }
    
    return mst;
}

UnidirectionalGraph computeMSTPrim(const UnidirectionalGraph& graph,
                                  const std::vector<std::vector<double>>& weights) {
    BoostGraph boost_graph = convertToBoostGraph(graph, weights);
    std::vector<int> vertices = graph.getVertices();
    
    if (vertices.empty()) {
        return UnidirectionalGraph();
    }
    
    // Run Prim's algorithm
    std::vector<BoostVertex> predecessors(boost::num_vertices(boost_graph));
    boost::prim_minimum_spanning_tree(boost_graph, &predecessors[0]);
    
    // Convert result back to UnidirectionalGraph
    UnidirectionalGraph mst;
    
    // Add all vertices
    for (int v : vertices) {
        mst.addVertex(v);
    }
    
    // Add MST edges
    for (size_t i = 0; i < predecessors.size(); ++i) {
        if (predecessors[i] != i) {  // Not the root
            int from = vertices[static_cast<size_t>(predecessors[i])];
            int to = vertices[i];
            mst.addEdge(from, to);
        }
    }
    
    return mst;
}

}  // namespace TSPBoost
