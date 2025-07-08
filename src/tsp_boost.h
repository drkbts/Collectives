#ifndef SRC_TSP_BOOST_H_
#define SRC_TSP_BOOST_H_

#include "graph.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/metric_tsp_approx.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <vector>

namespace TSPBoost {

// Define Boost Graph types
using BoostGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, 
                                        boost::no_property, boost::property<boost::edge_weight_t, double>>;
using BoostVertex = boost::graph_traits<BoostGraph>::vertex_descriptor;
using BoostEdge = boost::graph_traits<BoostGraph>::edge_descriptor;
using WeightMap = boost::property_map<BoostGraph, boost::edge_weight_t>::type;

/**
 * @brief Converts UnidirectionalGraph to Boost Graph format.
 * 
 * Creates an undirected Boost graph from our UnidirectionalGraph.
 * For directed graphs, treats edges as undirected for algorithms like TSP.
 * Assigns unit weights to all edges unless weights are provided.
 * 
 * @param graph Input UnidirectionalGraph
 * @param weights Optional edge weights (default: all edges have weight 1.0)
 * @return BoostGraph representation
 */
BoostGraph convertToBoostGraph(const UnidirectionalGraph& graph, 
                              const std::vector<std::vector<double>>& weights = {});

/**
 * @brief Converts Boost Graph back to UnidirectionalGraph.
 * 
 * @param boost_graph Input Boost graph
 * @return UnidirectionalGraph representation
 */
UnidirectionalGraph convertFromBoostGraph(const BoostGraph& boost_graph);

/**
 * @brief Finds an approximate Hamiltonian cycle using Boost's TSP approximation.
 * 
 * Uses Boost's metric_tsp_approx algorithm which provides a 2-approximation
 * for the traveling salesman problem on complete graphs with triangle inequality.
 * 
 * @param graph Input graph
 * @param weights Edge weights (if empty, uses unit weights)
 * @return Vector of vertices forming the approximate Hamiltonian cycle,
 *         or empty vector if no solution found
 */
std::vector<int> findApproximateHamiltonianCycle(const UnidirectionalGraph& graph,
                                                const std::vector<std::vector<double>>& weights = {});

/**
 * @brief Computes minimum spanning tree using Boost's Kruskal algorithm.
 * 
 * @param graph Input graph
 * @param weights Edge weights (if empty, uses unit weights)
 * @return UnidirectionalGraph representing the MST
 */
UnidirectionalGraph computeMSTKruskal(const UnidirectionalGraph& graph,
                                     const std::vector<std::vector<double>>& weights = {});

/**
 * @brief Computes minimum spanning tree using Boost's Prim algorithm.
 * 
 * @param graph Input graph
 * @param weights Edge weights (if empty, uses unit weights)
 * @return UnidirectionalGraph representing the MST
 */
UnidirectionalGraph computeMSTPrim(const UnidirectionalGraph& graph,
                                  const std::vector<std::vector<double>>& weights = {});

}  // namespace TSPBoost

#endif  // SRC_TSP_BOOST_H_
