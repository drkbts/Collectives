#ifndef SRC_HAMILTONIAN_H_
#define SRC_HAMILTONIAN_H_

#include "graph.h"
#include <vector>
#include <unordered_set>

namespace Hamiltonian {

/**
 * @brief Finds a Hamiltonian cycle in the graph if one exists.
 * 
 * A Hamiltonian cycle is a cycle that visits each vertex exactly once and 
 * returns to the starting vertex. This is an NP-complete problem, so the
 * algorithm uses backtracking which may be slow for large graphs.
 * 
 * For directed graphs, this function respects edge directions when finding
 * the cycle.
 * 
 * @param graph The input graph to search for a Hamiltonian cycle
 * @return std::vector<int> containing the vertices in the Hamiltonian cycle
 *         (including the starting vertex at the end to complete the cycle),
 *         or an empty vector if no Hamiltonian cycle exists
 */
std::vector<int> findExactHamiltonianCycle(const UnidirectionalGraph& graph);

/**
 * @brief Finds an approximate Hamiltonian cycle using Boost Graph Library's TSP approximation.
 * 
 * Uses Boost's metric_tsp_approx algorithm which provides a 2-approximation
 * for the traveling salesman problem. This is much faster than the exact
 * algorithm but only provides an approximation.
 * 
 * The algorithm works by:
 * 1. Converting the graph to a complete graph with appropriate weights
 * 2. Running Boost's TSP approximation algorithm
 * 3. Converting the result back to our graph format
 * 
 * @param graph The input graph to search for an approximate Hamiltonian cycle
 * @param weights Optional edge weights (default: unit weights for existing edges,
 *                penalty weights for non-existing edges)
 * @return std::vector<int> containing the vertices in the approximate Hamiltonian cycle
 *         (including the starting vertex at the end to complete the cycle),
 *         or an empty vector if no solution found
 */
std::vector<int> findApproximateHamiltonianCycle(const UnidirectionalGraph& graph,
                                                const std::vector<std::vector<double>>& weights = {});

/**
 * @brief Gets all vertices in the graph.
 * 
 * Helper function to extract all vertices from the graph.
 * 
 * @param graph The input graph
 * @return std::vector<int> containing all vertices in the graph
 */
std::vector<int> getAllVertices(const UnidirectionalGraph& graph);

}  // namespace Hamiltonian

#endif  // SRC_HAMILTONIAN_H_
