#ifndef SRC_GRAPH_CONNECTIVITY_H_
#define SRC_GRAPH_CONNECTIVITY_H_

#include "graph.h"

namespace GraphConnectivity {

/**
 * @brief Checks if the given graph is connected.
 * 
 * For directed graphs, this checks if the graph is weakly connected
 * (i.e., connected when treating edges as undirected).
 * 
 * @param graph The input graph to check connectivity
 * @return true if the graph is connected, false otherwise
 */
bool isConnected(const UnidirectionalGraph& graph);

/**
 * @brief Gets all vertices in the graph.
 * 
 * Helper function to extract all vertices from the graph.
 * 
 * @param graph The input graph
 * @return std::vector<int> containing all vertices in the graph
 */
std::vector<int> getAllVertices(const UnidirectionalGraph& graph);

}  // namespace GraphConnectivity

#endif  // SRC_GRAPH_CONNECTIVITY_H_
