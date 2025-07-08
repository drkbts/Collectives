#ifndef SRC_SPANNING_TREE_H_
#define SRC_SPANNING_TREE_H_

#include "graph.h"
#include <vector>
#include <unordered_set>

namespace SpanningTree {

/**
 * @brief Computes a spanning tree of the given graph using Depth-First Search.
 * 
 * A spanning tree is a subgraph that:
 * - Is connected and acyclic (forms a tree)
 * - Includes all vertices of the original graph
 * - Has exactly |V| - 1 edges where |V| is the number of vertices
 * 
 * For directed graphs, this function treats edges as undirected when building
 * the spanning tree. The resulting tree will have directed edges, but the
 * connectivity is based on the undirected version of the input graph.
 * 
 * @param graph The input graph to compute spanning tree for
 * @return UnidirectionalGraph representing the spanning tree
 * @throws std::invalid_argument if the graph is empty
 * @throws std::runtime_error if the graph is not connected (no spanning tree exists)
 */
UnidirectionalGraph computeSpanningTree(const UnidirectionalGraph& graph);

/**
 * @brief Computes a spanning forest of the given graph using Depth-First Search.
 * 
 * A spanning forest is a collection of spanning trees, one for each connected
 * component of the graph. If the graph is connected, this returns the same
 * result as computeSpanningTree.
 * 
 * @param graph The input graph to compute spanning forest for
 * @return UnidirectionalGraph representing the spanning forest
 * @throws std::invalid_argument if the graph is empty
 */
UnidirectionalGraph computeSpanningForest(const UnidirectionalGraph& graph);

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
 * Helper function to extract all vertices from the graph since the
 * UnidirectionalGraph class doesn't provide direct access to the vertex set.
 * 
 * @param graph The input graph
 * @return std::vector<int> containing all vertices in the graph
 */
std::vector<int> getAllVertices(const UnidirectionalGraph& graph);

}  // namespace SpanningTree

#endif  // SRC_SPANNING_TREE_H_
