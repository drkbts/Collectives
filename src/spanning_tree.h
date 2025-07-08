#ifndef SRC_SPANNING_TREE_H_
#define SRC_SPANNING_TREE_H_

#include "graph.h"
#include "graph_connectivity.h"
#include <vector>

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

}  // namespace SpanningTree

#endif  // SRC_SPANNING_TREE_H_
