#ifndef SRC_GRAPH_UTILS_H_
#define SRC_GRAPH_UTILS_H_

#include "graph.h"
#include <string>

namespace GraphUtils {

/**
 * @brief Creates a unidirectional ring of size N.
 * 
 * A unidirectional ring connects vertices 0, 1, 2, ..., N-1 in a cycle where
 * each vertex has exactly one outgoing edge to the next vertex:
 * 0 -> 1 -> 2 -> ... -> (N-1) -> 0
 * 
 * @param n The number of vertices in the ring (must be > 1)
 * @return UnidirectionalGraph containing the ring structure
 * @throws std::invalid_argument if n <= 1
 */
UnidirectionalGraph createRing(int n);

/**
 * @brief Creates a bi-directional ring of size N.
 * 
 * A bi-directional ring connects vertices 0, 1, 2, ..., N-1 in a cycle where
 * each vertex has exactly two outgoing edges (forward and backward):
 * 0 <-> 1 <-> 2 <-> ... <-> (N-1) <-> 0
 * 
 * @param n The number of vertices in the ring (must be > 1)
 * @return UnidirectionalGraph containing the bi-directional ring structure
 * @throws std::invalid_argument if n <= 1
 */
UnidirectionalGraph createBidirectionalRing(int n);

/**
 * @brief Creates the tensor product of two graphs.
 * 
 * The tensor product G ⊗ H of two graphs G and H is defined as:
 * - Vertices: For each vertex u in G and vertex v in H, there's a vertex (u,v) in G ⊗ H
 * - Edges: There's an edge between (u1,v1) and (u2,v2) in G ⊗ H if and only if
 *   there's an edge between u1 and u2 in G AND there's an edge between v1 and v2 in H
 * 
 * The vertices in the result graph are numbered as: u * |V_H| + v
 * where |V_H| is the number of vertices in graph H.
 * 
 * @param g1 The first input graph
 * @param g2 The second input graph
 * @return UnidirectionalGraph containing the tensor product G1 ⊗ G2
 * @throws std::invalid_argument if either graph is empty
 */
UnidirectionalGraph createTensorProduct(const UnidirectionalGraph& g1, const UnidirectionalGraph& g2);

/**
 * @brief Creates an N-dimensional torus graph.
 * 
 * A torus graph is the tensor product of N bidirectional rings, where each ring
 * has size D_i. The resulting graph represents an N-dimensional torus where
 * each dimension has size D_i.
 * 
 * For example, createTorus({3, 4}) creates a 2D torus (3×4 grid with wraparound)
 * by computing the tensor product of a 3-vertex bidirectional ring and a 
 * 4-vertex bidirectional ring.
 * 
 * @param dimensions A vector of positive integers representing the size of each dimension
 * @return UnidirectionalGraph containing the N-dimensional torus
 * @throws std::invalid_argument if dimensions is empty or any dimension <= 1
 */
UnidirectionalGraph createTorus(const std::vector<int>& dimensions);

/**
 * @brief Parses a graph expression string and creates the corresponding graph.
 * 
 * This function implements a simple DSL (Domain Specific Language) for graph creation.
 * Currently supported expressions:
 * - "uR[N]": Creates a unidirectional ring of size N
 * - "bR[N]": Creates a bidirectional ring of size N
 * - "E1 x E2": Creates the tensor product of expressions E1 and E2
 * - "(E)": Parentheses for grouping and resolving precedence
 * 
 * Examples:
 * - "uR[3] x bR[4]": Tensor product of 3-vertex unidirectional and 4-vertex bidirectional ring
 * - "bR[5] x bR[3]": Tensor product of two bidirectional rings
 * - "(uR[3] x bR[4]) x uR[5]": Explicit grouping for nested tensor products
 * - "uR[3] x (bR[4] x uR[5])": Alternative grouping
 * 
 * @param expression The graph expression string to parse
 * @return UnidirectionalGraph corresponding to the parsed expression
 * @throws std::invalid_argument if the expression is invalid or malformed
 * @throws std::out_of_range if numeric values are out of valid range
 */
UnidirectionalGraph parseGraphExpression(const std::string& expression);

}  // namespace GraphUtils

#endif  // SRC_GRAPH_UTILS_H_
