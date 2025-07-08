#ifndef SRC_UNIDIRECTIONAL_GRAPH_H_
#define SRC_UNIDIRECTIONAL_GRAPH_H_

#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/**
 * @brief Represents a uni-directional (directed) graph using an adjacency list.
 *
 * This class allows for the creation and manipulation of a directed graph where
 * edges have a specific direction from a source vertex to a destination vertex.
 */
class UnidirectionalGraph {
 public:
  /**
   * @brief Constructs an empty graph.
   */
  UnidirectionalGraph();

  /**
   * @brief Adds a vertex to the graph. If it already exists, no action is taken.
   * @param v The vertex to add.
   */
  void addVertex(int v);

  /**
   * @brief Adds a directed edge from vertex 'from' to vertex 'to'.
   * If the vertices do not exist, they are automatically created.
   * @param from The source vertex.
   * @param to The destination vertex.
   */
  void addEdge(int from, int to);

  bool hasVertex(int v) const;
  bool hasEdge(int from, int to) const;
  size_t getVertexCount() const;
  size_t getEdgeCount() const;

  /**
   * @brief Gets the neighbors of a given vertex.
   * @param v The vertex whose neighbors to retrieve.
   * @return A vector containing the neighbors.
   * @throws std::out_of_range if the vertex does not exist.
   */
  std::vector<int> getNeighbors(int v) const;

  /**
   * @brief Gets all vertices in the graph.
   * @return A vector containing all vertices in the graph.
   */
  std::vector<int> getVertices() const;

  void printGraph() const;

 private:
  // Using unordered_set for neighbors provides O(1) average time for edge
  // lookups and insertions.
  std::unordered_map<int, std::unordered_set<int>> adjList_;
  size_t edgeCount_;
};

#endif  // SRC_UNIDIRECTIONAL_GRAPH_H_