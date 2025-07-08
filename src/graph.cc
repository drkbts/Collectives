#include "graph.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <unordered_set>
#include <vector>

UnidirectionalGraph::UnidirectionalGraph() : edgeCount_(0) {}

void UnidirectionalGraph::addVertex(int v) {
  // Insert an empty set if the key does not exist.
  adjList_.emplace(v, std::unordered_set<int>());
}

void UnidirectionalGraph::addEdge(int from, int to) {
  // Ensure both vertices exist in the graph before adding the edge.
  addVertex(from);
  addVertex(to);

  // insert() returns a pair, with .second being a bool that is true if
  // insertion took place. This provides an efficient way to check for
  // duplicates and add the edge in one operation.
  if (adjList_.at(from).insert(to).second) {
    edgeCount_++;
  }
}

bool UnidirectionalGraph::hasVertex(int v) const {
  return adjList_.count(v) > 0;
}

bool UnidirectionalGraph::hasEdge(int from, int to) const {
  if (!hasVertex(from)) {
    return false;
  }
  // With an unordered_set, checking for edge existence is O(1) on average.
  const auto& neighbors = adjList_.at(from);
  return neighbors.count(to) > 0;
}

size_t UnidirectionalGraph::getVertexCount() const { return adjList_.size(); }

size_t UnidirectionalGraph::getEdgeCount() const { return edgeCount_; }

std::vector<int> UnidirectionalGraph::getNeighbors(int v) const {
  try {
    // Return a copy of the neighbors for API consistency.
    const auto& neighbors_set = adjList_.at(v);
    return std::vector<int>(neighbors_set.begin(), neighbors_set.end());
  } catch (const std::out_of_range& oor) {
    // Provide a more informative error message.
    throw std::out_of_range("Vertex not found in graph.");
  }
}

std::vector<int> UnidirectionalGraph::getVertices() const {
  std::vector<int> vertices;
  vertices.reserve(adjList_.size());
  
  for (const auto& pair : adjList_) {
    vertices.push_back(pair.first);
  }
  
  return vertices;
}

void UnidirectionalGraph::printGraph() const {
  for (const auto& pair : adjList_) {
    std::cout << pair.first << " -> { ";
    for (int neighbor : pair.second) {
      std::cout << neighbor << " ";
    }
    std::cout << "}\n";
  }
}