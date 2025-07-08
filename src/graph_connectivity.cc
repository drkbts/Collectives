#include "graph_connectivity.h"
#include <unordered_set>
#include <stack>

namespace GraphConnectivity {

std::vector<int> getAllVertices(const UnidirectionalGraph& graph) {
    return graph.getVertices();
}

bool isConnected(const UnidirectionalGraph& graph) {
    std::vector<int> vertices = getAllVertices(graph);
    
    if (vertices.empty()) {
        return true;  // Empty graph is considered connected
    }
    
    if (vertices.size() == 1) {
        return true;  // Single vertex is connected
    }
    
    // Use DFS to check if all vertices are reachable from the first vertex
    std::unordered_set<int> visited;
    std::stack<int> stack;
    
    stack.push(vertices[0]);
    visited.insert(vertices[0]);
    
    while (!stack.empty()) {
        int current = stack.top();
        stack.pop();
        
        // Check all neighbors (treating edges as undirected)
        for (int vertex : vertices) {
            if (visited.find(vertex) == visited.end()) {
                // Check if there's an edge in either direction
                if (graph.hasEdge(current, vertex) || graph.hasEdge(vertex, current)) {
                    visited.insert(vertex);
                    stack.push(vertex);
                }
            }
        }
    }
    
    return visited.size() == vertices.size();
}

}  // namespace GraphConnectivity
