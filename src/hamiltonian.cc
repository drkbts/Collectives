#include "hamiltonian.h"
#include "tsp_boost.h"
#include <stdexcept>
#include <algorithm>
#include <unordered_set>

namespace Hamiltonian {

std::vector<int> getAllVertices(const UnidirectionalGraph& graph) {
    return graph.getVertices();
}

// Helper function for Hamiltonian cycle backtracking
bool hamiltonianCycleHelper(const UnidirectionalGraph& graph, 
                           const std::vector<int>& vertices,
                           std::vector<int>& path, 
                           std::unordered_set<int>& visited,
                           int pos) {
    // Base case: if we've visited all vertices
    if (pos == vertices.size()) {
        // Check if there's an edge from the last vertex back to the start
        if (graph.hasEdge(path[pos - 1], path[0])) {
            path.push_back(path[0]); // Complete the cycle
            return true;
        }
        return false;
    }
    
    // Try all vertices as the next vertex in the path
    for (int vertex : vertices) {
        // Skip if already visited
        if (visited.find(vertex) != visited.end()) {
            continue;
        }
        
        // If this is not the first vertex, check if there's an edge from previous vertex
        if (pos > 0 && !graph.hasEdge(path[pos - 1], vertex)) {
            continue;
        }
        
        // Add vertex to path
        path[pos] = vertex;
        visited.insert(vertex);
        
        // Recursively solve for the rest
        if (hamiltonianCycleHelper(graph, vertices, path, visited, pos + 1)) {
            return true;
        }
        
        // Backtrack
        visited.erase(vertex);
    }
    
    return false;
}

std::vector<int> findExactHamiltonianCycle(const UnidirectionalGraph& graph) {
    std::vector<int> vertices = getAllVertices(graph);
    
    // Need at least 3 vertices for a cycle
    if (vertices.size() < 3) {
        return {};
    }
    
    // Sort vertices for consistent ordering
    std::sort(vertices.begin(), vertices.end());
    
    // Try starting from each vertex
    for (int startVertex : vertices) {
        std::vector<int> path(vertices.size());
        std::unordered_set<int> visited;
        
        // Start with the chosen vertex
        path[0] = startVertex;
        visited.insert(startVertex);
        
        // Try to find a Hamiltonian cycle starting from this vertex
        if (hamiltonianCycleHelper(graph, vertices, path, visited, 1)) {
            return path;
        }
    }
    
    // No Hamiltonian cycle found
    return {};
}

std::vector<int> findApproximateHamiltonianCycle(const UnidirectionalGraph& graph,
                                                const std::vector<std::vector<double>>& weights) {
    return TSPBoost::findApproximateHamiltonianCycle(graph, weights);
}

}  // namespace Hamiltonian
