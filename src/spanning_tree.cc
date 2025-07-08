#include "spanning_tree.h"
#include <stdexcept>
#include <stack>
#include <unordered_set>

namespace SpanningTree {

UnidirectionalGraph computeSpanningTree(const UnidirectionalGraph& graph) {
    std::vector<int> vertices = GraphConnectivity::getAllVertices(graph);
    
    if (vertices.empty()) {
        throw std::invalid_argument("Cannot compute spanning tree of empty graph");
    }
    
    if (!GraphConnectivity::isConnected(graph)) {
        throw std::runtime_error("Graph is not connected - no spanning tree exists");
    }
    
    UnidirectionalGraph spanningTree;
    
    // Add all vertices to the spanning tree
    for (int vertex : vertices) {
        spanningTree.addVertex(vertex);
    }
    
    if (vertices.size() == 1) {
        return spanningTree;  // Single vertex tree
    }
    
    // Use DFS to build spanning tree
    std::unordered_set<int> visited;
    std::stack<int> stack;
    
    // Start from the first vertex
    stack.push(vertices[0]);
    visited.insert(vertices[0]);
    
    while (!stack.empty() && visited.size() < vertices.size()) {
        int current = stack.top();
        stack.pop();
        
        // Find an unvisited neighbor
        for (int neighbor : vertices) {
            if (visited.find(neighbor) == visited.end()) {
                // Check if there's an edge in either direction
                if (graph.hasEdge(current, neighbor)) {
                    // Add edge from current to neighbor
                    spanningTree.addEdge(current, neighbor);
                    visited.insert(neighbor);
                    stack.push(neighbor);
                    break;  // Only add one edge per iteration
                } else if (graph.hasEdge(neighbor, current)) {
                    // Add edge from neighbor to current
                    spanningTree.addEdge(neighbor, current);
                    visited.insert(neighbor);
                    stack.push(neighbor);
                    break;  // Only add one edge per iteration
                }
            }
        }
    }
    
    return spanningTree;
}

UnidirectionalGraph computeSpanningForest(const UnidirectionalGraph& graph) {
    std::vector<int> vertices = GraphConnectivity::getAllVertices(graph);
    
    if (vertices.empty()) {
        throw std::invalid_argument("Cannot compute spanning forest of empty graph");
    }
    
    UnidirectionalGraph spanningForest;
    
    // Add all vertices to the spanning forest
    for (int vertex : vertices) {
        spanningForest.addVertex(vertex);
    }
    
    if (vertices.size() == 1) {
        return spanningForest;  // Single vertex forest
    }
    
    std::unordered_set<int> globalVisited;
    
    // Process each connected component
    for (int startVertex : vertices) {
        if (globalVisited.find(startVertex) != globalVisited.end()) {
            continue;  // Already processed this component
        }
        
        // DFS for this component
        std::unordered_set<int> componentVisited;
        std::stack<int> stack;
        
        stack.push(startVertex);
        componentVisited.insert(startVertex);
        globalVisited.insert(startVertex);
        
        while (!stack.empty()) {
            int current = stack.top();
            stack.pop();
            
            // Find an unvisited neighbor in this component
            for (int neighbor : vertices) {
                if (componentVisited.find(neighbor) == componentVisited.end()) {
                    // Check if there's an edge in either direction
                    if (graph.hasEdge(current, neighbor)) {
                        // Add edge from current to neighbor
                        spanningForest.addEdge(current, neighbor);
                        componentVisited.insert(neighbor);
                        globalVisited.insert(neighbor);
                        stack.push(neighbor);
                        break;  // Only add one edge per iteration
                    } else if (graph.hasEdge(neighbor, current)) {
                        // Add edge from neighbor to current
                        spanningForest.addEdge(neighbor, current);
                        componentVisited.insert(neighbor);
                        globalVisited.insert(neighbor);
                        stack.push(neighbor);
                        break;  // Only add one edge per iteration
                    }
                }
            }
        }
    }
    
    return spanningForest;
}

}  // namespace SpanningTree
