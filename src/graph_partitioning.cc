#include "graph_partitioning.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <bitset>

namespace GraphPartitioning {

void validateGraphForPartitioning(const UnidirectionalGraph& graph) {
    std::vector<int> vertices = GraphConnectivity::getAllVertices(graph);
    if (vertices.empty()) {
        throw std::invalid_argument("Cannot compute bisectional bandwidth of empty graph");
    }
}

int calculateCutSize(const UnidirectionalGraph& graph, 
                    const std::vector<int>& vertices, 
                    int partition_mask) {
    int cut_size = 0;
    int n = vertices.size();
    
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            int vertex_i = vertices[i];
            int vertex_j = vertices[j];
            
            // Check if vertices are in different partitions
            bool i_in_partition1 = (partition_mask & (1 << i)) != 0;
            bool j_in_partition1 = (partition_mask & (1 << j)) != 0;
            
            if (i_in_partition1 != j_in_partition1) {
                // Count edges between different partitions
                if (graph.hasEdge(vertex_i, vertex_j)) {
                    cut_size++;
                }
                if (graph.hasEdge(vertex_j, vertex_i)) {
                    cut_size++;
                }
            }
        }
    }
    
    return cut_size;
}

GraphPartition createPartitionFromMask(const UnidirectionalGraph& graph,
                                      const std::vector<int>& vertices,
                                      int partition_mask) {
    GraphPartition partition;
    int n = vertices.size();
    
    // Separate vertices into two partitions
    for (int i = 0; i < n; i++) {
        if (partition_mask & (1 << i)) {
            partition.partition1.push_back(vertices[i]);
        } else {
            partition.partition2.push_back(vertices[i]);
        }
    }
    
    // Find cut edges
    for (int v1 : partition.partition1) {
        for (int v2 : partition.partition2) {
            if (graph.hasEdge(v1, v2)) {
                partition.cut_edges.push_back({v1, v2});
            }
            if (graph.hasEdge(v2, v1)) {
                partition.cut_edges.push_back({v2, v1});
            }
        }
    }
    
    partition.bandwidth = partition.cut_edges.size();
    return partition;
}

int computeExactBisectionalBandwidth(const UnidirectionalGraph& graph) {
    validateGraphForPartitioning(graph);
    
    std::vector<int> vertices = GraphConnectivity::getAllVertices(graph);
    int n = vertices.size();
    
    // Handle trivial cases
    if (n <= 1) {
        return 0;
    }
    
    if (n == 2) {
        // For 2 vertices, bandwidth is the number of edges between them
        int edges = 0;
        if (graph.hasEdge(vertices[0], vertices[1])) edges++;
        if (graph.hasEdge(vertices[1], vertices[0])) edges++;
        return edges;
    }
    
    // Check if graph is too large for exact computation
    if (n > 20) {
        throw std::runtime_error("Graph too large for exact computation (>" + 
                                std::to_string(20) + " vertices). Use approximation methods.");
    }
    
    int min_cut = INT_MAX;
    int target_size = n / 2;
    
    // Try all possible partitions where one partition has exactly target_size vertices
    for (int mask = 0; mask < (1 << n); mask++) {
        if (__builtin_popcount(mask) == target_size) {
            int cut_size = calculateCutSize(graph, vertices, mask);
            min_cut = std::min(min_cut, cut_size);
        }
    }
    
    // If n is odd, also try partitions where one partition has target_size + 1 vertices
    if (n % 2 == 1) {
        for (int mask = 0; mask < (1 << n); mask++) {
            if (__builtin_popcount(mask) == target_size + 1) {
                int cut_size = calculateCutSize(graph, vertices, mask);
                min_cut = std::min(min_cut, cut_size);
            }
        }
    }
    
    return min_cut;
}

GraphPartition findBisectionalPartition(const UnidirectionalGraph& graph,
                                       PartitioningMethod method) {
    validateGraphForPartitioning(graph);
    
    std::vector<int> vertices = GraphConnectivity::getAllVertices(graph);
    int n = vertices.size();
    
    // Handle trivial cases
    if (n <= 1) {
        GraphPartition partition;
        if (n == 1) {
            partition.partition1.push_back(vertices[0]);
        }
        partition.bandwidth = 0;
        return partition;
    }
    
    if (n == 2) {
        GraphPartition partition;
        partition.partition1.push_back(vertices[0]);
        partition.partition2.push_back(vertices[1]);
        
        // Count edges between the two vertices
        if (graph.hasEdge(vertices[0], vertices[1])) {
            partition.cut_edges.push_back({vertices[0], vertices[1]});
        }
        if (graph.hasEdge(vertices[1], vertices[0])) {
            partition.cut_edges.push_back({vertices[1], vertices[0]});
        }
        
        partition.bandwidth = partition.cut_edges.size();
        return partition;
    }
    
    // For now, only exact brute force is implemented
    if (method == PartitioningMethod::AUTO) {
        method = (n <= 20) ? PartitioningMethod::EXACT_BRUTE_FORCE : PartitioningMethod::EXACT_BRUTE_FORCE;
    }
    
    if (method != PartitioningMethod::EXACT_BRUTE_FORCE) {
        throw std::runtime_error("Only EXACT_BRUTE_FORCE method is currently implemented");
    }
    
    if (n > 20) {
        throw std::runtime_error("Graph too large for exact computation");
    }
    
    GraphPartition best_partition;
    int min_cut = INT_MAX;
    int target_size = n / 2;
    
    // Try all possible partitions
    for (int mask = 0; mask < (1 << n); mask++) {
        int partition_size = __builtin_popcount(mask);
        if (partition_size == target_size || (n % 2 == 1 && partition_size == target_size + 1)) {
            GraphPartition current_partition = createPartitionFromMask(graph, vertices, mask);
            
            if (current_partition.bandwidth < min_cut) {
                min_cut = current_partition.bandwidth;
                best_partition = current_partition;
            }
        }
    }
    
    return best_partition;
}

int computeApproximateBisectionalBandwidth(const UnidirectionalGraph& graph, 
                                          PartitioningMethod method) {
    validateGraphForPartitioning(graph);
    
    if (method == PartitioningMethod::AUTO) {
        // For now, fall back to exact method for small graphs
        std::vector<int> vertices = GraphConnectivity::getAllVertices(graph);
        if (vertices.size() <= 20) {
            return computeExactBisectionalBandwidth(graph);
        } else {
            throw std::runtime_error("Approximation methods not yet implemented for large graphs");
        }
    }
    
    throw std::runtime_error("Approximation methods not yet implemented");
}

int computeBisectionalBandwidth(const UnidirectionalGraph& graph) {
    validateGraphForPartitioning(graph);
    
    std::vector<int> vertices = GraphConnectivity::getAllVertices(graph);
    int n = vertices.size();
    
    // Choose method based on graph size
    if (n <= 20) {
        return computeExactBisectionalBandwidth(graph);
    } else {
        // For now, throw error for large graphs until approximation is implemented
        throw std::runtime_error("Large graphs not yet supported. Approximation methods coming soon.");
    }
}

int computeKnownTopologyBandwidth(const UnidirectionalGraph& graph, 
                                 GraphTopology topology) {
    validateGraphForPartitioning(graph);
    
    std::vector<int> vertices = GraphConnectivity::getAllVertices(graph);
    int n = vertices.size();
    
    switch (topology) {
        case GraphTopology::COMPLETE_GRAPH:
            // For complete graph K_n, bisectional bandwidth is n²/4 (approximately)
            return (n * n) / 4;
            
        case GraphTopology::TREE:
            // For any tree, bisectional bandwidth is 1 (single edge removal disconnects)
            return (n > 1) ? 1 : 0;
            
        default:
            throw std::invalid_argument("Topology not yet supported");
    }
}

}  // namespace GraphPartitioning
