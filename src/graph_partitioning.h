#ifndef SRC_GRAPH_PARTITIONING_H_
#define SRC_GRAPH_PARTITIONING_H_

#include "graph.h"
#include "graph_connectivity.h"
#include <vector>
#include <climits>

namespace GraphPartitioning {

/**
 * @brief Partitioning methods for bisectional bandwidth computation.
 */
enum class PartitioningMethod {
    EXACT_BRUTE_FORCE,    // Try all possible partitions (small graphs only)
    EXACT_BRANCH_BOUND,   // Branch and bound optimization (future)
    SPECTRAL,             // Spectral partitioning using Fiedler vector (future)
    KERNIGHAN_LIN,        // Kernighan-Lin heuristic (future)
    KNOWN_TOPOLOGY,       // Use known formulas for special graphs (future)
    AUTO                  // Choose best method based on graph size
};

/**
 * @brief Known graph topologies for optimized computation.
 */
enum class GraphTopology {
    MESH_2D,
    MESH_3D,
    TORUS_2D,
    TORUS_3D,
    HYPERCUBE,
    COMPLETE_GRAPH,
    TREE,
    UNKNOWN
};

/**
 * @brief Represents a graph partition with associated metadata.
 */
struct GraphPartition {
    std::vector<int> partition1;              // Vertices in first partition
    std::vector<int> partition2;              // Vertices in second partition
    std::vector<std::pair<int, int>> cut_edges; // Edges between partitions
    int bandwidth;                            // Size of the cut (number of edges)
    
    /**
     * @brief Checks if the partition is valid (balanced).
     * @return true if partitions are balanced (differ by at most 1 vertex)
     */
    bool isBalanced() const {
        int size_diff = std::abs(static_cast<int>(partition1.size()) - 
                                static_cast<int>(partition2.size()));
        return size_diff <= 1;
    }
};

/**
 * @brief Computes the exact bisectional bandwidth of a graph.
 * 
 * Finds the minimum number of edges that must be removed to partition
 * the graph into two equal (or nearly equal) halves. Uses exact algorithms
 * suitable for small graphs only.
 * 
 * @param graph Input graph to analyze
 * @return int The bisectional bandwidth (minimum cut size)
 * @throws std::invalid_argument if graph is empty
 * @throws std::runtime_error if graph is too large for exact computation
 */
int computeExactBisectionalBandwidth(const UnidirectionalGraph& graph);

/**
 * @brief Computes an approximate bisectional bandwidth using heuristics.
 * 
 * Uses various approximation methods to find a good estimate of the
 * bisectional bandwidth. Suitable for larger graphs.
 * 
 * @param graph Input graph to analyze
 * @param method Partitioning method (defaults to AUTO)
 * @return int Approximate bisectional bandwidth
 * @throws std::invalid_argument if graph is empty
 * @throws std::runtime_error if method is not yet implemented
 */
int computeApproximateBisectionalBandwidth(const UnidirectionalGraph& graph, 
                                          PartitioningMethod method = PartitioningMethod::AUTO);

/**
 * @brief Computes bisectional bandwidth using automatic method selection.
 * 
 * Automatically chooses the best algorithm based on graph size and structure.
 * 
 * @param graph Input graph to analyze
 * @return int Bisectional bandwidth
 * @throws std::invalid_argument if graph is empty
 */
int computeBisectionalBandwidth(const UnidirectionalGraph& graph);

/**
 * @brief Finds the actual partition that achieves the bisectional bandwidth.
 * 
 * @param graph Input graph
 * @param method Partitioning method to use (defaults to AUTO)
 * @return GraphPartition Object containing the two partitions and cut edges
 * @throws std::invalid_argument if graph is empty
 * @throws std::runtime_error if method is not yet implemented
 */
GraphPartition findBisectionalPartition(const UnidirectionalGraph& graph,
                                       PartitioningMethod method = PartitioningMethod::AUTO);

/**
 * @brief Computes bisectional bandwidth for known graph topologies.
 * 
 * Provides O(1) or O(V) solutions for well-known graph types.
 * Currently supports basic topologies with plans for expansion.
 * 
 * @param graph Input graph
 * @param topology Known topology type
 * @return int Bisectional bandwidth for the topology
 * @throws std::invalid_argument if graph is empty or topology unsupported
 */
int computeKnownTopologyBandwidth(const UnidirectionalGraph& graph, 
                                 GraphTopology topology);

/**
 * @brief Helper function to calculate cut size for a given partition.
 * 
 * @param graph Input graph
 * @param vertices All vertices in the graph
 * @param partition_mask Bitmask representing the partition (1 = partition1, 0 = partition2)
 * @return int Number of edges crossing the partition
 */
int calculateCutSize(const UnidirectionalGraph& graph, 
                    const std::vector<int>& vertices, 
                    int partition_mask);

/**
 * @brief Validates that a graph is suitable for partitioning.
 * 
 * @param graph Input graph to validate
 * @throws std::invalid_argument if graph is empty or has issues
 */
void validateGraphForPartitioning(const UnidirectionalGraph& graph);

}  // namespace GraphPartitioning

#endif  // SRC_GRAPH_PARTITIONING_H_
