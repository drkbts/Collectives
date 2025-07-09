#ifndef SRC_GRAPH_COLLECTIVE_H_
#define SRC_GRAPH_COLLECTIVE_H_

#include "graph.h"
#include "graph_connectivity.h"
#include "spanning_tree.h"
#include <vector>
#include <map>
#include <memory>
#include <limits>
#include <algorithm>
#include <stdexcept>

namespace GraphCollective {

/**
 * @brief Communication patterns for all-reduce operations.
 */
enum class CommunicationPattern {
    TREE,           // Binary tree pattern
    HYPERCUBE,      // Hypercube communication
    BUTTERFLY,      // Butterfly network
    RING,           // Ring topology
    MESH_2D,        // 2D mesh
    CUSTOM_GRAPH,   // Use provided graph structure
    OPTIMAL         // Choose optimal pattern based on graph
};

/**
 * @brief Represents a single communication step in all-reduce.
 */
struct CommunicationStep {
    int source_node;
    int target_node;
    int round;              // Communication round number
    std::string operation;  // Description of operation performed
    
    CommunicationStep(int src, int tgt, int r, const std::string& op)
        : source_node(src), target_node(tgt), round(r), operation(op) {}
};

/**
 * @brief Result of an all-reduce operation with performance metrics.
 */
template<typename T>
struct AllReduceResult {
    T final_value;                                    // Final reduced value
    std::vector<CommunicationStep> communication_steps; // All communication steps
    int total_rounds;                                 // Number of communication rounds
    int total_messages;                               // Total messages sent
    double efficiency_metric;                         // Performance efficiency score
    
    AllReduceResult() : total_rounds(0), total_messages(0), efficiency_metric(0.0) {}
};

/**
 * @brief Abstract base class for reduction operations.
 */
template<typename T>
class ReductionOp {
public:
    virtual ~ReductionOp() = default;
    virtual T operator()(const T& a, const T& b) const = 0;
    virtual T identity() const = 0;  // Identity element for the operation
    virtual std::string name() const = 0;
};

/**
 * @brief Sum reduction operation.
 */
template<typename T>
class SumReduction : public ReductionOp<T> {
public:
    T operator()(const T& a, const T& b) const override { return a + b; }
    T identity() const override { return T{}; }
    std::string name() const override { return "SUM"; }
};

/**
 * @brief Maximum reduction operation.
 */
template<typename T>
class MaxReduction : public ReductionOp<T> {
public:
    T operator()(const T& a, const T& b) const override { return std::max(a, b); }
    T identity() const override { return std::numeric_limits<T>::lowest(); }
    std::string name() const override { return "MAX"; }
};

/**
 * @brief Minimum reduction operation.
 */
template<typename T>
class MinReduction : public ReductionOp<T> {
public:
    T operator()(const T& a, const T& b) const override { return std::min(a, b); }
    T identity() const override { return std::numeric_limits<T>::max(); }
    std::string name() const override { return "MIN"; }
};

/**
 * @brief Average reduction operation.
 */
template<typename T>
class AverageReduction : public ReductionOp<T> {
public:
    explicit AverageReduction(int count) : count_(count) {}
    T operator()(const T& a, const T& b) const override { return a + b; }
    T identity() const override { return T{}; }
    std::string name() const override { return "AVERAGE"; }
    
    T finalize(const T& sum) const { return sum / static_cast<T>(count_); }
    
private:
    int count_;
};

/**
 * @brief Node state management for all-reduce operations.
 */
template<typename T>
class NodeState {
public:
    NodeState(int node_id, T initial_value) : node_id_(node_id), value_(initial_value) {}
    
    void updateValue(const T& new_value) { value_ = new_value; }
    T getValue() const { return value_; }
    int getNodeId() const { return node_id_; }
    
    void addNeighbor(int neighbor_id) { neighbors_.push_back(neighbor_id); }
    const std::vector<int>& getNeighbors() const { return neighbors_; }
    
private:
    int node_id_;
    T value_;
    std::vector<int> neighbors_;
};

/**
 * @brief Main class for graph-based collective operations.
 */
template<typename T>
class GraphCollective {
public:
    /**
     * @brief Performs all-reduce operation on graph vertices.
     * 
     * Simulates distributed all-reduce where each vertex contributes a value,
     * and the reduction result is available to all vertices.
     * 
     * @param graph Communication topology graph
     * @param initial_values Initial value at each vertex
     * @param reduction_op Reduction operation to perform
     * @param pattern Communication pattern to use
     * @return AllReduceResult containing communication steps and final result
     */
    static AllReduceResult<T> allReduce(const UnidirectionalGraph& graph,
                                       const std::map<int, T>& initial_values,
                                       const ReductionOp<T>& reduction_op,
                                       CommunicationPattern pattern = CommunicationPattern::OPTIMAL);
    
    /**
     * @brief Performs all-reduce using tree-based communication.
     * 
     * Uses spanning tree for reduce phase, then broadcasts result.
     * 
     * @param graph Communication topology
     * @param initial_values Values at each vertex
     * @param reduction_op Reduction operation
     * @return AllReduceResult with communication pattern
     */
    static AllReduceResult<T> allReduceTree(const UnidirectionalGraph& graph,
                                           const std::map<int, T>& initial_values,
                                           const ReductionOp<T>& reduction_op);
    
    /**
     * @brief Performs all-reduce using hypercube communication pattern.
     * 
     * Requires graph to be a hypercube topology for optimal performance.
     * 
     * @param graph Hypercube topology graph
     * @param initial_values Values at each vertex
     * @param reduction_op Reduction operation
     * @return AllReduceResult with hypercube communication
     */
    static AllReduceResult<T> allReduceHypercube(const UnidirectionalGraph& graph,
                                                const std::map<int, T>& initial_values,
                                                const ReductionOp<T>& reduction_op);
    
    /**
     * @brief Performs all-reduce using ring topology.
     * 
     * Uses reduce-scatter followed by all-gather pattern.
     * 
     * @param graph Ring topology graph
     * @param initial_values Values at each vertex
     * @param reduction_op Reduction operation
     * @return AllReduceResult with ring communication
     */
    static AllReduceResult<T> allReduceRing(const UnidirectionalGraph& graph,
                                           const std::map<int, T>& initial_values,
                                           const ReductionOp<T>& reduction_op);

private:
    /**
     * @brief Performs tree reduction phase (leaf to root).
     */
    static void performTreeReduction(std::vector<NodeState<T>>& nodes,
                                    const UnidirectionalGraph& spanning_tree,
                                    const ReductionOp<T>& reduction_op,
                                    AllReduceResult<T>& result);
    
    /**
     * @brief Performs tree broadcast phase (root to leaves).
     */
    static void performTreeBroadcast(std::vector<NodeState<T>>& nodes,
                                    const UnidirectionalGraph& spanning_tree,
                                    const T& final_value,
                                    AllReduceResult<T>& result);
    
    /**
     * @brief Finds root node for tree-based operations.
     */
    static int findTreeRoot(const UnidirectionalGraph& tree);
    
    /**
     * @brief Calculates tree depth for a given root.
     */
    static int calculateTreeDepth(const UnidirectionalGraph& tree, int root);
    

};

/**
 * @brief Analyzes graph topology to determine optimal all-reduce pattern.
 */
CommunicationPattern analyzeOptimalPattern(const UnidirectionalGraph& graph);

/**
 * @brief Validates that graph topology supports the requested pattern.
 */
bool validateTopologyForPattern(const UnidirectionalGraph& graph, CommunicationPattern pattern);

/**
 * @brief Calculates efficiency metrics for all-reduce performance.
 */
template<typename T>
double calculateEfficiencyMetric(const AllReduceResult<T>& result, int num_nodes);

/**
 * @brief Validates initial values for all-reduce operation.
 */
template<typename T>
void validateInitialValues(const UnidirectionalGraph& graph, const std::map<int, T>& initial_values);

}  // namespace GraphCollective

#endif  // SRC_GRAPH_COLLECTIVE_H_
