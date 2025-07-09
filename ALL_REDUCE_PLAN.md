# All-Reduce Functionality Implementation Plan

## Background & Motivation

**All-Reduce** is a collective communication operation where:
1. Each process/node contributes data
2. A reduction operation (sum, max, min, etc.) is performed across all contributions
3. The result is distributed back to all processes/nodes

In the context of our graph library, this enables:
- **Distributed graph algorithms** (PageRank, consensus algorithms, etc.)
- **Graph-based communication patterns** (tree, hypercube, butterfly networks)
- **Collective operations on graph vertices** (global aggregation, synchronization)
- **Simulation of distributed systems** communication patterns

## New Module: graph_collective.h/.cc

### Purpose
Implement collective communication operations using graph topologies as communication patterns.

### Module Location
```
src/
├── graph_collective.h/.cc        # New module for collective operations
```

### Namespace
```cpp
namespace GraphCollective {
    // All collective communication functions
}
```

## Core Concepts & Data Structures

### Reduction Operations
```cpp
template<typename T>
class ReductionOp {
public:
    virtual ~ReductionOp() = default;
    virtual T operator()(const T& a, const T& b) const = 0;
    virtual T identity() const = 0;  // Identity element for the operation
};

// Predefined reduction operations
template<typename T>
class SumReduction : public ReductionOp<T> {
    T operator()(const T& a, const T& b) const override { return a + b; }
    T identity() const override { return T{}; }
};

template<typename T>
class MaxReduction : public ReductionOp<T> {
    T operator()(const T& a, const T& b) const override { return std::max(a, b); }
    T identity() const override { return std::numeric_limits<T>::lowest(); }
};

template<typename T>
class MinReduction : public ReductionOp<T> {
    T operator()(const T& a, const T& b) const override { return std::min(a, b); }
    T identity() const override { return std::numeric_limits<T>::max(); }
};
```

### Communication Patterns
```cpp
enum class CommunicationPattern {
    TREE,           // Binary tree pattern
    HYPERCUBE,      // Hypercube communication
    BUTTERFLY,      // Butterfly network
    RING,           // Ring topology
    MESH_2D,        // 2D mesh
    CUSTOM_GRAPH,   // Use provided graph structure
    OPTIMAL         // Choose optimal pattern based on graph
};

struct CommunicationStep {
    int source_node;
    int target_node;
    int round;              // Communication round number
    std::vector<int> data;  // Simulated data being communicated
};

struct AllReduceResult {
    std::vector<CommunicationStep> communication_steps;
    int total_rounds;
    int total_messages;
    double efficiency_metric;
};
```

### Node State Management
```cpp
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
```

## Core API Functions

### Primary All-Reduce Functions
```cpp
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
    static AllReduceResult allReduce(const UnidirectionalGraph& graph,
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
    static AllReduceResult allReduceTree(const UnidirectionalGraph& graph,
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
    static AllReduceResult allReduceHypercube(const UnidirectionalGraph& graph,
                                            const std::map<int, T>& initial_values,
                                            const ReductionOp<T>& reduction_op);
    
    /**
     * @brief Performs all-reduce using butterfly network pattern.
     * 
     * Efficient for power-of-2 number of nodes.
     * 
     * @param graph Butterfly network topology
     * @param initial_values Values at each vertex
     * @param reduction_op Reduction operation
     * @return AllReduceResult with butterfly communication
     */
    static AllReduceResult allReduceButterfly(const UnidirectionalGraph& graph,
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
    static AllReduceResult allReduceRing(const UnidirectionalGraph& graph,
                                       const std::map<int, T>& initial_values,
                                       const ReductionOp<T>& reduction_op);
};
```

### Utility Functions
```cpp
namespace GraphCollective {
    /**
     * @brief Analyzes graph topology to determine optimal all-reduce pattern.
     */
    CommunicationPattern analyzeOptimalPattern(const UnidirectionalGraph& graph);
    
    /**
     * @brief Generates spanning tree for tree-based all-reduce.
     */
    UnidirectionalGraph generateSpanningTreeForAllReduce(const UnidirectionalGraph& graph);
    
    /**
     * @brief Validates that graph topology supports the requested pattern.
     */
    bool validateTopologyForPattern(const UnidirectionalGraph& graph, CommunicationPattern pattern);
    
    /**
     * @brief Simulates one round of communication for all-reduce.
     */
    template<typename T>
    void simulateCommunicationRound(std::vector<NodeState<T>>& nodes,
                                   const std::vector<CommunicationStep>& steps,
                                   const ReductionOp<T>& reduction_op);
    
    /**
     * @brief Calculates efficiency metrics for all-reduce performance.
     */
    double calculateEfficiencyMetric(const AllReduceResult& result, int num_nodes);
}
```

## Implementation Strategy

### Phase 1: Core Infrastructure (Week 1)
```cpp
// Basic all-reduce with tree pattern
template<typename T>
AllReduceResult GraphCollective::allReduceTree(const UnidirectionalGraph& graph,
                                               const std::map<int, T>& initial_values,
                                               const ReductionOp<T>& reduction_op) {
    AllReduceResult result;
    
    // Phase 1: Build spanning tree
    UnidirectionalGraph spanning_tree = SpanningTree::computeSpanningTree(graph);
    
    // Phase 2: Reduce phase (leaf to root)
    std::vector<NodeState<T>> nodes;
    for (const auto& [node_id, value] : initial_values) {
        nodes.emplace_back(node_id, value);
    }
    
    // Perform reduction up the tree
    performTreeReduction(nodes, spanning_tree, reduction_op, result);
    
    // Phase 3: Broadcast phase (root to leaves)
    performTreeBroadcast(nodes, spanning_tree, result);
    
    return result;
}
```

### Phase 2: Specialized Patterns (Week 2)
```cpp
// Hypercube all-reduce implementation
template<typename T>
AllReduceResult GraphCollective::allReduceHypercube(const UnidirectionalGraph& graph,
                                                   const std::map<int, T>& initial_values,
                                                   const ReductionOp<T>& reduction_op) {
    AllReduceResult result;
    
    // Validate hypercube topology
    if (!validateHypercubeTopology(graph)) {
        throw std::invalid_argument("Graph is not a valid hypercube topology");
    }
    
    int dimensions = calculateHypercubeDimensions(graph);
    std::vector<NodeState<T>> nodes;
    
    // Initialize node states
    for (const auto& [node_id, value] : initial_values) {
        nodes.emplace_back(node_id, value);
    }
    
    // Perform log(n) rounds of communication
    for (int dim = 0; dim < dimensions; ++dim) {
        std::vector<CommunicationStep> round_steps;
        
        for (auto& node : nodes) {
            int partner = node.getNodeId() ^ (1 << dim);  // XOR to find partner
            
            // Exchange and reduce with partner
            CommunicationStep step;
            step.source_node = node.getNodeId();
            step.target_node = partner;
            step.round = dim;
            round_steps.push_back(step);
        }
        
        simulateCommunicationRound(nodes, round_steps, reduction_op);
        result.communication_steps.insert(result.communication_steps.end(),
                                         round_steps.begin(), round_steps.end());
    }
    
    result.total_rounds = dimensions;
    result.total_messages = nodes.size() * dimensions;
    result.efficiency_metric = calculateEfficiencyMetric(result, nodes.size());
    
    return result;
}
```

### Phase 3: Advanced Patterns (Week 3)
```cpp
// Butterfly network implementation
template<typename T>
AllReduceResult GraphCollective::allReduceButterfly(const UnidirectionalGraph& graph,
                                                   const std::map<int, T>& initial_values,
                                                   const ReductionOp<T>& reduction_op) {
    // Implement butterfly network all-reduce
    // More complex but potentially more efficient for certain topologies
}

// Ring-based all-reduce (reduce-scatter + all-gather)
template<typename T>
AllReduceResult GraphCollective::allReduceRing(const UnidirectionalGraph& graph,
                                               const std::map<int, T>& initial_values,
                                               const ReductionOp<T>& reduction_op) {
    // Implement ring-based all-reduce
    // Good for bandwidth-limited scenarios
}
```

## Integration with Existing Modules

### Dependencies
```cpp
#include "graph.h"                    // Core graph structure
#include "graph_connectivity.h"      // For connectivity analysis
#include "spanning_tree.h"           // For tree-based all-reduce
#include "graph_partitioning.h"      // For load balancing
```

### Module Relationships
```
graph.h/.cc
    ↓
graph_connectivity.h/.cc
    ↓
spanning_tree.h/.cc
    ↓
graph_collective.h/.cc (NEW)
    ↓
Applications (PageRank, distributed consensus, etc.)
```

## Use Cases & Applications

### 1. Distributed PageRank
```cpp
// Use all-reduce to compute PageRank in distributed manner
template<typename T>
std::vector<T> distributedPageRank(const UnidirectionalGraph& graph,
                                  const std::vector<T>& initial_ranks,
                                  double damping_factor = 0.85,
                                  int iterations = 100) {
    auto current_ranks = initial_ranks;
    SumReduction<T> sum_op;
    
    for (int iter = 0; iter < iterations; ++iter) {
        // Compute new ranks using all-reduce
        std::map<int, T> rank_contributions;
        for (int vertex : graph.getVertices()) {
            rank_contributions[vertex] = computePageRankContribution(vertex, current_ranks);
        }
        
        auto result = GraphCollective::allReduce(graph, rank_contributions, sum_op);
        // Update ranks based on result
        updatePageRanks(current_ranks, result, damping_factor);
    }
    
    return current_ranks;
}
```

### 2. Consensus Algorithms
```cpp
// Byzantine fault-tolerant consensus using all-reduce
template<typename T>
T distributedConsensus(const UnidirectionalGraph& graph,
                      const std::map<int, T>& initial_values,
                      int fault_tolerance = 1) {
    // Use all-reduce with custom reduction operation for consensus
    ConsensusReduction<T> consensus_op(fault_tolerance);
    auto result = GraphCollective::allReduce(graph, initial_values, consensus_op);
    return extractConsensusValue(result);
}
```

### 3. Distributed Optimization
```cpp
// Gradient descent using all-reduce for parameter updates
template<typename T>
std::vector<T> distributedGradientDescent(const UnidirectionalGraph& graph,
                                         const std::map<int, std::vector<T>>& local_gradients,
                                         double learning_rate = 0.01) {
    VectorSumReduction<T> gradient_sum;
    auto result = GraphCollective::allReduce(graph, local_gradients, gradient_sum);
    
    // Apply averaged gradient update
    std::vector<T> averaged_gradient = result.final_value;
    for (auto& grad : averaged_gradient) {
        grad /= graph.getVertexCount();
    }
    
    return averaged_gradient;
}
```

## Performance Characteristics

### Algorithm Complexity
| Pattern | Time Complexity | Message Complexity | Best For |
|---------|-----------------|-------------------|----------|
| **Tree** | O(log n) | O(n) | General purpose |
| **Hypercube** | O(log n) | O(n log n) | Power-of-2 nodes |
| **Butterfly** | O(log n) | O(n log n) | High bandwidth |
| **Ring** | O(n) | O(n) | Bandwidth-limited |

### Communication Efficiency
```cpp
// Metrics for evaluating all-reduce performance
struct PerformanceMetrics {
    int total_communication_rounds;
    int total_messages_sent;
    double bandwidth_utilization;
    double latency_overhead;
    double fault_tolerance_level;
};
```

## Testing Strategy

### Test Coverage
```cpp
// Test file: tests/graph_collective_test.cc

class GraphCollectiveTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create test topologies
        createHypercubeGraph();
        createButterflyGraph();
        createRingGraph();
        createTreeGraph();
    }
    
    UnidirectionalGraph hypercube_graph;
    UnidirectionalGraph butterfly_graph;
    UnidirectionalGraph ring_graph;
    UnidirectionalGraph tree_graph;
};

// Test cases
TEST_F(GraphCollectiveTest, AllReduceTreeSum);
TEST_F(GraphCollectiveTest, AllReduceHypercubeMax);
TEST_F(GraphCollectiveTest, AllReduceButterflyMin);
TEST_F(GraphCollectiveTest, AllReduceRingAverage);
TEST_F(GraphCollectiveTest, PerformanceComparison);
TEST_F(GraphCollectiveTest, FaultTolerance);
```

### Example Tests
```cpp
TEST_F(GraphCollectiveTest, AllReduceTreeSum) {
    // Test tree-based all-reduce with sum operation
    std::map<int, int> initial_values = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
    SumReduction<int> sum_op;
    
    auto result = GraphCollective::allReduceTree(tree_graph, initial_values, sum_op);
    
    EXPECT_EQ(result.final_value, 10);  // 1+2+3+4 = 10
    EXPECT_LE(result.total_rounds, 3);  // Should be log(n) rounds
    EXPECT_EQ(result.total_messages, 6); // 2*(n-1) messages for tree
}
```

## Future Extensions

### Advanced Features
- **Fault-tolerant all-reduce**: Handle node failures during communication
- **Adaptive patterns**: Dynamically choose optimal communication pattern
- **Bandwidth-aware scheduling**: Optimize for network bandwidth utilization
- **Heterogeneous networks**: Support for different node capabilities

### Integration Opportunities
- **Graph partitioning**: Use partitioning for load balancing in distributed algorithms
- **TSP optimization**: Distributed optimization using all-reduce
- **Hamiltonian cycle**: Distributed search algorithms
- **Spanning tree**: Multiple spanning trees for fault tolerance

## Build System Changes

### BUILD File Updates
```bazel
cc_library(
    name = "graph",
    srcs = [
        # ... existing files ...
        "src/graph_collective.cc",
    ],
    hdrs = [
        # ... existing headers ...
        "src/graph_collective.h",
    ],
    # ... existing deps ...
)

cc_test(
    name = "graph_collective_test",
    srcs = ["tests/graph_collective_test.cc"],
    deps = [
        ":graph",
        "@googletest//:gtest_main",
    ],
)
```

## Implementation Timeline

### Phase 1 (Week 1): Foundation
- [ ] Create module structure and basic types
- [ ] Implement tree-based all-reduce
- [ ] Add basic reduction operations (sum, max, min)
- [ ] Create initial test suite

### Phase 2 (Week 2): Core Patterns
- [ ] Implement hypercube all-reduce
- [ ] Add butterfly network pattern
- [ ] Implement ring-based all-reduce
- [ ] Performance benchmarking

### Phase 3 (Week 3): Advanced Features
- [ ] Add automatic pattern selection
- [ ] Implement fault tolerance mechanisms
- [ ] Add distributed algorithm examples
- [ ] Comprehensive testing

### Phase 4 (Week 4): Integration & Polish
- [ ] Integration with existing modules
- [ ] Performance optimization
- [ ] Documentation and examples
- [ ] Real-world application demos

This plan provides a comprehensive approach to adding all-reduce functionality while maintaining clean integration with the existing modular architecture and enabling powerful distributed computing capabilities.
