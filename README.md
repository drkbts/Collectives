# Graph Algorithms Library

A comprehensive C++ library for graph algorithms and analysis, featuring both exact and approximate solutions for fundamental graph problems.

## Overview

This library provides a modular collection of graph algorithms built around a `UnidirectionalGraph` data structure. It combines custom implementations with industry-standard algorithms from the Boost Graph Library to offer both precision and performance.

## Key Features

### 🔗 **Graph Data Structure**
- **UnidirectionalGraph**: Efficient directed graph representation using adjacency lists
- **O(1) vertex and edge operations** with unordered containers
- **Dynamic graph construction** with automatic vertex creation

### 🧭 **Connectivity Analysis**
- **Graph connectivity checking** (treats directed graphs as weakly connected)
- **Vertex enumeration** and graph property analysis
- **Connected component detection**

### 🌳 **Spanning Tree Algorithms**
- **Spanning Tree**: DFS-based spanning tree computation for connected graphs
- **Spanning Forest**: Multi-component spanning forest for disconnected graphs
- **Error handling** for invalid inputs (empty/disconnected graphs)

### 🔄 **Hamiltonian Cycle Algorithms**
- **Exact Solution**: Backtracking algorithm for small graphs (< 15 vertices)
- **Approximate Solution**: 2-approximation using Boost's TSP algorithm
- **Performance flexibility**: Choose speed vs. optimality based on requirements

### 🚀 **TSP & Optimization Algorithms**
- **TSP Approximation**: Boost Graph Library integration for larger instances
- **Minimum Spanning Trees**: Kruskal's and Prim's algorithms
- **Weighted graph support** with customizable edge weights

### 🔀 **Graph Partitioning**
- **Bisectional Bandwidth**: Exact computation for small graphs (≤20 vertices)
- **Balanced Partitioning**: Finds optimal graph cuts into equal halves
- **Multiple Algorithms**: Brute force exact method with approximation methods planned
- **Known Topology Support**: Optimized solutions for complete graphs and trees

## Architecture

The library follows a modular design with clear separation of concerns:

```
src/
├── graph.h/.cc                   # Core UnidirectionalGraph data structure
├── graph_utils.h/.cc             # Graph utility functions
├── graph_connectivity.h/.cc      # Connectivity analysis algorithms
├── spanning_tree.h/.cc           # Spanning tree and forest algorithms
├── hamiltonian.h/.cc             # Hamiltonian cycle algorithms
├── tsp_boost.h/.cc               # TSP algorithms & Boost integration
├── graph_partitioning.h/.cc      # Graph partitioning & bisectional bandwidth
└── graph_collective.h/.cc        # Collective operations & all-reduce algorithms
```

### Module Responsibilities

| Module | Purpose | Key Functions |
|--------|---------|---------------|
| **Core Graph** | Basic data structure | `addVertex()`, `addEdge()`, `hasEdge()` |
| **Connectivity** | Graph analysis | `isConnected()`, `getAllVertices()` |
| **Spanning Tree** | Tree algorithms | `computeSpanningTree()`, `computeSpanningForest()` |
| **Hamiltonian** | Cycle algorithms | `findExactHamiltonianCycle()`, `findApproximateHamiltonianCycle()` |
| **TSP Boost** | Optimization | `computeMSTKruskal()`, `computeMSTPrim()` |
| **Graph Partitioning** | Bisectional bandwidth | `computeBisectionalBandwidth()`, `findBisectionalPartition()` |

## Usage Examples

### Basic Graph Operations
```cpp
#include "graph.h"
#include "graph_connectivity.h"

// Create and populate a graph
UnidirectionalGraph graph;
graph.addEdge(0, 1);
graph.addEdge(1, 2);
graph.addEdge(2, 0);

// Check connectivity
bool connected = GraphConnectivity::isConnected(graph);
std::vector<int> vertices = GraphConnectivity::getAllVertices(graph);
```

### Spanning Tree Computation
```cpp
#include "spanning_tree.h"

// Compute spanning tree for connected graph
UnidirectionalGraph tree = SpanningTree::computeSpanningTree(graph);

// Compute spanning forest for any graph
UnidirectionalGraph forest = SpanningTree::computeSpanningForest(graph);
```

### Hamiltonian Cycle Finding
```cpp
#include "hamiltonian.h"

// Exact solution (good for small graphs)
std::vector<int> exact_cycle = Hamiltonian::findExactHamiltonianCycle(graph);

// Fast approximation (good for larger graphs)
std::vector<int> approx_cycle = Hamiltonian::findApproximateHamiltonianCycle(graph);
```

### Advanced TSP and MST Algorithms
```cpp
#include "tsp_boost.h"

// Minimum spanning trees
UnidirectionalGraph mst_kruskal = TSPBoost::computeMSTKruskal(graph);
UnidirectionalGraph mst_prim = TSPBoost::computeMSTPrim(graph);

// TSP with custom weights
std::vector<std::vector<double>> weights = {{0, 1, 2}, {1, 0, 1}, {2, 1, 0}};
std::vector<int> tsp_tour = TSPBoost::findApproximateHamiltonianCycle(graph, weights);
```

### Graph Partitioning and Bisectional Bandwidth
```cpp
#include "graph_partitioning.h"

// Compute bisectional bandwidth (automatic method selection)
int bandwidth = GraphPartitioning::computeBisectionalBandwidth(graph);

// Find the actual partition that achieves minimum bandwidth
GraphPartitioning::GraphPartition partition = GraphPartitioning::findBisectionalPartition(graph);
std::vector<int> partition1 = partition.partition1;
std::vector<int> partition2 = partition.partition2;
int cut_size = partition.bandwidth;

// Exact computation for small graphs
int exact_bandwidth = GraphPartitioning::computeExactBisectionalBandwidth(graph);

// Known topology optimizations
int tree_bandwidth = GraphPartitioning::computeKnownTopologyBandwidth(graph, 
    GraphPartitioning::GraphTopology::TREE);
```

### Collective Operations and All-Reduce
```cpp
#include "graph_collective.h"

// Basic all-reduce with automatic pattern selection
std::map<int, int> initial_values = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
GraphCollective::SumReduction<int> sum_op;

auto result = GraphCollective::GraphCollective<int>::allReduce(graph, initial_values, sum_op);
int final_value = result.final_value;  // Sum of all values
int total_rounds = result.total_rounds;
int total_messages = result.total_messages;

// Tree-based all-reduce (good for general graphs)
auto tree_result = GraphCollective::GraphCollective<int>::allReduceTree(graph, initial_values, sum_op);

// Hypercube all-reduce (optimal for power-of-2 nodes)
auto hypercube_result = GraphCollective::GraphCollective<int>::allReduceHypercube(graph, initial_values, sum_op);

// Ring-based all-reduce (good for bandwidth-limited scenarios)
auto ring_result = GraphCollective::GraphCollective<int>::allReduceRing(graph, initial_values, sum_op);

// Different reduction operations
GraphCollective::MaxReduction<int> max_op;
GraphCollective::MinReduction<int> min_op;
GraphCollective::AverageReduction<double> avg_op(4);

auto max_result = GraphCollective::GraphCollective<int>::allReduce(graph, initial_values, max_op);
auto min_result = GraphCollective::GraphCollective<int>::allReduce(graph, initial_values, min_op);

// Custom communication patterns
auto optimal_result = GraphCollective::GraphCollective<int>::allReduce(
    graph, initial_values, sum_op, GraphCollective::CommunicationPattern::OPTIMAL);
```

## Performance Characteristics

| Algorithm | Time Complexity | Best Use Case |
|-----------|----------------|---------------|
| **Connectivity Check** | O(V + E) | All graph sizes |
| **Spanning Tree** | O(V + E) | Connected graphs |
| **Exact Hamiltonian** | O(V!) | Small graphs (≤ 15 vertices) |
| **Approximate Hamiltonian** | O(V²) | Larger graphs, approximation acceptable |
| **MST (Kruskal/Prim)** | O(E log V) | Weighted graphs |
| **Bisectional Bandwidth** | O(2^V) | Small graphs (≤ 20 vertices) |

## Build System

This project uses **Bazel** with **Bzlmod** for dependency management.

### Prerequisites
- **Bazel 7.0+** with Bzlmod support
- **C++17 compatible compiler**
- **Boost Graph Library 1.87.0** (automatically managed)

### Building
```bash
# Build all targets
bazel build //...

# Build specific components
bazel build //:graph
```

### Testing
```bash
# Run all tests
bazel test //...

# Run specific test suites
bazel test //:spanning_tree_test
bazel test //:graph_test
bazel test //:graph_utils_test
bazel test //:graph_partitioning_test
```

### Test Coverage
- **107 comprehensive tests** covering all algorithms
- **Edge case handling**: empty graphs, single vertices, disconnected components
- **Performance validation**: algorithm correctness and efficiency
- **Integration testing**: module interaction verification

## Dependencies

- **Boost Graph Library 1.87.0**: TSP approximation and MST algorithms
- **Google Test**: Testing framework
- **Bazel Rules CC**: C++ build rules

All dependencies are managed automatically through Bazel's module system.

## Design Principles

### **Modularity**
Each algorithm family is in its own module with clear interfaces and minimal coupling.

### **Performance Flexibility**
Multiple algorithm implementations (exact vs. approximate) allow choosing the right trade-off for each use case.

### **Extensibility**
Clean module boundaries make it easy to add new algorithms without affecting existing code.

### **Robustness**
Comprehensive error handling and edge case coverage ensure reliable operation.

## Contributing

This library demonstrates modern C++ practices:
- **RAII resource management**
- **Clear module boundaries**
- **Comprehensive testing**
- **Modern build system (Bazel + Bzlmod)**
- **Industry-standard library integration (Boost)**

## Examples in Action

The library has been developed to handle various graph scenarios:

- **Small complete graphs**: Exact Hamiltonian cycle finding
- **Large sparse graphs**: Approximate TSP solutions
- **Disconnected networks**: Spanning forest computation
- **Weighted optimization**: MST algorithms with custom costs

This makes it suitable for research, education, and practical applications in network analysis, route optimization, and graph theory studies.
