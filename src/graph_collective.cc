#include "graph_collective.h"
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cmath>

namespace GraphCollective {

template<typename T>
void validateInitialValues(const UnidirectionalGraph& graph, const std::map<int, T>& initial_values) {
    std::vector<int> vertices = GraphConnectivity::getAllVertices(graph);
    if (vertices.empty()) {
        throw std::invalid_argument("Cannot perform all-reduce on empty graph");
    }
    
    if (initial_values.empty()) {
        throw std::invalid_argument("Initial values cannot be empty");
    }
    
    // Check that all vertices have initial values
    for (int vertex : vertices) {
        if (initial_values.find(vertex) == initial_values.end()) {
            throw std::invalid_argument("Missing initial value for vertex " + std::to_string(vertex));
        }
    }
}

template<typename T>
double calculateEfficiencyMetric(const AllReduceResult<T>& result, int num_nodes) {
    if (num_nodes <= 1) return 1.0;
    
    // Ideal metrics for comparison
    int ideal_rounds = static_cast<int>(std::ceil(std::log2(num_nodes)));
    int ideal_messages = num_nodes - 1;  // Minimum messages for any reduction
    
    // Calculate efficiency based on rounds and messages
    double round_efficiency = static_cast<double>(ideal_rounds) / std::max(1, result.total_rounds);
    double message_efficiency = static_cast<double>(ideal_messages) / std::max(1, result.total_messages);
    
    // Combined efficiency metric (weighted average)
    return 0.7 * round_efficiency + 0.3 * message_efficiency;
}

bool validateTopologyForPattern(const UnidirectionalGraph& graph, CommunicationPattern pattern) {
    std::vector<int> vertices = GraphConnectivity::getAllVertices(graph);
    int n = vertices.size();
    
    switch (pattern) {
        case CommunicationPattern::TREE:
            // Any connected graph can use tree pattern
            return GraphConnectivity::isConnected(graph);
            
        case CommunicationPattern::HYPERCUBE: {
            // Must be power of 2 vertices and form a hypercube
            if (n == 0 || ((n & (n - 1)) != 0)) {
                return false;  // Not power of 2
            }
            
            // Check if it's actually a hypercube topology
            int dimensions = static_cast<int>(std::log2(n));
            for (int i = 0; i < n; i++) {
                auto neighbors = graph.getNeighbors(i);
                if (neighbors.size() != dimensions) {
                    return false;  // Each node should have exactly 'dimensions' neighbors
                }
                
                // Check if neighbors are exactly the nodes that differ by one bit
                for (int dim = 0; dim < dimensions; dim++) {
                    int expected_neighbor = i ^ (1 << dim);
                    if (std::find(neighbors.begin(), neighbors.end(), expected_neighbor) == neighbors.end()) {
                        return false;
                    }
                }
            }
            return true;
        }
            
        case CommunicationPattern::RING:
            // Any connected graph can use ring pattern
            return GraphConnectivity::isConnected(graph);
            
        case CommunicationPattern::MESH_2D: {
            // Must be perfect square number of vertices
            int sqrt_n = static_cast<int>(std::sqrt(n));
            return sqrt_n * sqrt_n == n;
        }
            
        case CommunicationPattern::CUSTOM_GRAPH:
        case CommunicationPattern::OPTIMAL:
            return true;  // Always valid
            
        default:
            return false;
    }
}

CommunicationPattern analyzeOptimalPattern(const UnidirectionalGraph& graph) {
    std::vector<int> vertices = GraphConnectivity::getAllVertices(graph);
    int n = vertices.size();
    
    // Check for hypercube topology (power of 2)
    if (n > 0 && ((n & (n - 1)) == 0) && validateTopologyForPattern(graph, CommunicationPattern::HYPERCUBE)) {
        return CommunicationPattern::HYPERCUBE;
    }
    
    // Check for 2D mesh topology
    int sqrt_n = static_cast<int>(std::sqrt(n));
    if (sqrt_n * sqrt_n == n && validateTopologyForPattern(graph, CommunicationPattern::MESH_2D)) {
        return CommunicationPattern::MESH_2D;
    }
    
    // Default to tree pattern for general graphs
    return CommunicationPattern::TREE;
}

// Template implementations
template<typename T>
AllReduceResult<T> GraphCollective<T>::allReduce(const UnidirectionalGraph& graph,
                                                 const std::map<int, T>& initial_values,
                                                 const ReductionOp<T>& reduction_op,
                                                 CommunicationPattern pattern) {
    validateInitialValues(graph, initial_values);
    
    // Choose pattern if AUTO
    if (pattern == CommunicationPattern::OPTIMAL) {
        pattern = analyzeOptimalPattern(graph);
    }
    
    // Validate topology for chosen pattern
    if (!validateTopologyForPattern(graph, pattern)) {
        throw std::invalid_argument("Graph topology does not support requested communication pattern");
    }
    
    // Dispatch to appropriate implementation
    switch (pattern) {
        case CommunicationPattern::TREE:
            return allReduceTree(graph, initial_values, reduction_op);
            
        case CommunicationPattern::HYPERCUBE:
            return allReduceHypercube(graph, initial_values, reduction_op);
            
        case CommunicationPattern::RING:
            return allReduceRing(graph, initial_values, reduction_op);
            
        default:
            throw std::runtime_error("Communication pattern not yet implemented");
    }
}

template<typename T>
int GraphCollective<T>::findTreeRoot(const UnidirectionalGraph& tree) {
    std::vector<int> vertices = GraphConnectivity::getAllVertices(tree);
    if (vertices.empty()) {
        throw std::invalid_argument("Cannot find root of empty tree");
    }
    
    // Choose vertex with minimum eccentricity (center of tree)
    int best_root = vertices[0];
    int min_depth = INT_MAX;
    
    for (int candidate : vertices) {
        int depth = calculateTreeDepth(tree, candidate);
        if (depth < min_depth) {
            min_depth = depth;
            best_root = candidate;
        }
    }
    
    return best_root;
}

template<typename T>
int GraphCollective<T>::calculateTreeDepth(const UnidirectionalGraph& tree, int root) {
    std::unordered_map<int, int> depths;
    std::queue<int> queue;
    
    queue.push(root);
    depths[root] = 0;
    
    int max_depth = 0;
    while (!queue.empty()) {
        int current = queue.front();
        queue.pop();
        
        std::vector<int> neighbors = tree.getNeighbors(current);
        for (int neighbor : neighbors) {
            if (depths.find(neighbor) == depths.end()) {
                depths[neighbor] = depths[current] + 1;
                max_depth = std::max(max_depth, depths[neighbor]);
                queue.push(neighbor);
            }
        }
    }
    
    return max_depth;
}

template<typename T>
void GraphCollective<T>::performTreeReduction(std::vector<NodeState<T>>& nodes,
                                              const UnidirectionalGraph& spanning_tree,
                                              const ReductionOp<T>& reduction_op,
                                              AllReduceResult<T>& result) {
    std::vector<int> vertices = GraphConnectivity::getAllVertices(spanning_tree);
    
    if (vertices.empty()) {
        return;
    }
    
    // Create node value mapping
    std::unordered_map<int, T> node_values;
    for (const auto& node : nodes) {
        node_values[node.getNodeId()] = node.getValue();
    }
    
    // For 2-node case, just combine the values
    if (vertices.size() == 2) {
        result.total_rounds = 1;
        result.total_messages = 1;
        result.communication_steps.emplace_back(
            vertices[1], vertices[0], 1, "REDUCE " + reduction_op.name());
        
        result.final_value = reduction_op(node_values[vertices[0]], node_values[vertices[1]]);
        return;
    }
    
    // The spanning tree is already a directed graph from parent to child
    // We can use it directly to build the parent-child relationships
    std::unordered_map<int, int> parent;
    std::unordered_map<int, std::vector<int>> children;
    
    // Find root (node with no incoming edges in spanning tree)
    int root = -1;
    for (int vertex : vertices) {
        bool has_incoming = false;
        for (int other : vertices) {
            if (other != vertex && spanning_tree.hasEdge(other, vertex)) {
                has_incoming = true;
                break;
            }
        }
        if (!has_incoming) {
            root = vertex;
            break;
        }
    }
    
    if (root == -1) {
        root = vertices[0];  // Fallback
    }
    
    // Build parent-child relationships from spanning tree
    for (int vertex : vertices) {
        parent[vertex] = -1;
        for (int neighbor : spanning_tree.getNeighbors(vertex)) {
            parent[neighbor] = vertex;
            children[vertex].push_back(neighbor);
        }
    }
    
    // Perform post-order traversal (children to parent)
    std::function<void(int)> reduceSubtree = [&](int node) {
        // First, process all children
        for (int child : children[node]) {
            reduceSubtree(child);
            
            // Child sends its value to parent
            result.communication_steps.emplace_back(
                child, node, result.total_rounds + 1,
                "REDUCE " + reduction_op.name());
            result.total_messages++;
            
            // Combine child's value with parent's
            node_values[node] = reduction_op(node_values[node], node_values[child]);
        }
        
        if (!children[node].empty()) {
            result.total_rounds++;
        }
    };
    
    // Start reduction from root
    reduceSubtree(root);
    
    result.final_value = node_values[root];
}

template<typename T>
void GraphCollective<T>::performTreeBroadcast(std::vector<NodeState<T>>& nodes,
                                              const UnidirectionalGraph& spanning_tree,
                                              const T& final_value,
                                              AllReduceResult<T>& result) {
    // Build adjacency list for tree
    std::unordered_map<int, std::vector<int>> tree_adj;
    std::vector<int> vertices = GraphConnectivity::getAllVertices(spanning_tree);
    
    for (int vertex : vertices) {
        tree_adj[vertex] = spanning_tree.getNeighbors(vertex);
    }
    
    // Find root and establish parent-child relationships
    int root = findTreeRoot(spanning_tree);
    std::unordered_map<int, std::vector<int>> children;
    
    // BFS to establish parent-child relationships
    std::queue<int> queue;
    std::unordered_set<int> visited;
    
    queue.push(root);
    visited.insert(root);
    
    while (!queue.empty()) {
        int current = queue.front();
        queue.pop();
        
        for (int neighbor : tree_adj[current]) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                children[current].push_back(neighbor);
                queue.push(neighbor);
            }
        }
    }
    
    // Broadcast from root to leaves
    std::queue<int> broadcast_queue;
    broadcast_queue.push(root);
    
    while (!broadcast_queue.empty()) {
        int current = broadcast_queue.front();
        broadcast_queue.pop();
        
        if (!children[current].empty()) {
            result.total_rounds++;
            
            for (int child : children[current]) {
                result.communication_steps.emplace_back(
                    current, child, result.total_rounds, "BROADCAST");
                result.total_messages++;
                broadcast_queue.push(child);
            }
        }
    }
    
    // Update all node values to final result
    for (auto& node : nodes) {
        node.updateValue(final_value);
    }
}

template<typename T>
AllReduceResult<T> GraphCollective<T>::allReduceTree(const UnidirectionalGraph& graph,
                                                    const std::map<int, T>& initial_values,
                                                    const ReductionOp<T>& reduction_op) {
    validateInitialValues(graph, initial_values);
    
    AllReduceResult<T> result;
    
    // Handle single node case
    if (initial_values.size() == 1) {
        result.final_value = initial_values.begin()->second;
        result.efficiency_metric = 1.0;
        return result;
    }
    
    // Build spanning tree
    UnidirectionalGraph spanning_tree = SpanningTree::computeSpanningTree(graph);
    
    // Initialize node states
    std::vector<NodeState<T>> nodes;
    for (const auto& [node_id, value] : initial_values) {
        nodes.emplace_back(node_id, value);
    }
    
    // Phase 1: Reduce phase (leaf to root)
    performTreeReduction(nodes, spanning_tree, reduction_op, result);
    
    // Phase 2: Broadcast phase (root to leaves)
    performTreeBroadcast(nodes, spanning_tree, result.final_value, result);
    
    // Calculate efficiency metrics
    result.efficiency_metric = calculateEfficiencyMetric(result, initial_values.size());
    
    return result;
}

template<typename T>
AllReduceResult<T> GraphCollective<T>::allReduceHypercube(const UnidirectionalGraph& graph,
                                                         const std::map<int, T>& initial_values,
                                                         const ReductionOp<T>& reduction_op) {
    validateInitialValues(graph, initial_values);
    
    std::vector<int> vertices = GraphConnectivity::getAllVertices(graph);
    int n = vertices.size();
    
    // Validate hypercube topology
    if (n == 0 || ((n & (n - 1)) != 0)) {
        throw std::invalid_argument("Graph must have power-of-2 vertices for hypercube all-reduce");
    }
    
    AllReduceResult<T> result;
    int dimensions = static_cast<int>(std::log2(n));
    
    // Initialize node values
    std::unordered_map<int, T> node_values;
    for (const auto& [node_id, value] : initial_values) {
        node_values[node_id] = value;
    }
    
    // Perform log(n) rounds of communication
    for (int dim = 0; dim < dimensions; ++dim) {
        result.total_rounds++;
        
        // Store values for this round to avoid double-updates
        std::unordered_map<int, T> round_values = node_values;
        
        for (int node_id : vertices) {
            int partner = node_id ^ (1 << dim);  // XOR to find partner
            
            if (initial_values.find(partner) != initial_values.end()) {
                // Process each communication (both directions)
                result.communication_steps.emplace_back(
                    node_id, partner, result.total_rounds,
                    "HYPERCUBE " + reduction_op.name() + " dim=" + std::to_string(dim));
                result.total_messages++;
                
                // Combine values (but only update once per round)
                if (node_id < partner) {
                    T combined = reduction_op(round_values[node_id], round_values[partner]);
                    node_values[node_id] = combined;
                    node_values[partner] = combined;
                }
            }
        }
    }
    
    // All nodes should have the same final value
    if (!node_values.empty()) {
        result.final_value = node_values.begin()->second;
    }
    
    result.efficiency_metric = calculateEfficiencyMetric(result, n);
    
    return result;
}

template<typename T>
AllReduceResult<T> GraphCollective<T>::allReduceRing(const UnidirectionalGraph& graph,
                                                    const std::map<int, T>& initial_values,
                                                    const ReductionOp<T>& reduction_op) {
    validateInitialValues(graph, initial_values);
    
    std::vector<int> vertices = GraphConnectivity::getAllVertices(graph);
    int n = vertices.size();
    
    if (n <= 1) {
        AllReduceResult<T> result;
        if (n == 1) {
            result.final_value = initial_values.begin()->second;
        }
        result.efficiency_metric = 1.0;
        return result;
    }
    
    AllReduceResult<T> result;
    
    // Sort vertices to create ring order
    std::sort(vertices.begin(), vertices.end());
    
    // Initialize node values
    std::unordered_map<int, T> node_values;
    for (const auto& [node_id, value] : initial_values) {
        node_values[node_id] = value;
    }
    
    // Simplified ring all-reduce: just accumulate all values
    T final_value = reduction_op.identity();
    
    for (int i = 0; i < n; ++i) {
        int current = vertices[i];
        final_value = reduction_op(final_value, node_values[current]);
        
        if (i > 0) {
            result.communication_steps.emplace_back(
                current, vertices[0], 1, "RING_REDUCE " + reduction_op.name());
            result.total_messages++;
        }
    }
    
    result.total_rounds = 1;
    
    // Broadcast final value to all nodes
    for (int i = 1; i < n; ++i) {
        result.communication_steps.emplace_back(
            vertices[0], vertices[i], 2, "RING_BROADCAST");
        result.total_messages++;
        node_values[vertices[i]] = final_value;
    }
    
    if (n > 1) {
        result.total_rounds = 2;
    }
    
    node_values[vertices[0]] = final_value;
    
    // All nodes should have the same final value
    if (!node_values.empty()) {
        result.final_value = node_values.begin()->second;
    }
    
    result.efficiency_metric = calculateEfficiencyMetric(result, n);
    
    return result;
}

// Explicit template instantiations
template class GraphCollective<int>;
template class GraphCollective<double>;
template class GraphCollective<float>;
template class GraphCollective<long>;

template void validateInitialValues<int>(const UnidirectionalGraph& graph, const std::map<int, int>& initial_values);
template void validateInitialValues<double>(const UnidirectionalGraph& graph, const std::map<int, double>& initial_values);
template void validateInitialValues<float>(const UnidirectionalGraph& graph, const std::map<int, float>& initial_values);
template void validateInitialValues<long>(const UnidirectionalGraph& graph, const std::map<int, long>& initial_values);

template double calculateEfficiencyMetric<int>(const AllReduceResult<int>& result, int num_nodes);
template double calculateEfficiencyMetric<double>(const AllReduceResult<double>& result, int num_nodes);
template double calculateEfficiencyMetric<float>(const AllReduceResult<float>& result, int num_nodes);
template double calculateEfficiencyMetric<long>(const AllReduceResult<long>& result, int num_nodes);

}  // namespace GraphCollective
