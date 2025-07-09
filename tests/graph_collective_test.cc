#include <gtest/gtest.h>
#include "src/graph_collective.h"
#include "src/graph.h"
#include <map>
#include <vector>
#include <algorithm>

class GraphCollectiveTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple 4-vertex tree
        // 0 -- 1 -- 2 -- 3
        tree_graph.addEdge(0, 1);
        tree_graph.addEdge(1, 0);
        tree_graph.addEdge(1, 2);
        tree_graph.addEdge(2, 1);
        tree_graph.addEdge(2, 3);
        tree_graph.addEdge(3, 2);
        
        // Create a 4-vertex hypercube
        // 0 -- 1
        // |    |
        // 2 -- 3
        hypercube_graph.addEdge(0, 1);
        hypercube_graph.addEdge(1, 0);
        hypercube_graph.addEdge(0, 2);
        hypercube_graph.addEdge(2, 0);
        hypercube_graph.addEdge(1, 3);
        hypercube_graph.addEdge(3, 1);
        hypercube_graph.addEdge(2, 3);
        hypercube_graph.addEdge(3, 2);
        
        // Create a 4-vertex ring
        // 0 -- 1
        // |    |
        // 3 -- 2
        ring_graph.addEdge(0, 1);
        ring_graph.addEdge(1, 0);
        ring_graph.addEdge(1, 2);
        ring_graph.addEdge(2, 1);
        ring_graph.addEdge(2, 3);
        ring_graph.addEdge(3, 2);
        ring_graph.addEdge(3, 0);
        ring_graph.addEdge(0, 3);
        
        // Create a complete graph on 4 vertices
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i != j) {
                    complete_graph.addEdge(i, j);
                }
            }
        }
    }
    
    UnidirectionalGraph tree_graph;
    UnidirectionalGraph hypercube_graph;
    UnidirectionalGraph ring_graph;
    UnidirectionalGraph complete_graph;
};

// Test reduction operations
TEST_F(GraphCollectiveTest, SumReduction) {
    GraphCollective::SumReduction<int> sum_op;
    
    EXPECT_EQ(sum_op(5, 3), 8);
    EXPECT_EQ(sum_op.identity(), 0);
    EXPECT_EQ(sum_op.name(), "SUM");
}

TEST_F(GraphCollectiveTest, MaxReduction) {
    GraphCollective::MaxReduction<int> max_op;
    
    EXPECT_EQ(max_op(5, 3), 5);
    EXPECT_EQ(max_op(3, 5), 5);
    EXPECT_EQ(max_op.identity(), std::numeric_limits<int>::lowest());
    EXPECT_EQ(max_op.name(), "MAX");
}

TEST_F(GraphCollectiveTest, MinReduction) {
    GraphCollective::MinReduction<int> min_op;
    
    EXPECT_EQ(min_op(5, 3), 3);
    EXPECT_EQ(min_op(3, 5), 3);
    EXPECT_EQ(min_op.identity(), std::numeric_limits<int>::max());
    EXPECT_EQ(min_op.name(), "MIN");
}

TEST_F(GraphCollectiveTest, AverageReduction) {
    GraphCollective::AverageReduction<double> avg_op(4);
    
    EXPECT_EQ(avg_op(5.0, 3.0), 8.0);
    EXPECT_EQ(avg_op.identity(), 0.0);
    EXPECT_EQ(avg_op.name(), "AVERAGE");
    EXPECT_EQ(avg_op.finalize(20.0), 5.0);
}

// Test validation functions
TEST_F(GraphCollectiveTest, ValidateInitialValues) {
    std::map<int, int> valid_values = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
    std::map<int, int> invalid_values = {{0, 1}, {1, 2}};  // Missing vertices
    
    EXPECT_NO_THROW(GraphCollective::validateInitialValues(tree_graph, valid_values));
    EXPECT_THROW(GraphCollective::validateInitialValues(tree_graph, invalid_values), 
                 std::invalid_argument);
    
    UnidirectionalGraph empty_graph;
    EXPECT_THROW(GraphCollective::validateInitialValues(empty_graph, valid_values), 
                 std::invalid_argument);
    
    std::map<int, int> empty_values;
    EXPECT_THROW(GraphCollective::validateInitialValues(tree_graph, empty_values), 
                 std::invalid_argument);
}

TEST_F(GraphCollectiveTest, ValidateTopologyForPattern) {
    // Tree pattern should work with any connected graph
    EXPECT_TRUE(GraphCollective::validateTopologyForPattern(tree_graph, 
                GraphCollective::CommunicationPattern::TREE));
    EXPECT_TRUE(GraphCollective::validateTopologyForPattern(ring_graph, 
                GraphCollective::CommunicationPattern::TREE));
    
    // Hypercube requires power-of-2 vertices
    EXPECT_TRUE(GraphCollective::validateTopologyForPattern(hypercube_graph, 
                GraphCollective::CommunicationPattern::HYPERCUBE));
    
    // Ring works with any connected graph
    EXPECT_TRUE(GraphCollective::validateTopologyForPattern(ring_graph, 
                GraphCollective::CommunicationPattern::RING));
    
    // 2D mesh requires perfect square vertices
    EXPECT_TRUE(GraphCollective::validateTopologyForPattern(hypercube_graph, 
                GraphCollective::CommunicationPattern::MESH_2D));
}

TEST_F(GraphCollectiveTest, AnalyzeOptimalPattern) {
    // Power-of-2 connected graph should prefer hypercube
    auto pattern = GraphCollective::analyzeOptimalPattern(hypercube_graph);
    EXPECT_EQ(pattern, GraphCollective::CommunicationPattern::HYPERCUBE);
    
    // Non-power-of-2 should default to tree
    UnidirectionalGraph triangle;
    triangle.addEdge(0, 1);
    triangle.addEdge(1, 2);
    triangle.addEdge(2, 0);
    triangle.addEdge(1, 0);
    triangle.addEdge(2, 1);
    triangle.addEdge(0, 2);
    
    pattern = GraphCollective::analyzeOptimalPattern(triangle);
    EXPECT_EQ(pattern, GraphCollective::CommunicationPattern::TREE);
}

// Test tree-based all-reduce
TEST_F(GraphCollectiveTest, AllReduceTreeSum) {
    std::map<int, int> initial_values = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
    GraphCollective::SumReduction<int> sum_op;
    
    auto result = GraphCollective::GraphCollective<int>::allReduceTree(
        tree_graph, initial_values, sum_op);
    
    EXPECT_EQ(result.final_value, 10);  // 1+2+3+4 = 10
    EXPECT_GT(result.total_rounds, 0);
    EXPECT_GT(result.total_messages, 0);
    EXPECT_GT(result.efficiency_metric, 0.0);
    EXPECT_LE(result.efficiency_metric, 1.0);
}

TEST_F(GraphCollectiveTest, AllReduceTreeMax) {
    std::map<int, int> initial_values = {{0, 1}, {1, 8}, {2, 3}, {3, 4}};
    GraphCollective::MaxReduction<int> max_op;
    
    auto result = GraphCollective::GraphCollective<int>::allReduceTree(
        tree_graph, initial_values, max_op);
    
    EXPECT_EQ(result.final_value, 8);
    EXPECT_GT(result.total_rounds, 0);
    EXPECT_GT(result.total_messages, 0);
}

TEST_F(GraphCollectiveTest, AllReduceTreeMin) {
    std::map<int, int> initial_values = {{0, 5}, {1, 2}, {2, 7}, {3, 4}};
    GraphCollective::MinReduction<int> min_op;
    
    auto result = GraphCollective::GraphCollective<int>::allReduceTree(
        tree_graph, initial_values, min_op);
    
    EXPECT_EQ(result.final_value, 2);
    EXPECT_GT(result.total_rounds, 0);
    EXPECT_GT(result.total_messages, 0);
}

// Test hypercube all-reduce
TEST_F(GraphCollectiveTest, AllReduceHypercubeSum) {
    std::map<int, int> initial_values = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
    GraphCollective::SumReduction<int> sum_op;
    
    auto result = GraphCollective::GraphCollective<int>::allReduceHypercube(
        hypercube_graph, initial_values, sum_op);
    
    EXPECT_EQ(result.final_value, 10);  // 1+2+3+4 = 10
    EXPECT_EQ(result.total_rounds, 2);  // log2(4) = 2 rounds
    EXPECT_EQ(result.total_messages, 8);  // 4 nodes * 2 rounds = 8 messages
    EXPECT_GT(result.efficiency_metric, 0.0);
}

TEST_F(GraphCollectiveTest, AllReduceHypercubeMax) {
    std::map<int, int> initial_values = {{0, 1}, {1, 8}, {2, 3}, {3, 4}};
    GraphCollective::MaxReduction<int> max_op;
    
    auto result = GraphCollective::GraphCollective<int>::allReduceHypercube(
        hypercube_graph, initial_values, max_op);
    
    EXPECT_EQ(result.final_value, 8);
    EXPECT_EQ(result.total_rounds, 2);
}

TEST_F(GraphCollectiveTest, AllReduceHypercubeInvalidTopology) {
    std::map<int, int> initial_values = {{0, 1}, {1, 2}, {2, 3}};  // 3 vertices (not power of 2)
    GraphCollective::SumReduction<int> sum_op;
    
    UnidirectionalGraph triangle;
    triangle.addEdge(0, 1);
    triangle.addEdge(1, 2);
    triangle.addEdge(2, 0);
    
    EXPECT_THROW(GraphCollective::GraphCollective<int>::allReduceHypercube(
        triangle, initial_values, sum_op), std::invalid_argument);
}

// Test ring all-reduce
TEST_F(GraphCollectiveTest, AllReduceRingSum) {
    std::map<int, int> initial_values = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
    GraphCollective::SumReduction<int> sum_op;
    
    auto result = GraphCollective::GraphCollective<int>::allReduceRing(
        ring_graph, initial_values, sum_op);
    
    EXPECT_EQ(result.final_value, 10);  // 1+2+3+4 = 10
    EXPECT_EQ(result.total_rounds, 2);  // Simplified ring: reduce + broadcast
    EXPECT_EQ(result.total_messages, 6);  // 3 reduce + 3 broadcast = 6 messages
    EXPECT_GT(result.efficiency_metric, 0.0);
}

TEST_F(GraphCollectiveTest, AllReduceRingSingleVertex) {
    std::map<int, int> initial_values = {{0, 5}};
    GraphCollective::SumReduction<int> sum_op;
    
    UnidirectionalGraph single_vertex;
    single_vertex.addVertex(0);
    
    auto result = GraphCollective::GraphCollective<int>::allReduceRing(
        single_vertex, initial_values, sum_op);
    
    EXPECT_EQ(result.final_value, 5);
    EXPECT_EQ(result.total_rounds, 0);
    EXPECT_EQ(result.total_messages, 0);
    EXPECT_EQ(result.efficiency_metric, 1.0);
}

// Test automatic pattern selection
TEST_F(GraphCollectiveTest, AllReduceAutoPattern) {
    std::map<int, int> initial_values = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
    GraphCollective::SumReduction<int> sum_op;
    
    // Should automatically choose hypercube for power-of-2 graph
    auto result = GraphCollective::GraphCollective<int>::allReduce(
        hypercube_graph, initial_values, sum_op, 
        GraphCollective::CommunicationPattern::OPTIMAL);
    
    EXPECT_EQ(result.final_value, 10);
    EXPECT_EQ(result.total_rounds, 2);  // Hypercube rounds
    
    // Should choose tree for non-power-of-2 graph
    UnidirectionalGraph triangle;
    triangle.addEdge(0, 1);
    triangle.addEdge(1, 0);
    triangle.addEdge(1, 2);
    triangle.addEdge(2, 1);
    triangle.addEdge(2, 0);
    triangle.addEdge(0, 2);
    
    std::map<int, int> triangle_values = {{0, 1}, {1, 2}, {2, 3}};
    result = GraphCollective::GraphCollective<int>::allReduce(
        triangle, triangle_values, sum_op, 
        GraphCollective::CommunicationPattern::OPTIMAL);
    
    EXPECT_EQ(result.final_value, 6);  // 1+2+3 = 6
}

// Test communication steps tracking
TEST_F(GraphCollectiveTest, CommunicationStepsTracking) {
    std::map<int, int> initial_values = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
    GraphCollective::SumReduction<int> sum_op;
    
    auto result = GraphCollective::GraphCollective<int>::allReduceTree(
        tree_graph, initial_values, sum_op);
    
    // Check that communication steps are recorded
    EXPECT_GT(result.communication_steps.size(), 0);
    
    // Check that steps have valid structure
    for (const auto& step : result.communication_steps) {
        EXPECT_GE(step.source_node, 0);
        EXPECT_GE(step.target_node, 0);
        EXPECT_GE(step.round, 1);
        EXPECT_FALSE(step.operation.empty());
    }
}

// Test efficiency metric calculation
TEST_F(GraphCollectiveTest, EfficiencyMetricCalculation) {
    // Test with different sizes
    GraphCollective::AllReduceResult<int> result1;
    result1.total_rounds = 2;
    result1.total_messages = 4;
    
    double efficiency1 = GraphCollective::calculateEfficiencyMetric(result1, 4);
    EXPECT_GT(efficiency1, 0.0);
    EXPECT_LE(efficiency1, 1.0);
    
    // More efficient result should have higher metric
    GraphCollective::AllReduceResult<int> result2;
    result2.total_rounds = 1;
    result2.total_messages = 2;
    
    double efficiency2 = GraphCollective::calculateEfficiencyMetric(result2, 4);
    EXPECT_GT(efficiency2, efficiency1);
}

// Test node state management
TEST_F(GraphCollectiveTest, NodeStateManagement) {
    GraphCollective::NodeState<int> node(1, 42);
    
    EXPECT_EQ(node.getNodeId(), 1);
    EXPECT_EQ(node.getValue(), 42);
    
    node.updateValue(100);
    EXPECT_EQ(node.getValue(), 100);
    
    node.addNeighbor(2);
    node.addNeighbor(3);
    
    auto neighbors = node.getNeighbors();
    EXPECT_EQ(neighbors.size(), 2);
    EXPECT_EQ(neighbors[0], 2);
    EXPECT_EQ(neighbors[1], 3);
}

// Test with different data types
TEST_F(GraphCollectiveTest, DifferentDataTypes) {
    // Test with doubles
    std::map<int, double> double_values = {{0, 1.5}, {1, 2.5}, {2, 3.5}, {3, 4.5}};
    GraphCollective::SumReduction<double> double_sum_op;
    
    auto double_result = GraphCollective::GraphCollective<double>::allReduceTree(
        tree_graph, double_values, double_sum_op);
    
    EXPECT_DOUBLE_EQ(double_result.final_value, 12.0);
    
    // Test with floats
    std::map<int, float> float_values = {{0, 1.0f}, {1, 2.0f}, {2, 3.0f}, {3, 4.0f}};
    GraphCollective::SumReduction<float> float_sum_op;
    
    auto float_result = GraphCollective::GraphCollective<float>::allReduceTree(
        tree_graph, float_values, float_sum_op);
    
    EXPECT_FLOAT_EQ(float_result.final_value, 10.0f);
}

// Test error handling
TEST_F(GraphCollectiveTest, ErrorHandling) {
    std::map<int, int> initial_values = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
    GraphCollective::SumReduction<int> sum_op;
    
    // Test unsupported pattern
    EXPECT_THROW(GraphCollective::GraphCollective<int>::allReduce(
        tree_graph, initial_values, sum_op, 
        GraphCollective::CommunicationPattern::BUTTERFLY), std::invalid_argument);
    
    // Test invalid topology for pattern
    EXPECT_THROW(GraphCollective::GraphCollective<int>::allReduce(
        tree_graph, initial_values, sum_op, 
        GraphCollective::CommunicationPattern::HYPERCUBE), std::invalid_argument);
}

// Performance comparison test
TEST_F(GraphCollectiveTest, PerformanceComparison) {
    std::map<int, int> initial_values = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
    GraphCollective::SumReduction<int> sum_op;
    
    auto tree_result = GraphCollective::GraphCollective<int>::allReduceTree(
        complete_graph, initial_values, sum_op);
    
    auto hypercube_result = GraphCollective::GraphCollective<int>::allReduceHypercube(
        hypercube_graph, initial_values, sum_op);
    
    auto ring_result = GraphCollective::GraphCollective<int>::allReduceRing(
        ring_graph, initial_values, sum_op);
    
    // All should produce same result
    EXPECT_EQ(tree_result.final_value, 10);
    EXPECT_EQ(hypercube_result.final_value, 10);
    EXPECT_EQ(ring_result.final_value, 10);
    
    // Hypercube should be most efficient for power-of-2 nodes
    EXPECT_LE(hypercube_result.total_rounds, tree_result.total_rounds);
    EXPECT_LE(hypercube_result.total_rounds, ring_result.total_rounds);
}
