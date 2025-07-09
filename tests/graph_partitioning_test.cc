#include <gtest/gtest.h>
#include "src/graph_partitioning.h"
#include "src/graph.h"
#include <algorithm>

class GraphPartitioningTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple 4-vertex cycle for testing
        // 0 -- 1
        // |    |
        // 3 -- 2
        cycle_graph.addEdge(0, 1);
        cycle_graph.addEdge(1, 2);
        cycle_graph.addEdge(2, 3);
        cycle_graph.addEdge(3, 0);
    }

    UnidirectionalGraph cycle_graph;
};

TEST_F(GraphPartitioningTest, ValidateGraphForPartitioning) {
    UnidirectionalGraph empty_graph;
    EXPECT_THROW(GraphPartitioning::validateGraphForPartitioning(empty_graph), 
                 std::invalid_argument);
    
    UnidirectionalGraph single_vertex;
    single_vertex.addVertex(0);
    EXPECT_NO_THROW(GraphPartitioning::validateGraphForPartitioning(single_vertex));
    
    EXPECT_NO_THROW(GraphPartitioning::validateGraphForPartitioning(cycle_graph));
}

TEST_F(GraphPartitioningTest, CalculateCutSize) {
    std::vector<int> vertices = {0, 1, 2, 3};
    
    // Test partition {0, 1} vs {2, 3}
    // This should cut edges: 1->2 and 3->0
    int mask = 0b0011; // vertices 0 and 1 in partition 1
    int cut_size = GraphPartitioning::calculateCutSize(cycle_graph, vertices, mask);
    EXPECT_EQ(cut_size, 2);
    
    // Test partition {0, 2} vs {1, 3}  
    // This should cut edges: 0->1, 1->2, 2->3, 3->0
    mask = 0b0101; // vertices 0 and 2 in partition 1
    cut_size = GraphPartitioning::calculateCutSize(cycle_graph, vertices, mask);
    EXPECT_EQ(cut_size, 4);
}

TEST_F(GraphPartitioningTest, ComputeExactBisectionalBandwidthEmptyGraph) {
    UnidirectionalGraph empty_graph;
    EXPECT_THROW(GraphPartitioning::computeExactBisectionalBandwidth(empty_graph), 
                 std::invalid_argument);
}

TEST_F(GraphPartitioningTest, ComputeExactBisectionalBandwidthSingleVertex) {
    UnidirectionalGraph single_vertex;
    single_vertex.addVertex(0);
    
    int bandwidth = GraphPartitioning::computeExactBisectionalBandwidth(single_vertex);
    EXPECT_EQ(bandwidth, 0);
}

TEST_F(GraphPartitioningTest, ComputeExactBisectionalBandwidthTwoVertices) {
    UnidirectionalGraph two_vertices;
    two_vertices.addEdge(0, 1);
    
    int bandwidth = GraphPartitioning::computeExactBisectionalBandwidth(two_vertices);
    EXPECT_EQ(bandwidth, 1);
    
    // Test with bidirectional edge
    two_vertices.addEdge(1, 0);
    bandwidth = GraphPartitioning::computeExactBisectionalBandwidth(two_vertices);
    EXPECT_EQ(bandwidth, 2);
}

TEST_F(GraphPartitioningTest, ComputeExactBisectionalBandwidthTriangle) {
    UnidirectionalGraph triangle;
    triangle.addEdge(0, 1);
    triangle.addEdge(1, 2);
    triangle.addEdge(2, 0);
    
    int bandwidth = GraphPartitioning::computeExactBisectionalBandwidth(triangle);
    // Best partition for triangle: {0} vs {1, 2}
    // Cuts edges: 0->1 and 2->0, so bandwidth = 2
    EXPECT_EQ(bandwidth, 2);
}

TEST_F(GraphPartitioningTest, ComputeExactBisectionalBandwidthCycle) {
    // For 4-cycle, optimal partition is {0, 2} vs {1, 3}? No, that cuts 4 edges
    // Better partition is {0, 1} vs {2, 3}, which cuts 2 edges
    int bandwidth = GraphPartitioning::computeExactBisectionalBandwidth(cycle_graph);
    EXPECT_EQ(bandwidth, 2);
}

TEST_F(GraphPartitioningTest, ComputeExactBisectionalBandwidthCompleteGraph) {
    UnidirectionalGraph complete4;
    // Create complete graph on 4 vertices
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i != j) {
                complete4.addEdge(i, j);
            }
        }
    }
    
    int bandwidth = GraphPartitioning::computeExactBisectionalBandwidth(complete4);
    // For complete graph K4, best partition is {0,1} vs {2,3}
    // This cuts 2*2*2 = 8 edges (each vertex in one partition connects to each in the other, counting both directions)
    EXPECT_EQ(bandwidth, 8);
}

TEST_F(GraphPartitioningTest, ComputeExactBisectionalBandwidthLargeGraph) {
    UnidirectionalGraph large_graph;
    // Create a graph with > 20 vertices
    for (int i = 0; i < 25; i++) {
        large_graph.addVertex(i);
    }
    
    EXPECT_THROW(GraphPartitioning::computeExactBisectionalBandwidth(large_graph), 
                 std::runtime_error);
}

TEST_F(GraphPartitioningTest, FindBisectionalPartitionTwoVertices) {
    UnidirectionalGraph two_vertices;
    two_vertices.addEdge(0, 1);
    
    GraphPartitioning::GraphPartition partition = 
        GraphPartitioning::findBisectionalPartition(two_vertices);
    
    EXPECT_EQ(partition.partition1.size(), 1);
    EXPECT_EQ(partition.partition2.size(), 1);
    EXPECT_EQ(partition.bandwidth, 1);
    EXPECT_EQ(partition.cut_edges.size(), 1);
    EXPECT_TRUE(partition.isBalanced());
}

TEST_F(GraphPartitioningTest, FindBisectionalPartitionTriangle) {
    UnidirectionalGraph triangle;
    triangle.addEdge(0, 1);
    triangle.addEdge(1, 2);
    triangle.addEdge(2, 0);
    
    GraphPartitioning::GraphPartition partition = 
        GraphPartitioning::findBisectionalPartition(triangle);
    
    EXPECT_EQ(partition.bandwidth, 2);
    EXPECT_EQ(partition.cut_edges.size(), 2);
    EXPECT_TRUE(partition.isBalanced());
    
    // One partition should have 1 vertex, the other should have 2
    EXPECT_TRUE((partition.partition1.size() == 1 && partition.partition2.size() == 2) ||
                (partition.partition1.size() == 2 && partition.partition2.size() == 1));
}

TEST_F(GraphPartitioningTest, FindBisectionalPartitionCycle) {
    GraphPartitioning::GraphPartition partition = 
        GraphPartitioning::findBisectionalPartition(cycle_graph);
    
    EXPECT_EQ(partition.bandwidth, 2);
    EXPECT_EQ(partition.cut_edges.size(), 2);
    EXPECT_TRUE(partition.isBalanced());
    
    // Both partitions should have 2 vertices
    EXPECT_EQ(partition.partition1.size(), 2);
    EXPECT_EQ(partition.partition2.size(), 2);
}

TEST_F(GraphPartitioningTest, ComputeBisectionalBandwidth) {
    // Test automatic method selection
    int bandwidth = GraphPartitioning::computeBisectionalBandwidth(cycle_graph);
    EXPECT_EQ(bandwidth, 2);
    
    // Test single vertex
    UnidirectionalGraph single_vertex;
    single_vertex.addVertex(0);
    bandwidth = GraphPartitioning::computeBisectionalBandwidth(single_vertex);
    EXPECT_EQ(bandwidth, 0);
}

TEST_F(GraphPartitioningTest, ComputeKnownTopologyBandwidth) {
    // Test complete graph
    UnidirectionalGraph complete4;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i != j) {
                complete4.addEdge(i, j);
            }
        }
    }
    
    int bandwidth = GraphPartitioning::computeKnownTopologyBandwidth(
        complete4, GraphPartitioning::GraphTopology::COMPLETE_GRAPH);
    EXPECT_EQ(bandwidth, 4); // 4*4/4 = 4 (this is a theoretical approximation)
    
    // Test tree
    UnidirectionalGraph tree;
    tree.addEdge(0, 1);
    tree.addEdge(1, 2);
    tree.addEdge(2, 3);
    
    bandwidth = GraphPartitioning::computeKnownTopologyBandwidth(
        tree, GraphPartitioning::GraphTopology::TREE);
    EXPECT_EQ(bandwidth, 1);
}

TEST_F(GraphPartitioningTest, UnsupportedMethods) {
    // Test unsupported approximation methods
    EXPECT_THROW(GraphPartitioning::computeApproximateBisectionalBandwidth(
        cycle_graph, GraphPartitioning::PartitioningMethod::SPECTRAL), 
        std::runtime_error);
    
    EXPECT_THROW(GraphPartitioning::computeApproximateBisectionalBandwidth(
        cycle_graph, GraphPartitioning::PartitioningMethod::KERNIGHAN_LIN), 
        std::runtime_error);
}

TEST_F(GraphPartitioningTest, GraphPartitionIsBalanced) {
    GraphPartitioning::GraphPartition partition;
    
    // Balanced partitions
    partition.partition1 = {0, 1};
    partition.partition2 = {2, 3};
    EXPECT_TRUE(partition.isBalanced());
    
    partition.partition1 = {0};
    partition.partition2 = {1, 2};
    EXPECT_TRUE(partition.isBalanced());
    
    // Unbalanced partition
    partition.partition1 = {0};
    partition.partition2 = {1, 2, 3};
    EXPECT_FALSE(partition.isBalanced());
    
    // Empty partitions
    partition.partition1 = {};
    partition.partition2 = {0};
    EXPECT_TRUE(partition.isBalanced());
}

TEST_F(GraphPartitioningTest, DisconnectedGraph) {
    UnidirectionalGraph disconnected;
    disconnected.addEdge(0, 1);
    disconnected.addEdge(2, 3);
    
    // Should still work - bisectional bandwidth is about partitioning, not connectivity
    int bandwidth = GraphPartitioning::computeExactBisectionalBandwidth(disconnected);
    // Best partition: {0, 2} vs {1, 3} cuts 2 edges, or {0, 1} vs {2, 3} cuts 0 edges
    EXPECT_EQ(bandwidth, 0);
}

TEST_F(GraphPartitioningTest, PathGraph) {
    UnidirectionalGraph path;
    path.addEdge(0, 1);
    path.addEdge(1, 2);
    path.addEdge(2, 3);
    
    int bandwidth = GraphPartitioning::computeExactBisectionalBandwidth(path);
    // Best partition for path: {0, 1} vs {2, 3} cuts edge 1->2
    EXPECT_EQ(bandwidth, 1);
}
