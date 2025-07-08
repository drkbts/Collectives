#include <gtest/gtest.h>
#include "src/spanning_tree.h"
#include "src/graph.h"

class SpanningTreeTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple connected graph for testing
        // 0 -- 1
        // |    |
        // 2 -- 3
        connected_graph.addEdge(0, 1);
        connected_graph.addEdge(0, 2);
        connected_graph.addEdge(1, 3);
        connected_graph.addEdge(2, 3);
    }

    UnidirectionalGraph connected_graph;
};

TEST_F(SpanningTreeTest, IsConnectedTest) {
    EXPECT_TRUE(SpanningTree::isConnected(connected_graph));
    
    // Test disconnected graph
    UnidirectionalGraph disconnected_graph;
    disconnected_graph.addEdge(0, 1);
    disconnected_graph.addEdge(2, 3);
    EXPECT_FALSE(SpanningTree::isConnected(disconnected_graph));
    
    // Test single vertex
    UnidirectionalGraph single_vertex;
    single_vertex.addVertex(0);
    EXPECT_TRUE(SpanningTree::isConnected(single_vertex));
    
    // Test empty graph
    UnidirectionalGraph empty_graph;
    EXPECT_TRUE(SpanningTree::isConnected(empty_graph));
}

TEST_F(SpanningTreeTest, ComputeSpanningTreeConnectedGraph) {
    UnidirectionalGraph spanning_tree = SpanningTree::computeSpanningTree(connected_graph);
    
    // A spanning tree should have n-1 edges where n is number of vertices
    EXPECT_EQ(spanning_tree.getEdgeCount(), 3);
    
    // Verify that the spanning tree is actually connected
    EXPECT_TRUE(SpanningTree::isConnected(spanning_tree));
}

TEST_F(SpanningTreeTest, ComputeSpanningTreeDisconnectedGraph) {
    UnidirectionalGraph disconnected_graph;
    disconnected_graph.addEdge(0, 1);
    disconnected_graph.addEdge(2, 3);
    
    // Should throw exception for disconnected graph
    EXPECT_THROW(SpanningTree::computeSpanningTree(disconnected_graph), std::runtime_error);
}

TEST_F(SpanningTreeTest, ComputeSpanningTreeSingleVertex) {
    UnidirectionalGraph single_vertex;
    single_vertex.addVertex(0);
    
    UnidirectionalGraph spanning_tree = SpanningTree::computeSpanningTree(single_vertex);
    EXPECT_EQ(spanning_tree.getEdgeCount(), 0);
}

TEST_F(SpanningTreeTest, ComputeSpanningTreeEmptyGraph) {
    UnidirectionalGraph empty_graph;
    
    // Should throw exception for empty graph
    EXPECT_THROW(SpanningTree::computeSpanningTree(empty_graph), std::invalid_argument);
}

TEST_F(SpanningTreeTest, ComputeSpanningForestConnectedGraph) {
    UnidirectionalGraph spanning_forest = SpanningTree::computeSpanningForest(connected_graph);
    
    // For a connected graph, spanning forest should be same as spanning tree
    EXPECT_EQ(spanning_forest.getEdgeCount(), 3);
    
    // Verify connectivity
    EXPECT_TRUE(SpanningTree::isConnected(spanning_forest));
}

TEST_F(SpanningTreeTest, ComputeSpanningForestDisconnectedGraph) {
    UnidirectionalGraph disconnected_graph;
    // Component 1: 0-1
    disconnected_graph.addEdge(0, 1);
    // Component 2: 2-3-4
    disconnected_graph.addEdge(2, 3);
    disconnected_graph.addEdge(3, 4);
    
    UnidirectionalGraph spanning_forest = SpanningTree::computeSpanningForest(disconnected_graph);
    
    // Should have 3 edges total (1 for first component, 2 for second component)
    EXPECT_EQ(spanning_forest.getEdgeCount(), 3);
    
    // Verify that each component is represented by checking connectivity
    EXPECT_TRUE(spanning_forest.hasVertex(0));
    EXPECT_TRUE(spanning_forest.hasVertex(1));
    EXPECT_TRUE(spanning_forest.hasVertex(2));
    EXPECT_TRUE(spanning_forest.hasVertex(3));
    EXPECT_TRUE(spanning_forest.hasVertex(4));
}

TEST_F(SpanningTreeTest, ComputeSpanningForestMultipleComponents) {
    UnidirectionalGraph multi_component;
    // Component 1: single vertex
    multi_component.addVertex(0);
    // Component 2: 1-2
    multi_component.addEdge(1, 2);
    // Component 3: 3-4-5 (triangle)
    multi_component.addEdge(3, 4);
    multi_component.addEdge(4, 5);
    multi_component.addEdge(5, 3);
    
    UnidirectionalGraph spanning_forest = SpanningTree::computeSpanningForest(multi_component);
    
    // Should have 3 edges total (0 for component 1, 1 for component 2, 2 for component 3)
    EXPECT_EQ(spanning_forest.getEdgeCount(), 3);
}

TEST_F(SpanningTreeTest, ComputeSpanningForestEmptyGraph) {
    UnidirectionalGraph empty_graph;
    
    // Should throw exception for empty graph
    EXPECT_THROW(SpanningTree::computeSpanningForest(empty_graph), std::invalid_argument);
}

TEST_F(SpanningTreeTest, GetAllVerticesTest) {
    std::vector<int> vertices = SpanningTree::getAllVertices(connected_graph);
    
    // Should contain vertices 0, 1, 2, 3
    EXPECT_EQ(vertices.size(), 4);
    
    std::sort(vertices.begin(), vertices.end());
    std::vector<int> expected = {0, 1, 2, 3};
    EXPECT_EQ(vertices, expected);
}
