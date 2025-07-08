#include <gtest/gtest.h>
#include "src/spanning_tree.h"
#include "src/graph_connectivity.h"
#include "src/hamiltonian.h"
#include "src/graph.h"
#include "src/tsp_boost.h"

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
    EXPECT_TRUE(GraphConnectivity::isConnected(connected_graph));
    
    // Test disconnected graph
    UnidirectionalGraph disconnected_graph;
    disconnected_graph.addEdge(0, 1);
    disconnected_graph.addEdge(2, 3);
    EXPECT_FALSE(GraphConnectivity::isConnected(disconnected_graph));
    
    // Test single vertex
    UnidirectionalGraph single_vertex;
    single_vertex.addVertex(0);
    EXPECT_TRUE(GraphConnectivity::isConnected(single_vertex));
    
    // Test empty graph
    UnidirectionalGraph empty_graph;
    EXPECT_TRUE(GraphConnectivity::isConnected(empty_graph));
}

TEST_F(SpanningTreeTest, ComputeSpanningTreeConnectedGraph) {
    UnidirectionalGraph spanning_tree = SpanningTree::computeSpanningTree(connected_graph);
    
    // A spanning tree should have n-1 edges where n is number of vertices
    EXPECT_EQ(spanning_tree.getEdgeCount(), 3);
    
    // Verify that the spanning tree is actually connected
    EXPECT_TRUE(GraphConnectivity::isConnected(spanning_tree));
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
    EXPECT_TRUE(GraphConnectivity::isConnected(spanning_forest));
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
    std::vector<int> vertices = GraphConnectivity::getAllVertices(connected_graph);
    
    // Should contain vertices 0, 1, 2, 3
    EXPECT_EQ(vertices.size(), 4);
    
    std::sort(vertices.begin(), vertices.end());
    std::vector<int> expected = {0, 1, 2, 3};
    EXPECT_EQ(vertices, expected);
}

TEST_F(SpanningTreeTest, FindHamiltonianCycleSquareGraph) {
    // Create a square graph with a valid Hamiltonian cycle
    UnidirectionalGraph square;
    square.addEdge(0, 1);
    square.addEdge(1, 2);
    square.addEdge(2, 3);
    square.addEdge(3, 0);
    
    std::vector<int> cycle = Hamiltonian::findExactHamiltonianCycle(square);
    
    // Should find a valid Hamiltonian cycle
    EXPECT_EQ(cycle.size(), 5); // 4 vertices + return to start
    EXPECT_EQ(cycle[0], cycle[4]); // Should start and end with same vertex
    
    // Verify it's a valid cycle (each edge exists)
    for (size_t i = 0; i < cycle.size() - 1; i++) {
        EXPECT_TRUE(square.hasEdge(cycle[i], cycle[i + 1]));
    }
    
    // Verify all vertices are visited exactly once (excluding the return)
    std::vector<int> visited(cycle.begin(), cycle.end() - 1);
    std::sort(visited.begin(), visited.end());
    std::vector<int> expected = {0, 1, 2, 3};
    EXPECT_EQ(visited, expected);
}

TEST_F(SpanningTreeTest, FindHamiltonianCycleTriangle) {
    UnidirectionalGraph triangle;
    triangle.addEdge(0, 1);
    triangle.addEdge(1, 2);
    triangle.addEdge(2, 0);
    
    std::vector<int> cycle = Hamiltonian::findExactHamiltonianCycle(triangle);
    
    // Should find a valid Hamiltonian cycle
    EXPECT_EQ(cycle.size(), 4); // 3 vertices + return to start
    EXPECT_EQ(cycle[0], cycle[3]); // Should start and end with same vertex
    
    // Verify it's a valid cycle
    for (size_t i = 0; i < cycle.size() - 1; i++) {
        EXPECT_TRUE(triangle.hasEdge(cycle[i], cycle[i + 1]));
    }
}

TEST_F(SpanningTreeTest, FindHamiltonianCycleNoSolution) {
    UnidirectionalGraph no_cycle;
    // Create a graph with no Hamiltonian cycle
    // 0 -> 1, 1 -> 2, but no way back to 0
    no_cycle.addEdge(0, 1);
    no_cycle.addEdge(1, 2);
    
    std::vector<int> cycle = Hamiltonian::findExactHamiltonianCycle(no_cycle);
    
    // Should return empty vector (no Hamiltonian cycle)
    EXPECT_TRUE(cycle.empty());
}

TEST_F(SpanningTreeTest, FindHamiltonianCycleOriginalConnectedGraph) {
    // Test with the original connected_graph from setUp
    // This graph doesn't have a Hamiltonian cycle
    std::vector<int> cycle = Hamiltonian::findExactHamiltonianCycle(connected_graph);
    
    // Should return empty vector (no Hamiltonian cycle)
    EXPECT_TRUE(cycle.empty());
}

TEST_F(SpanningTreeTest, FindHamiltonianCycleCompleteGraph) {
    UnidirectionalGraph complete;
    // Complete graph on 4 vertices
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i != j) {
                complete.addEdge(i, j);
            }
        }
    }
    
    std::vector<int> cycle = Hamiltonian::findExactHamiltonianCycle(complete);
    
    // Should find a valid Hamiltonian cycle
    EXPECT_EQ(cycle.size(), 5); // 4 vertices + return to start
    EXPECT_EQ(cycle[0], cycle[4]); // Should start and end with same vertex
    
    // Verify it's a valid cycle
    for (size_t i = 0; i < cycle.size() - 1; i++) {
        EXPECT_TRUE(complete.hasEdge(cycle[i], cycle[i + 1]));
    }
    
    // Verify all vertices are visited exactly once (excluding the return)
    std::vector<int> visited(cycle.begin(), cycle.end() - 1);
    std::sort(visited.begin(), visited.end());
    std::vector<int> expected = {0, 1, 2, 3};
    EXPECT_EQ(visited, expected);
}

TEST_F(SpanningTreeTest, FindHamiltonianCycleSingleVertex) {
    UnidirectionalGraph single;
    single.addVertex(0);
    
    std::vector<int> cycle = Hamiltonian::findExactHamiltonianCycle(single);
    
    // Single vertex cannot form a cycle
    EXPECT_TRUE(cycle.empty());
}

TEST_F(SpanningTreeTest, FindHamiltonianCycleTwoVertices) {
    UnidirectionalGraph two_vertices;
    two_vertices.addEdge(0, 1);
    two_vertices.addEdge(1, 0);
    
    std::vector<int> cycle = Hamiltonian::findExactHamiltonianCycle(two_vertices);
    
    // Two vertices cannot form a Hamiltonian cycle (need at least 3)
    EXPECT_TRUE(cycle.empty());
}

TEST_F(SpanningTreeTest, FindHamiltonianCycleEmptyGraph) {
    UnidirectionalGraph empty;
    
    std::vector<int> cycle = Hamiltonian::findExactHamiltonianCycle(empty);
    
    // Empty graph has no Hamiltonian cycle
    EXPECT_TRUE(cycle.empty());
}

TEST_F(SpanningTreeTest, FindHamiltonianCycleDisconnectedGraph) {
    UnidirectionalGraph disconnected;
    // Two disconnected triangles
    disconnected.addEdge(0, 1);
    disconnected.addEdge(1, 2);
    disconnected.addEdge(2, 0);
    
    disconnected.addEdge(3, 4);
    disconnected.addEdge(4, 5);
    disconnected.addEdge(5, 3);
    
    std::vector<int> cycle = Hamiltonian::findExactHamiltonianCycle(disconnected);
    
    // Disconnected graph cannot have a Hamiltonian cycle
    EXPECT_TRUE(cycle.empty());
}

TEST_F(SpanningTreeTest, FindApproximateHamiltonianCycleSquareGraph) {
    // Create a square graph with a valid Hamiltonian cycle
    UnidirectionalGraph square;
    square.addEdge(0, 1);
    square.addEdge(1, 2);
    square.addEdge(2, 3);
    square.addEdge(3, 0);
    
    std::vector<int> cycle = Hamiltonian::findApproximateHamiltonianCycle(square);
    
    // Should find a valid approximate Hamiltonian cycle
    EXPECT_EQ(cycle.size(), 5); // 4 vertices + return to start
    EXPECT_EQ(cycle[0], cycle[4]); // Should start and end with same vertex
    
    // Verify all vertices are visited exactly once (excluding the return)
    std::vector<int> visited(cycle.begin(), cycle.end() - 1);
    std::sort(visited.begin(), visited.end());
    std::vector<int> expected = {0, 1, 2, 3};
    EXPECT_EQ(visited, expected);
}

TEST_F(SpanningTreeTest, FindApproximateHamiltonianCycleCompleteGraph) {
    UnidirectionalGraph complete;
    // Complete graph on 4 vertices
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i != j) {
                complete.addEdge(i, j);
            }
        }
    }
    
    std::vector<int> cycle = Hamiltonian::findApproximateHamiltonianCycle(complete);
    
    // Should find a valid approximate Hamiltonian cycle
    EXPECT_EQ(cycle.size(), 5); // 4 vertices + return to start
    EXPECT_EQ(cycle[0], cycle[4]); // Should start and end with same vertex
    
    // Verify all vertices are visited exactly once (excluding the return)
    std::vector<int> visited(cycle.begin(), cycle.end() - 1);
    std::sort(visited.begin(), visited.end());
    std::vector<int> expected = {0, 1, 2, 3};
    EXPECT_EQ(visited, expected);
}

TEST_F(SpanningTreeTest, FindApproximateHamiltonianCycleWithWeights) {
    UnidirectionalGraph weighted_graph;
    weighted_graph.addEdge(0, 1);
    weighted_graph.addEdge(1, 2);
    weighted_graph.addEdge(2, 0);
    
    // Create weight matrix
    std::vector<std::vector<double>> weights(3, std::vector<double>(3, 1.0));
    weights[0][1] = 1.0;
    weights[1][2] = 1.0;
    weights[2][0] = 1.0;
    
    std::vector<int> cycle = Hamiltonian::findApproximateHamiltonianCycle(weighted_graph, weights);
    
    // Should find a valid approximate Hamiltonian cycle
    EXPECT_EQ(cycle.size(), 4); // 3 vertices + return to start
    EXPECT_EQ(cycle[0], cycle[3]); // Should start and end with same vertex
    
    // Verify all vertices are visited exactly once (excluding the return)
    std::vector<int> visited(cycle.begin(), cycle.end() - 1);
    std::sort(visited.begin(), visited.end());
    std::vector<int> expected = {0, 1, 2};
    EXPECT_EQ(visited, expected);
}

TEST_F(SpanningTreeTest, FindApproximateHamiltonianCycleEmptyGraph) {
    UnidirectionalGraph empty;
    
    std::vector<int> cycle = Hamiltonian::findApproximateHamiltonianCycle(empty);
    
    // Empty graph has no Hamiltonian cycle
    EXPECT_TRUE(cycle.empty());
}

TEST_F(SpanningTreeTest, FindApproximateHamiltonianCycleTwoVertices) {
    UnidirectionalGraph two_vertices;
    two_vertices.addEdge(0, 1);
    two_vertices.addEdge(1, 0);
    
    std::vector<int> cycle = Hamiltonian::findApproximateHamiltonianCycle(two_vertices);
    
    // Two vertices cannot form a Hamiltonian cycle (need at least 3)
    EXPECT_TRUE(cycle.empty());
}

TEST_F(SpanningTreeTest, BoostMSTKruskal) {
    UnidirectionalGraph test_graph;
    test_graph.addEdge(0, 1);
    test_graph.addEdge(0, 2);
    test_graph.addEdge(1, 2);
    test_graph.addEdge(1, 3);
    test_graph.addEdge(2, 3);
    
    UnidirectionalGraph mst = TSPBoost::computeMSTKruskal(test_graph);
    
    // MST should have n-1 edges where n is number of vertices
    EXPECT_EQ(mst.getVertexCount(), 4);
    EXPECT_EQ(mst.getEdgeCount(), 3);
    
    // Verify MST is connected
    EXPECT_TRUE(GraphConnectivity::isConnected(mst));
}

TEST_F(SpanningTreeTest, BoostMSTPrim) {
    UnidirectionalGraph test_graph;
    test_graph.addEdge(0, 1);
    test_graph.addEdge(0, 2);
    test_graph.addEdge(1, 2);
    test_graph.addEdge(1, 3);
    test_graph.addEdge(2, 3);
    
    UnidirectionalGraph mst = TSPBoost::computeMSTPrim(test_graph);
    
    // MST should have n-1 edges where n is number of vertices
    EXPECT_EQ(mst.getVertexCount(), 4);
    EXPECT_EQ(mst.getEdgeCount(), 3);
    
    // Verify MST is connected
    EXPECT_TRUE(GraphConnectivity::isConnected(mst));
}
