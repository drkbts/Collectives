#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "graph.h"

// Test fixture for UnidirectionalGraph to reuse graph object
class UnidirectionalGraphTest : public ::testing::Test {
 protected:
  UnidirectionalGraph g;
};

// Test the constructor and initial state of an empty graph
TEST_F(UnidirectionalGraphTest, Constructor) {
  EXPECT_EQ(g.getVertexCount(), 0);
  EXPECT_EQ(g.getEdgeCount(), 0);
  EXPECT_FALSE(g.hasVertex(1));
}

// Test adding vertices
TEST_F(UnidirectionalGraphTest, AddVertex) {
  g.addVertex(1);
  EXPECT_TRUE(g.hasVertex(1));
  EXPECT_EQ(g.getVertexCount(), 1);
  EXPECT_EQ(g.getEdgeCount(), 0);

  // Adding an existing vertex should have no effect
  g.addVertex(1);
  EXPECT_TRUE(g.hasVertex(1));
  EXPECT_EQ(g.getVertexCount(), 1);

  // Add another vertex
  g.addVertex(2);
  EXPECT_TRUE(g.hasVertex(2));
  EXPECT_EQ(g.getVertexCount(), 2);
}

// Test hasVertex method
TEST_F(UnidirectionalGraphTest, HasVertex) {
  EXPECT_FALSE(g.hasVertex(1));
  g.addVertex(1);
  EXPECT_TRUE(g.hasVertex(1));
  EXPECT_FALSE(g.hasVertex(99));
}

// Test adding edges
TEST_F(UnidirectionalGraphTest, AddEdge) {
  g.addEdge(1, 2);
  EXPECT_TRUE(g.hasVertex(1));
  EXPECT_TRUE(g.hasVertex(2));
  EXPECT_TRUE(g.hasEdge(1, 2));
  EXPECT_FALSE(g.hasEdge(2, 1));  // It's a uni-directional graph
  EXPECT_EQ(g.getVertexCount(), 2);
  EXPECT_EQ(g.getEdgeCount(), 1);

  // Add another edge from the same source
  g.addEdge(1, 3);
  EXPECT_TRUE(g.hasVertex(3));
  EXPECT_TRUE(g.hasEdge(1, 3));
  EXPECT_EQ(g.getVertexCount(), 3);
  EXPECT_EQ(g.getEdgeCount(), 2);

  // Adding a duplicate edge should not change anything
  g.addEdge(1, 2);
  EXPECT_EQ(g.getVertexCount(), 3);
  EXPECT_EQ(g.getEdgeCount(), 2);
}

// Test hasEdge method
TEST_F(UnidirectionalGraphTest, HasEdge) {
  EXPECT_FALSE(g.hasEdge(1, 2));
  g.addEdge(1, 2);
  EXPECT_TRUE(g.hasEdge(1, 2));
  EXPECT_FALSE(g.hasEdge(2, 1));
  EXPECT_FALSE(g.hasEdge(1, 99));
  EXPECT_FALSE(g.hasEdge(99, 1));  // 'from' vertex doesn't exist
}

// Test getting neighbors of a vertex
TEST_F(UnidirectionalGraphTest, GetNeighbors) {
  g.addVertex(1);
  EXPECT_TRUE(g.getNeighbors(1).empty());

  g.addEdge(1, 2);
  g.addEdge(1, 3);
  g.addEdge(2, 3);

  const auto& neighbors1 = g.getNeighbors(1);
  // Order isn't guaranteed in the current implementation if it were to change.
  // So, we check for size and presence of each element.
  ASSERT_EQ(neighbors1.size(), 2);
  EXPECT_NE(std::find(neighbors1.begin(), neighbors1.end(), 2),
            neighbors1.end());
  EXPECT_NE(std::find(neighbors1.begin(), neighbors1.end(), 3),
            neighbors1.end());

  const auto& neighbors2 = g.getNeighbors(2);
  std::vector<int> expected2 = {3};
  EXPECT_EQ(neighbors2.size(), 1);
  EXPECT_NE(std::find(neighbors2.begin(), neighbors2.end(), 3),
            neighbors2.end());

  const auto& neighbors3 = g.getNeighbors(3);
  EXPECT_TRUE(neighbors3.empty());
}

// Test that getting neighbors for a non-existent vertex throws an exception
TEST_F(UnidirectionalGraphTest, GetNeighborsThrows) {
  EXPECT_THROW(g.getNeighbors(99), std::out_of_range);
  g.addVertex(1);
  EXPECT_NO_THROW(g.getNeighbors(1));
}

// Test getting all vertices
TEST_F(UnidirectionalGraphTest, GetVertices) {
  // Empty graph
  std::vector<int> vertices = g.getVertices();
  EXPECT_TRUE(vertices.empty());

  // Add vertices
  g.addVertex(1);
  g.addVertex(3);
  g.addVertex(2);
  
  vertices = g.getVertices();
  EXPECT_EQ(vertices.size(), 3);
  
  // Sort for consistent comparison
  std::sort(vertices.begin(), vertices.end());
  std::vector<int> expected = {1, 2, 3};
  EXPECT_EQ(vertices, expected);

  // Add vertex via edge
  g.addEdge(1, 5);
  vertices = g.getVertices();
  EXPECT_EQ(vertices.size(), 4);
  
  std::sort(vertices.begin(), vertices.end());
  expected = {1, 2, 3, 5};
  EXPECT_EQ(vertices, expected);
}

// Test vertex and edge counts during various operations
TEST_F(UnidirectionalGraphTest, Counts) {
  EXPECT_EQ(g.getVertexCount(), 0);
  EXPECT_EQ(g.getEdgeCount(), 0);

  g.addVertex(1);
  g.addVertex(2);
  EXPECT_EQ(g.getVertexCount(), 2);
  EXPECT_EQ(g.getEdgeCount(), 0);

  g.addEdge(1, 2);
  EXPECT_EQ(g.getVertexCount(), 2);
  EXPECT_EQ(g.getEdgeCount(), 1);

  g.addEdge(2, 3);  // Creates vertex 3
  EXPECT_EQ(g.getVertexCount(), 3);
  EXPECT_EQ(g.getEdgeCount(), 2);

  g.addEdge(1, 2);  // Duplicate edge
  EXPECT_EQ(g.getVertexCount(), 3);
  EXPECT_EQ(g.getEdgeCount(), 2);

  g.addVertex(1);  // Duplicate vertex
  EXPECT_EQ(g.getVertexCount(), 3);
  EXPECT_EQ(g.getEdgeCount(), 2);
}
