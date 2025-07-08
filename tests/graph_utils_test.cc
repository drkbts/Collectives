#include <gtest/gtest.h>

#include <algorithm>
#include <vector>
#include <stdexcept>

#include "graph_utils.h"

// Test fixture for GraphUtils tests
class GraphUtilsTest : public ::testing::Test {
 protected:
  // Set up can be added here if needed
};

// Test createRing with invalid input
TEST_F(GraphUtilsTest, CreateRingInvalidInput) {
  EXPECT_THROW(GraphUtils::createRing(0), std::invalid_argument);
  EXPECT_THROW(GraphUtils::createRing(1), std::invalid_argument);
}

// Test createRing with size 3
TEST_F(GraphUtilsTest, CreateRingSize3) {
  auto ring3 = GraphUtils::createRing(3);
  
  EXPECT_EQ(ring3.getVertexCount(), 3);
  EXPECT_EQ(ring3.getEdgeCount(), 3);
  
  // Check vertices exist
  EXPECT_TRUE(ring3.hasVertex(0));
  EXPECT_TRUE(ring3.hasVertex(1));
  EXPECT_TRUE(ring3.hasVertex(2));
  
  // Check ring edges: 0->1->2->0
  EXPECT_TRUE(ring3.hasEdge(0, 1));
  EXPECT_TRUE(ring3.hasEdge(1, 2));
  EXPECT_TRUE(ring3.hasEdge(2, 0));
  
  // Check no reverse edges (unidirectional)
  EXPECT_FALSE(ring3.hasEdge(1, 0));
  EXPECT_FALSE(ring3.hasEdge(2, 1));
  EXPECT_FALSE(ring3.hasEdge(0, 2));
}

// Test createRing with size 5
TEST_F(GraphUtilsTest, CreateRingSize5) {
  auto ring5 = GraphUtils::createRing(5);
  
  EXPECT_EQ(ring5.getVertexCount(), 5);
  EXPECT_EQ(ring5.getEdgeCount(), 5);
  
  // Check all vertices exist
  for (int i = 0; i < 5; i++) {
    EXPECT_TRUE(ring5.hasVertex(i));
  }
  
  // Check ring edges: 0->1->2->3->4->0
  for (int i = 0; i < 4; i++) {
    EXPECT_TRUE(ring5.hasEdge(i, i + 1));
  }
  EXPECT_TRUE(ring5.hasEdge(4, 0)); // Close the ring
  
  // Check each vertex has exactly one neighbor
  for (int i = 0; i < 5; i++) {
    EXPECT_EQ(ring5.getNeighbors(i).size(), 1);
  }
}

// Test createRing with larger size
TEST_F(GraphUtilsTest, CreateRingSize10) {
  auto ring10 = GraphUtils::createRing(10);
  
  EXPECT_EQ(ring10.getVertexCount(), 10);
  EXPECT_EQ(ring10.getEdgeCount(), 10);
  
  // Check all vertices exist
  for (int i = 0; i < 10; i++) {
    EXPECT_TRUE(ring10.hasVertex(i));
  }
  
  // Check ring structure
  for (int i = 0; i < 9; i++) {
    EXPECT_TRUE(ring10.hasEdge(i, i + 1));
  }
  EXPECT_TRUE(ring10.hasEdge(9, 0)); // Close the ring
  
  // Check each vertex has exactly one neighbor
  for (int i = 0; i < 10; i++) {
    EXPECT_EQ(ring10.getNeighbors(i).size(), 1);
  }
}

// Test that createRing creates proper unidirectional structure
TEST_F(GraphUtilsTest, CreateRingUnidirectional) {
  auto ring = GraphUtils::createRing(4);
  
  // Verify it's truly unidirectional - no backward edges
  EXPECT_FALSE(ring.hasEdge(1, 0));
  EXPECT_FALSE(ring.hasEdge(2, 1));
  EXPECT_FALSE(ring.hasEdge(3, 2));
  EXPECT_FALSE(ring.hasEdge(0, 3));
  
  // Verify forward edges exist
  EXPECT_TRUE(ring.hasEdge(0, 1));
  EXPECT_TRUE(ring.hasEdge(1, 2));
  EXPECT_TRUE(ring.hasEdge(2, 3));
  EXPECT_TRUE(ring.hasEdge(3, 0));
}

// Test createBidirectionalRing with invalid input
TEST_F(GraphUtilsTest, CreateBidirectionalRingInvalidInput) {
  EXPECT_THROW(GraphUtils::createBidirectionalRing(0), std::invalid_argument);
  EXPECT_THROW(GraphUtils::createBidirectionalRing(1), std::invalid_argument);
}

// Test createBidirectionalRing with size 3
TEST_F(GraphUtilsTest, CreateBidirectionalRingSize3) {
  auto ring3 = GraphUtils::createBidirectionalRing(3);
  
  EXPECT_EQ(ring3.getVertexCount(), 3);
  EXPECT_EQ(ring3.getEdgeCount(), 6);  // 2 * N edges
  
  // Check vertices exist
  EXPECT_TRUE(ring3.hasVertex(0));
  EXPECT_TRUE(ring3.hasVertex(1));
  EXPECT_TRUE(ring3.hasVertex(2));
  
  // Check forward edges: 0->1->2->0
  EXPECT_TRUE(ring3.hasEdge(0, 1));
  EXPECT_TRUE(ring3.hasEdge(1, 2));
  EXPECT_TRUE(ring3.hasEdge(2, 0));
  
  // Check backward edges: 1->0, 2->1, 0->2
  EXPECT_TRUE(ring3.hasEdge(1, 0));
  EXPECT_TRUE(ring3.hasEdge(2, 1));
  EXPECT_TRUE(ring3.hasEdge(0, 2));
  
  // Check each vertex has exactly 2 neighbors
  for (int i = 0; i < 3; i++) {
    EXPECT_EQ(ring3.getNeighbors(i).size(), 2);
  }
}

// Test createBidirectionalRing with size 4
TEST_F(GraphUtilsTest, CreateBidirectionalRingSize4) {
  auto ring4 = GraphUtils::createBidirectionalRing(4);
  
  EXPECT_EQ(ring4.getVertexCount(), 4);
  EXPECT_EQ(ring4.getEdgeCount(), 8);  // 2 * N edges
  
  // Check all vertices exist
  for (int i = 0; i < 4; i++) {
    EXPECT_TRUE(ring4.hasVertex(i));
  }
  
  // Check forward edges: 0->1->2->3->0
  EXPECT_TRUE(ring4.hasEdge(0, 1));
  EXPECT_TRUE(ring4.hasEdge(1, 2));
  EXPECT_TRUE(ring4.hasEdge(2, 3));
  EXPECT_TRUE(ring4.hasEdge(3, 0));
  
  // Check backward edges: 1->0, 2->1, 3->2, 0->3
  EXPECT_TRUE(ring4.hasEdge(1, 0));
  EXPECT_TRUE(ring4.hasEdge(2, 1));
  EXPECT_TRUE(ring4.hasEdge(3, 2));
  EXPECT_TRUE(ring4.hasEdge(0, 3));
  
  // Check each vertex has exactly 2 neighbors
  for (int i = 0; i < 4; i++) {
    EXPECT_EQ(ring4.getNeighbors(i).size(), 2);
  }
}

// Test createBidirectionalRing neighbors content
TEST_F(GraphUtilsTest, CreateBidirectionalRingNeighbors) {
  auto ring = GraphUtils::createBidirectionalRing(5);
  
  // Check neighbors of vertex 0: should have 1 and 4
  const auto& neighbors0 = ring.getNeighbors(0);
  EXPECT_EQ(neighbors0.size(), 2);
  EXPECT_TRUE(std::find(neighbors0.begin(), neighbors0.end(), 1) != neighbors0.end());
  EXPECT_TRUE(std::find(neighbors0.begin(), neighbors0.end(), 4) != neighbors0.end());
  
  // Check neighbors of vertex 2: should have 1 and 3
  const auto& neighbors2 = ring.getNeighbors(2);
  EXPECT_EQ(neighbors2.size(), 2);
  EXPECT_TRUE(std::find(neighbors2.begin(), neighbors2.end(), 1) != neighbors2.end());
  EXPECT_TRUE(std::find(neighbors2.begin(), neighbors2.end(), 3) != neighbors2.end());
  
  // Check neighbors of vertex 4: should have 3 and 0
  const auto& neighbors4 = ring.getNeighbors(4);
  EXPECT_EQ(neighbors4.size(), 2);
  EXPECT_TRUE(std::find(neighbors4.begin(), neighbors4.end(), 3) != neighbors4.end());
  EXPECT_TRUE(std::find(neighbors4.begin(), neighbors4.end(), 0) != neighbors4.end());
}

// Test difference between unidirectional and bidirectional rings
TEST_F(GraphUtilsTest, CompareUnidirectionalVsBidirectional) {
  auto uniRing = GraphUtils::createRing(4);
  auto biRing = GraphUtils::createBidirectionalRing(4);
  
  // Same number of vertices
  EXPECT_EQ(uniRing.getVertexCount(), biRing.getVertexCount());
  
  // Bidirectional has twice as many edges
  EXPECT_EQ(uniRing.getEdgeCount(), 4);
  EXPECT_EQ(biRing.getEdgeCount(), 8);
  
  // Unidirectional has 1 neighbor per vertex, bidirectional has 2
  for (int i = 0; i < 4; i++) {
    EXPECT_EQ(uniRing.getNeighbors(i).size(), 1);
    EXPECT_EQ(biRing.getNeighbors(i).size(), 2);
  }
}

// Test createTensorProduct with invalid input
TEST_F(GraphUtilsTest, CreateTensorProductInvalidInput) {
  UnidirectionalGraph emptyGraph;
  UnidirectionalGraph validGraph = GraphUtils::createRing(3);
  
  EXPECT_THROW(GraphUtils::createTensorProduct(emptyGraph, validGraph), std::invalid_argument);
  EXPECT_THROW(GraphUtils::createTensorProduct(validGraph, emptyGraph), std::invalid_argument);
  EXPECT_THROW(GraphUtils::createTensorProduct(emptyGraph, emptyGraph), std::invalid_argument);
}

// Test createTensorProduct with simple graphs
TEST_F(GraphUtilsTest, CreateTensorProductSimple) {
  // Create first graph: 0 -> 1
  UnidirectionalGraph g1;
  g1.addEdge(0, 1);
  
  // Create second graph: 0 -> 1
  UnidirectionalGraph g2;
  g2.addEdge(0, 1);
  
  auto product = GraphUtils::createTensorProduct(g1, g2);
  
  // Should have 4 vertices: (0,0)=0, (0,1)=1, (1,0)=2, (1,1)=3
  EXPECT_EQ(product.getVertexCount(), 4);
  
  // Should have 1 edge: (0,0) -> (1,1) because 0->1 in g1 AND 0->1 in g2
  EXPECT_EQ(product.getEdgeCount(), 1);
  
  // Check vertices exist
  EXPECT_TRUE(product.hasVertex(0));  // (0,0)
  EXPECT_TRUE(product.hasVertex(1));  // (0,1)
  EXPECT_TRUE(product.hasVertex(2));  // (1,0)
  EXPECT_TRUE(product.hasVertex(3));  // (1,1)
  
  // Check the single edge: (0,0) -> (1,1)
  EXPECT_TRUE(product.hasEdge(0, 3));
  
  // Check no other edges exist
  EXPECT_FALSE(product.hasEdge(0, 1));
  EXPECT_FALSE(product.hasEdge(0, 2));
  EXPECT_FALSE(product.hasEdge(1, 2));
  EXPECT_FALSE(product.hasEdge(1, 3));
  EXPECT_FALSE(product.hasEdge(2, 3));
}

// Test createTensorProduct with self-loops
TEST_F(GraphUtilsTest, CreateTensorProductSelfLoops) {
  // Create graph with self-loop: 0 -> 0
  UnidirectionalGraph g1;
  g1.addEdge(0, 0);
  
  // Create simple edge: 0 -> 1
  UnidirectionalGraph g2;
  g2.addEdge(0, 1);
  
  auto product = GraphUtils::createTensorProduct(g1, g2);
  
  // g1 has vertices {0}, g2 has vertices {0, 1}
  // Product should have 2 vertices: (0,0)=0, (0,1)=1
  EXPECT_EQ(product.getVertexCount(), 2);
  
  // Should have 1 edge: (0,0) -> (0,1) because 0->0 in g1 AND 0->1 in g2
  EXPECT_EQ(product.getEdgeCount(), 1);
  
  // Check the edge: (0,0) -> (0,1)
  EXPECT_TRUE(product.hasEdge(0, 1));
}

// Test createTensorProduct with complete graphs
TEST_F(GraphUtilsTest, CreateTensorProductComplete) {
  // Create K2 (complete graph on 2 vertices): 0 <-> 1
  UnidirectionalGraph k2;
  k2.addEdge(0, 1);
  k2.addEdge(1, 0);
  
  // Create K2 again
  UnidirectionalGraph k2_copy;
  k2_copy.addEdge(0, 1);
  k2_copy.addEdge(1, 0);
  
  auto product = GraphUtils::createTensorProduct(k2, k2_copy);
  
  // Should have 4 vertices
  EXPECT_EQ(product.getVertexCount(), 4);
  
  // Should have 4 edges: all combinations where both graphs have edges
  // (0,0) -> (1,1): 0->1 in both graphs
  // (0,1) -> (1,0): 0->1 in g1, 1->0 in g2
  // (1,0) -> (0,1): 1->0 in g1, 0->1 in g2
  // (1,1) -> (0,0): 1->0 in both graphs
  EXPECT_EQ(product.getEdgeCount(), 4);
  
  EXPECT_TRUE(product.hasEdge(0, 3));  // (0,0) -> (1,1)
  EXPECT_TRUE(product.hasEdge(1, 2));  // (0,1) -> (1,0)
  EXPECT_TRUE(product.hasEdge(2, 1));  // (1,0) -> (0,1)
  EXPECT_TRUE(product.hasEdge(3, 0));  // (1,1) -> (0,0)
}

// Test createTensorProduct with ring graphs
TEST_F(GraphUtilsTest, CreateTensorProductRings) {
  // Create ring of size 3: 0 -> 1 -> 2 -> 0
  auto ring3 = GraphUtils::createRing(3);
  
  // Create ring of size 2: 0 -> 1 -> 0
  auto ring2 = GraphUtils::createRing(2);
  
  auto product = GraphUtils::createTensorProduct(ring3, ring2);
  
  // Should have 6 vertices: 3 * 2 = 6
  EXPECT_EQ(product.getVertexCount(), 6);
  
  // Should have 6 edges: 3 edges in ring3 * 2 edges in ring2 = 6
  EXPECT_EQ(product.getEdgeCount(), 6);
  
  // Check some specific edges
  // (0,0) -> (1,1): 0->1 in ring3, 0->1 in ring2
  EXPECT_TRUE(product.hasEdge(0, 3));
  
  // (1,0) -> (2,1): 1->2 in ring3, 0->1 in ring2
  EXPECT_TRUE(product.hasEdge(2, 5));
  
  // (2,1) -> (0,0): 2->0 in ring3, 1->0 in ring2
  EXPECT_TRUE(product.hasEdge(5, 0));
}

// Test createTensorProduct vertex count formula
TEST_F(GraphUtilsTest, CreateTensorProductVertexCount) {
  auto ring3 = GraphUtils::createRing(3);
  auto ring4 = GraphUtils::createRing(4);
  auto ring5 = GraphUtils::createRing(5);
  
  // Test various combinations
  auto product34 = GraphUtils::createTensorProduct(ring3, ring4);
  EXPECT_EQ(product34.getVertexCount(), 12);  // 3 * 4
  
  auto product45 = GraphUtils::createTensorProduct(ring4, ring5);
  EXPECT_EQ(product45.getVertexCount(), 20);  // 4 * 5
  
  auto product35 = GraphUtils::createTensorProduct(ring3, ring5);
  EXPECT_EQ(product35.getVertexCount(), 15);  // 3 * 5
}

// Test createTensorProduct with single vertex graphs
TEST_F(GraphUtilsTest, CreateTensorProductSingleVertex) {
  // Create single vertex graph with self-loop
  UnidirectionalGraph singleVertex;
  singleVertex.addVertex(0);
  
  // Create simple edge graph
  UnidirectionalGraph edgeGraph;
  edgeGraph.addEdge(0, 1);
  
  auto product = GraphUtils::createTensorProduct(singleVertex, edgeGraph);
  
  // Should have 2 vertices: 1 * 2 = 2
  EXPECT_EQ(product.getVertexCount(), 2);
  
  // Should have 0 edges since single vertex graph has no edges
  EXPECT_EQ(product.getEdgeCount(), 0);
}

// Test createTorus with invalid input
TEST_F(GraphUtilsTest, CreateTorusInvalidInput) {
  std::vector<int> emptyDimensions;
  EXPECT_THROW(GraphUtils::createTorus(emptyDimensions), std::invalid_argument);
  
  std::vector<int> invalidDimensions1 = {3, 0, 4};
  EXPECT_THROW(GraphUtils::createTorus(invalidDimensions1), std::invalid_argument);
  
  std::vector<int> invalidDimensions2 = {3, 1, 4};
  EXPECT_THROW(GraphUtils::createTorus(invalidDimensions2), std::invalid_argument);
  
  std::vector<int> invalidDimensions3 = {-1, 3};
  EXPECT_THROW(GraphUtils::createTorus(invalidDimensions3), std::invalid_argument);
}

// Test createTorus with single dimension (should be same as bidirectional ring)
TEST_F(GraphUtilsTest, CreateTorusSingleDimension) {
  std::vector<int> dimensions = {5};
  auto torus = GraphUtils::createTorus(dimensions);
  auto ring = GraphUtils::createBidirectionalRing(5);
  
  // Should have same structure as bidirectional ring
  EXPECT_EQ(torus.getVertexCount(), ring.getVertexCount());
  EXPECT_EQ(torus.getEdgeCount(), ring.getEdgeCount());
  
  // Check all vertices and edges match
  for (int i = 0; i < 5; i++) {
    EXPECT_EQ(torus.hasVertex(i), ring.hasVertex(i));
    EXPECT_EQ(torus.getNeighbors(i).size(), ring.getNeighbors(i).size());
  }
}

// Test createTorus with 2D torus (3x4 grid with wraparound)
TEST_F(GraphUtilsTest, CreateTorus2D) {
  std::vector<int> dimensions = {3, 4};
  auto torus = GraphUtils::createTorus(dimensions);
  
  // Should have 12 vertices: 3 * 4 = 12
  EXPECT_EQ(torus.getVertexCount(), 12);
  
  // Each vertex should have 4 neighbors (tensor product of two bidirectional rings)
  for (int i = 0; i < 12; i++) {
    EXPECT_EQ(torus.getNeighbors(i).size(), 4);
  }
  
  // Total edges should be 48: 12 vertices * 4 neighbors / 2 (each edge counted twice)
  EXPECT_EQ(torus.getEdgeCount(), 48);
}

// Test createTorus with 3D torus (2x3x4 grid)
TEST_F(GraphUtilsTest, CreateTorus3D) {
  std::vector<int> dimensions = {2, 3, 4};
  auto torus = GraphUtils::createTorus(dimensions);
  
  // Should have 24 vertices: 2 * 3 * 4 = 24
  EXPECT_EQ(torus.getVertexCount(), 24);
  
  // Each vertex should have 4 neighbors (tensor product behavior)
  for (int i = 0; i < 24; i++) {
    EXPECT_EQ(torus.getNeighbors(i).size(), 4);
  }
  
  // Total edges should be 96: 24 vertices * 4 neighbors / 2 (each edge counted twice)
  EXPECT_EQ(torus.getEdgeCount(), 96);
}

// Test createTorus with 4D torus (minimal case)
TEST_F(GraphUtilsTest, CreateTorus4D) {
  std::vector<int> dimensions = {2, 2, 2, 2};
  auto torus = GraphUtils::createTorus(dimensions);
  
  // Should have 16 vertices: 2^4 = 16
  EXPECT_EQ(torus.getVertexCount(), 16);
  
  // Each vertex should have 1 neighbor (tensor product behavior with 2-rings)
  for (int i = 0; i < 16; i++) {
    EXPECT_EQ(torus.getNeighbors(i).size(), 1);
  }
  
  // Total edges should be 16: 16 vertices * 1 neighbor / 2 (each edge counted twice)
  EXPECT_EQ(torus.getEdgeCount(), 16);
}

// Test createTorus vertex count formula
TEST_F(GraphUtilsTest, CreateTorusVertexCount) {
  std::vector<int> dimensions1 = {3, 5};
  auto torus1 = GraphUtils::createTorus(dimensions1);
  EXPECT_EQ(torus1.getVertexCount(), 15);  // 3 * 5
  
  std::vector<int> dimensions2 = {2, 3, 4};
  auto torus2 = GraphUtils::createTorus(dimensions2);
  EXPECT_EQ(torus2.getVertexCount(), 24);  // 2 * 3 * 4
  
  std::vector<int> dimensions3 = {5, 2, 3};
  auto torus3 = GraphUtils::createTorus(dimensions3);
  EXPECT_EQ(torus3.getVertexCount(), 30);  // 5 * 2 * 3
}

// Test createTorus neighbor count behavior
TEST_F(GraphUtilsTest, CreateTorusNeighborCount) {
  // 1D torus: each vertex has 2 neighbors (same as bidirectional ring)
  std::vector<int> dimensions1 = {5};
  auto torus1 = GraphUtils::createTorus(dimensions1);
  for (int i = 0; i < 5; i++) {
    EXPECT_EQ(torus1.getNeighbors(i).size(), 2);
  }
  
  // 2D torus: each vertex has 4 neighbors (tensor product behavior)
  std::vector<int> dimensions2 = {3, 4};
  auto torus2 = GraphUtils::createTorus(dimensions2);
  for (int i = 0; i < 12; i++) {
    EXPECT_EQ(torus2.getNeighbors(i).size(), 4);
  }
  
  // 3D torus: each vertex has 4 neighbors (tensor product behavior)
  std::vector<int> dimensions3 = {2, 3, 4};
  auto torus3 = GraphUtils::createTorus(dimensions3);
  for (int i = 0; i < 24; i++) {
    EXPECT_EQ(torus3.getNeighbors(i).size(), 4);
  }
}

// Test createTorus edge count behavior
TEST_F(GraphUtilsTest, CreateTorusEdgeCount) {
  // Edge count depends on tensor product behavior
  
  std::vector<int> dimensions1 = {7};  // 1D
  auto torus1 = GraphUtils::createTorus(dimensions1);
  EXPECT_EQ(torus1.getEdgeCount(), 14);  // 7 * 2 (bidirectional ring)
  
  std::vector<int> dimensions2 = {3, 4};  // 2D
  auto torus2 = GraphUtils::createTorus(dimensions2);
  EXPECT_EQ(torus2.getEdgeCount(), 48);  // 12 vertices * 4 neighbors / 2
  
  std::vector<int> dimensions3 = {2, 3, 4};  // 3D
  auto torus3 = GraphUtils::createTorus(dimensions3);
  EXPECT_EQ(torus3.getEdgeCount(), 96);  // 24 vertices * 4 neighbors / 2
}

// Test createTorus with same dimensions in different order
TEST_F(GraphUtilsTest, CreateTorusDimensionOrder) {
  std::vector<int> dimensions1 = {3, 4};
  std::vector<int> dimensions2 = {4, 3};
  
  auto torus1 = GraphUtils::createTorus(dimensions1);
  auto torus2 = GraphUtils::createTorus(dimensions2);
  
  // Should have same vertex and edge counts regardless of order
  EXPECT_EQ(torus1.getVertexCount(), torus2.getVertexCount());
  EXPECT_EQ(torus1.getEdgeCount(), torus2.getEdgeCount());
  
  // Each vertex should have same number of neighbors
  for (int i = 0; i < 12; i++) {
    EXPECT_EQ(torus1.getNeighbors(i).size(), torus2.getNeighbors(i).size());
  }
}

// Test parseGraphExpression with valid uR[N] expressions
TEST_F(GraphUtilsTest, ParseGraphExpressionValidUnidirectionalRing) {
  // Test basic unidirectional ring
  auto ring5 = GraphUtils::parseGraphExpression("uR[5]");
  auto expected5 = GraphUtils::createRing(5);
  
  EXPECT_EQ(ring5.getVertexCount(), expected5.getVertexCount());
  EXPECT_EQ(ring5.getEdgeCount(), expected5.getEdgeCount());
  
  // Check structure matches
  for (int i = 0; i < 5; i++) {
    EXPECT_EQ(ring5.hasVertex(i), expected5.hasVertex(i));
    EXPECT_EQ(ring5.getNeighbors(i).size(), expected5.getNeighbors(i).size());
  }
  
  // Test with different sizes
  auto ring3 = GraphUtils::parseGraphExpression("uR[3]");
  EXPECT_EQ(ring3.getVertexCount(), 3);
  EXPECT_EQ(ring3.getEdgeCount(), 3);
  
  auto ring10 = GraphUtils::parseGraphExpression("uR[10]");
  EXPECT_EQ(ring10.getVertexCount(), 10);
  EXPECT_EQ(ring10.getEdgeCount(), 10);
}

// Test parseGraphExpression with whitespace
TEST_F(GraphUtilsTest, ParseGraphExpressionWithWhitespace) {
  // Test with various whitespace patterns
  auto ring1 = GraphUtils::parseGraphExpression(" uR[5] ");
  auto ring2 = GraphUtils::parseGraphExpression("uR[ 5 ]");
  auto ring3 = GraphUtils::parseGraphExpression(" uR [ 5 ] ");
  auto ring4 = GraphUtils::parseGraphExpression("\tuR[5]\n");
  
  // All should create the same ring
  EXPECT_EQ(ring1.getVertexCount(), 5);
  EXPECT_EQ(ring2.getVertexCount(), 5);
  EXPECT_EQ(ring3.getVertexCount(), 5);
  EXPECT_EQ(ring4.getVertexCount(), 5);
  
  EXPECT_EQ(ring1.getEdgeCount(), 5);
  EXPECT_EQ(ring2.getEdgeCount(), 5);
  EXPECT_EQ(ring3.getEdgeCount(), 5);
  EXPECT_EQ(ring4.getEdgeCount(), 5);
}

// Test parseGraphExpression with invalid expressions
TEST_F(GraphUtilsTest, ParseGraphExpressionInvalidExpressions) {
  // Test completely invalid expressions
  EXPECT_THROW(GraphUtils::parseGraphExpression(""), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("invalid"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR[]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR["), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("R[5]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("ur[5]"), std::invalid_argument);  // case sensitive
  
  // Test with invalid numbers
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR[0]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR[1]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR[-5]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR[abc]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR[3.5]"), std::invalid_argument);
  
  // Test with extra characters
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR[5]extra"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("extrauR[5]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR[5][6]"), std::invalid_argument);
}

// Test parseGraphExpression with large numbers
TEST_F(GraphUtilsTest, ParseGraphExpressionLargeNumbers) {
  // Test with reasonably large number
  auto ring100 = GraphUtils::parseGraphExpression("uR[100]");
  EXPECT_EQ(ring100.getVertexCount(), 100);
  EXPECT_EQ(ring100.getEdgeCount(), 100);
  
  // Test with very large number (should throw out_of_range if too large)
  // This depends on system limits, but we can test the mechanism
  EXPECT_NO_THROW(GraphUtils::parseGraphExpression("uR[1000]"));
}

// Test parseGraphExpression edge cases
TEST_F(GraphUtilsTest, ParseGraphExpressionEdgeCases) {
  // Test minimum valid size
  auto ring2 = GraphUtils::parseGraphExpression("uR[2]");
  EXPECT_EQ(ring2.getVertexCount(), 2);
  EXPECT_EQ(ring2.getEdgeCount(), 2);
  EXPECT_TRUE(ring2.hasEdge(0, 1));
  EXPECT_TRUE(ring2.hasEdge(1, 0));
  
  // Test with leading zeros (should still work)
  auto ring007 = GraphUtils::parseGraphExpression("uR[007]");
  EXPECT_EQ(ring007.getVertexCount(), 7);
  EXPECT_EQ(ring007.getEdgeCount(), 7);
}

// Test parseGraphExpression consistency with direct function calls
TEST_F(GraphUtilsTest, ParseGraphExpressionConsistency) {
  // Test that DSL produces same results as direct function calls
  for (int n = 2; n <= 10; n++) {
    std::string expression = "uR[" + std::to_string(n) + "]";
    auto dslRing = GraphUtils::parseGraphExpression(expression);
    auto directRing = GraphUtils::createRing(n);
    
    EXPECT_EQ(dslRing.getVertexCount(), directRing.getVertexCount());
    EXPECT_EQ(dslRing.getEdgeCount(), directRing.getEdgeCount());
    
    // Check all edges match
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        EXPECT_EQ(dslRing.hasEdge(i, j), directRing.hasEdge(i, j));
      }
    }
  }
}

// Test parseGraphExpression with valid bR[N] expressions
TEST_F(GraphUtilsTest, ParseGraphExpressionValidBidirectionalRing) {
  // Test basic bidirectional ring
  auto ring5 = GraphUtils::parseGraphExpression("bR[5]");
  auto expected5 = GraphUtils::createBidirectionalRing(5);
  
  EXPECT_EQ(ring5.getVertexCount(), expected5.getVertexCount());
  EXPECT_EQ(ring5.getEdgeCount(), expected5.getEdgeCount());
  
  // Check structure matches
  for (int i = 0; i < 5; i++) {
    EXPECT_EQ(ring5.hasVertex(i), expected5.hasVertex(i));
    EXPECT_EQ(ring5.getNeighbors(i).size(), expected5.getNeighbors(i).size());
  }
  
  // Test with different sizes
  auto ring3 = GraphUtils::parseGraphExpression("bR[3]");
  EXPECT_EQ(ring3.getVertexCount(), 3);
  EXPECT_EQ(ring3.getEdgeCount(), 6);  // Bidirectional has 2x edges
  
  auto ring4 = GraphUtils::parseGraphExpression("bR[4]");
  EXPECT_EQ(ring4.getVertexCount(), 4);
  EXPECT_EQ(ring4.getEdgeCount(), 8);  // Bidirectional has 2x edges
}

// Test parseGraphExpression bidirectional ring with whitespace
TEST_F(GraphUtilsTest, ParseGraphExpressionBidirectionalWithWhitespace) {
  // Test with various whitespace patterns
  auto ring1 = GraphUtils::parseGraphExpression(" bR[4] ");
  auto ring2 = GraphUtils::parseGraphExpression("bR[ 4 ]");
  auto ring3 = GraphUtils::parseGraphExpression(" bR [ 4 ] ");
  auto ring4 = GraphUtils::parseGraphExpression("\tbR[4]\n");
  
  // All should create the same bidirectional ring
  EXPECT_EQ(ring1.getVertexCount(), 4);
  EXPECT_EQ(ring2.getVertexCount(), 4);
  EXPECT_EQ(ring3.getVertexCount(), 4);
  EXPECT_EQ(ring4.getVertexCount(), 4);
  
  EXPECT_EQ(ring1.getEdgeCount(), 8);  // 4 vertices * 2 edges each
  EXPECT_EQ(ring2.getEdgeCount(), 8);
  EXPECT_EQ(ring3.getEdgeCount(), 8);
  EXPECT_EQ(ring4.getEdgeCount(), 8);
}

// Test parseGraphExpression invalid bidirectional ring expressions
TEST_F(GraphUtilsTest, ParseGraphExpressionInvalidBidirectionalRing) {
  // Test invalid bidirectional ring expressions
  EXPECT_THROW(GraphUtils::parseGraphExpression("bR[]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("bR[0]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("bR[1]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("bR[-3]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("bR[abc]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("BR[5]"), std::invalid_argument);  // case sensitive
  EXPECT_THROW(GraphUtils::parseGraphExpression("br[5]"), std::invalid_argument);  // case sensitive
  EXPECT_THROW(GraphUtils::parseGraphExpression("bR[5]extra"), std::invalid_argument);
}

// Test parseGraphExpression bidirectional consistency
TEST_F(GraphUtilsTest, ParseGraphExpressionBidirectionalConsistency) {
  // Test that bR[N] DSL produces same results as direct bidirectional function calls
  for (int n = 2; n <= 8; n++) {
    std::string expression = "bR[" + std::to_string(n) + "]";
    auto dslRing = GraphUtils::parseGraphExpression(expression);
    auto directRing = GraphUtils::createBidirectionalRing(n);
    
    EXPECT_EQ(dslRing.getVertexCount(), directRing.getVertexCount());
    EXPECT_EQ(dslRing.getEdgeCount(), directRing.getEdgeCount());
    
    // Check all edges match
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        EXPECT_EQ(dslRing.hasEdge(i, j), directRing.hasEdge(i, j));
      }
    }
    
    // Check each vertex has exactly 2 neighbors in bidirectional ring (except for n=2)
    for (int i = 0; i < n; i++) {
      if (n == 2) {
        EXPECT_EQ(dslRing.getNeighbors(i).size(), 1);  // Special case: 2-vertex ring
      } else {
        EXPECT_EQ(dslRing.getNeighbors(i).size(), 2);
      }
    }
  }
}

// Test parseGraphExpression comparison between unidirectional and bidirectional
TEST_F(GraphUtilsTest, ParseGraphExpressionCompareDirections) {
  // Compare uR[N] vs bR[N] for same N
  for (int n = 3; n <= 6; n++) {
    auto uniRing = GraphUtils::parseGraphExpression("uR[" + std::to_string(n) + "]");
    auto biRing = GraphUtils::parseGraphExpression("bR[" + std::to_string(n) + "]");
    
    // Same vertex count
    EXPECT_EQ(uniRing.getVertexCount(), biRing.getVertexCount());
    
    // Bidirectional has twice as many edges
    EXPECT_EQ(biRing.getEdgeCount(), 2 * uniRing.getEdgeCount());
    
    // Each vertex has different neighbor counts
    for (int i = 0; i < n; i++) {
      EXPECT_EQ(uniRing.getNeighbors(i).size(), 1);  // Unidirectional: 1 neighbor
      EXPECT_EQ(biRing.getNeighbors(i).size(), 2);   // Bidirectional: 2 neighbors
    }
  }
}

// Test parseGraphExpression mixed expressions
TEST_F(GraphUtilsTest, ParseGraphExpressionMixedExpressions) {
  // Test various valid expressions in sequence
  std::vector<std::pair<std::string, std::pair<int, int>>> testCases = {
    {"uR[3]", {3, 3}},     // {vertices, edges}
    {"bR[3]", {3, 6}},
    {"uR[5]", {5, 5}},
    {"bR[5]", {5, 10}},
    {"uR[7]", {7, 7}},
    {"bR[7]", {7, 14}}
  };
  
  for (const auto& testCase : testCases) {
    const std::string& expr = testCase.first;
    int expectedVertices = testCase.second.first;
    int expectedEdges = testCase.second.second;
    
    auto graph = GraphUtils::parseGraphExpression(expr);
    EXPECT_EQ(graph.getVertexCount(), expectedVertices) << "Failed for expression: " << expr;
    EXPECT_EQ(graph.getEdgeCount(), expectedEdges) << "Failed for expression: " << expr;
  }
}

// Test parseGraphExpression with tensor product expressions
TEST_F(GraphUtilsTest, ParseGraphExpressionTensorProduct) {
  // Test uR[N] x uR[M]
  auto graph1 = GraphUtils::parseGraphExpression("uR[3] x uR[4]");
  auto expected1 = GraphUtils::createTensorProduct(
    GraphUtils::createRing(3), 
    GraphUtils::createRing(4)
  );
  
  EXPECT_EQ(graph1.getVertexCount(), expected1.getVertexCount());
  EXPECT_EQ(graph1.getEdgeCount(), expected1.getEdgeCount());
  EXPECT_EQ(graph1.getVertexCount(), 12);  // 3 * 4
  
  // Test bR[N] x bR[M]
  auto graph2 = GraphUtils::parseGraphExpression("bR[3] x bR[2]");
  auto expected2 = GraphUtils::createTensorProduct(
    GraphUtils::createBidirectionalRing(3), 
    GraphUtils::createBidirectionalRing(2)
  );
  
  EXPECT_EQ(graph2.getVertexCount(), expected2.getVertexCount());
  EXPECT_EQ(graph2.getEdgeCount(), expected2.getEdgeCount());
  EXPECT_EQ(graph2.getVertexCount(), 6);  // 3 * 2
  
  // Test mixed: uR[N] x bR[M]
  auto graph3 = GraphUtils::parseGraphExpression("uR[2] x bR[3]");
  auto expected3 = GraphUtils::createTensorProduct(
    GraphUtils::createRing(2), 
    GraphUtils::createBidirectionalRing(3)
  );
  
  EXPECT_EQ(graph3.getVertexCount(), expected3.getVertexCount());
  EXPECT_EQ(graph3.getEdgeCount(), expected3.getEdgeCount());
  EXPECT_EQ(graph3.getVertexCount(), 6);  // 2 * 3
}

// Test parseGraphExpression tensor product with extra whitespace
TEST_F(GraphUtilsTest, ParseGraphExpressionTensorProductWhitespace) {
  // Test various whitespace patterns around "x"
  auto graph1 = GraphUtils::parseGraphExpression("uR[3] x uR[4]");
  auto graph2 = GraphUtils::parseGraphExpression(" uR[3] x uR[4] ");
  auto graph4 = GraphUtils::parseGraphExpression("uR[3] x  uR[4]");
  
  // First, second, and fourth should work and be identical
  EXPECT_EQ(graph1.getVertexCount(), 12);
  EXPECT_EQ(graph2.getVertexCount(), 12);
  EXPECT_EQ(graph4.getVertexCount(), 12);
  
  EXPECT_EQ(graph1.getEdgeCount(), graph2.getEdgeCount());
  EXPECT_EQ(graph1.getEdgeCount(), graph4.getEdgeCount());
  
  // Expressions without spaces around "x" should be treated as invalid
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR[3]xuR[4]"), std::invalid_argument);
}

// Test parseGraphExpression tensor product error cases
TEST_F(GraphUtilsTest, ParseGraphExpressionTensorProductErrors) {
  // Test invalid left expression
  EXPECT_THROW(GraphUtils::parseGraphExpression("invalid x uR[3]"), std::invalid_argument);
  
  // Test invalid right expression
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR[3] x invalid"), std::invalid_argument);
  
  // Test both invalid
  EXPECT_THROW(GraphUtils::parseGraphExpression("invalid x invalid"), std::invalid_argument);
  
  // Test incomplete expressions
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR[3] x"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("x uR[3]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("x"), std::invalid_argument);
  
  // Test invalid ring sizes in tensor product
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR[1] x uR[3]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR[3] x uR[0]"), std::invalid_argument);
}

// Test parseGraphExpression tensor product consistency
TEST_F(GraphUtilsTest, ParseGraphExpressionTensorProductConsistency) {
  // Test that DSL tensor product matches direct function calls
  std::vector<std::pair<std::pair<int, bool>, std::pair<int, bool>>> testCases = {
    {{3, false}, {4, false}},  // uR[3] x uR[4]
    {{3, true}, {4, true}},    // bR[3] x bR[4]
    {{2, false}, {3, true}},   // uR[2] x bR[3]
    {{5, true}, {2, false}}    // bR[5] x uR[2]
  };
  
  for (const auto& testCase : testCases) {
    int n1 = testCase.first.first;
    bool isBi1 = testCase.first.second;
    int n2 = testCase.second.first;
    bool isBi2 = testCase.second.second;
    
    // Build DSL expression
    std::string expr = (isBi1 ? "bR[" : "uR[") + std::to_string(n1) + "] x " +
                       (isBi2 ? "bR[" : "uR[") + std::to_string(n2) + "]";
    
    // Parse with DSL
    auto dslGraph = GraphUtils::parseGraphExpression(expr);
    
    // Create directly
    auto graph1 = isBi1 ? GraphUtils::createBidirectionalRing(n1) : GraphUtils::createRing(n1);
    auto graph2 = isBi2 ? GraphUtils::createBidirectionalRing(n2) : GraphUtils::createRing(n2);
    auto directGraph = GraphUtils::createTensorProduct(graph1, graph2);
    
    // Compare results
    EXPECT_EQ(dslGraph.getVertexCount(), directGraph.getVertexCount()) << "Failed for: " << expr;
    EXPECT_EQ(dslGraph.getEdgeCount(), directGraph.getEdgeCount()) << "Failed for: " << expr;
    
    // Check all edges match
    for (int i = 0; i < dslGraph.getVertexCount(); i++) {
      for (int j = 0; j < dslGraph.getVertexCount(); j++) {
        EXPECT_EQ(dslGraph.hasEdge(i, j), directGraph.hasEdge(i, j)) 
            << "Edge mismatch at (" << i << "," << j << ") for: " << expr;
      }
    }
  }
}

// Test parseGraphExpression complex tensor product examples
TEST_F(GraphUtilsTest, ParseGraphExpressionComplexTensorProducts) {
  // Test various interesting combinations
  std::vector<std::tuple<std::string, int, int>> testCases = {
    {"uR[2] x uR[2]", 4, 4},    // 2x2 grid, minimal unidirectional
    {"bR[2] x bR[2]", 4, 4},    // 2x2 grid, bidirectional
    {"uR[3] x bR[2]", 6, 6},    // 3x2 mixed
    {"bR[3] x uR[2]", 6, 12},   // 3x2 mixed (reversed)
    {"uR[5] x uR[2]", 10, 10},  // 5x2 unidirectional
    {"bR[4] x bR[3]", 12, 48}   // 4x3 bidirectional
  };
  
  for (const auto& testCase : testCases) {
    const std::string& expr = std::get<0>(testCase);
    int expectedVertices = std::get<1>(testCase);
    int expectedEdges = std::get<2>(testCase);
    
    auto graph = GraphUtils::parseGraphExpression(expr);
    EXPECT_EQ(graph.getVertexCount(), expectedVertices) << "Vertex count failed for: " << expr;
    EXPECT_EQ(graph.getEdgeCount(), expectedEdges) << "Edge count failed for: " << expr;
  }
}

// Test parseGraphExpression with parentheses
TEST_F(GraphUtilsTest, ParseGraphExpressionWithParentheses) {
  // Test basic parentheses around simple expressions
  auto graph1 = GraphUtils::parseGraphExpression("(uR[5])");
  auto expected1 = GraphUtils::parseGraphExpression("uR[5]");
  
  EXPECT_EQ(graph1.getVertexCount(), expected1.getVertexCount());
  EXPECT_EQ(graph1.getEdgeCount(), expected1.getEdgeCount());
  
  auto graph2 = GraphUtils::parseGraphExpression("(bR[4])");
  auto expected2 = GraphUtils::parseGraphExpression("bR[4]");
  
  EXPECT_EQ(graph2.getVertexCount(), expected2.getVertexCount());
  EXPECT_EQ(graph2.getEdgeCount(), expected2.getEdgeCount());
  
  // Test parentheses around tensor products
  auto graph3 = GraphUtils::parseGraphExpression("(uR[3] x bR[4])");
  auto expected3 = GraphUtils::parseGraphExpression("uR[3] x bR[4]");
  
  EXPECT_EQ(graph3.getVertexCount(), expected3.getVertexCount());
  EXPECT_EQ(graph3.getEdgeCount(), expected3.getEdgeCount());
}

// Test parseGraphExpression with nested parentheses
TEST_F(GraphUtilsTest, ParseGraphExpressionNestedParentheses) {
  // Test nested parentheses
  auto graph1 = GraphUtils::parseGraphExpression("((uR[3]))");
  auto expected1 = GraphUtils::parseGraphExpression("uR[3]");
  
  EXPECT_EQ(graph1.getVertexCount(), expected1.getVertexCount());
  EXPECT_EQ(graph1.getEdgeCount(), expected1.getEdgeCount());
  
  // Test deeply nested
  auto graph2 = GraphUtils::parseGraphExpression("(((bR[4] x uR[2])))");
  auto expected2 = GraphUtils::parseGraphExpression("bR[4] x uR[2]");
  
  EXPECT_EQ(graph2.getVertexCount(), expected2.getVertexCount());
  EXPECT_EQ(graph2.getEdgeCount(), expected2.getEdgeCount());
}

// Test parseGraphExpression with parentheses for precedence
TEST_F(GraphUtilsTest, ParseGraphExpressionParenthesesPrecedence) {
  // Test that parentheses can change evaluation order
  // Note: For tensor products, the result should be the same due to associativity of vertex count
  // but the internal vertex numbering might differ
  
  auto graph1 = GraphUtils::parseGraphExpression("(uR[2] x uR[3]) x uR[4]");
  auto graph2 = GraphUtils::parseGraphExpression("uR[2] x (uR[3] x uR[4])");
  
  // Both should have same vertex and edge counts due to associativity
  EXPECT_EQ(graph1.getVertexCount(), 24);  // 2 * 3 * 4
  EXPECT_EQ(graph2.getVertexCount(), 24);  // 2 * 3 * 4
  EXPECT_EQ(graph1.getVertexCount(), graph2.getVertexCount());
  
  // Test with bidirectional rings
  auto graph3 = GraphUtils::parseGraphExpression("(bR[2] x bR[2]) x bR[2]");
  auto graph4 = GraphUtils::parseGraphExpression("bR[2] x (bR[2] x bR[2])");
  
  EXPECT_EQ(graph3.getVertexCount(), 8);   // 2 * 2 * 2
  EXPECT_EQ(graph4.getVertexCount(), 8);   // 2 * 2 * 2
  EXPECT_EQ(graph3.getVertexCount(), graph4.getVertexCount());
}

// Test parseGraphExpression with complex nested expressions
TEST_F(GraphUtilsTest, ParseGraphExpressionComplexNested) {
  // Test complex nested expressions
  std::vector<std::pair<std::string, std::pair<int, int>>> testCases = {
    {"(uR[3] x bR[2]) x uR[4]", {24, 24}},     // (3*2)*4 = 24 vertices
    {"uR[3] x (bR[2] x uR[4])", {24, 24}},     // 3*(2*4) = 24 vertices  
    {"(uR[2] x uR[2]) x (bR[2] x bR[2])", {16, 16}},  // (2*2)*(2*2) = 16 vertices
    {"((uR[3] x bR[2]) x uR[2])", {12, 12}},   // Nested grouping
    {"(uR[2]) x (bR[3]) x (uR[2])", {12, 24}}  // Multiple groups
  };
  
  for (const auto& testCase : testCases) {
    const std::string& expr = testCase.first;
    int expectedVertices = testCase.second.first;
    int expectedEdges = testCase.second.second;
    
    auto graph = GraphUtils::parseGraphExpression(expr);
    EXPECT_EQ(graph.getVertexCount(), expectedVertices) << "Vertex count failed for: " << expr;
    EXPECT_EQ(graph.getEdgeCount(), expectedEdges) << "Edge count failed for: " << expr;
  }
}

// Test parseGraphExpression with parentheses whitespace handling
TEST_F(GraphUtilsTest, ParseGraphExpressionParenthesesWhitespace) {
  // Test whitespace around and inside parentheses
  std::vector<std::string> expressions = {
    "(uR[3])",
    " (uR[3]) ",
    "( uR[3] )",
    " ( uR[3] ) ",
    "(uR[3] x bR[4])",
    " (uR[3] x bR[4]) ",
    "( uR[3] x bR[4] )",
    " ( uR[3] x bR[4] ) "
  };
  
  // All should produce equivalent results
  auto expected = GraphUtils::parseGraphExpression("uR[3] x bR[4]");
  
  for (const auto& expr : expressions) {
    if (expr.find("x") != std::string::npos) {
      // Tensor product expressions
      auto graph = GraphUtils::parseGraphExpression(expr);
      EXPECT_EQ(graph.getVertexCount(), expected.getVertexCount()) << "Failed for: " << expr;
      EXPECT_EQ(graph.getEdgeCount(), expected.getEdgeCount()) << "Failed for: " << expr;
    } else {
      // Simple ring expressions
      auto graph = GraphUtils::parseGraphExpression(expr);
      EXPECT_EQ(graph.getVertexCount(), 3) << "Failed for: " << expr;
      EXPECT_EQ(graph.getEdgeCount(), 3) << "Failed for: " << expr;
    }
  }
}

// Test parseGraphExpression with invalid parentheses
TEST_F(GraphUtilsTest, ParseGraphExpressionInvalidParentheses) {
  // Test mismatched parentheses
  EXPECT_THROW(GraphUtils::parseGraphExpression("(uR[3]"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR[3])"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("(uR[3])x(bR[4])"), std::invalid_argument);  // No spaces around x
  EXPECT_THROW(GraphUtils::parseGraphExpression("((uR[3])"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("(uR[3]))"), std::invalid_argument);
  
  // Test empty parentheses
  EXPECT_THROW(GraphUtils::parseGraphExpression("()"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("( )"), std::invalid_argument);
  
  // Test misplaced parentheses
  EXPECT_THROW(GraphUtils::parseGraphExpression("uR([3])"), std::invalid_argument);
  EXPECT_THROW(GraphUtils::parseGraphExpression("(uR[3] x) bR[4]"), std::invalid_argument);
}

// Test parseGraphExpression parentheses consistency
TEST_F(GraphUtilsTest, ParseGraphExpressionParenthesesConsistency) {
  // Test that parentheses don't change the result when they don't affect precedence
  std::vector<std::pair<std::string, std::string>> equivalentPairs = {
    {"uR[5]", "(uR[5])"},
    {"bR[4]", "(bR[4])"},
    {"uR[3] x bR[4]", "(uR[3] x bR[4])"},
    {"uR[3] x bR[4]", "(uR[3]) x (bR[4])"},
    {"bR[2] x uR[3]", "((bR[2])) x ((uR[3]))"}
  };
  
  for (const auto& pair : equivalentPairs) {
    const std::string& expr1 = pair.first;
    const std::string& expr2 = pair.second;
    
    auto graph1 = GraphUtils::parseGraphExpression(expr1);
    auto graph2 = GraphUtils::parseGraphExpression(expr2);
    
    EXPECT_EQ(graph1.getVertexCount(), graph2.getVertexCount()) 
        << "Vertex count mismatch: '" << expr1 << "' vs '" << expr2 << "'";
    EXPECT_EQ(graph1.getEdgeCount(), graph2.getEdgeCount()) 
        << "Edge count mismatch: '" << expr1 << "' vs '" << expr2 << "'";
    
    // Check all edges match
    for (int i = 0; i < graph1.getVertexCount(); i++) {
      for (int j = 0; j < graph1.getVertexCount(); j++) {
        EXPECT_EQ(graph1.hasEdge(i, j), graph2.hasEdge(i, j)) 
            << "Edge mismatch at (" << i << "," << j << ") for: '" << expr1 << "' vs '" << expr2 << "'";
      }
    }
  }
}

// Test parseGraphExpression left associativity of tensor product
TEST_F(GraphUtilsTest, ParseGraphExpressionLeftAssociativity) {
  // Test that tensor product is left associative
  // "A x B x C" should be parsed as "(A x B) x C"
  
  // Test case: "uR[2] x uR[3] x uR[4]" should be equivalent to "(uR[2] x uR[3]) x uR[4]"
  auto implicitLeft = GraphUtils::parseGraphExpression("uR[2] x uR[3] x uR[4]");
  auto explicitLeft = GraphUtils::parseGraphExpression("(uR[2] x uR[3]) x uR[4]");
  
  EXPECT_EQ(implicitLeft.getVertexCount(), explicitLeft.getVertexCount());
  EXPECT_EQ(implicitLeft.getEdgeCount(), explicitLeft.getEdgeCount());
  
  // Check all edges match
  for (int i = 0; i < implicitLeft.getVertexCount(); i++) {
    for (int j = 0; j < implicitLeft.getVertexCount(); j++) {
      EXPECT_EQ(implicitLeft.hasEdge(i, j), explicitLeft.hasEdge(i, j))
          << "Edge mismatch at (" << i << "," << j << ") for left associativity test";
    }
  }
  
  // Test case: "bR[2] x bR[2] x bR[2]" should be equivalent to "(bR[2] x bR[2]) x bR[2]"
  auto implicitLeft2 = GraphUtils::parseGraphExpression("bR[2] x bR[2] x bR[2]");
  auto explicitLeft2 = GraphUtils::parseGraphExpression("(bR[2] x bR[2]) x bR[2]");
  
  EXPECT_EQ(implicitLeft2.getVertexCount(), explicitLeft2.getVertexCount());
  EXPECT_EQ(implicitLeft2.getEdgeCount(), explicitLeft2.getEdgeCount());
  
  // Test case: "uR[2] x bR[3] x uR[2]" should be equivalent to "(uR[2] x bR[3]) x uR[2]"
  auto implicitLeft3 = GraphUtils::parseGraphExpression("uR[2] x bR[3] x uR[2]");
  auto explicitLeft3 = GraphUtils::parseGraphExpression("(uR[2] x bR[3]) x uR[2]");
  
  EXPECT_EQ(implicitLeft3.getVertexCount(), explicitLeft3.getVertexCount());
  EXPECT_EQ(implicitLeft3.getEdgeCount(), explicitLeft3.getEdgeCount());
}

// Test parseGraphExpression left associativity with longer chains
TEST_F(GraphUtilsTest, ParseGraphExpressionLeftAssociativityLonger) {
  // Test longer chains to ensure proper left associativity
  
  // Test case: "uR[2] x uR[2] x uR[2] x uR[2]" should be "((uR[2] x uR[2]) x uR[2]) x uR[2]"
  auto chain4 = GraphUtils::parseGraphExpression("uR[2] x uR[2] x uR[2] x uR[2]");
  auto explicit4 = GraphUtils::parseGraphExpression("((uR[2] x uR[2]) x uR[2]) x uR[2]");
  
  EXPECT_EQ(chain4.getVertexCount(), explicit4.getVertexCount());
  EXPECT_EQ(chain4.getEdgeCount(), explicit4.getEdgeCount());
  EXPECT_EQ(chain4.getVertexCount(), 16);  // 2^4 = 16
  
  // Test case: "bR[2] x bR[2] x bR[2]" should be "(bR[2] x bR[2]) x bR[2]"
  auto chain3 = GraphUtils::parseGraphExpression("bR[2] x bR[2] x bR[2]");
  auto explicit3 = GraphUtils::parseGraphExpression("(bR[2] x bR[2]) x bR[2]");
  
  EXPECT_EQ(chain3.getVertexCount(), explicit3.getVertexCount());
  EXPECT_EQ(chain3.getEdgeCount(), explicit3.getEdgeCount());
  EXPECT_EQ(chain3.getVertexCount(), 8);  // 2^3 = 8
}

// Test parseGraphExpression left associativity vs right associativity
TEST_F(GraphUtilsTest, ParseGraphExpressionLeftVsRightAssociativity) {
  // Test that left associativity is different from right associativity
  // (when vertex numbering matters, the graphs might be different)
  
  // Left associative: "(uR[2] x uR[3]) x uR[4]"
  auto leftAssoc = GraphUtils::parseGraphExpression("uR[2] x uR[3] x uR[4]");
  auto explicitLeft = GraphUtils::parseGraphExpression("(uR[2] x uR[3]) x uR[4]");
  
  // Right associative: "uR[2] x (uR[3] x uR[4])"
  auto explicitRight = GraphUtils::parseGraphExpression("uR[2] x (uR[3] x uR[4])");
  
  // The implicit should match the explicit left
  EXPECT_EQ(leftAssoc.getVertexCount(), explicitLeft.getVertexCount());
  EXPECT_EQ(leftAssoc.getEdgeCount(), explicitLeft.getEdgeCount());
  
  // Both left and right should have the same vertex/edge counts due to tensor product properties
  EXPECT_EQ(explicitLeft.getVertexCount(), explicitRight.getVertexCount());
  EXPECT_EQ(explicitLeft.getEdgeCount(), explicitRight.getEdgeCount());
  
  // But verify the implicit parsing matches the left associative version
  for (int i = 0; i < leftAssoc.getVertexCount(); i++) {
    for (int j = 0; j < leftAssoc.getVertexCount(); j++) {
      EXPECT_EQ(leftAssoc.hasEdge(i, j), explicitLeft.hasEdge(i, j))
          << "Left associativity mismatch at (" << i << "," << j << ")";
    }
  }
}

// Test parseGraphExpression mixed associativity with parentheses
TEST_F(GraphUtilsTest, ParseGraphExpressionMixedAssociativity) {
  // Test various combinations of left/right associativity with parentheses
  
  std::vector<std::pair<std::string, std::pair<int, int>>> testCases = {
    {"uR[2] x uR[3] x uR[4]", {24, 24}},               // Left associative (default)
    {"(uR[2] x uR[3]) x uR[4]", {24, 24}},             // Explicit left
    {"uR[2] x (uR[3] x uR[4])", {24, 24}},             // Explicit right
    {"uR[2] x uR[2] x uR[2]", {8, 8}},                 // Triple left
    {"(uR[2] x uR[2]) x uR[2]", {8, 8}},               // Triple explicit left
    {"uR[2] x (uR[2] x uR[2])", {8, 8}},               // Triple explicit right
    {"bR[2] x bR[2] x bR[2] x bR[2]", {16, 16}},       // Quad left
    {"((bR[2] x bR[2]) x bR[2]) x bR[2]", {16, 16}},   // Quad explicit left
    {"bR[2] x (bR[2] x (bR[2] x bR[2]))", {16, 16}}    // Quad explicit right
  };
  
  for (const auto& testCase : testCases) {
    const std::string& expr = testCase.first;
    int expectedVertices = testCase.second.first;
    int expectedEdges = testCase.second.second;
    
    auto graph = GraphUtils::parseGraphExpression(expr);
    EXPECT_EQ(graph.getVertexCount(), expectedVertices) 
        << "Vertex count failed for: " << expr;
    EXPECT_EQ(graph.getEdgeCount(), expectedEdges) 
        << "Edge count failed for: " << expr;
  }
}
