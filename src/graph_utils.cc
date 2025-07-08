#include "graph_utils.h"
#include <stdexcept>
#include <regex>
#include <string>
#include <algorithm>
#include <cctype>

namespace GraphUtils {

UnidirectionalGraph createRing(int n) {
    if (n <= 1) {
        throw std::invalid_argument("Ring size must be greater than 1");
    }
    
    UnidirectionalGraph graph;
    
    // Create vertices 0, 1, 2, ..., n-1
    for (int i = 0; i < n; i++) {
        graph.addVertex(i);
    }
    
    // Create edges to form a ring: 0->1, 1->2, ..., (n-1)->0
    for (int i = 0; i < n - 1; i++) {
        graph.addEdge(i, i + 1);
    }
    
    // Close the ring: (n-1) -> 0
    graph.addEdge(n - 1, 0);
    
    return graph;
}

UnidirectionalGraph createBidirectionalRing(int n) {
    if (n <= 1) {
        throw std::invalid_argument("Ring size must be greater than 1");
    }
    
    UnidirectionalGraph graph;
    
    // Create vertices 0, 1, 2, ..., n-1
    for (int i = 0; i < n; i++) {
        graph.addVertex(i);
    }
    
    // Create forward edges: 0->1, 1->2, ..., (n-1)->0
    for (int i = 0; i < n - 1; i++) {
        graph.addEdge(i, i + 1);
    }
    graph.addEdge(n - 1, 0);
    
    // Create backward edges: 1->0, 2->1, ..., 0->(n-1)
    for (int i = 1; i < n; i++) {
        graph.addEdge(i, i - 1);
    }
    graph.addEdge(0, n - 1);
    
    return graph;
}

UnidirectionalGraph createTensorProduct(const UnidirectionalGraph& g1, const UnidirectionalGraph& g2) {
    if (g1.getVertexCount() == 0 || g2.getVertexCount() == 0) {
        throw std::invalid_argument("Both graphs must have at least one vertex");
    }
    
    UnidirectionalGraph result;
    
    // Get all vertices from both graphs
    std::vector<int> vertices1;
    std::vector<int> vertices2;
    
    // Find all vertices in g1 by checking for vertex existence
    // Since we don't have direct access to vertex list, we'll use a different approach
    // We'll iterate through possible vertex values and check if they exist
    int maxVertex2 = -1;
    
    // Find the maximum vertex in each graph by testing
    for (int v = 0; v < 1000; v++) {  // Reasonable upper bound
        if (g1.hasVertex(v)) {
            vertices1.push_back(v);
        }
        if (g2.hasVertex(v)) {
            vertices2.push_back(v);
            maxVertex2 = v;
        }
    }
    
    // Helper function to map vertex pair (u,v) to single integer
    auto pairToVertex = [&](int u, int v) -> int {
        return u * (maxVertex2 + 1) + v;
    };
    
    // Create all vertices in the tensor product
    for (int u : vertices1) {
        for (int v : vertices2) {
            result.addVertex(pairToVertex(u, v));
        }
    }
    
    // Create edges: (u1,v1) -> (u2,v2) exists iff u1->u2 in g1 AND v1->v2 in g2
    for (int u1 : vertices1) {
        for (int v1 : vertices2) {
            for (int u2 : vertices1) {
                for (int v2 : vertices2) {
                    if (g1.hasEdge(u1, u2) && g2.hasEdge(v1, v2)) {
                        result.addEdge(pairToVertex(u1, v1), pairToVertex(u2, v2));
                    }
                }
            }
        }
    }
    
    return result;
}

UnidirectionalGraph createTorus(const std::vector<int>& dimensions) {
    if (dimensions.empty()) {
        throw std::invalid_argument("Dimensions vector cannot be empty");
    }
    
    // Check that all dimensions are > 1
    for (int dim : dimensions) {
        if (dim <= 1) {
            throw std::invalid_argument("All dimensions must be greater than 1");
        }
    }
    
    // Start with the first bidirectional ring
    UnidirectionalGraph result = createBidirectionalRing(dimensions[0]);
    
    // Take tensor product with each subsequent bidirectional ring
    for (size_t i = 1; i < dimensions.size(); i++) {
        UnidirectionalGraph nextRing = createBidirectionalRing(dimensions[i]);
        result = createTensorProduct(result, nextRing);
    }
    
    return result;
}

UnidirectionalGraph parseGraphExpression(const std::string& expression) {
    // Trim leading and trailing whitespace
    std::string trimmed = expression;
    size_t start = trimmed.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) {
        throw std::invalid_argument("Empty or whitespace-only expression");
    }
    size_t end = trimmed.find_last_not_of(" \t\n\r");
    trimmed = trimmed.substr(start, end - start + 1);
    
    // Check for outermost parentheses and strip them if they wrap the entire expression
    if (trimmed.length() >= 2 && trimmed.front() == '(' && trimmed.back() == ')') {
        // Verify that these are actually the outermost parentheses
        int parenCount = 0;
        bool isOutermost = true;
        for (size_t i = 0; i < trimmed.length() - 1; i++) {
            if (trimmed[i] == '(') parenCount++;
            else if (trimmed[i] == ')') parenCount--;
            if (parenCount == 0) {
                isOutermost = false;
                break;
            }
        }
        if (isOutermost) {
            trimmed = trimmed.substr(1, trimmed.length() - 2);
            return parseGraphExpression(trimmed);
        }
    }
    
    // Find rightmost tensor product operator " x " at the top level for left associativity
    // For "A x B x C", we want to parse as "(A x B) x C", so we split at the rightmost "x"
    int parenDepth = 0;
    int rightmostXPos = -1;
    
    if (trimmed.length() >= 3) {
        for (size_t i = 0; i <= trimmed.length() - 3; i++) {
            if (trimmed[i] == '(') {
                parenDepth++;
            } else if (trimmed[i] == ')') {
                parenDepth--;
            } else if (parenDepth == 0 && i + 2 < trimmed.length() && trimmed.substr(i, 3) == " x ") {
                // Found top-level tensor product operator - remember the rightmost one
                rightmostXPos = static_cast<int>(i);
            }
        }
        
        // If we found a tensor product operator, split at the rightmost one
        if (rightmostXPos != -1) {
            std::string leftExpr = trimmed.substr(0, rightmostXPos);
            std::string rightExpr = trimmed.substr(rightmostXPos + 3);
            
            // Recursively parse both sides
            try {
                UnidirectionalGraph leftGraph = parseGraphExpression(leftExpr);
                UnidirectionalGraph rightGraph = parseGraphExpression(rightExpr);
                return createTensorProduct(leftGraph, rightGraph);
            } catch (const std::exception& e) {
                throw std::invalid_argument("Error in tensor product expression '" + expression + "': " + e.what());
            }
        }
    }
    
    // Remove whitespace for simple expressions (reuse trimmed variable)
    trimmed.erase(std::remove_if(trimmed.begin(), trimmed.end(), ::isspace), trimmed.end());
    
    // Regular expression for unidirectional ring: uR[N]
    std::regex uRingRegex(R"(^uR\[(\d+)\]$)");
    std::smatch matches;
    
    if (std::regex_match(trimmed, matches, uRingRegex)) {
        // Extract the number N
        std::string numStr = matches[1].str();
        
        try {
            int n = std::stoi(numStr);
            if (n <= 1) {
                throw std::invalid_argument("Ring size must be greater than 1");
            }
            return createRing(n);
        } catch (const std::invalid_argument& e) {
            throw std::invalid_argument("Invalid ring size: " + numStr);
        } catch (const std::out_of_range& e) {
            throw std::out_of_range("Ring size out of range: " + numStr);
        }
    }
    
    // Regular expression for bidirectional ring: bR[N]
    std::regex bRingRegex(R"(^bR\[(\d+)\]$)");
    
    if (std::regex_match(trimmed, matches, bRingRegex)) {
        // Extract the number N
        std::string numStr = matches[1].str();
        
        try {
            int n = std::stoi(numStr);
            if (n <= 1) {
                throw std::invalid_argument("Ring size must be greater than 1");
            }
            return createBidirectionalRing(n);
        } catch (const std::invalid_argument& e) {
            throw std::invalid_argument("Invalid ring size: " + numStr);
        } catch (const std::out_of_range& e) {
            throw std::out_of_range("Ring size out of range: " + numStr);
        }
    }
    
    // If no pattern matches, throw an error
    throw std::invalid_argument("Invalid graph expression: " + expression);
}

}  // namespace GraphUtils
