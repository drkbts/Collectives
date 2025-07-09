load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

package(default_visibility = ["//visibility:public"])

# Graph library
cc_library(
    name = "graph",
    srcs = [
        "src/graph.cc",
        "src/graph_utils.cc",
        "src/graph_connectivity.cc",
        "src/spanning_tree.cc",
        "src/hamiltonian.cc",
        "src/tsp_boost.cc",
        "src/graph_partitioning.cc",
        "src/graph_collective.cc",
    ],
    hdrs = [
        "src/graph.h",
        "src/graph_utils.h",
        "src/graph_connectivity.h",
        "src/spanning_tree.h",
        "src/hamiltonian.h",
        "src/tsp_boost.h",
        "src/graph_partitioning.h",
        "src/graph_collective.h",
    ],
    deps = [
        "@boost.graph",
    ],
    strip_include_prefix = "src",
)

# Test executable using Google Test
cc_test(
    name = "graph_test",
    srcs = ["tests/graph_test.cc"],
    deps = [
        ":graph",
        "@googletest//:gtest_main",
    ],
)

# Test executable for graph utilities
cc_test(
    name = "graph_utils_test",
    srcs = ["tests/graph_utils_test.cc"],
    deps = [
        ":graph",
        "@googletest//:gtest_main",
    ],
)

# Test executable for spanning tree
cc_test(
    name = "spanning_tree_test",
    srcs = ["tests/spanning_tree_test.cc"],
    deps = [
        ":graph",
        "@googletest//:gtest_main",
    ],
)

# Test executable for graph partitioning
cc_test(
    name = "graph_partitioning_test",
    srcs = ["tests/graph_partitioning_test.cc"],
    deps = [
        ":graph",
        "@googletest//:gtest_main",
    ],
)

# Test executable for graph collective operations
cc_test(
    name = "graph_collective_test",
    srcs = ["tests/graph_collective_test.cc"],
    deps = [
        ":graph",
        "@googletest//:gtest_main",
    ],
)
