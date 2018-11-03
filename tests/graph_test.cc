#include <memory>

#include "catch.hpp"

#include "graph.h"

// Check Graph class
TEST_CASE("Graph class and shortest-path is checked", "[graph][sp][od]") {
  // Tolerance
  const double Tolerance = 1.E-7;
  // Path
  const std::string path = "../";

  // Test directed graph
  SECTION("Test SSSP in directed graph") {
    // Set graph properties
    const bool directed = true;

    // Create graph object
    auto graph = std::make_unique<abm::Graph>(directed);
    // Create a simple example graph
    graph->generate_simple_graph();

    // Run Dijkstra Priority Queue
    abm::graph::vertex_t source = 1;
    abm::graph::vertex_t destination = 3;
    auto sp = graph->dijkstra_priority_queue(source, destination);
    auto distances = sp.distances;

    // Check distances
    REQUIRE(distances.size() == graph->nvertices());
    // Check shortest path
    REQUIRE(distances.at(3) == Approx(7.2).epsilon(Tolerance));

    // Check update edge
    SECTION("Check update edge") {
      // Update edge (1, 3) to weight 3.7
      graph->update_edge(1, 3, 3.7);
      // Run Dijkstra Priority Queue
      sp = graph->dijkstra_priority_queue(source, destination);
      // Get distances
      distances = sp.distances;
      // Check shortest path
      REQUIRE(distances.at(3) == Approx(3.7).epsilon(Tolerance));
    }

    // Check remove edge
    SECTION("Check remove edge") {
      // Remove edge (3, 1)
      graph->remove_edge(3, 1);
      // Run Dijkstra Priority Queue
      sp = graph->dijkstra_priority_queue(source, destination);
      // Get distances
      distances = sp.distances;
      // Check shortest path
      REQUIRE(distances.at(3) == Approx(7.2).epsilon(Tolerance));

      // Remove edge (2, 4)
      graph->remove_edge(2, 4);
      // Run Dijkstra Priority Queue
      sp = graph->dijkstra_priority_queue(source, destination);
      // Get distances
      distances = sp.distances;
      // Check shortest path
      REQUIRE(distances.at(3) == Approx(9.1).epsilon(Tolerance));
    }
  }

  // Test Single Source Multiple Destinations directed graph
  SECTION("Test Single Source Multiple Destinations in directed graph") {
    // Set graph properties
    const bool directed = true;

    // Create graph object
    auto graph = std::make_unique<abm::Graph>(directed);
    // Create a simple example graph
    graph->generate_simple_graph();

    // Run Dijkstra Priority Queue
    abm::graph::vertex_t source = 1;
    std::vector<abm::graph::vertex_t> destinations{2, 3};
    auto sp = graph->dijkstra_priority_queue(source, destinations);
    auto distances = sp.distances;

    // Check distances
    REQUIRE(distances.size() == graph->nvertices());
    // Check shortest path
    REQUIRE(distances.at(2) == Approx(1.5).epsilon(Tolerance));
    // Check shortest path
    REQUIRE(distances.at(3) == Approx(7.2).epsilon(Tolerance));
  }

  // Test undirected graph
  SECTION("Test SSSP in undirected graph") {
    const bool directed = false;
    // Create graph object
    auto graph = std::make_unique<abm::Graph>(directed);
    // Create a simple example graph
    graph->generate_simple_graph();

    // Run Dijkstra Priority Queue
    abm::graph::vertex_t source = 1;
    abm::graph::vertex_t destination = 3;
    auto sp = graph->dijkstra_priority_queue(source, destination);
    // Get distances
    auto distances = sp.distances;

    // Check distances
    REQUIRE(distances.size() == graph->nvertices());
    // Check shortest path
    REQUIRE(distances.at(3) == Approx(5.6).epsilon(Tolerance));

    // Check remove edge
    SECTION("Check remove edge") {
      // Remove edge (3, 1)
      graph->remove_edge(3, 1);
      // Run Dijkstra Priority Queue
      sp = graph->dijkstra_priority_queue(source, destination);
      // Get distances
      distances = sp.distances;
      // Check shortest path
      REQUIRE(distances.at(3) == Approx(9.1).epsilon(Tolerance));
    }
  }

  // Test Single Source Multiple Destinations undirected graph
  SECTION("Test Single Source Multiple Destinations in undirected graph") {
    // Set graph properties
    const bool directed = false;

    // Create graph object
    auto graph = std::make_unique<abm::Graph>(directed);
    // Create a simple example graph
    graph->generate_simple_graph();

    // Run Dijkstra Priority Queue
    abm::graph::vertex_t source = 1;
    std::vector<abm::graph::vertex_t> destinations{2, 3};
    auto sp = graph->dijkstra_priority_queue(source, destinations);
    auto distances = sp.distances;

    // Check distances
    REQUIRE(distances.size() == graph->nvertices());
    // Check shortest path
    REQUIRE(distances.at(2) == Approx(1.5).epsilon(Tolerance));
    // Check shortest path
    REQUIRE(distances.at(3) == Approx(5.6).epsilon(Tolerance));
  }

  SECTION("Test SSSP in directed graph from file") {
    // Set graph properties
    const bool directed = true;

    // Create graph object
    auto graph = std::make_unique<abm::Graph>(directed);
    // Read MatrixMarket file
    const std::string filename = path + "network.mtx";
    // Read file
    REQUIRE(graph->read_graph_matrix_market(filename) == true);

    // Run Dijkstra Priority Queue
    abm::graph::vertex_t source = 1020;
    abm::graph::vertex_t destination = 20;
    auto sp = graph->dijkstra_priority_queue(source, destination);
    // Get distances
    auto distances = sp.distances;

    // Check distances
    REQUIRE(distances.size() == graph->nvertices());
    // Check shortest path
    REQUIRE(distances.at(20) == Approx(12409.660000000002).epsilon(Tolerance));

    SECTION("Test non-existant file") {
      // Create graph object
      auto graph = std::make_unique<abm::Graph>(directed);
      // Read MatrixMarket file
      const std::string filename = "nofile.mtx";
      // Read file should fail
      REQUIRE(graph->read_graph_matrix_market(filename) == false);
    }
  }
}
