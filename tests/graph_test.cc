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

    auto path = graph->dijkstra(source, destination);
    // Check distances
    REQUIRE(path.size() == 150);

    auto edges = graph->dijkstra_edges(source, destination);
    // Check distances
    REQUIRE(edges.size() == 150);

    std::vector<abm::graph::vertex_t> route{
        1721,  75304, 1345,  69459, 1341,  1338,  1335,  1334,  72156, 6806,
        6763,  72153, 48665, 62070, 60701, 62063, 60689, 62051, 60679, 60677,
        60660, 60669, 60656, 60668, 62067, 79653, 79577, 67951, 49942, 75674,
        52097, 7103,  7090,  6980,  30489, 6857,  7143,  6610,  25420, 40280,
        48602, 48552, 66759, 78698, 49953, 40426, 9928,  40286, 40288, 70488,
        40573, 40579, 40302, 40307, 40312, 40310, 40315, 40767, 46485, 66033,
        65978, 66009, 65922, 58060, 25392, 43947, 25381, 76975, 25376, 65932,
        18039, 8806,  56885, 65963, 76977, 66012, 18064, 76971, 65947, 18531,
        78995, 4382,  67570, 40853, 67558, 67550, 16414, 16372, 43386, 30410,
        43377, 967,   63877, 60695, 77971, 73684, 17261, 54418, 76054, 22698,
        40752, 62481, 58062, 27544, 216,   79384, 220,   72534, 922,   53960,
        53745, 7054,  7020,  7019,  54036, 2891,  7029,  7409,  6987,  73837,
        69117, 56071, 6943,  27731, 67163, 76786, 67164, 71638, 66197, 67171,
        61076, 68098, 71778, 71779, 50152, 50019, 67652, 55748, 67656, 55749,
        3173,  21875, 35,    23,    24067, 24053, 22030, 22026, 22122, 65887};

    for (unsigned i = 0; i < route.size(); ++i)
      REQUIRE(edges.at(i) == route.at(i));
  }

  SECTION("Test SSSP in OSM directed graph from file") {
    // Set graph properties
    const bool directed = true;

    // Create graph object
    auto graph = std::make_unique<abm::Graph>(directed);
    // Read MatrixMarket file
    const std::string filename = path + "osm/edges.csv";
    // Read file
    REQUIRE(graph->read_osm_graph(filename) == true);

    // Run Dijkstra Priority Queue
    abm::graph::vertex_t source = 65358739;
    abm::graph::vertex_t destination = 2304626340;
    auto sp = graph->dijkstra_edges(source, destination);

    // Check distances
    REQUIRE(sp.size() == 58);

    SECTION("Test non-existant file") {
      // Create graph object
      auto graph = std::make_unique<abm::Graph>(directed);
      // Read MatrixMarket file
      const std::string filename = "osm/nofile.csv";
      // Read file should fail
      REQUIRE(graph->read_osm_graph(filename) == false);
    }
  }
}
