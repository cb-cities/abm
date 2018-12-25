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
    const auto path = graph->dijkstra_vertices(source, destination);

    // Check distances
    REQUIRE(path.size() == 3);
    // Check shortest path
    REQUIRE(graph->path_cost(path) == Approx(7.2).epsilon(Tolerance));

    // Check update edge
    SECTION("Check update edge") {
      // Update edge (1, 3) to weight 3.7
      graph->update_edge(1, 3, 3.7);
      // Run Dijkstra Priority Queue
      const auto path = graph->dijkstra_vertices(source, destination);
      // Check shortest path
      REQUIRE(graph->path_cost(path) == Approx(3.7).epsilon(Tolerance));
    }

    // Check remove edge
    SECTION("Check remove edge") {
      // Remove edge (3, 1)
      graph->remove_edge(3, 1);
      // Run Dijkstra Priority Queue
      const auto path = graph->dijkstra_vertices(source, destination);
      // Check shortest path
      REQUIRE(graph->path_cost(path) == Approx(7.2).epsilon(Tolerance));

      // Remove edge (2, 4)
      graph->remove_edge(2, 4);
      // Run Dijkstra Priority Queue
      const auto new_path = graph->dijkstra_vertices(source, destination);
      // Check shortest path
      REQUIRE(graph->path_cost(new_path) == Approx(9.1).epsilon(Tolerance));
    }
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
    const auto path = graph->dijkstra_vertices(source, destination);

    /*
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
    */
  }

  SECTION("Test SP in directed graph from file") {
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

    std::vector<abm::graph::vertex_t> vertices{
        1020,  39432, 794,   35644, 793,   792,   791,   790,   37527, 3802,
        3781,  37526, 24098, 31169, 30417, 31166, 30410, 31160, 30401, 30400,
        30389, 30396, 30387, 30395, 31168, 41904, 41863, 34668, 24734, 39642,
        25893, 3965,  3961,  3912,  14938, 3826,  3999,  3728,  12499, 19790,
        24057, 24018, 33932, 41303, 24741, 19849, 5376,  19794, 19795, 36365,
        19942, 19944, 19801, 19803, 19806, 19804, 19808, 20017, 22901, 33572,
        33536, 33556, 33502, 29020, 12487, 21616, 12483, 40398, 12481, 33509,
        9117,  4864,  28397, 33524, 40399, 33559, 9125,  40396, 33516, 9394,
        41476, 2497,  34451, 20062, 34447, 34443, 8366,  8334,  21305, 14905,
        21299, 512,   32121, 30412, 40884, 38447, 8752,  27099, 39857, 11226,
        20010, 31373, 29021, 13582, 130,   41772, 134,   37768, 486,   26861,
        26756, 3944,  3929,  3928,  26903, 1725,  3935,  4174,  3914,  38612,
        35459, 27994, 3892,  13726, 34186, 40267, 34187, 37158, 33669, 34193,
        30608, 34785, 37289, 37290, 24829, 24771, 34505, 27850, 34507, 27851,
        1897,  10798, 25,    19,    11913, 11907, 10899, 10897, 10953, 33466,
        20};

    const auto route_vertices = graph->dijkstra(source, destination);
    for (unsigned i = 0; i < vertices.size(); ++i)
      REQUIRE(route_vertices.at(i) == vertices.at(i));

    const auto path = graph->dijkstra_vertices(source, destination);

    // Check # of edges
    REQUIRE(path.size() == 150);
    REQUIRE(graph->path_cost(path) ==
            Approx(12409.660000000002).epsilon(Tolerance));

    // Check shortest path
    const auto edges = graph->dijkstra_edges(source, destination);

    // Check distances
    REQUIRE(edges.size() == 150);
    // REQUIRE(graph->path_cost(edges) ==
    // Approx(12409.660000000002).epsilon(Tolerance));

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

    SECTION("Test non-existant file") {
      // Create graph object
      auto graph = std::make_unique<abm::Graph>(directed);
      // Read MatrixMarket file
      const std::string filename = "nofile.mtx";
      // Read file should fail
      REQUIRE(graph->read_graph_matrix_market(filename) == false);
    }
  }

  SECTION("Test SSSP in OSM directed graph from file") {
    // Set graph properties
    const bool directed = true;

    // Create graph object
    auto graph = std::make_unique<abm::Graph>(directed);
    // Read MatrixMarket file
    const std::string filename = path + "osm/edges.csv";
    // Read file
    REQUIRE(graph->read_graph_osm(filename) == true);

    // Run Dijkstra Priority Queue
    abm::graph::vertex_t source = 65358739;
    abm::graph::vertex_t destination = 2304626340;

    const auto path = graph->dijkstra_edges(source, destination);

    // Check # of edges
    REQUIRE(path.size() == 58);

    SECTION("Test non-existant file") {
      // Create graph object
      auto graph = std::make_unique<abm::Graph>(directed);
      // Read MatrixMarket file
      const std::string filename = "osm/nofile.csv";
      // Read file should fail
      REQUIRE(graph->read_graph_osm(filename) == false);
    }
  }
}
