#include <memory>

#include "catch.hpp"
#include "mpi.h"

#include "router.h"

// Check Router class
TEST_CASE("Router class is checked", "[router]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  // Graph file
  const std::string path = "../";
  const std::string graph_file = path + "network.mtx";
  const std::string od_pairs = path + "sf-od-50k.csv";

  // Test read od pairs
  SECTION("Test read OD pairs") {
    // Set graph properties
    const bool directed = true;

    // Create graph object
    auto graph = std::make_shared<abm::Graph>(directed);
    // Create a simple example graph
    graph->read_graph_matrix_market(graph_file);

    auto router = std::make_unique<abm::Router>(graph);
    REQUIRE(router->read_od_pairs(od_pairs) == true);

    auto all_routes = router->od_pairs();
    REQUIRE(all_routes.size() == 49946);
    all_routes.resize(5000);
    REQUIRE(all_routes.size() == 5000);

    SECTION("Test non-existant file") {
      // Read OD pairs file
      const std::string filename = "nofile.csv";
      // Read file should fail
      REQUIRE(router->read_od_pairs(filename) == false);
    }

    SECTION("Test non-existant column in file") {
      // Read OD pairs file
      const std::string filename = path + "incorrect-od.csv";
      // Read file should fail
      REQUIRE(router->read_od_pairs(filename) == false);
    }

    SECTION("Compute routes") {
      // MPI ranks
      int mpi_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

      // Get number of MPI ranks
      int mpi_size;
      MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

      REQUIRE(router->read_od_pairs(od_pairs, 50) == true);
      const auto all_paths = router->compute_routes(mpi_rank, mpi_size);
      REQUIRE(all_paths.size() == 3564);
    }
  }
}
