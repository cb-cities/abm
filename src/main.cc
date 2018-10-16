#include "mpi.h"
#include "omp.h"

#include <chrono>
#include <ctime>
#include <memory>

#include "graph.h"
#include "router.h"

int main(int argc, char** argv) {
  // Initialise MPI
  MPI_Init(&argc, &argv);
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  const bool directed = true;
  auto graph = std::make_shared<abm::Graph>(directed);

  std::string od_file;
  if (argc == 3) {
    // Read MatrixMarket file
    const std::string filename = argv[1];
    od_file = argv[2];
    graph->read_graph_matrix_market(filename);
  } else {
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // On MPI rank 0 create a router and fetch all OD pairs
  std::shared_ptr<abm::Router> router;
  if (mpi_rank == 0) {
    router = std::make_shared<abm::Router>(graph);
    router->read_od_pairs(od_file, 5000);
  }

  const auto all_paths = router->compute_routes(mpi_rank);

  if (mpi_rank == 0)
    std::cout << "Collected paths: " << all_paths.size() << std::endl;

  MPI_Finalize();
}
