#include <memory>

#include "mpi.h"

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
  auto router = std::make_unique<abm::Router>(graph);
  router->read_od_pairs(od_file, 5000);

  // Get number of MPI ranks
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  const auto all_paths = router->compute_routes(mpi_rank, mpi_size);

  if (mpi_rank == 0)
    std::cout << "Collected paths: " << all_paths.size() << std::endl;

  MPI_Finalize();
}
