#include "mpi.h"
#include "omp.h"

#include <chrono>
#include <ctime>
#include <memory>
#include <numeric>

#include "graph.h"
#include "router.h"

int main(int argc, char** argv) {
  // Initialise MPI
  MPI_Init(&argc, &argv);
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  const bool directed = true;
  auto graph = std::make_unique<abm::Graph>(directed);
  if (argc > 1) {
    // Read MatrixMarket file
    const std::string filename = argv[1];
    graph->read_graph_matrix_market(filename);
  } else {
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  std::vector<std::array<abm::graph::vertex_t, 2>> all_routes;

  if (mpi_rank == 0) {
    auto router = std::make_unique<abm::Router>(10);
    router->read_od_pairs("../sf-od-50000.csv");

    all_routes = router->od_pairs();
    all_routes.resize(5000);
  }

  MPI_Datatype pair_t;
  MPI_Type_vector(2, 1, 1, MPI_INT, &pair_t);
  MPI_Type_commit(&pair_t);

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int chunk_size = all_routes.size() / mpi_size;
  MPI_Bcast(&chunk_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  std::vector<std::array<abm::graph::vertex_t, 2>> routes(chunk_size);

  //  std::cout << "routes size\t" << mpi_rank << '\t' << routes.size() << '\n';
  // std::cout << "all_routes size\t" << mpi_rank << '\t' << all_routes.size()
  // << '\n';

  MPI_Scatter(all_routes.data(), chunk_size, pair_t, routes.data(),
              routes.size(), pair_t, 0, MPI_COMM_WORLD);
  int chunk_remainder = all_routes.size() % mpi_size;
  if (mpi_rank == 0) {
    routes.insert(routes.begin(), all_routes.end() - chunk_remainder,
                  all_routes.end());
  }

  // Paths (vector of edges)
  std::vector<std::array<abm::graph::vertex_t, 2>> paths;
  paths.reserve(graph->nedges());
  unsigned i = 0;
#pragma omp parallel for schedule(dynamic)
  for (i = 0; i < routes.size(); ++i) {
    const auto sp = graph->dijkstra(routes[i][0], routes[i][1]);

#pragma omp critical
    paths.insert(std::end(paths), std::begin(sp), std::end(sp));
  }

#if 1
  unsigned path_size = paths.size();
  std::vector<int> paths_sizes(mpi_size);
  MPI_Gather(&path_size, 1, MPI_INT, paths_sizes.data(), 1, MPI_INT, 0,
             MPI_COMM_WORLD);

  std::vector<std::array<abm::graph::vertex_t, 2>> all_paths(
      std::accumulate(paths_sizes.begin(), paths_sizes.end(), 0));
  auto paths_scan = paths_sizes;
  std::partial_sum(paths_scan.begin(), paths_scan.end(), paths_scan.begin());
  paths_scan.insert(paths_scan.begin(), 0);

  MPI_Gatherv(paths.data(), paths.size(), pair_t, all_paths.data(),
              paths_sizes.data(), paths_scan.data(), pair_t, 0, MPI_COMM_WORLD);
#else
  auto all_paths = paths;
#endif

  if (mpi_rank == 0) {
    std::cout << all_paths.size() << std::endl;
  }

  MPI_Finalize();
}
