#include "mpi.h"
#include "omp.h"

#include <chrono>
#include <ctime>
#include <memory>
#include <numeric>

#include "graph.h"
#include "router.h"

template <typename T, size_t N>
std::vector<std::array<T, N>> thebestgather(std::vector<std::array<T, N>> const &x) {
  MPI_Datatype arr_t;
  MPI_Type_vector(N, 1, 1, MPI_INT, &arr_t);
  MPI_Type_commit(&arr_t);

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  unsigned x_size = x.size();
  std::vector<int> x_sizes(mpi_size);
  MPI_Gather(&x_size, 1, MPI_INT, x_sizes.data(), 1, MPI_INT, 0,
             MPI_COMM_WORLD);

  std::vector<std::array<abm::graph::vertex_t, N>> all_x(
      std::accumulate(x_sizes.begin(), x_sizes.end(), 0));
  auto x_scan = x_sizes;
  std::partial_sum(x_scan.begin(), x_scan.end(), x_scan.begin());
  x_scan.insert(x_scan.begin(), 0);

  MPI_Gatherv(x.data(), x.size(), arr_t, all_x.data(),
              x_sizes.data(), x_scan.data(), arr_t, 0, MPI_COMM_WORLD);

  MPI_Type_free(&arr_t);

  return all_x;
}

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

  // Vector of all origin destination pairs from router
  std::vector<std::array<abm::graph::vertex_t, 2>> all_od_pairs;

  // On MPI rank 0 create a router and fetch all OD pairs
  if (mpi_rank == 0) {
    auto router = std::make_unique<abm::Router>(5000, graph);
    router->read_od_pairs(od_file);

    all_od_pairs = router->od_pairs();
    //all_od_pairs.resize(5000);
  }

  // Create MPI pair type
  MPI_Datatype pair_t;
  MPI_Type_vector(2, 1, 1, MPI_INT, &pair_t);
  MPI_Type_commit(&pair_t);

  // Get number of MPI ranks
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // Calculate chunk size to split router
  int chunk_size = all_od_pairs.size() / mpi_size;
  MPI_Bcast(&chunk_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  std::vector<std::array<abm::graph::vertex_t, 2>> od_pairs(chunk_size);

  // Send route chunks to different compute nodes
  MPI_Scatter(all_od_pairs.data(), chunk_size, pair_t, od_pairs.data(),
              od_pairs.size(), pair_t, 0, MPI_COMM_WORLD);
  // Calculate the remaining chunk of od_pairs and add to rank 0
  int chunk_remainder = all_od_pairs.size() % mpi_size;
  if (mpi_rank == 0) {
    od_pairs.insert(od_pairs.begin(), all_od_pairs.end() - chunk_remainder,
                    all_od_pairs.end());
  }

  // Paths (vector of edges)
  std::vector<std::array<abm::graph::vertex_t, 2>> paths;
  paths.reserve(graph->nedges());

  std::vector<std::array<abm::graph::vertex_t, 3>> paths_idx;
  paths_idx.reserve(od_pairs.size());

#pragma omp parallel for schedule(dynamic)
  for (unsigned i = 0; i < od_pairs.size(); ++i) {
    const auto sp = graph->dijkstra(od_pairs[i][0], od_pairs[i][1]);

#pragma omp critical
    {
      paths_idx.emplace_back(std::array<abm::graph::vertex_t, 3>(
            {i, (int)paths.size(), (int)sp.size()}));
      paths.insert(std::end(paths), std::begin(sp), std::end(sp));
    }
  }

//// Get size of routes from each MPI rank
//unsigned path_size = paths.size();
//std::vector<int> paths_sizes(mpi_size);
//MPI_Gather(&path_size, 1, MPI_INT, paths_sizes.data(), 1, MPI_INT, 0,
//           MPI_COMM_WORLD);

//// Create a vector of cumulative size of path from each MPI rank
//std::vector<std::array<abm::graph::vertex_t, 2>> all_paths(
//    std::accumulate(paths_sizes.begin(), paths_sizes.end(), 0));
//auto paths_scan = paths_sizes;
//std::partial_sum(paths_scan.begin(), paths_scan.end(), paths_scan.begin());
//paths_scan.insert(paths_scan.begin(), 0);

//// Accumulate paths from all MPI ranks
//MPI_Gatherv(paths.data(), paths.size(), pair_t, all_paths.data(),
//            paths_sizes.data(), paths_scan.data(), pair_t, 0, MPI_COMM_WORLD);

  const auto all_paths = thebestgather(paths);
  auto const all_paths_idx = thebestgather(paths_idx);

  if (mpi_rank == 0)
    std::cout << "Collected paths: " << all_paths.size() << std::endl;

  MPI_Finalize();
}
