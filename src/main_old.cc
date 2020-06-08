#include <memory>
#include <unistd.h>

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "graph.h"
// #include "router.h"
#include "router_mpi.h"

int main(int argc, char** argv) {

  // Get number of threads and this rank
  int nproc, myrank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0) {
    std::cout << "MPI_Comm_size " << nproc << std::endl;
  }

  const bool directed = true;
  auto graph = std::make_shared<abm::Graph>(directed);
  graph->read_graph_csv("../osm/tokyo_edges.csv");
  if (myrank == 0) {
    std::cout << "graph has " << graph->nedges() << " edges, " << graph->nvertices() << " vertices" << std::endl;
  }
  auto ag = std::make_unique<abm::Router_mpi>(graph);
  int nagents = 5000;

  MPI_Barrier(MPI_COMM_WORLD);
  double time_start = MPI_Wtime();
  if (myrank == 0){
    ag->read_timed_od_pairs("../osm/tokyo_demands_0.csv", nagents);
    ag->master(nproc, myrank, nagents);
  } else {
    ag->worker(nproc, myrank, nagents);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  double time_end = MPI_Wtime();

  MPI_Finalize();
  if (myrank == 0) {
    std::cout << "run time is " << time_end - time_start << std::endl;
  }
  return 0;
}

// int main_2(int argc, char** argv) {

//   int mpi_rank = 0;
//   int mpi_size = 1;

// #ifdef USE_MPI
//   // Initialise MPI
//   MPI_Init(&argc, &argv);
//   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
//   MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
// #endif

//   const bool directed = true;
//   auto graph = std::make_shared<abm::Graph>(directed);

//   std::string od_file;
//   if (argc == 3) {
//     // Read MatrixMarket file
//     const std::string filename = argv[1];
//     od_file = argv[2];
//     graph->read_graph_osm(filename);
//   } else {
// #ifdef USE_MPI
//     MPI_Abort(MPI_COMM_WORLD, 1);
// #endif
//   }

//   // On MPI rank 0 create a router and fetch all OD pairs
//   auto router = std::make_unique<abm::Router>(graph);
//   router->read_od_pairs(od_file);

//   const auto all_paths = router->compute_routes(mpi_rank, mpi_size);

//   if (mpi_rank == 0)
//     std::cout << "Collected paths: " << all_paths.size() << std::endl;

// #ifdef USE_MPI
//   MPI_Finalize();
// #endif
// }
