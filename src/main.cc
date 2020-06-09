#include <memory>
#include <unistd.h>
#include <omp.h>

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "graph.h"
// #include "router.h"
#include "router_hybrid.h"

int main(int argc, char** argv) {

  // MPI
  int nproc, myrank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  // OpenMP export OMP_NUM_TRHEADS=36
  int s = omp_get_max_threads();
  if (myrank == 0) {
    std::cout << "MPI_Comm_size " << nproc << "; Available OpenMP threads " << s << std::endl;
  }

  // read graph
  const bool directed = true;
  auto graph = std::make_shared<abm::Graph>(directed);
  graph->read_graph_csv("../osm/tokyo_edges.csv");
  if (myrank == 0) {
    std::cout << "graph has " << graph->nedges() << " edges, " << graph->nvertices() << " vertices" << std::endl;
  }
  // read OD from rank 0 and scatter it into each rank
  auto ag = std::make_unique<abm::Router_hybrid>(graph);
  int nagents = 5000;
  int npagents = 1000;
  std::vector<std::array<abm::graph::vertex_t, 3>> all_od_pairs;
  if (myrank == 0) {
    ag->read_timed_od_pairs("../osm/tokyo_demands_0.csv", nagents);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  // int sub_avg=0;
  // int *sub_avgs = (int *)malloc(sizeof(int) * nproc);
  // MPI_Allgather(&sub_avg, 1, MPI_INT, sub_avgs, 1, MPI_INT, MPI_COMM_WORLD);
  double time_start = MPI_Wtime();

  ag->make_timed_od_map(0, npagents, nproc, myrank);

  for (int hour=3; hour!=4; ++hour) {
    for (int quarter=0; quarter!=2; ++quarter){
      ag->router(hour, quarter, npagents, myrank);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  } 
  std::cout << " finish all time steps " << myrank << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  // std::cout << " pass barrier " << myrank << std::endl;
  // int sub_avg2=0;
  // int *sub_avgs2 = (int *)malloc(sizeof(int) * nproc);
  // MPI_Allgather(&sub_avg2, 1, MPI_INT, sub_avgs2, 1, MPI_INT, MPI_COMM_WORLD);
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
