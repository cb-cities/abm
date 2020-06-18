#include <memory>
#include <unistd.h>
#include <omp.h>

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "graph.h"
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
  graph->read_graph_csv("/work/07427/bingyu/stampede2/abm/osm/tokyo_edges.csv");
  if (myrank == 0) {
    std::cout << "graph has " << graph->nedges() << " edges, " << graph->nvertices() << " vertices" << std::endl;
  }
  // read OD from rank 0 and scatter it into each rank
  auto ag = std::make_unique<abm::Router_hybrid>(graph);
  int nagents = 5000000; // simulation demand
  int subp_agents = 20000; // update graph after 5000 agents are assigned
  std::vector<std::array<abm::graph::vertex_t, 3>> all_od_pairs;
  if (myrank == 0) {
    std::vector<std::string> demand_input_files = {
      "/work/07427/bingyu/stampede2/abm/osm/tokyo_demands_0.csv", 
      "/work/07427/bingyu/stampede2/abm/osm/tokyo_demands_1.csv", 
      "/work/07427/bingyu/stampede2/abm/osm/tokyo_demands_2.csv" };
    ag->read_timed_od_pairs(demand_input_files, nagents);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double time_start = MPI_Wtime();

  // ag->make_timed_od_map(0, npagents, nproc, myrank);

  for (int hour=3; hour!=20; ++hour) {
    for (int quarter=0; quarter!=4; ++quarter){
      ag->quarter_router(hour, quarter, subp_agents, myrank, nproc);
    }
  } 
  std::cout << " finish all time steps " << myrank << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  double time_end = MPI_Wtime();

  MPI_Finalize();
  if (myrank == 0) {
    std::cout << "run time is " << time_end - time_start << std::endl;
  }
  return 0;
}
