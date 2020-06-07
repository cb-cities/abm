// http://mirror.olnevhost.net/pub/lam/tutorials/one-step/ezstart.php
#ifndef _ABM_ROUTER_MPI_H_
#define _ABM_ROUTER_MPI_H_

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "config.h"
#include "graph.h"

namespace abm {

class Router_mpi {
  public:
    explicit Router_mpi(const std::shared_ptr<abm::Graph>& graph) : graph_{graph} {};
    bool read_od_pairs(const std::string& filename,
                     int nagents = std::numeric_limits<int>::max());
    bool read_timed_od_pairs(const std::string& filename,
                     int nagents = std::numeric_limits<int>::max());
    void make_timed_od_map(std::vector<std::array<abm::graph::vertex_t, 3>>& timed_od_pairs, 
                     bool print_od_map, int nagents);
    void make_edge_vol_map (std::vector<abm::graph::vertex_t>& path);
    void output_edge_vol_map (const std::string& output_filename);
    
    void send_od_to_workers(int nproc, int myrank);
    void master(int nproc, int myrank, int nagents);
    void worker(int nproc, int myrank, int nagents);
  
  private:
   std::shared_ptr<abm::Graph> graph_;
   std::vector<std::array<abm::graph::vertex_t, 2>> all_od_pairs_;
   std::map<int, std::vector<std::array<abm::graph::vertex_t, 2>>> all_timed_od_pairs_;
   std::map<abm::graph::vertex_t, abm::graph::vertex_t> edge_vol_;
};
}

#endif  // _ABM_ROUTER_MPI_H_