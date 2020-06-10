// http://mirror.olnevhost.net/pub/lam/tutorials/one-step/ezstart.php
#ifndef _ABM_ROUTER_HYBRID_H_
#define _ABM_ROUTER_HYBRID_H_

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "config.h"
#include "graph.h"

namespace abm {

class Router_hybrid {
  public:
    explicit Router_hybrid(const std::shared_ptr<abm::Graph>& graph) : graph_{graph} {};
    bool read_timed_od_pairs(const std::string& filename,
                     int nagents = std::numeric_limits<int>::max());
    void make_timed_od_map(bool print_od_map, int npagents, int nproc, int myrank);
    // void make_edge_vol_map (std::vector<abm::graph::vertex_t>& path);
    // void output_edge_vol_map (const std::string& output_filename);
    
    // void send_od_to_workers(int nproc, int myrank);
    // void master(int nproc, int myrank, int nagents);
    // void worker(int nproc, int myrank, int nagents);
    void router(int hour, int quarter, int npagents, int myrank, int nproc);
  
  private:
   std::shared_ptr<abm::Graph> graph_;
   std::vector<std::array<abm::graph::vertex_t, 4>> input_ods_;
   std::vector<abm::graph::vertex_t> input_ods_test_;
   std::map<int, std::map<int, std::vector<std::array<abm::graph::vertex_t, 2>>>> partial_timed_od_pairs_;
   std::map<abm::graph::vertex_t, abm::graph::vertex_t> edge_vol_;
};
}

#endif  // _ABM_ROUTER_HYBRID_H_