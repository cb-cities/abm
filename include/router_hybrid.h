// http://mirror.olnevhost.net/pub/lam/tutorials/one-step/ezstart.php
#ifndef _ABM_ROUTER_HYBRID_H_
#define _ABM_ROUTER_HYBRID_H_

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "config.h"
#include "graph.h"

namespace abm {

struct Volume_and_Residual
{
  std::vector<std::array<abm::graph::vertex_t, 2>> Volume_Vector;
  std::vector<std::array<abm::graph::vertex_t, 2>> Residual_OD_Vector;
};

// struct Route_and_Weight
// {
//   std::vector<graph::vertex_t> Route_Vector;
//   double Weight;
// };

class Router_hybrid {
  public:
    explicit Router_hybrid(const std::shared_ptr<abm::Graph>& graph) : graph_{graph} {};
    // read od from csv file
    bool read_timed_od_pairs(const std::vector<std::string>& filename, int nagents = std::numeric_limits<int>::max());
    // quarterly routing
    void quarter_router (int hour, int quarter, int npagents, int myrank, int nproc, std::string output_suffix);
  
  private:
    // convert input OD (csv) into map with hour and quarter as the keys
    std::map<int, std::map<int, std::vector<std::array<abm::graph::vertex_t, 2>>>> make_timed_od_map (
                    bool print_od_map, int nagents, std::vector<std::array<abm::graph::vertex_t, 4>>& od_inputs);
    // called by quarter router for substep (OpenMP) path calculation
    Volume_and_Residual substep_router (
                    int myrank, std::vector<std::array<abm::graph::vertex_t, 2>>& partial_ods);
    // save quarterly edge volume to file
    void output_edge_vol_map (const std::string& output_filename);

    // private variables
    std::shared_ptr<abm::Graph> graph_;
    std::map<int, std::map<int, std::vector<std::array<abm::graph::vertex_t, 2>>>> input_ods_;
    std::map<abm::graph::vertex_t, abm::graph::vertex_t> edge_vol_;
};
}

#endif  // _ABM_ROUTER_HYBRID_H_