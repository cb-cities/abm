#ifndef _ABM_ROUTER_H_
#define _ABM_ROUTER_H_

#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <vector>

#include "csv.h"

#include "config.h"
#include "graph.h"
#include "mpi_wrapper.h"

namespace abm {
//! \brief Router of agents
class Router {
 public:
  //! Constructor with number of agents
  explicit Router(const std::shared_ptr<abm::Graph>& graph) : graph_{graph} {};

  //! Get OD pairs
  //! \param[in] filename containing origin destination pairs
  bool read_od_pairs(const std::string& filename,
                     int nagents = std::numeric_limits<int>::max());

  //! Return OD pairs
  std::vector<std::array<abm::graph::vertex_t, 2>> od_pairs() const {
    return this->all_od_pairs_;
  }

  //! Compute routes
  std::vector<std::array<abm::graph::vertex_t, 2>> compute_routes(int mpi_rank,
                                                                  int mpi_size);

 private:
  //! Graph
  std::shared_ptr<abm::Graph> graph_;
  //! All OD pairs
  std::vector<std::array<abm::graph::vertex_t, 2>> all_od_pairs_;
  //! All paths
  std::vector<std::array<abm::graph::vertex_t, 2>> all_paths_;
  //! All paths indices
  std::vector<std::array<abm::graph::vertex_t, 3>> all_paths_idx_;
};
}  // namespace abm
#endif  // _ABM_ROUTER_H_
