#ifndef _ABM_ROUTER_H_
#define _ABM_ROUTER_H_

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "csv.h"

#include "config.h"
#include "graph.h"

namespace abm {
//! \brief Router of agents
class Router {
 public:
  //! Constructor with number of agents
  explicit Router(unsigned nagents, const std::shared_ptr<abm::Graph>& graph)
      : nagents_{nagents}, graph_{graph} {};

  //! Get OD pairs
  //! \param[in] filename containing origin destination pairs
  bool read_od_pairs(const std::string& filename);

  //! Return OD pairs
  std::vector<std::array<abm::graph::vertex_t, 2>> od_pairs() const {
    return this->od_pairs_;
  }

 private:
  //! Number of agents
  unsigned nagents_;
  //! Graph
  std::shared_ptr<abm::Graph> graph_;
  //! OD pairs
  std::vector<std::array<abm::graph::vertex_t, 2>> od_pairs_;
};
}  // namespace abm
#endif  // _ABM_ROUTER_H_
