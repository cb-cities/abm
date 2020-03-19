#include <vector>

#include "agent.h"
#include "config.h"
#include "graph.h"

// Compute path from the current_node_ to destination_
void abm::Agent::compute_path() {
  // dijkstra function returns a list of vertices
  const auto path = graph_->dijkstra(current_node_, destination_);
  this->path_ = path;
};

// Print path in the format of (start_of_edge, end_of_edge, w=weight_of_edge),
// (s, e, w=), ...
void abm::Agent::print_path() {
  if (path_.size() > 0) {
    for (auto itr = path_.begin(); itr != path_.end() - 1; ++itr) {
      auto nitr = itr + 1;
      auto weight =
          graph_->get_edge_cost(static_cast<abm::graph::vertex_t>(*itr),
                                static_cast<abm::graph::vertex_t>(*nitr));
      std::cout << "(" << *itr << ", " << *nitr << ", w=" << weight << ") ";
    }
  } else {
    std::cout << "path size is zero" << std::endl;
  }
}

// Move agent
void abm::Agent::move_agent(graph::weight_t time_limit) {
  if (path_.size() > 0) {
    // Keep record of cumulated edge cost (cumulative travel time) not to exceed
    // the time_limit
    abm::graph::weight_t agent_time = 0;
    // Change/Keep agent status to 1 "enroute"
    this->status_ = 1;
    for (auto itr = path_.begin(); itr != path_.end() - 1; ++itr) {
      auto nitr = itr + 1;
      if (agent_time < time_limit) {
        // Add the cost of the next edge
        agent_time +=
            graph_->get_edge_cost(static_cast<abm::graph::vertex_t>(*itr),
                                  static_cast<abm::graph::vertex_t>(*nitr));
        this->current_node_ = static_cast<abm::graph::vertex_t>(*nitr);
        if (agent_time >= time_limit) {
          // If cumulative cost exceed the time_limit, make agent stops there
          // (even though it does not finishing the current link)
          this->current_node_ = static_cast<abm::graph::vertex_t>(*itr);
          // Update path, set the beginning of the path to the current node
          this->path_.erase(path_.begin(), itr);
          break;
        }
      }
    }
  }
  // When agent reachs destination, update its status to 2 "arrival" and delete
  // the path
  if (current_node_ == destination_) {
    this->status_ = 2;
    this->path_.clear();
  }
};