#include <memory>
#include <vector>

#include "config.h"
#include "graph.h"
#include "agent.h"

// Compute path
void abm::Agent::compute_agent_path () {
  const auto path = graph_ -> dijkstra(current_node_, destination_);
  this -> path_ = path;
};

// Print path
void abm::Agent::print_agent_path () {
  if (path_.size()>0) {
    for (auto itr = path_.begin(); itr != path_.end() - 1; ++itr) {
      auto nitr = itr + 1;
      auto weight = graph_ -> get_edge_cost(
          static_cast<abm::graph::vertex_t>(*itr), 
          static_cast<abm::graph::vertex_t>(*nitr));
      std::cout << "(" << *itr << ", " << *nitr << ", w=" << weight << ") ";
    }
  }
  else { std::cout << "path size is zero" << std::endl; }
}

// Move agent
void abm::Agent::move_agent(graph::weight_t time_limit) {
  if (path_.size() > 0) {
    abm::graph::weight_t agent_time = 0; // keep record of cumulated edge cost
    this -> status_ = 1; // agent is now enroute
    for (auto itr = path_.begin(); itr != path_.end() - 1; ++itr) {  
      auto nitr = itr + 1;
      if (agent_time < time_limit) {
        agent_time += graph_->get_edge_cost(
          static_cast<abm::graph::vertex_t>(*itr), 
          static_cast<abm::graph::vertex_t>(*nitr));
        this -> current_node_ = static_cast<abm::graph::vertex_t>(*nitr);
        if (agent_time >= time_limit) {
          this -> current_node_ = static_cast<abm::graph::vertex_t>(*itr);
          this -> path_.erase(path_.begin(), itr);
          break;
        }
      }
    }
  }
  // reach destination
  if (current_node_ == destination_) {
    this -> status_ = 2;
    this -> path_.clear();
  }
};