#ifndef _ABM_AGENT_H_
#define _ABM_AGENT_H_

#include <memory>
#include <vector>

#include "config.h"
#include "graph.h"
#include "agent.h"

namespace abm {

// Agent class with origin, destination, departure time. Computes current node and path.
class Agent {
 public:

  // Construct an agent class
  explicit Agent(const std::shared_ptr<abm::Graph>& graph) : graph_{graph} {};

  // Set and get origin
  graph::vertex_t get_origin() const { return origin_; }
  void set_origin(graph::vertex_t origin) { 
    this->origin_ = origin;
    this->current_node_ = origin;
  }

  // Set and get destination
  graph::vertex_t get_destination() const { return destination_; }
  void set_destination(graph::vertex_t destination) { this->destination_ = destination; }

  // Set and get departure_time
  graph::weight_t get_departure_time() const { return departure_time_; }
  void set_departure_time(graph::weight_t departure_time) { this->departure_time_ = departure_time; }

  // Get current node
  graph::vertex_t get_current_node() const { return current_node_; }

  // Get status
  int get_status() const { return status_;}

  // Compute path
  void compute_agent_path ();
  
  // Get path
  //std::vector<abm::graph::vertex_t> get_agent_path () { return path_; }
  
  // Print path
  void print_agent_path ();

  // Move agent
  void move_agent(graph::weight_t time_limit);

 private:
  // origin vertex
  graph::vertex_t origin_{0};

  // destination vertex
  graph::vertex_t destination_{0};

  // departure time
  graph::weight_t departure_time_{0};

  // departure node
  graph::vertex_t current_node_{origin_};

  // graph
  std::shared_ptr<abm::Graph> graph_;

  // path
  std::vector<abm::graph::vertex_t> path_;

  // status: 0: haven't started routing. 1: en_route. 2: arrived.
  int status_{0}; 
};

}

#endif  // _ABM_AGENT_H_