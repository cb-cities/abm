#ifndef _ABM_AGENT_H_
#define _ABM_AGENT_H_

#include <memory>
#include <vector>

#include "config.h"
#include "graph.h"
#include "agent.h"

namespace abm {

//! Agent class that has origin, destination, departure time, current node and compute path from current node to destination.
class Agent {
 public:
  //! Construct an agent class
  //! \param[in] graph Pass a reference of the graph that the agent does routing on
  explicit Agent(const std::shared_ptr<abm::Graph>& graph) : graph_{graph} {};

  //! Set origin and initialize current node to origin
  //! \param[in] origin Origin vertex assigned to agent
  void set_origin(graph::vertex_t origin) { 
    this->origin_ = origin;
    this->current_node_ = origin;
  }
  //! Get origin
  graph::vertex_t get_origin() const { return origin_; }

  //! Set destination
  //! \param[in] destination Destination vertex assigned to agent
  void set_destination(graph::vertex_t destination) { this->destination_ = destination; }
  //! Get destination
  graph::vertex_t get_destination() const { return destination_; }

  //! Set departure_time
  //! \param[in] departure_time Departure time of an agent
  void set_departure_time(graph::weight_t departure_time) { this->departure_time_ = departure_time; }
  //! Get departure time
  graph::weight_t get_departure_time() const { return departure_time_; }

  // Get current node
  graph::vertex_t get_current_node() const { return current_node_; }

  //! Get status
  //! \retval Status of the agent. 0: waiting to depart; 1: enroute; 2: reached destination.
  int get_status() const { return status_;}

  //! Compute path based on current node and destination
  void compute_agent_path ();
  
  //! Get path
  //! \retval Return path from current node to destination as a list of vertices
  std::vector<abm::graph::vertex_t> get_agent_path () { return path_; }
  
  //! Print path
  void print_agent_path ();

  //! Move agent, update current_node_, path_ and status_ of the agent object
  //! \param[in] time_limit The time that sets the limit on the total weights of a sub route that an agent can cover in the current time step 
  void move_agent(graph::weight_t time_limit);

 private:
  //! Agent origin vertex
  graph::vertex_t origin_{0};

  //! Agent destination vertex
  graph::vertex_t destination_{0};

  //! Agent departure time
  graph::weight_t departure_time_{0};

  //! Agent current node
  graph::vertex_t current_node_{origin_};

  //! Graph that agent travels on
  std::shared_ptr<abm::Graph> graph_;

  //! Path from current node to destination
  std::vector<abm::graph::vertex_t> path_;

  //! status: 0: haven't started routing. 1: en_route. 2: arrived.
  int status_{0}; 
};

} // namespace abm

#endif  // _ABM_AGENT_H_