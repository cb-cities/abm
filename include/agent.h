#ifndef _ABM_AGENT_H_
#define _ABM_AGENT_H_

#include <memory>
#include <vector>

#include "agent.h"
#include "config.h"
#include "graph.h"

namespace abm {

//! Agent class that has id, origin, destination, departure time, current node
//! and compute path from current node to destination.
class Agent {
 public:
  //! Construct an agent class from either an id, or id, origin and destination
  //! \param[in] id Agent id
  //! \param[in] origin Agent origin
  //! \param[in] destination Agent destination
  explicit Agent(const graph::vertex_t& id) : id_{id} {};
  explicit Agent(const graph::vertex_t& id, const graph::vertex_t& origin,
                 const graph::vertex_t& destination)
      : id_{id}, origin_{origin}, destination_{destination} {};

  //! Get and set id
  graph::vertex_t id() const { return id_; }

  //! Get and set origin
  graph::vertex_t origin() const { return origin_; }

  //! Get and set destination
  graph::vertex_t destination() const { return destination_; }

  //! Get and set departure time
  double departure_time() const { return departure_time_; }
  void departure_time(const double departure_time) {
    departure_time_ = departure_time;
  }

  // Get current node
  graph::vertex_t current_node() const { return current_node_; }

  //! status: 0: haven't started routing. 1: en_route. 2: arrived.
  enum class Status { PREDEPART, ENROUTE, ARRIVE };
  //! Get status
  //! \retval Status of the agent. 0: waiting to depart; 1: enroute; 2: reached
  //! destination.
  Status status() const { return status_; }

 private:
  //! Agent id
  graph::vertex_t id_{std::numeric_limits<graph::vertex_t>::max()};

  //! Agent origin vertex
  graph::vertex_t origin_{std::numeric_limits<graph::vertex_t>::max()};

  //! Agent destination vertex
  graph::vertex_t destination_{std::numeric_limits<graph::vertex_t>::max()};

  //! Agent departure time
  graph::weight_t departure_time_{std::numeric_limits<graph::weight_t>::max()};

  //! Agent current node
  graph::vertex_t current_node_{std::numeric_limits<graph::vertex_t>::max()};

  //! status: 0: haven't started routing. 1: en_route. 2: arrived.
  Status status_;
};

}  // namespace abm

#endif  // _ABM_AGENT_H_
