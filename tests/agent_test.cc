#include <memory>

#include "catch.hpp"

#include "agent.h"

// Check agent class
TEST_CASE("Agent class is checked", "[agent]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  SECTION("Test agent initialization, parameter assignment and routing") {
    // Set graph properties
    const bool directed = true;
    // Create simple graph object for the test agent
    auto graph = std::make_shared<abm::Graph>(directed);
    graph->generate_simple_graph();

    // Create test agent
    auto test_agent = std::make_unique<abm::Agent>(graph);

    // Parameters to be assigned to the test agent
    abm::graph::vertex_t agent_id = 100;
    abm::graph::vertex_t origin = 1;
    abm::graph::vertex_t destination = 3;
    abm::graph::weight_t departure_time = 0;

    // Assign agent parameters
    test_agent->set_id(agent_id);
    test_agent->set_origin(origin);
    test_agent->set_destination(destination);
    test_agent->set_departure_time(departure_time);

    // Check assigned values
    REQUIRE(test_agent->get_id() == agent_id);
    REQUIRE(test_agent->get_origin() == origin);
    REQUIRE(test_agent->get_destination() == destination);
    REQUIRE(test_agent->get_departure_time() == departure_time);

    // Check current nodes is initialized to the origin
    REQUIRE(test_agent->get_current_node() == origin);

    // Test agent routing
    // Compute shortest path between curret node and destination node
    REQUIRE(test_agent->get_current_node() == origin);  // At first, the current
                                                        // node is the origin
    REQUIRE(test_agent->get_status() == 0);  // At first, the agent status is 0
                                             // "waiting to depart"
    test_agent->compute_agent_path();
    // test_agent -> print_agent_path(); // (1, 2, w=1.5) (2, 4, w=5.5) (4, 3,
    // w=0.2)

    // Check agent position after 2 units of cost
    test_agent->move_agent(2);
    REQUIRE(test_agent->get_current_node() == 2);  // After 2 units of cost,
                                                   // agent can reach vertex 2
    REQUIRE(test_agent->get_status() == 1);  // At this time, agent status is 1
                                             // "enroute"

    // Check agent position after 10 units of cost
    test_agent->move_agent(10);
    REQUIRE(test_agent->get_current_node() == 3);  // After 10 unit of cost,
                                                   // agent can reach vertex 3
    REQUIRE(test_agent->get_status() == 2);  // At this time, agent status is 2
                                             // "reached destination"
  }
}