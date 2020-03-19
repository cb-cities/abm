#include "catch.hpp"

#include "agent.h"

// Check agent class
TEST_CASE("Agent class is checked", "[agent]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  SECTION("Test agent initialization") {

    // Parameters to be assigned to the test agent
    abm::graph::vertex_t agent_id = 100;
    abm::graph::vertex_t origin = 1;
    abm::graph::vertex_t destination = 3;
    abm::graph::weight_t departure_time = 0;

    SECTION("Initialization with id") {
      // Create test agent
      auto test_agent = std::make_unique<abm::Agent>(agent_id);

      // Check assigned values
      REQUIRE(test_agent->id() == agent_id);
    }

    SECTION("Initialization with id, origin and destination") {
      // Create test agent
      auto test_agent =
          std::make_unique<abm::Agent>(agent_id, origin, destination);

      // Check assigned values
      REQUIRE(test_agent->id() == agent_id);
      REQUIRE(test_agent->origin() == origin);
      REQUIRE(test_agent->destination() == destination);
    }
  }
}
