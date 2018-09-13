#include <chrono>
#include <ctime>
#include <memory>

#include "graph.h"
#include "router.h"

int main(int argc, char** argv) {
  const bool directed = true;
  auto graph = std::make_unique<abm::Graph>(directed);
  if (argc > 1) {
    // Read MatrixMarket file
    const std::string filename = argv[1];
    graph->read_graph_matrix_market(filename);
  } else {
    // Generate a simple graph
    graph->generate_simple_graph();
  }

  /*
  auto router = std::make_unique<abm::Router>(10);
  router->read_od_pairs("../sf-graph-od-50000.csv");

  const auto routes = router->od_pairs();

  unsigned i = 0;
  for (const auto& route : routes) {
    // const auto distances = graph->dijkstra_priority_queue(1, -1);
    const auto sp = graph->dijkstra_priority_queue(route.first, route.second);
    std::cout << ++i << "\n";
    if (i == 1000) break;
  }
  //*/

  auto start = std::chrono::system_clock::now();
  const auto path = graph->dijkstra(1020, 20);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end - start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);

  std::cout << "finished computation at " << std::ctime(&end_time)
            << "elapsed time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   elapsed_seconds)
                   .count()
            << "ms\n";
  std::cout << "Dijkstra PriorityQueue\n";
  for (const auto& vertex : path) std::cout << vertex << "\t";
  /*
  unsigned i = 0;
  for (const auto& distance : sp.distances) {
    std::cout << i << "\t" << distance << "\n";
    ++i;
  }
  */
}
