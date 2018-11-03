#ifndef ABM_GRAPH_H_
#define ABM_GRAPH_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <sstream>
#include <tuple>
#include <vector>

#include <tsl/robin_map.h>

#include "config.h"

namespace abm {
//! \brief ShortestPath struct to return source, distance and parents
struct ShortestPath {
  //! Get path from source to j using parent array
  //! \param[in] parent Map of vertex to its parent id
  //! \param[in] destination Destination vertex id to get path
  //! \param[in] source Source vertex id to get path (SSSP set as -1)
  //! \retval path Path from source to destination
  std::vector<graph::vertex_t> get_path(graph::vertex_t source,
                                        graph::vertex_t destination) {
    // Create an empty path
    std::vector<graph::vertex_t> path;
    // Iterate until source has been reached
    while (destination != source) {
      destination = parent.at(destination);
      if (destination != source) path.emplace_back(destination);
    }
    // Reverse to arrange from source to destination
    std::reverse(path.begin(), path.end());
    return path;
  }

  //! Source
  graph::vertex_t source;
  //! Destination
  // std::vector<vertex_t> destinations;
  //! Distances
  std::vector<graph::weight_t> distances;
  //! Parent array to store shortest path tree
  std::vector<graph::vertex_t> parent;
};

//! \brief Graph class to store vertices and edge and compute shortest path
//! \details Graph class has Priority Queue Dijkstra algorithm for SSSP
class Graph {
 public:
  //! Edge {{v1, v2}, weight}
  using Edge =
      std::pair<std::pair<graph::vertex_t, graph::vertex_t>, graph::weight_t>;

  //! Construct directed / undirected graph
  //! \param[in] directed Defines if the graph is directed or not
  explicit Graph(bool directed) : directed_{directed} {};

  //! Return number of vertices
  unsigned nvertices() const { return nvertices_; }

  //! Number of edges
  graph::vertex_t nedges() const { return edges_.size(); }

  //! Add edge to graph
  //! \param[in] vertex1 ID of vertex1
  //! \param[in] vertex2 ID of vertex2
  //! \param[in] weight Weight of edge connecting vertex 1 and 2
  void add_edge(graph::vertex_t vertex1, graph::vertex_t vertex2,
                graph::weight_t weight);

  //! Update edge of a graph
  //! \param[in] vertex1 ID of vertex1
  //! \param[in] vertex2 ID of vertex2
  //! \param[in] weight Weight of edge connecting vertex 1 and 2
  void update_edge(graph::vertex_t vertex1, graph::vertex_t vertex2,
                   graph::weight_t weight);

  //! Remove edge from graph
  //! \param[in] vertex1 ID of vertex1
  //! \param[in] vertex2 ID of vertex2
  void remove_edge(graph::vertex_t vertex1, graph::vertex_t vertex2);

  //! Generate a simple graph
  void generate_simple_graph();

  //! Read MatrixMarket graph file format
  //! \param[in] filename Name of input MatrixMarket file
  //! \retval status File read status
  bool read_graph_matrix_market(const std::string& filename);

  //! Compute the shortest path using priority queue
  //! \param[in] source ID of source vertex1
  //! \param[in] destination ID of destination vertex
  //! \retval route_edges Edges of the route from source to destination
  std::vector<std::array<graph::vertex_t, 2>> dijkstra(
      graph::vertex_t source, graph::vertex_t destination);

  //! Compute the shortest path using priority queue
  //! \param[in] source ID of source vertex1
  //! \param[in] destination ID of destination vertex (default is -1 for SSSP)
  //! \retval sp Shortest path and distances
  ShortestPath dijkstra_priority_queue(graph::vertex_t source,
                                       graph::vertex_t destination = -1);
 private:
  //! Assign number of vertices
  //! \param[in] nvertices Number of vertices in graph
  void assign_nvertices(unsigned nvertices) { this->nvertices_ = nvertices; }

  // Directed / undirected
  bool directed_{false};
  // Number of graph vertices
  unsigned nvertices_{std::numeric_limits<unsigned>::max()};
  // Edges
  std::map<std::tuple<graph::vertex_t, graph::vertex_t>, std::shared_ptr<Edge>>
      edges_;
  // adjacency list with iteration over each edge
  tsl::robin_map<graph::vertex_t, std::vector<std::shared_ptr<Edge>>>
      vertex_edges_;
};

}  // namespace abm
#endif  // ABM_GRAPH_H_
