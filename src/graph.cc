#include "graph.h"

// Add edge
inline void abm::Graph::add_edge(abm::graph::vertex_t vertex1,
                                 abm::graph::vertex_t vertex2,
                                 abm::graph::weight_t weight = 1) {

  if (!this->directed_)
    if (vertex1 > vertex2) std::swap(vertex1, vertex2);

  // Create an edge
  auto edge = std::make_shared<Graph::Edge>(
      std::make_pair(std::make_pair(vertex1, vertex2), weight));
  edges_[std::make_tuple(vertex1, vertex2)] = edge;

  // Vertex 1
  auto vertex1_edges = vertex_edges_[vertex1];
  vertex1_edges.emplace_back(edge);
  vertex_edges_[vertex1] =
      std::vector<std::shared_ptr<Graph::Edge>>(vertex1_edges);

  if (!this->directed_) {
    // Vertex 2
    auto vertex2_edges = vertex_edges_[vertex2];
    vertex2_edges.emplace_back(edge);
    vertex_edges_[vertex2] =
        std::vector<std::shared_ptr<Graph::Edge>>(vertex2_edges);
  }
}

// Update edge
void abm::Graph::update_edge(abm::graph::vertex_t vertex1,
                             abm::graph::vertex_t vertex2,
                             abm::graph::weight_t weight) {
  // Get pointer to specified edge connecting vertex 1 and 2
  auto edge = edges_.at(std::make_tuple(vertex1, vertex2));
  // Update edge weight
  edge->second = weight;
}

// Remove edge
void abm::Graph::remove_edge(abm::graph::vertex_t vertex1,
                             abm::graph::vertex_t vertex2) {
  auto edge = edges_[std::make_tuple(vertex1, vertex2)];
  edges_.erase(edges_.find(std::make_tuple(vertex1, vertex2)));

  auto v1edge = vertex_edges_.at(vertex1);
  auto v2edge = vertex_edges_.at(vertex2);

  v1edge.erase(std::remove(v1edge.begin(), v1edge.end(), edge));
  v2edge.erase(std::remove(v2edge.begin(), v2edge.end(), edge));

  vertex_edges_[vertex1] = v1edge;
  vertex_edges_[vertex2] = v2edge;
}

// Read MatrixMarket graph file format
bool abm::Graph::read_graph_matrix_market(const std::string& filename) {
  bool status = true;
  try {
    std::fstream file;
    file.open(filename.c_str(), std::ios::in);
    if (file.is_open() && file.good()) {
      // Line
      std::string line;
      bool header = true;
      double ignore;
      while (std::getline(file, line)) {
        std::istringstream istream(line);
        int v1, v2;
        double weight;
        unsigned nvertices;
        // ignore comment lines (# or !) or blank lines
        if ((line.find('#') == std::string::npos) &&
            (line.find('%') == std::string::npos) && (line != "")) {
          if (header) {
            // Ignore header
            istream >> nvertices;
            while (istream.good()) istream >> ignore;
            header = false;
            this->assign_nvertices(nvertices + 1);
          }
          while (istream.good()) {
            // Read vertices edges and weights
            istream >> v1 >> v2 >> weight;
            this->add_edge(v1, v2, weight);
          }
        }
      }
      std::cout << "Graph summary #edges: " << this->edges_.size()
                << " #vertices: " << this->nvertices_ << "\n";
    } else {
      throw std::runtime_error("Input file not found");
    }
  } catch (std::exception& exception) {
    std::cout << "Read matrix market file: " << exception.what() << "\n";
    status = false;
  }
  return status;
}

void abm::Graph::generate_simple_graph() {
  this->assign_nvertices(7);
  // set up a simple graph
  this->add_edge(1, 2, 1.5);
  this->add_edge(1, 3, 9.1);
  this->add_edge(1, 6, 14.3);
  this->add_edge(2, 3, 15.9);
  this->add_edge(2, 4, 5.5);
  this->add_edge(3, 1, 5.6);
  this->add_edge(3, 4, 11.6);
  this->add_edge(3, 6, 2.4);
  this->add_edge(4, 3, 0.2);
  this->add_edge(4, 5, 6.2);
  this->add_edge(5, 6, 9.7);
}

// Dijktra shortest paths from src to all other vertices
std::vector<abm::graph::vertex_t> abm::Graph::dijkstra(
    abm::graph::vertex_t source, abm::graph::vertex_t destination) {

  // Using lambda to compare elements.
  auto compare =
      [](std::pair<abm::graph::weight_t, abm::graph::vertex_t> left,
         std::pair<abm::graph::weight_t, abm::graph::vertex_t> right) {
        return left.first > right.first;
      };

  // Create a priority queue to store weights and vertices
  std::priority_queue<
      std::pair<abm::graph::weight_t, abm::graph::vertex_t>,
      std::vector<std::pair<abm::graph::weight_t, abm::graph::vertex_t>>,
      decltype(compare)>
      priority_queue(compare);

  // Create a vector for distances and initialize all to max
  std::vector<graph::weight_t> distances;
  distances.resize(nvertices_,
                   std::numeric_limits<abm::graph::weight_t>::max());

  // Parent array to store shortest path tree
  std::vector<graph::vertex_t> parent;
  parent.resize(nvertices_, -1);

  // Insert source itself in priority queue & initialize its distance as 0.
  priority_queue.push(std::make_pair(0., source));
  distances.at(source) = 0.;

  // Looping till priority queue becomes empty (or all
  // distances are not finalized)
  while (!priority_queue.empty()) {
    // {min_weight, vertex} sorted based on weights (distance)
    abm::graph::vertex_t u = priority_queue.top().second;
    priority_queue.pop();

    // Break if destination is reached
    if (u == destination) break;

    // Get all adjacent vertices of a vertex
    for (const auto& edge : vertex_edges_[u]) {
      // Get vertex label and weight of neighbours of u.
      const abm::graph::vertex_t neighbour = edge->first.second;
      const abm::graph::weight_t weight = edge->second;

      // Distance from source to neighbour
      // distance_u = distance to current node + weight of edge u to
      // neighbour
      const abm::graph::weight_t distance_u = distances.at(u) + weight;
      // If there is shorted path to neighbour vertex through u.
      if (distances.at(neighbour) > distance_u) {
        parent[neighbour] = u;
        // Update distance of the vertex
        distances.at(neighbour) = distance_u;
        priority_queue.push(std::make_pair(distance_u, neighbour));
      }
    }
  }

  std::vector<abm::graph::vertex_t> path;
  // Iterate until source has been reached
  while (destination != source) {
    destination = parent.at(destination);
    if (destination != source) path.emplace_back(destination);
  }
  // Reverse to arrange from source to destination
  std::reverse(path.begin(), path.end());
  return path;
}

// Dijktra shortest paths from src to all other vertices
abm::ShortestPath abm::Graph::dijkstra_priority_queue(
    abm::graph::vertex_t source, abm::graph::vertex_t destination) {

  // Using lambda to compare elements.
  auto compare =
      [](std::pair<abm::graph::weight_t, abm::graph::vertex_t> left,
         std::pair<abm::graph::weight_t, abm::graph::vertex_t> right) {
        return left.first > right.first;
      };

  // Create a priority queue to store weights and vertices
  std::priority_queue<
      std::pair<abm::graph::weight_t, abm::graph::vertex_t>,
      std::vector<std::pair<abm::graph::weight_t, abm::graph::vertex_t>>,
      decltype(compare)>
      priority_queue(compare);

  // Create a shortest path object.
  ShortestPath sp;
  // Create a vector for distances and initialize all to max
  sp.distances.clear();
  sp.distances.resize(nvertices_,
                      std::numeric_limits<abm::graph::weight_t>::max());

  // Parent array to store shortest path tree
  sp.parent.clear();
  sp.parent.resize(nvertices_, -1);

  // Insert source itself in priority queue and initialize its distance as
  // 0.
  priority_queue.push(std::make_pair(0., source));
  sp.distances.at(source) = 0.;

  // Looping till priority queue becomes empty (or all
  // distances are not finalized)
  while (!priority_queue.empty()) {
    // {min_weight, vertex} sorted based on weights (distance)
    abm::graph::vertex_t u = priority_queue.top().second;
    priority_queue.pop();

    // Break if destination is reached
    if (u == destination) break;

    // Get all adjacent vertices of a vertex
    for (const auto& edge : vertex_edges_[u]) {
      // Get vertex label and weight of neighbours of u.
      const abm::graph::vertex_t neighbour = edge->first.second;
      const abm::graph::weight_t weight = edge->second;

      // Distance from source to neighbour
      // distance_u = distance to current node + weight of edge u to
      // neighbour
      const abm::graph::weight_t distance_u = sp.distances.at(u) + weight;
      // If there is shorted path to neighbour vertex through u.
      if (sp.distances.at(neighbour) > distance_u) {
        sp.parent.at(neighbour) = u;
        // Update distance of the vertex
        sp.distances.at(neighbour) = distance_u;
        priority_queue.push(std::make_pair(distance_u, neighbour));
      }
    }
  }
  // print the path
  if (destination != -1) {
    std::cout << "\nSource: " << source << " destination: " << destination
              << " distance: " << sp.distances.at(destination);
    const auto&& path = sp.get_path(source, destination);
    std::cout << " path: ";
    for (const auto& item : path) std::cout << item << "->";
    std::cout << destination << "\n";
  }
  return sp;
}

// Dijktra shortest paths from src to all other vertices
abm::ShortestPath abm::Graph::dijkstra_priority_queue(
    abm::graph::vertex_t source,
    const std::vector<abm::graph::vertex_t>& destinations) {

  // Using lambda to compare elements.
  auto compare =
      [](std::pair<abm::graph::weight_t, abm::graph::vertex_t> left,
         std::pair<abm::graph::weight_t, abm::graph::vertex_t> right) {
        return left.first > right.first;
      };

  // Create a priority queue to store weights and vertices
  std::priority_queue<
      std::pair<abm::graph::weight_t, abm::graph::vertex_t>,
      std::vector<std::pair<abm::graph::weight_t, abm::graph::vertex_t>>,
      decltype(compare)>
      priority_queue(compare);

  // Check if all destinations have been reached, set all to false
  std::unordered_map<abm::graph::vertex_t, bool> sp_destinations;
  for (const auto& v : destinations) sp_destinations.insert({v, false});

  // Create a shortest path object.
  ShortestPath sp;
  // Create a vector for distances and initialize all to max
  sp.distances.clear();
  sp.distances.resize(nvertices_,
                      std::numeric_limits<abm::graph::weight_t>::max());

  // Parent array to store shortest path tree
  sp.parent.clear();
  sp.parent.resize(nvertices_, -1);

  // Insert source itself in priority queue and initialize its distance as
  // 0.
  priority_queue.push(std::make_pair(0., source));
  sp.distances.at(source) = 0.;

  // Looping till priority queue becomes empty (or all
  // distances are not finalized)
  while (!priority_queue.empty()) {
    // {min_weight, vertex} sorted based on weights (distance)
    abm::graph::vertex_t u = priority_queue.top().second;
    priority_queue.pop();

    // Set when a current destination is reached
    auto itr = sp_destinations.find(u);
    if (itr != sp_destinations.end()) itr->second = true;
    // Break if all destinations are reached
    if (std::all_of(
            sp_destinations.begin(), sp_destinations.end(),
            [&](std::pair<abm::graph::vertex_t, bool> const& destination) {
              return destination.second;
            }))
      break;

    // Get all adjacent vertices of a vertex
    for (const auto& edge : vertex_edges_[u]) {
      // Get vertex label and weight of neighbours of u.
      const abm::graph::vertex_t neighbour = edge->first.second;
      const abm::graph::weight_t weight = edge->second;

      // Distance from source to neighbour
      // distance_u = distance to current node + weight of edge u to
      // neighbour
      const abm::graph::weight_t distance_u = sp.distances.at(u) + weight;
      // If there is shorted path to neighbour vertex through u.
      if (sp.distances.at(neighbour) > distance_u) {
        sp.parent.at(neighbour) = u;
        // Update distance of the vertex
        sp.distances.at(neighbour) = distance_u;
        priority_queue.push(std::make_pair(distance_u, neighbour));
      }
    }
  }
  return sp;
}
