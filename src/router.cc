#include "router.h"

// Read OD pairs file format
bool abm::Router::read_od_pairs(const std::string& filename, int nagents) {
  bool status = true;
  try {
    io::CSVReader<2> in(filename);
    in.read_header(io::ignore_extra_column, "origin", "destination");
    int v1, v2;
    double weight;
    while (in.read_row(v1, v2)) {
      std::array<abm::graph::vertex_t, 2> od = {v1, v2};
      this->all_od_pairs_.emplace_back(od);
    }
    if (nagents != std::numeric_limits<int>::max())
      all_od_pairs_.resize(nagents);
  } catch (std::exception& exception) {
    std::cout << "Read OD file: " << exception.what() << "\n";
    status = false;
  }
  return status;
}

std::vector<std::array<abm::graph::vertex_t, 2>> abm::Router::compute_routes(
    int mpi_rank) {
  // Create MPI pair type
  MPI_Datatype pair_t;
  MPI_Type_vector(2, 1, 1, MPI_INT, &pair_t);
  MPI_Type_commit(&pair_t);

  // Get number of MPI ranks
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // Calculate chunk size to split router
  int chunk_size = this->all_od_pairs_.size() / mpi_size;
  MPI_Bcast(&chunk_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  std::vector<std::array<abm::graph::vertex_t, 2>> od_pairs(chunk_size);

  // Send route chunks to different compute nodes
  MPI_Scatter(all_od_pairs_.data(), chunk_size, pair_t, od_pairs.data(),
              od_pairs.size(), pair_t, 0, MPI_COMM_WORLD);
  // Calculate the remaining chunk of od_pairs and add to rank 0
  int chunk_remainder = all_od_pairs_.size() % mpi_size;
  if (mpi_rank == 0) {
    od_pairs.insert(od_pairs.begin(), all_od_pairs_.end() - chunk_remainder,
                    all_od_pairs_.end());
  }

  // Paths (vector of edges)
  std::vector<std::array<abm::graph::vertex_t, 2>> paths;
  paths.reserve(graph_->nedges());

  // Indices of start of path and length for each agent
  std::vector<std::array<abm::graph::vertex_t, 3>> paths_idx;
  paths_idx.reserve(od_pairs.size());

#pragma omp parallel for schedule(dynamic)
  for (abm::graph::vertex_t i = 0; i < od_pairs.size(); ++i) {
    const auto sp = graph_->dijkstra(od_pairs[i][0], od_pairs[i][1]);
#pragma omp critical
    {
      paths_idx.emplace_back(std::array<abm::graph::vertex_t, 3>(
          {i, static_cast<abm::graph::vertex_t>(paths.size()),
           static_cast<abm::graph::vertex_t>(sp.size())}));
      paths.insert(std::end(paths), std::begin(sp), std::end(sp));
    }
  }

  if (mpi_rank == 0) {
    // Get all paths and indices
    this->all_paths_ = abm::gather_vector_arrays(paths);
    const auto all_paths_idx = abm::gather_vector_arrays(paths_idx);
  }
  return all_paths_;
}
