#include "router.h"

// Read OD pairs file format
bool abm::Router::read_od_pairs(const std::string& filename) {
  bool status = true;
  try {
    io::CSVReader<2> in(filename);
    in.read_header(io::ignore_extra_column, "origin", "destination");
    int v1, v2;
    double weight;
    while (in.read_row(v1, v2)) {
      std::array<abm::graph::vertex_t, 2> od = {v1, v2};
      this->od_pairs_.emplace_back(od);
    }
  } catch (std::exception& exception) {
    std::cout << "Read OD file: " << exception.what() << "\n";
    status = false;
  }
  return status;
}
