#include "router.h"

// Read OD pairs file format
bool abm::Router::read_od_pairs(const std::string& filename) {
  bool status = true;
  try {
    std::fstream file;
    file.open(filename.c_str(), std::ios::in);
    if (file.is_open() && file.good()) {
      // Line
      std::string line;
      double ignore;
      while (std::getline(file, line)) {
        std::istringstream istream(line);
        int v1, v2;
        double weight;
        // ignore comment lines (# or !) or blank lines
        if ((line.find('#') == std::string::npos) &&
            (line.find('%') == std::string::npos) && (line != "")) {
          while (istream.good()) {
            // Read vertices edges and weights
            istream >> ignore >> v1 >> v2 >> weight;
            this->od_pairs_.emplace_back(std::make_pair(v1, v2));
          }
        }
      }
    } else {
      throw std::runtime_error("Input file not found");
    }
  } catch (std::exception& exception) {
    std::cout << "Read OD file: " << exception.what() << "\n";
    status = false;
  }
  return status;
}
