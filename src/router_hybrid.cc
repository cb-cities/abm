// check all reserves
// update edge weights
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include "router_hybrid.h"

// Read OD pairs file format
bool abm::Router_hybrid::read_od_pairs(const std::string& filename, int nagents) {
  bool status = true;
  try {
    io::CSVReader<2> in(filename);
    // in.read_header(io::ignore_extra_column, "origin", "destination");
    in.read_header(io::ignore_extra_column, "node_id_igraph_O", "node_id_igraph_D");
    abm::graph::vertex_t v1, v2;
    abm::graph::weight_t weight;
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
  std::cout << "Read " << all_od_pairs_.size() << " OD pairs" << std::endl;
  return status;
}

// Read OD pairs file format
bool abm::Router_hybrid::read_timed_od_pairs(const std::string& filename, int nagents) {
  bool status = true;
  int od_count = 0;
  int random_quarter, tag;
  std::vector<std::array<abm::graph::vertex_t, 3>> timed_od_pairs;
  timed_od_pairs.reserve(nagents);
  srand(time(NULL));
  try {
    io::CSVReader<3> in(filename);
    // in.read_header(io::ignore_extra_column, "origin", "destination");
    in.read_header(io::ignore_extra_column, "node_id_igraph_O", "node_id_igraph_D", "hour");
    abm::graph::vertex_t v1, v2;
    int hour;
    while (in.read_row(v1, v2, hour) && od_count<nagents) {
      if (v1 != v2) {
        random_quarter = rand() % 4;
        tag = hour*4 + random_quarter;
        std::array<abm::graph::vertex_t, 3> timed_od = {tag, v1, v2};
        timed_od_pairs.emplace_back(timed_od);
        od_count++;
      }
    }
  } catch (std::exception& exception) {
    std::cout << "Read OD file: " << exception.what() << "\n";
    status = false;
  }
  std::cout << "Read " << od_count << " OD pairs" << std::endl;

  this->make_timed_od_map(timed_od_pairs, 0, nagents);
  std::cout << "\n" << std::endl;
  return status;
}

void abm::Router_hybrid::make_timed_od_map (
  std::vector<std::array<abm::graph::vertex_t, 3>>& timed_od_pairs, bool print_od_map, int nagents) {
  
  for (auto itr=timed_od_pairs.begin(); itr != timed_od_pairs.end(); ++itr) {
    std::array<abm::graph::vertex_t, 2> od = {(*itr)[1], (*itr)[2]};
    if (this->all_timed_od_pairs_.find((*itr)[0]) == all_timed_od_pairs_.end()) {
      std::vector<std::array<abm::graph::vertex_t, 2>> od_values;
      od_values.reserve(nagents);
      this->all_timed_od_pairs_[(*itr)[0]] = od_values;
    }
    this->all_timed_od_pairs_[(*itr)[0]].emplace_back(od);
  }

  if (print_od_map){
    for (int tag = 0; tag != 30*4; ++tag){
      if (all_timed_od_pairs_[tag].size()>0) {
        std::cout << "H" << (int)tag/4 << " QT" << tag%4 << " OD" << all_timed_od_pairs_[tag].size() << "; ";
      }
    }
  }
}

void abm::Router_hybrid::make_edge_vol_map (
  std::vector<abm::graph::vertex_t>& path) {
  
  for (auto itr=path.begin(); itr != path.end(); ++itr) {
    if (this->edge_vol_.find(*itr) == edge_vol_.end()) {
      this->edge_vol_[(*itr)] = 0;
    }
    this->edge_vol_[(*itr)] ++;
  }
}

void abm::Router_hybrid::output_edge_vol_map (const std::string& output_filename) {

  if (output_filename != "no") {
    std::ofstream outputfile;
    outputfile.open (output_filename);
    outputfile << "edgeid, vol\n";
    for (auto itr = edge_vol_.begin(); itr != edge_vol_.end(); ++itr){
      outputfile << itr->first << ", " << itr->second << "\n";
    }
    outputfile.close();
  }
}

void abm::Router_hybrid::send_od_to_workers (int nproc, int myrank) {
  // Create MPI pair type
  MPI_Datatype triple_t;
  MPI_Type_vector(3, 1, 1, MPI_LONG_LONG_INT, &triple_t);
  MPI_Type_commit(&triple_t);
  int rank=1;
  // Sent od related
  std::array<abm::graph::vertex_t, 3> timed_od;
  // int timed_od_cnt=0;
  // std::vector<int> timed_od_sent;

  for (int tag=0; tag<30*4; ++tag) {
    for (auto itr=all_timed_od_pairs_[tag].begin(); itr<all_timed_od_pairs_[tag].end(); ++itr) {
      timed_od = {static_cast<abm::graph::vertex_t>(tag), (*itr)[0], (*itr)[1]};
      // timed_od_cnt ++;
      // timed_od_sent.push_back(rank);
      MPI_Send(&timed_od, 1, triple_t, rank, tag, MPI_COMM_WORLD);
      // std::cout << "send " << timed_od[0] << " " << timed_od[1] << " " << timed_od[2] << std::endl;
      rank++;
      if (rank==nproc) {rank=1;}
    }
  }
  // std::cout << "send " << timed_od_cnt << " ODs" << std::endl;
  // for (auto i = timed_od_sent.begin(); i != timed_od_sent.end(); ++i) {
  //   std::cout << *i << " ";
  // }
  this->all_timed_od_pairs_.clear();
}

void abm::Router_hybrid::router (int hour, int quarter) {

  // quarter_od = rank_timed_od_pairs_[make_pair(hour,quarter)];
  std::vector<std::array<abm::graph::vertex_t, 2>> quarter_od;
  quarter_od.reserve(all_od_pairs_.size());
  quarter_od = this->all_timed_od_pairs_[12];
  std::array<abm::graph::vertex_t, 2> od;
  std::cout << "router " << quarter_od.size() << std::endl;

  int quarter_od_routed, mm;
  quarter_od_routed = 0;
  int s = omp_get_num_threads();
  std::map<int, int> count;
  // std::map<int, std::vector<abm::graph::vertex_t>> path_collection; // {thread_id: [edge1, edge2, ...]}
  // std::map<int, std::vector<std::array<abm::graph::vertex_t, 2>>> residual_od; // {thread_id: [(o1,d1), (o2,d2), ...]}
  // std::vector<std::array<abm::graph::vertex_t, 2>> residual_od_list;
  // std::map<abm::graph::vertex_t, abm::graph::vertex_t> path_volume; //{edge1: edge1_count, edge2: edge2_count, ...}
  MPI_Request sum_count;
  for (int t=0; t<s; t++) {
    count[t]=0;
  //   path_collection[t].reserve(graph_->nedges()*5);
  //   residual_od[t].reserve(all_od_pairs_.size());
  }
  // residual_od_list.reserve(all_od_pairs_.size());

  // while(quarter_od_routed < quarter_od.size())
  while (quarter_od_routed < std::min(500,static_cast<int>(quarter_od.size()*0.5)))
  {
    #pragma omp parallel for num_threads(6)
    for (unsigned int i=0; i<std::min(1000, static_cast<int>(quarter_od.size()-quarter_od_routed)); i++)
    {
      int t = omp_get_thread_num();
      od = quarter_od[quarter_od_routed+i];
      auto sp = graph_->dijkstra_edges_with_limit(od[0], od[1], 15*60);
      // total_cost += graph_->path_cost(sp);
      // path_collection[t].insert(std::end(path_collection[t]), std::begin(sp), std::end(sp));
      // sp_last_node = graph_->get_edge_ends(sp.back())[1];
      // std::cout << " end at " << sp_last_node << " " << sp.back() << std::endl;
      // if (sp_last_node != od[2]){
      //   one_residual_od = {od[0] + 1, sp_last_node, od[2]};
      //   residual_od[t].emplace_back(one_residual_od);
      // } else {
      //   arrival[t]++;
      // }
      count[t]++;
    } 

    for (int t=0; t<s; t++) {
      // for (auto itr=path_collection[t].begin(); itr != path_collection[t].end(); ++itr) {
      //   if (this->path_volume.find((*itr)) == path_volume.end()) {
      //     this->path_volume[(*itr)] = 1;
      //   }
      //   this->path_volume[(*itr)] += 1;
      // }
      // path_collection[t].clear();

      // for (auto itr=residual_od[t].begin(); itr != residual_od[t].end(); ++itr) {
      //   residual_od_list.emplace_back({(*itr)[0], (*itr)[1]});
      // }
      quarter_od_routed += count[t];
    }

    // Ireduce
    // If Ireduce flag == True: update
    MPI_Iallreduce(&quarter_od_routed, &mm, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, &sum_count);
    std::cout << " all reduce " << mm << " hour " << hour << " quarter " << quarter  << std::endl;
  }
};