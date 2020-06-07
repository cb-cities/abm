// check all reserves
// update edge weights
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include "router_mpi.h"

// Read OD pairs file format
bool abm::Router_mpi::read_od_pairs(const std::string& filename, int nagents) {
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
bool abm::Router_mpi::read_timed_od_pairs(const std::string& filename, int nagents) {
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

  this->make_timed_od_map(timed_od_pairs, 1, nagents);
  std::cout << "\n" << std::endl;
  return status;
}

void abm::Router_mpi::make_timed_od_map (
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

void abm::Router_mpi::make_edge_vol_map (
  std::vector<abm::graph::vertex_t>& path) {
  
  for (auto itr=path.begin(); itr != path.end(); ++itr) {
    if (this->edge_vol_.find(*itr) == edge_vol_.end()) {
      this->edge_vol_[(*itr)] = 0;
    }
    this->edge_vol_[(*itr)] ++;
  }
}

void abm::Router_mpi::output_edge_vol_map (const std::string& output_filename) {

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

void abm::Router_mpi::send_od_to_workers (int nproc, int myrank) {
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

void abm::Router_mpi::master (int nproc, int myrank, int nagents) {
  
  // Create MPI pair type
  MPI_Datatype triple_t;
  MPI_Type_vector(3, 1, 1, MPI_LONG_LONG_INT, &triple_t);
  MPI_Type_commit(&triple_t);
  // Status objects
  MPI_Status status;
  // rank iterator
  int rank;
  int tag;
  // Retrieve paths related
  int received_path_cnt=0, received_od_cnt=0;
  int update_path_cnt=0, update_od_cnt=0;
  int received_path_current;
  int receive_path_flag, receive_od_flag;
  int received_edges_size, received_od_size;
  std::vector<abm::graph::vertex_t> received_path_collection;
  std::vector<std::array<abm::graph::vertex_t, 3>> received_od_collection, received_od_single_rank;
  received_path_collection.reserve(10*graph_->nedges());
  received_od_collection.reserve(nagents);
  MPI_Status receive_path_status, receive_od_status;
  // Arrival
  int received_arrival_cnt = 0, received_arrival_single_rank;
  int receive_arrival_flag;
  MPI_Status receive_arrival_status;
  // Output
  int output_file_no = 0;
  
  // Send OD
  send_od_to_workers(nproc, myrank);
  std::cout << "Simulation has " << nagents << " OD " << std::endl;

  // Retrieve paths
  while (1) {
    for (rank = 1; rank < nproc; ++rank) {
      MPI_Iprobe(rank, 2000, MPI_COMM_WORLD, &receive_path_flag, &receive_path_status);
      MPI_Iprobe(rank, 5000, MPI_COMM_WORLD, &receive_od_flag, &receive_od_status);
      MPI_Iprobe(rank, 6000, MPI_COMM_WORLD, &receive_arrival_flag, &receive_arrival_status);
      if (receive_path_flag && receive_od_flag && receive_arrival_flag) {
        // receive residual od
        MPI_Get_count(&receive_od_status, triple_t, &received_od_size);
        received_od_single_rank.resize(received_od_size);
        MPI_Recv(received_od_single_rank.data(), received_od_size, triple_t, rank, 5000, MPI_COMM_WORLD, &status);
        received_od_collection.insert(std::end(received_od_collection), std::begin(received_od_single_rank), std::end(received_od_single_rank));
        received_od_cnt += received_od_size;
        update_od_cnt += received_od_size;
        // receive arrival
        MPI_Recv(&received_arrival_single_rank, 1, MPI_INT, rank, 6000, MPI_COMM_WORLD, &status);
        received_arrival_cnt += received_arrival_single_rank;
        // receive path
        MPI_Get_count(&receive_path_status, MPI_LONG_LONG_INT, &received_edges_size);
        received_path_collection.resize(received_edges_size);
        MPI_Recv(received_path_collection.data(), received_edges_size, MPI_LONG_LONG_INT, rank, 2000, MPI_COMM_WORLD, &status);
        received_path_current = received_path_collection.back();
        received_path_collection.pop_back();
        // reduce to edge volume counts
        make_edge_vol_map(received_path_collection);
        received_path_cnt += received_path_current;
        update_path_cnt += received_path_current;
      }
    }
    // terminate condition
    if (received_arrival_cnt >= nagents*0.9) {
      break;
    }

    // add an update function here. Update when receives certain counts
    if (update_path_cnt > 1000) {
      std::cout << "master receives " << received_path_cnt << " paths, " << received_od_cnt << " residual ods, find " << received_arrival_cnt << " arrival " << std::endl;
    //   reduce_edges;
    //   update_graph; // unblocking send
      make_timed_od_map(received_od_collection, 0, nagents);
      received_od_collection.clear();
      send_od_to_workers(nproc, myrank);
      update_path_cnt = 0;
      update_od_cnt = 0;
    }
  }

  // termination workers
  int stop_value = 1;
  for (rank = 1; rank < nproc; ++rank) {
    MPI_Ssend(&stop_value, 1, MPI_INT, rank, 10000, MPI_COMM_WORLD);
    // std::cout << "sending stop to worker " << rank << std::endl;
  }
  std::cout << "Finishes sending termination requests " << std::endl;

  std::cout << "Before termination, master receives " << received_path_cnt << " paths, " << received_od_cnt << " residual_ods " << received_arrival_cnt << " arrivals" << std::endl;
  for (rank = 1; rank < nproc; ++rank) {
    MPI_Iprobe(rank, 2000, MPI_COMM_WORLD, &receive_path_flag, &receive_path_status);
    MPI_Iprobe(rank, 5000, MPI_COMM_WORLD, &receive_od_flag, &receive_od_status);
    if (receive_path_flag && receive_od_flag) {
      // receive residual od
      MPI_Get_count(&receive_od_status, triple_t, &received_od_size);
      received_od_single_rank.resize(received_od_size);
      MPI_Recv(received_od_single_rank.data(), received_od_size, triple_t, rank, 5000, MPI_COMM_WORLD, &status);
      received_od_cnt += received_od_size;
      update_od_cnt += received_od_size;
      // receive path
      MPI_Get_count(&receive_path_status, MPI_LONG_LONG_INT, &received_edges_size);
      received_path_collection.resize(received_edges_size);
      MPI_Recv(received_path_collection.data(), received_edges_size, MPI_LONG_LONG_INT, rank, 2000, MPI_COMM_WORLD, &status);
      received_path_current = received_path_collection.back();
      received_path_collection.pop_back();
      // reduce to edge volume counts
      make_edge_vol_map(received_path_collection);
      received_path_cnt += received_path_current;
      update_path_cnt += received_path_current;
    }
  }
  std::cout << "After termination, master receives " << received_path_cnt << " paths, " << received_od_cnt << " residual_ods " << received_arrival_cnt << " arrivals" << std::endl;

  output_edge_vol_map("../osm/tokyo_edge_vol.csv");
  return;
};

void abm::Router_mpi::worker (int nproc, int myrank, int nagents) {

  // Create MPI pair type
  MPI_Datatype triple_t;
  MPI_Type_vector(3, 1, 1, MPI_LONG_LONG_INT, &triple_t);
  MPI_Type_commit(&triple_t);
  // Status objects
  MPI_Status status;
  // Related to the numbers of ODs that a worker has received
  int received_od = 0;
  float total_cost = 0;
  int new_od_flag = 0;
  MPI_Request new_od_request;
  // Residual od
  std::vector<std::array<abm::graph::vertex_t, 3>> residual_od, sent_residual_od;
  residual_od.reserve(nagents*5);
  sent_residual_od.reserve(nagents*5);
  int sent_od = 0;
  int sent_od_flag = 1;
  long long int sp_last_node;
  std::array<abm::graph::vertex_t, 3> od, one_residual_od;
  MPI_Request sent_od_request;
  // Related to path collections (and numbers of ODs associated with this collection) sent to worker
  std::vector<abm::graph::vertex_t> worker_path_collection, sent_path_collection;
  worker_path_collection.reserve(graph_->nedges()*5);
  sent_path_collection.reserve(graph_->nedges()*5);
  int sent_path_flag = 1;
  MPI_Request sent_path_request;
  // Related to master issued termination command
  int stop_value = 0;
  int stop_flag = 0;
  // int stop_probe_flag;
  MPI_Request stop_request;
  // arrival counts
  int arrival = 0, sent_arrival;
  int sent_arrival_flag = 1;
  MPI_Request sent_arrival_request;
  
  sleep(1);
  while (1) {

    if (worker_path_collection.size()>50000) {
      std::cout << "rank " << myrank << " path size " << worker_path_collection.size() << std::endl;
    }

    for (int tag = 0; tag<30*4; ++tag) {
      MPI_Irecv(od.data(), 1, triple_t, 0, tag, MPI_COMM_WORLD, &new_od_request); // tag 0 for OD
      // MPI_Wait(&new_od_request, &status);
      MPI_Test(&new_od_request, &new_od_flag, &status);
      if (new_od_flag) {
        // std::cout << od[0] << " " << od[1] << " " << od[2] << std::endl;
        received_od++;
        sent_od++;
        // get path (vector of edges) in 15 min time limit
        auto sp = graph_->dijkstra_edges_with_limit(od[1], od[2], 15*60);
        total_cost += graph_->path_cost(sp);
        worker_path_collection.insert(std::end(worker_path_collection), std::begin(sp), std::end(sp));
        sp_last_node = graph_->get_edge_ends(sp.back())[1];
        // std::cout << " end at " << sp_last_node << " " << sp.back() << std::endl;
        if (sp_last_node != od[2]){
          one_residual_od = {od[0] + 1, sp_last_node, od[2]};
          residual_od.emplace_back(one_residual_od);
        } else {
          arrival ++;
        }
        break;
      }
    }

    if (sent_path_flag && sent_od_flag && sent_arrival_flag && sent_od>0) {
      sent_path_collection.clear();
      sent_path_collection = worker_path_collection;
      sent_residual_od.clear();
      sent_residual_od = residual_od;
      sent_arrival = arrival;
      sent_path_collection.emplace_back(static_cast<abm::graph::vertex_t>(sent_od));
      MPI_Issend(sent_path_collection.data(), sent_path_collection.size(), MPI_LONG_LONG_INT, 0, 2000, MPI_COMM_WORLD, &sent_path_request); 
      MPI_Test(&sent_path_request, &sent_path_flag, &status);
      MPI_Issend(sent_residual_od.data(), sent_residual_od.size(), triple_t, 0, 5000, MPI_COMM_WORLD, &sent_od_request); 
      MPI_Test(&sent_od_request, &sent_od_flag, &status);
      MPI_Issend(&arrival, 1, MPI_INT, 0, 6000, MPI_COMM_WORLD, &sent_arrival_request); 
      MPI_Test(&sent_arrival_request, &sent_arrival_flag, &status);
      
      // std::cout << "worker " << myrank << " posts send " << sent_od << std::endl;
      worker_path_collection.clear();
      residual_od.clear();
      sent_od = 0;
      arrival = 0;
    }
    MPI_Test(&sent_path_request, &sent_path_flag, &status);
    MPI_Test(&sent_od_request, &sent_od_flag, &status);
    MPI_Test(&sent_arrival_request, &sent_arrival_flag, &status);

    // stop if receive master command
    MPI_Irecv(&stop_value, 1, MPI_INT, 0, 10000, MPI_COMM_WORLD, &stop_request); // tag 100 for stop
    MPI_Test(&stop_request, &stop_flag, &status);
    // sleep(1);
    // std::cout << "worker stop value " << stop_value << std::endl;
    if (stop_value) {
      std::cout << "worker " << myrank << " stops by master command, processed " << received_od << " od pairs, average weight " << total_cost/received_od << "  exits" << std::endl;
      return;
    }
  }
};