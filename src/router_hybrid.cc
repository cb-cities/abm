// check all reserves
// update edge weights
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <random>
#include "router_hybrid.h"

// Read OD pairs file format
bool abm::Router_hybrid::read_timed_od_pairs(const std::string& filename, int nagents) {
  bool status = true;
  int od_count = 0;
  int hour, quarter;
  this->input_ods_.reserve(nagents);
  this->input_ods_test_.reserve(nagents);
  srand(time(NULL));
  try {
    io::CSVReader<3> in(filename);
    // in.read_header(io::ignore_extra_column, "origin", "destination");
    in.read_header(io::ignore_extra_column, "node_id_igraph_O", "node_id_igraph_D", "hour");
    abm::graph::vertex_t v1, v2;
    int hour;
    while (in.read_row(v1, v2, hour) && od_count<nagents) {
      if (v1 != v2) {
        quarter = rand() % 4;
        std::array<abm::graph::vertex_t, 4> timed_od = {hour, quarter, v1, v2};
        this->input_ods_.emplace_back(timed_od);
        this->input_ods_test_.emplace_back(od_count);
        od_count++;
      }
    }
  } catch (std::exception& exception) {
    std::cout << "Read OD file: " << exception.what() << "\n";
    status = false;
  }
  std::cout << "Read " << od_count << " OD pairs" << std::endl;
  
  auto rng = std::default_random_engine {};
  std::shuffle(std::begin(this->input_ods_), std::end(this->input_ods_), rng);
  return status;
}

void abm::Router_hybrid::make_timed_od_map (bool print_od_map, int npagents, int nproc, int myrank) {
  
  std::vector<std::array<abm::graph::vertex_t, 4>> partial_ods;
  std::vector<int> partial_ods_test;
  partial_ods_test.resize(npagents);
  partial_ods.resize(npagents);
  int sendcounts[nproc], displs[nproc];
  for (int i = 0; i < nproc; i++) {
    sendcounts[i] = npagents;
    displs[i] = npagents*i;
  }
  // std::cout << "rank " << myrank << " sendcounts " << sendcounts[myrank] << " displs " << displs[myrank] << std::endl;

  // MPI_Scatterv((this->input_ods_test_).data(), sendcounts, displs, MPI_INT, partial_ods_test.data(), 1000, MPI_INT, 0, MPI_COMM_WORLD );
  // partial_ods_test.shrink_to_fit();
  // std::cout << "Rank " << myrank << " assigned " << partial_ods_test.size() << " od pairs out of " << (this->input_ods_test_).size() << " 100th element " << partial_ods_test.at(100) << " 2000th element " << partial_ods_test.at(2000) << std::endl;
  
  MPI_Datatype mpi_od;
  MPI_Type_vector(4, 1, 1, MPI_LONG_LONG_INT, &mpi_od);
  MPI_Type_commit(&mpi_od);
  MPI_Scatterv((this->input_ods_).data(), sendcounts, displs, mpi_od, partial_ods.data(), npagents, mpi_od, 0, MPI_COMM_WORLD );
  MPI_Type_free( &mpi_od );
  std::cout << "Rank " << myrank << " assigned " << partial_ods.size() << " od pairs out of " << (this->input_ods_).size() << " 100th element " << partial_ods.at(100)[2] << " 2000th element " << partial_ods.at(999)[2] << std::endl;

  for (auto itr=partial_ods.begin(); itr != partial_ods.end(); ++itr) {
    std::array<abm::graph::vertex_t, 2> od = {(*itr)[2], (*itr)[3]};
    if (this->partial_timed_od_pairs_.find((*itr)[0]) == partial_timed_od_pairs_.end()) {
      std::vector<std::array<abm::graph::vertex_t, 2>> od_values;
      od_values.reserve(npagents);
      this->partial_timed_od_pairs_[(*itr)[0]][(*itr)[1]] = od_values;
    }
    if (this->partial_timed_od_pairs_[(*itr)[0]].find((*itr)[1]) == partial_timed_od_pairs_[(*itr)[0]].end()) {
      std::vector<std::array<abm::graph::vertex_t, 2>> od_values;
      od_values.reserve(npagents);
      this->partial_timed_od_pairs_[(*itr)[0]][(*itr)[1]] = od_values;
    }
    this->partial_timed_od_pairs_[(*itr)[0]][(*itr)[1]].emplace_back(od);
  }
  std::cout << "Rank " << myrank << " assigned " << partial_ods.size() << " od pairs." << std::endl;

  if (print_od_map){
    for (int hour = 0; hour != 20; ++hour){
      if (partial_timed_od_pairs_[hour].size()>0) {
        for (int quarter = 0; quarter != 4; ++quarter){
          if (partial_timed_od_pairs_[hour][quarter].size()>0) {
            std::cout << "H" << hour << " QT" << quarter << " OD" << partial_timed_od_pairs_[hour][quarter].size() << "; ";
          }
        }
      }
    }
    std::cout << "\n" << std::endl;
  }
}

// void abm::Router_hybrid::make_edge_vol_map (
//   std::vector<abm::graph::vertex_t>& path) {
  
//   for (auto itr=path.begin(); itr != path.end(); ++itr) {
//     if (this->edge_vol_.find(*itr) == edge_vol_.end()) {
//       this->edge_vol_[(*itr)] = 0;
//     }
//     this->edge_vol_[(*itr)] ++;
//   }
// }

// void abm::Router_hybrid::output_edge_vol_map (const std::string& output_filename) {

//   if (output_filename != "no") {
//     std::ofstream outputfile;
//     outputfile.open (output_filename);
//     outputfile << "edgeid, vol\n";
//     for (auto itr = edge_vol_.begin(); itr != edge_vol_.end(); ++itr){
//       outputfile << itr->first << ", " << itr->second << "\n";
//     }
//     outputfile.close();
//   }
// }

// void abm::Router_hybrid::send_od_to_workers (int nproc, int myrank) {
//   // Create MPI pair type
//   MPI_Datatype triple_t;
//   MPI_Type_vector(3, 1, 1, MPI_LONG_LONG_INT, &triple_t);
//   MPI_Type_commit(&triple_t);
//   int rank=1;
//   // Sent od related
//   std::array<abm::graph::vertex_t, 3> timed_od;
//   // int timed_od_cnt=0;
//   // std::vector<int> timed_od_sent;

//   for (int tag=0; tag<30*4; ++tag) {
//     for (auto itr=all_timed_od_pairs_[tag].begin(); itr<all_timed_od_pairs_[tag].end(); ++itr) {
//       timed_od = {static_cast<abm::graph::vertex_t>(tag), (*itr)[0], (*itr)[1]};
//       // timed_od_cnt ++;
//       // timed_od_sent.push_back(rank);
//       MPI_Send(&timed_od, 1, triple_t, rank, tag, MPI_COMM_WORLD);
//       // std::cout << "send " << timed_od[0] << " " << timed_od[1] << " " << timed_od[2] << std::endl;
//       rank++;
//       if (rank==nproc) {rank=1;}
//     }
//   }
//   // std::cout << "send " << timed_od_cnt << " ODs" << std::endl;
//   // for (auto i = timed_od_sent.begin(); i != timed_od_sent.end(); ++i) {
//   //   std::cout << *i << " ";
//   // }
//   this->partial_timed_od_pairs_.clear();
// }

void reduce_edge_vol(
  std::array<abm::graph::vertex_t, 2> *one_path_volume_for_reduce,
  std::vector<std::array<abm::graph::vertex_t, 2>> *path_volume_collected,
  int *len, MPI_Datatype *dtype) {

    (*path_volume_collected).emplace_back(*one_path_volume_for_reduce);

  }

void abm::Router_hybrid::router (int hour, int quarter, int npagents, int myrank, int nproc) {

  std::vector<std::array<abm::graph::vertex_t, 2>> quarter_od;
  if (myrank==0) {
    std::cout << "hour " << hour << " quarter " << quarter << " on rank " << myrank << std::endl;
  }
  quarter_od.reserve(npagents);
  quarter_od = this->partial_timed_od_pairs_[hour][quarter];
  std::array<abm::graph::vertex_t, 2> od;
  if (myrank==0) {
    std::cout << "hour " << hour << " quarter " << quarter << " routing " << quarter_od.size() << " on rank " << myrank << std::endl;
  }

  int quarter_od_routed=0, quarter_od_routed_for_reduce=0, imm=1, mm=0;
  int t, s = omp_get_max_threads();
  abm::graph::vertex_t sp_last_node;
  std::array<abm::graph::vertex_t, 4> one_residual_od;
  std::map<int, int> count, arrival;
  std::map<int, std::vector<abm::graph::vertex_t>> path_collection; // {thread_id: [edge1, edge2, ...]}
  std::map<int, std::vector<std::array<abm::graph::vertex_t, 4>>> residual_od; // {thread_id: [(o1,d1), (o2,d2), ...]}
  std::vector<std::array<abm::graph::vertex_t, 2>> residual_od_list, path_volume_for_reduce, path_volume_collected;
  path_volume_for_reduce.reserve(npagents*1000);
  path_volume_collected.reserve(npagents*1000);
  path_volume_collected.resize(1);
  std::array<abm::graph::vertex_t, 2> one_edge_vol;
  std::map<abm::graph::vertex_t, abm::graph::vertex_t> path_volume, path_volume_for_update; //{edge1: edge1_count, edge2: edge2_count, ...}
  // path_volume.reserve((this->graph_)->nedges());
  // path_volume_for_reduce.reserve((this->graph_)->nedges());
  // path_volume_collected.reserve((this->graph_)->nedges());
  // path_volume_for_reduce.resize(0);
  // path_volume_collected.resize(0);
  int cum_displs, recvcounts[nproc], displs[nproc];

  MPI_Request sum_request, edge_vol_request;
  int sum_flag=1, edge_vol_flag=1;
  MPI_Status  sum_status, edge_vol_status;

  MPI_Datatype mpi_edge_vol;
  MPI_Type_vector(2, 1, 1, MPI_LONG_LONG_INT, &mpi_edge_vol);
  MPI_Type_commit(&mpi_edge_vol);
  // struct edge_vol {abm::graph::vertex_t e; abm::graph::vertex_t v;} e_v;
  // int ev_lengths[2] = {1, 1};
  // MPI_Aint ev_offsets[2] = {offsetof(edge_vol, e), offsetof(edge_vol, v)};
  // MPI_Datatype ev_types[2] = {MPI_LONG_LONG_INT, MPI_LONG_LONG_INT};
  // MPI_Datatype mpi_edge_vol;
  // MPI_Type_struct(2, ev_lengths, ev_offsets, ev_types, &mpi_edge_vol);
  // MPI_Type_commit(&mpi_edge_vol);

  MPI_Op mpi_reduce_edge_vol_op;
  MPI_Op_create((MPI_User_function *)reduce_edge_vol, 1, &mpi_reduce_edge_vol_op);

  for (int t=0; t<s; t++) {
    count[t]=0;
    arrival[t]=0;
    path_collection[t].reserve((this->graph_)->nedges()*5);
    residual_od[t].reserve(npagents);
  }
  residual_od_list.reserve(npagents*10);

  if (myrank==0) {
    std::cout << " Before while loop: rank " << myrank << " quarter od routed " << quarter_od_routed << " limit " << quarter_od.size() << " quarter " << quarter << std::endl;
  }
  while (quarter_od_routed < static_cast<int>(quarter_od.size()))
  // while (quarter_od_routed < 100000)
  {
    if (myrank==0) {
      std::cout << " While loop: rank " << myrank << " quarter od routed " << quarter_od_routed << " limit " << quarter_od.size() << " quarter " << quarter << std::endl;
    }

    #pragma omp parallel for private(t, od, sp_last_node, one_residual_od)
    for (unsigned int i=0; i<std::min(1000, static_cast<int>(quarter_od.size()-quarter_od_routed)); i++)
    // for (unsigned int i=0; i<1000; i++)
    {
      t = omp_get_thread_num();
      od = quarter_od[quarter_od_routed+i];
      auto sp = graph_->dijkstra_edges_with_limit(od[0], od[1], 15*60);
      // total_cost += graph_->path_cost(sp);
      path_collection[t].insert(std::end(path_collection[t]), std::begin(sp), std::end(sp));
      sp_last_node = graph_->get_edge_ends(sp.back())[1];
      // std::cout << " end at " << sp_last_node << " " << sp.back() << std::endl;
      if (sp_last_node != od[1]){
        one_residual_od = {(int)(hour+quarter+1)/4, (hour+quarter+1)%4, sp_last_node, od[1]};
        residual_od[t].emplace_back(one_residual_od);
      } else {
        arrival[t]++;
      }
      count[t]++;
    } 
    if (myrank==0) {
      std::cout << " After OMP: rank " << myrank << " quarter od routed " << quarter_od_routed << " limit " << quarter_od.size() << " quarter " << quarter << std::endl;
    }

    for (int t=0; t<s; t++) {
      for (auto itr=path_collection[t].begin(); itr != path_collection[t].end(); ++itr) {
        if (path_volume.find((*itr)) == path_volume.end()) {
          path_volume[(*itr)] = 1;
        }
        path_volume[(*itr)] += 1;
      }
      path_collection[t].clear();

      // for (auto itr=residual_od[t].begin(); itr != residual_od[t].end(); ++itr) {
      //   residual_od_list.emplace_back({(*itr)[0], (*itr)[1]});
      // }
      quarter_od_routed += count[t];
      count[t] = 0;
    }
    if (myrank==0) {
      std::cout << " Collect OMP: rank " << myrank << " quarter od routed " << quarter_od_routed << " limit " << quarter_od.size() << " quarter " << quarter << std::endl;
    }

    // path_volume_for_reduce.assign(path_volume.begin(), path_volume.end());
    for( auto const& x : path_volume ) {
      one_edge_vol = {x.first, x.second};
      path_volume_for_reduce.emplace_back(one_edge_vol);
    }
    int path_vol_length = path_volume_for_reduce.size();
    MPI_Allgather(&path_vol_length, 1, MPI_INT, &recvcounts, 1, MPI_INT, MPI_COMM_WORLD);
    if (myrank==0) {
      std::cout << " Gather length: rank " << myrank << " quarter od routed " << quarter_od_routed << " limit " << quarter_od.size() << " quarter " << quarter << std::endl;
    }
    displs[0] = 0;              // offsets into the global array
    for (int i=1; i<nproc; i++) {
      displs[i] = displs[i-1] + recvcounts[i-1];
    }
    MPI_Allgatherv(path_volume_for_reduce.data(), path_volume_for_reduce.size(), mpi_edge_vol, path_volume_collected.data(), recvcounts, displs, mpi_edge_vol, MPI_COMM_WORLD);
    if (myrank==0) {
      std::cout << " Gather volume: rank " << myrank << " quarter od routed " << quarter_od_routed << " limit " << quarter_od.size() << " quarter " << quarter << std::endl;
    }

    // // Ireduce
    // if (sum_flag && edge_vol_flag) {
    //   if (myrank==0) {
    //     std::cout << " rank " << myrank << " test true quarter " << quarter << std::endl;
    //   }
    //   quarter_od_routed_for_reduce = quarter_od_routed;
    //   MPI_Iallreduce(&quarter_od_routed_for_reduce, &mm, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, &sum_request);
    //   // path_volume_for_reduce = path_volume;
    //   for( auto const& x : path_volume ) {
    //     one_edge_vol = {x.first, x.second};
    //     path_volume_for_reduce.emplace_back(one_edge_vol);
    //   }
    //   // MPI_Iallreduce(path_volume_for_reduce.data(), path_volume_collected.data(), path_volume_for_reduce.size(), mpi_edge_vol, mpi_reduce_edge_vol_op, MPI_COMM_WORLD, &edge_vol_request);
    //   for (int i = 0; i < nproc; i++) {
    //     if (myrank==i) {
    //       recvcounts[i] = path_volume_for_reduce.size();
    //     }
    //   }
    //   displs[0] = 0;              // offsets into the global array
    //   for (int i=1; i<nproc; i++) {
    //     displs[i] = displs[i-1] + recvcounts[i-1];
    //   }
    //   MPI_Iallgatherv(path_volume_for_reduce.data(), path_volume_for_reduce.size(), mpi_edge_vol, path_volume_collected.data(), recvcounts, displs, mpi_edge_vol, MPI_COMM_WORLD, &edge_vol_request);
    // }
    // MPI_Test(&sum_request, &sum_flag, &sum_status);
    // MPI_Test(&edge_vol_request, &edge_vol_flag, &edge_vol_status);
  }
  // blocking synchronize
  // quarter_od_routed_for_reduce = quarter_od_routed;
  MPI_Allreduce(&quarter_od_routed, &mm, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Type_free( &mpi_edge_vol );
  MPI_Op_free( &mpi_reduce_edge_vol_op );
  if (myrank==0) {
    std::cout << " time step rank " << myrank << " hour " << hour << " quarter " << quarter << " mm  " << mm << std::endl;
  }
};