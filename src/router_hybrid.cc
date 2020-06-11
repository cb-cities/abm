// check all reserves
// update edge weights
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
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
  partial_ods.resize(npagents);
  int sendcounts[nproc], displs[nproc];
  for (int i = 0; i < nproc; i++) {
    sendcounts[i] = npagents;
    displs[i] = npagents*i;
  }
  
  // receive OD from rank 0
  MPI_Datatype mpi_od;
  MPI_Type_vector(4, 1, 1, MPI_LONG_LONG_INT, &mpi_od);
  MPI_Type_commit(&mpi_od);
  MPI_Scatterv((this->input_ods_).data(), sendcounts, displs, mpi_od, partial_ods.data(), npagents, mpi_od, 0, MPI_COMM_WORLD );
  MPI_Type_free( &mpi_od );
  std::cout << "Rank " << myrank << " assigned " << partial_ods.size() << " od pairs out of " << (this->input_ods_).size() << " 100th element " << partial_ods.at(100)[2] << " 2000th element " << partial_ods.at(999)[2] << std::endl;

  // make a map with (hour, quarter) as the key and ODs as the value
  for (auto itr=partial_ods.begin(); itr != partial_ods.end(); ++itr) {
    std::array<abm::graph::vertex_t, 2> od = {(*itr)[2], (*itr)[3]};
    // create first key: hour
    if (this->partial_timed_od_pairs_.find((*itr)[0]) == partial_timed_od_pairs_.end()) {
      std::vector<std::array<abm::graph::vertex_t, 2>> od_values;
      od_values.reserve(npagents);
      this->partial_timed_od_pairs_[(*itr)[0]][(*itr)[1]] = od_values;
    }
    // create second key: quarter
    if (this->partial_timed_od_pairs_[(*itr)[0]].find((*itr)[1]) == partial_timed_od_pairs_[(*itr)[0]].end()) {
      std::vector<std::array<abm::graph::vertex_t, 2>> od_values;
      od_values.reserve(npagents);
      this->partial_timed_od_pairs_[(*itr)[0]][(*itr)[1]] = od_values;
    }
    // append OD to values
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

void abm::Router_hybrid::output_edge_vol_map (const std::string& output_filename) {

  if (output_filename != "no") {
    std::ofstream outputfile;
    outputfile.open (output_filename);
    outputfile << "edgeid,vol\n";
    for (auto itr = edge_vol_.begin(); itr != edge_vol_.end(); ++itr){
      outputfile << itr->first << "," << itr->second << "\n";
    }
    outputfile.close();
  }
}

void abm::Router_hybrid::router (int hour, int quarter, int npagents, int myrank, int nproc) {

  // These variables are used across the whole quarter time step
  int s = omp_get_max_threads(); // maximum OpenMP threads.
  int quarter_od_routed=0, quarter_od_arrival=0; // track total number of ODs whose shortest paths have been calculated.
  std::vector<std::array<abm::graph::vertex_t, 2>> quarter_od = this->partial_timed_od_pairs_[hour][quarter]; // list of ODs to be routed in this time step
  if (myrank==0) {
    std::cout << "hour:" << hour << " quarter:" << quarter << " rank:" << myrank << " routing:" << quarter_od.size() << std::endl;
  }

  // ?? // std::map<int, std::vector<std::array<abm::graph::vertex_t, 4>>> residual_od; // {thread_id: [(o1,d1), (o2,d2), ...]}
  // ?? // std::vector<std::array<abm::graph::vertex_t, 2>> quarter_residual_od_list;
  // for (int t=0; t<s; t++) {
  //   residual_od[t].reserve(npagents);
  // }
  // residual_od_list.reserve(npagents*10);

  // Test shortest path
  // auto itPos = quarter_od.begin() + 1;
  // std::array<abm::graph::vertex_t, 2> prob_od = {401047,299770};
  // quarter_od.insert(itPos, prob_od);
  // auto sp = graph_->dijkstra_edges_with_limit(401047, 299770, 15*60);
  // for (int i=0; i < sp.size(); i++)
  //   std::cout << sp.at(i) << ' ';
  // exit(0);

  // Non-blocking MPI command related
  // MPI_Request sum_request, edge_vol_request;
  // int sum_flag=1, edge_vol_flag=1;
  // MPI_Status  sum_status, edge_vol_status;

  MPI_Datatype mpi_edge_vol;
  MPI_Type_vector(2, 1, 1, MPI_LONG_LONG_INT, &mpi_edge_vol);
  MPI_Type_commit(&mpi_edge_vol);

  while (quarter_od_routed < static_cast<int>(quarter_od.size()))
  // while (quarter_od_routed < 100000)
  {
    if (myrank==0) {
      std::cout << " While loop, rank:" << myrank << " quarter od routed:" << quarter_od_routed << " total:" << quarter_od.size() << std::endl;
    }

    // keep results for this substep per thread
    std::map<int, std::vector<abm::graph::vertex_t>> substep_allthread_edges; // {thread_id: [edge1, edge2, ...]}
    std::map<int, int> substep_allthread_routed, substep_allthread_arrival;
    for (int t=0; t<s; t++) {
      substep_allthread_edges[t].reserve((this->graph_)->nedges()*5);
      substep_allthread_routed[t] = 0;
      substep_allthread_arrival[t] = 0;
    }
    #pragma omp parallel for
    for (unsigned int i=0; i<std::min(1000, (int)(quarter_od.size()-quarter_od_routed)); i++)
    {
      int t = omp_get_thread_num();
      std::array<abm::graph::vertex_t, 2> od = quarter_od[quarter_od_routed+i];
      auto sp = graph_->dijkstra_edges_with_limit(od[0], od[1], 15*60);
      substep_allthread_edges[t].insert(std::end(substep_allthread_edges[t]), std::begin(sp), std::end(sp));
      substep_allthread_routed[t]++;
      abm::graph::vertex_t sp_last_node = graph_->get_edge_ends(sp.back())[1];
      if (sp_last_node != od[1]){
        std::array<abm::graph::vertex_t, 4> one_residual_od = {(int)(hour+quarter+1)/4, (hour+quarter+1)%4, sp_last_node, od[1]};
        // residual_od[t].emplace_back(one_residual_od);
      } else {
        substep_allthread_arrival[t]++;
      }
    } 
    if (myrank==0) {
      std::cout << " After OMP, rank:" << myrank << " routed:" << substep_allthread_arrival[0] << std::endl;
    }

    // Collect results from each thread
    std::map<abm::graph::vertex_t, abm::graph::vertex_t> substep_volume, substep_volume_for_update; //{edge1: edge1_count, edge2: edge2_count, ...}
    for (int t=0; t<s; t++) {
      for (auto & x: substep_allthread_edges[t]) {
        if (substep_volume.find(x) == substep_volume.end())
          substep_volume[x] = 1;
        else
          substep_volume[x]++;      
      }
      substep_allthread_edges[t].clear();
      quarter_od_routed += substep_allthread_routed[t];
      quarter_od_arrival += substep_allthread_arrival[t];
    }
    if (myrank==0) {
      std::cout << " Collect OMP, rank:" << myrank << " so far routed:" << quarter_od_routed << std::endl;
    }

    // Convert map of edge-volume to a vector (easier for MPI to send)
    std::vector<std::array<abm::graph::vertex_t, 2>> substep_volume_vector;
    substep_volume_vector.reserve(npagents*1000);
    for( auto & x : substep_volume ) {
      std::array<abm::graph::vertex_t, 2> one_edge_vol = {x.first, x.second};
      substep_volume_vector.emplace_back(one_edge_vol);
    }

    // Determine total numbers of edge-volume pairs received from each rank
    int recvcounts[nproc], displs[nproc]; // related to gatherv of edge volume from each rank
    int substep_volume_vector_size = substep_volume_vector.size();
    MPI_Allgather(&substep_volume_vector_size, 1, MPI_INT, &recvcounts, 1, MPI_INT, MPI_COMM_WORLD);
    if (myrank==0) {
      std::cout << " Gather length: rank " << myrank << " sending " << substep_volume_vector.size() << " edge-volume pairs" << std::endl;
    }

    // Gather edge-volume pairs from each rank
    std::vector<std::array<abm::graph::vertex_t, 2>> substep_volume_gathered;
    substep_volume_gathered.reserve(npagents*1000);
    int cum_recvcounts = recvcounts[0]; 
    displs[0] = 0;              // offsets into the global array
    for (int i=1; i<nproc; i++) {
      displs[i] = displs[i-1] + recvcounts[i-1];
      cum_recvcounts += recvcounts[i];
    }
    substep_volume_gathered.resize(cum_recvcounts);
    MPI_Allgatherv(substep_volume_vector.data(), substep_volume_vector.size(), mpi_edge_vol, substep_volume_gathered.data(), recvcounts, displs, mpi_edge_vol, MPI_COMM_WORLD);
    if (myrank==0) {
      std::cout << " Gather length: rank " << myrank << " gathering " << substep_volume_gathered.size() << " edge-volume pairs" << std::endl;
    }

    // Combine edge-volume pairs based on edge-id
    for (auto & e_v : substep_volume_gathered) {
      abm::graph::vertex_t edge_id=e_v.at(0);
      if ((this->edge_vol_).find(edge_id) == (this->edge_vol_).end()) {
        (this->edge_vol_)[edge_id] = 1;
      }
      (this->edge_vol_)[edge_id] += 1;
    }
    substep_volume_gathered.clear();
    // Update graph
    for (auto & edge_update : (this->edge_vol_)) {
      abm::graph::vertex_t vol = edge_update.second;
      abm::graph::weight_t fft = graph_->get_edge_fft(edge_update.first);
      abm::graph::weight_t new_weight = fft*(1+pow((vol/2000),2));
      (this->graph_)->update_edge_by_id(edge_update.first, new_weight);
    }
  }

  // Output quarterly edge volume results
  // if (myrank==0) {
  //   output_edge_vol_map("/scratch/07427/bingyu/abm/simulation_outputs/edge_vol/edge_vol_h"+std::to_string(hour)+"_q"+std::to_string(quarter)+".csv");
  // }
  (this->edge_vol_).clear();
  MPI_Type_free( &mpi_edge_vol );
  // MPI_Op_free( &mpi_reduce_edge_vol_op );
  if (myrank==0) {
    std::cout << " Finish hour " << hour << " quarter " << quarter << " on rank " << myrank << std::endl;
  }
};