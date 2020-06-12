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
  std::vector<std::array<abm::graph::vertex_t, 4>> od_inputs;
  od_inputs.reserve(nagents);
  srand(time(NULL));
  try {
    io::CSVReader<3> in(filename);
    in.read_header(io::ignore_extra_column, "node_id_igraph_O", "node_id_igraph_D", "hour");
    abm::graph::vertex_t v1, v2;
    int hour;
    while (in.read_row(v1, v2, hour) && od_count<nagents) {
      if (v1 != v2) {
        int quarter = rand() % 4;
        std::array<abm::graph::vertex_t, 4> timed_od = {hour, quarter, v1, v2};
        od_inputs.emplace_back(timed_od);
        od_count++;
      }
    }
  } catch (std::exception& exception) {
    std::cout << "Read OD file: " << exception.what() << "\n";
    status = false;
  }
  std::cout << "Rank 0 reads " << od_count << " OD pairs" << std::endl;
  
  // give it as an instance variable to the router object at rank 0
  auto rng = std::default_random_engine {};
  std::shuffle(std::begin(od_inputs), std::end(od_inputs), rng);
  this->input_ods_ = make_timed_od_map(0, nagents, od_inputs);
  
  return status;
}

std::map<int, std::map<int, std::vector<std::array<abm::graph::vertex_t, 2>>>> abm::Router_hybrid::make_timed_od_map (
  bool print_od_map, int nagents, std::vector<std::array<abm::graph::vertex_t, 4>>& od_inputs) {

  std::map<int, std::map<int, std::vector<std::array<abm::graph::vertex_t, 2>>>> timed_od_map;
  // make a map with (hour, quarter) as the key and ODs as the value
  for (auto & x: od_inputs) {
    // create first key: hour
    if (timed_od_map.find(x[0]) == timed_od_map.end()) {
      std::vector<std::array<abm::graph::vertex_t, 2>> od_values;
      od_values.reserve(nagents);
      timed_od_map[x[0]][x[1]] = od_values;
    }
    // create second key: quarter
    if (timed_od_map[x[0]].find(x[1]) == timed_od_map[x[0]].end()) {
      std::vector<std::array<abm::graph::vertex_t, 2>> od_values;
      od_values.reserve(nagents);
      timed_od_map[x[0]][x[1]] = od_values;
    }
    // append OD to values
    timed_od_map[x[0]][x[1]].push_back({x[2], x[3]});
  }

  if (print_od_map)
  {
    for (int hour = 0; hour != 20; ++hour)
      if (timed_od_map[hour].size()>0) 
        for (int quarter = 0; quarter != 4; ++quarter)
          if (timed_od_map[hour][quarter].size()>0) 
            std::cout << "H" << hour << " QT" << quarter << " OD" << timed_od_map[hour][quarter].size() << "; ";
    std::cout << "\n" << std::endl;
  }

  return timed_od_map;
}

void abm::Router_hybrid::quarter_router (int hour, int quarter, int subp_agents, int myrank, int nproc) {

  std::vector<std::array<abm::graph::vertex_t, 2>> quarter_ods; // quarter_ods are from rank 0 and are to be scattered to each rank.
  if (myrank==0)
    quarter_ods = (this->input_ods_)[hour][quarter];

  int quarter_od_routed = 0, quarter_od_total;
  MPI_Bcast(&quarter_od_total, 1, MPI_INT, 0, MPI_COMM_WORLD);

  while (quarter_od_routed < quarter_od_total) {
    
    // each substep needs to process among all ranks
    int substep_od_size = std::min(subp_agents, quarter_od_total-quarter_od_routed);
    if ((quarter_od_total-quarter_od_routed) < (quarter_od_routed + nproc*10)) {
      substep_od_size = quarter_od_total-quarter_od_routed;
    }
    quarter_od_routed += substep_od_size;
    // scatter to each rank
    int substep_od_size_per_rank = (int)(substep_od_size/nproc);
    int sendcounts[nproc], senddispls[nproc];
    for (int i = 0; i < nproc; i++) {
      sendcounts[i] = substep_od_size_per_rank;
      senddispls[i] = substep_od_size_per_rank*i;
    }
    // if cannot be distributed evenly among ranks
    sendcounts[nproc] = substep_od_size - sendcounts[nproc-1];
    
    // each rank receives 
    std::vector<std::array<abm::graph::vertex_t, 2>> partial_ods;
    partial_ods.resize(sendcounts[myrank]);
    MPI_Datatype mpi_od;
    MPI_Type_vector(2, 1, 1, MPI_LONG_LONG_INT, &mpi_od);
    MPI_Type_commit(&mpi_od);
    MPI_Scatterv(quarter_ods.data(), sendcounts, senddispls, mpi_od, partial_ods.data(), sendcounts[nproc], mpi_od, 0, MPI_COMM_WORLD);
    std::cout << "Rank " << myrank << " assigned " << partial_ods.size() << " od pairs out of " << quarter_od_total << " 1st element " << partial_ods[0][0] << std::endl;

    // route
    std::vector<std::array<abm::graph::vertex_t, 2>> substep_volume_vector = substep_router (myrank, partial_ods);

    // Determine total numbers of edge-volume pairs received from each rank
    int recvcounts[nproc], recvdispls[nproc]; // related to gatherv of edge volume from each rank
    int substep_volume_vector_size = substep_volume_vector.size();
    MPI_Allgather(&substep_volume_vector_size, 1, MPI_INT, &recvcounts, 1, MPI_INT, MPI_COMM_WORLD);
    if (myrank==0) {
      std::cout << " Gather length: rank " << myrank << " sending " << substep_volume_vector.size() << " edge-volume pairs" << std::endl;
    }

    // Gather edge-volume pairs from each rank
    std::vector<std::array<abm::graph::vertex_t, 2>> substep_volume_gathered;
    int cum_recvcounts = recvcounts[0]; 
    recvdispls[0] = 0;              // offsets into the global array
    for (int i=1; i<nproc; i++) {
      recvdispls[i] = recvdispls[i-1] + recvcounts[i-1];
      cum_recvcounts += recvcounts[i];
    }
    substep_volume_gathered.resize(cum_recvcounts);
    MPI_Allgatherv(substep_volume_vector.data(), substep_volume_vector.size(), mpi_od, substep_volume_gathered.data(), recvcounts, recvdispls, mpi_od, MPI_COMM_WORLD);
    if (myrank==0) {
      std::cout << " Gather length: rank " << myrank << " gathering " << substep_volume_gathered.size() << " edge-volume pairs" << std::endl;
    }
    MPI_Type_free( &mpi_od );

    // Combine edge-volume pairs based on edge-id
    for (auto & e_v : substep_volume_gathered) {
      abm::graph::vertex_t edge_id=e_v[0];
      if ((this->edge_vol_).find(edge_id) == (this->edge_vol_).end())
        (this->edge_vol_)[edge_id] = 1;
      else
        (this->edge_vol_)[edge_id]++;
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

  // clear results for next quarter
  (this->edge_vol_).clear();
  if (myrank==0) {
    std::cout << " Finish hour " << hour << " quarter " << quarter << " on rank " << myrank << std::endl;
  }
}

std::vector<std::array<abm::graph::vertex_t, 2>> abm::Router_hybrid::substep_router (
  int myrank, std::vector<std::array<abm::graph::vertex_t, 2>>& partial_ods) {
  
  // OpenMP calculate shortest path
  std::map<abm::graph::vertex_t, abm::graph::vertex_t> substep_volume_map; // {edge_id: vol}
  int substep_arrival;

  #pragma omp parallel
  {
    std::vector<abm::graph::vertex_t> substep_onethread_edges;
    substep_onethread_edges.reserve((this->graph_)->nedges()*5);
    int substep_onethread_arrival=0;
    #pragma omp for
    for (unsigned int i=0; i<partial_ods.size(); i++) {
      int t = omp_get_thread_num();
      std::array<abm::graph::vertex_t, 2> od = partial_ods[i];
      auto sp = graph_->dijkstra_edges_with_limit(od[0], od[1], 15*60);
      substep_onethread_edges.insert(std::end(substep_onethread_edges), std::begin(sp), std::end(sp));
      // residual demand
      abm::graph::vertex_t sp_last_node = graph_->get_edge_ends(sp.back())[1];
      if (sp_last_node != od[1]){
        std::array<abm::graph::vertex_t, 2> one_residual_od = {sp_last_node, od[1]};
        // residual_od[t].emplace_back(one_residual_od);
      } else {
        substep_onethread_arrival++;
      }
    }
    // join results from each thread
    #pragma omp critical
    {
      substep_arrival += substep_onethread_arrival;
      for (auto & x: substep_onethread_edges) {
        if (substep_volume_map.find(x) == substep_volume_map.end())
          substep_volume_map[x] = 1;
        else
          substep_volume_map[x]++;      
      }
    }
  }
  if (myrank==0) {
    std::cout << " After OMP, rank:" << myrank << std::endl;
  }

  // Convert the std::map<edge, volume> to a vector (easier to send using MPI)
  std::vector<std::array<abm::graph::vertex_t, 2>> substep_volume_vector; // [(edge_id: vol), ...] for gathering
  substep_volume_vector.reserve(partial_ods.size()*1000);
  for( auto & x : substep_volume_map ) {
    substep_volume_vector.push_back({x.first, x.second});
  }

  return substep_volume_vector;
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

// void abm::Router_hybrid::router (int hour, int quarter, int npagents, int myrank, int nproc) {
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
// };