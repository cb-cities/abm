import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse as ssparse
from ctypes import *
import scipy.io as sio
import subprocess

sys.path.insert(0, '/home/bingyu')
from sp import interface

def process_observations(links_df0):

    ### read validation data: traffic flow by quarter
    quarterly_measures = pd.read_csv('quarterly_measure.csv')

    ### get observation groups
    group_measures = quarterly_measures.groupby('obs_grp_id').agg({'start': 'first', 'end': 'first', 'dir': 'first'}).reset_index()

    ### this network prioritize motorway and motorway link
    g = interface.readgraph(bytes('network_sparse_mex_analysis.mtx', encoding='utf-8'))

    ### find the node id in graph
    nodes_df = pd.read_csv('../osm/tokyo_nodes.csv')
    group_measures = group_measures.merge(nodes_df[['node_osmid', 'node_id_igraph']], how='left', left_on='start', right_on='node_osmid')
    group_measures = group_measures.merge(nodes_df[['node_osmid', 'node_id_igraph']], how='left', left_on='end', right_on='node_osmid', suffixes=['_start', '_end'])
    group_measures = group_measures.dropna(subset=['node_id_igraph_start', 'node_id_igraph_end'])
    group_measures['start_igraph'] = group_measures['node_id_igraph_start'].astype(int)
    group_measures['end_igraph'] = group_measures['node_id_igraph_end'].astype(int)

    ### build edge-obs_grp dataframe
    obs_grp_edge_list = []
    for row in group_measures.itertuples():
        obs_grp_id = getattr(row, 'obs_grp_id')
        dir = getattr(row, 'dir')
        if obs_grp_id in [110, 116]:
            continue
        elif obs_grp_id in range(43,88):
            if dir == 2:
                start_igraph = getattr(row, 'end_igraph')
                end_igraph = getattr(row, 'start_igraph')
            elif dir == 1:
                start_igraph = getattr(row, 'start_igraph')
                end_igraph = getattr(row, 'end_igraph')
            else:
                print('invalid direction')
                print(dir)
        else:
            if dir == 1:
                start_igraph = getattr(row, 'end_igraph')
                end_igraph = getattr(row, 'start_igraph')
            elif dir == 2:
                start_igraph = getattr(row, 'start_igraph')
                end_igraph = getattr(row, 'end_igraph')
            else:
                print('invalid direction')
                print(dir)

        try:
            sp = g.dijkstra(start_igraph+1, end_igraph+1)
        except ArgumentError:
            print(end_igraph)
        sp_dist = sp.distance(end_igraph+1)
        if sp_dist > 10e7:
            print('route not found')
            sp.clear()
        else:
            sp_route = sp.route(end_igraph+1)
            route_igraph = [(start_sp-1, end_sp-1) for (start_sp, end_sp) in sp_route]
            if len(route_igraph)>20:
                pass
            else:
                obs_grp_edge_list += [(start, end, obs_grp_id) for (start, end) in route_igraph]
            sp.clear()

    obs_grp_edge_df = pd.DataFrame(obs_grp_edge_list, columns=['start_igraph', 'end_igraph', 'obs_grp_id'])
    obs_grp_edge_df = pd.merge(obs_grp_edge_df, links_df0[['start_igraph', 'end_igraph', 'edge_id_igraph', 'length']], how='left', on=['start_igraph', 'end_igraph'])
    return quarterly_measures, obs_grp_edge_df

def main():

    ### read edges_df file
    links_df0 = pd.read_csv('../osm/tokyo_edges.csv')
    links_df0['fft'] = links_df0['length']/links_df0['maxmph']*2.237

    ### read observation file
    quarterly_measures, obs_grp_edge_df = process_observations(links_df0)

    ### daily measures
    # daily_measures = quarterly_measures.groupby('obs_grp_id').agg({'Q': np.sum}).reset_index()
    daily_measures = quarterly_measures[quarterly_measures['start_quarter']==12].reset_index()

    # for highway_discount in [0.001, 0.01] + list(range(0.1, 2.2, 0.3)):
    for highway_discount in [2]:

        ### adjust edge weights
        links_df0['fft_adj'] = np.where(links_df0['type'].isin(['motorway', 'motorway_link']), links_df0['fft']*highway_discount, links_df0['fft'])
        links_df0[["edge_id_igraph", "start_igraph", "end_igraph", "fft_adj"]].to_csv('../osm/tokyo_edges_discount.csv', index=False)

        ### run simulation
        # subprocess.run(["export", "OMP_NUM_THREADS=10"])
        os.environ["OMP_NUM_THREADS"] = "10"
        subprocess.run(["mpirun", "-n", "1", "../build/abm", "{}".format(highway_discount)])

        ### process simulationr results
        daily_sim_vol = obs_grp_edge_df.copy()
        daily_sim_vol = daily_sim_vol.drop_duplicates(subset='edge_id_igraph')
        daily_sim_vol['daily_vol'] = 0
        for hour in range(3,4):
            for quarter in range(1):
                quarter_sim_vol = pd.read_csv('../simulation_outputs/edge_vol_discount/edge_vol_h{}_q{}_adj{}.csv'.format(hour, quarter, highway_discount))
                quarter_sim_vol = quarter_sim_vol.loc[quarter_sim_vol['edgeid'].isin(daily_sim_vol['edge_id_igraph'])]
                daily_sim_vol = pd.merge(daily_sim_vol, quarter_sim_vol, how='left', left_on='edge_id_igraph', right_on='edgeid')
                daily_sim_vol['daily_vol'] += daily_sim_vol.fillna(value={'vol': 0})['vol']
                daily_sim_vol = daily_sim_vol[['obs_grp_id', 'edge_id_igraph', 'length', 'daily_vol']]
                # break
            # break
        # daily_sim_vol.sort_values(by='daily_vol', ascending=False).head()

        ### compare
        daily_sim_vol['Q_sim'] = daily_sim_vol['daily_vol']*24*4 ### tmp
        compare_df = daily_measures.merge(daily_sim_vol[['obs_grp_id', 'Q_sim']], how='left', on='obs_grp_id')
        compare_df['Q_sim'] = compare_df['Q_sim'].fillna(0)
        compare_df['diff'] = compare_df['Q_sim'] - compare_df['Q']
        print(compare_df[['Q', 'Q_sim', 'diff']].describe())

if __name__ == "__main__":
    main()