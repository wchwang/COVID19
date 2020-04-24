# Created by woochanghwang at 22/04/2020

import sys
sys.path.insert(0,'/home/wch23/Project/LifeArc/General')
from toolbox.network_utilities import randomize_graph, calc_RWR_from_Graph
import networkx as nx
import random
import time
import multiprocessing as mp
import toolbox.data_handlers as dh


def centrality_dict_to_list(centriality_result_dict):

    centrality_result_list = []
    for key,value in centriality_result_dict.items():
        centrality_result_list.append([key, value])

    return centrality_result_list

def analysis_centrality_RWR(input):


    original_G, process_num = input
    # eigen_permutation_result = []
    # bc_permutation_result = []
    # dc_permutation_result = []
    # rwr_permutation_result = []

    start_time = time.time()

    # dep_num = 164  # 24hr
    dep_num = 64  # 6hr
    # dep_num = 117   # RNA_24hr
    print("DEP num", dep_num)
    random_G = randomize_graph(original_G, randomization_type="preserve_degree_distribution")
    # random_G = randomize_graph(original_G, randomization_type="random")
    print("Random node, edges", len(random_G.nodes()), len(random_G.edges()))

    random_start_genes = random.choices(list(random_G.nodes()), k=334)
    random_target_genes = random.choices(list(random_G.nodes()), k=dep_num)
    between_centrality_dict = nx.betweenness_centrality_subset(random_G, sources=random_start_genes,
                                                               targets=random_target_genes)

    degree_centrality_dict = nx.degree_centrality(random_G)
    eigenvector_centrality_dict = nx.eigenvector_centrality(random_G)

    # eigen_permutation_result += centrality_dict_to_list(eigenvector_centrality_dict)
    # dc_permutation_result += centrality_dict_to_list(degree_centrality_dict)
    # bc_permutation_result += centrality_dict_to_list(between_centrality_dict)

    random_start_genes = random.choices(list(random_G.nodes()), k=1)

    random_edge_BC = nx.edge_betweenness_centrality_subset(random_G, sources=random_start_genes,
                                                           targets=random_target_genes)

    random_edge_BC_noZero = dict()

    for k, v in random_edge_BC.items():
        random_edge_BC_noZero[k] = v + 0.1
    nx.set_edge_attributes(random_G, random_edge_BC_noZero, 'weight')

    start_genes_for_PR = {random_start_genes[0] : 1}
    rwr_dict = nx.pagerank(random_G, personalization=start_genes_for_PR)
    # print(rwr_dict)
    # rwr_permutation_result += centrality_dict_to_list(rwr_dict)


    print("{} times training Runtime: {:2f} Minutes".format(process_num, ((time.time() - start_time) / 60)))

    return {'eigen' : eigenvector_centrality_dict, 'between':between_centrality_dict, 'degee': degree_centrality_dict, 'rwr': rwr_dict}
    # return eigen_permutation_result

def main_COVID19():
    '''
    create 2020-04-01
    :return:
    '''
    # round_num = 'round1'
    # hour = "RNA_24hr"
    # data_type = "daniel_rnaseq_24h"

    round_num = "round2"
    hour = "6hr"
    data_type = "6hr"
    original_G = dh.load_obj(
        "/mnt/raid0_data/wch23/Project/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_graph_{}".format(
            data_type))
    # original_G = nx.gnp_random_graph(n=50, p=0.5)

    permutation_num = list(range(1,501))
    graph_list = [original_G] * len(permutation_num)
    print(permutation_num)

    pool = mp.Pool(processes=20)

    results = pool.map(analysis_centrality_RWR, zip(graph_list, permutation_num))

    pool.close()

    dh.save_obj(results,"/mnt/raid0_data/wch23/Project/COVID-19/result/Permutation_NA/{}/{}/COVID_network_analysis_permutation_{}".format(hour, round_num, data_type) )



if __name__ == '__main__':
    main_COVID19()
    # test_main()