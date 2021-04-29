# Created by woochanghwang at 18/11/2020

import pandas as pd
import multiprocessing as mp
import itertools
import networkx as nx
import sys
import pickle


def save_obj(obj, file_addr ):
    with open(file_addr + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def make_shortest_path_sources_to_targets(pair):
    '''
    To make network between covid ppi to covid DEG genes
    All pairs are connected on STRING
    :return:
    '''


    string_addr = "../Data/STRING/9606.protein.links.v11.0.400.symbol.v2.txt"
    with open(string_addr) as network_f:
        string_network_edges = [x.strip().split('\t') for x in network_f.readlines()]

    string_G = nx.Graph(string_network_edges)

    all_path_list = []

    try:
        for p in nx.all_shortest_paths(string_G, source=pair[0], target=pair[1]):
            pairs_from_p = pairwise(p)
            all_path_list += pairs_from_p
    except:
        print(pair, "No Path")


    all_path_G = nx.Graph(all_path_list)
    all_edges = all_path_G.edges()

    all_edges_w = ['\t'.join(edge) for edge in all_edges]

    # with open("/home/wch23/Project/LifeArc/COVID-19/result/network/ppi_to_deg_paths.tsv",'w') as ppi_to_deg_f:
    #     ppi_to_deg_f.write('\n'.join(all_edges_w))

    return all_edges

def main():
    '''
    sys.argv[1] : Hour of DEP (6 or 24)
    sys.argv[2] : The number of proccess for multiprocess
    :return:
    '''
    if len(sys.argv) < 1:
        sys.exit(0)

    hour = sys.argv[1] #    "6"  # It is for 6hr SIP network
    # hour = "24" # It is for 24hr SIP network

    DIP_addr = "../Data/DIP/DIP_MS_result.txt"
    DEP_addr = "../Data/DIP/DEP/DEP_MS_{}h.p0.05.fc0.5.txt".format(hour)

    DIP_df = pd.read_csv(DIP_addr, sep='\t')
    print("Create {}hour SIP network".format(hour))

    DEP_df = pd.read_csv(DEP_addr, sep='\t', names=['Gene'])
    DEP_list = DEP_df['Gene'].to_list()

    # COVID19 - sub structures
    covid_sub_structures = list(set(DIP_df['Bait'].to_list()))

    sub_structure_result = dict()
    # pool = mp.Pool(processes=20)
    pool = mp.Pool(processes=sys.argv[2])
    for sub_structure in covid_sub_structures:
        DIP_sub_list = DIP_df[DIP_df['Bait'] == sub_structure]['StringGene'].to_list()
        print("SUB : {}".format(sub_structure))
        print("DIP: ", len(DIP_sub_list))
        covid_ppi_deg_pair_list = list(itertools.product(DIP_sub_list, DEP_list))

        # SP_result_addr = "/home/wch23/Project/LifeArc/COVID-19/result/network/SP_based/COVID_{}_all_shortest_paths.tsv".format(
        #     sub_structure)

        SP_result = pool.map(make_shortest_path_sources_to_targets, covid_ppi_deg_pair_list)

        sub_structure_result[sub_structure] = SP_result
    pool.close()

    SIG_G = nx.Graph(sub_structure_result)

    save_obj(SIG_G,
                "../Result/Network/COVID_All_Structure_All_Shortest_Paths_graph_{}h".format(hour))
    print("Finish to create SIP network")


if __name__ == '__main__':
    main()