# Created by woochanghwang at 21/04/2020

import sys
sys.path.append('/home/wch23/Project/LifeArc/General')
import pandas as pd
import networkx as nx
import multiprocessing as mp
import toolbox.data_handlers as dh
import itertools


def make_network_based_RWR(graph_G, DIP_genes, DEP_genes, RWR_result_addr):
    # graph_G, DIP_genes, DEP_genes, RWR_result_addr = input_tuple
    print("Start Calculating BC")
    # edge_BC = nx.edge_betweenness_centrality(graph_G)
    edge_BC = nx.edge_betweenness_centrality_subset(graph_G, sources=DIP_genes, targets=DEP_genes)
    start_genes_for_PR = dict()
    for gene in DIP_genes:
        start_genes_for_PR[gene] = 1
    PR_score = nx.pagerank(graph_G, personalization=start_genes_for_PR)
    PR_score_df = pd.DataFrame.from_dict(PR_score, orient='index', columns=['PageRank'])

    PR_score_df.to_csv(RWR_result_addr, sep='\t')
    return PR_score_df

# def run_pagerank(graph_G, DIP_genes, DEP_genes, RWR_result_addr):
#
#     print("Start Calculating BC")
#     # edge_BC = nx.edge_betweenness_centrality(graph_G)
#     edge_BC = nx.edge_betweenness_centrality_subset(graph_G, sources=DIP_genes, targets=DEP_genes)
#     start_genes_for_PR = dict()
#     for gene in DIP_genes:
#         start_genes_for_PR[gene] = 1
#     PR_score = nx.pagerank(graph_G, personalization=start_genes_for_PR)
#     PR_score_df = pd.DataFrame.from_dict(PR_score, orient='index', columns=['PageRank'])
#
#     PR_score_df.to_csv(RWR_result_addr, sep='\t')
#     return PR_score_df

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


    string_addr = "/home/wch23/Project/LifeArc/General/data/STRING/9606.protein.links.v11.0.400.symbol.v2.txt"
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


def get_DIP_to_DEP_paths(RWR_score_df, DEP_list):
    '''
    genes over minimum RWR score of DEP
    :param RWR_score_df:
    :param DEP_list:
    :return:
    '''

    RWR_DEP_df = RWR_score_df[RWR_score_df.index.isin(DEP_list)]
    print(RWR_DEP_df)
    # print(RWR_DEP_df.min(level="PageRank"))
    min_RWR = float(RWR_DEP_df.min())

    DIP_to_DEG_df = RWR_score_df[RWR_score_df['PageRank']>=min_RWR]
    print(DIP_to_DEG_df)

    return RWR_DEP_df
def main_RWR():
    hour = "24"
    DIP_addr = "/home/wch23/Project/LifeArc/COVID-19/data/DIP/covid_krogan_ms_result.txt"
    DEP_addr = "/home/wch23/Project/LifeArc/COVID-19/data/DEP/christian_ms_{}h.p0.05.fc0.5.txt".format(hour)
    string_addr = "/home/wch23/Project/LifeArc/General/data/STRING/9606.protein.links.v11.0.400.symbol.v2.txt"  #Backbone

    with open(string_addr) as network_f:
        string_network_edges = [x.strip().split('\t') for x in network_f.readlines()]

    DIP_df = pd.read_csv(DIP_addr, sep='\t')
    print(DIP_df.columns)

    DEP_df = pd.read_csv(DEP_addr, sep='\t',names=['Gene'])
    DEP_list = DEP_df['Gene'].to_list()
    print(len(DEP_list))

    # COVID19 - sub structures
    covid_sub_structures = list(set(DIP_df['Bait'].to_list()))
    print(covid_sub_structures)
    pool = mp.Pool(mp.cpu_count())

    sub_structure_result = dict()
    for sub_structure in covid_sub_structures:
        DIP_sub_list = DIP_df[DIP_df['Bait'] == sub_structure]['StringGene'].to_list()
        string_G = nx.Graph(string_network_edges)
        # edge_BC = nx.edge_betweenness_centrality_subset(string_G, sources=DIP_sub_list, targets=DEP_list)
        # # for k,v in edge_BC.items():
        # #     print(k,v)
        #
        # edge_BC_no_Zero = dict()
        #
        # for k, v in edge_BC.items():
        #     edge_BC_no_Zero[k] = v + 0.1
        #
        # nx.set_edge_attributes(string_G, edge_BC_no_Zero, 'weight')

        print("SUB : {}".format(sub_structure))
        print("DIP: ", len(DIP_sub_list) )
        RWR_result_addr = "/home/wch23/Project/LifeArc/COVID-19/result/network/RWR_based/RWR_results/COVID_{}_RWR_BC_weighted_result_RWR.tsv".format(
            sub_structure)
        # RWR_result_df = run_pagerank(string_G, DIP_sub_list, DEP_list, RWR_result_addr )
        # DIP_to_DEP_df = get_DIP_to_DEP_paths(RWR_result_df, DEP_list)
        rwr_result = pool.apply(make_network_based_RWR, args=(string_G, DIP_sub_list, DEP_list, sub_structure))

        sub_structure_result[sub_structure] = rwr_result
    pool.close()

    dh.save_obj(sub_structure_result ,"/home/wch23/Project/LifeArc/COVID-19/result/network/RWR_based/RWR_results/COVID_All_Structure_RWR" )
    print(sub_structure_result)

def main_ShortestPath():
    hour = "24"
    data_set = "daniel_rnaseq"
    # data_set = "christian_ms"
    DIP_addr = "/home/wch23/Project/LifeArc/COVID-19/data/DIP/covid_krogan_ms_result.txt"
    DEP_addr = "/home/wch23/Project/LifeArc/COVID-19/data/DEP/{}_{}h.p0.05.fc0.5.txt".format(data_set,hour)




    DIP_df = pd.read_csv(DIP_addr, sep='\t')
    print(hour)
    print(DIP_df.columns)

    DEP_df = pd.read_csv(DEP_addr, sep='\t', names=['Gene'])
    DEP_list = DEP_df['Gene'].to_list()
    print(len(DEP_list))

    # COVID19 - sub structures
    covid_sub_structures = list(set(DIP_df['Bait'].to_list()))
    print(covid_sub_structures)

    sub_structure_result = dict()
    pool = mp.Pool(processes=20)
    for sub_structure in covid_sub_structures:

        DIP_sub_list = DIP_df[DIP_df['Bait'] == sub_structure]['StringGene'].to_list()
        print("SUB : {}".format(sub_structure))
        print("DIP: ", len(DIP_sub_list))
        covid_ppi_deg_pair_list = list(itertools.product(DIP_sub_list, DEP_list))

        # SP_result_addr = "/home/wch23/Project/LifeArc/COVID-19/result/network/SP_based/COVID_{}_all_shortest_paths.tsv".format(
        #     sub_structure)

        SP_result = pool.map(make_shortest_path_sources_to_targets, covid_ppi_deg_pair_list)
        print(len(SP_result))
        sub_structure_result[sub_structure] = SP_result
    pool.close()

    dh.save_obj(sub_structure_result,
                "/home/wch23/Project/LifeArc/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_{}_{}h".format(data_set,hour))
    print(sub_structure_result)

def test():
    test_dict = {'SARS-CoV2 nsp11': {'A': 1}, 'SARS-CoV2 orf8': {'A': 1}, 'SARS-CoV2 nsp2': {'A': 1}, 'SARS-CoV2 nsp4': {'A': 1},
     'SARS-CoV2 N': {'A': 1}, 'SARS-CoV2 nsp13': {'A': 1}, 'SARS-CoV2 E': {'A': 1}, 'SARS-CoV2 Spike': {'A': 1},
     'SARS-CoV2 orf3a': {'A': 1}, 'SARS-CoV2 nsp5': {'A': 1}, 'SARS-CoV2 nsp15': {'A': 1}, 'SARS-CoV2 nsp14': {'A': 1},
     'SARS-CoV2 orf7a': {'A': 1}, 'SARS-CoV2 orf9b': {'A': 1}, 'SARS-CoV2 nsp1': {'A': 1}, 'SARS-CoV2 M': {'A': 1},
     'SARS-CoV2 nsp7': {'A': 1}, 'SARS-CoV2 orf3b': {'A': 1}, 'SARS-CoV2 orf6': {'A': 1}, 'SARS-CoV2': {'A': 1},
     'SARS-CoV2 orf9c': {'A': 1}, 'SARS-CoV2 nsp8': {'A': 1}, 'SARS-CoV2 nsp10': {'A': 1}, 'SARS-CoV2 nsp12': {'A': 1},
     'SARS-CoV2 nsp5_C145A': {'A': 1}, 'SARS-CoV2 orf10': {'A': 1}, 'SARS-CoV2 nsp9': {'A': 1},
     'SARS-CoV2 nsp6': {'A': 1}}

    test_df = pd.DataFrame.from_dict(test_dict)
    print(test_df)
    print(test_df.columns)

if __name__ == '__main__':
    main_ShortestPath()
    # test()