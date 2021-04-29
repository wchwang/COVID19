# Created by woochanghwang at 18/11/2020

'''
 I used Personlaised PageRank as RWR in this code.
'''
import sys
import pandas as pd
import networkx as nx
import pickle

def load_obj(file_addr):
    with open(file_addr+ '.pkl', 'rb') as f:
        return pickle.load(f)

def calc_network_centrality_RWR(network, DIP_list, DEP_list, result_save_dir):
    network_nodes = network.nodes()

    try:
        eigenvector_centrality_dict = nx.eigenvector_centrality(network)
    except:
        eigenvector_centrality_dict = dict()
        for node in network_nodes:
            eigenvector_centrality_dict[node] = 0.0
    try:
        degree_centrality_dict = nx.degree_centrality(network)
    except:
        degree_centrality_dict = dict()
        for node in network_nodes:
            degree_centrality_dict[node] = 0.0
    try:
        between_centrality_dict = nx.betweenness_centrality_subset(network, sources=DIP_list, targets=DEP_list)
    except:
        between_centrality_dict = dict()
        for node in network_nodes:
            between_centrality_dict[node] = 0.0

    edge_BC_no_Zero = dict()
    start_genes = ['COVID19']
    edge_BC = nx.edge_betweenness_centrality_subset(network, sources=start_genes, targets=DEP_list)

    for k,v in edge_BC.items():
        edge_BC_no_Zero[k] = v+0.1

    nx.set_edge_attributes(network, edge_BC_no_Zero, 'weight')

    start_genes_for_PR = {'COVID19':1}
    PR_score= nx.pagerank(network, personalization=start_genes_for_PR)

    network_property_df = pd.DataFrame(
        columns=['Eigen', 'Degree', 'Between', 'RWR'])
    for node in network_nodes:
        network_property_df.loc[node] = [eigenvector_centrality_dict[node], degree_centrality_dict[node],
                                         between_centrality_dict[node], PR_score[node]]

    network_property_df.to_csv(result_save_dir)


def main():
    '''
    sys.argv[1] : Hour of DEP (6 or 24)
    :return:
    '''

    hour = sys.argv[1]

    data_type = "{}hr".format(hour)

    print("Read {}hr SIP network".format(hour))
    SIG_G = load_obj("../Result/Network/COVID_All_Structure_All_Shortest_Paths_graph_{}".format(data_type))


    centrality_result_addr= "../Result/Network_analysis/COVID19_centrality_RWR_result_{}.csv".format(data_type)

    DIP_addr = "../Data/DIP/DIP_MS_result.txt"
    DEP_addr = "../Data/DIP/DEP/DEP_MS_{}h.p0.05.fc0.5.txt".format(hour)

    DIP_df = pd.read_csv(DIP_addr, sep='\t')
    DEP_df = pd.read_csv(DEP_addr, sep='\t', names=['Gene'])
    DEP_list = DEP_df['Gene'].to_list()
    DIP_list = DIP_df['StringGene'].to_list()

    print("Start to analyze SIP network using multiple centrality methods")
    calc_network_centrality_RWR(SIG_G, DIP_list, DEP_list, centrality_result_addr)
    print("Finish")


if __name__ == '__main__':
    main()