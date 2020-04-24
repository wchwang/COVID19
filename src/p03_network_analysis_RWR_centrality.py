# Created by woochanghwang at 22/04/2020
'''
 I used Personlaised PageRank as RWR in this code.

'''
import sys
sys.path.append('/home/wch23/Project/LifeArc/General')
import pandas as pd
import networkx as nx
import toolbox.data_handlers as dh

def calc_network_centrality_RWR(network, DIP_list, DEP_list, result_save_dir):
    network_nodes = network.nodes()

    # print(degseq_result_df.head())
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
        # between_centrality_dict = nx.betweenness_centrality(network)
        # covid_sig_gene_addr = "/home/wch23/Project/LifeArc/COVID-19/data/covid_ppi_ms_24h.p0.05.fc0.5.txt"
        # covid_sig_gene_df = pd.read_csv(covid_sig_gene_addr, sep='\t', names=['Type', 'Gene'])
        # deg_genes = covid_sig_gene_df[covid_sig_gene_df['Type'] == 'MS']['Gene'].to_list()
        # ppi_genes = covid_sig_gene_df[covid_sig_gene_df['Type'] == 'PPI']['Gene'].to_list()
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
    # PR_score_df = pd.DataFrame.from_dict(PR_score, orient='index', columns=['PageRank'])
    # RWR_result_addr = "/home/wch23/Project/LifeArc/COVID-19/result/RWR/test_PR/COVID_{}_RWR_BC_weighted_result_Pagerank_t1.tsv".format(
    #     '_'.join(start_genes))
    # PR_score_df.to_csv(RWR_result_addr, sep='\t')


    network_property_df = pd.DataFrame(
        columns=['Eigen', 'Degree', 'Between', 'RWR'])
    for node in network_nodes:
        network_property_df.loc[node] = [eigenvector_centrality_dict[node], degree_centrality_dict[node],
                                         between_centrality_dict[node], PR_score[node]]

    network_property_df.to_csv(result_save_dir)
    print(network_property_df.head())



def main_COVID19():

    # sig_interaction_addr = "/home/wch23/Project/LifeArc/COVID-19/result/network/COVID_core_interactions_both_with_covid.tsv"  # Network1
    # sig_interaction_addr = "/home/wch23/Project/LifeArc/COVID-19/result/network/COVID_core_interactions_shortest_both_with_PPI_DEG.tsv" #Network2
    # sig_interaction_addr = "/home/wch23/Project/LifeArc/COVID-19/result/network/COVID_shortest_both_with_PPI_DEG.tsv"  # Network3
    # sig_interaction_addr = "/home/wch23/Project/LifeArc/COVID-19/result/network/COVID_core_interactions_shortest_both_with_PPI_DEG_n4.tsv" #Network4
    # sig_interaction_addr = "/home/wch23/Project/LifeArc/COVID-19/result/network/COVID_core_interactions_shortest_both_with_PPI_DEG_n5.tsv"  # 24hr Network5
    # sig_interaction_addr = "/home/wch23/Project/LifeArc/COVID-19/result/network/COVID_core_interactions_shortest_both_with_PPI_DEG_6hr_n1.tsv" # 6hr Network1
    round_num = 'round6'
    # hour = "RNA_24"
    hour = "24"
    # data_type = "6hr"
    data_type = "24hr"
    # data_type = "daniel_rnaseq_24h"

    covid19_G = dh.load_obj("/home/wch23/Project/LifeArc/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_graph_{}".format(data_type))


    centrality_result_addr= "/home/wch23/Project/LifeArc/COVID-19/result/Centrality/{}hr/{}/COVID19_centrality_RWR_result_{}_{}.csv".format(
        hour, round_num, data_type, round_num
    )

    DIP_addr = "/home/wch23/Project/LifeArc/COVID-19/data/DIP/covid_krogan_ms_result.txt"
    DEP_addr = "/home/wch23/Project/LifeArc/COVID-19/data/DEP/christian_ms_{}h.p0.05.fc0.5.txt".format(hour)
    # DEP_addr = "/home/wch23/Project/LifeArc/COVID-19/data/DEP/daniel_rnaseq_24h.p0.05.fc0.5.txt"

    DIP_df = pd.read_csv(DIP_addr, sep='\t')
    DEP_df = pd.read_csv(DEP_addr, sep='\t', names=['Gene'])
    DEP_list = DEP_df['Gene'].to_list()
    DIP_list = DIP_df['StringGene'].to_list()
    calc_network_centrality_RWR(covid19_G, DIP_list, DEP_list, centrality_result_addr)



if __name__ == '__main__':
    main_COVID19()