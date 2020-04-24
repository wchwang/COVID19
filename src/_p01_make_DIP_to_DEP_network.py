# Created by woochanghwang at 21/04/2020

import pandas as pd
import networkx as nx
import itertools
def run_pagerank(graph_G, DIP_genes, DEP_genes, RWR_result_addr):

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


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def make_shortest_path_sources_to_targets(covid_sig_gene_addr):
    '''
    To make network between covid ppi to covid DEG genes
    All pairs are connected on STRING
    :return:
    '''

    # covid_sig_gene_addr = "/home/wch23/Project/LifeArc/COVID-19/data/covid_ppi_ms_24h.p0.01.txt"

    covid_sig_gene_df = pd.read_csv(covid_sig_gene_addr, sep='\t', names=['Type', 'Gene'])
    deg_genes = covid_sig_gene_df[covid_sig_gene_df['Type'] == 'MS']['Gene'].to_list()
    ppi_genes = covid_sig_gene_df[covid_sig_gene_df['Type'] == 'PPI']['Gene'].to_list()

    string_addr = "/home/wch23/Project/LifeArc/General/data/STRING/9606.protein.links.v11.0.400.symbol.v2.txt"
    with open(string_addr) as network_f:
        string_network_edges = [x.strip().split('\t') for x in network_f.readlines()]

    string_G = nx.Graph(string_network_edges)

    covid_ppi_deg_pair_list = list(itertools.product(ppi_genes, deg_genes))
    print(len(covid_ppi_deg_pair_list))

    pair_shortest_path=[]
    all_path_list = []
    for pair in covid_ppi_deg_pair_list:
        # all_shorttest_paths = [p for p in nx.all_shortest_paths(string_G, source=pair[0], target=pair[1])]
        # print(all_shorttest_paths)
        for p in nx.all_shortest_paths(string_G, source=pair[0], target=pair[1]):
            pairs_from_p = pairwise(p)
            all_path_list += pairs_from_p

    print(all_path_list)

    all_path_G = nx.Graph(all_path_list)
    all_edges = all_path_G.edges()

    print(all_edges)
    all_edges_w = ['\t'.join(edge) for edge in all_edges]

    with open("/home/wch23/Project/LifeArc/COVID-19/result/network/ppi_to_deg_paths.tsv",'w') as ppi_to_deg_f:
        ppi_to_deg_f.write('\n'.join(all_edges_w))

    return all_edges
def main():
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

    for sub_structure in covid_sub_structures[:2]:

        string_G = nx.Graph(string_network_edges[1:])
        DIP_sub_list = DIP_df[DIP_df['Bait']== sub_structure]['StringGene'].to_list()
        print("SUB : {}".format(sub_structure))
        print("DIP: ", len(DIP_sub_list) )
        RWR_result_addr = "/home/wch23/Project/LifeArc/COVID-19/result/network/RWR_based/RWR_results/COVID_{}_RWR_BC_weighted_result_RWR.tsv".format(
            sub_structure)
        RWR_result_df = run_pagerank(string_G, DIP_sub_list, DEP_list, RWR_result_addr )
        DIP_to_DEP_df = get_DIP_to_DEP_paths(RWR_result_df, DEP_list)

if __name__ == '__main__':
    main()