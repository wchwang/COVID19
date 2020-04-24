# Created by woochanghwang at 02/04/2020

import pandas as pd
import networkx as nx
import itertools

def main():
    covid_network_addr = "/home/wch23/Project/LifeArc/COVID-19/result/network/COVID_core_interactions_both_with_covid.tsv"
    covid_sig_gene_addr = "/home/wch23/Project/LifeArc/COVID-19/result/Sig.Genes/covid_ppi_ms_24h.p0.01.txt"

    with open(covid_network_addr) as network_f:
        covid_network_edges = [x.strip().split('\t') for x in network_f.readlines()]

    covid_sig_gene_df = pd.read_csv(covid_sig_gene_addr, sep='\t', names=['Type','Gene'])
    deg_genes = covid_sig_gene_df[covid_sig_gene_df['Type']=='MS']['Gene'].to_list()
    ppi_genes = covid_sig_gene_df[covid_sig_gene_df['Type']=='PPI']['Gene'].to_list()

    covid_G = nx.Graph(covid_network_edges)
    covid_G_nodes = covid_G.nodes()

    # print("{} deg genes {} in network".format(len(deg_genes),len(set(deg_genes)&set(covid_G_nodes))))

    string_addr = "/home/wch23/Project/LifeArc/General/data/STRING/9606.protein.links.v11.0.400.symbol.v2.txt"
    with open(string_addr) as network_f:
        string_network_edges = [x.strip().split('\t') for x in network_f.readlines()]

    string_G = nx.Graph(string_network_edges)
    string_G_nodes= string_G.nodes()

    covid_ppi_deg_pair_list = list(itertools.product(ppi_genes,deg_genes))
    print(len(covid_ppi_deg_pair_list))

    non_reachable_pair = []
    for pair in covid_ppi_deg_pair_list:
        if not nx.has_path(string_G,pair[0],pair[1]):
            non_reachable_pair.append(pair)

    print("{} are not reachalb".format(len(non_reachable_pair)))

    # print("{} deg genes {} in network".format(len(deg_genes),len(set(deg_genes)&set(string_G_nodes))))
    # print("{} ppi genes {} in network".format(len(ppi_genes), len(set(ppi_genes) & set(string_G_nodes))))
    #
    # print(set(deg_genes)-set(string_G_nodes))
    # print(set(ppi_genes) - set(string_G_nodes))

if __name__ == '__main__':
    main()