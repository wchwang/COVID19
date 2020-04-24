# Created by woochanghwang at 22/04/2020

import toolbox.data_handlers as dh
import pandas as pd
import networkx as nx

def save_RWR_based():
    hour = "24"
    COVID_substructure_RWR_dict = dh.load_obj("/home/wch23/Project/LifeArc/COVID-19/result/network/RWR_based/RWR_results/COVID_All_Structure_RWR")

    COVID_sub_RWR_result_dict = dict()
    for sub, RWR_result in COVID_substructure_RWR_dict.items():
        print(sub)
        # print(RWR_result.head())
        RWR_dict= RWR_result.to_dict()
        COVID_substructure_RWR_dict[sub] = RWR_dict['PageRank']

    COVID_sub_RWR_result_df = pd.DataFrame.from_dict(COVID_substructure_RWR_dict)

    # print(COVID_sub_RWR_result_df.head())
    # print(COVID_sub_RWR_result_df.columns)

    COVID_sub_RWR_result_df['COVID19'] = COVID_sub_RWR_result_df.mean(axis=1)
    print(COVID_sub_RWR_result_df['COVID19'])

    COVID_sub_RWR_result_df.to_csv("/home/wch23/Project/LifeArc/COVID-19/result/network/RWR_based/RWR_results/COVID_All_Structure_RWR_24hr.csv")

    DEP_addr = "/home/wch23/Project/LifeArc/COVID-19/data/DEP/christian_ms_{}h.p0.05.fc0.5.txt".format(hour)
    DEP_df = pd.read_csv(DEP_addr, sep='\t', names=['Gene'])
    DEP_list = DEP_df['Gene'].to_list()

    DIP_to_DEP_df = get_DIP_to_DEP_paths_from_RWR(COVID_sub_RWR_result_df[['COVID19']], DEP_list)
    # print(DIP_to_DEP_df)
    # print(DIP_to_DEP_df.describe())
    DIP_to_DEP_df.to_csv("/home/wch23/Project/LifeArc/COVID-19/result/network/RWR_based/RWR_results/COVID_All_Structure_DIP_to_DEP_24hr.csv")


def get_DIP_to_DEP_paths_from_RWR(RWR_score_df, DEP_list):
    '''
    genes over minimum RWR score of DEP
    :param RWR_score_df:
    :param DEP_list:
    :return:
    '''

    RWR_DEP_df = RWR_score_df[RWR_score_df.index.isin(DEP_list)]
    print(RWR_DEP_df.describe())
    print("Min:",RWR_DEP_df.min())
    print("Mean:",RWR_DEP_df.mean())
    print("quantile:", RWR_DEP_df.quantile(q=0.25))
    min_RWR = float(RWR_DEP_df.min())
    mean_RWR = float(RWR_DEP_df.mean())
    quantile_RWR = float(RWR_DEP_df.quantile(q=0.25))

    DIP_to_DEG_df = RWR_score_df[RWR_score_df['COVID19']>= quantile_RWR]
    return DIP_to_DEG_df

def save_Shortest_path_based():
    from networkx.classes.reportviews import EdgeView
    # data_type = "6hr"
    data_type = "daniel_rnaseq_24h"
    COVID_sub_SP_dict = dh.load_obj("/home/wch23/Project/LifeArc/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_{}".format(data_type))

    sub_sp_edges = []
    for sub, SP_result in COVID_sub_SP_dict.items():
        print(sub)
        for edges in SP_result:
            sub_sp_edges += list(edges)

    #############
    # Add virus to human interaction , Krogan
    ###########
    with open("/home/wch23/Project/LifeArc/COVID-19/data/covid_protein_interactions.tsv") as covid_dip_f:
        covid_dip = [x.strip().split('\t') for x in covid_dip_f.readlines()[1:]]


    sub_sp_edges += covid_dip

    COVID_Graph = nx.Graph(sub_sp_edges)

    print("NODE:" , len(COVID_Graph.nodes()))
    print("EDGES: ", len(COVID_Graph.edges()))

    dh.save_obj(COVID_Graph, "/home/wch23/Project/LifeArc/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_graph_{}".format(data_type))


def main():
    # save_RWR_based()
    save_Shortest_path_based()

if __name__ == '__main__':
    main()