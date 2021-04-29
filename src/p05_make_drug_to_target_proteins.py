
# Created by woochanghwang at 12/05/2020

import pandas as pd
import networkx as nx
import pickle

def load_obj(file_addr):
    with open(file_addr+ '.pkl', 'rb') as f:
        return pickle.load(f)

def save_obj(obj, file_addr ):
    with open(file_addr + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def make_union_network():
    network_6hr = load_obj(
        "../Result/Network/COVID_All_Structure_All_Shortest_Paths_graph_6hr")
    network_24hr = load_obj(
        "../Result/Network/COVID_All_Structure_All_Shortest_Paths_graph_24hr")
    edges_6hr = list(network_6hr.edges())
    edges_24hr = list(network_24hr.edges())

    edges_union = edges_6hr+edges_24hr
    network_union = nx.Graph(edges_union)

    save_obj(network_union,"../Result/Network/COVID_All_Structure_All_Shortest_Paths_graph_6hr_24hr")

def make_extended_target_list():


    network_union = load_obj("../Result/Network/COVID_All_Structure_All_Shortest_Paths_graph_6hr_24hr")

    candidate_drug_excel = "../Result/Drugs/COVID19_candidate_approved_drugs_combined.xlsx"

    candidate_drug_df = pd.read_excel(candidate_drug_excel,index_col=0)
    print(candidate_drug_df.head())

    candidatge_drug_target_df = candidate_drug_df[['Name','Target Proteins in Network']]

    candidate_drug_extended_target_list = []
    for index, row in candidatge_drug_target_df.iterrows():
        drug_targets = row['Target Proteins in Network'].split(',')
        print(drug_targets)
        print(len(drug_targets))
        if len(drug_targets) <= 5:
            extended_targets = []
            for target in drug_targets:
                extended_targets.append(target)
                extended_targets +=[n for n in network_union.neighbors(target)]


            extended_targets = list(set(extended_targets))
            candidate_drug_extended_target_list.append([index, row['Name'], ','.join(extended_targets)])
        else:
            candidate_drug_extended_target_list.append([index, row['Name'],','.join(drug_targets)])

    candidate_drug_extended_target_df = pd.DataFrame(candidate_drug_extended_target_list, columns=['DrugID','Name','Targets'])
    print(candidate_drug_extended_target_df)

    candidate_drug_extended_target_df.to_csv("../Result/Drugs/COVID19_candidate_approved_drugs_extended_targets.csv",index=False)
        # print(len(list(row['Targets'])))


def main():
    '''
    find enriched pathway
    :return:
    '''
    make_union_network()
    make_extended_target_list()


if __name__ == '__main__':
    main()