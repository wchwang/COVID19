# Created by woochanghwang at 19/11/2020

# Created by woochanghwang at 27/07/2020

import pandas as pd
import random
import toolbox.data_handler as dh

def target_in_network(drug_targets):
    cip_nodes = pd.read_csv(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/network/SP_based/COVID_6hr_network_nodes.csv")[
        'ID'].to_list()
    drug_targets_nospace = [x.strip()for x in drug_targets]
    # print(drug_targets_nospace)
    return list(set(drug_targets_nospace).intersection(set(cip_nodes)))

def permutation_pvalue_for_drugs():
    covid_cancidiate_drug_df = pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/COVID19_candidate_approved_drugs_combined.v7.xlsx")
    # print(covid_cancidiate_drug_df)
    chembl_drugbank_target_df = pd.read_csv("/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Drug/Drugbank_MTI/drugbank_mti_chembl_tareget_dataframe.csv")
    # print(chembl_drugbank_target_df)

    all_drug_ids = chembl_drugbank_target_df['DrugID'].to_list()

    # print(candidate_drug_ids)
    # candidate_drug_targets_df = chembl_drugbank_target_df[chembl_drugbank_target_df['DrugID'].isin(candidate_drug_ids)]
    # print(candidate_drug_targets_df)
     # chembl_drugbank_target_df['Drug_Target_length']= chembl_drugbank_target_df['Drug_Target'].str.split(',').str.len()

    # chembl_drugbank_target_df = chembl_drugbank_target_df.head()
    chembl_drugbank_target_df['Drug_Targetlist'] = chembl_drugbank_target_df['Drug_Target'].str.split(',')
    chembl_drugbank_target_df['Drug_Target_CIP'] = chembl_drugbank_target_df['Drug_Targetlist'].apply(target_in_network)
    chembl_drugbank_target_df['Drug_Target_CIP_length'] = chembl_drugbank_target_df['Drug_Target_CIP'].str.len()
    print(chembl_drugbank_target_df)
    #
    random_selected_drugs = pd.DataFrame()

    for i in range(1000):
        random_drugs = random.sample(all_drug_ids,k=200)
        # print(len(random_drugs))
        # print(random_drugs)
        # print(len(chembl_drugbank_target_df))
        random_drugs_target_df = chembl_drugbank_target_df[chembl_drugbank_target_df['DrugID'].isin(random_drugs)]
        # print(len(random_drugs_target_df))
        random_selected_drugs = pd.concat([random_selected_drugs,random_drugs_target_df],ignore_index=True)

    # print(len(random_selected_drugs))
    count_drugs = random_selected_drugs.shape[0]
    random_selected_drugs = random_selected_drugs.sort_values(by='Drug_Target_CIP_length',ascending=False)
    random_selected_drugs['Rank'] = random_selected_drugs['Drug_Target_CIP_length'].rank(method='max', ascending=False)
    random_selected_drugs['Target_pvalue'] = random_selected_drugs['Rank']/count_drugs

    print(random_selected_drugs)
    random_selected_drugs.to_csv("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/numberOfTarget_in_CIP_pvalue.csv",index=False)


def make_background_drugs_in_CIP():
    backgroud_drug_df = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/data/CT_n_simulated_drug/COVID_19_simulated_drug_phase2over_v2.xlsx")
    # print(covid_cancidiate_drug_df)
    chembl_drugbank_target_df = pd.read_csv(
        "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Drug/Drugbank_MTI/drugbank_mti_chembl_tareget_dataframe.csv")

    backgroud_drugs = backgroud_drug_df['source'].to_list()

    chembl_drugbank_target_df = chembl_drugbank_target_df[chembl_drugbank_target_df['DrugID'].isin(backgroud_drugs)]
    # print(chembl_drugbank_target_df)
    all_drug_ids = chembl_drugbank_target_df['DrugID'].to_list()

    # # print(candidate_drug_ids)
    # # candidate_drug_targets_df = chembl_drugbank_target_df[chembl_drugbank_target_df['DrugID'].isin(candidate_drug_ids)]
    # # print(candidate_drug_targets_df)
    # # chembl_drugbank_target_df['Drug_Target_length']= chembl_drugbank_target_df['Drug_Target'].str.split(',').str.len()
    #
    # # chembl_drugbank_target_df = chembl_drugbank_target_df.head()
    # #############################################
    chembl_drugbank_target_df['Drug_Targetlist'] = chembl_drugbank_target_df['Drug_Target'].str.split(',')
    chembl_drugbank_target_df['Drug_Target_CIP'] = chembl_drugbank_target_df['Drug_Targetlist'].apply(target_in_network)
    chembl_drugbank_target_df['Drug_Target_CIP_length'] = chembl_drugbank_target_df['Drug_Target_CIP'].str.len()
    print(chembl_drugbank_target_df)
    chembl_drugbank_target_df.to_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/data/CT_n_simulated_drug/COVID_19_simulated_drug_phase2over_CIP_Targets.xlsx",index=False)

def permutation_pvalue_for_targets():
    # covid_cancidiate_drug_df = pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/COVID19_candidate_approved_drugs_combined.v7.xlsx")
    # # print(covid_cancidiate_drug_df)
    # chembl_drugbank_target_df = pd.read_csv("/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Drug/Drugbank_MTI/drugbank_mti_chembl_tareget_dataframe.csv")
    # # print(chembl_drugbank_target_df)
    #
    # all_drug_ids = chembl_drugbank_target_df['DrugID'].to_list()
    #
    # # print(candidate_drug_ids)
    # # candidate_drug_targets_df = chembl_drugbank_target_df[chembl_drugbank_target_df['DrugID'].isin(candidate_drug_ids)]
    # # print(candidate_drug_targets_df)
    #  # chembl_drugbank_target_df['Drug_Target_length']= chembl_drugbank_target_df['Drug_Target'].str.split(',').str.len()
    #
    # # chembl_drugbank_target_df = chembl_drugbank_target_df.head()
    # #############################################
    # chembl_drugbank_target_df['Drug_Targetlist'] = chembl_drugbank_target_df['Drug_Target'].str.split(',')
    # chembl_drugbank_target_df['Drug_Target_CIP'] = chembl_drugbank_target_df['Drug_Targetlist'].apply(target_in_network)
    # chembl_drugbank_target_df['Drug_Target_CIP_length'] = chembl_drugbank_target_df['Drug_Target_CIP'].str.len()
    # print(chembl_drugbank_target_df)

    # make_background_drugs_in_CIP()
    backgound_drugs_df = pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/data/CT_n_simulated_drug/COVID_19_simulated_drug_phase2over_CIP_Targets.xlsx")
    backgound_drugs_df = backgound_drugs_df[backgound_drugs_df['Drug_Target_CIP_length']>0]
    backgound_drugs_df['Drug_Target_CIP'] = backgound_drugs_df['Drug_Target_CIP'].str.replace("[", "")
    backgound_drugs_df['Drug_Target_CIP'] = backgound_drugs_df['Drug_Target_CIP'].str.replace("]", "")
    all_drug_ids = backgound_drugs_df['DrugID'].to_list()
    random_selected_drugs = pd.DataFrame()


    #########################
    # random_drugs = random.sample(all_drug_ids, k=200)
    # random_drugs_target_df = backgound_drugs_df[backgound_drugs_df['DrugID'].isin(random_drugs)]
    # print(list(random_drugs_target_df))
    # random_drugs_target_df = random_drugs_target_df[['DrugID','Drug_Target_CIP']]
    # random_drugs_target_df['Drug_Target_CIP'] = random_drugs_target_df['Drug_Target_CIP'].str.split(',')
    # random_drugs_target_df = random_drugs_target_df.Drug_Target_CIP.apply(pd.Series)\
    #     .merge(random_drugs_target_df, right_index = True, left_index = True)\
    #     .drop('Drug_Target_CIP',axis=1).melt(id_vars = ['DrugID'],value_name = "Drug_Target_CIP")\
    #     .drop("variable",axis=1).dropna()
    #
    # random_drugs_groupby_target = random_drugs_target_df.groupby('Drug_Target_CIP').agg(list)
    # random_drugs_groupby_target['Drug_count'] = random_drugs_groupby_target['DrugID'].str.len()
    # random_drugs_groupby_target = random_drugs_groupby_target.reset_index()
    # print(random_drugs_target_df)
    # print(random_drugs_groupby_target)
    #########################

    for i in range(100):
        random_drugs = random.sample(all_drug_ids, k=200)
        random_drugs_target_df = backgound_drugs_df[backgound_drugs_df['DrugID'].isin(random_drugs)]
        # print(list(random_drugs_target_df))
        random_drugs_target_df = random_drugs_target_df[['DrugID','Drug_Target_CIP']]
        random_drugs_target_df['Drug_Target_CIP'] = random_drugs_target_df['Drug_Target_CIP'].str.split(',')
        random_drugs_target_df = random_drugs_target_df.Drug_Target_CIP.apply(pd.Series)\
            .merge(random_drugs_target_df, right_index = True, left_index = True)\
            .drop('Drug_Target_CIP',axis=1).melt(id_vars = ['DrugID'],value_name = "Drug_Target_CIP")\
            .drop("variable",axis=1).dropna()

        random_drugs_groupby_target = random_drugs_target_df.groupby('Drug_Target_CIP').agg(list)
        random_drugs_groupby_target['Drug_count'] = random_drugs_groupby_target['DrugID'].str.len()
        random_drugs_groupby_target = random_drugs_groupby_target.reset_index()
        random_selected_drugs = pd.concat([random_selected_drugs,random_drugs_groupby_target],ignore_index=True)

    # print(len(random_selected_drugs))
    count_drugs = random_selected_drugs.shape[0]
    random_selected_drugs = random_selected_drugs.sort_values(by='Drug_count',ascending=False)
    random_selected_drugs['Rank'] = random_selected_drugs['Drug_count'].rank(method='max', ascending=False)
    random_selected_drugs['Target_pvalue'] = random_selected_drugs['Rank']/count_drugs

    print(random_selected_drugs)
    random_selected_drugs.to_csv("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/target_numberOfdrugs_in_CIP_pvalue.csv",index=False)

def get_one_neighbor( proteins):
    cip_network = dh.load_obj(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_graph_6hr")
    one_neighbors = []

    proteins_nospace = [x.strip() for x in proteins]
    try:
        for protein in proteins_nospace:
            one_neighbors += [n for n in cip_network.neighbors(protein)]
    except :
        print("Not in this network")
    return len(set(one_neighbors))

def check_network_coverage():
    covid_cancidiate_drug_df = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/COVID19_candidate_approved_drugs_combined.v7.xlsx")

    candidate_drugs_df = covid_cancidiate_drug_df[['Drug_ID','Name','Target Proteins in Network']]
    candidate_drugs_df['Target_Proteins'] = candidate_drugs_df['Target Proteins in Network'].str.split(',')
    candidate_drugs_df['Target_Proteins_length'] = candidate_drugs_df['Target_Proteins'].str.len()
    # candidate_drugs_df['Network_coverage'] = candidate_drugs_df['Target_Proteins'].apply(get_one_neighbor)
    #
    # print(candidate_drugs_df)
    # candidate_drugs_df.to_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/target_neighbor_network_coverage.xlsx",index=False)
    cip_network = dh.load_obj(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_graph_6hr")
    print(len(cip_network.nodes()))

def hypergeomatric_test_for_drug_numberOfTargets():
    pass


def main():
    # permutation_pvalue_for_drugs()
    permutation_pvalue_for_targets()
    # check_network_coverage()

    hypergeomatric_test_for_drug_numberOfTargets()
if __name__ == '__main__':
    main()