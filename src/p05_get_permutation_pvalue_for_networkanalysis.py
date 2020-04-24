# Created by woochanghwang at 23/04/2020

import pandas as pd
import glob
import os
import toolbox.data_handlers as dh


def get_permutation_result(permutation_result, key_method):
    # 'eigen': eigenvector_centrality_dict
    # 'between': between_centrality_dict,
    # 'degee': degree_centrality_dict,
    # 'rwr': rwr_dict
    score_list = []
    for permutation_process in permutation_result:
        network_analisys_score_dict = permutation_process[key_method]
        score_list += list(network_analisys_score_dict.values())

    print(score_list[:5])
    print(len(score_list))

    permutation_df= pd.DataFrame(score_list, columns=['Score'])
    permutation_df= permutation_df.sort_values('Score',ascending=False)
    return permutation_df


def get_permutation_result_df(permutation_result, key_method):
    # 'eigen': eigenvector_centrality_dict
    # 'between': between_centrality_dict,
    # 'degee': degree_centrality_dict,
    # 'rwr': rwr_dict
    score_list = []
    for permutation_process in permutation_result:
        network_analisys_score_dict = permutation_process[key_method]
        score_list += list(network_analisys_score_dict.values())

    print(score_list[:5])
    print(len(score_list))

    permutation_df= pd.DataFrame(score_list, columns=[key_method])
    # permutation_df[:5].apply(lambda row: print(row[key_method]), axis=1)
    # permutation_df= permutation_df.sort_values(key_method,ascending=False)
    # permutation_df[:5].apply(lambda row: print(row[key_method]), axis=1)
    return permutation_df


def get_pvalue(df_permutation_results_scaling_sorted,p_value):


    # p_value = 0.01
    p_value_index = int(len(df_permutation_results_scaling_sorted.index)*p_value)
    significant_score = df_permutation_results_scaling_sorted.iloc[[p_value_index], [0]]
    # print("pivalue index:", p_value_index)
    print("sifgnicant score of path:", df_permutation_results_scaling_sorted.iloc[[p_value_index],[0]])
    return significant_score


def main():
    cancer_type = "COVID"
    start_gene = "COVID19"

    # round_num = "round6"
    # data_type = "24hr"
    round_num = "round2"
    data_type = "6hr"

    permutation_result_all = dh.load_obj("/mnt/raid0_data/wch23/Project/COVID-19/result/Permutation_NA/{}/{}/COVID_network_analysis_permutation_{}".format(data_type, round_num, data_type))

    permutation_eigen_result_all = get_permutation_result_df(permutation_result_all, 'eigen')
    permutation_degree_result_all = get_permutation_result_df(permutation_result_all, 'degee')
    permutation_between_result_all = get_permutation_result_df(permutation_result_all, 'between')
    permutation_rwr_result_all = get_permutation_result_df(permutation_result_all, 'rwr')
    print(permutation_eigen_result_all.head())
    centrality_addr = "/mnt/raid0_data/wch23/Project/COVID-19/result/Centrality/{}/{}/COVID19_centrality_RWR_result_{}_{}.csv".format(
        data_type, round_num, data_type, round_num)
    df_centrality_result = pd.read_csv(centrality_addr, index_col=0)
    print(df_centrality_result.head())
    centrality_pvalue_addr = "/mnt/raid0_data/wch23/Project/COVID-19/result/Centrality/{}/{}/{}_{}_gene_score_by_centrality_pvalue.csv".format(
        data_type, round_num, cancer_type, data_type)

    df_final_result = pd.concat([permutation_eigen_result_all, permutation_degree_result_all, permutation_between_result_all, permutation_rwr_result_all], axis=1)
    # print(len(df_final_result[df_final_result.between > 1.445]))


    df_centrality_result['eigen_pvalue'] = df_centrality_result.apply(lambda row: len(df_final_result[df_final_result.eigen > row.eigen]) / len(
            df_final_result.index),axis=1)
    df_centrality_result['degree_pvalue'] = df_centrality_result.apply(
        lambda row: len(df_final_result[df_final_result.degee > row.degre]) / len(
            df_final_result.index), axis=1)

    df_centrality_result['between_pvalue'] = df_centrality_result.apply(
        lambda row: len(df_final_result[df_final_result.between > row['between']]) / len(
            df_final_result.index), axis=1)

    df_centrality_result['rwr_pvalue'] = df_centrality_result.apply(
        lambda row: len(df_final_result[df_final_result.rwr > row.RWR]) / len(
            df_final_result.index), axis=1)
    print(df_centrality_result)
    df_centrality_result.to_csv(centrality_pvalue_addr, index=False)


    # df_final_result = pd.DataFrame(
    #     columns=['Gene', 'Eigen', 'Degree', 'Between', 'RWR', 'Eigen_pvalue', 'Degree_pvalue', 'Between_plvaue',
    #              'RWR_pvalue'])

    # a=0
    # for gene in df_centrality_result.index:
    #     # print(len(permutation_result_all[permutation_result_all.score >= df_rwr_result.loc[gene,'RWR']]))
    #     # print(permutation_result_all[permutation_result_all.score >= df_rwr_result.loc[gene,'RWR']])
    #     # print(len(permutation_result_all.index))
    #     # print(df_rwr_result.loc[gene,'RWR'])
    #     a += 1
    #     if a % 100 ==0 : print(a)
    #     beween_pvalue_gene = len(permutation_between_result_all[permutation_between_result_all.Score > df_centrality_result.loc[gene, 'between']]) / len(
    #         permutation_between_result_all.index)
    #     degree_pvalue_gene = len(permutation_degree_result_all[
    #                                  permutation_degree_result_all.Score > df_centrality_result.loc[
    #                                      gene, 'degre']]) / len(permutation_degree_result_all.index)
    #     eigen_pvalue_gene = len(permutation_eigen_result_all[
    #                                  permutation_eigen_result_all.Score > df_centrality_result.loc[
    #                                      gene, 'eigen']]) / len(permutation_eigen_result_all.index)
    #     rwr_pvalue_gene = len(permutation_rwr_result_all[
    #                                 permutation_rwr_result_all.Score > df_centrality_result.loc[
    #                                     gene, 'RWR']]) / len(permutation_rwr_result_all.index)
    #
    #     df_final_result = df_final_result.append(
    #         pd.Series([gene,
    #                    df_centrality_result.loc[gene, 'eigen'],
    #                    df_centrality_result.loc[gene, 'degre'],
    #                    df_centrality_result.loc[gene, 'between'],
    #                    df_centrality_result.loc[gene, 'RWR'],
    #                    eigen_pvalue_gene,
    #                    degree_pvalue_gene,
    #                    beween_pvalue_gene,
    #                    rwr_pvalue_gene],
    #                   index=['Gene', 'Eigen', 'Degree', 'Between', 'RWR', 'Eigen_pvalue', 'Degree_pvalue', 'Between_plvaue','RWR_pvalue']), ignore_index=True)
    #
    # print(df_final_result)
    # df_final_result.to_csv(centrality_pvalue_addr)


def main_COVID19():
    cancer_type = "COVID"
    start_gene = "COVID19"

    round_num = "round6"
    data_type = "24hr"

    permutation_result_all = dh.load_obj("/mnt/raid0_data/wch23/Project/COVID-19/result/Permutation_NA/{}/{}/COVID_network_analysis_permutation_{}".format(data_type, round_num, data_type))

    permutation_eigen_result_all = get_permutation_result(permutation_result_all, 'eigen')
    permutation_degree_result_all = get_permutation_result(permutation_result_all, 'degee')
    permutation_between_result_all = get_permutation_result(permutation_result_all, 'between')
    permutation_rwr_result_all = get_permutation_result(permutation_result_all, 'rwr')
    print(permutation_eigen_result_all.head())
    centrality_addr = "/mnt/raid0_data/wch23/Project/COVID-19/result/Centrality/{}/{}/COVID19_centrality_RWR_result_{}_{}.csv".format(
        data_type, round_num, data_type, round_num)
    df_centrality_result = pd.read_csv(centrality_addr, index_col=0)
    print(df_centrality_result.head())
    centrality_pvalue_addr = "/mnt/raid0_data/wch23/Project/COVID-19/result/Centrality/{}/{}/{}_{}_gene_score_by_centrality_pvalue.csv".format(
        data_type, round_num, cancer_type, data_type)
    df_final_result = pd.DataFrame(
        columns=['Gene', 'Eigen', 'Degree', 'Between', 'RWR', 'Eigen_pvalue', 'Degree_pvalue', 'Between_plvaue',
                 'RWR_pvalue'])

    a=0
    for gene in df_centrality_result.index:
        # print(len(permutation_result_all[permutation_result_all.score >= df_rwr_result.loc[gene,'RWR']]))
        # print(permutation_result_all[permutation_result_all.score >= df_rwr_result.loc[gene,'RWR']])
        # print(len(permutation_result_all.index))
        # print(df_rwr_result.loc[gene,'RWR'])
        a += 1
        if a % 100 ==0 : print(a)
        beween_pvalue_gene = len(permutation_between_result_all[permutation_between_result_all.Score > df_centrality_result.loc[gene, 'between']]) / len(
            permutation_between_result_all.index)
        degree_pvalue_gene = len(permutation_degree_result_all[
                                     permutation_degree_result_all.Score > df_centrality_result.loc[
                                         gene, 'degre']]) / len(permutation_degree_result_all.index)
        eigen_pvalue_gene = len(permutation_eigen_result_all[
                                     permutation_eigen_result_all.Score > df_centrality_result.loc[
                                         gene, 'eigen']]) / len(permutation_eigen_result_all.index)
        rwr_pvalue_gene = len(permutation_rwr_result_all[
                                    permutation_rwr_result_all.Score > df_centrality_result.loc[
                                        gene, 'RWR']]) / len(permutation_rwr_result_all.index)

        df_final_result = df_final_result.append(
            pd.Series([gene,
                       df_centrality_result.loc[gene, 'eigen'],
                       df_centrality_result.loc[gene, 'degre'],
                       df_centrality_result.loc[gene, 'between'],
                       df_centrality_result.loc[gene, 'RWR'],
                       eigen_pvalue_gene,
                       degree_pvalue_gene,
                       beween_pvalue_gene,
                       rwr_pvalue_gene],
                      index=['Gene', 'Eigen', 'Degree', 'Between', 'RWR', 'Eigen_pvalue', 'Degree_pvalue', 'Between_plvaue','RWR_pvalue']), ignore_index=True)

    print(df_final_result)
    df_final_result.to_csv(centrality_pvalue_addr)



if __name__ == '__main__':
    # main_COVID19()
    main()