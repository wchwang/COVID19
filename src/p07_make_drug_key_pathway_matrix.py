# Created by woochanghwang at 17/07/2020
'''
Make Drug - enrichmemt pathway score ( F1 score) matrix
all pathways: 206
- Remvoe Top level pathways(Metabolism of lipids, Apoptosis et al)
'''
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def make_low_level_pathways_matrix(low_level_pathways):
    # 200 drug list
    candidate_drug_addr = "/COVID-19/result/Drug/COVID19_candidate_approved_drugs_combined.v5.multimoa.category.subcategory.xlsx"
    # pathway enrichment result
    pathway_addr_prefix = "/COVID-19/result/Drug/6hr_24hr/6-2_24-6/Drug_pathways/enriched_reactome_pathways_"

    candidate_drug_df = pd.read_excel(candidate_drug_addr)



    sorted_pathway = sorted(low_level_pathways)
    print(len(sorted_pathway))
    print(sorted_pathway[:3])

    drug_names = candidate_drug_df['Name'].to_list()

    drug_pathway_matrix = []

    for drug in drug_names:
        drug_pathway_enriched_df = pd.read_csv(pathway_addr_prefix + drug + '.csv')
        drug_pathway_values = [drug]
        # print(drug)
        # print(drug_pathway_enriched_df.head())
        for path in sorted_pathway:
            path_df = drug_pathway_enriched_df[drug_pathway_enriched_df['term_name'] == path][['precision', 'recall']]
            if len(path_df) > 0:
                recall = path_df.iloc[0]['recall']
                precision = path_df.iloc[0]['precision']
                f1_score = 2 * ((precision * recall) / (precision + recall))
                drug_pathway_values.append(f1_score)
            else:
                drug_pathway_values.append(0.0)
        drug_pathway_matrix.append(drug_pathway_values)

    drug_pathway_df = pd.DataFrame(drug_pathway_matrix, columns=['Drug_name'] + sorted_pathway)
    print(drug_pathway_df)

    # sns.clustermap(drug_pathway_df)
    #
    # plt.show()

    drug_pathway_df.to_csv(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/candidate_drug_to_pathways_matrix_low_level.csv",
        index=False)


def find_low_level_pathways(moa_to_pathway_df):

    top2_groupby = moa_to_pathway_df.groupby('Top2')['Pathways'].agg(list).reset_index(name='Pathways')
    top2_groupby_size = moa_to_pathway_df.groupby('Top2')['Pathways'].size().reset_index(name='counts')
    top2_pathways_size_df = top2_groupby.merge(top2_groupby_size)

    high_level_top2_df = top2_pathways_size_df[top2_pathways_size_df['counts']>1]
    high_level_top2 = high_level_top2_df['Top2'].to_list()
    high_level_top2 = [x.strip() for x in high_level_top2]
    # print(high_level_top2)

    moa_to_pathway_df = moa_to_pathway_df[['Pathways','Parents','Top2']]
    moa_to_pathway_df['Parents_counts'] = moa_to_pathway_df['Parents'].str.split('|').str.len()

    low_level_pathways_df = moa_to_pathway_df.dropna()


    low_level_pathways_df = low_level_pathways_df[~low_level_pathways_df['Pathways'].isin(high_level_top2)]


    return low_level_pathways_df['Pathways'].to_list()


def make_drug_enriched_pathways(candidate_drug_moa_df):
    '''
    Pandas : splitting a colum into multipul rows
    :param candidate_drug_moa_df:
    :return:
    '''

    print(candidate_drug_moa_df)
    print(candidate_drug_moa_df.columns)
    drug_pathway_df = candidate_drug_moa_df[['Name','Pathway']]
    drug_MoA_df = candidate_drug_moa_df[['Name','Multi_MoA_Sub_Category']]

    print(drug_MoA_df)

    drug_pathway_split_df = pd.DataFrame(drug_pathway_df.Pathway.str.split('|').tolist(),index=drug_pathway_df.Name).stack()

    drug_pathway_split_df = drug_pathway_split_df.reset_index([0,'Name'])
    drug_pathway_split_df.columns = ['Name','Pathways']


    # drug_pathway_split_df = drug_pathway_split_df.set_index('Name')
    print(drug_pathway_split_df)

    drug_pathway_moa_df = drug_pathway_split_df.merge(drug_MoA_df, how='left', left_on='Name', right_on='Name')
    print(drug_pathway_moa_df)
    print(drug_pathway_moa_df.columns)

    drug_pathway_moa_df.to_excel("/COVID-19/result/Drug/Drug_enriched_pathways.xlsx", index=False)

def main():

    candidate_drug_df = pd.read_excel("/COVID-19/result/Drug/COVID19_candidate_approved_drugs_combined.xlsx")
    make_drug_enriched_pathways(candidate_drug_df)
    moa_to_pathway_df = pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/Drug/MoA/MoA_Reactome_pathways_category.xlsx")
    low_level_pathways = find_low_level_pathways(moa_to_pathway_df)

    make_low_level_pathways_matrix(low_level_pathways)


if __name__ == '__main__':
    main()
