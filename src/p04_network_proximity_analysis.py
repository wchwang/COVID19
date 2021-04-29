# Created by woochanghwang at 18/11/2020


import sys, time
import pandas as pd
import numpy
import networkx, random
import mygene
import multiprocessing as mp


def check_has_path(G,node_from,node_to):
    return(networkx.has_path(G,node_from,node_to))

def get_shortest_path_length_between(G, source_id, target_id):
    return networkx.shortest_path_length(G, source_id, target_id)

def calculate_closest_distance(network, nodes_from, nodes_to, lengths=None):
    values_outer = []
    if lengths is None:
        for node_from in nodes_from:
            values = []
            for node_to in nodes_to:
                # print("from - to", node_from, node_to)
                if not check_has_path(network,node_from,node_to): continue
                val = get_shortest_path_length_between(network, node_from, node_to)
                values.append(val)
            if len(values) == 0:    continue
            d = min(values)
            # print (d)
            values_outer.append(d)
    else:
        for node_from in nodes_from:
            values = []
            vals = lengths[node_from]
            for node_to in nodes_to:
                val = vals[node_to]
                values.append(val)
            d = min(values)
            values_outer.append(d)
    d = numpy.mean(values_outer)
    # print d
    return d

def get_degree_binning(g, bin_size, lengths=None):
    '''

    It tried to make number of gene list( with same degree) to bin size.
    If number of gene list with some degree , it combine with ohter genes with another degree to meet bin size.
    '''
    degree_to_nodes = {}
    for node, degree in g.degree():  # .items(): # iterator in networkx 2.0
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)
    values = degree_to_nodes.keys()
    # values.sort() #pytyon 2.x
    values = sorted(values)
    bins = []
    i = 0
    while i < len(values):
        low = values[i]
        val = degree_to_nodes[values[i]]
        while len(val) < bin_size:
            i += 1
            if i == len(values):
                break
            val.extend(degree_to_nodes[values[i]])
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1
        # print i, low, high, len(val)
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
    return bins

def get_degree_equivalents(seeds, bins, g):
    seed_to_nodes = {}
    for seed in seeds:
        d = g.degree(seed)
        for l, h, nodes in bins:
            if l <= d and h >= d:
                mod_nodes = list(nodes)
                mod_nodes.remove(seed)
                seed_to_nodes[seed] = mod_nodes
                break
    return seed_to_nodes

def pick_random_nodes_matching_selected(network, bins, nodes_selected, n_random, degree_aware=True, connected=False,
                                        seed=None):
    """
    Use get_degree_binning to get bins
    """
    if seed is not None:
        random.seed(seed)
    values = []
    nodes = network.nodes()
    for i in range(n_random):
        if degree_aware:
            if connected:
                raise ValueError("Not implemented!")
            nodes_random = set()
            node_to_equivalent_nodes = get_degree_equivalents(nodes_selected, bins, network)
            for node, equivalent_nodes in node_to_equivalent_nodes.items():
                # nodes_random.append(random.choice(equivalent_nodes))
                chosen = random.choice(equivalent_nodes)
                for k in range(20):  # Try to find a distinct node (at most 20 times)
                    if chosen in nodes_random:
                        chosen = random.choice(equivalent_nodes)
                nodes_random.add(chosen)
            nodes_random = list(nodes_random)
        else:
            if connected:
                nodes_random = [random.choice(nodes)]
                k = 1
                while True:
                    if k == len(nodes_selected):
                        break
                    node_random = random.choice(nodes_random)
                    node_selected = random.choice(network.neighbors(node_random))
                    if node_selected in nodes_random:
                        continue
                    nodes_random.append(node_selected)
                    k += 1
            else:
                nodes_random = random.sample(nodes, len(nodes_selected))
        values.append(nodes_random)
    return values

def get_random_nodes(nodes, network, bins=None, n_random=1000, min_bin_size=100, degree_aware=True, seed=None):
    if bins is None:
        # Get degree bins of the network
        bins = get_degree_binning(network, min_bin_size)
    nodes_random = pick_random_nodes_matching_selected(network, bins, nodes, n_random, degree_aware,
                                                                         seed=seed)
    return nodes_random

def calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=None,
                        n_random=1000, min_bin_size=100, seed=452456, lengths=None):
    """
    Calculate proximity from nodes_from to nodes_to
    If degree binning or random nodes are not given, they are generated
    lengths: precalculated shortest path length dictionary
    """

    nodes_network = set(network.nodes())
    nodes_from = set(nodes_from) & nodes_network
    nodes_to = set(nodes_to) & nodes_network
    if len(nodes_from) == 0 or len(nodes_to) == 0:
        return None  # At least one of the node group not in network
    d = calculate_closest_distance(network, nodes_from, nodes_to, lengths)
    if bins is None and (nodes_from_random is None or nodes_to_random is None):
        bins = get_degree_binning(network, min_bin_size,
                                                    lengths)  # if lengths is given, it will only use those nodes
    if nodes_from_random is None:
        nodes_from_random = get_random_nodes(nodes_from, network, bins=bins, n_random=n_random,
                                             min_bin_size=min_bin_size, seed=seed)
    if nodes_to_random is None:
        nodes_to_random = get_random_nodes(nodes_to, network, bins=bins, n_random=n_random, min_bin_size=min_bin_size,
                                           seed=seed)
    random_values_list = zip(nodes_from_random, nodes_to_random)
    values = numpy.empty(len(nodes_from_random))  # n_random
    for i, values_random in enumerate(random_values_list):
        nodes_from, nodes_to = values_random
        # values[i] = network_utilities.get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
        values[i] = calculate_closest_distance(network, nodes_from, nodes_to, lengths)
    # pval = float(sum(values <= d)) / len(values) # needs high number of n_random
    m, s = numpy.mean(values), numpy.std(values)
    if s == 0:
        z = 0.0
    else:
        z = (d - m) / s
    return d, z, (m, s)  # (z, pval)

def calculate_proximity_multiple_multiprocessing_return(input_tuple):
    """
    Run proximity on each entries of from and to files in a pairwise manner
    output is saved in out_file (e.g., output.txt)
    disease_mode : whole - from disease_gene.tsv
    disease_mode : disease_name - from our analysis
    """
    network, disease_to_genes, drug_to_targets = input_tuple
    drug = drug_to_targets[0]
    nodes_from = drug_to_targets[1]

    start_time = time.time()

    n_random = 1000
    min_bin_size = 100
    seed = 452456
    lengths = None

    # Get degree binning
    bins = get_degree_binning(network, min_bin_size)

    proximity_result = []

    for disease, nodes_to in disease_to_genes.items():
        # print(drug, disease)
        d, z, (m, s) = calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None,
                                           nodes_to_random=None, bins=bins, n_random=n_random,
                                           min_bin_size=min_bin_size, seed=seed, lengths=lengths)
        proximity_result.append((drug, disease, len(nodes_from), len(nodes_to), d, z))


    print("{} times training Runtime: {:2f} Minutes".format(drug, ((time.time() - start_time) / 60)))
    # f.close()
    return  proximity_result[0]

def get_nodes_and_edges_from_sif_file(file_name, store_edge_type=False, delim=None, data_to_float=True):
    """
    Parse sif file into node and edge sets and dictionaries
    returns setNode, setEdge, dictNode, dictEdge
    store_edge_type: if True, dictEdge[(u,v)] = edge_value
    delim: delimiter between elements in sif file, if None all whitespaces between letters are considered as delim
    """
    setNode = set()
    setEdge = set()
    dictNode = {}
    dictEdge = {}
    flag = False
    f = open(file_name)
    for line in f:
        if delim is None:
            words = line.rstrip("\n").split()
        else:
            words = line.rstrip("\n").split(delim)
        id1 = words[0]
        setNode.add(id1)
        if len(words) == 2:
            if data_to_float:
                score = float(words[1])
            else:
                score = words[1]
            dictNode[id1] = score
        elif len(words) >= 3:
            if len(words) > 3:
                flag = True
            id2 = words[2]
            setNode.add(id2)
            setEdge.add((id1, id2))
            if store_edge_type:
                if data_to_float:
                    dictEdge[(id1, id2)] = float(words[1])
                else:
                    dictEdge[(id1, id2)] = words[1]
    f.close()
    if len(setEdge) == 0:
        setEdge = None
    if len(dictNode) == 0:
        dictNode = None
    if len(dictEdge) == 0:
        dictEdge = None
    if flag:
        print("Warning: Ignored extra columns in the file!")
    return setNode, setEdge, dictNode, dictEdge

def create_graph(directed=False):
    """
        Creates & returns a graph
    """
    if directed:
        g = networkx.DiGraph()
    else:
        g = networkx.Graph()
    return g

def create_network_from_sif_file(network_file_in_sif, use_edge_data=False, delim=None, include_unconnected=True):
    setNode, setEdge, dictDummy, dictEdge = get_nodes_and_edges_from_sif_file(network_file_in_sif,
                                                                              store_edge_type=use_edge_data,
                                                                              delim=delim)
    g = create_graph()
    if include_unconnected:
        g.add_nodes_from(setNode)
    if use_edge_data:
        for e, w in dictEdge.items():
            u, v = e
            g.add_edge(u, v, w=w)  # ,{'w':w})
    else:
        g.add_edges_from(setEdge)
    return g

def get_connected_components(G, return_as_graph_list=True):
    """
        Finds (strongly in the case of directed network) connected components of graph
        returnAsGraphList: returns list of graph objects corresponding to connected components (from larger to smaller)
        otherwise returns list of node list corresponding nodes in connected components
    """
    result_list = []

    if return_as_graph_list:
        result_list = networkx.connected_component_subgraphs(G)
    else:
        result_list = [c for c in sorted(networkx.connected_components(G), key=len, reverse=True)]

    return result_list

def get_drug_target_drugbank(drug_target_file, nodes=None, network=None):
    """
    If nodes is not None, keep only nodes in the network
    If network is not None, keep only LCC
    drugbank file = [[drugid,entrez],[]...]
    groupby drugid
    return = drug = gene_set
    """

    drug_target_df = pd.read_table(drug_target_file,sep='\t')
    drug_target_df = drug_target_df.applymap(str)
    drug_targets_df = drug_target_df.groupby('DrugID').agg({'Drug_Target':list})
    print(drug_targets_df.head())
    drug_targets_dict = drug_targets_df.to_dict()['Drug_Target']

    drug_to_genes_dict = {}
    ## Check weather targets are in network

    for drug, genes in drug_targets_dict.items():
        genes = set(genes)
        if nodes is not None:
            genes &= nodes
            if len(genes) == 0:
                continue
        if network is not None:
            network_sub = network.subgraph(genes)
            genes = get_connected_components(network_sub, False)[0]
        drug_to_genes_dict[drug] = genes


    return drug_to_genes_dict

def mapping_genes_id_to(gene_list,id_from, id_to, species = 9606, as_dataframe = False):
    mg = mygene.MyGeneInfo()
    return mg.querymany(gene_list, scopes=id_from ,fields= id_to ,species=species,as_dataframe=as_dataframe)

def get_diseasome_genes_from_selectedGenes(gene_list_file, disease_name, nodes=None, network=None ):
    """
    It is exactly same function as get_diseasome_genes
    from gene list from my analysis
    one disease type
    :param gene_list_file: result from my anlysis (RWR)
    :param disease_name: pre difined
    :param nodes: nodes in network
    :param network:
    :return: dict(disease_to_gene[disease]={genes})
    """
    gene_df = pd.read_table(gene_list_file, sep='\t')
    gene_list = list(gene_df['Gene'])
    print("gene list",gene_list)
    id_from='symbol'
    id_to = 'entrezgene'
    geneID_df = mapping_genes_id_to(gene_list,id_from=id_from, id_to=id_to, species = 9606, as_dataframe = True)
    print(geneID_df.head())

    geneID_list = list(geneID_df[id_to])

    disease_to_genes = dict()
    disease_to_category = dict()

    # consider genes in only network
    genes = set(geneID_list)
    if nodes is not None:
        genes &= nodes
        # if len(genes) == 0:
        #     continue
    if network is not None:
        network_sub = network.subgraph(genes)
        genes = get_connected_components(network_sub, False)[0]
    disease_to_genes[disease_name] = genes


    return disease_to_genes

def mp_drug_simulation(disease_name, network_file, gene_list_file, drug_target_file, output_file, num_process):
    network = create_network_from_sif_file(network_file, use_edge_data=False, delim=None,
                                              include_unconnected=True)
    nodes = set(network.nodes())
    drug_to_targets = get_drug_target_drugbank(drug_target_file, nodes=nodes)

    disease_to_genes = get_diseasome_genes_from_selectedGenes(gene_list_file,disease_name,nodes=nodes)

    # # Calculate proximity values
    print(len(drug_to_targets), len(disease_to_genes))
    print("diseas_genes:", disease_to_genes)

    drug_to_targets_list = []
    for drug, targets in drug_to_targets.items():
        drug_to_targets_list.append([drug,targets])

    print(len(drug_to_targets_list))
    print(len(drug_to_targets))

    network_list = [network] * len(drug_to_targets_list)
    disease_to_genes_list = [disease_to_genes] * len(drug_to_targets_list)

    pool = mp.Pool(processes=num_process)

    start_time = time.time()

    results = pool.imap(calculate_proximity_multiple_multiprocessing_return,
                        zip(network_list, disease_to_genes_list, drug_to_targets_list))

    pool.close()

    results = list(results)
    # results.insert(0,["Drug","Disease","n.source","n.target","d","z"])
    results_df = pd.DataFrame(columns=["Drug", "Disease", "n.source", "n.target", "d", "z"], data=results)
    print(results_df)
    results_df.to_csv(output_file, index=False)
    print("{} times training Runtime: {:2f} Minutes".format(disease_name, ((time.time() - start_time) / 60)))

def main():

    '''
    sys.argv[1] : Hour of DEP (6 or 24)
    sys.argv[2] : The number of proccess for multiprocess
    :return:
    '''

    print("Start network-based proxivity analysis multiprocessing")
    disease_name = "COVID"

    ex_time = sys.argv[1]
    data_type = "{}hr".format(ex_time)



    network_file = "../Data/Human_network/human_protein_interactome.txt"
    gene_list_file = "../Result/Key_proteins/COVID_{}_key_protein.txt".format(data_type)
    drug_target_file = "../Data/Drug_Target/stitch_drugbank_target.v5.1.5_interactions_900_th_400.onlyTarget.Approved.tsv"
    output_file = "../Result/Drug/COVID_drug_network_proximity_{}.csv".format(data_type) # Clinicaltrial TRUE

    num_process = sys.argv[2]
    mp_drug_simulation(disease_name, network_file, gene_list_file, drug_target_file, output_file, num_process)




if __name__ == '__main__':
    main()