# Created by woochanghwang at 15/04/2020

import networkx as nx
import collections
import matplotlib.pyplot as plt


def degree_histogram(G, type):

    # G = nx.gnp_random_graph(100, 0.02)

    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    # print "Degree sequence", degree_sequence
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())

    fig, ax = plt.subplots()
    plt.bar(deg, cnt, width=0.80, color='b')

    plt.title("Degree Histogram")
    plt.ylabel("Count")
    plt.xlabel("Degree")
    ax.set_xticks([d + 0.4 for d in deg])
    ax.set_xticklabels(deg)

    # draw graph in inset
    plt.axes([0.4, 0.4, 0.5, 0.5])
    Gcc = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])
    pos = nx.spring_layout(G)
    plt.axis('off')
    nx.draw_networkx_nodes(G, pos, node_size=20)
    nx.draw_networkx_edges(G, pos, alpha=0.4)

    plt.savefig("/home/wch23/Project/LifeArc/COVID-19/result/network/{}_degree_histogram.pdf".format(type))

    plt.show()

def main():
    background_addr = "/home/wch23/Project/LifeArc/General/src_drug/Data/human_protein_interactome.txt"
    corona_addr = "/home/wch23/Project/LifeArc/COVID-19/result/network/COVID_core_interactions_shortest_both_with_PPI_DEG_n5.tsv"

    target_nodes = ["EGFR", "TP53", "MYC"]
    target_entrez =["1956", "7157", "4609"]

    with open(background_addr) as back_f:
        back_edges = [x.strip().split('\t') for x in back_f.readlines()]

    with open(corona_addr) as corona_f:
        corona_edges = [x.strip().split('\t') for x in corona_f.readlines()]


    back_edges = [[x[0],x[2]] for x in back_edges[1:]]
    corona_edges = corona_edges[1:]

    back_G = nx.Graph(back_edges)
    corona_G = nx.Graph(corona_edges)

    for target in target_entrez:
        print(target, back_G.degree[target])

    # back_degree = sorted([d for n, d in back_G.degree()], reverse=True)
    # print("max: ", back_degree[:10])
    # print("min: ", back_degree[-2:])
    # print("average:", sum(back_degree)/len(back_degree))

    for target in target_nodes:
        print(target, corona_G.degree[target])

    back_degree = sorted([d for n, d in corona_G.degree()], reverse=True)
    print("max: ", back_degree[:10])
    print("min: ", back_degree[-2:])
    print("average:", sum(back_degree)/len(back_degree))
    # degree_histogram(back_G, "Background")
    # degree_histogram(corona_G, "Corona")






if __name__ == '__main__':
    main()