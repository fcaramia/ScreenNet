__author__ = 'fcaramia'
import networkx as nx
import csv
from pylab import *
import operator


def print_results(graph, paths, mir_graph, candidates, network_scores, mir_scores, mirs, config, out_file,
                  raw_si_scores, raw_net_scores, raw_mir_scores, interest_genes):
    total_scores = {}
    with open(out_file, 'w') as csvfile:
        out_writer = csv.writer(csvfile, delimiter='\t')
        # print run info
        line = ['#']
        for k in config:
            line += [k + ':' + str(config[k]) + ' ']

        out_writer.writerow(line)

        # print header
        line = ['Gene', 'Total_Score', 'siRNA_Score', 'ScreenNet_Score', 'miRNA_Score', 'siRNA_Raw_Score',
                'ScreenNet_Raw_Score', 'miRNA_Raw_Score', 'ScreenNet_Detail', 'miRNA_Detail']

        for g in interest_genes:
            line += [g]

        out_writer.writerow(line)

        for c in candidates:
            line = [c]
            mir_score = 0
            network_score = 0
            raw_mir_score = 0
            raw_net_score = 0

            if c in mir_scores:
                mir_score = mir_scores[c]
                raw_mir_score = raw_mir_scores[c]

            if c in network_scores:
                network_score = network_scores[c]
                raw_net_score = raw_net_scores[c]

            line += [candidates[c] + mir_score - network_score]
            total_scores[c] = candidates[c] + mir_score - network_score
            line += [str(candidates[c])]
            line += [str(network_score)]
            line += [str(mir_score)]
            line += [str(raw_si_scores[c])]
            line += [str(raw_net_score)]
            line += [str(raw_mir_score)]

            net_detail = ""
            if c in paths:
                for t in paths[c]:
                    p = nx.shortest_path(graph, c, t)
                    for i in range(len(p) - 1):
                        net_detail += p[i] + "->" + p[i + 1] + "(" + str(graph.edge[p[i]][p[i + 1]]['score']) + ") "
                    net_detail += "; "

            line += [net_detail]

            mir_detail = ""
            if c in mir_graph:
                for mir in mir_graph.predecessors(c):
                    if mir in mirs:
                        mir_detail += mir + ":" + str(mirs[mir]) + " "

            line += [mir_detail]

            for g in interest_genes:
                net_detail = ""
                if c in graph:
                    try:
                        p = nx.shortest_path(graph, c, g)
                        for i in range(len(p) - 1):
                            net_detail += p[i] + "->" + p[i + 1] + "(" + str(graph.edge[p[i]][p[i + 1]]['score']) + ") "

                        net_detail += "; "
                        line += [net_detail]
                    except nx.exception.NetworkXNoPath:
                        line += [net_detail]
                        continue






            out_writer.writerow(line)
    return total_scores


def print_stats(norm_si_scores, norm_network_scores, norm_mir_scores, total_scores):
    sorted_si = sorted(norm_si_scores.items(), key=operator.itemgetter(1), reverse=True)
    vals = []
    index = []
    for [x, y] in sorted_si:
        vals.append(y)
        index.append(x)
    plot(vals)
    vals = []
    for x in index:
        if x in norm_network_scores:
            vals.append(norm_network_scores[x])
        else:
            vals.append(0.0)
    plot(vals)
    vals = []
    for x in index:
        if x in total_scores:
            vals.append(total_scores[x])
        else:
            vals.append(0.0)
    plot(vals)
    vals = []
    for x in index:
        if x in norm_mir_scores:
            vals.append(norm_mir_scores[x])
        else:
            vals.append(0.0)
    plot(vals)

    show()
