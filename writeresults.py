__author__ = 'fcaramia'
import networkx as nx
import csv
from pylab import *
import os


def print_results(out_file, graph, paths, mir_graph, total_scores, si_norm, network_norm_neg, network_norm_pos, mir_norm, config,
                  raw_si_scores, raw_net_scores_neg, raw_net_scores_pos, raw_mir_scores, miRNAs, interactome):

    interest_genes = config['interest_genes']
    sorted_total = sorted(total_scores.items(), key=lambda x: (x[1], x[0]), reverse=True)
    sorted_si = sorted(si_norm.items(), key=lambda x: (x[1], x[0]), reverse=True)

    indexes_total_scores = {}
    indexes_si_scores = {}
    for i in range(len(sorted_total)):
        indexes_total_scores[sorted_total[i][0]] = i
        indexes_si_scores[sorted_si[i][0]] = i

    with open(out_file, 'w') as csvfile:
        out_writer = csv.writer(csvfile, delimiter='\t')
        # print run info
        line = ['#']
        for k in config:
            line += [k + ':' + str(config[k]) + ' ']

        out_writer.writerow(line)

        for k in config:
            line += [k + ':' + str(config[k]) + ' ']

        # print header
        line = ['New_Rank', 'Original_Rank', 'Ranking_difference', 'Gene', 'Total_Score', 'siRNA_Score',
                'ScreenNet_Score_Neg', 'ScreenNet_Score_Pos', 'miRNA_Score', 'siRNA_Raw_Score',
                'ScreenNet_Raw_Score_Neg', 'ScreenNet_Raw_Score_Pos', 'miRNA_Raw_Score',
                'interactome_info', 'ScreenNet_Detail', 'miRNA_Detail']
        # interest genes in header
        for g in interest_genes:
            line += [g]

        out_writer.writerow(line)

        for c in total_scores:
            line = [str(indexes_total_scores[c]), str(indexes_si_scores[c]),
                    str(indexes_si_scores[c]-indexes_total_scores[c]), c]

            line += [total_scores[c]]
            line += [str(si_norm[c])]

            network_score_neg = 0
            raw_net_score_neg = 0
            network_score_pos = 0
            raw_net_score_pos = 0
            if c in network_norm_neg:
                network_score_neg = network_norm_neg[c]
                raw_net_score_neg = raw_net_scores_neg[c]

            if c in network_norm_pos:
                network_score_pos = network_norm_pos[c]
                raw_net_score_pos = raw_net_scores_pos[c]

            mir_score = ''
            raw_mir_score = ''
            if c in mir_norm:
                mir_score = mir_norm[c]
                raw_mir_score = raw_mir_scores[c]

            inter = ''
            if c in interactome:
                inter = interactome[c]

            line += [str(network_score_neg)]
            line += [str(network_score_pos)]
            line += [str(mir_score)]
            line += [str(raw_si_scores[c])]
            line += [str(raw_net_score_neg)]
            line += [str(raw_net_score_pos)]
            line += [str(raw_mir_score)]
            line += [str(inter)]

            net_detail = ""
            if c in paths:
                for t in paths[c]:
                    p = nx.shortest_path(graph, c, t)
                    for i in range(len(p) - 1):
                        net_detail += p[i] + "->" + p[i + 1] + "(" + str(graph.edge[p[i]][p[i + 1]]['score']) + ', ' + \
                            str(graph.edge[p[i]][p[i + 1]]['interaction']) + ") "
                    net_detail += "; "

            line += [net_detail]

            mir_detail = ""
            if c in mir_graph:
                for mir in mir_graph.predecessors(c):
                    if mir in miRNAs:
                        mir_detail += mir + ":" + str(miRNAs[mir]) + " "

            line += [mir_detail]

            # Path to interest genes
            for g in interest_genes:
                net_detail = ""
                if c in graph:
                    try:
                        p = nx.shortest_path(graph, c, g)
                        for i in range(len(p) - 1):
                            net_detail += p[i] + "->" + p[i + 1] + "(" + str(graph.edge[p[i]][p[i + 1]]['score']) + ', ' + \
                                str(graph.edge[p[i]][p[i + 1]]['interaction']) + ") "
                        net_detail += "; "
                        line += [net_detail]
                    except nx.exception.NetworkXNoPath:
                        line += [net_detail]
                        continue

            out_writer.writerow(line)
        # print re-ranking
        ranks = [50, 100, 150, 200, 500, 1000]
        # line = ['# Re-Rankings']
        sorted_si = sorted(raw_si_scores.items(), key=lambda x: (x[1], x[0]), reverse=True)
        total_scores_sorted = sorted(total_scores.items(), key=lambda x: (x[1], x[0]), reverse=True)
        for r in ranks:
            re_ranked = 0
            d = dict(total_scores_sorted[0:r])

            for i in sorted_si[0:r]:
                if i[0] in d:
                    re_ranked += 1
            # line += [str(r)+":"+str(re_ranked)+"; "]
            print(str(r)+":"+str(re_ranked)+"; ")
        # out_writer.writerow(line)


def print_stats(norm_si_scores, norm_network_scores_neg, norm_network_scores_pos, norm_mir_scores, total_scores, out):
    sorted_si = sorted(norm_si_scores.items(), key=lambda w: (w[1], w[0]), reverse=True)
    vals = []
    index = []
    for [x, y] in sorted_si:
        vals.append(y)
        index.append(x)
    plot(vals)
    vals = []
    for x in index:
        if x in norm_network_scores_neg:
            vals.append(norm_network_scores_neg[x])
        else:
            vals.append(0.0)
    plot(vals)
    vals = []
    for x in index:
        if x in norm_network_scores_pos:
            vals.append(norm_network_scores_pos[x])
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

    savefig(os.path.splitext(out)[0]+'.png')
