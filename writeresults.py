__author__ = 'fcaramia'
import networkx as nx
import csv


def print_results(graph, paths, mir_graph, candidates, network_scores, mir_scores, mirs, config, out_file):
    with open(out_file, 'w') as csvfile:
        out_writer = csv.writer(csvfile, delimiter='\t')
        # print run info
        line = ['#']
        for k in config:
            line += [k + ':' + str(config[k]) + ' ']

        out_writer.writerow(line)

        # print header
        line = ['Gene', 'Total_Score', 'siRNA_Score', 'ScreenNet_score', 'miRNA_Score', 'ScreenNet_Detail',
                'miRNA_detail']

        out_writer.writerow(line)

        for c in candidates:
            line = [c]
            mir_score = 0
            network_score = 0

            if c in mir_scores:
                mir_score = mir_scores[c]

            if c in network_scores:
                network_score = network_scores[c]

            line += [config['siRNA'] * candidates[c] + config['miRNA'] * mir_score - config['network'] * network_score]

            line += [str(candidates[c])]
            line += [str(network_score)]
            line += [str(mir_score)]

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

            out_writer.writerow(line)