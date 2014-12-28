__author__ = 'fcaramia'
import networkx as nx
import csv


def print_results(graph, candidates, network_scores, mir_scores, config, out_file):
    with open(out_file, 'wb') as csvfile:
        out_writer = csv.writer(csvfile, delimiter='\t')
        # print run info
        line = '#'
        for k in config:
            line += k + ':' + config[k] + ' '
        # print header
        out_writer.writerow(['Gene', 'Total_Score', 'siRNA_Score', 'ScreenNet_score', 'miRNA_Score', 'ScreenNet_Detail',
                             'miRNA_detail'])

        for c in candidates:
            line = [c]

            out_writer.writerow(line)