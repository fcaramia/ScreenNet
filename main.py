__author__ = 'fcaramia'
import optparse
from readinputs import *
from utils import *
import networkx as nx


def main():

    parser = optparse.OptionParser()

    parser.add_option('-c', '--config', type='string', dest="config_file", help="config file")
    parser.add_option('-s', '--siRNA', type='string', dest="siRNA_file", help="siRNA candidates file")
    parser.add_option('-o', '--outfile', dest='out', help="output file", default='out.txt', type='string')
    parser.add_option('-m', '--miRNA', dest='miRNA_file', help="miRNA file", type='string')
    (options, args) = parser.parse_args()

    if not options.config_file:
        parser.error("config file is needed")

    if not options.siRNA_file:
        parser.error("siRNA file is needed")

    #Read Config
    config = {}
    config = read_config(options.config_file, config)

    #print config
    for k in config:
        print(k, config[k])

    if not validate_config(config):
        parser.error("Invalid config file")

    #Read Candidates
    candidates = {}
    candidates = read_candidates(options.siRNA_file, candidates)

    #Read miRNA
    if options.miRNA_file:
        miRNAs = {}
        miRNAs = read_candidates(options.miRNA_file, miRNAs)

    #Check sign of candidates
    candidates = check_sign_of_candidates(candidates, parser)

    #set graph
    graph = nx.DiGraph()
    mir_graph = nx.DiGraph()

    #Read DBs
    marks = {}

    if config["mirTarBase"] == 'yes':
        print("Loading TransFac")
        mir_graph = read_reg_db(config["gene_db_dir"] + "hsa_MTI.csv", mir_graph, 'mirTarBase', 'miRNA interaction')

    if config["transFac"] == 'yes':
        print("Loading TransFac")
        graph = read_reg_db(config["gene_db_dir"] + "transfac_interactions.csv", graph, 'transfac', 'transcription factor')

    if config["phosphoSite"] == 'yes':
        print("Loading Phosphosite")
        graph = read_reg_db(config["gene_db_dir"] + "kinase_curated_db.csv", graph, "phosphosite", "phosphorylation")

    if config['string'] == 'yes':
        print("Loading String Actions")
        graph = read_string_action_db(config["gene_db_dir"] + "actions_curated.tsv", graph, config['directed'])

    if config["expand"] == 0:
        #Check only for direct regulation
        marks = check_reg_graph(list(candidates.keys()), graph, marks)
    else:
        if config['directed'] == 'no':
            print("Warning: expanding without directionality could add false positives")
        print("Expanding Search")
        print(len(candidates), "siRNA candidates")
        print(len(graph.nodes()), "nodes")
        print(len(graph.edges()), "edges")

        marks = check_graph(list(candidates.keys()), graph, marks, options.max_expand)

    i = 0
    for r in marks:
        i += len(marks[r])

    f = open(options.out, 'w+')
    f.write("Input file: " + options.candidate_file + '\n')

    not_in_db = 0
    nodes = list(graph.nodes())
    for c in candidates:
        if not c in nodes:
            not_in_db += 1

    f.write("Candidate genes with no evidence for selected score: " + str(not_in_db))
    print("Candidate genes with no evidence for selected score: " + str(not_in_db))

    ##DISCARD GENES
    if len(marks) > 0:
        discarded = discard_candidates(candidates, marks, options.std_dev, options.zscore)

        ##PRUNING
        pruned = []
        if options.prune:
            for d in discarded:
                for m in marks:
                    if d in marks[m] and m not in discarded and candidates[m] <= options.zscore:
                        discarded.append(m)
                        pruned.append([m, d])

        f.write("\n\nDiscarded Genes: "+str(len(discarded)) + "\n")
        print("Discarded Genes: "+str(len(discarded)) + "\n")
        print("Pruned Genes: "+str(len(pruned)) + "\n")

        for r in discarded:
            f.write(r + ' ' + str(candidates[r]) + "\n")

        f.write("\n\nDiscarded Genes detailed (empty for pruned nodes):\n")
        for r in discarded:
            f.write(r + ' ' + str(candidates[r]) + " regulates: \n")
            for g in marks[r]:
                if candidates[r] < candidates[g]:
                    f.write("\t" + g + " " + str(candidates[g]) + " " + str(nx.shortest_path(graph, r, g)) + '\n')

        f.write("\n\nPruned detailed:\n")
        for [p, d] in pruned:
            f.write(p + ' pruned by ' + d + "\n")

    i = 0
    f.write("\n\nGenes not dicarded:\n")
    for c in candidates:
        if c not in discarded:
            f.write(c + ' ' + str(candidates[c]) + '\n')
            i += 1
    print("Genes not discarded "+str(i))

if __name__ == '__main__':
    main()