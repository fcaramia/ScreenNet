__author__ = 'fcaramia'
import optparse
from readinputs import *
from utils import *
from writeresults import *
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

    # Read Config
    config = {}
    config = read_config(options.config_file, config)

    # Print config
    for k in config:
        print(k, config[k])

    if not validate_config(config):
        parser.error("Invalid config file")

    # Read Candidates
    candidates = {}
    candidates = read_candidates(options.siRNA_file, candidates)

    miRNAs = {}
    # Read miRNA
    if options.miRNA_file:
        miRNAs = read_candidates(options.miRNA_file, miRNAs)
        # Check sign of miRNAs
        miRNAs = check_sign_of_candidates(miRNAs, parser)

    # Check sign of candidates
    candidates = check_sign_of_candidates(candidates, parser)

    # set graph
    graph = nx.DiGraph()
    mir_graph = nx.DiGraph()

    # Read DBs
    marks = {}

    if config["mirTarBase"] == 'yes':
        print("Loading mirTarBase")
        mir_graph = read_reg_db(config["gene_db_dir"] + "hsa_MTI.csv", mir_graph, 'mirTarBase', 'miRNA interaction')

    if config["transFac"] == 'yes':
        print("Loading TransFac")
        graph = read_reg_db(config["gene_db_dir"] + "transfac_interactions.csv", graph, 'transfac',
                            "transcription factor", 1000)

    if config["phosphoSite"] == 'yes':
        print("Loading Phosphosite")
        graph = read_reg_db(config["gene_db_dir"] + "kinase_curated_db.csv", graph, "phosphosite",
                            "phosphorylation", 1000)

    if config['string'] == 'yes':
        print("Loading String Actions")
        graph = read_string_action_db(config["gene_db_dir"] + "actions_curated.tsv", graph, config['directed'],
                                      config["min_score"])

    # Mark genes that are connected to each other
    if config["expand"] == 0:
        # Check only for direct regulation
        marks = check_reg_graph(list(candidates.keys()), graph, marks)
    else:
        if config['directed'] == 'no':
            print("Warning: creating a non directed network could add false positives")
        print("Expanding Search")
        print(len(candidates), "siRNA candidates")
        print(len(graph.nodes()), "nodes")
        print(len(graph.edges()), "edges")

        marks = check_graph(list(candidates.keys()), graph, marks, config['expand'])

    i = 0
    for r in marks:
        i += len(marks[r])

    not_in_db = 0
    nodes = list(graph.nodes())
    for c in candidates:
        if c not in nodes:
            not_in_db += 1

    print("Candidate genes with no evidence for selected score: " + str(not_in_db))

    # Generate Network score for marked genes
    norm_network_scores = {}
    paths_for_scoring = {}
    if len(marks) > 0:

        paths_for_scoring = mark_for_scoring(candidates, marks, config['std_dev_val'])
        # Generate network scores
        network_scores = get_network_scores(paths_for_scoring, graph, config['score_reduce_fun'],
                                            config['network_score_select'])

        norm_network_scores = normalize_scores(network_scores)

    # Generate miRNA scores
    norm_mir_scores = {}
    if len(mir_graph) > 0 and len(miRNAs) > 0:
        # Search predecessors of genes and compute mir scores
        mir_scores = get_mir_scores(candidates, mir_graph, miRNAs, config['mir_score_select'])
        norm_mir_scores = normalize_scores(mir_scores)

    norm_si_scores = normalize_scores(candidates)

    print_results(graph, paths_for_scoring, mir_graph, norm_si_scores, norm_network_scores, norm_mir_scores, miRNAs, config, options.out)


if __name__ == '__main__':
    main()