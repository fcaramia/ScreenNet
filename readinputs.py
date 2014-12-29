__author__ = 'fcaramia'
import configparser


# read config
def read_config(file, config):
    config_parser = configparser.ConfigParser()
    config_parser.read_file(open(file))
    config["transFac"] = config_parser.get("databases", "transFac")
    config["phosphoSite"] = config_parser.get("databases", "phosphoSite")
    config["string"] = config_parser.get("databases", "string")
    config["mirTarBase"] = config_parser.get("databases", "mirTarBase")
    config["gene_db_dir"] = config_parser.get("databases", "gene_db_dir")
    config["directed"] = config_parser.get("network", "directed")
    config["expand"] = config_parser.getint("network", "expand")
    config["siRNA"] = config_parser.getfloat("weights", "siRNA")
    config["network"] = config_parser.getfloat("weights", "network")
    config["miRNA"] = config_parser.getfloat("weights", "miRNA")
    config["std_dev_val"] = config_parser.getfloat("scores", "std_dev_val")
    config["score_reduce_fun"] = config_parser.get("network", "score_reduce_fun")
    config["network_score_select"] = config_parser.get("network", "network_score_select")
    config["mir_score_select"] = config_parser.get("miRNA", "mir_score_select")
    config['min_score'] = config_parser.get('network', 'min_score')
    config['filter_by_si'] = config_parser.get('network', 'filter_by_si')

    return config


# read String db by certain filters
def read_string_action_db(db_file, graph, direction, score_val=0, mode=None, action=None):
    db_obj = open(db_file, 'r+')
    # Skip header
    next(db_obj)
    for line in db_obj:
        [a, b, mode, action, a_is_acting, score, source, source2] = line.rstrip('\n').split('\t')
        # Apply filters
        if a == ' ' or b == ' ':
            continue
        if a == '' or b == '':
            continue
        if a == "NA" or b == "NA":
            continue
        if direction == 'yes' and a_is_acting == '0':
            continue
        if score < score_val:
            continue

        if direction == 'yes':
            graph.add_edge(a, b, source=source, interaction=mode, score=float(score))
        else:
            if a_is_acting == "0":
                graph.add_edge(a, b, source=source, interaction=mode, score=float(score))
                graph.add_edge(b, a, source=source, interaction=mode, score=float(score))
            else:
                graph.add_edge(a, b, source=source, interaction=mode, score=float(score))

    return graph


# Same function for transfac and phosphosite
def read_reg_db(db_file, graph, src, inter, score=0):
    db_obj = open(db_file, 'r+')
    # Skip header
    next(db_obj)
    for line in db_obj:
        [reg, gene] = line.rstrip('\n').split(',')
        if reg == ' ' or gene == ' ':
            continue
        if gene == '' or reg == '':
            continue
        if gene == "NA" or reg == "NA":
            continue

        graph.add_edge(reg, gene, source=src, interaction=inter, score=float(score))

    return graph


def read_candidates(file, candidates):
    candidates_obj = open(file)
    # Skip header
    next(candidates_obj)
    for line in candidates_obj:
        [val, score] = line.rstrip('\n').split(',')

        candidates[val] = float(score)

    return candidates
