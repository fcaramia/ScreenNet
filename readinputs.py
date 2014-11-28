__author__ = 'fcaramia'


#read String db by certain filters
def read_string_action_db(db_file, graph, min_score, direction, mode=None, action=None):
    db_obj = open(db_file, 'r+')
    #Skip header
    next(db_obj)
    for line in db_obj:
        [a, b, mode, action, a_is_acting, score, source, source2] = line.rstrip('\n').split('\t')
        #Apply filters
        if min_score > int(score):
            continue
        if a == ' ' or b == ' ':
            continue
        if a == '' or b == '':
            continue
        if a == "NA" or b == "NA":
            continue
        if direction and a_is_acting == '0':
            continue

        if direction:
            graph.add_edge(a, b, source=source, interaction=mode, score=score)
        else:
            if a_is_acting == "0":
                graph.add_edge(a, b, source=source, interaction=mode, score=score)
                graph.add_edge(b, a, source=source, interaction=mode, score=score)
            else:
                graph.add_edge(a, b, source=source, interaction=mode, score=score)

    return graph


#Same function for transfac and phosphosite
def read_reg_db(db_file, graph, src, inter):
    db_obj = open(db_file, 'r+')
    #Skip header
    next(db_obj)
    for line in db_obj:
        [reg, gene] = line.rstrip('\n').split(',')
        if reg == ' ' or gene == ' ':
            continue
        if gene == '' or reg == '':
            continue
        if gene == "NA" or reg == "NA":
            continue

        graph.add_edge(reg, gene, source=src, interaction=inter)

    return graph


def read_candidates(candidate_file, candidates):
    candidates_obj = open(candidate_file)
    #Skip header
    next(candidates_obj)
    for line in candidates_obj:
        [gene, score] = line.rstrip('\n').split(',')

        candidates[gene] = float(score)

    return candidates
