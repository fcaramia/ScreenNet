__author__ = 'fcaramia'
import numpy
import math
import networkx as nx
import os


def validate_config(config):
    valid_set = ["yes", "no"]
    bool_keys = ["transFac", "phosphoSite", "string", "mirTarBase", "directed", "filter_by_si"]

    if not all(config[k] in valid_set for k in bool_keys):
        return False

    if not os.path.isdir(config["gene_db_dir"]):
        print("Invalid DB directory")
        return False

    if config["expand"] < 0:
        print("Expand must be bigger than zero")
        return False

    if config["std_dev_val"] < 0:
        print("std_dev_val must be bigger than zero")
        return False

    if config["score_reduce_fun"] not in ["exponential", "lineal"]:
        print("Invalid Score reduce function")
        return False

    if config["network_score_select"] not in ["best", "average"]:
        print("Invalid network score selection")
        return False

    if config["mir_score_select"] not in ["best", "average", "sum"]:
        print("Invalid mir score selection")
        return False

    return True


def check_sign_of_candidates(candidates, parser):
    scores = []
    for c in candidates:
        scores.append(candidates[c])

    sign = numpy.mean(scores)

    # Validate candidates
    for c in candidates:
        if sign > 0 > candidates[c]:
            parser.error("z-scores should have same sign")

        if sign < 0 < candidates[c]:
            parser.error("z-scores should have same sign")

    # Change direction of regulation if negative
    if sign < 0:
        for c in candidates:
            candidates[c] *= -1

    return candidates


# Check if regulation interaction is present in candidates
# Check for repeats and self-regulation
def check_reg_graph(candidates, graph, marks):
    candidates.sort()
    nodes = graph.nodes()
    for c in candidates:
        if c in nodes:
            succ = graph.successors(c)
            for s in succ:
                if s != c and s in candidates:
                    if c in marks and s not in marks[c]:
                        marks[c].append(s)
                    else:
                        marks[c] = [s]
    return marks


def mark_for_scoring(candidates, marks, std_dev_filter):
    scores = []
    ret = {}
    for m in marks:
        for g in marks[m]:
            if candidates[m] < candidates[g]:
                scores.append(math.fabs(candidates[m] - candidates[g]))

    std = numpy.std(scores)

    for m in marks:
        for g in marks[m]:
            if candidates[m] < candidates[g]:
                if math.fabs(candidates[m] - candidates[g]) >= std_dev_filter * std:
                    if m not in ret:
                        ret[m] = [g]
                    else:
                        ret[m].append(g)

    return ret


def check_graph(candidates, graph, marks, max_depth):
    candidates.sort()
    for i in range(len(candidates)):
        for j in range(i + 1, len(candidates)):
            c1 = candidates[i]
            c2 = candidates[j]

            if c1 in marks and c2 in marks[c1]:
                continue
            if c2 in marks and c1 in marks[c2]:
                continue

            if c1 not in graph or c2 not in graph:
                continue

            # Check one way
            try:
                path = nx.shortest_path(graph, c1, c2)
            except nx.exception.NetworkXNoPath:
                continue

            if 0 < len(path) <= max_depth:

                if c1 in marks:
                    marks[c1].append(c2)
                else:
                    marks[c1] = [c2]
            # Check other way
            try:
                path = nx.shortest_path(graph, c2, c1)
            except nx.exception.NetworkXNoPath:
                continue

            if 0 < len(path) <= max_depth:
                if c2 in marks:
                    marks[c2].append(c1)
                else:
                    marks[c2] = [c1]

    return marks


def get_network_scores(paths, graph, reduce_fun, select_score):

    ret = {}

    for s in paths:
        score_sum = 0
        for t in paths[s]:
            p = nx.shortest_path(graph, s, t)
            n = len(p) - 1
            if n <= 0:
                continue
            path_sum = 0
            for i in range(len(p)-1):
                path_sum += graph.edge[p[i]][p[i+1]]['score']

            score = float(path_sum) / float(n)

            if reduce_fun == 'exponential':
                score *= (1.0/2**(n-1))
            else:
                score *= (1.0/2*(n-1))

            score_sum += score
            if s not in ret:
                ret[s] = score
            else:
                if ret[s] < score:
                    ret[s] = score

        if select_score == 'average':
            ret[s] = score_sum/len(paths[s])

    return ret


def get_mir_scores(candidates, mir_graph, mirs, mir_score_select):

    ret = {}

    for c in candidates:
            if c in mir_graph:
                mir_sum = 0
                for mir in mir_graph.predecessors(c):
                    if mir in mirs:
                        mir_sum += mirs[mir]
                        if c not in ret:
                            ret[c] = mirs[mir]
                        else:
                            if ret[c] < mirs[mir]:
                                ret[c] = mirs[mir]

                if mir_score_select is 'sum':
                    ret[c] = mir_sum

                elif mir_score_select is 'average':
                    ret[c] = mir_sum/len(mir_graph.predecessors(c))

    return ret


def normalize_scores(scores, a, b, capping):
    ret = {}
    values = list(scores.values())
    median = numpy.median(values)
    std = numpy.std(values)
    if capping != 0.0:
        for k in scores:
            if(std * capping) < (scores[k] - median):
                scores[k] = median + (std * capping)

    values = list(scores.values())
    min_val = 0
    max_val = max(values)
    # (b-a) (x - min) / max - min

    for k in scores:
        ret[k] = ((1.0-a)*(scores[k] - min_val)) / (max_val - min_val)
        ret[k] = numpy.log2(ret[k] + 1.0) * b
    return ret