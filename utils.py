__author__ = 'fcaramia'
import numpy
import math
import networkx as nx


def check_sign_of_candidates(candidates, parser):
    scores = []
    for c in candidates:
        scores.append(candidates[c])

    sign = numpy.mean(scores)

    #Validate candidates
    for c in candidates:
        if sign > 0 > candidates[c]:
            parser.error("z-scores should have same sign")

        if sign < 0 < candidates[c]:
            parser.error("z-scores should have same sign")

    ##Change direction of regulation if negative
    if sign < 0:
        for c in candidates:
            candidates[c] *= -1

    return candidates


# Check if regulation interaction is present in candidates
#Check for repeats and self-regulation
def check_reg_graph(candidates, graph, marks):
    for c in candidates:
        if c in graph:
            for s in graph.successors(c):
                if s != c and s in candidates:
                    if c in marks and s not in marks[c]:
                        marks[c].append(s)
                    else:
                        marks[c] = [s]
    return marks


def discard_candidates(candidates, marks, std_dev_filter, zscore):
    scores = []
    discarded = []
    for m in marks:
        for g in marks[m]:
            if candidates[m] > candidates[g]:
                scores.append(math.fabs(candidates[m] - candidates[g]))

    std = numpy.std(scores)

    for m in marks:
        for g in marks[m]:
            if candidates[m] < candidates[g]:
                if math.fabs(candidates[m] - candidates[g]) >= std_dev_filter * std:
                    if m not in discarded and candidates[m] <= zscore:
                        discarded.append(m)

    return discarded


def check_graph(candidates, graph, marks, max_depth):
    part = (len(candidates) ** 2) / 100
    count = 0
    for i in range(len(candidates)):
        for j in range(i + 1, len(candidates)):
            count += 1

            if count % part == 0:
                print
                (count / part)*2, "%"

            c1 = candidates[i]
            c2 = candidates[j]
            # print "Evaluate candidates: ",c1,c2
            if c1 in marks and c2 in marks[c1]:
                continue
            if c2 in marks and c1 in marks[c2]:
                continue
            if not c1 in graph or c2 not in graph:
                continue

            try:
                path = nx.shortest_path_length(graph, c1, c2)
            except nx.exception.NetworkXNoPath:
                continue

            if 0 < path < max_depth:

                if c1 in marks:
                    marks[c1].append(c2)
                else:
                    marks[c1] = [c2]

            else:
                try:
                    path = nx.shortest_path_length(graph, c2, c1)
                except nx.exception.NetworkXNoPath:
                    continue

                if 0 < path < max_depth:

                    if c2 in marks:
                        marks[c2].append(c1)
                    else:
                        marks[c2] = [c1]

                    #if path:
                    #	print path
    return marks