__author__ = 'fcaramia'
import numpy
import math


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