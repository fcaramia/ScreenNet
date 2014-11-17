__author__ = 'fcaramia'
import networkx as nx
import optparse
from io import *
from utils import *


def have_path(graph, start, end, max_depth, path=[], visited=[], depth=0):
    if depth > max_depth:
        return None
    path = path + [start]
    if start == end:
        return path
    if start not in graph:
        return None
    visited.append(start)
    for node in graph[start]:
        if node in visited:
            continue
        ret = have_path(graph, node, end, max_depth, path, visited, depth + 1)
        if ret: return ret

    return None


def check_graph(candidates, graph, marks, max_depth):
    part = (len(candidates) ** 2) / 100
    count = 0
    for i in range(len(candidates)):
        for j in range(i + 1, len(candidates)):
            count = count + 1

            if count % part == 0:
                print
                count / part, "%"

            c1 = candidates[i]
            c2 = candidates[j]
            # print "Evaluate candidates: ",c1,c2
            if c1 in marks and c2 in marks[c1]:
                continue
            if c2 in marks and c1 in marks[c2]:
                continue

            path = have_path(graph, c1, c2, max_depth, visited=[])
            if path:

                if c1 in marks:
                    marks[c1].append(c2)
                else:
                    marks[c1] = [c2]

            else:
                path = have_path(graph, c2, c1, max_depth, visited=[])

                if path:

                    if c2 in marks:
                        marks[c2].append(c1)
                    else:
                        marks[c2] = [c1]

                    #if path:
                    #	print path
    return marks


def main():
    gene_db_dir = "../gene_db/"
    parser = optparse.OptionParser()
    parser.add_option('-t', '--transFac', action='store_true', dest="trans_fac", help="use TransFac DB")
    parser.add_option('-p', '--phosphoSite', action='store_true', dest="phospho_site", help="use Phosphosite DB")
    parser.add_option('-s', '--string', action='store_true', dest="string_db", help="use String DB")
    parser.add_option('-c', '--candidates', type='string', dest="candidate_file", help="candidates file")
    parser.add_option('-r', '--score', type='int', dest='score', help="minimum evidence score, default: 400",
                      default=400)
    parser.add_option('-d', '--direction', action='store_true', dest='direction',
                      help="evidence of direction for interaction", default=False)
    parser.add_option('-x', '--standardDev', dest='std_dev',
                      help="number of standard deviations to mark candidate for removal, default: 2", default=2.0,
                      type='float')
    parser.add_option('-z', '--maxZScore', dest='zscore', help="max zscore of discarded gene", default=3.0,
                      type='float')
    parser.add_option('-o', '--outfile', dest='out', help="output file", default='out.txt', type='string')
    parser.add_option('-e', '--expand', dest='expand', help="create graph to expand regulation assumption",
                      action="store_true", default=False)
    parser.add_option('-m', '--maxExpand', dest='max_expand', help="max size of path, default: 5", default=5, type=int)

    (options, args) = parser.parse_args()

    if not options.candidate_file:
        parser.error("candidate file is needed")

    #Read Candidates
    candidates = {}
    candidates = read_candidates(options.candidate_file, candidates)

    #Check sign of candidates
    candidates = check_sign_of_candidates(candidates, parser)

    #set graph
    graph = nx.DiGraph()

    #Read DBs
    marks = {}
    if options.trans_fac:
        print
        "Loading TransFac"
        graph = read_reg_db(gene_db_dir + "transfac_interactions.csv", graph, 'transfac', 'transcription factor')

    if options.phospho_site:
        print
        "Loading Phosphosite"
        graph = read_reg_db(gene_db_dir + "kinase_curated_db.csv", graph, "phosphosite", "phosphorylation")

    if options.string_db:
        print
        "Loading String Actions"
        graph = read_string_action_db(gene_db_dir + "actions_curated.tsv", graph, options.score, options.direction)

    #Check for direct regulation
    marks = check_reg_graph(candidates, graph, marks)

    if options.expand:
        print("Warning: expanding increases running time considerably")
        if options.direction is False:
            print("Warning: expanding without directionality could add false positives")
        print("Expanding Search")
        marks = check_graph(candidates.keys(), graph, marks, options.max_expand)

    i = 0
    for r in marks:
        i += len(marks[r])
    print
    i, "interactions marked"

    f = open(options.out, 'w+')
    f.write("Input file: " + options.candidate_file + '\n')
    f.write("DB Interactions marked: " + str(i) + '\n')
    if len(marks) > 0:
        discarded = discard_candidates(candidates, marks, reg_db, string_db, options.std_dev, options.zscore)
        f.write("Discarded Genes:\n")
        for r in discarded:
            f.write(r + ' ' + str(candidates[r]) + "\n")

        f.write("\n\nDiscarded Genes detailed:\n")
        for r in discarded:
            f.write(r + ' ' + str(candidates[r]) + " regulates: \n")
            for g in marks[r]:
                f.write("\t" + g + " " + str(candidates[g]) + '\n')
        print
        "Discarded Genes: ", len(discarded)

    i = 0
    print
    "Direct Regulators:"
    f.write("\n\nGenes not dicarded:\n")
    for c in candidates:
        if c not in discarded:
            f.write(c + ' ' + str(candidates[c]) + '\n')
            i += 1
    print
    i


if __name__ == '__main__':
    main()