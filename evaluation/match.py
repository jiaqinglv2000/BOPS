# -*- coding: utf-8 -*-

"""\
%prog [options] reference_complexes predicted_complexes

Calculates matching scores between a set of reference and predicted complexes.
The input files must contain the reference and the predicted complexes, one
complex per line.
"""

from __future__ import division
from textwrap import dedent
import optparse
import sys
import xlwt
import pandas as pd
from mwmatching import maxWeightMatching

###########################################################################

def canonical_protein_name(name):
    return name.strip().upper()

def is_numeric(x):
    try:
        float(x)
        return True
    except:
        return False

def matching_score(set1, set2):
    return len(set1.intersection(set2))**2 / (float(len(set1)) * len(set2))

def fraction_matched(reference, predicted, score_threshold=0.25):
    recall=recall_matched(reference, predicted,score_threshold);
    precision=predicted_matched(reference, predicted, score_threshold);
    if recall*precision==0 :
        return 0
    return 2*recall*precision/(recall+precision)

def accuracy(reference, predicted):
    return (clusteringwise_sensitivity(reference, predicted) * positive_predictive_value(reference, predicted)) ** 0.5

def clusteringwise_sensitivity(reference, predicted):
    if len(predicted)==0 :
        return 0
    num, den = 0., 0.
    for complex in reference:
        den += len(complex)
        num += max(len(complex.intersection(cluster)) for cluster in predicted)
    if den == 0.:
        return 0.
    return num / den

def recall_matched(reference, predicted, score_threshold=0.25):
    result=0
    for id1, c1 in enumerate(reference):
        for id2, c2 in enumerate(predicted):
            score = matching_score(c1, c2)
            if score > score_threshold:
                result += 1
                break
    if len(reference)==0:
        return 0
    return result/len(reference)

def predicted_matched(reference, predicted, score_threshold=0.25):
    result = 0
    for id1, c1 in enumerate(predicted):
        for id2, c2 in enumerate(reference):
            score = matching_score(c1, c2)
            if score > score_threshold:
                result += 1
                break
    if len(predicted)==0:
        return 0
    return result/len(predicted)

def positive_predictive_value(reference, predicted):
    num, den = 0., 0.
    for cluster in predicted:
        isects = [len(cluster.intersection(complex)) for complex in reference]
        isects.append(0.)
        num += max(isects)
        den += sum(isects)
    if den == 0.:
        return 0.
    return num / den

def maximum_matching_ratio(reference, predicted, score_threshold=0.25):
    scores = {}

    n = len(reference)
    for id1, c1 in enumerate(reference):  #遍历参考集的所有元素
        for id2, c2 in enumerate(predicted):  #遍历预测集的所有元素
            score = matching_score(c1, c2)
            if score <= score_threshold:
                continue

            scores[id1, id2+n] = score

    input = [(v1, v2, w) for (v1, v2), w in scores.items()]  #items
    mates = maxWeightMatching(input)
    score = sum(scores[i, mate] for i, mate in enumerate(mates) if i < mate)
    return score / n

###########################################################################

class MatchApplication(object):
    def __init__(self, run_print=True,):
        self.__run_print=True
        self.measures = dict(
                ACC=accuracy,
                Sn=clusteringwise_sensitivity,
                PPV=positive_predictive_value,
                Fscore=fraction_matched,
                Recall=recall_matched,
                Predicted=predicted_matched,
                MMR=maximum_matching_ratio,
        )
        self.parser = self.create_parser()

    def create_parser(self):
        parser = optparse.OptionParser(usage=dedent(sys.modules[__name__].__doc__).strip())

        parser.add_option("-m", "--measure", action="append", dest="measures", default=[],
                metavar="MEASURE", help="calculate the quality measure given by MEASURE. "
                "Possible values are: %s. May be given multiple times." %
                ", ".join(sorted(self.measures.keys())))
        parser.add_option("-n", "--network", metavar="FILE", dest="network",
                help="read the PPI network from FILE and assume that only these complexes "
                "were known to the clustering algorithm")
        parser.add_option("-q", "--quiet", action="store_true", dest="quiet",
                default=False, help="be quiet")

        return parser

    def log(self, msg):
        if self.options.quiet:
            return
        print >>sys.stderr, msg

    def read_complexes(self, fname, known_proteins=None, strictness=0.5,
                       min_size=3, max_size=100,min_num=0,max_num=10):
        result = []
        for line in open(fname):
            ps = set(canonical_protein_name(x) for x in line.strip().split() if x)
            if known_proteins is not None:
                isect = ps.intersection(known_proteins)
                ps = list(x for x in ps if x[0] != 'C')
                if len(isect) < max(min_size, len(ps) * strictness):
                    continue
                if len(isect) > max_size:
                    continue
                if len(ps) < min_num:
                    continue
                if len(ps) > max_num:
                    continue
                ps = isect
            result.append(set(ps))
        to_delete = set()
        for idx, cluster in enumerate(result):
            for idx2, cluster2 in enumerate(result):
                if idx == idx2 or idx2 in to_delete:
                    continue
                if cluster == cluster2:
                    to_delete.add(idx2)

        result = [r for i, r in enumerate(result) if i not in to_delete]
        return result
    def read_network(self, fname):
        known_proteins = set()
        for line in open(fname):
            parts = [canonical_protein_name(part) for part in line.strip().split()
                    if not is_numeric(part)]
            known_proteins.update(parts)
        return known_proteins
    def run(self):
        self.options, self.args = self.parser.parse_args()
        if len(self.args) != 2:
            print self.options,self.args
            self.parser.print_help()
            return 1

        if not self.options.measures:
            self.options.measures = sorted(self.measures.keys())

        if self.options.network:
            known_proteins = self.read_network(self.options.network)
            self.log("%d known proteins found in network" % len(known_proteins))
        else:
            known_proteins = None

        reference_complexes = self.read_complexes(self.args[0], known_proteins)
        predicted_complexes = self.read_complexes(self.args[1], known_proteins)
        self.log("%d reference complexes, %d predicted complexes" %
                (len(reference_complexes), len(predicted_complexes)))
        ppiresult={'referP':len(reference_complexes),'preP':len(predicted_complexes),}
        for measure in self.options.measures:
            if measure not in self.measures:
                self.log("Ignoring unknown measure: %s" % measure)
                continue
            result = self.measures[measure](reference_complexes, predicted_complexes)
            if self.__run_print==True:
                print "%s = %.4f" % (measure, result)
            ppiresult.update({measure:result})
        if self.__run_print==True:
            return 0
        else:
            return ppiresult

def main():
    return MatchApplication().run()

if __name__ == "__main__":
    sys.exit(main())