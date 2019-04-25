#!/usr/env python2.7

__author__ = 'newuser'

import csv
from Bio import Align
from Bio import AlignIO
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import pyfaidx
import re
from os.path import basename
from sys import argv
from collections import Counter
from dendropy import popgenstat
import dendropy
from dendropy.interop import raxml
import shutil
import matplotlib.pyplot as plt
from ete2 import Tree
script = argv[0]
skip_locus = False
tax_table = {"pol_bear": ["pob", "uma", "pb", "bgi"],
             "bro_bear": ["br", "br", "abc", "sw", "uar"],
             "amb_bear": ["amb", "american", "^7_", "uam"],
             "asb_bear": ["asb", "uth", "^16_", "asiatic", "asian"],
             "sun_bear": ["sun", "hma", "^13_fr", "^15_fr"],
             "slo_bear": ["slo", "mur", "^14_fr"],
             "spc_bear": ["spec", "tor", "^17_", "^18_"]}

result_table = []


def load_matrix():
    matrix = []
    with open(options.matrix) as matrix_infile:
        print "[{}] Info: Loading matrix {} with {} lines".format(script, options.matrix,
                                                                  len(matrix_infile.readlines()))
        matrix_infile.seek(0)  # put cursor in file back to zero
        matrix_csv = csv.reader(matrix_infile, delimiter="\t")
        for marker_file in matrix_csv:
            seqname, start, stop, taxa = [e.strip() for e in marker_file[0:4]]
            taxa = [associate_taxon(t) for t in taxa.split(",")]
            marker = {"seqname": seqname, "start": int(start), "stop": int(stop), "taxa": taxa}
            matrix.append(marker)
    return matrix


def associate_taxon(sample_name):
    m = []
    for k, v in tax_table.iteritems():
        re_match = [re.search(e, sample_name.lower()) for e in v]
        m.extend([i for i in re_match if i])
    if m:
        taxon = [k for k, v in tax_table.items() if m[0].group() in [i.strip("^") for i in v]]
        if len(set(taxon)) > 1:
            print "[{}] Error: Ambiguos sample name {}".format(script, sample_name)
            return 0
        elif len(taxon) == 1:
            taxon = [k for k, v in tax_table.items() if m[0].group() in [i.strip("^") for i in v]][0]
            print "[{}] Info: identified {} as {}".format(script, sample_name.lower(), taxon)
            return taxon
        else:
            print "Something else happend {}".format(sample_name)
    else:
        print "[{}] Error: Sample name {sample} could not be associated to known taxa.".format(script,
                                                                                               sample=sample_name)
    return False


def concat_sequences(sequences):
    concatenated = Seq("")
    for s in sequences:
        concatenated += s

    return concatenated


def identify_ST_SNPs(column):
    ST_congruent_splits = {"pol_bro": ["pol_bear", "bro_bear"],
                           "american": ["pol_bear", "bro_bear", "amb_bear"],
                           "asian": ["asb_bear", "sun_bear", "slo_bear"],
                           "sun_slo": ["sun_bear", "slo_bear"],
                           "ursinae": ["pol_bear", "bro_bear", "amb_bear", "asb_bear", "sun_bear", "slo_bear"]
                           }

    for splitname, splittaxa in ST_congruent_splits.iteritems():
        rows_state1 = [i for i in xrange(len(column)) if column[i].id in splittaxa]
        rows_state2 = [i for i in xrange(len(column)) if column[i].id not in splittaxa]

        state1 = set(concat_sequences(column[row, :].seq for row in rows_state1))
        state2 = set(concat_sequences(column[row, :].seq for row in rows_state2))

        if len(rows_state1) == len(splittaxa) and len(state1) == 1 and list(state1)[0] not in state2:
            # check that state is unambiguos and is different than that absent
            return splitname
        else:
            continue
def get_popgen_stats(msa_fname, ins_present):

    seqs = dendropy.DnaCharacterMatrix.get_from_path(msa_fname,"fasta")

    try:
        pi = popgenstat.nucleotide_diversity(seqs)
        seg_sites = popgenstat.num_segregating_sites(seqs)
        avg_pd = popgenstat.average_number_of_pairwise_differences(seqs)
    except ZeroDivisionError:
        pi = "NA"
        seg_sites = "NA"
        avg_pd = "NA"
    stats = {   "seg_sites": seg_sites,
                "nuc_div":  pi,
                "avg_pd":  avg_pd}


    return stats
def process_alignment(msa, ins_present):
    SNPcounter = {"allSNP": 0,
                  "infSNP": 0,
                  "othSNP": 0,
                  "insSNP": 0,
                  "noSNP": 0,
                  "sptreeSNP": 0,
                  "baseSNP": 0,
                  "insDiv":0,
                  "absDiv": 0 }
    rows_ins_present = [i for i in xrange(len(msa)) if msa[i].id in ins_present]
    rows_ins_absent = set(xrange(len(msa))) - set(rows_ins_present)
    nuc_ins_present = ""
    nuc_ins_absent = ""

    for i in xrange(0, len(msa[1])):
        column = msa[:, i:i + 1]
        badchar = False

        for nuc in column:
            if str(nuc.seq).upper() in "YRKMWSBDHVN-":  # really exclude het sites?
                badchar = True
                break

        if not badchar:
            # TODO: Check whether observed SNPs are in taxa with MEI present
            # or with MEI absent
            if 2 <= len(set(column[:, 0])) < len(column[:, 0]):  # is site phyologentically informative
                SNPcounter["allSNP"] += 1
                if len(rows_ins_present) == 0:
                    continue
                if len([e for e in Counter(column[:, 0]).values() if e >= 2]) >= 2:
                    SNPcounter["infSNP"] += 1  # informative site

                    spTree = identify_ST_SNPs(column)
                    if spTree and spTree != "ursinae":
                        print "Species Tree SNP: %s" % spTree
                        SNPcounter["sptreeSNP"] += 1
                    elif spTree and spTree == "ursinae":
                        SNPcounter["baseSNP"] += 1

                    state_ins_present = set(concat_sequences(
                        column[row, :].seq for row in rows_ins_present))  # collect states of ins_present samples
                    state_ins_absent = set(concat_sequences(
                        column[row, :].seq for row in rows_ins_absent))  # states for ins_absent samples

                    # if nuc_ins_present and nuc_ins_present == "C":
                    #     if list(state_ins_present)[0] == "G":
                    #         SNPcounter["insCpG"] += 1
                    #
                    # if nuc_ins_absent and nuc_ins_absent == "C":
                    #     if list(state_ins_absent)[0] == "G":
                    #         SNPcounter["absCpG"] += 1
                    # nuc_ins_present = list(state_ins_present)[0]
                    # nuc_ins_absent = list(state_ins_absent)[0]
                    if len(state_ins_present) == 1 and list(state_ins_present)[0] not in state_ins_absent:
                        # check that state is unambiguos and is different than that absent
                        SNPcounter["insSNP"] += 1
                        # SNPcounter["insDiv"] += 1.0/(len(rows_ins_present))
                        pass
                    # Record genetic diversity within sequence with MEI present
                    # require this site to be different then sequence without MEI
                    if len(state_ins_present) >= 2:
                        SNPcounter["insDiv"] += 1.0/(len(rows_ins_present))
                    # Record genetic diversity within sequence without MEI
                    # require this site to be different then sequence with MEI
                    if len(state_ins_absent) >= 2:
                        SNPcounter["absDiv"] += 1.0/(len(rows_ins_present))
                else:
                    SNPcounter["othSNP"] += 1  # uninformative site
            else:
                SNPcounter["noSNP"] += 1

    return SNPcounter


def check_monophyly(msa_fname, ins_present):
    print msa_fname
    d = dendropy.DnaCharacterMatrix.get_from_path(msa_fname, "fasta")
    for taxon, map in d.items():
        taxon.label =  taxon.label.split(" ")[0]
    # initiate raxml (make sure that raxmlHPC is executable)
    rx = raxml.RaxmlRunner()
    # calculate tree
    tree  = rx.estimate_tree(d, ['-m', 'GTRCAT', '-V','-N', '100'])
    #### the following is only for drawing the tree
    # tree.write_to_path(msa_fname.replace(".fa", ".nexus"), 'nexus')
    # my_tree = Phylo.read(msa_fname.replace(".fa", ".nexus", 'nexus')
    # fig = plt.figure(figsize=(16,18))
    # ax = fig.add_subplot(1,1,1)
    # Phylo.draw(my_tree, axes=ax)

   # t = Tree(tree.compose_newick().replace("'","") + ";") # load tree into ETE
    t = Tree(str(tree.extract_tree()).replace("'", "") + ";")  # load tree into ETE

    t.set_outgroup("spc_bear")
    print set(ins_present)
    print t
    is_monophyletic = t.check_monophyly(values=ins_present, target_attr="name")[0]
    return {"ins_monophyletic": is_monophyletic}


def get_sequence(infile, seqname, start, stop):
    print "[{}] Info: Load Sequence from {} ...".format(script, infile),
    fasta = pyfaidx.Fasta(infile)  # , filt_function=lambda x: x[0] == seqname)
    print "OK."
    fasta_id = associate_taxon(basename(infile))

    if stop < start:
        start, stop = stop, start

    seq = fasta[seqname][start:stop]
    record = SeqRecord(Seq(seq.seq), id=fasta_id, description=seq.longname)
    if set(record) == set("N"):
        print "Error: Sequence %s consists only of Ns, skipping locus." %record.id
        return False

    return record


def create_interval(start, stop, side):
    start = start - 500  # added insert size lengthed slop to account for potential SNP calling errors
    stop = stop + 500
    if side == "5prime":
        return {"start": start - options.flank,
                "stop": start - (options.flank - options.windowsize)}
    elif side == "3prime":
        return {"start": stop + (options.flank - options.windowsize),
                "stop": stop + options.flank}
    else:
        return False


def iterate_loci(matrix):

    for i_locus, locus in enumerate(matrix):
        skip_locus = False
        print "[{}] Info: Processing locus no {} - {}:{}-{}".format(script,i_locus, locus["seqname"], locus["start"], locus["stop"],
                                                            locus["taxa"])
        # fetch sequences
        sequence_collection = {"5prime": [],
                               "3prime": [],
                               "combined": []}  # Align.MultipleSeqAlignment([])

        for fasta_in in options.fasta:
            interval_5prime = create_interval(locus["start"], locus["stop"], "5prime")

            seq_5_prime = get_sequence(fasta_in, locus["seqname"],
                                       interval_5prime["start"], interval_5prime["stop"])

            interval_3prime = create_interval(locus["start"], locus["stop"], "3prime")

            seq_3_prime = get_sequence(fasta_in, locus["seqname"],
                                       interval_3prime["start"], interval_3prime["stop"])

            for coordinate in interval_3prime.values()+ interval_5prime.values():
                if coordinate < 0:
                    skip_locus = True

            for sequence in [seq_5_prime, seq_3_prime]:
                if not sequence:
                    skip_locus = True

            if skip_locus:
                break
            sequence_collection["5prime"].append(seq_5_prime)
            sequence_collection["3prime"].append(seq_3_prime)
            sequence_collection["combined"].append(seq_5_prime + seq_3_prime)





            print "[{}] Info: Total sequence length {} ({}-{} + {}-{})".format(script,
                                                                               len(seq_5_prime + seq_3_prime),
                                                                               interval_5prime["start"],
                                                                               interval_5prime["stop"],
                                                                               interval_3prime["start"],
                                                                               interval_3prime["stop"])
        if skip_locus:
            continue

        for flanktype, seqs in sequence_collection.iteritems():

            min_alignment_length = min([len(sequence) for sequence in seqs])
            msa = Align.MultipleSeqAlignment([sequence[:min_alignment_length] for sequence in seqs])  # build alignment
            # process MSA
            sites = process_alignment(msa, locus["taxa"])
            # save alignment for dendropy-popgen stats
            msa_fname = "msa_%s.fa" %flanktype
            AlignIO.write(msa, "msa_%s.fa" %flanktype, "fasta")
            monophyly = check_monophyly(msa_fname, locus["taxa"])

            for k, v in locus.iteritems():
                if k == "taxa" and len(v) >= 2:
                    sites[k] = ",".join(v)
                    monophyly[k] = ",".join(v)
                else:
                    sites[k] = v
                    monophyly[k] = v
                sites["flank"] = flanktype
            sites_stats = sites.copy()
            sites_stats.update(monophyly)
            result_table.append(sites_stats)
            print sites_stats


    with open('sites.csv', 'w') as csvfile:
        fieldnames = result_table[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, dialect="excel-tab")

        writer.writeheader()
        for line in result_table:
            writer.writerow(line)
    return 1


def main():
    matrix = load_matrix()
    iterate_loci(matrix)

    return 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This is {}".format(script))
    parser.add_argument('-i', '--fasta', required=True, action="append", help='Input files')
    parser.add_argument('-f', '--flank', required=True, type=int, help='Size of flanking region')
    parser.add_argument('-o', '--overlapping', required=False, type=bool, help='extend windows or being separate',
                        default=False)
    parser.add_argument('-w', '--windowsize', required=False, type=int, help='window size', default=1000)
    parser.add_argument('-m', '--matrix', required=True, help='matrix TSV file')
    options = parser.parse_args()

    main()
