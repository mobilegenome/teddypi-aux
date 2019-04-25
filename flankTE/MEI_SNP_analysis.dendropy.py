#!/usr/env python2.7

__author__ = 'newuser'

import csv
from Bio import Align
from Bio import AlignIO
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
script = argv[0]

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


def get_sequence(infile, seqname, start, stop):
    print "[{}] Info: Load Sequence from {} ...".format(script, infile),
    fasta = pyfaidx.Fasta(infile)  # , filt_function=lambda x: x[0] == seqname)
    print "OK."
    fasta_id = associate_taxon(basename(infile))
    seq = fasta[seqname][start:stop]
    record = SeqRecord(Seq(seq.seq), id=fasta_id, description=seq.longname)

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
    for locus in matrix:
        print "[{}] Info: Processing locus {}:{}-{}".format(script, locus["seqname"], locus["start"], locus["stop"],
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

            sequence_collection["5prime"].append(seq_5_prime)
            sequence_collection["3prime"].append(seq_3_prime)
            sequence_collection["combined"].append(seq_5_prime + seq_3_prime)

            print "[{}] Info: Total sequence length {} ({}-{} + {}-{})".format(script,
                                                                               len(seq_5_prime + seq_3_prime),
                                                                               interval_5prime["start"],
                                                                               interval_5prime["stop"],
                                                                               interval_3prime["start"],
                                                                               interval_3prime["stop"])

        for flanktype, seqs in sequence_collection.iteritems():

            min_alignment_length = min([len(sequence) for sequence in seqs])


            msa = Align.MultipleSeqAlignment([sequence[:min_alignment_length] for sequence in seqs])  # build alignment
            msa_fname = "msa_%s.fa" %flanktype
            AlignIO.write(msa, "msa_%s.fa" %flanktype, "fasta")
            # process MSA
            sites = process_alignment(msa_fname, locus["taxa"])
            for k, v in locus.iteritems():
                if k == "taxa" and len(v) >= 2:
                    sites[k] = ",".join(v)
                else:
                    sites[k] = v
                sites["flank"] = flanktype
            print sites
            result_table.append(sites)

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
