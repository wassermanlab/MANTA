#!/usr/bin/env python2.7
#*-* coding: utf-8 *-*

"""
Search MANTA for all potential SNVs in given set of regions.

Usage: fetch_tfbs.py -b bed -o out_file

Where:
    -b bed      - bed file containing the regions where to search for
                  variants
    -o file     - ouput file listing the SNV and TFBS information

"""

import argparse
import pymongo
import bson

from pymongo import MongoClient

MANTA_DB_HOST = 'manta.cmmt.ubc.ca'
MANTA_DB_NAME = 'manta'
MANTA_DB_USER = 'manta_r'
MANTA_DB_PASS = 'mantapw'


def fetch_tfbs(db, chrom, start, end):
    """ Search the database for TFBS in the given range. """

    query = {'chrom': chrom, 'end': {'$gte': start},
             'start': {'$lte': end}}

    projection = {'_id': 0}

    #
    # Sorting also slowed down query somewhat and is probably not
    # necessary.
    #
    sortfields = [
        ("start", pymongo.ASCENDING),
        ("jaspar_tf_id", pymongo.ASCENDING)
    ]

    tfbs_list = []
    for tfbs in db.tfbs_snvs.find(query, projection).sort(sortfields):
        #
        # Get the JASPAR TF NAME from the first experiment
        #
        exp_ids = tfbs['experiment_ids']
        experiment = db.experiments.find_one(
            {'_id': bson.ObjectId(exp_ids[0])}
        )
        tfbs['jaspar_tf_name'] = experiment['jaspar_tf_name']

        tfbs_list.append(tfbs)

    return tfbs_list


def fetch_variants(db, bedfile):
    """ Fetch the variants in the regions from the BED file. """

    var = []
    with open(bedfile) as stream:
        for line in stream:
            spl = line.split('\t')
            chrom = spl[0]
            start = eval(spl[1]) + 1
            end = eval(spl[2])
            var.extend(fetch_tfbs(db, chrom, start, end))
    return var


def write_variants(filename, tfbs):
    """ Write the variants. """

    fh = open(filename, 'w')
    # Write full TFBS/SNV detailed information in the same format as the
    # website results table.
    #
    for t in tfbs:
        for snv in t['snvs']:
            pos = snv['pos']
            ref_allele = snv['ref_allele']

            for alt_allele in ['A', 'C', 'G', 'T']:
                if alt_allele in snv:
                    impact = snv[alt_allele]

                    fh.write("chr{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.1f}%\t{}\t{}\t{}\t{:.3f}\t{:.1f}%\t{:.3f}\n".format(
                        t['chrom'],
                        pos,
                        ref_allele,
                        alt_allele,
                        t['jaspar_tf_name'],
                        t['jaspar_tf_id'],
                        t['start'],
                        t['end'],
                        '-' if t['strand'] == -1 else '+',
                        t['abs_score'],
                        t['rel_score'] * 100,
                        impact['start'],
                        impact['end'],
                        '-' if impact['strand'] == -1 else '+',
                        impact['abs_score'],
                        impact['rel_score'] * 100,
                        impact['impact']))
    fh.close()


###############################################################################
#                               MAIN
###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Fetch SNVs in the given ranges.'
    )

    parser.add_argument('-b', '--bed', nargs='?', required=True,
                        help='BED file with regions coord')

    parser.set_defaults(write_var=False)

    parser.add_argument('-o', '--out_file', nargs='?', required=True,
                        help='Output tab delimited file containing SNVs')

    args = parser.parse_args()

    bed = args.bed
    out_file = args.out_file

    uri = "mongodb://{0}:{1}@{2}/{3}".format(MANTA_DB_USER, MANTA_DB_PASS,
                                             MANTA_DB_HOST, MANTA_DB_NAME)
    client = MongoClient(uri)
    database = client.manta

    variants = fetch_variants(database, bed)

    client.close()

    write_variants(out_file, variants)
