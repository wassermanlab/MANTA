#!/usr/bin/env python2.7
#*-* coding: utf-8 *-*

import sys
import argparse
import pymongo
import bson

from pymongo import MongoClient

MANTA_DB_HOST = 'manta.cmmt.ubc.ca'
MANTA_DB_NAME = 'manta'
MANTA_DB_USER = 'manta_r'
MANTA_DB_PASS = 'mantapw'

def fetch_tfbs(db, chrom, start, end, get_exp):
    """Search the database for TFBS in the given range.
    """

    query = {
        'chrom' : chrom,
        'end'   : {'$gte' : start},
        'start' : {'$lte' : end}
    }

    projection = {
        '_id'               : 0,
        #'jaspar_tf_id'      : 1,
        #'chrom'             : 1,
        #'start'             : 1,
        #'end'               : 1,
        #'strand'            : 1,
        #'rel_score'         : 1,
        #'experiment_ids'    : 1
    }

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
        if get_exp:
            #
            # XXX
            # If we have a lot of TFBS from different experiments it would
            # be more efficient to pre-select all the experiments from the
            # DB and build a hash keyed on experiment ID.
            # XXX
            #
            exp_names = []
            experiments = []
            exp_ids = tfbs['experiment_ids']
            for exp_id in exp_ids:
                experiment = db.experiments.find_one(
                    {'_id' : bson.ObjectId(exp_id)}
                )

                experiments.append(experiment)
                exp_names.append(experiment['name'])

            tfbs['experiment_names'] = exp_names
            tfbs['jaspar_tf_name'] = experiments[0]['jaspar_tf_name']
        else:
            #
            # Get the JASPAR TF NAME from the first experiment 
            #
            exp_ids = tfbs['experiment_ids']
            experiment = db.experiments.find_one(
                {'_id' : bson.ObjectId(exp_ids[0])}
            )
            tfbs['jaspar_tf_name'] = experiment['jaspar_tf_name']

        tfbs_list.append(tfbs)

    return tfbs_list


def write_tfbs(filename, tfbs, have_experiments, write_var):
    """Write the TFBS.
    """

    fh = open(filename, 'w')

    if write_var:
        #
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
    else:
        #
        # Just write the TFBS information in BED format
        #
        for t in tfbs:
            fh.write("chr{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(t['chrom'], t['start'] - 1, t['end'], t['jaspar_tf_name'], '%.0f'%(float(t['rel_score']) * 1000), '-' if t['strand'] == -1 else '+'))

    fh.close()

    return


def write_tfbs_experiments(filename, tfbs):
    """Write the TFBS / experiment relationships.
    """

    fh = open(filename, 'w')

    for t in tfbs:
        exp_names = t['experiment_names']

        fh.write("{0}\t{1}\t{2}\n".format(t['jaspar_tf_id'], t['jaspar_tf_name'], ','.join(exp_names)))

    fh.close()

    return


###############################################################################
#                               MAIN
###############################################################################
if __name__ == '__main__':
    """Search the CRV database for all impacted TFBSs for a given set of
    variants.

    Usage: fetch_tfbs.py -c chrom -s start -e end -o out_file

    Where:
        -c chrom    - chromosome
        -s start    - start position
        -e end      - end position
        -x exp_file - file to which the relationship between the TFBS and
                      experiments is output
        -v          - flag indicating that variant detail should also be
                      written
        -o FILE     - ouput BED file listing the TFBSs
    """

    parser = argparse.ArgumentParser(
        description='Fetch TFBSs in the given range. Output as BED'
    )

    parser.add_argument(
        '-c', '--chromosome', nargs='?', required=True, help='Chromosome name'
    )

    parser.add_argument(
        '-s', '--start', nargs='?', help='Start position on chromosome'
    )

    parser.add_argument(
        '-e', '--end', nargs='?', help='End position on chromosome'
    )

    parser.add_argument(
        '-x', '--exp_file', nargs='?', help='Optionally output a separate file showing the relationship between TFBS and the experiment to which this binding site relates to (can be a one-to-many relationship)'
    )

    parser.add_argument(
        '-v', '--write_variant_info', dest='write_var', action='store_true', help='Optionally include variant details (variant alleles/impact scores) in the output file'
    )

    parser.set_defaults(write_var=False)

    parser.add_argument(
        '-o', '--out_file', nargs='?', required=True, help='Output tab delimited BED file containing TFBSs'
    )

    args = parser.parse_args()

    chrom = args.chromosome
    start = args.start
    end = args.end
    out_file = args.out_file
    exp_file = args.exp_file
    write_var = args.write_var

    uri = "mongodb://{0}:{1}@{2}/{3}".format(
        MANTA_DB_USER, MANTA_DB_PASS, MANTA_DB_HOST, MANTA_DB_NAME
    )
    client = MongoClient(uri)
    db = client.manta

    get_experiments = 0
    if exp_file:
        get_experiments = 1

    tfbs = fetch_tfbs(db, chrom, int(start), int(end), get_experiments)

    client.close()

    write_tfbs(out_file, tfbs, get_experiments, write_var)

    if exp_file:
        write_tfbs_experiments(exp_file, tfbs)
