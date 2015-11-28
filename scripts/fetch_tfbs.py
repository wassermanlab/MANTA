#!/usr/bin/env python2.7
#*-* coding: utf-8 *-*

import sys
import argparse
import pymongo

from pymongo import MongoClient

MANTA_DB_HOST = 'manta.cmmt.ubc.ca'
MANTA_DB_NAME = 'manta'
MANTA_DB_USER = 'manta_r'
MANTA_DB_PASS = 'mantapw'

def fetch_tfbs(db, chrom, start, end):
    '''
    Search the database for TFBS in the given range.
    '''

    query = {
        'chrom' : chrom,
        'end'   : {'$gte' : start},
        'start' : {'$lte' : end}
    }

    projection = {
        '_id'               : 0,
        'jaspar_tf_id'      : 1,
        'chrom'             : 1,
        'start'             : 1,
        'end'               : 1,
        'strand'            : 1
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
        tfbs_list.append(tfbs)

    return tfbs_list


def write_tfbs(filename, tfbs):
    '''
    Write the TFBS.
    '''

    fh = open(filename, 'w')

    for t in tfbs:
        fh.write("chr{0}\t{1}\t{2}\t{3}\n".format(t['chrom'], t['start'] - 1, t['end'], t['jaspar_tf_id']))

    fh.close()

    return


###############################################################################
#                               MAIN
###############################################################################
if __name__ == '__main__':
    '''
    Search the CRV database for all impacted TFBSs for a given set of
    variants.

    Usage: fetch_tfbs.py -c chrom -s start -e end -o out_file

    Where:
        -c chrom    - chromosome
        -s start    - start position
        -e end      - end position
        -o FILE     - Ouput BED file listing the TFBSs
    '''

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
        '-o', '--out_file', nargs='?', required=True, help='Output tab delimited BED file containing TFBSs'
    )

    args = parser.parse_args()

    chrom = args.chromosome
    start = args.start
    end = args.end
    out_file = args.out_file

    uri = "mongodb://{0}:{1}@{2}/{3}".format(
        MANTA_DB_USER, MANTA_DB_PASS, MANTA_DB_HOST, MANTA_DB_NAME
    )
    client = MongoClient(uri)
    db = client.manta

    tfbs = fetch_tfbs(db, chrom, int(start), int(end))

    client.close()

    write_tfbs(out_file, tfbs)
