#!/usr/bin/env python2.7
#*-* coding: utf-8 *-*

import argparse
import pymongo

from pymongo import MongoClient

def clean_database(db):
    '''
    Clean the database by dropping all collections.
    Then recreate the secondary indexes.
    '''

    # Drop the collections
    db.experiments.drop()
    db.peaks.drop()
    db.tfbs_snvs.drop()

    # Re-create the indexes
    db.experiments.ensure_index('jaspar_tf_id')
    db.experiments.ensure_index('tf_name')
    db.experiments.ensure_index('name')

    db.peaks.ensure_index('experiment_id')
    db.peaks.ensure_index([('chrom', pymongo.ASCENDING),
                           ('start', pymongo.ASCENDING),
                           ('end', pymongo.DESCENDING)])

    db.tfbs_snvs.ensure_index([('jaspar_tf_id', pymongo.ASCENDING),
                               ('chrom', pymongo.ASCENDING),
                               ('start', pymongo.ASCENDING),
                               ('strand', pymongo.ASCENDING)])
    db.tfbs_snvs.ensure_index([('chrom', pymongo.ASCENDING),
                               ('snvs.pos', pymongo.ASCENDING)])


###############################################################################
#                               MAIN
###############################################################################
if __name__ == '__main__':
    '''
    Clean the CRV database by removing all documents from the experiments,
    peaks and tfbs_snvs collections

    Usage: clean_database.py -h host -d name -u user -p pass

    Where:
        -h host - The manta database host
        -d name - The manta database name
        -u user - The manta database user name
        -p pass - The manta database password
    '''

    parser = argparse.ArgumentParser(
        description='Clean the CRV database by removing all documents from the experiments, peaks and tfbs_snvs collections'
    )

    parser.add_argument(
        '-h', '--host', nargs='?', required=True, help='Manta database host'
    )

    parser.add_argument(
        '-d', '--dbname', nargs='?', required=True, help='Manta database name'
    )

    parser.add_argument(
        '-u', '--user', nargs='?', required=True, help='Manta database user (with read/write privileges)'
    )

    parser.add_argument(
        '-p', '--password', nargs='?', required=True, help='Manta database password'
    )

    args = parser.parse_args()

    manta_db_host = args.host
    manta_db_name = args.dbname
    manta_db_user = args.user
    manta_db_pass = args.password

    uri = "mongodb://{0}:{1}@{2}/{3}".format(
        manta_db_user, manta_db_pass, manta_db_host, manta_db_name
    )
    client = MongoClient(uri)
    db = client.manta

    clean_database(db)

    client.close()
