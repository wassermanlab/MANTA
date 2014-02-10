#!/raid2/local/python2.7.3/bin/python2.7
#*-* coding: utf-8 *-*

import argparse
import pymongo
import bson

from pymongo import MongoClient

def remove_documents(db, jaspar_tf_id):
    '''
    Remove all documents from the experiments, peaks and tfbs_snvs collections
    for the given JASPAR TF ID
    '''

    #
    # Find all experiments with this JASPAR TF ID
    #
    for exp in db.experiments.find({'jaspar_tf_id' : jaspar_tf_id}):
        exp_id = exp['_id']

        # Delete the experiment
        print "Deleting experiment {0}\n\n".format(exp_id)
        db.experiments.remove({'_id' : bson.ObjectId(exp_id)})

        # Delete all peaks associated with this experiment
        print "Deleting all peaks with experiment ID {0}\n\n".format(exp_id)
        db.peaks.remove({'experiment_id' : bson.ObjectId(exp_id)})

    # Delete all tfbs_snvs with this JASPAR TF ID
    print "Deleting all TFBS SNVs with JASPAR TF ID {0}\n\n".format(jaspar_tf_id)
    db.tfbs_snvs.remove({'jaspar_tf_id' : jaspar_tf_id})

###############################################################################
#                               MAIN
###############################################################################
if __name__ == '__main__':
    '''
    Remove all experiments, peaks, and tfbs_snvs documents related to a
    specific JASPAR TF ID from the CRV database

    Usage: delete_by_jaspar_tf_id.py -i ID

    Where:
        -i ID   - The JASPAR TF ID for which all documents are removed
    '''

    parser = argparse.ArgumentParser(
        description='Parse the input JASPAR TF ID and delete all related documents in the experiments, peaks and tfbs_snvs collections.'
    )

    parser.add_argument(
        '-i', '--id', nargs='?', required=True, help='The JASPAR TF ID for which all experiments, peaks and tfbs_snvs records should be removed'
    )

    args = parser.parse_args()

    jaspar_tf_id = args.id

    client = MongoClient('mongodb://manta.cmmt.ubc.ca/')
    db = client.manta

    remove_documents(db, jaspar_tf_id)
