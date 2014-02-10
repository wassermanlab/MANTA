#!/usr/bin/env python2.7
#*-* coding: utf-8 *-*

import sys
import argparse
import os
import pymongo
import logging
import datetime

from pymongo import MongoClient

LOGFILE = 'load_database.log'

# Base path for the TFBS data files. This is not actually used as all the TFBS
# information is also contained in the SNV_impacts files.
TFBS_BASE_PATH = '/export/home/amathelier/CRV_TFBS/CRV_TFBS_JASPAR2014/predicted_TFBSs'

# Base path to the TFBS SNV impacts files
#SNV_BASE_PATH = '/raid2/amathelier/SNV_impacts_disruptalt'
SNV_BASE_PATH = '/raid2/amathelier/SNV_impacts_allalt'

def clean_database(db):
    logging.info("Cleaning database\n")

    # Drop the collections
    logging.info("Dropping experiments collection\n")
    db.experiments.drop()
    logging.info("Dropping peaks collection\n")
    db.peaks.drop()
    logging.info("Dropping tfbs_snvs collection\n")
    db.tfbs_snvs.drop()

    # Re-create the indexes
    logging.info("Creating experiments indexes\n")
    db.experiments.ensure_index('jaspar_tf_id')
    db.experiments.ensure_index('tf_name')
    db.experiments.ensure_index('name')

    logging.info("Creating peaks indexes\n")
    db.peaks.ensure_index('experiment_id')
    db.peaks.ensure_index([('chrom', pymongo.ASCENDING),
                           ('start', pymongo.ASCENDING),
                           ('end', pymongo.ASCENDING)])

    logging.info("Creating tfbs_snvs indexes\n")
    db.tfbs_snvs.ensure_index([('jaspar_tf_id', pymongo.ASCENDING),
                               ('chrom', pymongo.ASCENDING),
                               ('start', pymongo.ASCENDING),
                               ('strand', pymongo.ASCENDING)])
    db.tfbs_snvs.ensure_index([('chrom', pymongo.ASCENDING),
                               ('snvs.pos', pymongo.ASCENDING)])


def load_experiments(db, master_file):
    experiments = db.experiments

    with open(master_file, 'r') as f:
        for line in f:
            line = line.rstrip()

            tf, peaks_file_path, jaspar_tf_id, jaspar_tf_name = line.split(',')

            #print "{0}\t{1}\t{2}\t{3}".format(
            #    tf, peaks_file_path, jaspar_tf_id, jaspar_tf_name
            #)

            if tf == 'TF':
                # skip header
                continue

            logging.info("Processing experiment {0}\n".format(line))

            peaks_file_name = os.path.basename(peaks_file_path)

            # Use the base name and extension of the file as the experiment
            # name and type respectively, e.g.:
            #     name = wgEncodeHaibTfbsGm12878BatfPcr1xPkRep1
            #     type = broadPeak
            exp_name, exp_type = os.path.splitext(peaks_file_name)
            exp_type = exp_type.lstrip('.')

            experiment = {
                "name"           : exp_name,
                "type"           : exp_type,
                "tf_name"        : tf,
                "jaspar_tf_id"   : jaspar_tf_id,
                "jaspar_tf_name" : jaspar_tf_name
            }

            exp_id = experiments.insert(experiment)

            experiment['id'] = exp_id

            #
            # We explicitly associate peaks and TFBSs/SNV impacts to the
            # experiment via the experiment ID. Originally, it was also
            # intended to associate the TFBSs/SNV impacts to the peaks via
            # the peak ID, but the organization of the data files makes this
            # somewhat difficult to do efficiently. It is not clear we really
            # need this explicit association anyway (they are still associated
            # via their coordinates).
            # If necessary this can done in a post processing step after all
            # the data is loaded into the database.
            #
            load_peaks(db, experiment, peaks_file_path)
            load_tfbs_snvs(db, experiment) 


def load_peaks(db, experiment, peaks_file):
    logging.info("Processing peaks file {0}\n".format(peaks_file))

    peaks = db.peaks

    exp_type = experiment['type']
    exp_id   = experiment['id']

    if exp_type == 'broadPeak':
        with open(peaks_file, 'r') as f:
            for line in f:
                line = line.rstrip()

                chrom, start, end, name, score, strand, signal, pvalue, qvalue = line.split('\t')

                chrom = chrom.lstrip('chr')

                #print "{0}\t{1}\t{2}\t{3}\t{4}".format(
                #    chrom, start, end, name, score
                #)

                # Convert start to 1-based
                start = int(start) + 1

                peak = {
                    "experiment_id" : exp_id,
                    "chrom"         : chrom,
                    "start"         : start,
                    "end"           : int(end),
                    "score"         : int(score),
                    "signal"        : float(signal)
                }

                # Store if not indicated as missing values
                if name != '.':
                    peak['name'] = name
                if pvalue != '-1':
                    peak['pvalue'] = float(pvalue)
                if qvalue != '-1':
                    peak['qvalue'] = float(qvalue)

                peak_id = peaks.insert(peak)

    elif exp_type == 'narrowPeak':
        with open(peaks_file, 'r') as f:
            for line in f:
                line = line.rstrip()

                chrom, start, end, name, score, strand, signal, pvalue, qvalue, peak_pos = line.split('\t')

                chrom = chrom.lstrip('chr')

                #print "{0}\t{1}\t{2}\t{3}\t{4}".format(
                #    chrom, start, end, score, peak_max_pos
                #)

                # Convert start to 1-based
                start = int(start) + 1
                # Convert to peak position to chromosomal coord (specified
                # in file as 0-based offset).
                peak_pos = start + int(peak_pos)

                peak = {
                    "experiment_id" : exp_id,
                    "chrom"         : chrom,
                    "start"         : start,
                    "end"           : int(end),
                    "position"      : peak_pos,
                    "score"         : int(score),
                    "signal"        : float(signal)
                }

                # Store if not indicated as missing values
                if name != '.':
                    peak['name'] = name
                if pvalue != '-1':
                    peak['pvalue'] = float(pvalue)
                if qvalue != '-1':
                    peak['qvalue'] = float(qvalue)

                peak_id = peaks.insert(peak)

    elif exp_type == 'paz':
        with open(peaks_file, 'r') as f:
            for line in f:
                line = line.rstrip()

                chrom, start, end, peak_pos, peak_center, score = line.split('\t')

                chrom = chrom.lstrip('chr')

                #print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(
                #    chrom, start, end, score, height, peak_center
                #)

                # Convert start to 1-based
                start = int(start) + 1
                peak = {
                    "experiment_id" : exp_id,
                    "chrom"         : chrom,
                    "start"         : start,
                    "end"           : int(end),
                    "center"        : peak_center
                }

                if peak_pos != '.':
                    peak['position'] = int(peak_pos)
                if score != '.':
                    peak['score'] = float(score)

                peak_id = peaks.insert(peak)
    else:
        sys.exit("Unknown peaks file format {0}".format(exp_type))


def load_tfbs_snvs(db, experiment):
    tfbs_snvs = db.tfbs_snvs

    exp_id = experiment['id']
    exp_name = experiment['name']
    exp_type = experiment['type']
    jaspar_tf_id = experiment['jaspar_tf_id']

    #
    # For computational efficiency, check if there are multiple experiments
    # for this JASPAR matrix. If not, we can safely insert tfbs_snvs without
    # checking to see if duplicate ones already exists at the same position
    # from another experiment. Otherwise for each tfbs_snv we have to check
    # and if a record already exists at the same postions we have to update
    # it by adding the current experiment ID to it rather than inserting a
    # new tfbs_snv record which is slow.
    #
    cur = db.experiments.find({'jaspar_tf_id' : jaspar_tf_id})
    exp_count = cur.count()
    check_for_update = False
    if exp_count > 1:
        check_for_update = True

    #
    # The TFBS file naming convention is:
    # <TFBS_BASE_PATH>/<jaspar_tf_id>/<exp_name>.<exp_type>_TFBSs.bed
    #
    # Note, there is also a related .txt file.
    #
    # For broad peak and narrow peak experiments, the TFBS SNV impacts file
    # naming convention is:
    # <SNV_BASE_PATH>/<jaspar_tf_id>/<exp_name>.<exp_type>_TFBSs_SNV_impacts.txt]
    #
    # For PAZAR experiments, the TFBS SNV impacts file naming convention is:
    # <SNV_BASE_PATH>/<jaspar_tf_id>/<exp_name>_TFBSs_SNV_impacts.txt]
    #
    # We do not even need the TFBS file as the SNV impacts file seems to
    # contain all the information we need
    #
    #tfbs_file = "{0}.{1}_TFBSs.bed".format(exp_name, exp_type)
    #tfbs_path = os.path.join(TFBS_BASE_PATH, jaspar_tf_id, tfbs_file)

    if exp_type == 'paz':
        snv_file = "{0}_TFBSs_SNV_impacts.txt".format(exp_name)
    else:
        snv_file = "{0}.{1}_TFBSs_SNV_impacts.txt".format(exp_name, exp_type)

    snv_path = os.path.join(SNV_BASE_PATH, jaspar_tf_id, snv_file)

    logging.info("Processing TFBS SNV file {0}\n".format(snv_path))

    tfbs_snv = {}
    first = True
    line_num = 0

    try:
        f = open(snv_path, 'r')
    except:
        errstr = "TFBS SNV impacts file {0} does not exist!".format(snv_path)
        logging.error(errstr)
    else:
        for line in f:
            line = line.rstrip()
            line_num += 1

            # NOTE: not a BED file, coords are already 1-based
            chrom, position, ignore, ref_allele, alt_allele, ref_tfbs_start, ref_tfbs_end, ref_tfbs_abs_score, ref_tfbs_rel_score, ref_tfbs_strand, snv_jaspar_tf_id, alt_tfbs_abs_score, alt_tfbs_rel_score, alt_tfbs_strand, alt_tfbs_start, alt_tfbs_end, impact, num_alt_considered, species = line.split('\t')

            if snv_jaspar_tf_id != jaspar_tf_id:
                errstr = "JASPAR TF ID '{0}' from experiment file does not match ID '{1}' from TFBS_SNV file {2} line {3}".format(jaspar_tf_id, snv_jaspar_tf_id, snv_path, line_num)
                logging.error(errstr)
                sys.exit(errstr)

            #
            # Explicitly cast elements to required types for efficient
            # storage in DB.
            #
            position = int(position)
            ref_tfbs_start = int(ref_tfbs_start)
            ref_tfbs_end = int(ref_tfbs_end)
            ref_tfbs_abs_score = float(ref_tfbs_abs_score)
            ref_tfbs_rel_score = float(ref_tfbs_rel_score)
            alt_tfbs_abs_score = float(alt_tfbs_abs_score)
            alt_tfbs_rel_score = float(alt_tfbs_rel_score)
            ref_tfbs_strand = int(ref_tfbs_strand)
            alt_tfbs_strand = int(alt_tfbs_strand)
            alt_tfbs_start = int(alt_tfbs_start)
            alt_tfbs_end = int(alt_tfbs_end)
            num_alt_considered = int(num_alt_considered)

            if first:
                # The first SNV of the file. Create a new tfbs_snv object
                tfbs_snv = {
                    'jaspar_tf_id'  : jaspar_tf_id,
                    'chrom'         : chrom,
                    'start'         : ref_tfbs_start,
                    'end'           : ref_tfbs_end,
                    'abs_score'     : ref_tfbs_abs_score,
                    'rel_score'     : ref_tfbs_rel_score,
                    'strand'        : ref_tfbs_strand,
                    'alt_considered': num_alt_considered,
                    'snvs' : []
                }

                snv = {
                    'pos'           : position,
                    'ref_allele'    : ref_allele,
                    alt_allele : {
                        'start'     : alt_tfbs_start,
                        'end'       : alt_tfbs_end,
                        'strand'    : alt_tfbs_strand,
                        'abs_score' : alt_tfbs_abs_score,
                        'rel_score' : alt_tfbs_rel_score
                    }
                }

                if impact != 'N/A':
                    snv[alt_allele]['impact'] = float(impact)

                tfbs_snv['snvs'].append(snv)

                first = False

            else:
                if chrom == tfbs_snv['chrom'] and ref_tfbs_start == tfbs_snv['start'] and ref_tfbs_strand == tfbs_snv['strand']:
                    # Same binding site. Add SNV data to existing site.

                    position_exists = False
                    for snv in tfbs_snv['snvs']:
                        if snv['pos'] == position:
                            # Same position. Add new alt. allele data.
                            position_exists = True

                            if alt_allele in snv:
                                #
                                # XXX
                                # This should not happen.
                                # Throw exception once debugged.
                                # XXX
                                #
                                errstr = "Duplicate alt. allele {0}, position {1}, ref. allele {2} from SNV impacts file {3}, line {4}:\n{5}\nexisting: {6}".format(alt_allele, position, ref_allele, snv_file, line_num, line, snv)
                                logging.error(errstr)
                                sys.exit(errstr)

                            snv[alt_allele] = {
                                'start'     : alt_tfbs_start,
                                'end'       : alt_tfbs_end,
                                'strand'    : alt_tfbs_strand,
                                'abs_score' : alt_tfbs_abs_score,
                                'rel_score' : alt_tfbs_rel_score
                            }

                            if impact != 'N/A':
                                snv[alt_allele]['impact'] = float(impact)

                            break

                    if not position_exists:
                        # New position within the binding site.
                        snv = {
                            'pos'           : position,
                            'ref_allele'    : ref_allele,
                            alt_allele : {
                                'start'     : alt_tfbs_start,
                                'end'       : alt_tfbs_end,
                                'strand'    : alt_tfbs_strand,
                                'abs_score' : alt_tfbs_abs_score,
                                'rel_score' : alt_tfbs_rel_score,
                            }
                        }

                        if impact != 'N/A':
                            snv[alt_allele]['impact'] = float(impact)

                        tfbs_snv['snvs'].append(snv)
                else:
                    # New TFBS, save previous and initialize new one
                    if check_for_update:
                        insert_or_update_tfbs_snv(db, exp_id, tfbs_snv)
                    else:
                        tfbs_snv['experiment_ids'] = [exp_id]
                        tfbs_snvs.insert(tfbs_snv)

                    tfbs_snv = {
                        'jaspar_tf_id'  : jaspar_tf_id,
                        'chrom'         : chrom,
                        'start'         : ref_tfbs_start,
                        'end'           : ref_tfbs_end,
                        'abs_score'     : ref_tfbs_abs_score,
                        'rel_score'     : ref_tfbs_rel_score,
                        'strand'        : ref_tfbs_strand,
                        'alt_considered': num_alt_considered,
                        'snvs' : []
                    }

                    snv = {
                        'pos'           : position,
                        'ref_allele'    : ref_allele,
                        alt_allele : {
                            'start'     : alt_tfbs_start,
                            'end'       : alt_tfbs_end,
                            'strand'    : alt_tfbs_strand,
                            'abs_score' : alt_tfbs_abs_score,
                            'rel_score' : alt_tfbs_rel_score
                        }
                    }

                    if impact != 'N/A':
                        snv[alt_allele]['impact'] = float(impact)

                    tfbs_snv['snvs'].append(snv)

        # Save last TFBS_SNV
        if line_num > 0 and tfbs_snv['jaspar_tf_id']:
            if check_for_update:
                insert_or_update_tfbs_snv(db, exp_id, tfbs_snv)
            else:
                tfbs_snv['experiment_ids'] = [exp_id]
                tfbs_snvs.insert(tfbs_snv)
        

def insert_or_update_tfbs_snv(db, exp_id, tfbs_snv):
    tfbs_snvs = db.tfbs_snvs

    query_tfbs_snv = {
        'jaspar_tf_id'  : tfbs_snv['jaspar_tf_id'],
        'chrom'         : tfbs_snv['chrom'],
        'start'         : tfbs_snv['start'],
        'strand'        : tfbs_snv['strand']
    }

    existing_tfbs_snv = tfbs_snvs.find_one(query_tfbs_snv)

    if existing_tfbs_snv:
        # Append this experiment ID to the tfbs_snv experiment IDs and
        # update database record.
        exp_ids = existing_tfbs_snv['experiment_ids']
        exp_ids.append(exp_id)
        tfbs_snvs.update(query_tfbs_snv, {'$set': {'experiment_ids' : exp_ids}})
    else:
        # Insert a new tfbs_snv. The tfbs_snv may refer to more than one
        # experiments, so store experiments IDs as a list.
        tfbs_snv['experiment_ids'] = [exp_id]
        tfbs_snvs.insert(tfbs_snv)


###############################################################################
#                               MAIN
###############################################################################
if __name__ == '__main__':
    '''
    Load the CRV database.

    Usage: load_database.py -i data_set.csv -u user -p pass [-l log_file] [-c]

    Where:
        -i FILE - Specifies the input data sets (experiments) file
        -l FILE - Specifies a file to which warning and error messages
                  are logged. Default = load_database.log
        -c      - If specified, the database is completely cleaned of all
                  existing documents.
        -u user - The manta database user name
        -p pass - The manta database password
    '''

    parser = argparse.ArgumentParser(
        description='Parse the various files associated with the CRVs and load the information into the database'
    )

    parser.add_argument(
        '-i', '--data_sets_file', nargs='?', required=True, help='Master CSV file containing the TF name, full path to the peaks BED file, JASPAR TF ID and JASPAR TF name'
    )

    parser.add_argument(
        '-l', '--log', nargs='?', const=LOGFILE, default=LOGFILE, help='Name of file to which logging messages are written'
    )

    parser.add_argument(
        '-c', '--clean_db', nargs='?', const=True, help='If specified, clean the database before update'
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

    data_sets_file = args.data_sets_file
    logfile = args.log
    clean_db = args.clean_db
    manta_db_host = args.host
    manta_db_name = args.dbname
    manta_db_user = args.user
    manta_db_pass = args.password

    logging.basicConfig(filename=logfile, filemode='w', level=logging.INFO)

    start_datetime = datetime.datetime.now()
    logging.info("CRV DB build started on {0}\n".format(start_datetime))
    logging.info("Data sets file: {0}\n".format(data_sets_file))

    uri = "mongodb://{0}:{1}@{2}/{3}".format(
        manta_db_user, manta_db_pass, manta_db_host, manta_db_name)
    client = MongoClient(uri)
    db = client.manta

    #
    # If -c option specified, clean database (remove all documents from all
    # collections) before insert/update of new data.
    #
    if clean_db:
        clean_database(db)

    load_experiments(db, data_sets_file)

    client.close()

    end_datetime = datetime.datetime.now()
    elapsed_time = end_datetime - start_datetime
    logging.info("CRV DB completed successfully on {0}\n".format(end_datetime))
    logging.info("Elapsed time: {0}\n".format(elapsed_time))
