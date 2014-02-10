#!/usr/bin/env python2.7
#*-* coding: utf-8 *-*

import sys
import argparse
import re
import pymongo
import bson

from pymongo import MongoClient

MANTA_DB_HOST = 'manta.cmmt.ubc.ca'
MANTA_DB_NAME = 'manta'
MANTA_DB_USER = 'manta_r'
MANTA_DB_PASS = 'mantapw'

REF_ALLELE_REGEX = re.compile('[ACGTNactgn]')
ALT_ALLELE_REGEX = re.compile('[ACGTactg]')


def read_variants_file(filename, filetype):
    '''
    Read the input variants file.
    Return the variants as a list.
    '''

    #sys.stderr.write("Reading variants file {0} of type {1}\n".format(filename, filetype))

    variants = []
    with open(filename, 'r') as f:
        if filetype == 'vcf':
            for line in f:
                if line.startswith("#"):
                    continue

                line = line.rstrip()

                cols = line.split('\t')

                if len(cols) < 5:
                    sys.exit("VCF file contains less than 5 columns")

                chrom       = cols[0]
                position    = cols[1]
                id          = cols[2]
                ref_allele  = cols[3]
                alt_allele  = cols[4]

                if not check_alleles(ref_allele, alt_allele):
                    continue

                chrom = chrom.lstrip('chr')

                var = {
                    'id'    : id,
                    'chrom' : chrom,
                    'position' : int(position),
                    'ref_allele' : ref_allele,
                    'alt_allele' : alt_allele
                }

                variants.append(var)

        elif filetype == 'gff':
            for line in f:
                if line.startswith("#"):
                    continue

                line = line.rstrip()

                cols = line.split('\t')

                if len(cols) < 9:
                    sys.exit("GFF file contains no attributes column")

                chrom, data_source, feature_type, start, end, score, strand, frame, attribute_list = cols

                if start != end:
                    # In this case, assume we have an indel rather than
                    # an SNV and ignore it
                    continue
                
                ref_allele = None
                alt_allele = None

                if attribute_list:
                    attributes = attribute_list.split('; ')
                    for attr in attributes:
                        attr_name, attr_val = attr.split('=')
                        if attr_name == 'ref_allele' or attr_name == 'reference_allele':
                            ref_allele = attr_val
                        if attr_name == 'alt_allele' or attr_name == 'alt_allele':
                            alt_allele = attr_val
                else:
                    sys.exit("GFF file has blank attributes field")

                if not ref_allele or not alt_allele:
                    sys.exit("GFF attributes field contains no ref_allele or alt_allele information")

                if not check_alleles(ref_allele, alt_allele):
                    continue

                chrom = chrom.lstrip('chr')

                #
                # XXX
                # If strand is given and negative, should we reverse
                # complement the alleles?
                # XXX
                #

                #
                # Assume the name field holds the alt. allele. We don't
                # know the ref. allele, so set it to 'N'
                #
                var = {
                    'id'    : '.',
                    'chrom' : chrom,
                    'position' : int(start),
                    'ref_allele' : ref_allele,
                    'alt_allele' : alt_allele
                }

                variants.append(var)

        elif filetype == 'bed':
            for line in f:
                if line.startswith("#"):
                    continue

                line = line.rstrip()

                #sys.stderr.write("BED line: {0}\n".format(line))

                cols = line.split('\t')
                
                if len(cols) < 4:
                    sys.exit("BED file has less than 4 columns")
                chrom = cols[0]
                start = int(cols[1])
                end   = int(cols[2])
                # NOTE: using the name field as alt. allele
                alt_allele = cols[3]
                # Don't have ref. allele info so set to 'N'
                ref_allele = 'N'

                if start < end - 1:
                    # In this case, assume we have an indel rather than
                    # an SNV and ignore it
                    continue

                if len(cols) >= 6:
                    score  = cols[4]
                    strand = cols[5]

                #
                # XXX
                # If strand is given and negative, should we reverse
                # complement the alleles?
                # XXX
                #

                if not check_alleles(ref_allele, alt_allele):
                    continue

                chrom = chrom.lstrip('chr')

                #
                # Assume the name field holds the alt. allele. We don't
                # know the ref. allele, so set it to 'N'
                #
                var = {
                    'id'    : '.',
                    'chrom' : chrom,
                    'position' : int(end),
                    'ref_allele' : ref_allele,
                    'alt_allele' : alt_allele
                }

                variants.append(var)

        elif filetype == 'simple':
            #
            # My own test file type. Not used...
            # Format:
            #   chrom  position  ref_allele  alt_allele  snv_id(optional)
            #
            for line in f:
                if line.startswith("#"):
                    continue

                line = line.rstrip()

                # NOTE: assuming 1-based coordinates for position
                cols = line.split('\t')

                if len(cols) < 4:
                    sys.exit("SNV 'simple' file contains less than 4 columns")

                chrom      = cols[0]
                position   = cols[1]
                ref_allele = cols[2]
                alt_allele = cols[3]

                id = '.'
                if len(cols) >= 5:
                    id = cols[4]

                chrom = chrom.lstrip('chr')

                var = {
                    'id'    : id,
                    'chrom' : chrom,
                    'position' : int(position),
                    'ref_allele' : ref_allele,
                    'alt_allele' : alt_allele
                }

                variants.append(var)

        else:
            sys.exit("Unknown SNV file type {0}".format(filetype))

    return variants


def search_database(db, variants):
    '''
    Search the database for TFBS impacts at each position in the list of
    variants. A list of the tfbs_snv records is appended to each of the
    variants records and returned.
    '''

    snv_impacts = []
    for var in variants:
        snv_id = var['id']
        chrom = var['chrom']
        pos = var['position']
        ref_allele = var['ref_allele']
        alt_allele = var['alt_allele']

        #query = {
        #    'chrom' : chrom,
        #    'start' : {'$lte': pos},
        #    'end' : {'$gte': pos}
        #}

        query = {
            'chrom' : chrom,
            'snvs.pos' : pos
        }

        for tfbs_snv in db.tfbs_snvs.find(query):
            #
            # Get the TF name from the experiment(s) associated with the
            # tfbs_snv. A tfbs_snv may have multiple experiments associated
            # with it. All experiments for a given tfbs_snv are for the same
            # JASPAR TF ID and therefore *should* have the same TF name, so
            # we only need to look at the first experiment.
            #
            exp_ids = tfbs_snv['experiment_ids']
            experiment = db.experiments.find_one(
                {'_id' : bson.ObjectId(exp_ids[0])}
            )

            snvs = tfbs_snv['snvs']
            for snv in snvs:
                if snv['pos'] == pos:
                    if snv['ref_allele'] != ref_allele:
                        sys.stderr.write("Ref allele mismatch {0} vs. {1} for TFBS {2} chr{3}:{4}-{5} at SNV position {6}\n".format(snv['ref_allele'], ref_allele, tfbs_snv['jaspar_tf_id'], chrom, tfbs_snv['start'], tfbs_snv['end'], pos))

                        break

                    if alt_allele in snv:
                        impact = snv[alt_allele]

                        snv_impacts.append(
                            {
                                'tf_name'       : experiment['tf_name'],
                                'snv_id'        : snv_id,
                                'chrom'         : chrom,
                                'position'      : pos,
                                'ref_allele'    : ref_allele,
                                'alt_allele'    : alt_allele,
                                'jaspar_tf_id'  : tfbs_snv['jaspar_tf_id'],
                                'start1'        : tfbs_snv['start'],
                                'end1'          : tfbs_snv['end'],
                                'strand1'       : tfbs_snv['strand'],
                                'abs_score1'    : tfbs_snv['abs_score'],
                                'rel_score1'    : tfbs_snv['rel_score'],
                                'start2'        : impact['start'],
                                'end2'          : impact['end'],
                                'strand2'       : impact['strand'],
                                'abs_score2'    : impact['abs_score'],
                                'rel_score2'    : impact['rel_score'],
                                'impact'        : impact['impact']
                            }
                        )

                        break
                    else:
                        sys.stderr.write("Alt. allele {0} not found in TFBS {1} chr{2}:{3}-{4} at position {5}\n".format(alt_allele, tfbs_snv['jaspar_tf_id'], chrom, tfbs_snv['start'], tfbs_snv['end'], pos))
                        break

    return snv_impacts


def write_snv_impacts(filename, snv_impacts):
    '''
    Write the TFBS impacts for each variant to the given output file.
    '''

    fh = open(filename, 'w')

    for si in snv_impacts:
        fh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10:0.3f}\t{11:0.3f}\t{12}\t{13}\t{14}\t{15:0.3f}\t{16:0.3f}\t{17:0.3f}\n".format(si['chrom'], si['position'], si['ref_allele'], si['alt_allele'], si['snv_id'], si['tf_name'], si['jaspar_tf_id'], si['start1'], si['end1'], si['strand1'], si['abs_score1'], si['rel_score1'], si['start2'], si['end2'], si['strand2'], si['abs_score2'], si['rel_score2'], si['impact']))

    fh.close()

    return


def check_alleles(ref_allele, alt_allele):
    '''
    Check alleles to make sure this looks like an SNV
    '''

    if len(ref_allele) != 1:
        return False
    if len(alt_allele) != 1:
        return False
    if ref_allele.upper() == alt_allele.upper():
        return False
    if not REF_ALLELE_REGEX.match(ref_allele):
        return False
    if not ALT_ALLELE_REGEX.match(alt_allele):
        return False

    return True



###############################################################################
#                               MAIN
###############################################################################
if __name__ == '__main__':
    '''
    Search the CRV database for all impacted TFBSs for a given set of
    variants.

    Usage: search_impacts.py -i variants_file -o out_file

    Where:
        -i FILE     - Input file containing a list of variants. The file
                      should contain the following tab delimited columns:
                          chromosome  position  ref_allele  alt_allele

        -o FILE     - Ouput file listing the TFBSs impacted by the given
                      variants.
    '''

    parser = argparse.ArgumentParser(
        description='Parse the input variant position file to search the database for damaging impacts to TFBSs'
    )

    parser.add_argument(
        '-i', '--variant_file', nargs='?', required=True, help='Input tab delimited file containing variants to search the database'
    )

    parser.add_argument(
        '-t', '--file_type', nargs='?', help='Input file type - one of: vcf, bed, gff or simple. Default = vcf'
    )

    parser.add_argument(
        '-o', '--out_file', nargs='?', required=True, help='Output tab delimited file containing results of TFBSs disrupted by variants'
    )

    args = parser.parse_args()

    variant_file = args.variant_file
    out_file = args.out_file
    file_type = args.file_type

    if file_type == None:
        file_type = 'vcf'

    file_type = file_type.lower()

    variants = read_variants_file(variant_file, file_type)

    uri = "mongodb://{0}:{1}@{2}/{3}".format(
        MANTA_DB_USER, MANTA_DB_PASS, MANTA_DB_HOST, MANTA_DB_NAME
    )
    client = MongoClient(uri)
    db = client.manta

    variant_impacts = search_database(db, variants)

    client.close()

    write_snv_impacts(out_file, variant_impacts)
