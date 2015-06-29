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

                (allele_type, indel_len) = check_alleles(ref_allele, alt_allele)
                if allele_type is None:
                    continue

                #
                # NOTE the length returned by check_alleles is the true length
                # of the indel not including the base before the indel.
                # The start refers to the position of the base before the indel.
                # The end is the last base of the indel.
                #
                start = int(position)
                end = start + indel_len

                chrom = chrom.lstrip('chr')

                var = {
                    'id'    : id,
                    'chrom' : chrom,
                    'start' : start,
                    'end'   : end,
                    'allele_type' : allele_type,
                    'ref_allele' : ref_allele,
                    'alt_allele' : alt_allele
                }

                variants.append(var)

        elif filetype == 'gff':
            sys.exit("Sorry, GFF input file format is not yet supported\n")

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
            sys.exit("Sorry, BED input file format is not yet supported\n")

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

                # Don't have ref. allele info so set to '.'
                # This is used later to determine if the input file was BED
                # and therefore the ref. allele is actually unknown.
                ref_allele = '.'

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
            sys.exit("Sorry, 'simple' input file format is not yet supported\n")

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
    Search the database for TFBS impacted by the given insertion and deletion
    events. A list of the tfbs_snv records is appended to each of the variants
    records and returned.
    '''

    impacted_tfbs = []
    for var in variants:
        indel_id = var['id']
        chrom = var['chrom']
        var_start = var['start']
        var_end = var['end']
        ref_allele = var['ref_allele']
        alt_allele = var['alt_allele']
        allele_type = var['allele_type']

        query = None
        if allele_type == 'insertion':
            print "checking TFBS overlapping insertion chr{0}:{1}-{2}\n".format(chrom, var_start, var_end)
            query = {
                'chrom' : chrom,
                'start' : {'$lte': var_start},
                'end' : {'$gt': var_start}
            }
        elif allele_type == 'deletion':
            print "checking TFBS overlapping deletion chr{0}:{1}-{2}\n".format(chrom, var_start, var_end)
            query = {
                'chrom' : chrom,
                'start' : {'$lte': var_end},
                'end' : {'$gt': var_start}
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

            impacted_tfbs.append(
                {
                    'tf_name'       : experiment['tf_name'],
                    'indel_id'      : indel_id,
                    'chrom'         : chrom,
                    'var_start'     : var_start,
                    'var_end'       : var_end,
                    'ref_allele'    : ref_allele,
                    'alt_allele'    : alt_allele,
                    'jaspar_tf_id'  : tfbs_snv['jaspar_tf_id'],
                    'tfbs_start'    : tfbs_snv['start'],
                    'tfbs_end'      : tfbs_snv['end'],
                    'strand'        : tfbs_snv['strand'],
                    'abs_score'     : tfbs_snv['abs_score'],
                    'rel_score'     : tfbs_snv['rel_score']
                }
            )

    return impacted_tfbs


def write_impacted_tfbs(filename, impacted_tfbs):
    '''
    Write out the impacted TFBS information for each insertion / deletion
    event to the given output file.
    '''

    fh = open(filename, 'w')

    fh.write("Chrom\tPosition\tRef_allele\tAlt_allele\tIndel_ID\tTF_name\tJASPAR_ID\tTFBS_start\tTFBS_end\tTFBS_strand\tTFBS_abs_score\tTFBS_rel_score\tImpact_score\n")

    for t in impacted_tfbs:
        fh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10:0.3f}\t{11:0.3f}\n".format(t['chrom'], t['var_start'], t['ref_allele'], t['alt_allele'], t['indel_id'], t['tf_name'], t['jaspar_tf_id'], t['tfbs_start'], t['tfbs_end'], '+' if t['strand'] == 1 else '-', t['abs_score'], t['rel_score']))

    fh.close()

    return


def check_alleles(ref_allele, alt_allele):
    '''
    Check alleles to make sure this looks like an insertion or deletion.
    Returns 'insertion' or 'deletion' and the length of the insertion /
    deletion as a tuple or None if this doesn't look like a simple insertion
    or deletion.
    '''

    # These are the lengths including the base before the actual start of
    # the indel
    ref_len = len(ref_allele)
    alt_len = len(alt_allele)

    allele_type = None
    indel_len = 0

    if alt_len > ref_len:
        if ref_len > 1:
            sys.stderr.write("Alleles {0} and {1} appear to indicate an insertion but the reference allele length is greater than 1.\n".format(ref_allele, alt_allele))
        else:
            allele_type = 'insertion'
            # The true length not including the base before the actual start
            # of the indel
            indel_len = alt_len - 1
    elif ref_len > alt_len:
        if alt_len > 1:
            sys.stderr.write("Alleles {0} and {1} appear to indicate a deletion but the alternate allele length is greater than 1.\n".format(ref_allele, alt_allele))
        else:
            allele_type = 'deletion'
            # The true length not including the base before the actual start
            # of the indel
            indel_len = ref_len - 1
    else:
        sys.stderr.write("Alleles {0} and {1} have the same length. Could not determine if this is an insertion or deletion.\n".format(ref_allele, alt_allele))

    return (allele_type, indel_len)



###############################################################################
#                               MAIN
###############################################################################
if __name__ == '__main__':
    '''
    Search the MANTA database for all impacted TFBSs for a given set of
    insertions / deletions.

    Usage: search_impacts.py -i variants_file [-t file_type] -o out_file

    Where:
        -i FILE     - Input file containing a list of variants. The file
                      should be in BED, VCF or GFF format.

        -t type     - The input file type: 'VCF', 'BED' or 'GFF'.
                      Default = 'VCF'
        -o FILE     - Ouput file listing the TFBSs impacted by the given
                      insertions / deletions.
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

    impacted_tfbs = search_database(db, variants)

    client.close()

    write_impacted_tfbs(out_file, impacted_tfbs)
