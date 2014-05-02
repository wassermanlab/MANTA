#!/usr/bin/env python2.7

import sys

sys.path.append("/usr/local/src/Template-Python")
#sys.path.append("/apps")

import os
import string
import tempfile
import re
import glob
import shutil
import time
import pymongo
import bson

from pymongo import MongoClient

#
# XXX Overriding standard Bio with Bio_dev to have access to latest JASPAR
# modules. Once the JASPAR modules are incorporated into the standard Biopython
# release, remove this override.
#
#from Bio.motifs.jaspar.db import JASPAR5

from cgi_app import CGI_Application
from template import Template
from subprocess import Popen
from constants import (
    ADMIN_EMAIL, HTDOCS_ABS_PATH, HTDOCS_REL_PATH, CGI_BIN_ABS_PATH,
    CGI_BIN_REL_PATH, RESULTS_ABS_PATH, RESULTS_REL_PATH,
    HTDOCS_ABS_TEMPLATE_PATH,
    MANTA_DB_HOST, MANTA_DB_NAME, MANTA_DB_USER, MANTA_DB_PASS,
    JASPAR_URL, JASPAR_DB_HOST, JASPAR_DB_NAME, JASPAR_DB_USER, JASPAR_DB_PASS,
    REMOVE_RESULTFILES_OLDER_THAN, REF_ALLELE_REGEX, ALT_ALLELE_REGEX
)


class MantaWebapp(CGI_Application):
    def setup(self):
        self.start_mode = 'start'

    def start(self):
        '''
        Display the start (home) page
        '''
        #vars = {
        #    'htdocs_rel_path' : HTDOCS_REL_PATH,
        #    'cgi_bin_rel_path' : CGI_BIN_REL_PATH,
        #    'header' : 'MANTA Home',
        #    'title' : 'Welcome to MANTA'
        #}

        #return self.process_template("start.html", vars)
        return self.snv_input()

    def snv_input(self):
        '''
        Display the page for inputting SNVs
        '''

        #jdb = JASPAR5(
        #    host=JASPAR_DB_HOST,
        #    name=JASPAR_DB_NAME,
        #    user=JASPAR_DB_USER,
        #    password=JASPAR_DB_PASS
        #)

        #motifs = jdb.fetch_motif_set(collection='CORE')
        #motifs = jdb.fetch_motifs(collection='CORE')

        vars = {
            'htdocs_rel_path' : HTDOCS_REL_PATH,
            'cgi_bin_rel_path' : CGI_BIN_REL_PATH,
            'header' : 'Search MANTA',
            'title' : 'Search for SNVs Impacting TFBSs',
            #'motifs' : sorted(motifs, compare_motifs)
        }

        return self.process_template("snv_input.html", vars)


    def search_snvs(self):
        '''
        Retrieve user entered SNVs
        Search the Manta database for TFBSs impacted by the given SNVs
        and display the results.
        '''

        snv_paste = self.param("snv_paste")
        snv_filetype = self.param("snv_filetype")

        snv_upload = self.form["snv_upload"]

        #
        # XXX
        # Need to check if we already created a temp directory and if so,
        # use that. This information has to be stored in a persistent
        # manner.
        # XXX
        #
        result_dir = tempfile.mkdtemp(dir = RESULTS_ABS_PATH) 

        snv_file = None
        if snv_paste:
            snv_file = os.path.join(result_dir, "snvs.txt")

            out = open(snv_file, "w")

            out.write(snv_paste)
            out.close()

        elif snv_upload.filename:
            snv_file = upload_named_file(
                snv_upload.file, result_dir, "snvs.txt"
            )

        else:
            return self.html_error("No SNVs were pasted or uploaded")

        snvs = self.read_variants_file(snv_file, snv_filetype)

        snv_impacts = self.search_database(snvs)

        snv_impacts_text_file = os.path.join(
            result_dir, "snv_impact_results.txt"
        )

        self.write_snv_impacts(snv_impacts_text_file, snv_impacts)

        # Convert text file path to URL for linking on results page
        snv_impacts_file_url = snv_impacts_text_file.replace(
            HTDOCS_ABS_PATH, HTDOCS_REL_PATH
        )

        vars = {
            'htdocs_rel_path' : HTDOCS_REL_PATH,
            'cgi_bin_rel_path' : CGI_BIN_REL_PATH,
            'header' : 'MANTA Results',
            'title' : 'SNVs Impacting TFBSs',
            'snv_impacts' : snv_impacts,
            'snv_impacts_file_url' : snv_impacts_file_url
        }

        return self.process_template("snv_impact_results.html", vars)


    def read_variants_file(self, filename, filetype):
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
                        return self.html_error("VCF file contains less than 5 columns")

                    chrom       = cols[0]
                    position    = cols[1]
                    id          = cols[2]
                    ref_allele  = cols[3]
                    alt_allele  = cols[4]

                    if not self.check_alleles(ref_allele, alt_allele):
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
                        return self.html_error("GFF file contains no attributes column")

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
                        return self.html_error("GFF file has blank attributes field")

                    if not ref_allele or not alt_allele:
                        return self.html_error("GFF attributes field contains no ref_allele or alt_allele information")

                    if not self.check_alleles(ref_allele, alt_allele):
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
                        'chrom' : chrom,
                        'id'    : '.',
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
                        return self.html_error("BED file has less than 4 columns")
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

                    if not self.check_alleles(ref_allele, alt_allele):
                        continue

                    chrom = chrom.lstrip('chr')

                    #
                    # Assume the name field holds the alt. allele. We don't
                    # know the ref. allele, so set it to 'N'
                    #
                    var = {
                        'chrom' : chrom,
                        'id'    : '.',
                        'position' : int(end),
                        'ref_allele' : ref_allele,
                        'alt_allele' : alt_allele
                    }

                    variants.append(var)

            elif filetype == 'simple':
                #
                # My own test file type. Not used...
                #
                for line in f:
                    if line.startswith("#"):
                        continue

                    line = line.rstrip()

                    # NOTE: assuming 1-based coordinates for position
                    cols = line.split('\t')

                    if len(cols) < 4:
                        return self.html_error("SNV file contains less than 4 columns")

                    chrom      = cols[0]
                    position   = cols[1]
                    ref_allele = cols[2]
                    alt_allele = cols[3]

                    chrom = chrom.lstrip('chr')

                    var = {
                        'chrom' : chrom,
                        'id'    : '.',
                        'position' : int(position),
                        'ref_allele' : ref_allele,
                        'alt_allele' : alt_allele
                    }

                    variants.append(var)

            else:
                return self.html_error("Unknown SNV file type {0}".format(filetype))

        return variants


    def search_database(self, variants):
        '''
        Search the database for TFBS impacts at each position in the list of
        variants. A list of the tfbs_snv records is appended to each of the
        variants records and returned.
        '''

        uri = "mongodb://{0}:{1}@{2}/{3}".format(
            MANTA_DB_USER, MANTA_DB_PASS, MANTA_DB_HOST, MANTA_DB_NAME
        )
        client = MongoClient(uri)
        db = client.manta

        snv_impacts = []
        for var in variants:
            snv_id = var['id']
            chrom = var['chrom']
            pos = var['position']
            ref_allele = var['ref_allele']
            alt_allele = var['alt_allele']

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


    def write_snv_impacts(self, filename, snv_impacts):
        '''
        Write the TFBS impacts for each variant to the given output file.
        '''

        fh = open(filename, 'w')

        for si in snv_impacts:
            fh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10:0.3f}\t{11:0.3f}\t{12}\t{13}\t{14}\t{15:0.3f}\t{16:0.3f}\t{17:0.3f}\n".format(si['chrom'], si['position'], si['ref_allele'], si['alt_allele'], si['snv_id'], si['tf_name'], si['jaspar_tf_id'], si['start1'], si['end1'], si['strand1'], si['abs_score1'], si['rel_score1'], si['start2'], si['end2'], si['strand2'], si['abs_score2'], si['rel_score2'], si['impact']))

        fh.close()

        return


    def check_alleles(self, ref_allele, alt_allele):
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


    def html_error(self, error):
        '''
        Display error as an html page
        '''

        sys.stderr.write(error)

        vars = {
            'htdocs_rel_path' : HTDOCS_REL_PATH,
            'cgi_bin_rel_path' : CGI_BIN_REL_PATH,
            'header' : 'MANTA Error',
            'title' : 'Error',
            'error' : error
        }

        return self.process_template("error.html", vars)


    def process_template(self, template_file, vars):
        config = {
            'ABSOLUTE'      : 1,
            'INCLUDE_PATH'  : HTDOCS_ABS_TEMPLATE_PATH,
            'INTERPOLATE'   : 1,
            'POST_CHOMP'    : 1,
            'EVAL_PYTHON'   : 1,
            #'PRE_PROCESS'   : 'header'
        }

        template = Template(config)

        return template.process(template_file, vars)


    def teardown(self):
        #sys.stderr.write("teardown called")
        self.clean_resultfiles()
    
        return


    def clean_resultfiles(self):
        #sys.stderr.write("clean_resultfiles called")

        rm_time = time.time() - REMOVE_RESULTFILES_OLDER_THAN * 24 * 60 * 60

        for path in glob.iglob(os.path.join(RESULTS_ABS_PATH, '*')):
            if os.path.isdir(path):
                st = os.stat(path)
                mtime = st.st_mtime
                if mtime < rm_time:
                    shutil.rmtree(path, True)

        return
        

#
# Upload/copy a user supplied file to the local destination directory with
# the supplied file name.
#
def upload_named_file(fh, destdir, out_file):
    out_file = os.path.join(destdir, out_file)

    out_fh = open(out_file, "w")

    while 1:
        chunk = fh.read(100000)
        if not chunk: break
        out_fh.write(chunk)

    out_fh.close()

    return out_file


#
# Upload/copy a user supplied file to the local destination directory creating
# a temporary file name. Return the name of the locally created file.
#
def upload_temp_file(fh, destdir, prefix, suffix):
    out_file = tempfile.mkstemp(
        prefix = prefix, suffix = suffix, dir = destdir
    )[1]

    out_fh = open(out_file, "w")

    #while 1:
    #    chunk = fh.read(100000)
    #    if not chunk: break
    #    out_fh.write(chunk)
    out_fh.write(fh.read())

    out_fh.close()

    return out_file


def compare_motifs(a, b):
    return cmp(string.upper(a.name), string.upper(b.name))


# this could be run as a mod_python handler:
def handler(req):
    e = MantaWebapp(req)
    return e.run()

# or as a regular CGI app...
if __name__ == "__main__":
    e = MantaWebapp()
    e.run()
