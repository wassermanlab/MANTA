"""
Define various constants needed by the Manta web application (including
standalone scripts).
"""
import re

ADMIN_EMAIL         = 'dave@cmmt.ubc.ca'
ERROR_LOG           = '/var/log/httpd/error_log'
SERVER_URL          = 'http://manta.cmmt.ubc.ca'
WWW_ROOT            = '/var/www'
HTDOCS_ROOT         = WWW_ROOT + '/html'
CGI_BIN_ROOT        = WWW_ROOT + '/cgi-bin'
HTDOCS_REL_PATH     = '/manta'
CGI_BIN_REL_PATH    = '/cgi-bin/manta'
HTDOCS_ABS_PATH     = HTDOCS_ROOT + HTDOCS_REL_PATH
CGI_BIN_ABS_PATH    = WWW_ROOT + CGI_BIN_REL_PATH
RESULTS_ABS_PATH    = HTDOCS_ABS_PATH + '/results'
RESULTS_REL_PATH    = HTDOCS_REL_PATH + '/results'
HTDOCS_ABS_TEMPLATE_PATH = HTDOCS_ABS_PATH + '/templates'
JASPAR_URL          = 'http://jaspar.genereg.net'
REMOVE_RESULTFILES_OLDER_THAN   = 7 # days
#JASPAR_MATRIX_FILE  = '/devel/TFFM/data/jaspar_pfms.txt'
MANTA_DB_HOST       = 'localhost'
MANTA_DB_NAME       = 'manta'
MANTA_DB_USER       = 'manta_r'
MANTA_DB_PASS       = 'mantapw'
REF_ALLELE_REGEX    = re.compile('[ACGTNactgn]')
ALT_ALLELE_REGEX    = re.compile('[ACGTactg]')

JASPAR_DB_HOST      = 'vm5.cmmt.ubc.ca'
JASPAR_DB_NAME      = 'JASPAR_2014'
JASPAR_DB_USER      = 'jaspar_r'
JASPAR_DB_PASS      = ''
