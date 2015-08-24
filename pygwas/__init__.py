
__author__ = 'Eric Torstenson'
__version__ = '0.5'

import subprocess
import sys
import exceptions

"""PYthon GWAS (pygwas) library

This library contains all of the classes required to parse standard and
transposed pedigree data using most of PLINK style enhancements (such as
tolerating bases and missing/additional columns). In addition to standard
text pedigrees, there is also support for bed files.

Each parser class will be able to parse the various sister files that are
necessary for complete data representation.


Threshold variables represent values which, when surpassed, will
trigger an exclusion (either SNP or individual)

Module Data Members:
missing_representation  Input file representation of missing genotypes
missing_storage Internal representation of missing genotypes
min_maf         Minimum threshold for minor allele frequency
max_maf         Maximum threshold for minor allele frequency
snp_miss_tol    Minimum threshold for perc. of missing ind
ind_miss_tol    Minimum threshold for perc. of missing loci
boundary        User configuration for genomic traversal
ind_exclusions  List of IDs to be ignored
has_sex         Std Pedigree col 5 is present
has_parents     Std pedigree cols 3 and 4 are present
has_fid         Std pedigree col 1 is present
has_pheno       Std pedigree col 6 is present
has_liability   Non std pedigree col 7 is present

"""


def sys_call(cmd):
    """Execute cmd and capture stdout and stderr"""
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    return p.stdout.readlines(), p.stderr.readlines()


def Exit(msg, code=1):
    """Exit execution with return code and message"""
    print >> sys.stderr, msg
    sys.exit(code)

def ExitIf(msg, do_exit, code=1):
    if do_exit:
        Exit(msg, code)

def BuildReportLine(key, value):
    """Prepare key/value for reporting in configuration report"""

    return "# " + key.ljust(20) + " : " + str(value)

