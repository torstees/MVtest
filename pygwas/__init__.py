
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

"""


def sys_call(cmd):
    """Execute cmd and capture stdout and stderr

    :param cmd: command to be executed
    :return: (stdout, stderr)
    """
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    return p.stdout.readlines(), p.stderr.readlines()


def Exit(msg, code=1):
    """Exit execution with return code and message
    :param msg: Message displayed prior to exit
    :param code: code returned upon exiting
    """
    print >> sys.stderr, msg
    sys.exit(code)

def ExitIf(msg, do_exit, code=1):
    """Exit if do_exit is true

    :param msg: Message displayed prior to exit
    :param do_exit: exit when true
    :param code: application's return code upon exit
    """
    if do_exit:
        Exit(msg, code)

def BuildReportLine(key, value):
    """Prepare key/value for reporting in configuration report

    :param key: configuration 'keyword'
    :param value: value reported to be associated with keyword

    :return: formatted line starting with a comment
    """

    return "# " + key.ljust(20) + " : " + str(value)

