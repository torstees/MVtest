#!/usr/bin/env python

"""mv-many is a helper script designed to help the user split large mv-test jobs across many nodes with minimal
   effort. The approach is simple: The user provides most of the arguments required for running MVTest, along with
   some basic control information such as how many jobs to split the task into, etc and then the user provides a
   template which contains all of the job submission details necessary to run on their cluster. Example templates
   can be found at scripts/templates.

   This will generate the scripts. The user is required to submit those jobs and make sure that they completed
   successfully.
"""

__author__ = 'Eric Torstenson'
__version__ = 1.0

import argparse
from pygwas import sys_call
from pygwas import ExitIf
import os
import sys
import math
from string import Template

def check_and_append(args, flag, new_args):
    try:
        v = eval("args.%s" % (flag.replace("-", "_")))
        if v is not None:
            if v == True:
                new_args.append("--%s" % (flag))
            elif v != False:
                new_args.append("--%s %s" % (flag, str(v)))
    except:
        pass

# Example template:
"""
#SBATCH --job-name=$jobname
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
#SBATCH --time=3:00:00
#SBATCH --error $logpath/$jobname.e
#SBATCH --output $respath/$jobname.txt

$body
"""

def generate_jobs(args, job_list, argument_string):
    mvtest_path = args.mvpath
    template = "\n".join(args.template.readlines())
    logpath = os.path.abspath(args.logpath)
    respath = os.path.abspath(args.res_path)
    scriptpath = os.path.abspath(args.script_path)

    for jobname in job_list.keys():
        filename = "%s/%s.sh" % (scriptpath, jobname)
        job_body = mvtest_path + " " + job_list[jobname]
        contents = Template(template).safe_substitute(logpath=logpath, respath=respath, body=job_body, jobname=jobname)

        file = open(filename, "w")
        print >> file,contents
def mkdir(dirname):
    try:
        os.mkdir(dirname)
    except:
        pass

def split_mach_jobs(args, filename):
    """Parse the IMPUTE file and generate the list of jobs """
    max_snp_count = args.snps_per_job
    file_count = args.impute_count

    job_list = {}
    cur = None
    last_pos = None
    job_string = ""
    job_name = ""

    ExitIf("mv-many doesn't support splitting mach jobs into pieces at this time", args.max_snp_count > 0)

    file_count = sys_call("wc -l %s" % (filename))
    job_count = int(math.ceil(float(file_count) / args.mach_count))

    for job_idx in range(job_count):
        job_string = "--mach-offset %d" % (job_idx)
        job_list[str(job_idx)] = "--mach-offset %d" % (job_idx)


    return job_list

def split_impute_jobs(args, filename):
    """Parse the IMPUTE file and generate the list of jobs """
    max_snp_count = args.snps_per_job
    file_count = args.impute_count

    ExitIf("mv-many doesn't support splitting IMPUTE jobs into pieces at this time", args.max_snp_count > 0)

    job_list = {}

    file_count = sys_call("wc -l %s" % (filename))
    job_count = int(math.ceil(float(file_count) / args.mach_count))

    for job_idx in range(job_count):
        job_string = "--mach-offset %d" % (job_idx)
        job_list[str(job_idx)] = "--mach-offset %d" % (job_idx)

    return job_list

    # For now, let's not deal with the complexity of splitting chromosomes in IMPUTE
    poscol = 2
    cur = None
    last_pos = None
    job_string = ""
    job_name = ""

    file_index = 0
    for line in open(filename):
        chr, genfile = line.strip().split()

        if max_snp_count > 0:
            locus_index = 0
            last_pos = 1
            for locus in open(genfile):
                if locus_index >= max_snp_count - 1:
                    rsid, pos = locus.split()[1:2]
                    job_name = "chr%d_%d" % (chr, last_pos)
                    job_string = "--chr %s --from-bp %d --to-bp %d" % (chr, last_pos, pos)
                    last_pos = pos + 1
                    job_list[job_name] = job_string
                    locus_index = 0

                if cur is None:
                    cur = pos



    for line in sys_call("cut -f 1,%d %s" % (poscol, chrom_file)):
        chrom, pos = [int(x) for x in line.split()]
        if cur is None:     # First line observed
            cur = chrom
            job_string = "--chr %d --from-bp %d" % (chrom, pos)
            job_name = "Chr%d_%d-" % (chrom, pos)
            snp_count = 0
        elif cur != cur:    # Changed chromosome
            job_string += " --to-bp %d" % (last_pos)
            job_name += str(last_pos)
            job_list[job_name] = job_string
            cur = chrom
            job_string = "--chr %d --from-bp %d" % (chrom, pos)
            job_name = "Chr%d_%d-" % (chrom, pos)
            snp_count = 0
                            # create new job based on snp count
        elif snp_count < max_snp_count:
            snp_count += 1
        else:
            job_string += " --to-bp %d" % (last_pos)
            job_name += str(last_pos)
            job_list[job_name] = job_string
            job_string = "--chr %d --from-bp" % (chrom, pos)
            job_name = "Chr%d_%d-" % (chrom, pos)
            snp_count = 0

        last_pos = pos
    if job_string != "":
        job_string += " --to-bp %d" % (last_pos)
        job_name += str(last_pos)
        job_list[job_name] = job_string

    return job_list

def split_chrom_jobs(args, chrom_file):
    max_snp_count = args.snps_per_job

    poscol = 3
    if args.map3:
        poscol = 2

    job_list = {}
    cur = None
    last_pos = None
    job_string = ""
    job_name = ""

    for line in sys_call("cut -f 1,%d %s" % (poscol, chrom_file)):
        pos = -1

        values = line.split()
        if len(values) > 0:
            chrom, pos = [int(x) for x in values]

        if cur is None:     # First line observed
            cur = chrom
            job_string = "--chr %d --from-bp %d" % (chrom, pos)
            job_name = "Chr%d_%d-" % (chrom, pos)
            snp_count = 0
        elif cur != cur:    # Changed chromosome
            job_string += " --to-bp %d" % (last_pos)
            job_name += str(last_pos)
            job_list[job_name] = job_string
            cur = chrom
            job_string = "--chr %d --from-bp %d" % (chrom, pos)
            job_name = "Chr%d_%d-" % (chrom, pos)
            snp_count = 0
                            # create new job based on snp count
        elif snp_count < max_snp_count:
            snp_count += 1
        else:
            job_string += " --to-bp %d" % (last_pos)
            job_name += str(last_pos)
            job_list[job_name] = job_string
            job_string = "--chr %d --from-bp" % (chrom, pos)
            job_name = "Chr%d_%d-" % (chrom, pos)
            snp_count = 0

        last_pos = pos
    if job_string != "":
        job_string += " --to-bp %d" % (last_pos)
        job_name += str(last_pos)
        job_list[job_name] = job_string

    return job_list

def main(print_cfg):
    parser = argparse.ArgumentParser(description="mv-many.py -- MVTEST helper script.")
    parser.add_argument("-v", action="store_true", help="Print version number")

    parser.add_argument("--mvpath", type=str, default="mvtest.py", help="The path to mvtest.py if it's not found in PATH")
    parser.add_argument("--logpath", type=str, default="log", help="Where to write the logs to")
    parser.add_argument("--res-path", type=str, default="result", help="Where to write the results to")
    parser.add_argument("--script-path", type=str, default="scripts", help="Where to write the various scripts to")
    parser.add_argument("--template", type=argparse.FileType('r'), help="File containing cluster related information that will be written as header to the various files")
    parser.add_argument("--snps-per-job", type=float, default=1.0, help="Allow the user to .")

    # Most of these will just be passed on directly to the various MVTest scripts.
    parser.add_argument("--exclude", type=str, default="", help="Comma-delimited list of rsids to be excluded")

    # For now, I'm not implementing keep, since we don't have any real meaningful need for analyzing individuals
    # PLINK does, but we don't do the QC stuff they do.
    parser.add_argument("--keep", type=str, default="", help="Comma-delimited list of individuals to be analyzed")
    parser.add_argument("--remove", type=str, default="", help="Comma-delimited list of individuals to be removed from analysis")

    parser.add_argument("--file", type=str, help="Prefix for .ped and .map files")
    parser.add_argument("--ped", type=argparse.FileType('r'), help="PLINK compatible .ped file")
    parser.add_argument("--map", type=argparse.FileType('r'), help="PLINK compatible .map file")
    parser.add_argument("--map3", action='store_true', help="MAP file has only 3 columns")
    parser.add_argument("--no-sex", action='store_true', help="Pedigree file doesn't have column 5 (sex)")
    parser.add_argument("--no-parents", action="store_true", help="Pedigree file doesn't have columns 3 and 4 (parents)")
    parser.add_argument("--no-fid", action="store_true", help="Pedigree file doesn't have column 1 (family ID)")
    parser.add_argument("--no-pheno", action="store_true", help="Pedigree file doesn't have column 6 (phenotype")
    parser.add_argument("--liability", action="store_true", help="Pedigree file has column 7 (liability)")

    parser.add_argument("--bfile", type=str, help="Prefix for .bed, .bim and .fam files")
    parser.add_argument("--bed", type=argparse.FileType('r'), help="Binary Ped file (.bed)")
    parser.add_argument("--bim", type=argparse.FileType('r'), help="Binary ped marker file (.bim)")
    parser.add_argument("--fam", type=argparse.FileType('r'), help="Binary ped family file (.fam)")

    parser.add_argument("--tfile", type=str, help="Prefix for .tped and .tfam files")
    parser.add_argument("--tped", type=argparse.FileType('r'), help="Transposed Pedigree file (.tped)")
    parser.add_argument("--tfam", type=argparse.FileType('r'), help="Transposed pedigre Family file (.tfam)")
    parser.add_argument("--compressed", action="store_true", help="Ped/TPed compressed with gzip (named .ped.tgz or .tped.tgz)")


    parser.add_argument("--impute", type=argparse.FileType('r'), help="File containing list of impute output for analysis")
    parser.add_argument("--impute-fam", type=argparse.FileType('r'), help="File containing family details for impute data")
    parser.add_argument("--impute-uncompressed", action="store_true", help="Indicate that the impute input is not gzipped, but plain text")
    parser.add_argument("--impute-encoding", type=str, choices=['additive', 'dominant', 'recessive', 'genotype'], default='additive', help='Genetic model to be used')
    parser.add_argument("--impute-info-ext", type=str, default='info', help="Portion of filename denotes info filename")
    parser.add_argument("--impute-gen-ext", type=str, default='gen.gz', help="Portion of filename that denotes gen file")
    parser.add_argument("--impute-info-thresh", type=float, default=0.4, help="Threshold for filtering imputed SNPs with poor 'info' values")
    parser.add_argument("--impute-files-per-job", type=integer, default=1, help="How many gen files to run per job")

    parser.add_argument("--mach", type=argparse.FileType('r'), help="File containing list of MACH output for analysis")
    parser.add_argument("--mach-uncompressed", action="store_true", help="Indicate that the mach input is not gzipped")
    parser.add_argument("--mach-info-ext", type=str, default="info.gz", help="Portion of filename denotes info filenames")
    parser.add_argument("--mach-dose-ext", type=str, default="dose.gz", help="Portion of filename that denotes dose files")
    parser.add_argument("--mach-min-rsquared", type=float, default=0.3, help="Filter out loci with RSquared < this value")
    parser.add_argument("--mach-files-per-job", type=integer, default=1, help="How many dosage files to run per job")

    parser.add_argument("--pheno", type=argparse.FileType('r'), help="File containing phenotypes")
    parser.add_argument("--sample-pheno", type=argparse.FileType('r'), help="(Mach) Sample file containing phenotypes")
    parser.add_argument("--mphenos", type=str, default="", help="Column number(s) for phenotype to be analyzed if number of columns > 1")
    parser.add_argument("--pheno-names", type=str, default="", help="Name for phenotype(s) to be analyzed (must be in --pheno file)")
    parser.add_argument("--all-pheno", action="store_true", help="Analyze all columns from the phenotype file")

    parser.add_argument("--covar", type=argparse.FileType('r'), help="File containing covariates")
    parser.add_argument("--sample-covar", type=argparse.FileType('r'), help="(Mach) Sample file containing covariates")
    parser.add_argument("--covar-numbers", type=str, default="", help="Comma-separated list of covariate indices")
    parser.add_argument("--covar-names", type=str, default="", help="Comma-separated list of covariate names")
    parser.add_argument("--sex", action='store_true', help="Use sex from the pedigree file as a covariate")
    parser.add_argument("--missing-phenotype", type=float, default=-9.0, help="Encoding for missing phenotypes")

    parser.add_argument("--maf", type=float, default=0.0, help="Minimum MAF allowed for analysis")
    parser.add_argument("--max-maf", type=float, default=1.0, help="MAX MAF allowed for analysis")
    parser.add_argument("--geno", type=float, default=1.0, help="MAX per-SNP missing for analysis")
    parser.add_argument("--mind", type=float, default=1.0, help="MAX per-person missing")

    parser.set_defaults(all_pheno=False, sex=False)
    args = parser.parse_args(args)

    # Report version, if requested, and exit
    if args.v:
        print >> sys.stderr, "%s: %s" % (os.path.basename(__file__), __version__)
        sys.exit(0)
    mkdir(args.script_path)
    mkdir(args.res_path)
    mkdir(args.logpath)
    if args.template is None:
        args.template = open(os.path.dirname(os.path.realpath(__file__)) + "/template.txt")

    general_arguments = []
    for flag in "exclude,keep,remove,file,ped,map,map3,no-sex,no-parents,no-fid,no-pheno,liability," + \
                        "bfile,bed,bim,fam,tfile,tped,tfam,compressed," + \
                        "impute,impute-fam,impute-uncompresed,impute-encoding,impute-info-ext,impute-gen-ext,impute-info-thresh," +\
                        "mach,mach-uncompressed,mach-info-ext,mach-dose-ext,mach-min-rsquared," +\
                        "pheno,sample-pheno,pheno-names,mphenos,all-pheno," + \
                        "covar,sample-covar,covar-numbers,covar-names,sex,missing-phenotype,maf,max-maf,gen,mind".split(","):
        check_and_append(args, flag, general_arguments)


    job_list = None
    map_file = None
    impute_file_list = None
    if args.file:
        map_file = "%s.map" % (args.file)
    elif args.map:
        map_file = args.map
    elif args.bfile:
        map_file = "%s.bim" % (args.bfile)
    elif args.bim:
        map_file = args.bim
    elif args.tfile:
        map_file = "%s.tped"
    elif args.tped:
        map_file = args.tped
    elif args.impute:
        impute_file_list = args.impute

    max_snp_count = args.snps_per_job

    if map_file is not None:
        job_list = split_chrom_jobs(args, map_file)

    elif args.impute is not None:
        job_list = split_impute_jobs(args, args.impute)

    elif args.mach is not None:
        job_list = split_mach_jobs(args, args.mach)

    generate_jobs(args, job_list, " ".join(general_arguments))
    sys.exit(0)

if __name__ == "__main__":
    main(print_cfg=True)

