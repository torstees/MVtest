# MVTEST

## What is MVTEST?
*TODO: Write some background information about the application and it's
scientific basis.*



## Installation
MVTest requires python 2.7 or later (but 3.0 and later) as well as the
following packages:

* NumPy (version 1.7.2 or later)   www.numpy.org
* SciPY (version 0.13.2 or later)  www.scipy.org

If these aren't already installed, and you don't have root access to the
machine, please see the section, Miniconda for easy instructions on
installing a compact version of python suitable for running MVTest. MVTest's
installer can install these for you, however, it assumes that you have
write access to your python library, which will not be the case on shared
systems.

Download the package at: TODO: URL

To install the software, run the setup script as shown below:
> `python setup.py install`

If no errors are reported, it should be installed and ready to use.

## System Requirements
Aside from the library dependencies, MVTest's requirements depend largely on
the number of SNPs and individuals being analyzed as well as the data format
being used. In general, GWAS sized datasets will require several gigabytes of
memory when using the traditional pedigree format, however, even 10s of
thousands of subjects can be analyzed with less than 1 gigabyte of RAM when
the data is formatted as transposed pedigree or PLINK's default bed format.

Otherwise, it is recommended that the system be run on a unix-like system
such as Linux or OS X, but it should work under windows as well (this is
not a fully supported platform).


## Running Unit Tests
MVTest comes with a unit test suite which can be run prior to installation.
To run the tests, simply run the following command from within the root
directory of the extracted archive's contents:
python setup.py test

If no errors are reported, then mvtest should run correctly on your system.


## Miniconda
Miniconda is a minimal version of the package manager used by the Anaconda
python distribution. It makes it easy to create local installations of python
with the latest versions of the common scientific libraries for users who don't
have root access to their target machines.

Firstly, download and install the appropriate version of miniconda at the
project website. Please be sure to choose the Python 2 version:
http://conda.pydata.org/miniconda.html

While it is doing the installation, please allow it to update your PATH
information. Also, be sure to follow directions such as starting a new
shell to allow those changes to take effect.

Once those changes have taken effect, install setuptools and scipy:
conda install setuptools scipy

Installing SciPy will also force the installation of NumPy, whish is
also required for running mvtest. (setuptools includes easy_install).

Once that has been completed successfully, you should be ready to follow
the standard instructions for installing mvtest.



# Documentation
Documentation for mvtest is still under construction. However, the application
provides reasonable inline help using standard unix help arguments:

> `mvtest.py -h`

or

> `mvtest.py --help`

In general, overlapping functionality should mimic that of PLINK.

## Command-Line Arguments

Command line arguments used by mvtest often mimick those used by PLINK, except
where there is no matching functionality (or the functionality differs
significantly.)

For the parameters listed below, when a parameter requires a value, the value
must follow the argument with a single space separating the two (no '=' signs.)
For flags with no specified value, passing the flag indicates that condition
is to be "activated".

### Getting help
Flag(s)          | Type     | Description
---------------- | -------- | -----------
  -h, --help     |          | show this help message and exit
  -v             |          | Print version number

### Input Data
MVTest attempts to mimic the interface for PLINK where appropriate.

All input files should be whitespace delimited. For text based allelic
annotations, 1|2 and A|C|G|T annotation is sufficient. All data must
be expressed as alleles, not as genotypes (except for IMPUTE output,
which is a specialized format that is very different from the other
forms).

For Pedigree, Transposed Pedigree and PLINK binary pedigree files, the using
the prefix arguments is sufficient and recommended if your files follow
the standard naming conventions.

#### Pedigree Data
Pedigree data is fully supported, however it is not recommended. When loading
pedigree data, mvtest must load the entire dataset into memory prior to
analysis, which can result in a substantial amount of memory overhead that is
unnecessary.

Flags like --no-pheno and --no-sex can be used in any combination creating
MAP files with highly flexible header structures.

Flag(s)          | Type        | Description
---------------- | ----------- | -----------
  --file FILE    | file prefix | Prefix for .ped and .map files
  --ped PED      | filename    | PLINK compatible .ped file
  --map MAP      | filename    | PLINK compatible .map file
  --map3         |             | MAP file has only 3 columns
  --no-sex       |             | Pedigree file doesn't have column 5 (sex)
  --no-parents   |             | Pedigree file doesn't have columns 3 and 4 (parents)
  --no-fid       |             | Pedigree file doesn't have column 1 (family ID)
  --no-pheno     |             | Pedigree file doesn't have column 6 (phenotype
  --liability    |             | Pedigree file has column 7 (liability)

#### PLINK Binary Pedigree
This format represents the most efficient storage for large GWAS datasets,
and can be used directly by mvtest. In addition to a minimal overhead, plink
style bed files will also run very quickly, due to the efficient disk layout.

Flag(s)          | Type        | Description
---------------- | ----------- | -----------
  --bfile BFILE  | file prefix |  Prefix for .bed, .bim and .fam files
  --bed BED      | filename    | Binary Ped file (.bed)
  --bim BIM      | filename    | Binary ped marker file (.bim)
  --fam FAM      | filename    | Binary ped family file (.fam)

#### Transposed Pedigree Data
Transposed Pedigree data is similar to standard pedigree except that the data
is arranged such that the data is organized as SNPs as rows, instead of
individuals. This allows mvtest to run it's analysis without loading the
entire dataset into memory.

Flag(s)          | Type        | Description
---------------- | ----------- | -----------
  --tfile TFILE  | file prefix | Prefix for .tped and .tfam files
  --tped TPED    | filename    | Transposed Pedigree file (.tped)
  --tfam TFAM    | filename    | Transposed pedigree Family file (.tfam)

#### Pedigree/Transposed Pedigree Common Flags
By default, Pedigree and Transposed Pedigree data is assumed to be uncompressed.
However, mvtest can directly use gzipped data files if they have the extension
.tgz with the addition of the --compressed argument.

Flag(s)        | Description
-------------- | -----------
  --compressed | Ped/TPed compressed with gzip (named .ped.tgz or .tped.tgz)

#### IMPUTE output
MVTest doesn't call genotypes when performing analysis, and allows users to
define which model to use when analyzing the data. Due to the fact that there
is no specific location for chromosome within the input files, mvtest requires
that users provide chromosome, impute input file and the corresponding .info
file for each imputed output.

Due to the huge number of expected loci, mvtest allows users to specify an
offset and file count for analysis. This is to allow users to run multiple
jobs simultaneously on a cluster and work individually on separate impute
region files. Users can segment those regions even further using standard
mvtest region selection as well.

By default, all imputed data is assumed to be compressed using gzip.

Default naming convention is for impute data files to end in .gen.gz and
the info files to have the same name except for the end being replaced
by .info.

Flag(s)                         | Type     | Default | Description
-----------------               | -------- | --------| -----------
  --impute IMPUTE               | filename |         | File containing list of impute output for analysis
  --impute-fam IMPUTE_FAM       | filename |         | File containing family details for impute data
  --impute-offset IMPUTE_OFFSET | int      | -1      | Impute file index (1 based) to begin analysis
  --impute-count IMPUTE_COUNT   | int      | -1      | Number of impute files to process (for this node). Defaults to all remaining.
  --impute-uncompressed         |          |         | Indicate that the impute input is not gzipped, but plain text
  --impute-encoding {additive,dominant,recessive} | one of the three options | additive | Genetic model to be used when analyzing imputed data.
  --impute-info-ext IMPUTE_INFO_EXT | file prefix | info | Portion of filename denotes info filename
  --impute-gen-ext IMPUTE_GEN_EXT | file suffix | gen.gz | Portion of filename that denotes gen file
  --impute-info-thresh IMPUTE_INFO_THRESH | floating point value | 0.4 | Threshold for filtering imputed SNPs with poor 'info' values


#### Phenotype/Covariate Data
Phenotypes and Covariate data can be found inside either the standard pedigree headers or within special PLINK style
covariate files. Users can specify phenotypes and covariates using either header names (if a header exists in the file)
or by 1 based column indices.

Flag(s)           | Type     | Description
----------------- | -------- | -----------
  --pheno PHENO   | filename | File containing phenotypes
  --mphenos MPHENOS | numbers| Column number(s) for phenotype to be analyzed if number of columns > 1. Comma separated list if more than one is to be used.
  --pheno-names PHENO_NAMES  | string | Name for phenotype(s) to be analyzed (must be in --pheno file). Comma separated list if more than one is to be used.
  --covar COVAR  | filename  | File containing covariates
  --covar-numbers COVAR_NUMBERS | numbers | Comma-separated list of covariate indices
  --covar-names COVAR_NAMES | | Comma-separated list of covariate names
  --sex          |          | Use sex from the pedigree file as a covariate
  --missing-phenotype MISSING_PHENOTYPE | character string | Encoding for missing phenotypes as can be found in the data.


### Restricting regions for analysis
When specifying a range of positions for analysis, a chromosome must be present.
If a chromosome is specified but is not accompanied by a range, the entire
chromosome will be used. Only one range can be specified per run.

Flag(s)           | Type     | Description
----------------- | -------- | -----------
  --snps SNPS     | string   | Comma-delimited list of SNP(s): rs1,rs2,rs3-rs6
  --chr N         | int      | Select Chromosome. If not selected, all chromosomes are to be analyzed.
  --from-bp START | int      | SNP range start
  --to-bp END     | int      | SNP range end
  --from-kb START | int      | SNP range start
  --to-kb END     | int      | SNP range end
  --from-mb START | int      | SNP range start
  --to-mb END     | int      | SNP range end
  --exclude EXCLUDE | string | Comma-delimited list of rsids to be excluded
  --remove REMOVE   | string | Comma-delimited list of individuals to be removed from analysis. This must be in the form of family_id:individual_id
  --maf MAF      | float    | Minimum MAF allowed for analysis
  --max-maf MAX_MAF | float | MAX MAF allowed for analysis
  --geno GENO    | int      | MAX per-SNP missing for analysis
  --mind MIND    | int      | MAX per-person missing
  --verbose      |          | Output additional data details


#Change History
* Version 1.2.1
    * Updated destandardization function to correct the standard error for variances after Chun studied the formula more carefully.
* Version 1.2.0
    * Initial modularization of data parsers to allow sharing with other gwas programs
    * Switched to setuptools for easier installation.
    * Bug fixes for IMPUTED data and improved unit testing integration.
* Version 1.1.1
    * Added support for negative positions to be interpreted as missing in order to remain consistent with PLINK's behavior.
    * Revised genotype parsing for pedigree data to ensure usability with large GWAS datasets
    * Added looping feature to Mean Var test with different starting values in order to help with convergence when using with large numbers of covariates.
    * Enabled gzip compression for standard and Transposed Pedigree pedigree data. (Bed files are not affected).