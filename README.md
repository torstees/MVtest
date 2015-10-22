What is MVTest?
===============

*TODO: Write some background information about the application and it's scientific basis.*

Documentation
-------------

Documentation for mvtest is still under construction. However, the application provides reasonable inline help using standard unix help arguments:

\> mvtest.py -h

or

\> mvtest.py --help

In general, overlapping functionality should mimic that of PLINK.

Command-Line Arguments
----------------------

Command line arguments used by mvtest often mimick those used by PLINK, except where there is no matching functionality (or the functionality differs significantly.)

For the parameters listed below, when a parameter requires a value, the value must follow the argument with a single space separating the two (no '=' signs.) For flags with no specified value, passing the flag indicates that condition is to be "activated".

When there is no value listed in the "Type" column, the arguments are *off* by default and *on* when the argument is present (i.e. by default, compression is turned off except when the flag, --compression, has been provided.)

### Getting help

| Flag(s)    | Description                     |
|------------|---------------------------------|
| -h, --help | show this help message and exit |
| -v         | Print version number            |

\newpage
### Input Data

MVTest attempts to mimic the interface for PLINK where appropriate.

All input files should be whitespace delimited. For text based allelic annotations, 1|2 and A|C|G|T annotation is sufficient. All data must be expressed as alleles, not as genotypes (except for IMPUTE output, which is a specialized format that is very different from the other forms).

For Pedigree, Transposed Pedigree and PLINK binary pedigree files, the using the prefix arguments is sufficient and recommended if your files follow the standard naming conventions.

#### Pedigree Data

Pedigree data is fully supported, however it is not recommended. When loading pedigree data, mvtest must load the entire dataset into memory prior to analysis, which can result in a substantial amount of memory overhead that is unnecessary.

Flags like --no-pheno and --no-sex can be used in any combination creating MAP files with highly flexible header structures.

| Flag(s)      | Type        | Description                                          |
|--------------|-------------|------------------------------------------------------|
| --file FILE  | file prefix | Prefix for .ped and .map files                       |
| --ped PED    | filename    | PLINK compatible .ped file                           |
| --map MAP    | filename    | PLINK compatible .map file                           |
| --map3       |             | MAP file has only 3 columns                          |
| --no-sex     |             | Pedigree file doesn't have column 5 (sex)            |
| --no-parents |             | Pedigree file doesn't have columns 3 and 4 (parents) |
| --no-fid     |             | Pedigree file doesn't have column 1 (family ID)      |
| --no-pheno   |             | Pedigree file doesn't have column 6 (phenotype       |
| --liability  |             | Pedigree file has column 7 (liability)               |

#### PLINK Binary Pedigree

This format represents the most efficient storage for large GWAS datasets, and can be used directly by mvtest. In addition to a minimal overhead, plink style bed files will also run very quickly, due to the efficient disk layout.

| Flag(s)      | Type        | Description                          |
|--------------|-------------|--------------------------------------|
| --bfile FILE | file prefix | Prefix for .bed, .bim and .fam files |
| --bed BED    | filename    | Binary Ped file (.bed)               |
| --bim MAP    | filename    | Binary ped marker file (.bim)        |
| --fam FAM    | filename    | Binary ped family file (.fam)        |

#### Transposed Pedigree Data

Transposed Pedigree data is similar to standard pedigree except that the data is arranged such that the data is organized as SNPs as rows, instead of individuals. This allows mvtest to run it's analysis without loading the entire dataset into memory.

| Flag(s)      | Type        | Description                             |
|--------------|-------------|-----------------------------------------|
| --tfile FILE | file prefix | Prefix for .tped and .tfam files        |
| --tped BED   | filename    | Transposed Pedigree file (.tped)        |
| --tfim MAP   | filename    | Transposed pedigree Family file (.tfam) |

#### Pedigree/Transposed Pedigree Common Flags

By default, Pedigree and Transposed Pedigree data is assumed to be uncompressed. However, mvtest can directly use gzipped data files if they have the extension .tgz with the addition of the --compressed argument.

| Flag(s)      | Type     | Description                             |
|--------------|----------|-----------------------------------------|
| --compressed | Ped/TPed | compressed with gzip (named .ped.tgz or |
|              |          | .tped.tgz)                              |

#### IMPUTE output

MVTest doesn't call genotypes when performing analysis and allows users to define which model to use when analyzing the data. Due to the fact that there is no specific location for chromosome within the input files, mvtest requires that users provide chromosome, impute input file and the corresponding .info file for each imputed output.

Due to the huge number of expected loci, mvtest allows users to specify an offset and file count for analysis. This is to allow users to run multiple jobs simultaneously on a cluster and work individually on separate impute region files. Users can segment those regions even further using standard mvtest region selection as well.

By default, all imputed data is assumed to be compressed using gzip.

Default naming convention is for impute data files to end in .gen.gz and the info files to have the same name except for the end being replaced by .info.

| Flag(s)                                         | Type        | Description                                                                   |
|-------------------------------------------------|-------------|-------------------------------------------------------------------------------|
| --impute IMPUTE                                 | filename    | File containing list of impute output for analysis                            |
| --impute-fam IMPUTE\_FAM                        | filename    | File containing family details for impute data                                |
| --impute-offset IMPUTE\_OFFSET                  | int         | Impute file index (1 based) to begin analysis                                 |
| --impute-count IMPUTE\_COUNT                    | int         | Number of impute files to process (for this node). Defaults to all remaining. |
| --impute-uncompressed                           |             | Indicate that the impute input is not gzipped, but plain text                 |
| --impute-encoding {additive,dominant,recessive} | selection   | Genetic model to be used when analyzing imputed data.                         |
| --impute-info-ext IMPUTE\_INFO\_EXT             | file prefix | Portion of filename denotes info filename                                     |
| --impute-gen-ext IMPUTE\_GEN\_EXT               | file suffix | Portion of filename that denotes gen file                                     |
| --impute-info-thresh IMPUTE\_INFO\_THRESH       | float       | Threshold for filtering imputed SNPs with poor 'info' values                  |

#### IMPUTE File Input

When performing an analysis on IMPUTE output, users must provide a single file which lists each of the gen files to be analyzed. This plain text file contains 2 (or optionally 3) columns for each gen file:

| Col 1 (chromosome) | Col 2 (gen file) | Col 3 (optional .info filename) |
|--------------------|------------------|---------------------------------|
| N (chromosome \#)  | filename         | filename                        |
| ...                | ...              | ...                             |

The 3rd column is only required if your .info files and .gen files are not the same except for the extension.

#### MACH output

Users can analyze data imputed with MACH. Because most situations require many files, the format is a single file which contains either pairs of dosage/info files, or, if the two files share the same filename except for extensions, one dosage file per line.

There is one caveat when using MACH output for analysis: MV-Test requires Chromosome and Position for consistency in reporting. As such, the IDs inside .info files must be of the form: chrom:pos

If RSIDs or solely positions are found, MVTest will exit with an error.

When running mvtest using MACH dosage on a cluster, users can instruct a given job to anlyze data from a portion of the files contained within the MACH dosage file list by changing the --mach-offset and --mach-count arguments. By default, the offset starts with 1 (the first file in the dosage list) and runs all it finds. However, if one were to want to split the jobs up to analyze three dosage files per job, they might set those values to --mach-offset 1 --mach-count 3 or --mach-offset 4 --mach-count 3 depending on which job is being defined.

In order to minimize memory requirements, MACH dosage files can be loaded incrementally such that only N loci are stored in memory at a time. This can be controlled using the --mach-chunk-size argument. The larger this number is, the faster MVTest will run (fewer times reading from file) but the more memory is required.

| Flag(s)                       | Type     | Description                                                                                                                                                                         |
|-------------------------------|----------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| --mach MACH                   | filename | File containing list of dosages, one per line. Optionally, lines may contain the info names as well (separated by whitespace) if the two filenames do not share a common base name. |
| --mach-offset OFFSET          | number   | Index into the MACH file to begin analyzing                                                                                                                                         |
| --mach-count COUNT            | number   | Number of dosage files to analyze                                                                                                                                                   |
| --mach-uncompressed           |          | By default, MACH input is expected to be gzip compressed. If data is plain text, add this flag                                                                                      |
| --mach-chunk-size CHUNK\_SIZE | number   | Due to the individual orientation of the data, large dosage files are parsed in chunks in order to minimize excessive memory during loading                                         |
| --mach-info-ext EXT           | string   | Indicate the extension used by the mach info files                                                                                                                                  |
| --mach-dose-ext EXT           | string   | Indicate the extension used by the mach dosage files                                                                                                                                |
| --mach-min-rsquared MIN       | float    | Indicate the minimum threshold for the rsqured value from the .info files required for analysis.                                                                                    |

#### MACH File Input

When running an analysis on MACH output, users must provide a single file which lists of each dosage file and (optionally) the matching .info file. This file is a simple text file with either 1 column (the dosage filename) or 2 (dosage filename followed by the info filename separated by whitespace).

The 2nd column is only required if the filenames aren't identical except for the extension.

| Col 1 (dosage filename) | Col 2 (optional info filename) |
|-------------------------|--------------------------------|
| filename.dose           | filename.info                  |
| ...                     | ...                            |

\newpage
#### Phenotype/Covariate Data

Phenotypes and Covariate data can be found inside either the standard pedigree headers or within special PLINK style covariate files. Users can specify phenotypes and covariates using either header names (if a header exists in the file) or by 1 based column indices.

| Flag(s)                                | Type      | Description                                                                                                                                |
|----------------------------------------|-----------|--------------------------------------------------------------------------------------------------------------------------------------------|
| --pheno PHENO                          | filename  | File containing phenotypes. Unless --all-pheno is present, user must provide either index(s) or label(s) of the phenotypes to be analyzed. |
| --mphenos MPHENOS                      | numbers   | Column number(s) for phenotype to be analyzed if number of columns \> 1. Comma separated list if more than one is to be used.              |
| --pheno-names PHENO\_NAMES             | string    | Name for phenotype(s) to be analyzed (must be in --pheno file). Comma separated list if more than one is to be used.                       |
| --covar COVAR                          | filename  | File containing covariates                                                                                                                 |
| --covar-numbers COVAR\_NUMBERS         | numbers   | Comma-separated list of covariate indices                                                                                                  |
| --covar-names COVAR\_NAMES             |           | Comma-separated list of covariate names                                                                                                    |
| --sex                                  |           | Use sex from the pedigree file as a covariate                                                                                              |
| --missing-phenotype MISSING\_PHENOTYPE | character | Encoding for missing phenotypes as can be found in the data.                                                                               |
| --all-pheno                            |           | When present, mv-test will run each phenotypes found inside the phenotype file.                                                            |

### Restricting regions for analysis

When specifying a range of positions for analysis, a chromosome must be present. If a chromosome is specified but is not accompanied by a range, the entire chromosome will be used. Only one range can be specified per run.

| Flag(s)            | Type   | Description                                                                                                            |
|--------------------|--------|------------------------------------------------------------------------------------------------------------------------|
| --snps SNPS        | string | Comma-delimited list of SNP(s): rs1,rs2,rs3-rs6                                                                        |
| --chr N            | int    | Select Chromosome. If not selected, all chromosomes are to be analyzed.                                                |
| --from-bp START    | int    | SNP range start                                                                                                        |
| --to-bp END        | int    | SNP range end                                                                                                          |
| --from-kb START    | int    | SNP range start                                                                                                        |
| --to-kb END        | int    | SNP range end                                                                                                          |
| --from-mb START    | int    | SNP range start                                                                                                        |
| --to-mb END        | int    | SNP range end                                                                                                          |
| --exclude EXCLUDE  | string | Comma-delimited list of rsids to be excluded                                                                           |
| --remove REMOVE    | string | Comma-delimited list of individuals to be removed from analysis. This must be in the form of family\_id:individual\_id |
| --maf MAF          | float  | Minimum MAF allowed for analysis                                                                                       |
| --max-maf MAX\_MAF | float  | MAX MAF allowed for analysis                                                                                           |
| --geno GENO        | int    | MAX per-SNP missing for analysis                                                                                       |
| --mind MIND        | int    | MAX per-person missing                                                                                                 |
| --verbose          |        | Output additional data details                                                                                         |

Installation
============

MVTest requires python 2.7.x as well as the following libraries:

-   NumPy (version 1.7.2 or later) www.numpy.org
-   SciPY (version 0.13.2 or later) www.scipy.org

MVTest's installation will attempt to install these required components for you, however, it requires that you have write permission to the installation directory. If you are using a shared system and lack the necessary privileges to install libraries and software yourself, you should please see one of the sections, Miniconda or virtual-env for instructions on different options for setting up your own python environement which will exist entirely under your own control.

Download the package at: TODO: URL

To install the software, run the setup script as shown below: $ python setup.py install

If no errors are reported, it should be installed and ready to use.

**Regarding PYTHON 3** I began the process of updating the code to work with both python versions 2 and 3, however, there are some real issues with some library support of version 3 that is discouraging. So, until those have been resolved, I have no plans to invest further time toward support for python 3.

System Requirements
-------------------

Aside from the library dependencies, MVTest's requirements depend largely on the number of SNPs and individuals being analyzed as well as the data format being used. In general, GWAS sized datasets will require several gigabytes of memory when using the traditional pedigree format, however, even 10s of thousands of subjects can be analyzed with less than 1 gigabyte of RAM when the data is formatted as transposed pedigree or PLINK's default bed format.

Otherwise, it is recommended that the system be run on a unix-like system such as Linux or OS X, but it should work under windows as well (we can't offer support for running MVTest under windows).

Running Unit Tests
------------------

MVTest comes with a unit test suite which can be run prior to installation. To run the tests, simply run the following command from within the root directory of the extracted archive's contents:

$ python setup.py test

If no errors are reported, then mvtest should run correctly on your system.

Virtual Env
-----------

Virtual ENV is a powerful too for python programmers and end users alike as it allows for users to deploy different versions of python applications without the need for root access to the machine.

Because MVTest requires version 2.7, you'll need to ensure that your machine's python version is in compliance. Virtual Env basically uses the the system version of python, but creates a user owned environment wrapper allowing users to install libraries easily without administrative rights to the machine.

For a helpful introduction to VirtualEnv, please have a look at the tutorial: <http://www.simononsoftware.com/virtualenv-tutorial/>

Miniconda
---------

Miniconda is a minimal version of the package manager used by the Anaconda python distribution. It makes it easy to create local installations of python with the latest versions of the common scientific libraries for users who don't have root access to their target machines. Basically, when you use miniconda, you'll be installing your own version of Python into a directory under your control which allows you to install anything else you need without having to submit a helpdesk ticket for administrative assistance.

Unlike pip, the folks behind the conda distributions provide binary downloads of it's selected library components. As such, only the most popular libraries, such as pip, NumPY and SciPy, are supported by conda itself. However, these do not require compilation and may be easier to get installed than when using pip alone. I have experienced difficulty installing SciPy through pip and setup tools on our cluster here at vanderbilt due to non-standard paths for certain required components, but mini-conda always comes through.

Firstly, download and install the appropriate version of miniconda at the project website. Please be sure to choose the Python 2 version: <http://conda.pydata.org/miniconda.html>

While it is doing the installation, please allow it to update your PATH information. If you prefer not to always use this version of python in the future, simple tell it not to update your .bashrc file and note the instructions for loading and unloading your new python environment. Please note that even if you chose to update your .bashrc file, you will need to follow directions for loading the changes into your current shell.

Once those changes have taken effect, install setuptools and scipy: $ conda install pip scipy

Installing SciPy will also force the installation of NumPy, which is also required for running mvtest. (setuptools includes easy\_install).

Once that has been completed successfully, you should be ready to follow the standard instructions for installing mvtest. :tocdepth: 2

Change Log
==========
