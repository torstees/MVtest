What is MVTest?
===============

*TODO: Write some background information about the application and it's
scientific basis.*

Documentation
+++++++++++++
Documentation for mvtest is still under construction. However, the application
provides reasonable inline help using standard unix help arguments:

> `mvtest.py -h`

or

> `mvtest.py --help`

In general, overlapping functionality should mimic that of PLINK.

Due to automatic conversion of two dashes, "--", into an emdash (single long
dash) when writing to PDF, you may need to

Command-Line Arguments
++++++++++++++++++++++
Command line arguments used by mvtest often mimick those used by PLINK, except
where there is no matching functionality (or the functionality differs
significantly.)

For the parameters listed below, when a parameter requires a value, the value
must follow the argument with a single space separating the two (no '=' signs.)
For flags with no specified value, passing the flag indicates that condition
is to be "activated".

Getting help
------------

====================  =========  ============================================
 Flag(s)              Type       Description
====================  =========  ============================================
  -h, --help                     show this help message and exit
  -v                             Print version number
====================  =========  ============================================

.. raw:: latex

    \newpage

Input Data
----------
MVTest attempts to mimic the interface for PLINK where appropriate.

All input files should be whitespace delimited. For text based allelic
annotations, 1|2 and A|C|G|T annotation is sufficient. All data must
be expressed as alleles, not as genotypes (except for IMPUTE output,
which is a specialized format that is very different from the other
forms).

For Pedigree, Transposed Pedigree and PLINK binary pedigree files, the using
the prefix arguments is sufficient and recommended if your files follow
the standard naming conventions.

Pedigree Data
^^^^^^^^^^^^^
Pedigree data is fully supported, however it is not recommended. When loading
pedigree data, mvtest must load the entire dataset into memory prior to
analysis, which can result in a substantial amount of memory overhead that is
unnecessary.

Flags like --no-pheno and --no-sex can be used in any combination creating
MAP files with highly flexible header structures.

====================  ===========  ============================================
 Flag(s)              Type         Description
====================  ===========  ============================================
  --file FILE         file prefix  Prefix for .ped and .map files
  --ped PED           filename     PLINK compatible .ped file
  --map MAP           filename     PLINK compatible .map file
  --map3                           MAP file has only 3 columns
  --no-sex                         Pedigree file doesn't have column 5 (sex)
  --no-parents                     Pedigree file doesn't have columns 3 and 4 (parents)
  --no-fid                         Pedigree file doesn't have column 1 (family ID)
  --no-pheno                       Pedigree file doesn't have column 6 (phenotype
  --liability                      Pedigree file has column 7 (liability)
====================  ===========  ============================================


PLINK Binary Pedigree
^^^^^^^^^^^^^^^^^^^^^
This format represents the most efficient storage for large GWAS datasets,
and can be used directly by mvtest. In addition to a minimal overhead, plink
style bed files will also run very quickly, due to the efficient disk layout.

====================  ===========  ============================================
 Flag(s)              Type         Description
====================  ===========  ============================================
  --bfile FILE        file prefix  Prefix for .bed, .bim and .fam files
  --bed BED           filename     Binary Ped file (.bed)
  --bim MAP           filename     Binary ped marker file (.bim)
  --fam FAM           filename     Binary ped family file (.fam)
====================  ===========  ============================================

Transposed Pedigree Data
^^^^^^^^^^^^^^^^^^^^^^^^
Transposed Pedigree data is similar to standard pedigree except that the data
is arranged such that the data is organized as SNPs as rows, instead of
individuals. This allows mvtest to run it's analysis without loading the
entire dataset into memory.

====================  ===========  ============================================
 Flag(s)              Type         Description
====================  ===========  ============================================
  --tfile FILE        file prefix  Prefix for .tped and .tfam files
  --tped BED          filename     Transposed Pedigree file (.tped)
  --tfim MAP          filename     Transposed pedigree Family file (.tfam)
====================  ===========  ============================================

Pedigree/Transposed Pedigree Common Flags
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
By default, Pedigree and Transposed Pedigree data is assumed to be uncompressed.
However, mvtest can directly use gzipped data files if they have the extension
.tgz with the addition of the --compressed argument.

====================  ===========  ============================================
 Flag(s)              Type         Description
====================  ===========  ============================================
  --compressed        Ped/TPed     compressed with gzip (named .ped.tgz or \
                                   .tped.tgz)
====================  ===========  ============================================

IMPUTE output
^^^^^^^^^^^^^
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

===================================================  ===========  =================================
 Flag(s)                                             Type         Description
===================================================  ===========  =================================
  --impute IMPUTE                                    filename     File containing list of impute output for analysis
  --impute-fam IMPUTE_FAM                            filename     File containing family details for impute data
  --impute-offset IMPUTE_OFFSET                      int          Impute file index (1 based) to begin analysis
  --impute-count IMPUTE_COUNT                        int          Number of impute files to process (for this node). Defaults to all remaining.
  --impute-uncompressed                                           Indicate that the impute input is not gzipped, but plain text
  --impute-encoding {additive,dominant,recessive}    selection    Genetic model to be used when analyzing imputed data.
  --impute-info-ext IMPUTE_INFO_EXT                  file prefix  Portion of filename denotes info filename
  --impute-gen-ext IMPUTE_GEN_EXT                    file suffix  Portion of filename that denotes gen file
  --impute-info-thresh IMPUTE_INFO_THRESH            float        Threshold for filtering imputed SNPs with poor 'info' values
===================================================  ===========  =================================


.. raw:: latex

    \newpage

Phenotype/Covariate Data
^^^^^^^^^^^^^^^^^^^^^^^^
Phenotypes and Covariate data can be found inside either the standard pedigree headers or within special PLINK style
covariate files. Users can specify phenotypes and covariates using either header names (if a header exists in the file)
or by 1 based column indices.

========================================  ===========  =================================
 Flag(s)                                  Type         Description
========================================  ===========  =================================
  --pheno PHENO                           filename     File containing phenotypes. Unless --all-pheno is present, user must provide either index(s) or label(s) of the phenotypes to be analyzed.
  --mphenos MPHENOS                       numbers      Column number(s) for phenotype to be analyzed if number of columns > 1. Comma separated list if more than one is to be used.
  --pheno-names PHENO_NAMES               string       Name for phenotype(s) to be analyzed (must be in --pheno file). Comma separated list if more than one is to be used.
  --covar COVAR                           filename     File containing covariates
  --covar-numbers COVAR_NUMBERS           numbers      Comma-separated list of covariate indices
  --covar-names COVAR_NAMES                            Comma-separated list of covariate names
  --sex                                                Use sex from the pedigree file as a covariate
  --missing-phenotype MISSING_PHENOTYPE   character    Encoding for missing phenotypes as can be found in the data.
  --all-pheno                                          When present, mv-test will run each phenotypes found inside the phenotype file.
========================================  ===========  =================================


Restricting regions for analysis
--------------------------------
When specifying a range of positions for analysis, a chromosome must be present.
If a chromosome is specified but is not accompanied by a range, the entire
chromosome will be used. Only one range can be specified per run.

========================  ===========  =================================
 Flag(s)                  Type         Description
========================  ===========  =================================
  --snps SNPS             string       Comma-delimited list of SNP(s): rs1,rs2,rs3-rs6
  --chr N                 int          Select Chromosome. If not selected, all chromosomes are to be analyzed.
  --from-bp START         int          SNP range start
  --to-bp END             int          SNP range end
  --from-kb START         int          SNP range start
  --to-kb END             int          SNP range end
  --from-mb START         int          SNP range start
  --to-mb END             int          SNP range end
  --exclude EXCLUDE       string       Comma-delimited list of rsids to be excluded
  --remove REMOVE         string       Comma-delimited list of individuals to be removed from analysis. This must be in the form of family_id:individual_id
  --maf MAF               float        Minimum MAF allowed for analysis
  --max-maf MAX_MAF       float        MAX MAF allowed for analysis
  --geno GENO             int          MAX per-SNP missing for analysis
  --mind MIND             int          MAX per-person missing
  --verbose                            Output additional data details
========================  ===========  =================================
