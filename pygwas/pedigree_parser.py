import gzip

import numpy

from pygwas.data_parser import DataParser
from pygwas.exceptions import MalformedInputFile
from pygwas.exceptions import TooManyAlleles
from pygwas import sys_call
from . import BuildReportLine


class Parser(DataParser):
    """Parse standard pedigree dataset

    Data should follow standard format for pedigree data, except alleles
    be either numerical (1 and 2) or as bases (A, C, T and G). All loci
    must have 2 alleles to be returned.

    Class Members:
    mapfile         Filename for the marker information
    datasource      Filename for the actual pedigree information
    genotypes       Matrix of genotype data
    invalid_loci    Loci that are being ignored due to filtration
    individual_mask Mask used to remove excluded and filtered calls
                    from the genotype data (each position represents
                    an individual)

    These are only available after load_genotypes() call
    locus_count     Number of valid loci
    markers         List of valid SNPs
    alleles         List of both alleles for each valid locus
    rsids           List of all SNP names for valid loci
    markers_maf     List of MAF at each locus




    """
    def __init__(self, mapfile, datasource):
        """Very basic initialization. """

        self.mapfile = mapfile
        self.datasource = datasource
        self.markers = []
        self.genotypes = []
        self.invalid_loci = []
        self.individual_mask = 0                         # Track which individuals were dropped due to missingness

    def ReportConfiguration(self, file):
        print >> file, BuildReportLine("PED FILE", self.datasource)
        print >> file, BuildReportLine("MAP FILE", self.mapfile)

    def load_mapfile(self, map3=False):
        """Load the marker data

        Builds up the marker list according to the boundary

        map3 indicate whether the file contains 3 or four columns.

        """
        global boundary
        cols            = [0, 1, 3]
        if map3:
            cols        = [0, 1, 2]
        markers         = numpy.loadtxt(self.mapfile, dtype=str, usecols=cols)

        self.snp_mask  = numpy.ones(markers.shape[0]*2, dtype=numpy.int8).reshape(-1, 2)

        if DataParser.boundary.NoExclusions():
            self.markers = numpy.zeros((markers.shape[0], 2), dtype=int)
            mask = markers[:, 2].astype(int) >= 0    # Check for plink's "off" mode
            self.snp_mask[:,0] = ~mask               # Turn them all 'on'
            self.snp_mask[:,1] = ~mask
            snp_count = numpy.sum(self.snp_mask[:, 0] == 0)

            self.markers[0:snp_count, 0] = markers[:, 0].astype(int)[mask]
            self.markers[0:snp_count, 1] = markers[:, 2].astype(int)[mask]
            self.rsids = markers[:, 1][mask]
            self.markers = self.markers[0:snp_count]
        else:
            idx = 0
            self.markers = []
            self.rsids   = []
            for locus in markers:
                if DataParser.boundary.TestBoundary(int(locus[0]), int(locus[2]), locus[1]):
                    self.markers.append([locus[0], locus[2]])
                    self.rsids.append(locus[1])
                    self.snp_mask[idx] = 0
                idx += 1
            self.markers = numpy.array(self.markers, dtype=numpy.int)
            self.rsids   = numpy.array(self.rsids)
        DataParser.boundary.beyond_upper_bound = False # We don't follow these rules here
        self.locus_count    = len(self.markers)


    def load_genotypes(self, pheno_covar):
        global has_sex, has_parents,has_pheno,has_fid, has_liability, ind_miss_tol
        global compressed_pedigree, ind_exclusions
        global missing_representation, snp_miss_tol, missing_storage
        first_genotype = 6
        pheno_col      = 5
        if not DataParser.has_sex:
            first_genotype -= 1
            pheno_col -= 1
        if not DataParser.has_parents:
            first_genotype -= 2
            pheno_col -= 2
        if not DataParser.has_pheno:
            first_genotype -= 1
        if not DataParser.has_fid:
            first_genotype -= 1
            pheno_col -= 1
        if DataParser.has_liability:
            first_genotype += 1

        sex_col = pheno_col - 1
        individual_mask = []
        self.individual_mask = []
        dropped_individuals = []

        # number of missing SNPs we can tolerate before dropping an individual
        max_missing_for_individual = numpy.sum(self.snp_mask[:, 0]==0) * DataParser.ind_miss_tol

        if DataParser.compressed_pedigree:
            ind_count, err = sys_call("gzip -cd %s | wc -l" % ("%s.gz" % (self.datasource)))
        else:
            ind_count, err = sys_call("wc -l %s" % (self.datasource))
        ind_count = int(ind_count[0].split()[0]) + 1

        snp_count = numpy.sum(self.snp_mask[:, 0] == 0)

        allelic_data = numpy.empty((ind_count, snp_count, 2), dtype='S1')

        valid_allele_count = 0
        if DataParser.compressed_pedigree:
            input_file = gzip.open("%s.gz" % self.datasource, 'rb')
        else:
            input_file = open(self.datasource)

        for line in input_file:
            line = line.strip()
            if len(line) > 0:
                raw_data = line.strip().split()
                alleles = numpy.ma.MaskedArray(numpy.array(raw_data[first_genotype:]).reshape(-1, 2), self.snp_mask).compressed().reshape(-1, 2)

                # Convert the alleles into genotypes

                indid = ":".join(raw_data[0:2])
                if not DataParser.has_fid:
                    indid = raw_data[0]

                # Ignore any subjects that are to be excluded and remove those that have too much missingness
                if DataParser.valid_indid(indid):
                    missing = numpy.sum(alleles[:, 0] == DataParser.missing_representation)

                    if missing > max_missing_for_individual:
                        individual_mask += [1, 1]
                        self.individual_mask.append(1)
                        dropped_individuals.append(indid)
                    else:
                        sex = None
                        phenotype = None
                        if DataParser.has_pheno:
                            phenotype = float(raw_data[pheno_col])
                        if DataParser.has_sex:
                            sex = int(raw_data[sex_col])
                        pheno_covar.add_subject(indid, sex, phenotype)
                        individual_mask += [0, 0]
                        self.individual_mask.append(0)
                        allelic_data[valid_allele_count] = alleles
                        valid_allele_count += 1

                else:
                    individual_mask += [1, 1]
                    self.individual_mask.append(1)
        self.ind_count = valid_allele_count
        allelic_data = allelic_data[0:valid_allele_count]
        self.genotypes = numpy.empty((snp_count, valid_allele_count))
        max_missing_individuals = DataParser.snp_miss_tol * ind_count
        dropped_loci = []
        valid_snps = 0
        valid_markers = []
        valid_rsids   = []
        valid_maf     = []
        valid_allele_list = []
        allele_count2s = []

        for i in xrange(0, snp_count):
            snp_geno = allelic_data[:,i]
            alleles = list(set(numpy.unique(snp_geno)) - set([DataParser.missing_representation]))

            if len(alleles) > 2:
                raise TooManyAlleles(chr=self.markers[i][0], rsid=self.rsids[i], alleles=alleles)
            #if len(alleles) == 1:
            #    print alleles, snp_geno, self.markers[i]
            #    raise TooFewAlleles(chr=self.markers[i][0], rsid=self.rsids[i], alleles=alleles)

            allele_count1 = numpy.sum(snp_geno==alleles[0])
            allele_count2 = 0
            maf = 0

            if len(alleles) > 1:
                allele_count2 = numpy.sum(snp_geno==alleles[1])
                real_allele_count2 = allele_count2

                if allele_count2 > allele_count1:
                    sorted_alleles = [alleles[1], alleles[0]]
                    alleles = sorted_alleles
                    allele_count = allele_count1
                    allele_count1 = allele_count2
                    allele_count2 = allele_count
                maf = allele_count2 / float(allele_count1 + allele_count2)
                allele_count2s.append(allele_count2)
                #genotypes = []
                major_allele       = alleles[0]
                minor_allele       = alleles[1]

                genotype_data = numpy.sum(snp_geno==alleles[1], axis=1)
                genotype_data[snp_geno[:, 0]==DataParser.missing_representation] = DataParser.missing_storage
            else:
                major_allele = alleles[0]
                minor_allele = '?'

            missing = numpy.sum(genotype_data==DataParser.missing_storage)
            if maf == 0 or maf < DataParser.min_maf or maf > DataParser.max_maf or max_missing_individuals < missing:
                locus_details = self.markers[i]
                DataParser.boundary.dropped_snps[locus_details[0]].add(locus_details[1])
                dropped_loci.append("%s:%s" % (locus_details[0], locus_details[1]))
                self.invalid_loci.append(i)
            else:
                self.genotypes[valid_snps, :] = genotype_data
                valid_snps += 1
                valid_markers.append(list(self.markers[i]))
                valid_rsids.append(self.rsids[i])
                valid_allele_list.append([major_allele, minor_allele])
                valid_maf.append(maf)

        self.markers = valid_markers
        self.alleles = valid_allele_list
        self.rsids   = valid_rsids
        self.locus_count = valid_snps
        self.genotypes = self.genotypes[0:self.locus_count, :]
        self.allele_count2s = allele_count2s

    def get_loci(self):
        return self.markers

    def populate_iteration(self, iteration):
        """Pour the current data into the iteration object"""

        cur_idx = iteration.cur_idx
        if cur_idx < self.locus_count:
            iteration.chr = self.markers[cur_idx][0]
            iteration.pos = self.markers[cur_idx][1]
            iteration.rsid = self.rsids[cur_idx]
            iteration.major_allele = self.alleles[cur_idx][0]
            iteration.minor_allele = self.alleles[cur_idx][1]
            iteration.allele_count2 = self.allele_count2s[cur_idx]
            iteration.genotype_data = self.genotypes[cur_idx, :]
            hetero = numpy.sum(iteration.genotype_data==1)
            iteration.min_allele_count = numpy.sum(iteration.genotype_data==2)*2 + hetero
            iteration.maj_allele_count = numpy.sum(iteration.genotype_data==0)*2 + hetero
            return True
        else:
            raise StopIteration

