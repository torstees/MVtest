
class Locus(object):
    def __init__(self, other=None):
        if other:
            self.chr = other.chr
            self.pos = other.pos
            self.rsid = other.rsid
            self.alleles = list(other.alleles)
            self.hetero_count = other.hetero_count
            self.min_allele_count = other.min_allele_count
            self.maj_allele_count = other.maj_allele_count
            self.missing_allele_count = other.missing_allele_count
            self._maf = other._maf
        else:
            self._maf = None
            self.chr = -1
            self.pos = -1
            self.rsid = ""
            self.alleles = ["",""]
            self.hetero_count = -1
            self.min_allele_count = -1
            self.maj_allele_count = -1
            self.missing_allele_count = -1

    def flip(self):
        """This will switch major/minor around, regardless of frequency truth.

        This is intended for forcing one of two populations to relate correctly to the same genotype definitions. When
        flipped, Ps and Qs will be backward, and the maf will no longer relate to the "minor" allele frequency. However,
        it does allow clients to use the same calls for each population without having to perform checks during
        those calculations. """

        maj_count = self.maj_allele_count
        self.maj_allele_count = self.min_allele_count
        self.min_allele_count = maj_count

        alleles = self.alleles
        self.alleles = [alleles[1], alleles[0]]

    @property
    def sample_size(self):
        return self.total_allele_count * 0.5

    @property
    def total_allele_count(self):
        return self.min_allele_count + self.maj_allele_count

    @property
    def hetero_freq(self):
        return self.hetero_count / (0.5 * (self.min_allele_count + self.maj_allele_count))
    @property
    def exp_hetero_freq(self):
        return 2 * self.p * self.q

    @property
    def p(self):
        return self.maj_allele_count / float(self.total_allele_count)

    @property
    def q(self):
        return self.min_allele_count / float(self.total_allele_count)

    @property
    def major_allele(self):
        return self.alleles[0]
    @major_allele.setter
    def major_allele(self, allele):
        self.alleles[0] = allele


    @property
    def minor_allele(self):
        return self.alleles[1]
    @minor_allele.setter
    def minor_allele(self, allele):
        self.alleles[1] = allele

    @property
    def maf(self):
        if self._maf:
            return self._maf
        return self.q

    def __cmp__(self, other):
        if self.chr == other.chr:
            return self.pos.__cmp__(other.pos)
        return self.chr.__cmp__(other.chr)

    def __str__(self):
        return "%d\t%d:%d %s %s %s %0.4f %0.4f %s" % (self.cur_idx, self.chr, self.pos, self.rsid, self.major_allele, self.minor_allele, self.maf, self.hetero_freq, self.genotype_data)
