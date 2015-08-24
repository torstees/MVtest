from locus import Locus

class ParsedLocus(Locus):
    """Locus data representing current iteration from a dataset

    Provide an iterator interface for all dataset types.

    Class Members:
        * chr             Chromosome of current locus
        * pos             BP offset from beginning of current locus ( >=1 )
        * rsid            SNP name of current locus
        * major_allele    Letter representation of the most common allele
        * minor_allele    Representation of the less common allele
        * maf             Frequency of the minor_allele
        * genotype_data   array of genotypes. Genotypes are represented as \
                            the number of minor alleles

    """

    def __init__(self, datasource, index=-1):
        """Basic initialization (nothing is currently valid)"""

        super(ParsedLocus, self).__init__()
        self.__datasource       = datasource
        self.cur_idx            = index
        self.genotype_data      = None

    def next(self):
        """Move to the next valid locus.

        Will only return valid loci or exit via StopIteration exception

        """
        while True:
            self.cur_idx += 1
            if self.__datasource.populate_iteration(self):
                return self

        raise StopIteration


    def __iter__(self):
        """Basic iterator functionality. """

        return self
