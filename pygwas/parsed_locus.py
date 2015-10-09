from locus import Locus

class ParsedLocus(Locus):
    """Locus data representing current iteration from a dataset

    Provide an iterator interface for all dataset types.

    """

    def __init__(self, datasource, index=-1):
        """Basic initialization (nothing is currently valid)"""

        super(ParsedLocus, self).__init__()
        #: Reference back to the parser that generated this object
        self.__datasource       = datasource
        #: Index within the list of loci being analyzed
        self.cur_idx            = index
        #: Actual genotype data for this locus
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
