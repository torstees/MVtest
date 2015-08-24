class ReportableException(Exception):
    """Simple exeception with message"""
    def __init__(self, msg):
        self.msg = msg


class UnsolvedLocus(ReportableException):
    def __init__(self, msg):
        super(UnsolvedLocus, self).__init__(msg)

class TooManyAlleles(ReportableException):
    """Indicate locus found with more than 2 alleles"""
    def __init__(self, chr=None, rsid=None, pos=None, alleles=None, index=None, prefix="Too many alleles: "):
        self.chr = chr
        self.pos = pos
        self.rsid = rsid
        self.alleles = alleles
        self.index = index
        super(TooManyAlleles, self).__init__("%s %s:%s (%s) %s" % (prefix, self.chr, self.pos, self.rsid, self.alleles))

class NanInResult(ReportableException):
    def __init__(self, msg = ""):
        super(NanInResult, self).__init__(msg)

class NoMatchedPhenoCovars(ReportableException):
    def __init__(self, msg = ""):
        super(NoMatchedPhenoCovars, self).__init__(msg)

class InvariantVar(ReportableException):
    def __init__(self, msg=""):
        super(InvariantVar, self).__init__(msg)

class TooFewAlleles(TooManyAlleles):
    """Indicate fixed allele was found"""
    def __init__(self, chr=None, rsid=None, pos=None, alleles=None, index=None):
        super(TooFewAlleles, self).__init__(chr, rsid, pos, alleles, index, "Too few alleles: ")


class InvalidBoundarySpec(ReportableException):
    """Indicate boundary specification was malformed or non-sensical"""
    def __init__(self, malformed_boundary):
        self.malformed_boundary = malformed_boundary
        super(InvalidBoundarySpec, self).__init__("%s" % (malformed_boundary))


class MalformedInputFile(ReportableException):
    """Error encountered in data from an input file"""
    def __init__(self, msg):
        super(MalformedInputFile, self).__init__(msg)


class InvalidSelection(MalformedInputFile):
    """Indicate that the user provided input that is meaningless.

    This is likely a situation where the user provided an invalid name
    for a phenotype or covariate. Probably a misspelling.

    """
    pass