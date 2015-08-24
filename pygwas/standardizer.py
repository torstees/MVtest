
import numpy
import pheno_covar

from exceptions import InvariantVar

def get_standardizer():
    global _standardizer
    return _standardizer
def set_standardizer(std):
    global _standardizer
    _standardizer = std

class StandardizedVariable(object):
    """Optional plugin object that can be used to standardize covariate and
       phenotype data.

       Many algorithms require that input be standardized in some way in order
       to work properly, however, rescaling the results is algorithm specific.
       In order to facilitate this situation, application authors can
       write up application specific Standardization objects for use with
       the data parsers.

       self.covariates_raw are unstandardized covariates with missingness over entire model
       self.covariates are standardized covariates
       self.phenotype_raw is unstandardized phenotype with missingness over entire model
       self.phenotype is standardized phenotype
       self.nonmissing is the mask associated with non-missing data (for entire group)
       """
    def __init__(self, pc):
        self.missing = []
        self.covar_count = len(pc.covariate_data)
        self.pheno_count = len(pc.phenotype_data)
        self.covariates = None
        self.phenotypes = None

        for pheno in pc.phenotype_data:
            missing = pheno == pheno_covar.PhenoCovar.missing_encoding
            for idx in range(0, self.covar_count):
                missing = missing | (pc.covariate_data[idx] == pheno_covar.PhenoCovar.missing_encoding)
            self.missing.append(missing)

        self.idx = 0
        self.datasource = pc


    def get_variables(self, missing_in_geno=None):
        count = 0
        mismatch = 0

        if missing_in_geno is None:
            nonmissing = numpy.invert(self.missing[self.idx])
        else:
            nonmissing = numpy.invert(self.missing[self.idx] | missing_in_geno)
        nmcount = sum(nonmissing)
        covars = numpy.zeros((self.covar_count, nmcount))
        for idx in range(0, self.covar_count):
            covars[idx] = self.covariates[idx][nonmissing]
            min = covars[idx][covars[idx] != pheno_covar.PhenoCovar.missing_encoding].min()
            max = covars[idx][covars[idx] != pheno_covar.PhenoCovar.missing_encoding].max()
            if min == max:
                raise InvariantVar("Covar %s doesn't have enough variation to continue" % (self.datasource.covariate_labels[idx]))
        min = self.phenotypes[self.idx][nonmissing].min()
        max = self.phenotypes[self.idx][nonmissing].max()
        if min == max:
            raise InvariantVar("Phenotype %s doesn't have enough variation to continue" % (self.datasource.phenotype_names[self.idx]))
        return (self.phenotypes[self.idx][nonmissing], covars, nonmissing)

    def get_phenotype_name(self):
        return self.datasource.phenotype_names[self.idx]

    def get_covariate_name(self, idx):
        return self.datasource.covariate_labels[idx]

    def get_covariate_names(self):
        return self.datasource.covariate_labels

    def standardize(self):
        pass

    def destandardize(self):
        pass


class NoStandardization(StandardizedVariable):
    def __init__(self, pc):
        super(NoStandardization, self).__init__(pc)
    def standardize(self):
        """Standardize the variables within a range [-1.0 and 1.0]

        This replaces the local copies of this data. When it's time to
        scale back, use destandardize from the datasource for that.

        """
        self.covariates = self.datasource.covariate_data
        self.phenotypes = self.datasource.phenotype_data

    def destandardize(self, estimates, se, **kwargs):
        """When the pheno/covar data has been standardized, this can be
        used to rescale the betas back to a meaningful value using the
        original data.

        For the "Un-standardized" data, we do no conversion.

        """

        return estimates, se, kwargs["pvalues"]

_standardizer = NoStandardization
