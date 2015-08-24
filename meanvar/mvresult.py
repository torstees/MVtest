import sys
import math
__author__ = 'torstees'
__version__ = "0.9.1"

class MVResult(object):
    """Result associated with a single locus/phenotype execution

    Class Members:
    chr         Chromosome Number
    pos         BP offset from beginning of chromosome (starts with 1)
    rsid        SNP name
    maj_allele  Most common allele representation found at this locus
    min_allele  Least common allele representation found at this locus
    non_miss    Number of individuals who weren't missing any of the
                components (genotype, outcome or variate values)
    fit         Fitness value returned from search
    fit0        Fitness returned from null distribution
    p_mvtest    p-value associated with the evaluation
    p_variance
    beta_lin    Betas associated with the mean
    beta_var    Betas associated with the variance
    ph_label    Phenotype label associated with this result
    maf         frequency of min_allele
    runtime     number seconds required to perform calculation

    """
    verbose = False

    #def __init__(self, chr, pos, rsid, maj, min, non_miss_count, ph_label, fit, fit0, p_mvtest, p_variance, blin, bvar, maf, runtime=-1):
    def __init__(self, chr, pos, rsid, maj, min, eff_alcount, non_miss_count, p_mvtest, ph_label, beta_values , pvalues, stderrors, maf, covar_labels=[], lm=-1,runtime=-1):
        self.chr = chr
        self.pos = pos
        self.rsid = rsid
        self.maj_allele = maj
        self.min_allele = min
        self.eff_alcount = eff_alcount
        self.non_miss = non_miss_count
        self.p_mvtest = p_mvtest
        self.betas = beta_values
        self.beta_pvalues = pvalues
        self.beta_stderr = stderrors
        self.ph_label   = ph_label
        self.maf = maf
        self.runtime = runtime
        self.lmpv = lm
        self.covar_labels = covar_labels
        #self.stdy = stdy
        #self.meanAVA = meanAVA
        #self.varAVA = varAVA

    @classmethod
    def set_verbose(cls, val):
        cls.verbose = val

    @property
    def p_variance(self):
        return self.beta_pvalues[3]

    def print_header(self, f=sys.stdout, verbose=False):
        self.var_count = 2 + len(self.covar_labels)

        if verbose:
            header = [
                "Chr",
                "Pos",
                "RSID",
                "Phenotype",
                "N",
                "Ref_allele",
                "Eff_allele",
                "Eff_Allele_Freq",
                "P-Value",
                "LM_PValue"
            ]
            for var in ["Intercept","Geno"] + self.covar_labels:
                for t in ["mean", "mean_stder", "mean_pval", "var", "var_stder", "var_pval"]:
                    header.append("%s_%s" % (var.lower(), t))
            print >> f, "\t".join(header)
        else:
            print >> f, "\t".join([
                "Chr",
                "Pos",
                "RSID",
                "Phenotype",
                "N",
                "Ref_allele",
                "Eff_allele",
                "Eff_Allele_Freq",
                "P-Value",
                "geno_mean",
                "geno_mean_stderr",
                "geno_mean_pval",
                "geno_var",
                "geno_var_stderr",
                "geno_var_pval",
            ])


    def stringify(self, value):
        try:
            return "%.3e" % (value)
        except:
            return value

    def print_result(self, f=sys.stdout, verbose=False):
        var_count = len(self.betas)/2
        if verbose:
            results = [str(x) for x in [
                self.chr,
                self.pos,
                self.rsid,
                self.ph_label,
                self.non_miss,
                self.maj_allele,
                self.min_allele,
                "%.4e" % self.eff_alcount,
                "%.4e" % self.p_mvtest,
                "%.4e" % self.lmpv
            ]]

            for i in range(var_count):
                results.append("%.3e" % self.betas[i])
                results.append(self.stringify(self.beta_stderr[i]))
                results.append(self.stringify(self.beta_pvalues[i]))
                results.append(self.stringify(self.betas[i+var_count]))
                results.append(self.stringify(self.beta_stderr[i+var_count]))
                results.append(self.stringify(self.beta_pvalues[i+var_count]))
            print >> f, "\t".join([str(x) for x in results])

        else:
            print >> f, "\t".join(str(x) for x in [
                self.chr,
                self.pos,
                self.rsid,
                self.ph_label,
                self.non_miss,
                self.maj_allele,
                self.min_allele,
                "%.3e" % self.eff_alcount,
                "%.3e" % self.p_mvtest,
                "%.3e" % self.betas[1],
                self.stringify(self.beta_stderr[1]),
                "%.3e" % self.beta_pvalues[1],
                "%.3e" % self.betas[1+var_count],
                self.stringify(self.beta_stderr[1+var_count]),
                "%.3e" % self.beta_pvalues[1+var_count],
            ])


