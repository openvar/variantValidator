import re
from . import vvFunctions as fn


class Variant(object):
    """
    This Variant object will contain the original input, the processed variant description and any other data that's
    relevant to what kind of variant it is.
    """

    def __init__(self, original, quibble=None, warnings='', write=True, primary_assembly=False, order=False):
        self.original = original
        if quibble is None:
            self.quibble = original
        else:
            self.quibble = quibble
        self.hgvs_formatted = None
        self.hgvs_genomic = None
        self.stashed = None
        self.trapped = None
        self.input_parses = None

        self.warnings = warnings
        self.description = ''  # hgnc_gene_info variable
        self.coding = ''
        self.coding_g = ''
        self.genomic_r = ''
        self.genomic_g = ''
        self.protein = ''
        self.write = write
        self.primary_assembly = primary_assembly
        self.order = order
        self.output_type_flag = 'warning'

        self.test_stash_tx_left = None
        self.test_stash_tx_right = None

        self.timing = {}

        self.refsource = None
        self.reftype = None

        self.hn = None
        self.reverse_normalizer = None
        self.evm = None
        self.no_norm_evm = None
        self.min_evm = None
        self.lose_vm = None

    def is_ascii(self):
        """
        Instead of the previous test for unicode rich text characters.
        Now going to test that all characters are within the ascii alphabet
        """
        try:
            self.quibble.encode('ascii')
            return True
        except UnicodeEncodeError or UnicodeDecodeError:
            # Will catch errors raised by python 2 and python 3
            return False

    def get_non_ascii(self):
        """
        Will return non ascii character positions within variant description
        :return:
        """
        chars = []
        positions = []

        for i, c in enumerate(self.quibble):
            try:
                c.encode('ascii')
            except UnicodeEncodeError or UnicodeDecodeError:
                chars.append(c)
                positions.append(i+1)

        return chars, positions

    def remove_whitespace(self):
        """
        Will remove all whitespace from quibble
        :return:
        """
        self.quibble = ''.join(self.quibble.split())

    def format_quibble(self):
        """
        Removes whitespace from the ends of the string
        Removes anything in brackets
        Identifies variant type (p. c. etc)
        Accepts c, g, n, r currently. And now P also 15.07.15
        """
        # Set regular expressions for if statements
        pat_gene = re.compile(r'\(.+?\)')  # Pattern looks for (....)

        if pat_gene.search(self.quibble):
            self.quibble = pat_gene.sub('', self.quibble)

        try:
            self.set_refsource()
        except fn.VariantValidatorError:
            return True

        try:
            self.set_reftype()
        except fn.VariantValidatorError:
            return True

        return False

    def set_reftype(self):
        """
        Method will set the reftype based on the quibble
        :return:
        """
        pat_est = re.compile(r'\d\:\d')

        if ':g.' in self.quibble:
            self.reftype = ':g.'
        elif ':r.' in self.quibble:
            self.reftype = ':r.'
        elif ':n.' in self.quibble:
            self.reftype = ':n.'
        elif ':c.' in self.quibble:
            self.reftype = ':c.'
        elif ':p.' in self.quibble:
            self.reftype = ':p.'
        elif ':m.' in self.quibble:
            self.reftype = ':m.'
        elif pat_est.search(self.quibble):
            self.reftype = 'est'
        else:
            raise fn.VariantValidatorError("Unable to identity reference type from %s" % self.quibble)

    def set_refsource(self):
        """
        Method will set the refsource based on the quibble
        :return:
        """
        if self.quibble.startswith('LRG'):
            self.refsource = 'LRG'
        elif self.quibble.startswith('ENS'):
            self.refsource = 'ENS'
        elif self.quibble.startswith('N'):
            self.refsource = 'RefSeq'
        else:
            raise fn.VariantValidatorError("Unable to identify reference source from %s" % self.quibble)

    def set_quibble(self, newval):
        """
        Method will set the quibble and reset the refsource and reftype
        :param newval:
        :return:
        """
        self.quibble = newval
        self.set_refsource()
        self.set_reftype()