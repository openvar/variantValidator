import re
import string


class Variant(object):
    """
    This Variant object will contain the original input, the processed variant description and any other data that's
    relevant to what kind of variant it is.
    """

    def __init__(self, original):
        self.original = original
        self.quibble = original
        self.hgvs_formatted = original

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

