#! /usr/bin/env python

import argparse
import sys
from VariantValidator import Validator


def output_results(valoutput, outformat, with_meta):
    if outformat == 'dict':
        return str(valoutput.format_as_dict(with_meta=with_meta))
    elif outformat == 'json':
        return str(valoutput.format_as_json(with_meta=with_meta))
    else:
        # table format
        table = valoutput.format_as_table(with_meta=with_meta)
        newtable = []
        for row in table:
            if isinstance(row, list):
                newrow = '\t'.join(row)
            else:
                newrow = str(row)
            newtable.append(newrow)
        return '\n'.join(newtable)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--variant', required=True, nargs='+', help="Variant(s) to validate")
    parser.add_argument('-g', '--genome', nargs='?', default='GRCh37', choices=['GRCh37', 'GRCh38', 'hg19', 'hg38'],
                        help="Genome assembly (default: %(default)s)")
    parser.add_argument('-t', '--transcripts',  nargs='?', default='all',
                        help='Transcripts to output results for (default: %(default)s)')
    parser.add_argument('-s', '--submission', choices=['individual', 'batch'], default='individual',
                        help='Submit variants individually or as a single batch validation (default: %(default)s)')
    parser.add_argument('-f', '--output_format', choices=['dict', 'table', 'json'], default='dict',
                        help='Output validations as a list or as a dictionary (default: %(default)s)')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default='-',
                        help='Specifies the output file (default: stdout)')
    parser.add_argument('-m', '--meta', action='store_true', default=False,
                        help='Also output metadata (default: %(default)s)')

    args = parser.parse_args()

    validator = Validator()

    if args.submission == 'individual':
        for variant in args.variant:
            output = validator.validate(variant, args.genome, args.transcripts)
            args.output.write(output_results(output, args.output_format, args.meta) + '\n')
    else:
        batch = '|'.join(args.variant)
        sys.stderr.write("Submitting batch query: %s\n" % batch)
        output = validator.validate(batch, args.genome, args.transcripts)
        args.output.write(output_results(output, args.output_format, args.meta) + '\n')

# <LICENSE>
# Copyright (C) 2016-2022 VariantValidator Contributors
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>
