#! /usr/bin/env python

import argparse
import json
import sys
from VariantValidator import Validator


def output_results(valoutput, outformat):
    if outformat == 'dict':
        return str(valoutput.format_as_dict())
    elif outformat == 'json':
        return json.dumps(valoutput.format_as_dict())
    else:
        return str(valoutput.format_as_dict())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--variant', required=True, nargs='+', help="Variant(s) to validate")
    parser.add_argument('-g', '--genome', nargs='?', default='GRCh37', choices=['GRCh37', 'GRCh38', 'hg19', 'hg38'],
                        help="Genome assembly (default: %(default)s)")
    parser.add_argument('-t', '--transcripts',  nargs='?', default='all',
                        help='Transcripts to output results for (default: %(default)s)')
    parser.add_argument('-s', '--submission', choices=['individual', 'batch'], default='individual',
                        help='Submit variants individually or as a single batch validation (default: %(default)s')
    parser.add_argument('-f', '--output_format', choices=['dict', 'list', 'json'], default='dict',
                        help='Output validations as a list or as a dictionary (default: %(default)s')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default='-',
                        help='Specifies the output file (default: stdout)')

    args = parser.parse_args()

    validator = Validator()

    if args.submission == 'individual':
        for variant in args.variant:
            output = validator.validate(variant, args.genome, args.transcripts)
            args.output.write(output_results(output, args.output_format) + '\n')
    else:
        batch = '|'.join(args.variant)
        sys.stderr.write("Submitting batch query: %s\n" % batch)
        output = validator.validate(batch, args.genome, args.transcripts)
        args.output.write(output_results(output, args.output_format) + '\n')
