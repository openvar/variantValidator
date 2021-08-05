import os
import sys
import VariantValidator
vval = VariantValidator.Validator()
cwd = os.path.dirname(os.path.abspath(__file__))


# Check Args
if len(sys.argv) != 3:
    print('Too few arguments. The command required is: python bin/batch_validator.py genome_build select_transcripts')
    exit()

# Check genome_build
genome_build = sys.argv[1]
genomes = ['GRCh38', 'hg38', 'GRCh37', 'hg19']
if genome_build not in genomes:
    warn = '%s is not a supported genome build' % genome_build
    print(warn)
    exit()

select_transcripts = sys.argv[2]

# Loop through file and detect fails
filepath = cwd.replace('bin', 'batch')
infile = filepath + '/input.txt'
outfile = open(filepath + "/output.txt", "w")
counter = 0
with open(infile) as fp:
    variants = fp.readlines()
    for variant in variants:
        variant = variant.strip()
        print(variant)
        try:
            validate = vval.validate(variant, genome_build, select_transcripts)
            validation = validate.format_as_table(with_meta=True)
            if counter == 0:
                line_counter = 0
                for line in validation:
                    if line_counter == 0:
                        outfile.write(line + '\n')
                    elif line_counter == 1:
                        ln_cat = '\t'.join(line)
                        outfile.write(ln_cat + '\n')
                    else:
                        copy_line = []
                        for element in line:
                            if element == '':
                                copy_line.append('None')
                            elif element is None:
                                copy_line.append('None')
                            else:
                                copy_line.append(element)
                        ln_cat = '\t'.join(copy_line)
                        outfile.write(ln_cat + '\n')
                    line_counter = line_counter + 1
            else:
                line_counter = 0
                for line in validation:
                    if line_counter == 0:
                        pass
                    elif line_counter == 1:
                        ln_cat = '\t'.join(line)
                        pass
                    else:
                        copy_line = []
                        for element in line:
                            if element == '':
                                copy_line.append('None')
                            elif element is None:
                                copy_line.append('None')
                            else:
                                copy_line.append(element)
                        ln_cat = '\t'.join(copy_line)
                        outfile.write(ln_cat + '\n')
                    line_counter = line_counter + 1
        except VariantValidator.modules.utils.VariantValidatorError as e:
            outfile.close()
            print(variant)
            print(e)
            exit()
        counter = counter + 1

print('Processing complete')
outfile.close()

# <LICENSE>
# Copyright (C) 2016-2021 VariantValidator Contributors
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
