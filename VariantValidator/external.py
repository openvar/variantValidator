# -*- coding: utf-8 -*-

def get_urls(dict_out):
    # Add urls
    report_urls = {}
    if 'NM_' in dict_out['hgvs_transcript_variant'] or 'NR_' in dict_out['hgvs_transcript_variant']:
        report_urls['transcript'] = 'https://www.ncbi.nlm.nih.gov' \
                                    '/nuccore/%s' % dict_out['hgvs_transcript_variant'].split(':')[0]
    if 'NP_' in str(predicted_protein_variant):
        report_urls['protein'] = 'https://www.ncbi.nlm.nih.gov' \
                                 '/nuccore/%s' % str(format_p).split(':')[0]
    if 'NG_' in dict_out['hgvs_refseqgene_variant']:
        report_urls['refseqgene'] = 'https://www.ncbi.nlm.nih.gov' \
                                    '/nuccore/%s' % dict_out['hgvs_refseqgene_variant'].split(':')[0]
    if 'LRG' in dict_out['hgvs_lrg_variant']:
        lrg_id = dict_out['hgvs_lrg_variant'].split(':')[0]
        lrg_data = va_dbCrl.data.get_LRG_data_from_LRGid(lrg_id)
        lrg_status = str(lrg_data[4])
        if lrg_status == 'public':
            report_urls['lrg'] = 'http://ftp.ebi.ac.uk/pub' \
                                 '/databases/lrgex/%s.xml' % dict_out['hgvs_lrg_variant'].split(':')[0]
        else:
            report_urls['lrg'] = 'http://ftp.ebi.ac.uk' \
                                 '/pub/databases/lrgex' \
                                 '/pending/%s.xml' % dict_out['hgvs_lrg_variant'].split(':')[0]
    # Ensembl needs to be added at a later data
    # "http://www.ensembl.org/id/" ? What about historic versions?????

    # Varsome and UCSC
    report_urls['external'] = {}
    for primary_assembly, data_set in dict_out['primary_assembly_loci'][]

        if not 'hg' in primary_assembly:
            continue

        if primary_assembly == 'hg19':
            varsome = "https://varsome.com/variant/hg19/"  # %s" %(coding.replace('dup', 'ins'))
        if primary_assembly == 'hg38':
            varsome = "https://varsome.com/variant/hg38/"  # %s" %(coding.replace('dup', 'ins'))

        # Report VCF from hgvs
        rp_vcf_component_list = [str(dataset['vcf']['chr']),
                                 str(dataset['vcf']['pos']),
                                 str(dataset['vcf']['ref']),
                                 str(dataset['vcf']['alt'])
                                 ]
        vcf_varsome = '-'.join(rp_vcf_component_list)
        varsome = varsome + vcf_varsome

        # Create a link to UCSC genome browser
        ucsc_assembly = primary_assembly
        vcf_components = vcf_varsome.split('-')
        vcf_components[0] = ucsc_chromosome
        vcf_varsome = '-'.join(vcf_components)

        browser_start = str(final_hgvs_genomic.posedit.pos.start.base - 11)
        browser_end = str(final_hgvs_genomic.posedit.pos.end.base + 11)
        ucsc_browser_position = '%s:%s-%s' % (ucsc_chromosome, browser_start, browser_end)
        coding = 'intergenic'
        warnings.warn(coding)
        ucsc_link = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s&hgt.customText=https://variantvalidator.org/bed/?variant=%s|%s|GRCh37|%s|%s' % (
        ucsc_assembly, ucsc_browser_position, coding, final_hgvs_genomic.ac, valstr(final_hgvs_genomic), vcf_varsome)
        # add links
        vcf_varsome = primary_assembly + ':' + vcf_varsome.replace('-', ':')
        report_urls['external'][vcf_varsome]['varsome'] = varsome
        report_urls['external'][vcf_varsome]['ucsc'] = ucsc_link


# <LICENSE>
# Copyright (C) 2018  Peter Causey-Freeman, University of Leicester
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