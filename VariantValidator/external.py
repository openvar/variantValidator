


if primary_assembly == 'GRCh37':
    varsome = "https://varsome.com/variant/hg19/"  # %s" %(coding.replace('dup', 'ins'))
if primary_assembly == 'GRCh38':
    varsome = "https://varsome.com/variant/hg38/"  # %s" %(coding.replace('dup', 'ins'))

# Report VCF from hgvs
reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(final_hgvs_genomic)
rp_vcf_dict = variantanalyser.hgvs2vcf.report_hgvs2vcf(reverse_normalized_hgvs_genomic, primary_assembly)
rp_vcf_component_list = [str(rp_vcf_dict['chr']),
                         str(rp_vcf_dict['pos']),
                         str(rp_vcf_dict['ref']),
                         str(rp_vcf_dict['alt']),
                         ]
vcf_varsome = '-'.join(rp_vcf_component_list)
varsome = varsome + vcf_varsome
if primary_assembly == 'GRCh37':
    varsome_api = ['hg19', vcf_varsome]
if primary_assembly == 'GRCh38':
    varsome_api = ['hg38', vcf_varsome]

# Create a link to UCSC genome browser
if primary_assembly == 'GRCh37':
    ucsc_assembly = 'hg19'
if primary_assembly == 'GRCh38':
    ucsc_assembly = 'hg38'
ucsc_chromosome = variantanalyser.supported_chromosome_builds.to_chr_num_ucsc(final_hgvs_genomic.ac, ucsc_assembly)
if ucsc_chromosome is not None:
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
