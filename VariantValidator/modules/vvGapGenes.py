"""
GapGenes contains code which handles instances where there is a discrepancy in the transcipt exon and the equivalent
alighed genome exon. This is only relevant to RefSeq and only in instances of the following blacklisted genes
"""

# -*- coding: utf-8 -*-

# Import python modules
from distutils.version import StrictVersion
import re
import copy
import hgvs.exceptions
import hgvs.assemblymapper
import hgvs.variantmapper

# These will need to be fixed
# import hgvs2vcf as va_H2V
# import variantanalyser.functions as va_func

# from gap_genes
def gap_black_list(symbol):
    """
    Lists of genes for GRCh37 and GRCh38 which require a gap to be inserted into either the
    transcript or the genome to maintain a perfect alignment
    """
    gapGene = {
                    "LPP": "",
                    "VPS13D": "",
                    "SSPO": "",
                    "HTT": "",
                    "PRKDC": "",
                    "RNA45SN4": "",
                    "RNA45SN1": "",
                    "RNA45SN2": "",
                    "RNA45SN3": "",
                    "ALMS1": "",
                    "ZNF141": "",
                    "PRLR": "",
                    "NBPF10": "",
                    "ACACA": "",
                    "ZMYM2": "",
                    "MIAT": "",
                    "WDFY4": "",
                    "CECR2": "",
                    "FAM30A": "",
                    "MYO15B": "",
                    "CELF2": "",
                    "JRK": "",
                    "PTEN": "",
                    "ZNF714": "",
                    "MGAT4C": "",
                    "SLITRK4": "",
                    "ZAN": "",
                    "COL19A1": "",
                    "CCDC144B": "",
                    "RAB11FIP4": "",
                    "ZNF516": "",
                    "ZNF518A": "",
                    "PROX1": "",
                    "HCG18": "",
                    "SON": "",
                    "ARMC9": "",
                    "CAMK1D": "",
                    "GRIP2": "",
                    "KLHL5": "",
                    "PPIP5K2": "",
                    "PKD1L2": "",
                    "SLC7A2": "",
                    "DGKK": "",
                    "IQSEC1": "",
                    "SYNM": "",
                    "SARM1": "",
                    "SMAD5": "",
                    "MAML3": "",
                    "CXorf40A": "",
                    "MAPT": "",
                    "ITIH5": "",
                    "NOTCH4": "",
                    "FER1L4": "",
                    "CNTNAP4": "",
                    "NLRC3": "",
                    "COL18A1": "",
                    "SLC6A6": "",
                    "DDX52": "",
                    "CDH4": "",
                    "SLC46A1": "",
                    "SLC35E2B": "",
                    "OCLN": "",
                    "DCAF7": "",
                    "SCAMP1": "",
                    "ATG13": "",
                    "SMAD3": "",
                    "DDX6": "",
                    "SLC25A53": "",
                    "ALG9": "",
                    "DCP1A": "",
                    "NCAM1": "",
                    "LINC00869": "",
                    "MYH7": "",
                    "DIXDC1": "",
                    "ZBTB4": "",
                    "RABEP1": "",
                    "PVR": "",
                    "POM121C": "",
                    "HOOK1": "",
                    "MAPK8IP2": "",
                    "ZNF280B": "",
                    "WASF2": "",
                    "PLEKHA2": "",
                    "PPP4R3B": "",
                    "FAM83H": "",
                    "SALL3": "",
                    "PHKG2": "",
                    "C18orf25": "",
                    "ZNF229": "",
                    "ZNF765-ZNF761": "",
                    "KANSL1": "",
                    "FAM102B": "",
                    "NOTCH2NL": "",
                    "YTHDF3": "",
                    "DPCR1": "",
                    "DACH1": "",
                    "PKD1L3": "",
                    "GRIA3": "",
                    "CYP1B1": "",
                    "LTBP4": "",
                    "SPON1": "",
                    "RNA28SN4": "",
                    "RNA28SN1": "",
                    "TRIL": "",
                    "RNA28SN3": "",
                    "RNA28SN2": "",
                    "XKR5": "",
                    "RBM8A": "",
                    "SALL2": "",
                    "JADE3": "",
                    "DHX57": "",
                    "PIGN": "",
                    "CPNE3": "",
                    "ANO1": "",
                    "NATD1": "",
                    "DKFZP434A062": "",
                    "TDRD9": "",
                    "BDNF": "",
                    "IVD": "",
                    "STIMATE": "",
                    "KCP": "",
                    "PRAG1": "",
                    "KLHL18": "",
                    "LYNX1": "",
                    "HYOU1": "",
                    "HLA-L": "",
                    "ATG9B": "",
                    "SLC6A14": "",
                    "PCSK6": "",
                    "MIR99AHG": "",
                    "TOX4": "",
                    "GABBR1": "",
                    "RABGEF1": "",
                    "PRR36": "",
                    "MAP3K14": "",
                    "PCDHB9": "",
                    "LOC102723753": "",
                    "MYO19": "",
                    "SRSF8": "",
                    "CTPS2": "",
                    "AHCYL1": "",
                    "UHRF1": "",
                    "MARCKS": "",
                    "ZMYM1": "",
                    "SENP3-EIF4A1": "",
                    "SEC14L2": "",
                    "RAPGEFL1": "",
                    "ZNF761": "",
                    "CNTROB": "",
                    "SSTR3": "",
                    "PAX2": "",
                    "GGA3": "",
                    "MCL1": "",
                    "EPS8": "",
                    "LINC02210": "",
                    "KRBA1": "",
                    "MSH5-SAPCD1": "",
                    "HLA-DPB1": "",
                    "PPP1R9B": "",
                    "OPLAH": "",
                    "UBXN4": "",
                    "ZNF2": "",
                    "EPHB6": "",
                    "LIX1L": "",
                    "RAPGEF4": "",
                    "MED22": "",
                    "POLR3C": "",
                    "DDR1": "",
                    "SIGLEC16": "",
                    "NEFL": "",
                    "ABCG4": "",
                    "BAG6": "",
                    "RECQL4": "",
                    "SPPL2B": "",
                    "RETREG3": "",
                    "FZD6": "",
                    "SCRT1": "",
                    "LSM14A": "",
                    "TAPBP": "",
                    "TWSG1": "",
                    "FRMD8": "",
                    "VPS26C": "",
                    "PNMA3": "",
                    "ZNF282": "",
                    "SP8": "",
                    "SRRM3": "",
                    "CCDC125": "",
                    "NPIPB3": "",
                    "FAM13C": "",
                    "GTF2IP1": "",
                    "ANKRD34A": "",
                    "PPP1R2": "",
                    "PHYHIPL": "",
                    "USH1G": "",
                    "LINC00461": "",
                    "ZNRD1ASP": "",
                    "TRIM10": "",
                    "SPIB": "",
                    "BCL6B": "",
                    "SCARF2": "",
                    "KIR3DX1": "",
                    "LOC400682": "",
                    "HLA-DOA": "",
                    "PLCD3": "",
                    "VPS11": "",
                    "FAM231D": "",
                    "TRIM52": "",
                    "ABCF1": "",
                    "ANP32E": "",
                    "COPG2IT1": "",
                    "TGIF2": "",
                    "LHX1": "",
                    "PIK3R6": "",
                    "APOL4": "",
                    "ZNF502": "",
                    "FGD5P1": "",
                    "LINC00624": "",
                    "ADRA2B": "",
                    "ZNF598": "",
                    "GNAZ": "",
                    "TMEM106A": "",
                    "SLC12A9": "",
                    "TCF19": "",
                    "CCDC3": "",
                    "EFHC2": "",
                    "KCNE1B": "",
                    "PBX2": "",
                    "PAMR1": "",
                    "GJA5": "",
                    "TYW1B": "",
                    "PLP1": "",
                    "ANKDD1A": "",
                    "GBE1": "",
                    "MAMDC2": "",
                    "PIGW": "",
                    "MOCOS": "",
                    "GRIPAP1": "",
                    "COL26A1": "",
                    "MAPT-IT1": "",
                    "SRRT": "",
                    "ZNF595": "",
                    "SEMA3B": "",
                    "C21orf58": "",
                    "RHBDF1": "",
                    "EGR2": "",
                    "ABRAXAS2": "",
                    "NPRL3": "",
                    "TXNIP": "",
                    "RYK": "",
                    "RXRB": "",
                    "LILRB2": "",
                    "SYT3": "",
                    "TRPV6": "",
                    "PARG": "",
                    "CSNK1G2": "",
                    "ARHGEF16": "",
                    "HSH2D": "",
                    "ALDH3B1": "",
                    "ZNF274": "",
                    "MUC13": "",
                    "LINC00842": "",
                    "AKT1": "",
                    "CHM": "",
                    "ZSCAN26": "",
                    "MAL2": "",
                    "PTH2R": "",
                    "GPANK1": "",
                    "LINC01623": "",
                    "CD86": "",
                    "RHBG": "",
                    "TMSB15B": "",
                    "ZCCHC3": "",
                    "TUBB": "",
                    "POLDIP2": "",
                    "PRMT3": "",
                    "PPT2-EGFL8": "",
                    "LINC02210-CRHR1": "",
                    "KIFC1": "",
                    "USP27X": "",
                    "HDGFL2": "",
                    "FOXI3": "",
                    "PAH": "",
                    "P3H3": "",
                    "CRHR1": "",
                    "LOC101927759": "",
                    "ARFRP1": "",
                    "C3orf38": "",
                    "DAXX": "",
                    "SLC37A4": "",
                    "IQCA1L": "",
                    "MMP28": "",
                    "LINC02197": "",
                    "NECAP1": "",
                    "CDSN": "",
                    "LOC440570": "",
                    "B3GNT6": "",
                    "AOAH": "",
                    "GAS2L1": "",
                    "MPIG6B": "",
                    "CDK11B": "",
                    "ASPN": "",
                    "HSPA1B": "",
                    "LOC100508631": "",
                    "MICB": "",
                    "LOC102724580": "",
                    "SENP3": "",
                    "RBM38": "",
                    "TMC4": "",
                    "LILRB5": "",
                    "C6orf47": "",
                    "RIOX1": "",
                    "BHLHE40-AS1": "",
                    "SRD5A2": "",
                    "TSEN34": "",
                    "EI24": "",
                    "PADI6": "",
                    "LINC00893": "",
                    "CYP2D7": "",
                    "LINC01622": "",
                    "LINC01879": "",
                    "REC8": "",
                    "UNC93B1": "",
                    "POU5F1": "",
                    "GPIHBP1": "",
                    "FOXD1": "",
                    "GPSM1": "",
                    "MICA": "",
                    "UGT2B15": "",
                    "KIZ": "",
                    "ARL17A": "",
                    "PRAMEF36P": "",
                    "HCG22": "",
                    "RNF39": "",
                    "BECN1": "",
                    "MOG": "",
                    "PROSER3": "",
                    "LINC01149": "",
                    "CYP21A2": "",
                    "PRAMEF18": "",
                    "TBC1D3G": "",
                    "NR2E3": "",
                    "NR1H2": "",
                    "VEGFC": "",
                    "TBC1D3F": "",
                    "C18orf65": "",
                    "HOXC11": "",
                    "TRY2P": "",
                    "LINC01138": "",
                    "LINC00243": "",
                    "HCG4": "",
                    "GBAP1": "",
                    "LYPD4": "",
                    "FAM226A": "",
                    "ZNF787": "",
                    "CYP11A1": "",
                    "EEF1A2": "",
                    "SLC38A5": "",
                    "MICB-DT": "",
                    "ZNF852": "",
                    "LOC441242": "",
                    "RNF115": "",
                    "SMA4": "",
                    "TAZ": "",
                    "LENG9": "",
                    "STRAP": "",
                    "CYP4F8": "",
                    "TSPAN10": "",
                    "KIR3DL1": "",
                    "HCP5B": "",
                    "MMP12": "",
                    "STAG3L2": "",
                    "GOLGA6L17P": "",
                    "ZBTB12": "",
                    "TREH": "",
                    "PMCHL2": "",
                    "LAGE3": "",
                    "ATRNL1": "",
                    "CEACAM20": "",
                    "ZG16": "",
                    "MIR3936HG": "",
                    "LOC102724562": "",
                    "INTS4P2": "",
                    "LINC00221": "",
                    "DHRS3": "",
                    "HCG27": "",
                    "CLTB": "",
                    "KLK6": "",
                    "HLA-H": "",
                    "SPANXA2-OT1": "",
                    "PRAMEF11": "",
                    "PPP1R11": "",
                    "NDUFA6-AS1": "",
                    "ECHDC3": "",
                    "HLA-DQB1": "",
                    "KIR2DS4": "",
                    "HLA-B": "",
                    "LOC102725121": "",
                    "CIB2": "",
                    "KIR2DL1": "",
                    "KIR2DL2": "",
                    "HLA-C": "",
                    "ABO": "",
                    "KRTAP10-7": "",
                    "HLA-G": "",
                    "CWC15": "",
                    "C17orf100": "",
                    "HLA-J": "",
                    "OR4K3": "",
                    "HLA-DQA1": "",
                    "LOC105379550": "",
                    "MRPS21": "",
                    "SIGLEC17P": "",
                    "LINC01115": "",
                    "NUDT18": "",
                    "ORAI1": "",
                    "PNLIPRP2": "",
                    "KLF14": "",
                    "SSX2B": "",
                    "CCL15-CCL14": "",
                    "UBXN8": "",
                    "IGFBP2": "",
                    "TMEM44-AS1": "",
                    "TEX13A": "",
                    "LCA10": "",
                    "SPANXN2": "",
                    "SYCE1": "",
                    "LILRA5": "",
                    "KRTAP5-4": "",
                    "FAM228B": "",
                    "OR12D1": "",
                    "SPC25": "",
                    "FCGR1CP": "",
                    "OR52E1": "",
                    "NOP16": "",
                    "EGFL8": "",
                    "PRAF2": "",
                    "LOC388282": "",
                    "CCNQ": "",
                    "VN1R3": "",
                    "HLA-V": "",
                    "SBK3": "",
                    "LOC100128594": "",
                    "KLRF1": "",
                    "EMG1": "",
                    "TARM1": "",
                    "UBE2NL": "",
                    "OR5AL1": "",
                    "TPSB2": "",
                    "PSORS1C2": "",
                    "HLA-DQA2": "",
                    "OR10AC1": "",
                    "OR2J1": "",
                    "OR10J4": "",
                    "CSNK2B": "",
                    "OR4Q2": "",
                    "LOC100507547": "",
                    "ZNF630-AS1": "",
                    "HLA-DMA": "",
                    "OR4E1": "",
                    "PRB3": "",
                    "CCL15": "",
                    "C8orf59": "",
                    "PSMB9": "",
                    "LINC01719": "",
                    "CT45A1": "",
                    "BST2": "",
                    "NCF4-AS1": "",
                    "FOLR3": "",
                    "KRTAP9-9": "",
                    "COPZ2": "",
                    "LYNX1-SLURP2": "",
                    "SAPCD1": "",
                    "PSORS1C1": "",
                    "ZNF793-AS1": "",
                    "ZNRD1": "",
                    "FRG1CP": "",
                    "LINC02362": "",
                    "KRTAP4-1": "",
                    "PICSAR": "",
                    "TWIST2": "",
                    "LINC01796": "",
                    "HCG25": "",
                    "KRTAP7-1": "",
                    "CRLF2": "",
                    "MDH2": "",
                    "HCG9": "",
                    "ATP5MC1": "",
                    "TTTY14": "",
                    "LOC100507384": "",
                    "PMS2P2": "",
                    "HCG23": "",
                    "LINC00226": "",
                    "RPP21": "",
                    "GPHB5": "",
                    "GAGE8": "",
                    "GAGE2E": "",
                    "LOC101928087": "",
                    "GAGE12B": "",
                    "GRIFIN": "",
                    "LOC102725193": "",
                    "HCG14": "",
                    "IFITM4P": "",
                    "SNORD48": "",
                    "MUC22": "",
                    "PTPRQ": "",
                    "HERC2": "",
                    "OTUD7A": "",
                    "LOC646214": "",
                    "TJP1": "",
                    "WDR81": "",
                    "KLF13": "",
                    "POLR2A": "",
                    "LOC100288637": "",
                    "GOLGA8N": "",
                    "GOLGA8J": "",
                    "GOLGA8K": "",
                    "GOLGA8R": "",
                    "MTMR10": "",
                    "SMIM10L1": "",
                    "KLLN": "",
                    "LINC02249": "",
                    "APBA2": "",
                    "CHRNA7": "",
                    "DBET": "",
                    "WNT3": "",
                    "GOLGA2P10": "",
                    "CHRFAM7A": "",
                    "RPH3AL": "",
                    "SORD2P": "",
                    "LINC00552": "",
                    "MPV17L": "",
                    "SLC22A18AS": "",
                    "C16orf45": "",
                    "GRK1": "",
                    "FRG2": "",
                    "LOC143666": "",
                    "FRG2EP": "",
                    "LOC105373100": "",
                    "GOLGA8Q": "",
                    "HERC2P7": "",
                    "SLC22A18": "",
                    "METRNL": "",
                    "BTNL2": "",
                    "ADAM18": "",
                    "PRSS22": "",
                    "C2orf27B": "",
                    "C2orf27A": "",
                    "LOC283710": "",
                    "LOC101928804": "",
                    "IFI27": "",
                    "ABCC6": "",
                    "LOC692247": ""
                }
    is_it_gapped = gapGene.get(symbol)
    if is_it_gapped == '':
        return True
    else:
        return False


"""
Head function for g_to_t mapping
Runs quick tests to see if gap compensation is required

RefSeq only
hgvs <= 1.1.3
"""


def compensate_g_to_t(hgvs_tx, hgvs_genomic, un_norm_hgvs_genomic, vm,
                      hn, reverse_normalizer, primary_assembly, hdp, hp, sf, hgvs_version):
    # Not required in these instances
    if re.match('ENST', hgvs_tx.ac) or (StrictVersion(str(hgvs_version)) >
                                        StrictVersion('1.1.3') is True):
        # Push to absolute position
        normalized_tx = fully_normalize(hgvs_tx, hgvs_genomic, hn, reverse_normalizer,
                                        hdp)
        hgvs_tx_returns = [normalized_tx, False, None, None, None]

    else:
        gene_symbol = hdp.get_tx_identity_info(hgvs_tx.ac)[6]
        # Check the blaccklist
        gap_compensation = gap_black_list(gene_symbol)
        if gap_compensation is False:
            normalized_tx = fully_normalize(hgvs_tx, hgvs_genomic, hn,
                                            reverse_normalizer, hdp)
            hgvs_tx_returns = [normalized_tx, False, None, None, None]
        else:
            # At this stage, we know that:
            # the gene is on the gap list
            # the hgvs version is <= 1.1.3
            # The requested transcript set is RefSeq
            gap_compensated_tx = g_to_t_compensation_code(hgvs_tx, hgvs_genomic,
                                                          un_norm_hgvs_genomic,
                                                          vm, hn,
                                                          reverse_normalizer,
                                                          primary_assembly,
                                                          hdp, hp, sf)
            if gap_compensated_tx[1] is False:
                hgvs_tx_returns = fully_normalize(hgvs_tx, hgvs_genomic, hn,
                                                  reverse_normalizer, hdp)
            else:
                hgvs_tx_returns = gap_compensated_tx

    hgvs_tx_dict = {'hgvs_transcript': hgvs_tx_returns[0],
                    'position_lock': hgvs_tx_returns[1],
                    'gapped_alignment_warning': hgvs_tx_returns[3],
                    'corrective_action': hgvs_tx_returns[2],
                    'gap_position': hgvs_tx_returns[-1],
                    'transcript_accession': hgvs_tx_returns[0].ac
                    }
    return hgvs_tx_dict


"""
Fully normalizes the hgvs_tx variant from the hgvs_genomic variant

Is only activated if the g_to_t_compensation_code IS NOT USED!
"""


def fully_normalize(hgvs_tx, hgvs_genomic, hn, reverse_normalizer, hdp):
    # set required variables
    tx_id = hgvs_tx.ac
    if re.match('ENST', tx_id):
        alt_aln_method = 'genebuild'
    else:
        alt_aln_method = 'splign'
    rhn = reverse_normalizer

    # Obtain the orientation of the transcript wrt selected genomic accession
    exon_alignments = hdp.get_tx_exons(tx_id, hgvs_genomic.ac, alt_aln_method)
    orientation = int(exon_alignments[0]['alt_strand'])
    # Normalize the genomic variant 5 prime if antisense or 3 prime if sense
    if orientation == -1:
        hgvs_genomic = rhn.normalize(hgvs_genomic)
    else:
        hgvs_genomic = hn.normalize(hgvs_genomic)
    try:
        hgvs_tx = hn.normalize(hgvs_tx)
    except hgvs.exceptions.HGVSError:
        pass

    return hgvs_tx


"""
Gap compensation code from genome to transcript

Source is VariantValidator

Requires an un-normalized genomic variant for stashing if available

Also requres a hgvs_genomic and tx id
"""


def g_to_t_compensation_code(hgvs_tx, hgvs_genomic, un_norm_hgvs_genomic, vm, hn,
                             reverse_normalizer, primary_assembly, hdp, hp, sf):
    """
    Gap aware projection from g. to c.
    """

    # Create mappers: has to be inline because requires primary_assembly and alt_aln_method
    alt_aln_method = 'splign'
    no_norm_evm = hgvs.assemblymapper.AssemblyMapper(hdp,
                                                     assembly_name=primary_assembly,
                                                     alt_aln_method=alt_aln_method,  # Only RefSeq should be here!!!
                                                     normalize=False,
                                                     replace_reference=True
                                                     )
    evm = hgvs.assemblymapper.AssemblyMapper(hdp,
                                             assembly_name=primary_assembly,
                                             alt_aln_method=alt_aln_method,  # Only RefSeq should be here!!!
                                             normalize=True,
                                             replace_reference=True
                                             )

    nr_vm = hgvs.variantmapper.VariantMapper(hdp, replace_reference=False)

    utilise_gap_code = True
    automap = ''

    # Set variables for problem specific warnings
    gapped_alignment_warning = None
    corrective_action_taken = None
    gapped_transcripts = ''
    auto_info = None

    # Set variables
    stash_input = un_norm_hgvs_genomic
    hgvs_genomic_variant = hgvs_genomic
    # Reverse normalize hgvs_genomic_variant: NOTE will replace ref
    reverse_normalized_hgvs_genomic = reverse_normalizer.normalize(hgvs_genomic_variant)
    hgvs_genomic_5pr = copy.deepcopy(reverse_normalized_hgvs_genomic)

    # Create a pseudo VCF so that normalization can be applied and a delins can be generated
    vcf_dict = va_H2V.hgvs2vcf(reverse_normalized_hgvs_genomic, primary_assembly, reverse_normalizer, sf)
    chr = vcf_dict['chr']
    pos = vcf_dict['pos']
    ref = vcf_dict['ref']
    alt = vcf_dict['alt']

    # Generate an end position
    end = str(int(pos) + len(ref) - 1)
    pos = str(pos)

    # take a look at the input genomic variant for potential base salvage
    stash_ac = vcf_dict['chr']
    stash_pos = int(vcf_dict['pos'])
    stash_ref = vcf_dict['ref']
    stash_alt = vcf_dict['alt']
    stash_end = end
    # Re-Analyse genomic positions
    if re.match('NG_', str(stash_input)):
        c = hgvs_tx
        try:
            c.posedit.edit.ref = c.posedit.edit.ref.upper()
        except Exception:
            pass
        try:
            c.posedit.edit.alt = c.posedit.edit.alt.upper()
        except Exception:
            pass
        stash_input = va_func.myevm_t_to_g(c, hdp, no_norm_evm, primary_assembly, vm, hp, hn, sf, nr_vm,
                                           utilise_gap_code)
    if re.match('NC_', str(stash_input)) or re.match('NT_', str(stash_input)) or re.match('NW_',
                                                                                          str(
                                                                                              stash_input)):
        try:
            hgvs_stash = hp.parse_hgvs_variant(stash_input)
        except:
            hgvs_stash = stash_input
        try:
            hgvs_stash.posedit.edit.ref = hgvs_stash.posedit.edit.ref.upper()
        except Exception:
            pass
        try:
            hgvs_stash.posedit.edit.alt = hgvs_stash.posedit.edit.alt.upper()
        except Exception:
            pass

        stash_ac = hgvs_stash.ac
        # MAKE A NO NORM HGVS2VCF
        stash_dict = va_H2V.pos_lock_hgvs2vcf(hgvs_stash, primary_assembly, reverse_normalizer, sf)

        stash_ac = hgvs_stash.ac
        stash_pos = int(stash_dict['pos'])
        stash_ref = stash_dict['ref']
        stash_alt = stash_dict['alt']
        # Generate an end position
        stash_end = str(stash_pos + len(stash_ref) - 1)

    # Store a not real deletion insertion
    stored_hgvs_not_delins = hp.parse_hgvs_variant(str(
        hgvs_genomic_5pr.ac) + ':' + hgvs_genomic_5pr.type + '.' + pos + '_' + end + 'del' + ref + 'ins' + alt)
    stash_hgvs_not_delins = hp.parse_hgvs_variant(
        stash_ac + ':' + hgvs_genomic_5pr.type + '.' + str(
            stash_pos) + '_' + stash_end + 'del' + stash_ref + 'ins' + stash_alt)

    # Set non-valid caution to false
    non_valid_caution = 'false'

    # Store the current hgvs:c. description
    saved_hgvs_coding = hgvs_tx

    # Get orientation of the gene wrt genome and a list of exons mapped to the genome
    ori = va_func.tx_exons(tx_ac=saved_hgvs_coding.ac, alt_ac=hgvs_genomic_5pr.ac,
                           alt_aln_method=alt_aln_method, hdp=hdp)
    orientation = int(ori[0]['alt_strand'])
    intronic_variant = 'false'

    if orientation == -1:
        # position genomic at its most 5 prime position
        try:
            query_genomic = reverse_normalizer.normalize(hgvs_genomic)
        except:
            query_genomic = hgvs_genomic
        # Map to the transcript ant test for movement
        try:
            hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
        except hgvs.exceptions.HGVSError as e:
            hgvs_seek_var = saved_hgvs_coding
        else:
            seek_var = str(hgvs_seek_var)
            seek_ac = str(hgvs_seek_var.ac)
        if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                saved_hgvs_coding.posedit.pos.start.base + saved_hgvs_coding.posedit.pos.start.offset) and (
                hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                saved_hgvs_coding.posedit.pos.end.base + saved_hgvs_coding.posedit.pos.end.offset):
            pass
        else:
            hgvs_seek_var = saved_hgvs_coding

    elif orientation != -1:
        # position genomic at its most 3 prime position
        try:
            query_genomic = hn.normalize(hgvs_genomic)
        except:
            query_genomic = hgvs_genomic
    # Map to the transcript and test for movement
    try:
        hgvs_seek_var = evm.g_to_t(query_genomic, saved_hgvs_coding.ac)
    except hgvs.exceptions.HGVSError as e:
        hgvs_seek_var = saved_hgvs_coding
    else:
        # seek_var = str(hgvs_seek_var)
        # seek_ac = str(hgvs_seek_var.ac)
        if (hgvs_seek_var.posedit.pos.start.base + hgvs_seek_var.posedit.pos.start.offset) > (
                saved_hgvs_coding.posedit.pos.start.base + saved_hgvs_coding.posedit.pos.start.offset) and (
                hgvs_seek_var.posedit.pos.end.base + hgvs_seek_var.posedit.pos.end.offset) > (
                saved_hgvs_coding.posedit.pos.end.base + saved_hgvs_coding.posedit.pos.end.offset):
            pass
        else:
            hgvs_seek_var = saved_hgvs_coding

    try:
        intron_test = hn.normalize(hgvs_seek_var)
    except hgvs.exceptions.HGVSUnsupportedOperationError as e:
        error = str(e)
        if re.match('Normalization of intronic variants is not supported', error) or re.match(
                'Unsupported normalization of variants spanning the exon-intron boundary',
                error):
            if re.match(
                    'Unsupported normalization of variants spanning the exon-intron boundary',
                    error):
                intronic_variant = 'hard_fail'
            else:
                # Double check to see whether the variant is actually intronic?
                for exon in ori:
                    genomic_start = int(exon['alt_start_i'])
                    genomic_end = int(exon['alt_end_i'])
                    if (
                            hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                            hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                        intronic_variant = 'false'
                        break
                    else:
                        intronic_variant = 'true'

    if intronic_variant != 'hard_fail':
        if re.search('\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\d+\-', str(
                hgvs_seek_var.posedit.pos)) or re.search('\*\d+\+', str(
            hgvs_seek_var.posedit.pos)) or re.search('\*\d+\-',
                                                     str(hgvs_seek_var.posedit.pos)):
            # Double check to see whether the variant is actually intronic?
            for exon in ori:
                genomic_start = int(exon['alt_start_i'])
                genomic_end = int(exon['alt_end_i'])
                if (
                        hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                        hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                    intronic_variant = 'false'
                    break
                else:
                    intronic_variant = 'true'

    if re.search('\d+\+', str(hgvs_seek_var.posedit.pos)) or re.search('\d+\-', str(
            hgvs_seek_var.posedit.pos)) or re.search('\*\d+\+', str(
        hgvs_seek_var.posedit.pos)) or re.search('\*\d+\-', str(hgvs_seek_var.posedit.pos)):
        # Double check to see whether the variant is actually intronic?
        for exon in ori:
            genomic_start = int(exon['alt_start_i'])
            genomic_end = int(exon['alt_end_i'])
            if (
                    hgvs_genomic_5pr.posedit.pos.start.base > genomic_start and hgvs_genomic_5pr.posedit.pos.start.base <= genomic_end) and (
                    hgvs_genomic_5pr.posedit.pos.end.base > genomic_start and hgvs_genomic_5pr.posedit.pos.end.base <= genomic_end):
                intronic_variant = 'false'
                break
            else:
                intronic_variant = 'true'

    # If exonic, process
    if intronic_variant != 'true':
        # map form reverse normalized g. to c.
        hgvs_from_5n_g = no_norm_evm.g_to_t(hgvs_genomic_5pr, saved_hgvs_coding.ac)

        # Attempt to find gaps in reference sequence by catching disparity in genome length and overlapping transcript lengths
        disparity_deletion_in = ['false', 'false']
        if stored_hgvs_not_delins != '':
            # Refresh hgvs_not_delins from stored_hgvs_not_delins
            hgvs_not_delins = copy.deepcopy(stored_hgvs_not_delins)
            # This test will only occur in dup of single base, insertion or substitution
            if not re.search('_', str(hgvs_not_delins.posedit.pos)):
                if re.search('dup', hgvs_genomic_5pr.posedit.edit.type) or re.search('ins',
                                                                                     hgvs_genomic_5pr.posedit.edit.type):
                    # For gap in chr, map to t. - but becaouse we have pushed to 5 prime by norm, add 1 to end pos
                    plussed_hgvs_not_delins = copy.deepcopy(hgvs_not_delins)
                    plussed_hgvs_not_delins.posedit.pos.end.base = plussed_hgvs_not_delins.posedit.pos.end.base + 1
                    plussed_hgvs_not_delins.posedit.edit.ref = ''
                    transcript_variant = no_norm_evm.g_to_t(plussed_hgvs_not_delins,
                                                            str(saved_hgvs_coding.ac))
                    if ((
                            transcript_variant.posedit.pos.end.base - transcript_variant.posedit.pos.start.base) > (
                            hgvs_genomic_5pr.posedit.pos.end.base - hgvs_genomic_5pr.posedit.pos.start.base)):
                        if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                            start = hgvs_not_delins.posedit.pos.start.base - 1
                            end = hgvs_not_delins.posedit.pos.end.base
                            ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                            hgvs_not_delins.posedit.edit.ref = ref_bases
                            hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                               :1] + hgvs_not_delins.posedit.edit.alt[
                                                                     1:] + ref_bases[1:]
                        elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search(
                                'del', str(hgvs_genomic_5pr.posedit.edit)):
                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                        elif re.search('ins',
                                       str(hgvs_genomic_5pr.posedit.edit)) and not re.search(
                            'del', str(hgvs_genomic_5pr.posedit.edit)):
                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                            start = hgvs_not_delins.posedit.pos.start.base - 1
                            end = hgvs_not_delins.posedit.pos.end.base
                            ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                            hgvs_not_delins.posedit.edit.ref = ref_bases
                            hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                               :1] + hgvs_not_delins.posedit.edit.alt[
                                                                     1:] + ref_bases[1:]
                    else:
                        if re.search('dup', str(hgvs_genomic_5pr.posedit.edit)):
                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                            start = hgvs_not_delins.posedit.pos.start.base - 1
                            end = hgvs_not_delins.posedit.pos.end.base
                            ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                            hgvs_not_delins.posedit.edit.ref = ref_bases
                            hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                               :1] + hgvs_not_delins.posedit.edit.alt[
                                                                     1:] + ref_bases[1:]
                        elif re.search('ins', str(hgvs_genomic_5pr.posedit.edit)) and re.search(
                                'del', str(hgvs_genomic_5pr.posedit.edit)):
                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                        elif re.search('ins',
                                       str(hgvs_genomic_5pr.posedit.edit)) and not re.search(
                            'del', str(hgvs_genomic_5pr.posedit.edit)):
                            hgvs_not_delins.posedit.pos.end.base = hgvs_not_delins.posedit.pos.start.base + 1
                            start = hgvs_not_delins.posedit.pos.start.base - 1
                            end = hgvs_not_delins.posedit.pos.end.base
                            ref_bases = sf.fetch_seq(str(hgvs_not_delins.ac), start, end)
                            hgvs_not_delins.posedit.edit.ref = ref_bases
                            hgvs_not_delins.posedit.edit.alt = ref_bases[
                                                               :1] + hgvs_not_delins.posedit.edit.alt[
                                                                     1:] + ref_bases[1:]
                else:
                    pass
            else:
                pass

            try:
                tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins, saved_hgvs_coding.ac)
            except hgvs.exceptions.HGVSInvalidIntervalError:
                tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_genomic_5pr, saved_hgvs_coding.ac)
            except hgvs.exceptions.HGVSError:
                if str(e) == 'start or end or both are beyond the bounds of transcript record':
                    tx_hgvs_not_delins = saved_hgvs_coding

            # Create normalized version of tx_hgvs_not_delins
            rn_tx_hgvs_not_delins = copy.deepcopy(tx_hgvs_not_delins)
            # Check for +ve base and adjust
            if (re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)) or re.search('\-',
                                                                                           str(
                                                                                               rn_tx_hgvs_not_delins.posedit.pos.start))) and (
                    re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)) or re.search(
                '\-', str(rn_tx_hgvs_not_delins.posedit.pos.end))):
                # Remove offsetting to span the gap
                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                try:
                    rn_tx_hgvs_not_delins.posedit.edit.alt = ''
                except:
                    pass
            elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                # move tx end base to next available non-offset base
                rn_tx_hgvs_not_delins.posedit.pos.end.base = tx_hgvs_not_delins.posedit.pos.end.base + 1
                rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                    test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                else:
                    test_tx_var = rn_tx_hgvs_not_delins
                # re-make genomic and tx
                hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm, primary_assembly, vm, hp, hn, sf,
                                                       nr_vm, utilise_gap_code)
                rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                           str(saved_hgvs_coding.ac))
            elif re.search('\+', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                # move tx start base to previous available non-offset base
                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                    test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                else:
                    test_tx_var = rn_tx_hgvs_not_delins
                # re-make genomic and tx
                hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm, primary_assembly, vm, hp, hn, sf,
                                                       nr_vm, utilise_gap_code)
                rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                           str(saved_hgvs_coding.ac))
                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
            #                                     else:
            #                                         pass

            # Check for -ve base and adjust
            elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)) and re.search('\-',
                                                                                           str(
                                                                                               rn_tx_hgvs_not_delins.posedit.pos.start)):
                # Remove offsetting to span the gap
                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                rn_tx_hgvs_not_delins.posedit.pos.end.base = rn_tx_hgvs_not_delins.posedit.pos.end.base + 1
                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                try:
                    rn_tx_hgvs_not_delins.posedit.edit.alt = ''
                except:
                    pass
            elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.end)):
                # move tx end base back to next available non-offset base
                rn_tx_hgvs_not_delins.posedit.pos.end.offset = 0
                # Delete the ref
                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                # Add the additional base to the ALT
                start = rn_tx_hgvs_not_delins.posedit.pos.end.base - 1
                end = rn_tx_hgvs_not_delins.posedit.pos.end.base
                ref_bases = sf.fetch_seq(str(tx_hgvs_not_delins.ac), start, end)
                rn_tx_hgvs_not_delins.posedit.edit.alt = rn_tx_hgvs_not_delins.posedit.edit.alt + ref_bases
                if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                    test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                else:
                    test_tx_var = rn_tx_hgvs_not_delins
                # re-make genomic and tx
                hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm, primary_assembly, vm, hp, hn, sf,
                                                       nr_vm, utilise_gap_code)
                rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                           str(saved_hgvs_coding.ac))
            elif re.search('\-', str(rn_tx_hgvs_not_delins.posedit.pos.start)):
                # move tx start base to previous available non-offset base
                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
                rn_tx_hgvs_not_delins.posedit.pos.start.base = rn_tx_hgvs_not_delins.posedit.pos.start.base - 1
                rn_tx_hgvs_not_delins.posedit.edit.ref = ''
                if re.match('NM_', str(rn_tx_hgvs_not_delins)):
                    test_tx_var = no_norm_evm.n_to_c(rn_tx_hgvs_not_delins)
                else:
                    test_tx_var = rn_tx_hgvs_not_delins
                # re-make genomic and tx
                hgvs_not_delins = va_func.myevm_t_to_g(test_tx_var, hdp, no_norm_evm, primary_assembly, vm, hp, hn, sf,
                                                       nr_vm, utilise_gap_code)
                rn_tx_hgvs_not_delins = no_norm_evm.g_to_n(hgvs_not_delins,
                                                           str(saved_hgvs_coding.ac))
                rn_tx_hgvs_not_delins.posedit.pos.start.offset = 0
            else:
                pass

            # Logic
            if len(hgvs_not_delins.posedit.edit.ref) < len(
                    rn_tx_hgvs_not_delins.posedit.edit.ref):
                gap_length = len(rn_tx_hgvs_not_delins.posedit.edit.ref) - len(
                    hgvs_not_delins.posedit.edit.ref)
                disparity_deletion_in = ['chromosome', gap_length]
            elif len(hgvs_not_delins.posedit.edit.ref) > len(
                    rn_tx_hgvs_not_delins.posedit.edit.ref):
                gap_length = len(hgvs_not_delins.posedit.edit.ref) - len(
                    rn_tx_hgvs_not_delins.posedit.edit.ref)
                disparity_deletion_in = ['transcript', gap_length]
            else:
                hgvs_stash_t = vm.g_to_t(stash_hgvs_not_delins, saved_hgvs_coding.ac)
                if len(stash_hgvs_not_delins.posedit.edit.ref) > len(
                        hgvs_stash_t.posedit.edit.ref):
                    try:
                        hn.normalize(hgvs_stash_t)
                    except:
                        pass
                    else:
                        gap_length = len(stash_hgvs_not_delins.posedit.edit.ref) - len(
                            hgvs_stash_t.posedit.edit.ref)
                        disparity_deletion_in = ['transcript', gap_length]
                        try:
                            tx_hgvs_not_delins = vm.c_to_n(hgvs_stash_t)
                        except:
                            tx_hgvs_not_delins = hgvs_stash_t
                        hgvs_not_delins = stash_hgvs_not_delins
                elif hgvs_stash_t.posedit.pos.start.offset != 0 or hgvs_stash_t.posedit.pos.end.offset != 0:
                    disparity_deletion_in = ['transcript', 'Requires Analysis']
                    try:
                        tx_hgvs_not_delins = vm.c_to_n(hgvs_stash_t)
                    except:
                        tx_hgvs_not_delins = hgvs_stash_t
                    hgvs_not_delins = stash_hgvs_not_delins
                    hgvs_genomic_5pr = stash_hgvs_not_delins
                else:
                    pass

        # Final sanity checks
        try:
            vm.g_to_t(hgvs_not_delins, tx_hgvs_not_delins.ac)
        except Exception as e:
            if str(e) == 'start or end or both are beyond the bounds of transcript record':
                hgvs_not_delins = saved_hgvs_coding
                disparity_deletion_in = ['false', 'false']
        try:
            hn.normalize(tx_hgvs_not_delins)
        except hgvs.exceptions.HGVSUnsupportedOperationError as e:
            error = str(e)
            if re.match('Normalization of intronic variants is not supported',
                        error) or re.match(
                'Unsupported normalization of variants spanning the exon-intron boundary',
                error):
                if re.match(
                        'Unsupported normalization of variants spanning the exon-intron boundary',
                        error):
                    hgvs_not_delins = saved_hgvs_coding
                    disparity_deletion_in = ['false', 'false']
                elif re.match('Normalization of intronic variants is not supported', error):
                    # We know that this cannot be because of an intronic variant, so must be aligned to tx gap
                    disparity_deletion_in = ['transcript', 'Requires Analysis']

        # Pre-processing of tx_hgvs_not_delins
        try:
            if tx_hgvs_not_delins.posedit.edit.alt is None:
                tx_hgvs_not_delins.posedit.edit.alt = ''
        except Exception as e:
            if str(e) == "'Dup' object has no attribute 'alt'":
                tx_hgvs_not_delins_delins_from_dup = tx_hgvs_not_delins.ac + ':' + tx_hgvs_not_delins.type + '.' + str(
                    tx_hgvs_not_delins.posedit.pos.start) + '_' + str(
                    tx_hgvs_not_delins.posedit.pos.end) + 'del' + tx_hgvs_not_delins.posedit.edit.ref + 'ins' + tx_hgvs_not_delins.posedit.edit.ref + tx_hgvs_not_delins.posedit.edit.ref
                tx_hgvs_not_delins = hp.parse_hgvs_variant(tx_hgvs_not_delins_delins_from_dup)

        # GAP IN THE TRANSCRIPT DISPARITY DETECTED
        if disparity_deletion_in[0] == 'transcript':
            if disparity_deletion_in[1] == 'Requires Analysis':
                analyse_gap = copy.deepcopy(tx_hgvs_not_delins)
                try:
                    analyse_gap = vm.n_to_c(analyse_gap)
                except hgvs.exceptions.HGVSError:
                    pass
                analyse_gap.posedit.pos.start.offset = 0
                analyse_gap.posedit.pos.end.offset = 0
                try:
                    analyse_gap.posedit.edit.ref = ''
                except AttributeError:
                    pass
                try:
                    analyse_gap.posedit.edit.alt = ''
                except AttributeError:
                    pass
                my_g_gap_v = vm.t_to_g(analyse_gap, hgvs_genomic_5pr.ac)
                g_gap = my_g_gap_v.posedit.pos.end.base - my_g_gap_v.posedit.pos.start.base + 1
                my_t_gap_v = vm.g_to_t(my_g_gap_v, analyse_gap.ac)
                t_gap = my_t_gap_v.posedit.pos.end.base - my_t_gap_v.posedit.pos.start.base
                disparity_deletion_in[1] = str(g_gap - t_gap - 1)

            gap_position = ''
            gapped_alignment_warning = str(
                hgvs_genomic_5pr) + ' may be an artefact of aligning ' + tx_hgvs_not_delins.ac + ' with genome build ' + primary_assembly

            # ANY VARIANT WHOLLY WITHIN THE GAP
            if (re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) or re.search('\-',
                                                                                        str(
                                                                                            tx_hgvs_not_delins.posedit.pos.start))) and (
                    re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) or re.search('\-',
                                                                                          str(
                                                                                              tx_hgvs_not_delins.posedit.pos.end))):
                gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
                # Copy the current variant
                tx_gap_fill_variant = copy.deepcopy(tx_hgvs_not_delins)
                try:
                    if tx_gap_fill_variant.posedit.edit.alt is None:
                        tx_gap_fill_variant.posedit.edit.alt = ''
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        tx_gap_fill_variant_delins_from_dup = tx_gap_fill_variant.ac + ':' + tx_gap_fill_variant.type + '.' + str(
                            tx_gap_fill_variant.posedit.pos.start) + '_' + str(
                            tx_gap_fill_variant.posedit.pos.end) + 'del' + tx_gap_fill_variant.posedit.edit.ref + 'ins' + tx_gap_fill_variant.posedit.edit.ref + tx_gap_fill_variant.posedit.edit.ref
                        tx_gap_fill_variant = hp.parse_hgvs_variant(
                            tx_gap_fill_variant_delins_from_dup)

                # Identify which half of the NOT-intron the start position of the variant is in
                if re.search('\-', str(tx_gap_fill_variant.posedit.pos.start)):
                    tx_gap_fill_variant.posedit.pos.start.base = tx_gap_fill_variant.posedit.pos.start.base - 1
                    tx_gap_fill_variant.posedit.pos.start.offset = int('0')  # int('+1')
                    tx_gap_fill_variant.posedit.pos.end.offset = int('0')  # int('-1')
                    tx_gap_fill_variant.posedit.edit.alt = ''
                    tx_gap_fill_variant.posedit.edit.ref = ''
                elif re.search('\+', str(tx_gap_fill_variant.posedit.pos.start)):
                    tx_gap_fill_variant.posedit.pos.start.offset = int('0')  # int('+1')
                    tx_gap_fill_variant.posedit.pos.end.base = tx_gap_fill_variant.posedit.pos.end.base + 1
                    tx_gap_fill_variant.posedit.pos.end.offset = int('0')  # int('-1')
                    tx_gap_fill_variant.posedit.edit.alt = ''
                    tx_gap_fill_variant.posedit.edit.ref = ''

                try:
                    tx_gap_fill_variant = vm.n_to_c(tx_gap_fill_variant)
                except:
                    pass
                genomic_gap_fill_variant = vm.t_to_g(tx_gap_fill_variant,
                                                     reverse_normalized_hgvs_genomic.ac)
                genomic_gap_fill_variant.posedit.edit.alt = genomic_gap_fill_variant.posedit.edit.ref

                try:
                    c_tx_hgvs_not_delins = vm.n_to_c(tx_hgvs_not_delins)
                except Exception:
                    c_tx_hgvs_not_delins = copy.copy(tx_hgvs_not_delins)
                genomic_gap_fill_variant_alt = vm.t_to_g(c_tx_hgvs_not_delins,
                                                         hgvs_genomic_5pr.ac)

                # Ensure an ALT exists
                try:
                    if genomic_gap_fill_variant_alt.posedit.edit.alt is None:
                        genomic_gap_fill_variant_alt.posedit.edit.alt = 'X'
                except Exception as e:
                    if str(e) == "'Dup' object has no attribute 'alt'":
                        genomic_gap_fill_variant_delins_from_dup = genomic_gap_fill_variant.ac + ':' + genomic_gap_fill_variant.type + '.' + str(
                            genomic_gap_fill_variant.posedit.pos.start.base) + '_' + str(
                            genomic_gap_fill_variant.posedit.pos.end.base) + 'del' + genomic_gap_fill_variant.posedit.edit.ref + 'ins' + genomic_gap_fill_variant.posedit.edit.ref + genomic_gap_fill_variant.posedit.edit.ref
                        genomic_gap_fill_variant = hp.parse_hgvs_variant(
                            genomic_gap_fill_variant_delins_from_dup)
                        genomic_gap_fill_variant_alt_delins_from_dup = genomic_gap_fill_variant_alt.ac + ':' + genomic_gap_fill_variant_alt.type + '.' + str(
                            genomic_gap_fill_variant_alt.posedit.pos.start.base) + '_' + str(
                            genomic_gap_fill_variant_alt.posedit.pos.end.base) + 'del' + genomic_gap_fill_variant_alt.posedit.edit.ref + 'ins' + genomic_gap_fill_variant_alt.posedit.edit.ref + genomic_gap_fill_variant_alt.posedit.edit.ref
                        genomic_gap_fill_variant_alt = hp.parse_hgvs_variant(
                            genomic_gap_fill_variant_alt_delins_from_dup)

                # Correct insertion alts
                if genomic_gap_fill_variant_alt.posedit.edit.type == 'ins':
                    append_ref = sf.fetch_seq(genomic_gap_fill_variant_alt.ac,
                                              genomic_gap_fill_variant_alt.posedit.pos.start.base - 1,
                                              genomic_gap_fill_variant_alt.posedit.pos.end.base)
                    genomic_gap_fill_variant_alt.posedit.edit.alt = append_ref[
                                                                        0] + genomic_gap_fill_variant_alt.posedit.edit.alt + \
                                                                    append_ref[1]

                # Split the reference and replacing alt sequence into a dictionary
                reference_bases = list(genomic_gap_fill_variant.posedit.edit.ref)
                if genomic_gap_fill_variant_alt.posedit.edit.alt is not None:
                    alternate_bases = list(genomic_gap_fill_variant_alt.posedit.edit.alt)
                else:
                    # Deletions with no ins
                    pre_alternate_bases = list(genomic_gap_fill_variant_alt.posedit.edit.ref)
                    alternate_bases = []
                    for base in pre_alternate_bases:
                        alternate_bases.append('X')

                # Create the dictionaries
                ref_start = genomic_gap_fill_variant.posedit.pos.start.base
                alt_start = genomic_gap_fill_variant_alt.posedit.pos.start.base
                ref_base_dict = {}
                for base in reference_bases:
                    ref_base_dict[ref_start] = str(base)
                    ref_start = ref_start + 1

                alt_base_dict = {}

                # NEED TO SEARCH FOR RANGE = and replace with interval_range
                # Need to search for int and replace with integer

                # Note, all variants will be forced into the format delete insert
                # Deleted bases in the ALT will be substituted for X
                for integer in range(genomic_gap_fill_variant_alt.posedit.pos.start.base,
                                     genomic_gap_fill_variant_alt.posedit.pos.end.base + 1, 1):
                    if integer == alt_start:
                        alt_base_dict[integer] = str(''.join(alternate_bases))
                    else:
                        alt_base_dict[integer] = 'X'

                # Generate the alt sequence
                alternate_sequence_bases = []
                for integer in range(genomic_gap_fill_variant.posedit.pos.start.base,
                                     genomic_gap_fill_variant.posedit.pos.end.base + 1, 1):
                    if integer in alt_base_dict.keys():
                        alternate_sequence_bases.append(alt_base_dict[integer])
                    else:
                        alternate_sequence_bases.append(ref_base_dict[integer])
                alternate_sequence = ''.join(alternate_sequence_bases)
                alternate_sequence = alternate_sequence.replace('X', '')

                # Add the new alt to the gap fill variant and generate transcript variant
                genomic_gap_fill_variant.posedit.edit.alt = alternate_sequence
                hgvs_refreshed_variant = vm.g_to_t(genomic_gap_fill_variant,
                                                   tx_gap_fill_variant.ac)

                # Set warning
                gap_size = str(len(genomic_gap_fill_variant.posedit.edit.ref) - 2)
                disparity_deletion_in[1] = [gap_size]
                auto_info = str(stored_hgvs_not_delins.ac) + ':g.' + str(
                    stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + gap_size + ' genomic base(s) that fail to align to transcript ' + str(
                    tx_hgvs_not_delins.ac)
                non_valid_caution = 'true'

                # Alignment position
                for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                if re.match('NM_', str(for_location_c)):
                    for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                if re.match('\-', str(for_location_c.posedit.pos.start.offset)):
                    gps = for_location_c.posedit.pos.start.base - 1
                    gpe = for_location_c.posedit.pos.start.base
                else:
                    gps = for_location_c.posedit.pos.start.base
                    gpe = for_location_c.posedit.pos.start.base + 1
                gap_position = ' between positions c.' + str(gps) + '_' + str(gpe)
                auto_info = '%s' % (gap_position)
            else:
                if tx_hgvs_not_delins.posedit.pos.start.offset == 0 and tx_hgvs_not_delins.posedit.pos.end.offset == 0:
                    # In this instance, we have identified a transcript gap but the n. version of
                    # the transcript variant but do not have a position which actually hits the gap,
                    # so the variant likely spans the gap, and is not picked up by an offset.
                    try:
                        c1 = vm.n_to_c(tx_hgvs_not_delins)
                    except:
                        c1 = tx_hgvs_not_delins
                    g1 = nr_vm.t_to_g(c1, hgvs_genomic.ac)
                    g3 = nr_vm.t_to_g(c1, hgvs_genomic.ac)
                    g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                    ng2 = hn.normalize(g2)
                    g3.posedit.pos.end.base = g3.posedit.pos.start.base + (
                            len(g3.posedit.edit.ref) - 1)
                    try:
                        c2 = vm.g_to_t(g3, c1.ac)
                        if c2.posedit.pos.start.offset == 0 and c2.posedit.pos.end.offset == 0:
                            pass
                        else:
                            tx_hgvs_not_delins = c2
                            try:
                                tx_hgvs_not_delins = vm.c_to_n(tx_hgvs_not_delins)
                            except hgvs.exceptions.HGVSError:
                                pass
                    except hgvs.exceptions.HGVSInvalidVariantError:
                        pass

                if re.search('\+', str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                        '\+', str(tx_hgvs_not_delins.posedit.pos.end)):
                    auto_info = str(stored_hgvs_not_delins.ac) + ':g.' + str(
                        stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                        disparity_deletion_in[
                            1]) + ' genomic base(s) that fail to align to transcript ' + str(
                        tx_hgvs_not_delins.ac)
                    non_valid_caution = 'true'
                    try:
                        c2 = vm.n_to_c(tx_hgvs_not_delins)
                    except:
                        c2 = tx_hgvs_not_delins
                    c1 = copy.deepcopy(c2)
                    c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                    c1.posedit.pos.start.offset = 0
                    c1.posedit.pos.end = c2.posedit.pos.start
                    c1.posedit.edit.ref = ''
                    c1.posedit.edit.alt = ''
                    if orientation != -1:
                        g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                        g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                        g1.posedit.edit.alt = g1.posedit.edit.ref
                    else:
                        g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                        g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                        g2.posedit.edit.alt = g2.posedit.edit.ref
                    reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                    alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                    g3 = copy.deepcopy(g1)
                    g3.posedit.pos.end.base = g2.posedit.pos.end.base
                    g3.posedit.edit.ref = reference
                    g3.posedit.edit.alt = alternate
                    c3 = vm.g_to_t(g3, c1.ac)
                    hgvs_refreshed_variant = c3
                    # Alignment position
                    for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                    if re.match('NM_', str(for_location_c)):
                        for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                        gps = for_location_c.posedit.pos.start.base
                        gpe = for_location_c.posedit.pos.start.base + 1
                    gap_position = ' between positions c.' + str(gps) + '_' + str(gpe)
                    # Warn update
                    auto_info = '%s' % (gap_position)
                elif re.search('\+', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search(
                        '\+', str(tx_hgvs_not_delins.posedit.pos.start)):
                    auto_info = 'Genome position ' + str(
                        stored_hgvs_not_delins.ac) + ':g.' + str(
                        stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(
                        disparity_deletion_in[1]) + '-bp gap in transcript ' + str(
                        tx_hgvs_not_delins.ac)
                    gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
                    non_valid_caution = 'true'
                    try:
                        c1 = vm.n_to_c(tx_hgvs_not_delins)
                    except:
                        c1 = tx_hgvs_not_delins
                    c2 = copy.deepcopy(c1)
                    c2.posedit.pos.start = c1.posedit.pos.end
                    c2.posedit.pos.end.base = c1.posedit.pos.end.base + 1
                    c2.posedit.pos.end.offset = 0
                    c2.posedit.edit.ref = ''
                    c2.posedit.edit.alt = ''
                    if orientation != -1:
                        g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                        g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                        g2.posedit.edit.alt = g2.posedit.edit.ref
                    else:
                        g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                        g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                        g1.posedit.edit.alt = g1.posedit.edit.ref
                    reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                    alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                    g3 = copy.deepcopy(g1)
                    g3.posedit.pos.end.base = g2.posedit.pos.end.base
                    g3.posedit.edit.ref = reference
                    g3.posedit.edit.alt = alternate
                    c3 = vm.g_to_t(g3, c1.ac)
                    hgvs_refreshed_variant = c3
                    # Alignment position
                    for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                    if re.match('NM_', str(for_location_c)):
                        for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                    gps = for_location_c.posedit.pos.end.base
                    gpe = for_location_c.posedit.pos.end.base + 1
                    gap_position = ' between positions c.' + str(gps) + '_' + str(gpe)
                    # Warn update
                    auto_info = '%s' % (gap_position)
                elif re.search('\-',
                               str(tx_hgvs_not_delins.posedit.pos.start)) and not re.search(
                    '\-', str(tx_hgvs_not_delins.posedit.pos.end)):
                    auto_info = str(stored_hgvs_not_delins.ac) + ':g.' + str(
                        stored_hgvs_not_delins.posedit.pos.start.base) + ' is one of ' + str(
                        disparity_deletion_in[
                            1]) + ' genomic base(s) that fail to align to transcript ' + str(
                        tx_hgvs_not_delins.ac)
                    non_valid_caution = 'true'
                    try:
                        c2 = vm.n_to_c(tx_hgvs_not_delins)
                    except:
                        c2 = tx_hgvs_not_delins
                    c1 = copy.deepcopy(c2)
                    c1.posedit.pos.start.base = c2.posedit.pos.start.base - 1
                    c1.posedit.pos.start.offset = 0
                    c1.posedit.pos.end = c2.posedit.pos.start
                    c1.posedit.edit.ref = ''
                    c1.posedit.edit.alt = ''
                    if orientation != -1:
                        g1 = vm.t_to_g(c1, hgvs_genomic.ac)
                        g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                        g1.posedit.edit.alt = g1.posedit.edit.ref
                    else:
                        g1 = vm.t_to_g(c2, hgvs_genomic.ac)
                        g2 = vm.t_to_g(c1, hgvs_genomic.ac)
                        g2.posedit.edit.alt = g2.posedit.edit.ref
                    reference = g1.posedit.edit.ref + g2.posedit.edit.ref[1:]
                    alternate = g1.posedit.edit.alt + g2.posedit.edit.alt[1:]
                    g3 = copy.deepcopy(g1)
                    g3.posedit.pos.end.base = g2.posedit.pos.end.base
                    g3.posedit.edit.ref = reference
                    g3.posedit.edit.alt = alternate
                    c3 = vm.g_to_t(g3, c1.ac)
                    hgvs_refreshed_variant = c3
                    # Alignment position
                    for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                    if re.match('NM_', str(for_location_c)):
                        for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                    gps = for_location_c.posedit.pos.start.base - 1
                    gpe = for_location_c.posedit.pos.start.base
                    gap_position = ' between positions c.' + str(gps) + '_' + str(gpe)
                    # Warn update
                    auto_info = '%s' % (gap_position)
                elif re.search('\-', str(tx_hgvs_not_delins.posedit.pos.end)) and not re.search(
                        '\-', str(tx_hgvs_not_delins.posedit.pos.start)):
                    auto_info = 'Genome position ' + str(
                        stored_hgvs_not_delins.ac) + ':g.' + str(
                        stored_hgvs_not_delins.posedit.pos.end.base + 1) + ' aligns within a ' + str(
                        disparity_deletion_in[1]) + '-bp gap in transcript ' + str(
                        tx_hgvs_not_delins.ac)
                    gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)
                    non_valid_caution = 'true'
                    try:
                        c1 = vm.n_to_c(tx_hgvs_not_delins)
                    except:
                        c1 = tx_hgvs_not_delins
                    c2 = copy.deepcopy(c1)
                    c2.posedit.pos.start = c1.posedit.pos.end
                    c2.posedit.pos.end.base = c1.posedit.pos.end.base
                    c2.posedit.pos.end.offset = 0
                    c2.posedit.edit.ref = ''
                    c2.posedit.edit.alt = ''
                    g2 = vm.t_to_g(c2, hgvs_genomic.ac)
                    c2 = vm.g_to_t(g2, c2.ac)
                    reference = c1.posedit.edit.ref + c2.posedit.edit.ref[1:]
                    alternate = c1.posedit.edit.alt + c2.posedit.edit.ref[1:]
                    c3 = copy.deepcopy(c1)
                    c3.posedit.pos.end = c2.posedit.pos.end
                    c3.posedit.edit.ref = ''  # reference
                    c3.posedit.edit.alt = alternate
                    hgvs_refreshed_variant = c3
                    # Alignment position
                    for_location_c = copy.deepcopy(hgvs_refreshed_variant)
                    if re.match('NM_', str(for_location_c)):
                        for_location_c = no_norm_evm.n_to_c(tx_hgvs_not_delins)
                    gps = for_location_c.posedit.pos.end.base - 1
                    gpe = for_location_c.posedit.pos.end.base
                    gap_position = ' between positions c.' + str(gps) + '_' + str(gpe)
                    # Warn update
                    auto_info = '%s' % (gap_position)
                else:
                    auto_info = str(stored_hgvs_not_delins.ac) + ':g.' + str(
                        stored_hgvs_not_delins.posedit.pos) + ' contains ' + str(
                        disparity_deletion_in[
                            1]) + ' genomic base(s) that fail to align to transcript ' + str(
                        tx_hgvs_not_delins.ac) + '\n'
                    hgvs_refreshed_variant = tx_hgvs_not_delins
                    gapped_transcripts = gapped_transcripts + ' ' + str(tx_hgvs_not_delins.ac)

        # GAP IN THE CHROMOSOME
        elif disparity_deletion_in[0] == 'chromosome':
            # Set warning variables
            gap_position = ''
            gapped_alignment_warning = str(
                hgvs_genomic_5pr) + ' may be an artefact of aligning ' + tx_hgvs_not_delins.ac + ' with genome build ' + primary_assembly
            hgvs_refreshed_variant = tx_hgvs_not_delins
            # Warn
            auto_info = str(hgvs_refreshed_variant.ac) + ':c.' + str(
                hgvs_refreshed_variant.posedit.pos) + ' contains ' + str(disparity_deletion_in[
                                                                             1]) + ' transcript base(s) that fail to align to chromosome ' + str(
                hgvs_genomic.ac) + '\n'
            gapped_transcripts = gapped_transcripts + str(hgvs_refreshed_variant.ac) + ' '
        else:
            # Try the push
            hgvs_stash = copy.deepcopy(stash_hgvs_not_delins)
            stash_ac = hgvs_stash.ac
            # Make a hard left and hard right not delins g.
            stash_dict_right = va_H2V.hard_right_hgvs2vcf(hgvs_stash, primary_assembly, hn, sf)
            stash_pos_right = int(stash_dict_right['pos'])
            stash_ref_right = stash_dict_right['ref']
            stash_alt_right = stash_dict_right['alt']
            stash_end_right = str(stash_pos_right + len(stash_ref_right) - 1)
            stash_hgvs_not_delins_right = hp.parse_hgvs_variant(stash_ac + ':' + hgvs_stash.type + '.' + str(
                stash_pos_right) + '_' + stash_end_right + 'del' + stash_ref_right + 'ins' + stash_alt_right)
            stash_dict_left = va_H2V.hard_left_hgvs2vcf(hgvs_stash, primary_assembly, reverse_normalizer, sf)
            stash_pos_left = int(stash_dict_left['pos'])
            stash_ref_left = stash_dict_left['ref']
            stash_alt_left = stash_dict_left['alt']
            stash_end_left = str(stash_pos_left + len(stash_ref_left) - 1)
            stash_hgvs_not_delins_left = hp.parse_hgvs_variant(stash_ac + ':' + hgvs_stash.type + '.' + str(
                stash_pos_left) + '_' + stash_end_left + 'del' + stash_ref_left + 'ins' + stash_alt_left)
            # Map in-situ to the transcript left and right
            try:
                tx_hard_right = vm.g_to_t(stash_hgvs_not_delins_right, saved_hgvs_coding.ac)
            except Exception as e:
                tx_hard_right = saved_hgvs_coding
            else:
                normalize_stash_right = hn.normalize(stash_hgvs_not_delins_right)
                if str(normalize_stash_right.posedit) == str(stash_hgvs_not_delins.posedit):
                    tx_hard_right = saved_hgvs_coding
            try:
                tx_hard_left = vm.g_to_t(stash_hgvs_not_delins_left, saved_hgvs_coding.ac)
            except Exception as e:
                tx_hard_left = saved_hgvs_coding
            else:
                normalize_stash_left = hn.normalize(stash_hgvs_not_delins_left)
                if str(normalize_stash_left.posedit) == str(stash_hgvs_not_delins.posedit):
                    tx_hard_left = saved_hgvs_coding
            # The Logic - Currently limited to genome gaps
            if len(stash_hgvs_not_delins_right.posedit.edit.ref) < len(
                    tx_hard_right.posedit.edit.ref):
                tx_hard_right = hn.normalize(tx_hard_right)
                gap_position = ''
                gapped_alignment_warning = str(
                    hgvs_genomic_5pr) + ' may be an artefact of aligning ' + tx_hgvs_not_delins.ac + ' with genome build ' + primary_assembly
                hgvs_refreshed_variant = tx_hard_right
                gapped_transcripts = gapped_transcripts + str(tx_hard_right.ac) + ' '
            elif len(stash_hgvs_not_delins_left.posedit.edit.ref) < len(
                    tx_hard_left.posedit.edit.ref):
                tx_hard_left = hn.normalize(tx_hard_left)
                gap_position = ''
                gapped_alignment_warning = str(
                    hgvs_genomic_5pr) + ' may be an artefact of aligning ' + tx_hgvs_not_delins.ac + ' with genome build ' + primary_assembly
                hgvs_refreshed_variant = tx_hard_left
                gapped_transcripts = gapped_transcripts + str(tx_hard_left.ac) + ' '
            else:
                # Keep the same by re-setting rel_var
                hgvs_refreshed_variant = saved_hgvs_coding

        # Edit the output
        if re.match('NM_', str(hgvs_refreshed_variant.ac)) and not re.search('c', str(
                hgvs_refreshed_variant.type)):
            hgvs_refreshed_variant = evm.n_to_c(hgvs_refreshed_variant)
        else:
            pass

        # Set pos_lock
        pos_lock = True

        try:
            hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)
            if hgvs_refreshed_variant.posedit.edit.type == 'delins' and \
                    hgvs_refreshed_variant.posedit.edit.ref[-1] == \
                    hgvs_refreshed_variant.posedit.edit.alt[-1]:
                hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[
                                                          0:-1]
                hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[
                                                          0:-1]
                hgvs_refreshed_variant.posedit.pos.end.base = hgvs_refreshed_variant.posedit.pos.end.base - 1
                hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)
            elif hgvs_refreshed_variant.posedit.edit.type == 'delins' and \
                    hgvs_refreshed_variant.posedit.edit.ref[0] == \
                    hgvs_refreshed_variant.posedit.edit.alt[0]:
                hgvs_refreshed_variant.posedit.edit.ref = hgvs_refreshed_variant.posedit.edit.ref[
                                                          1:]
                hgvs_refreshed_variant.posedit.edit.alt = hgvs_refreshed_variant.posedit.edit.alt[
                                                          1:]
                hgvs_refreshed_variant.posedit.pos.start.base = hgvs_refreshed_variant.posedit.pos.start.base + 1
                hgvs_refreshed_variant = hn.normalize(hgvs_refreshed_variant)

        except Exception as e:
            error = str(e)
            # Ensure the final variant is not intronic nor does it cross exon boundaries
            if re.match('Normalization of intronic variants is not supported',
                        error) or re.match(
                'Unsupported normalization of variants spanning the exon-intron boundary',
                error):
                hgvs_refreshed_variant = saved_hgvs_coding
                corrective_action_taken = None
                gapped_alignment_warning = None
                auto_info = None
                pos_lock = False
            else:
                pass

    # Otherwise these variants need to be set
    else:
        corrective_action_taken = None
        gapped_alignment_warning = None
        auto_info = None
        pos_lock = False
        hgvs_refreshed_variant = saved_hgvs_coding

    # Warn the user that the g. description is not valid
    if gapped_alignment_warning is not None:
        if disparity_deletion_in[0] == 'transcript':
            corrective_action_taken = 'Automap has deleted ' + str(
                disparity_deletion_in[1]) + ' bp from chromosomal reference sequence ' + str(
                hgvs_genomic.ac) + ' to ensure perfect alignment with transcript reference_sequence'
        if disparity_deletion_in[0] == 'chromosome':
            corrective_action_taken = 'Automap has added ' + str(
                disparity_deletion_in[1]) + ' bp to chromosomal reference sequence ' + str(
                hgvs_genomic.ac) + ' to ensure perfect alignment with transcript reference_sequence '

    # Add additional data to the front of automap
    if auto_info is not None:
        automap = auto_info + '\n' + automap

    # Make the return
    gap_report = [hgvs_refreshed_variant, pos_lock, corrective_action_taken,
                  gapped_alignment_warning, auto_info]

    return gap_report

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