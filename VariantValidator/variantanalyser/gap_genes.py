"""
Lists of genes for GRCh37 and GRCh38 which require a gap to be inserted into either the 
transcript or the genome to maintain a perfect alignment
"""
def gap_black_list(symbol):
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