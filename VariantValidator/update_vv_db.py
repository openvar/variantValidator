# -*- coding: utf-8 -*-
import requests
import copy
import gzip
import re
import logging
from configparser import ConfigParser
from .modules import vvDatabase
from . import configure
import VariantValidator
validator = VariantValidator.Validator()

logger = logging.getLogger(__name__)


def connect():
    config = ConfigParser()
    config.read(configure.CONFIG_DIR)

    db_config = {
        'user': config["mysql"]["user"],
        'password': config["mysql"]["password"],
        'host': config["mysql"]["host"],
        'port': int(config["mysql"]["port"]),
        'database': config["mysql"]["database"],
        'raise_on_warnings': True
    }
    logger.debug("Connecting to database with config %s", db_config)
    # Create database access objects
    db = vvDatabase.Database(db_config)
    return db


def delete():

    db = connect()

    db.execute('DELETE FROM transcript_info')
    db.execute('DELETE FROM refSeqGene_loci')
    db.execute('DELETE FROM LRG_transcripts')
    db.execute('DELETE FROM LRG_proteins')
    db.execute('DELETE FROM LRG_RSG_lookup')
    db.conn.commit()
    logger.debug("Deleted data from all tables including transcript_info")


def update():

    db = connect()

    update_refseq(db)
    update_lrg(db)


def update_refseq(dbcnx):
    logger.debug('Updating RefSeqGene no Missmatch MySQL data')

    # Download data from RefSeqGene
    # Download data
    rsg = requests.get('http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/gene_RefSeqGene')
    rsg_data = rsg.text.strip().split('\n')

    # Download data for RefSeqGene to GRCh38 alignments, we no longer have access to the smaller
    # *_refseqgene_alignments.gff3 files so have to use GCF*_GRCh38.*_genomic.gff instead.
    # We don't know what the latest build will be called exactly so navigate towards it first.
    base_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/' + \
            'latest_assembly_versions/'
    ncbi_web_resp = requests.get(base_url)
    ncbi_web_resp.raise_for_status()
    ncbi_web_resp = re.findall(r"href=[\"'][^\"']+[\"']", ncbi_web_resp.text)
    base_url = base_url +list(filter(lambda x: 'GCF_' in x and 'GRCh' in x, ncbi_web_resp))[0][6:-1]
    ncbi_web_resp = requests.get(base_url)
    ncbi_web_resp.raise_for_status()
    ncbi_web_resp = re.findall(r"href=[\"'][^\"']+[\"']", ncbi_web_resp.text)
    ncbi_web_resp = list(filter(
        lambda x: 'GCF_' in x and 'GRCh' in x and '_genomic.gff.gz' in x,
        ncbi_web_resp))
    grch38 = requests.get( base_url + ncbi_web_resp[0][6:-1])
    grch38.raise_for_status()
    grch38 = gzip.decompress(grch38.content)
    grch38 = grch38.decode()
    grch38_align_data = grch38.strip().split('\n')

    # Download data for GRCh37 alignments, this is optional at best, and lacking a current good
    # data-source is hard disabled for now
    #grch37 = requests.get(
    #    'http://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/GCF_000001405.25_refseqgene_alignments.gff3')
    #grch37_align_data = grch37.text.strip().split('\n')
    grch37_align_data = []

    # Open Lists
    # rsg_data = open(os.path.join(ROOT, 'gene_RefSeqGene'), 'r')
    rsg_id_info = []
    # grch37_align_data = open(os.path.join(ROOT, 'GCF_000001405.25_refseqgene_alignments.gff3'), 'r')
    grch37_align = []
    # grch38_align_data = open(os.path.join(ROOT, 'GCF_000001405.28_refseqgene_alignments.gff3'), 'r')
    grch38_align = []

    # Place the required data from each file into a dictionary

    for line in rsg_data:
        if line.startswith('#'):
            continue
        info = line.strip().split()
        if len(info) > 0:
            entry = {'symbol': info[2], 'rsg_id': info[3], 'gene_id': info[1]}
            rsg_id_info.append(entry)

    # Create dictionary to store RefSeqGene and gene symbol data NOTE RefSeqGene ID stored without version number!
    rsg_to_symbol = {}
    # Collect the data
    for ent in rsg_id_info:
        rsg_id = ent['rsg_id'].split('.')[0]
        rsg_to_symbol[rsg_id] = {'symbol': ent['symbol'], 'gene_id': ent['gene_id']}

    # Count total number of NG to NC mappings
    total_rsg_to_nc = 0
    total_rsg_to_nc_rejected = 0
    for line in grch37_align_data:
        ng_nc = count_ng_nc(line)
        if ng_nc is not None:
            total_rsg_to_nc = total_rsg_to_nc + 1
            if ng_nc == 'rejected':
                total_rsg_to_nc_rejected = total_rsg_to_nc_rejected + 1
            elif ng_nc != 'failed':
                grch37_align.append(ng_nc)

    for line in grch38_align_data:
        ng_nc = count_ng_nc(line)
        if ng_nc is not None:
            total_rsg_to_nc = total_rsg_to_nc + 1
            if ng_nc == 'rejected':
                total_rsg_to_nc_rejected = total_rsg_to_nc_rejected + 1
            elif ng_nc != 'failed':
                grch38_align.append(ng_nc)

    # Create a data array containing the database
    db = []
    # map line
    for line in grch37_align:
        ml = map_line(line, 'GRCh37', rsg_id_info)
        db.append(ml)

    for line in grch38_align:
        ml = map_line(line, 'GRCh38', rsg_id_info)
        db.append(ml)

    missing = []

    # Identify lines with missing data e.g. gene symbols
    of_these = len(db)
    counter = 1
    for line in db:
        print("\n")
        print(counter, "of", of_these)
        counter += 1
        update_success = False
        try:
            line[6]
        except IndexError:
            try:
                identifier = copy.deepcopy(line[0])
                identifier = identifier.split('.')[0]
                line.append(rsg_to_symbol[identifier]['symbol'])
                line.append(rsg_to_symbol[identifier]['gene_id'])
            except KeyError:
                logger.info("Can't identify gene symbol for %s", line[0])
                update_success = False
                try:
                    record = validator.entrez_efetch(db="nucleotide", id=line[0], rettype="gb", retmode="text")
                except Exception:
                    pass
                else:
                    gene_found = False
                    gene_id_found = False
                    for ft in record.features:
                        if "gene" in ft.qualifiers.keys():
                            gene_found = True
                            try:
                                # Line format is e.g. ['NG_030321.1', 'NC_000019.9', 'GRCh37', '7791843', '7802057', '-',
                                # 'CLEC4G', '339390'] We updated to the HGNC ID
                                line.append(ft.qualifiers["gene"][0])
                                line.append(ft.qualifiers["db_xref"][2].split(":")[-1])
                                update_success = True
                            except IndexError:
                                if "HGNC:" in str(ft.qualifiers["db_xref"]):
                                    for xref in ft.qualifiers["db_xref"]:
                                        if "HGNC" in xref:
                                            gene_id_found = True
                                            line.append(xref.split(":")[-1])
                                            update_success = True
                                else:
                                    # Will capture control elements like promoters, but captures the Entrez ID not the
                                    # HGNC ID
                                    for xref in ft.qualifiers["db_xref"]:
                                        if "GeneID" in xref:
                                            gene_id_found = True
                                            line.append(xref.split(":")[-1])
                                            update_success = True

                            except KeyError:
                                logger.info("Can't identify an HGNC ID for gene symbol for %s", line[0])
                            except Exception:
                                pass

                    if gene_found is False:
                        logger.info("Can't create an update gene symbol for %s", line[0])
                        for ft in record.features:
                            try:
                                if ft.type == "regulatory" and "GeneID:" in str(ft.qualifiers["db_xref"]):
                                    gene_id_found = True
                                    match = re.search(r'\(LOC\d+\)', record.description)
                                    if match:
                                        loc_var = re.sub(r'[()]', '', match.group())
                                        line.append(loc_var)
                                    else:
                                        continue
                                    line.append(ft.qualifiers["db_xref"][0].split(":")[-1])
                                    if "GRCh37" in ft.qualifiers["note"][0]:
                                        line[2] = "GRCh37"
                                        chr_num = line[1].split("NC_")[1].split(".")
                                        chr_num = int(chr_num[0])
                                        chr_acc = VariantValidator.modules.seq_data.to_accession(str(chr_num), "GRCh37")
                                        line[1] = chr_acc

                                    update_success = True
                                    break
                            except Exception:
                                continue
                    if gene_id_found is False:
                        logger.info("Can't create an update gene ID for %s", line[0])

        if update_success is False:
            print("Failed")
            print("Gene info found = ", gene_found)
            print("Gene ID found = ", gene_id_found)
            print(line)
            missing.append(line[0])
        else:
            print("Success")
            print(line)

    # Create a list of data to write to the database
    print("MAKE DATABASE ENTRIES")
    to_mysql = []
    for line in db:
        if line[0] in missing:
            continue
        # Only gap-less RefSeqGenes will have passed. The rest will be alternatively curated
        # Take the mapping data
        write = copy.deepcopy(line[0:6])
        # add RSG ranges
        end_rsg = int(line[4]) - int(line[3]) + 1
        end_rsg = str(end_rsg)
        write.append(end_rsg)
        # Create block data chr then rsg
        chr_block = str(line[3]) + '-' + str(line[4])
        write.append(chr_block)
        rsg_block = '1-' + str(write[6])
        write.append(rsg_block)
        # Add gene ID and Gene symbol(s)
        write.append(line[7])
        write.append(line[6])

        to_mysql.append(write)

    # Set up code to write to database
    print("WRITE TO DATABASE")
    for line in to_mysql:
        print("\nWriting to database:")
        print(line)
        current_symbol = dbcnx.get_gene_symbol_from_refseq_id(line[0])
        if line[10] != current_symbol:
            if current_symbol != 'none':
                line[10] = current_symbol
        try:
            dbcnx.update_refseqgene_loci(line)
        except Exception as e:
            # Badly formatted data due to RefSeqGene record errors so fail on our QC
            print("Failed to write to database: ", str(e))
            # import traceback
            # traceback.print_exc()
            missing.append(line[0])
            continue

    logger.info('Total NG_ to NC_ alignments = ' + str(total_rsg_to_nc))
    logger.info('Gaps within NG_ to NC_ alignments = ' + str(total_rsg_to_nc_rejected))
    print('Total NG_ to NC_ alignments = ' + str(total_rsg_to_nc))
    print('Gaps within NG_ to NC_ alignments = ' + str(total_rsg_to_nc_rejected))
    print('Total of badly formatted records = ' + str(len(missing)))
    return


def update_lrg(dbcnx):
    logger.debug('Updating LRG lookup tables')

    lr2rs_download = requests.get('http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_transcripts_xrefs.txt')
    lr2rs = lr2rs_download.text.strip().split('\n')

    # Download
    lrg_status_download = requests.get('http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_GRCh38.txt')
    lrg_status = lrg_status_download.text.strip().split('\n')

    # Download LRG transcript (_t) to LRG Protein (__p) data file
    lr_t2p_downloaded = requests.get('http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_proteins_RefSeq.txt')
    lr_t2p = lr_t2p_downloaded.text.strip().split('\n')

    # Dictionary the status by LRG_ID
    lrg_status_dict = {}
    # Compile dictionary
    for line in lrg_status:
        if line.startswith('#'):
            continue
        data = line.split()
        if len(data) < 3:
            continue
        lrg_status_dict[data[0]] = data[2]

    # Required lookup tables
    # LRG_ID	GeneSymbol	RefSeqGeneID	status
    # LRG_ID	RefSeqTranscriptID
    # LRG_T2LRG_P

    logger.debug('Update LRG and LRG_transcript lookup tables')
    # Populate lists lrg_rs_lookup (LRG to RefSeqGene) and lrg_t2nm_ (LRG Transcript to RefSeq Transcript)
    for line in lr2rs:
        if line.startswith('#'):
            continue
        data = line.split()
        if len(data) < 5:
            continue
        # Assign objects
        lrg_id = data[0]
        symbol = data[1]
        rsgid = data[2]
        lrg_tx = str(data[0]) + str(data[3])
        rstid = data[4]
        status = lrg_status_dict[lrg_id]
        # pass data to relevant lists
        # lrg_rs_lookup
        lrg_rs_lookup = [lrg_id, symbol, rsgid, status]

        # update LRG to RefSeqGene database
        dbcnx.update_lrg_rs_lookup(lrg_rs_lookup)

        # lrg_t2nm_
        lrgtx_to_rst_id = [lrg_tx, rstid]
        # update database
        dbcnx.update_lrgt_rst(lrgtx_to_rst_id)

    logger.debug('Update LRG protein lookup table')
    # Populate LRG protein RefSeqProtein lokup table
    for line in lr_t2p:
        if line.startswith('#'):
            continue
        data = line.split()
        # Assign objects
        lrg_p = data[0]
        rs_p = data[1]
        # update LRG to RefSeqGene database
        dbcnx.update_lrg_p_rs_p_lookup(lrg_p, rs_p)

    logger.info('LRG lookup tables updated')
    return


def count_ng_nc(line):
    # Count NG_ to NC_ and remove the entries we don't care about!
    if 'NC_' in line and 'NG_' in line:
        pass
    else:
        return None
    if '#' in line:
        return 'failed'

    if 'gap_count=0' not in line:
        return 'rejected'

    else:
        line = line.strip()
        info = line.split('\t')
        if len(info) != 9:
            return 'failed'

        metrics = info[8].split(';')
        id_ori = metrics[1].replace('Target=', '')
        id_ori_list = id_ori.split()
        entry = {'rsg_id': id_ori_list[0], 'chr_id': info[0], 'rsg_start': info[3], 'rsg_end': info[4],
                 'ori': id_ori_list[3]}
        return entry


def map_line(line, genome, rsg_id_info):
    ml = []
    link = line['rsg_id']
    ml.append(link)
    ml.append(line['chr_id'])
    ml.append(genome)
    ml.append(line['rsg_start'])
    ml.append(line['rsg_end'])
    ml.append(line['ori'])
    # Add the additional data from rsg_id_info
    for data in rsg_id_info:
        if link == data['rsg_id']:
            ml.append(data['symbol'])
            ml.append(data['gene_id'])
    # Create the entry and append to db
    return ml

# <LICENSE>
# Copyright (C) 2016-2024 VariantValidator Contributors
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
