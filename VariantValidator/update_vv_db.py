# -*- coding: utf-8 -*-

import os
import urllib.request, urllib.error, urllib.parse
import copy
from configparser import ConfigParser
from .modules import vvDatabase
from . import variantValidator
from . import configure


def update():

    config = ConfigParser()
    config.read(configure.CONFIG_DIR)

    dbConfig = {
        'user': config["mysql"]["user"],
        'password': config["mysql"]["password"],
        'host': config["mysql"]["host"],
        'database': config["mysql"]["database"],
        'raise_on_warnings': True
    }
    # Create database access objects
    db = vvDatabase.vvDatabase(variantValidator.Validator(), dbConfig)

    update_refseq(db)
    update_lrg(db)


def update_refseq(dbcnx):
    print('Updating RefSeqGene no Missmatch MySQL data')
    # Set os path
    # Set up os paths data and log folders
    ROOT = os.path.dirname(os.path.abspath(__file__))

    # Download data from RefSeqGene
    # Download data
    rsg = urllib.request.Request('http://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/gene_RefSeqGene')
    response = urllib.request.urlopen(rsg)
    rsg_file = response.read()
    rsg_data_line = rsg_file.split(b'\n')
    rsg_data = []
    for data in rsg_data_line:
        rsg_data.append(data.decode())

    # Download data
    grch37 = urllib.request.Request(
        'http://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/GCF_000001405.25_refseqgene_alignments.gff3')
    response = urllib.request.urlopen(grch37)
    grch37_file = response.read()
    grch37_data_line = grch37_file.split(b'\n')
    grch37_align_data = []
    for data in grch37_data_line:
        grch37_align_data.append(data.decode())

    # Download data
    grch38 = urllib.request.Request(
        'http://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/GCF_000001405.28_refseqgene_alignments.gff3')
    response = urllib.request.urlopen(grch38)
    grch38_file = response.read()
    grch38_data_line = grch38_file.split(b'\n')
    grch38_align_data = []
    for data in grch38_data_line:
        grch38_align_data.append(data.decode())

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
    for line in db:
        try:
            line[6]
        except IndexError:
            try:
                identifier = copy.deepcopy(line[0])
                identifier = identifier.split('.')[0]
                line.append(rsg_to_symbol[identifier]['symbol'])
                line.append(rsg_to_symbol[identifier]['gene_id'])
            except KeyError:
                print("Can't identify gene symbol for %s" % line[0])
                missing.append(line[0])

    # Open a text file to be used as a simple database and write the database
    # rsg_db = open(os.path.join(ROOT, 'rsg_chr_db.txt'), 'w')

    to_mysql = []
    for line in db:
        if line[0] in missing:
            continue
        # Only gap-less RefSeqGenes will have passed. The rest will be alternatively curated
        # Take the mapping data
        write = copy.deepcopy(line[0:6])
        # add RSG ranges
        write.append('1')
        end_rsg = int(line[4]) - int(line[3]) + 1
        end_rsg = str(end_rsg)
        write.append(end_rsg)
        # Create block data chr then rsg
        chr_block = str(line[3]) + '-' + str(line[4])
        write.append(chr_block)
        rsg_block = str(write[6]) + '-' + str(write[7])
        write.append(rsg_block)
        # Add gene ID and Gene symbol(s)
        write.append(line[7])
        write.append(line[6])

        del write[6]
        to_mysql.append(write)

    # Set up code to write to database
    for line in to_mysql:
        current_symbol = dbcnx.get_gene_symbol_from_refSeqGeneID(line[0])
        if line[10] != current_symbol:
            if current_symbol != 'none':
                line[10] = current_symbol
        dbcnx.update_refSeqGene_loci(line)

    print('Total NG_ to NC_ alignments = ' + str(total_rsg_to_nc))
    print('Gaps within NG_ to NC_ alignments = ' + str(total_rsg_to_nc_rejected))

    return

def update_lrg(dbcnx):
    print('Updating LRG lookup tables')

    lr2rs_download = urllib.request.Request('http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_transcripts_xrefs.txt')
    # Open and read
    lr2rs_data = urllib.request.urlopen(lr2rs_download)
    lr2rs = lr2rs_data.read()
    # List the data
    lr2rs = lr2rs.strip().decode()
    lr2rs = lr2rs.split('\n')

    # Download
    lrg_status_download = urllib.request.Request('http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_GRCh38.txt')
    # Open and read
    lrg_status_data = urllib.request.urlopen(lrg_status_download)
    lrg_status = lrg_status_data.read()
    # List the data
    lrg_status = lrg_status.strip().decode()
    lrg_status = lrg_status.split('\n')

    # Download LRG transcript (_t) to LRG Protein (__p) data file
    lr_t2p_downloaded = urllib.request.Request('http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_proteins_RefSeq.txt')
    # Open and read
    lr_t2p_data = urllib.request.urlopen(lr_t2p_downloaded)
    lr_t2p = lr_t2p_data.read()
    # List the data
    lr_t2p = lr_t2p.strip().decode()
    lr_t2p = lr_t2p.split('\n')

    # Dictionary the status by LRG_ID
    lrg_status_dict = {}
    # Compile dictionary
    for line in lrg_status:
        if line.startswith('#'):
            continue
        data = line.split()
        lrg_status_dict[data[0]] = data[2]

    # Required lookup tables
    # LRG_ID	GeneSymbol	RefSeqGeneID	status
    # LRG_ID	RefSeqTranscriptID
    # LRG_T2LRG_P

    print('Update LRG and LRG_transcript lookup tables')
    # Populate lists lrg_rs_lookup (LRG to RefSeqGene) and lrg_t2nm_ (LRG Transcript to RefSeq Transcript)
    for line in lr2rs:
        if line.startswith('#'):
            continue
        data = line.split()
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
        lrgtx_to_rstID = [lrg_tx, rstid]
        # update database
        dbcnx.update_lrgt_rst(lrgtx_to_rstID)

    print('Update LRG protein lookup table')
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

    print('LRG lookup tables updated')
    return


def count_ng_nc(line):
    # Count NG_ to NC_ and remove the entries we don't care about!
    if 'NC_' in line and 'NG_' in line:
        # print(line)
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
