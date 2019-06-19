from .logger import Logger
from . import vvFunctions as fn
from .vvFunctions import handleCursor
#from vvDBInsert import vvDBInsert
#from vvDBGet import vvDBGet
from . import vvDBInsert
import urllib.request, urllib.error, urllib.parse
import copy

import re
import os

class vvDatabase(vvDBInsert.Mixin):
    '''
    This class contains and handles the mysql connections for the variant validator database.

    It now uses mixins, and the order of inheritance is
    vvDBInit.Mixin
       v
    vvDBGet.Mixin
       v
    vvDBInsert.Mixin
       v
    vvDatabase
    '''
    # from dbquery
    @handleCursor
    def query_with_fetchone(self,entry, table):
        #if table == 'transcript_info':
        query = "SELECT refSeqID, description, transcriptVariant, currentVersion, hgncSymbol, utaSymbol, updated, IF(updated < NOW() - INTERVAL 3 MONTH , 'true', 'false') FROM transcript_info WHERE refSeqID = '%s'" %(entry)
        self.cursor.execute(query)
        row = self.cursor.fetchone()
        if row is None:
            row = ['none', 'No data']
            Logger.debug("No data returned from query " + str(query))
        return row
    # From data
    def data_add(self, accession, validator):
        '''
        # Add accurate transcript descriptions to the database
        :param accession:
        :return:
        '''
        self.update_transcript_info_record(accession, validator)
        entry = self.in_entries(accession, 'transcript_info')
        return entry

    def in_entries(self,entry, table):
        '''
        Retrieve transcript information
        :param entry:
        :param table:
        :return:
        '''
        data={}
        if table == 'transcript_info':
            row = self.query_with_fetchone(entry, table)
            if row[0] == 'error':
                data['error'] = row[0]
                data['description'] = row[1]
            elif row[0] == 'none':
                data['none'] = row[0]
                data['description'] = row[1]
            else:
                data['accession'] = row[0]
                data['description'] = row[1]
                data['variant'] = row[2]
                data['version'] = row[3]
                data['hgnc_symbol'] = row[4]
                data['uta_symbol'] = row[5]
                data['updated'] = row[6]
                data['expiry'] = row[7]
        return data
    def update_transcript_info_record(self,accession, validator):
        '''
        # Search Entrez for corresponding record for the RefSeq ID
        '''
        # Prime these entries, just in case.
        previous_entry = self.in_entries(accession, 'transcript_info')
        accession = accession
        if 'none' not in previous_entry.keys():
            description = previous_entry['description']
            variant = previous_entry['variant']
            version = previous_entry['version']
            hgnc_symbol = previous_entry['hgnc_symbol']
            uta_symbol = previous_entry['uta_symbol']
        try:
            record = validator.entrez_efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            version = record.id
            description = record.description
            variant = '0'

            if 'transcript variant' in description:
                tv = re.search('transcript variant \w+', description)
                tv = str(tv.group(0))
                tv = tv.replace('transcript variant', '')
                variant = tv.strip()
                variant = variant.upper()  # Some tv descriptions are a or A
            else:
                variant = '0'

            # Get information from UTA
            try:
                uta_info = validator.hdp.get_tx_identity_info(version)
            except:
                version_ac_ver = version.split('.')
                version = version_ac_ver[0] + '.' + str(int(version_ac_ver[1]) - 1)
                uta_info = validator.hdp.get_tx_identity_info(version)

            uta_symbol = str(uta_info[6])

            # First perform a search against the input gene symbol or the symbol inferred from UTA
            initial = fn.hgnc_rest(path = "/fetch/symbol/" + uta_symbol)
            # Check for a record
            if str(initial['record']['response']['numFound']) != '0':
                hgnc_symbol = uta_symbol
            # No record found, is it a previous symbol?
            else:
                # Search hgnc rest to see if symbol is out of date
                rest_data = fn.hgnc_rest(path = "/search/prev_symbol/" + uta_symbol)
                # If the name is correct no record will be found
                if rest_data['error'] == 'false':
                    if int(rest_data['record']['response']['numFound']) == 0:
                        hgnc_symbol = uta_info[6]
                    else:
                        hgnc_symbol = rest_data['record']['response']['docs'][0]['symbol']
                else:
                    hgnc_symbol = 'unassigned'

        # List of connection error types. May need to be expanded.
        # Outcome - Put off update for 3 months!
        except Exception as e:
            if not str(e) == '<urlopen error [Errno -2] Name or service not known>':
                # Issues with DNSSEC for the nih.gov
                raise

        # Query information
        # query_info = [accession, description, variant, version, hgnc_symbol, uta_symbol]
        query_info = [version, description, variant, version, hgnc_symbol, uta_symbol]
        table='transcript_info'

        # Update the transcript_info table (needs plugging in)
        returned_data = self.in_entries(version, table)
        # If the entry is not in the database add it
        if 'none' in returned_data:
            self.add_entry(version, query_info, table)
        # If the data in the entry has changed, update it
        else:
            self.update_entry(version, query_info, table)
        return

    def update_refSeqGene_loci(self,rsg_data):
        # First query the database
        entry_exists = self.get_refSeqGene_data_by_refSeqGeneID(rsg_data[0], rsg_data[2])
        if entry_exists[0] == 'none':
            self.insert_refSeqGene_data(rsg_data)
        else:
            self.update_refSeqGene_data(rsg_data)
    def update_lrg_rs_lookup(self,lrg_rs_lookup):
        # First query the database
        rsgID = self.get_RefSeqGeneID_from_lrgID(lrg_rs_lookup[0])
        if rsgID == 'none':
            self.insert_RefSeqGeneID_from_lrgID(lrg_rs_lookup)
    def update_lrgt_rst(self,lrgtx_to_rstID):
        # First query the database
        rstID = self.get_RefSeqTranscriptID_from_lrgTranscriptID(lrgtx_to_rstID[0])
        if rstID == 'none':
            self.insert_LRG_transcript_data(lrgtx_to_rstID)
    def update_lrg_p_rs_p_lookup(self,lrg_p, rs_p):
        # First query the database
        rspID = self.get_RefSeqProteinID_from_lrgProteinID(lrg_p)
        if rspID == 'none':
            self.insert_LRG_protein_data(lrg_p, rs_p)
    # From variantValidator.py
    def update_vv_data(self):
        # Update refSeqGene Primary assembly alignment data
        self.update_rsg()
        # Update LRG records
        self.update_lrg()
    # From update_refseqgene_nomissmatch.py
    def update_rsg(self):
        Logger.info('Updating RefSeqGene no Missmatch MySQL data')
        # Set os path
        # Set up os paths data and log folders
        ROOT = os.path.dirname(os.path.abspath(__file__))

        # Download data from RefSeqGene
        # Download data
        rsg = urllib.request.Request('http://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/gene_RefSeqGene')
        response = urllib.request.urlopen(rsg)
        rsg_file = response.read()
        rsg_data_line = rsg_file.split('\n')
        rsg_data = []
        for data in rsg_data_line:
            rsg_data.append(data)

        # Download data
        grch37 = urllib.request.Request(
            'http://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/GCF_000001405.25_refseqgene_alignments.gff3')
        response = urllib.request.urlopen(grch37)
        grch37_file = response.read()
        grch37_data_line = grch37_file.split('\n')
        grch37_align_data = []
        for data in grch37_data_line:
            grch37_align_data.append(data)

        # Download data
        grch38 = urllib.request.Request(
            'http://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/GCF_000001405.28_refseqgene_alignments.gff3')
        response = urllib.request.urlopen(grch38)
        grch38_file = response.read()
        grch38_data_line = grch38_file.split('\n')
        grch38_align_data = []
        for data in grch38_data_line:
            grch38_align_data.append(data)

        # Open Lists
        # rsg_data = open(os.path.join(ROOT, 'gene_RefSeqGene'), 'r')
        rsg_id_info = []
        # grch37_align_data = open(os.path.join(ROOT, 'GCF_000001405.25_refseqgene_alignments.gff3'), 'r')
        grch37_align = []
        # grch38_align_data = open(os.path.join(ROOT, 'GCF_000001405.28_refseqgene_alignments.gff3'), 'r')
        grch38_align = []

        # Place the required data from each file into a dictionary
        hash = re.compile('#')
        for line in rsg_data:
            if hash.search(line):
                pass
            else:
                line = line.strip()
                info = line.split()
                if len(info) == 0:
                    pass
                else:
                    dict = {'symbol': info[2], 'rsg_id': info[3], 'gene_id': info[1]}
                    rsg_id_info.append(dict)

        # Create dictionary to store RefSeqGene and gene symbol data NOTE RefSeqGene ID stored without version number!
        rsg_to_symbol = {}
        # Collect the data
        for ent in rsg_id_info:
            rsg_id = copy.deepcopy(ent['rsg_id'])
            rsg_id = rsg_id.split('.')[0]
            rsg_to_symbol[rsg_id] = {'symbol': ent['symbol'], 'gene_id': ent['gene_id']}

        # Count total number of NG to NC mappings
        total_rsg_to_nc = 0
        total_rsg_to_nc_rejected = 0
        for line in grch37_align_data:
            # Count NG_ to NC_ and remove the entries we don't care about!
            if re.search('NC_', line) and re.search('NG_', line):
                total_rsg_to_nc = total_rsg_to_nc + 1
            else:
                continue
            if hash.search(line):
                pass
            elif not re.search('gap_count=0', line):
                if re.search('NC_', line) and re.search('NG_', line):
                    total_rsg_to_nc_rejected = total_rsg_to_nc_rejected + 1
                # print line
                pass
            else:
                line = line.strip()
                info = line.split('\t')
                if len(info) != 9:
                    pass
                else:
                    metrics = info[8].split(';')
                    id_ori = metrics[1].replace('Target=', '')
                    id_ori_list = id_ori.split()
                    dict = {'rsg_id': id_ori_list[0], 'chr_id': info[0], 'rsg_start': info[3], 'rsg_end': info[4],
                            'ori': id_ori_list[3]}
                    grch37_align.append(dict)

        for line in grch38_align_data:
            if re.search('NC_', line) and re.search('NG_', line):
                total_rsg_to_nc = total_rsg_to_nc + 1
            else:
                continue
            if hash.search(line):
                pass
            elif not re.search('gap_count=0', line):
                if re.search('NC_', line) and re.search('NG_', line):
                    total_rsg_to_nc_rejected = total_rsg_to_nc_rejected + 1
                # print line
                pass
            else:
                line = line.strip()
                info = line.split('\t')
                if len(info) != 9:
                    pass
                else:
                    metrics = info[8].split(';')
                    id_ori = metrics[1].replace('Target=', '')
                    id_ori_list = id_ori.split()
                    dict = {'rsg_id': id_ori_list[0], 'chr_id': info[0], 'rsg_start': info[3], 'rsg_end': info[4],
                            'ori': id_ori_list[3]}
                    grch38_align.append(dict)

        # Create a data array containing the database
        db = []
        # map line
        for line in grch37_align:
            ml = []
            link = line['rsg_id']
            ml.append(link)
            ml.append(line['chr_id'])
            ml.append('GRCh37')
            ml.append(line['rsg_start'])
            ml.append(line['rsg_end'])
            ml.append(line['ori'])
            # Add the additional data from rsg_id_info
            for data in rsg_id_info:
                if link == data['rsg_id']:
                    ml.append(data['symbol'])
                    ml.append(data['gene_id'])
                else:
                    continue
            # Create the entry and append to db
            db.append(ml)

        for line in grch38_align:
            ml = []
            link = line['rsg_id']
            ml.append(link)
            ml.append(line['chr_id'])
            ml.append('GRCh38')
            ml.append(line['rsg_start'])
            ml.append(line['rsg_end'])
            ml.append(line['ori'])
            # Add the additional data from rsg_id_info
            for data in rsg_id_info:
                if link == data['rsg_id']:
                    ml.append(data['symbol'])
                    ml.append(data['gene_id'])
                else:
                    continue
            # Create the entry and append to db
            db.append(ml)

        # Known missing identifiers
        known = {
                'NG_021289.1' : {'symbol' : 'CFAP47', 'gene_id' : '286464'},
                'NG_027707.1' : {'symbol' : 'DUX4L1', 'gene_id' : '22947'},
                'NG_033266.1' : {'symbol' : 'DSE', 'gene_id': '29940'},
                'NG_061543.1' : {'symbol' : 'CYP1A2', 'gene_id': '1544'},
                'NG_061374.1' : {'symbol' : 'CYP1A1', 'gene_id': '1543'},
                'NG_059281.1' : {'symbol' : 'HBB', 'gene_id': '3043'},
                'NG_012639.1' : {'symbol' : 'VHLL', 'gene_id': '391104'},
                'NG_059186.1' : {'symbol' : 'HBA1', 'gene_id': '3040'},
                'NG_059271.1' : {'symbol' : 'HBA2', 'gene_id': '3040'}
                }

        # Known Obsolete identifiers
        obsolete = {
            'NG_016553.1': 'OBSOLETE',
            'NG_012639.1': 'Removed due to questionable status'
        }

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
                    try:
                        line.append(known[line[0]]['symbol'])
                        line.append(known[line[0]]['gene_id'])
                    except KeyError:
                        check = obsolete[line[0]]
                        Logger.info(str(line[0]) + ' : ' + check)

        # Open a text file to be used as a simple database and write the database
        # rsg_db = open(os.path.join(ROOT, 'rsg_chr_db.txt'), 'w')

        to_mysql = []
        for line in db:
            if line[0] in list(obsolete.keys()):
                continue
            # Only gap-less RefSeqGenes will have passed. The rest will be alternatively curated
            write = []
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
            # write_me = '\t'.join(write)
            # rsg_db.write(write_me + '\n')
            del write[6]
            to_mysql.append(write)

        # Set up code to write to database
        for line in to_mysql:
            current_symbol = self.get_gene_symbol_from_refSeqGeneID(line[0])
            if line[10] == current_symbol:
                pass
            else:
                if current_symbol != 'none':
                    line[10] = current_symbol
                else:
                    pass
            self.update_refSeqGene_loci(line)

        # Close database
        # rsg_db.close()

        Logger.info('Total NG_ to NC_ alignments = ' + str(total_rsg_to_nc))
        Logger.info('Gapps within NG_ to NC_ alignments = ' + str(total_rsg_to_nc_rejected))

        Logger.info('complete')
        return
    #from compile_lrg_data, this function was originally just called "update"
    def update_lrg(self):
        Logger.info('Updating LRG lookup tables')
        lr2rs_download = urllib.request.Request('http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_transcripts_xrefs.txt')
        # Open and read
        lr2rs_data = urllib.request.urlopen(lr2rs_download)
        lr2rs = lr2rs_data.read()
        # List the data
        lr2rs = lr2rs.strip()
        lr2rs = lr2rs.split('\n')

        # Download
        lrg_status_download = urllib.request.Request('http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_GRCh38.txt')
        # Open and read
        lrg_status_data = urllib.request.urlopen(lrg_status_download)
        lrg_status = lrg_status_data.read()
        # List the data
        lrg_status = lrg_status.strip()
        lrg_status = lrg_status.split('\n')

        # Download
        rs2lr_download = urllib.request.Request('http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene')
        # Open and read
        rs2lr_data = urllib.request.urlopen(rs2lr_download)
        rs2lr = rs2lr_data.read()
        # List the data
        rs2lr = rs2lr.strip()
        rs2lr = rs2lr.split('\n')

        # Download LRG transcript (_t) to LRG Protein (__p) data file
        lr_t2p_downloaded = urllib.request.Request('http://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_proteins_RefSeq.txt')
        # Open and read
        lr_t2p_data = urllib.request.urlopen(lr_t2p_downloaded)
        lr_t2p = lr_t2p_data.read()
        # List the data
        lr_t2p = lr_t2p.strip()
        lr_t2p = lr_t2p.split('\n')

        # Dictionary the status by LRG_ID
        lrg_status_dict = {}
        # Compile dictionary
        for line in lrg_status:
            if re.search('^#', line):
                continue
            else:
                list = line.split()
                lrgID = list[0]
                stat = list[2]
                lrg_status_dict[lrgID] = stat

        # Required lookup tables
        # LRG_ID	GeneSymbol	RefSeqGeneID	status
        # LRG_ID	RefSeqTranscriptID
        # LRG_T2LRG_P

        Logger.info('Update LRG and LRG_transcript lookup tables')
        # Populate lists lrg_rs_lookup (LRG to RefSeqGene) and lrg_t2nm_ (LRG Transcript to RefSeq Transcript)
        for line in lr2rs:
            if re.search('^#', line):
                continue
            else:
                list = line.split()
                # Assign objects
                lrg_id = list[0]
                symbol = list[1]
                rsgid = list[2]
                lrg_tx = str(list[0]) + str(list[3])
                rstid = list[4]
                status = lrg_status_dict[lrg_id]
                # pass data to relevant lists
                # lrg_rs_lookup
                lrg_rs_lookup = [lrg_id, symbol, rsgid, status]

                # update LRG to RefSeqGene database
                self.update_lrg_rs_lookup(lrg_rs_lookup)

                # lrg_t2nm_
                lrgtx_to_rstID = [lrg_tx, rstid]
                # update database
                self.update_lrgt_rst(lrgtx_to_rstID)

        Logger.info('Update LRG protein lookup table')
        # Populate LRG protein RefSeqProtein lokup table
        for line in lr_t2p:
            if re.search('^#', line):
                continue
            else:
                list = line.split()
                # Assign objects
                lrg_p = list[0]
                rs_p = list[1]
                # update LRG to RefSeqGene database
                self.update_lrg_p_rs_p_lookup(lrg_p, rs_p)

        Logger.info('LRG lookup tables updated')
        return
    #From ref_seq_type
    def ref_type_assign(self,accession):
        if 'NC_' in accession or 'NG_' in accession or 'NT_' in accession or 'NW_' in accession:
            ref_type = ':g.'
        elif re.match('NM_', accession):
            ref_type = ':c.'
        elif re.match('NR_', accession):
            ref_type = ':n.'
        elif re.match('NP_', accession):
            ref_type = ':p.'
        elif re.match('LRG_', accession):
            if re.search('t', accession):
                refseqtranscript_reference = self.get_RefSeqTranscriptID_from_lrgTranscriptID(accession)
                if re.match('NM_', refseqtranscript_reference):
                    ref_type = ':c.'
                else:
                    ref_type = ':n.'
            elif re.search('_p', accession):
                ref_type = ':p.'
            else:
                ref_type = ':g.'
        return ref_type


