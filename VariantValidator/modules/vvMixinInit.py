import os
from configparser import ConfigParser,RawConfigParser
import hgvs
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.dataproviders.seqfetcher
import hgvs.assemblymapper
import hgvs.variantmapper
import hgvs.sequencevariant
import hgvs.validator
import hgvs.exceptions
import hgvs.location
import hgvs.posedit
import hgvs.edit
import hgvs.normalizer
from Bio.Seq import Seq

import re
import copy
#import io
from vvDatabase import vvDatabase
from vvLogging import logger
import vvFunctions as fn


'''
This file contains the validator object, which is instantiated in order to perform validator functions.
The validator contains configuration information and permanent copies of database links and the like.
Much of the validator's inner workings are stored in special one-off function container objects:
validator.db : The validator's MySQL database access functions

The validator configuration is stored in ~/.config/VariantValidator/config.ini . This is loaded
when the validator object is initialized.

Running variant validator should hopefully be as simple as writing a script like this:
import VariantValidator

val=Validator()
val.validate("some kind of gene situation","the transcripts to use")

'''

'''        
        Renaming of variables :
        'seqrepo_directory': HGVS_SEQREPO_DIR,           #self.seqrepoPath
        'uta_url': UTA_DB_URL,                           #self.utaPath
        'py_liftover_directory': PYLIFTOVER_DIR,         #self.liftoverPath
        'variantvalidator_data_url': VALIDATOR_DB_URL,   #self.db.path
        'entrez_id': ENTREZ_ID,                          #self.entrezID
        'variantvalidator_version': VERSION,             #self.version
        'variantvalidator_hgvs_version': hgvs_version,   #self.hgvsVersion
        'uta_schema': str(hdp.data_version()),           #self.uta_schema
        'seqrepo_db': HGVS_SEQREPO_DIR.split('/')[-1]    #self.seqrepoVersion
'''



class Mixin():
    # This object contains configuration options for the validator, but it inherits the mixin
    # class in vvCore that contains the enormous validator function.
    def __init__(self):
        # First load from the configuration file, if it exists.
        configName="config.ini"
        homePath=os.path.expanduser("~")
        configPath=os.path.join(homePath,".config","VariantValidator")
        if not os.path.isdir(configPath):
            os.makedirs(configPath)
        # Now configpath points to the config file itself.
        configPath=os.path.join(configPath,configName)
        # Does the file exist?
        if not os.path.exists(configPath):
            self.createConfig(configPath)

        # Load the configuration file.
        with open(configPath) as file:
            lines=file.read()
        config=RawConfigParser(allow_no_value=True)
        config.read(configPath)
        # The custom vvLogging module will set itself up using the VALDIATOR_DEBUG environment variable.
        logString = config["logging"]['string']
        os.environ["VALIDATOR_DEBUG"] = logString

        # Handle databases
        self.entrezID=config["EntrezID"]["entrezID"]
        if config["seqrepo"]["location"]!=None:
            self.seqrepoVersion=config["seqrepo"]["version"]
            self.seqrepoPath=config["seqrepo"]["location"]+self.seqrepoVersion
            os.environ['HGVS_SEQREPO_DIR']=self.seqrepoPath
        else:
            raise ValueError("The seqrepo location has not been set in ~/.config/VariantValidator/config.ini")
        os.environ['UTA_DB_URL']=config["uta"]["location"]+config["uta"]["version"]
        self.utaPath=config["uta"]["location"]+config["uta"]["version"]
        self.dbConfig = {
            'user':    config["mysql"]["user"],
            'password':config["mysql"]["password"],
            'host':    config["mysql"]["host"],
            'database':config["mysql"]["database"],
    	    'raise_on_warnings': True
        }
        self.db=vvDatabase(self,self.dbConfig)
        # Set up versions
        __version__ = config["variantValidator"]['version']
        self.version=__version__
        if re.match('^\d+\.\d+\.\d+$', __version__) is not None:
            self.releasedVersion=True
            _is_released_version = True
        else:
            self.releasedVersion=False
        self.hgvsVersion=hgvs.__version__

        # Set up other configuration variables
        self.liftoverPath=config["liftover"]["location"]
        if not self.liftoverPath==None:
            os.environ['PYLIFTOVER_DIR']=self.liftoverPath
        self.entrezID=config["EntrezID"]['entrezid']

        # Set up HGVS
        # Configure hgvs package global settings
        hgvs.global_config.uta.pool_max = 25
        hgvs.global_config.formatting.max_ref_length = 1000000
        # Create HGVS objects
        self.hdp = hgvs.dataproviders.uta.connect(pooling=True)
        self.hp = hgvs.parser.Parser() #Parser
        self.vr = hgvs.validator.Validator(self.hdp) # Validator
        self.vm = hgvs.variantmapper.VariantMapper(self.hdp) # Variant mapper
        # Create a lose vm instance
        self.lose_vm = hgvs.variantmapper.VariantMapper(self.hdp,
                                                   replace_reference=True,
                                                   prevalidation_level=None
                                                   )
        self.nr_vm = hgvs.variantmapper.VariantMapper(self.hdp, replace_reference=False) #No reverse variant mapper
        self.sf = hgvs.dataproviders.seqfetcher.SeqFetcher() # Seqfetcher
        # Set standard genome builds
        self.genome_builds = ['GRCh37', 'hg19', 'GRCh38']
        self.utaSchema = str(self.hdp.data_version())

        # Create normalizer
        self.reverse_hn = hgvs.normalizer.Normalizer(self.hdp,
                                                cross_boundaries=False,
                                                shuffle_direction=5,
                                                alt_aln_method='splign'
                                                )

        # Create normalizer
        self.merge_normalizer = hgvs.normalizer.Normalizer(self.hdp,
                                                      cross_boundaries=False,
                                                      shuffle_direction=hgvs.global_config.normalizer.shuffle_direction,
                                                      alt_aln_method='splign',
                                                      validate=False
                                                      )
        self.reverse_merge_normalizer = hgvs.normalizer.Normalizer(self.hdp,
                                                              cross_boundaries=False,
                                                              shuffle_direction=hgvs.global_config.normalizer.shuffle_direction,
                                                              alt_aln_method='splign',
                                                              validate=False
                                                              )
        #create no_norm_evm
        self.no_norm_evm_38 = hgvs.assemblymapper.AssemblyMapper(self.hdp,
                                                            assembly_name='GRCh38',
                                                            alt_aln_method='splign',
                                                            normalize=False,
                                                            replace_reference=True
                                                            )

        self.no_norm_evm_37 = hgvs.assemblymapper.AssemblyMapper(self.hdp,
                                                            assembly_name='GRCh37',
                                                            alt_aln_method='splign',
                                                            normalize=False,
                                                            replace_reference=True
                                                            )



    def myConfig(self):
        #Returns configuration:
        #version, hgvs version, uta schema, seqrepo db.
        return {
            "Seqrepo path":self.seqrepoPath
        }
    def createConfig(self,outPath):
        # This function reads from the default configuration file stored in the same folder as this module,
        # and transfers it to outPath.
        # Outpath should include a filename.
        lines=[]
        inPath=os.path.join(os.path.dirname(os.path.realpath(__file__)),"defaultConfig.ini")
#        print(os.path.join(inPath,"defaultConfig.ini"))
        with open(inPath) as file:
            for l in file:
                lines.append(l)
        with open(outPath, "w") as file:
            for l in lines:
                file.write(l)
    def protein(self,variant, evm, hpUnused):
        # Set regular expressions for if statements
        pat_c = re.compile("\:c\.")  # Pattern looks for :c. Note (gene) has been removed

        # If the :c. pattern is present in the input variant
        if pat_c.search(variant):
            # convert the input string into a hgvs object
            var_c = self.hp.parse_hgvs_variant(variant)
            # Does the edit affect the start codon?
            if ((var_c.posedit.pos.start.base >= 1 and var_c.posedit.pos.start.base <= 3 and var_c.posedit.pos.start.offset == 0) or (
                    var_c.posedit.pos.end.base >= 1 and var_c.posedit.pos.end.base <= 3 and var_c.posedit.pos.end.offset == 0)) and not re.search('\*', str(
                    var_c.posedit.pos)):
                ass_prot = self.hdp.get_pro_ac_for_tx_ac(var_c.ac)
                if str(ass_prot) == 'None':
                    cod = str(var_c)
                    cod = cod.replace('inv', 'del')
                    cod = self.hp.parse_hgvs_variant(cod)
                    p = evm.c_to_p(cod)
                    ass_prot = p.ac
                var_p = hgvs.sequencevariant.SequenceVariant(ac=ass_prot, type='p', posedit='(Met1?)')
            else:
                var_p = evm.c_to_p(var_c)
            return var_p
        if re.search(':n.', variant):
            var_p = self.hp.parse_hgvs_variant(variant)
            var_p.ac = 'Non-coding transcript'
            var_p.posedit = ''
            return var_p
    def myc_to_p(self,hgvs_transcript, evm, re_to_p):
        # Create dictionary to store the information
        hgvs_transcript_to_hgvs_protein = {'error': '', 'hgvs_protein': '', 'ref_residues': ''}

        # Collect the associated protein
        if hgvs_transcript.type == 'c':
            associated_protein_accession = self.hdp.get_pro_ac_for_tx_ac(hgvs_transcript.ac)
            # This method sometimes fails
            if str(associated_protein_accession) == 'None':
                cod = str(hgvs_transcript)
                cod = cod.replace('inv', 'del')
                cod = self.hp.parse_hgvs_variant(cod)
                p = evm.c_to_p(cod)
                associated_protein_accession = p.ac
        else:
            pass

        # Check for non-coding transcripts
        if hgvs_transcript.type == 'c':
            # Handle non inversions with simple c_to_p mapping

            if (hgvs_transcript.posedit.edit.type != 'inv') and (hgvs_transcript.posedit.edit.type != 'delins') and (re_to_p is False):
                # Does the edit affect the start codon?
                if ((hgvs_transcript.posedit.pos.start.base >= 1 and hgvs_transcript.posedit.pos.start.base <= 3 and hgvs_transcript.posedit.pos.start.offset == 0) or (
                        hgvs_transcript.posedit.pos.end.base >= 1 and hgvs_transcript.posedit.pos.end.base <= 3 and hgvs_transcript.posedit.pos.end.offset == 0)) \
                        and not re.search('\*', str(
                    hgvs_transcript.posedit.pos)):
                    hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                        type='p', posedit='(Met1?)')
                else:
                    try:
                        hgvs_protein = evm.c_to_p(hgvs_transcript)
                    except IndexError as e:
                        error = str(e)
                        if re.search('string index out of range', error) and re.search('dup', str(hgvs_transcript)):
                            hgvs_ins = self.hp.parse_hgvs_variant(str(hgvs_transcript))
                            hgvs_ins = hn.normalize(hgvs_ins)
                            inst = hgvs_ins.ac + ':c.' + str(hgvs_ins.posedit.pos.start.base - 1) + '_' + str(hgvs_ins.posedit.pos.start.base) + 'ins' + hgvs_ins.posedit.edit.ref
                            hgvs_transcript = self.hp.parse_hgvs_variant(inst)
                            hgvs_protein = evm.c_to_p(hgvs_transcript)

                try:
                    hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                    return hgvs_transcript_to_hgvs_protein
                except UnboundLocalError:
                    hgvs_transcript_to_hgvs_protein = self.myc_to_p(hgvs_transcript, evm, re_to_p = True)
                    return hgvs_transcript_to_hgvs_protein

            else:
                # Additional code required to process inversions
                # Note, this code was developed for VariantValidator and is not native to the biocommons hgvs Python package
                # Convert positions to n. position
                hgvs_naughty = self.vm.c_to_n(hgvs_transcript)

                # Collect the deleted sequence using fetch_seq
                del_seq = self.sf.fetch_seq(str(hgvs_naughty.ac), start_i=hgvs_naughty.posedit.pos.start.base - 1, end_i=hgvs_naughty.posedit.pos.end.base)

                # Make the inverted sequence
                my_seq = Seq(del_seq)

                if hgvs_transcript.posedit.edit.type == 'inv':
                    inv_seq = my_seq.reverse_complement()
                else:
                    inv_seq = hgvs_transcript.posedit.edit.alt
                    if inv_seq is None:
                        inv_seq = ''

                # Look for p. delins or del
                not_delins = True
                if hgvs_transcript.posedit.edit.type != 'inv':
                    try:
                        shifts = evm.c_to_p(hgvs_transcript)
                        if re.search('del', shifts.posedit.edit.type):
                            not_delins = False
                    except Exception:
                        not_delins = False
                else:
                    not_delins = False

                # Use inv delins code?
                if not_delins == False:
                    # Collect the associated protein
                    associated_protein_accession = self.hdp.get_pro_ac_for_tx_ac(hgvs_transcript.ac)

                    # Intronic inversions are marked as uncertain i.e. p.?
                    if re.search('\d+\-', str(hgvs_transcript.posedit.pos)) or re.search('\d+\+', str(hgvs_transcript.posedit.pos)) or re.search('\*', str(hgvs_transcript.posedit.pos)) or re.search('[cn].\-', str(hgvs_transcript)):
                        if ((
                                    hgvs_transcript.posedit.pos.start.base >= 1 and hgvs_transcript.posedit.pos.start.base <= 3 and hgvs_transcript.posedit.pos.start.offset == 0)
                            or
                            (hgvs_transcript.posedit.pos.end.base >= 1 and hgvs_transcript.posedit.pos.end.base <= 3 and hgvs_transcript.posedit.pos.end.offset == 0)) \
                                and not re.search('\*', str(hgvs_transcript.posedit.pos)):
                            hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession, type='p',
                                                                                posedit='(Met1?)')
                        else:
                            # Make the variant
                            hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession, type='p', posedit='?')
                        hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                        return hgvs_transcript_to_hgvs_protein
                    else:
                        # Need to obtain the cds_start
                        inf = self.hdp.get_tx_identity_info(hgvs_transcript.ac)
                        cds_start = inf[3]

                        # Extract the reference coding sequence from SeqRepo
                        try:
                            ref_seq = self.sf.fetch_seq(str(hgvs_naughty.ac))
                        except Exception as e:
                            error = str(e)
                            hgvs_transcript_to_hgvs_protein['error'] = error
                            return hgvs_transcript_to_hgvs_protein

                        # Create the variant coding sequence
                        var_seq = fn.n_inversion(ref_seq, del_seq, inv_seq,
                                                    hgvs_naughty.posedit.pos.start.base,
                                                    hgvs_naughty.posedit.pos.end.base)
                        # Translate the reference and variant proteins
                        prot_ref_seq = fn.translate(ref_seq, cds_start)

                        try:
                            prot_var_seq = fn.translate(var_seq, cds_start)
                        except IndexError:
                            hgvs_transcript_to_hgvs_protein['error'] = 'Cannot identify an in-frame Termination codon in the variant mRNA sequence'
                            hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession, type='p',
                                                                                    posedit='?')
                            hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                            return hgvs_transcript_to_hgvs_protein

                        if prot_ref_seq == 'error':
                            error = 'Unable to generate protein variant description'
                            hgvs_transcript_to_hgvs_protein['error'] = error
                            return hgvs_transcript_to_hgvs_protein
                        elif prot_var_seq == 'error':
                            # Does the edit affect the start codon?
                            if ((
                                        hgvs_transcript.posedit.pos.start.base >= 1 and hgvs_transcript.posedit.pos.start.base <= 3 and hgvs_transcript.posedit.pos.start.offset == 0)
                                or
                                (hgvs_transcript.posedit.pos.end.base >= 1 and hgvs_transcript.posedit.pos.end.base <= 3 and hgvs_transcript.posedit.pos.end.offset == 0)) \
                                    and not re.search('\*', str(hgvs_transcript.posedit.pos)):
                                hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession, type='p',
                                                                                    posedit='(Met1?)')

                                hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                                return hgvs_transcript_to_hgvs_protein
                            else:
                                error = 'Unable to generate protein variant description'
                                hgvs_transcript_to_hgvs_protein['error'] = error
                                return hgvs_transcript_to_hgvs_protein
                        else:
                            # Gather the required information regarding variant interval and sequences
                            if hgvs_transcript.posedit.edit.type != 'delins':
                                pro_inv_info = fn.pro_inv_info(prot_ref_seq, prot_var_seq)
                            else:
                                pro_inv_info = fn.pro_delins_info(prot_ref_seq, prot_var_seq)

                            # Error has occurred
                            if pro_inv_info['error'] == 'true':
                                error = 'Translation error occurred, please contact admin'
                                hgvs_transcript_to_hgvs_protein['error'] = error
                                return hgvs_transcript_to_hgvs_protein

                            # The Nucleotide variant has not affected the protein sequence i.e. synonymous
                            elif pro_inv_info['variant'] != 'true':
                                # Make the variant
                                hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession, type='p',
                                                                                    posedit='=')
                                hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                                return hgvs_transcript_to_hgvs_protein

                            else:
                                # Early termination i.e. stop gained
                                # if pro_inv_info['terminate'] == 'true':
                                #     end = 'Ter' + str(pro_inv_info['ter_pos'])
                                #     pro_inv_info['prot_ins_seq'].replace('*', end)

                                # Complete variant description
                                # Recode the single letter del and ins sequences into three letter amino acid codes
                                del_thr = fn.one_to_three(pro_inv_info['prot_del_seq'])
                                ins_thr = fn.one_to_three(pro_inv_info['prot_ins_seq'])

                                # Write the HGVS position and edit
                                del_len = len(del_thr)
                                from_aa = del_thr[0:3]
                                to_aa = del_thr[del_len - 3:]

                                # Handle a range of amino acids
                                if pro_inv_info['edit_start'] != pro_inv_info['edit_end']:
                                    if len(ins_thr) > 0:
                                        if re.search('Ter', del_thr) and ins_thr[-3:] != 'Ter':
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + str(
                                                pro_inv_info['edit_end']) + 'delins' + ins_thr + '?)'
                                        else:
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + str(
                                                pro_inv_info['edit_end']) + 'delins' + ins_thr + ')'
                                    else:
                                        if re.search('Ter', del_thr) and ins_thr[-3:] != 'Ter':
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + str(
                                                pro_inv_info['edit_end']) + 'del?)'
                                        else:
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + str(
                                                pro_inv_info['edit_end']) + 'del)'
                                else:
                                    # Handle extended proteins i.e. stop_lost
                                    if del_thr == 'Ter' and (len(ins_thr) > len(del_thr)):
                                        # Nucleotide variant range aligns to the Termination codon
                                        if ins_thr[-3:] == 'Ter':
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + str(
                                                ins_thr[:3]) + 'ext' + str(ins_thr[-3:]) + str((len(ins_thr) / 3) - 1) + ')'
                                        # Nucleotide variant range spans the Termination codon
                                        else:
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + str(
                                                ins_thr[:3]) + 'ext?)'

                                    # Nucleotide variation has not affected the length of the protein thus substitution or del
                                    else:
                                        if len(ins_thr) == 3:
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + ins_thr + ')'
                                        elif len(ins_thr) == 0:
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + 'del)'
                                        else:
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + 'delins' + ins_thr + ')'

                                # Complete the variant
                                hgvs_protein = hgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession, type='p',
                                                                                    posedit=posedit)

                                hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein

                else:
                    hgvs_transcript_to_hgvs_protein['hgvs_protein'] = shifts

                # Return
                return hgvs_transcript_to_hgvs_protein


        # Handle non-coding transcript and non transcript descriptions
        elif hgvs_transcript.type == 'n':
            # non-coding transcripts
            hgvs_protein = copy.deepcopy(hgvs_transcript)
            hgvs_protein.ac = 'Non-coding '
            hgvs_protein.posedit = ''
            hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
            return hgvs_transcript_to_hgvs_protein
        else:
            hgvs_transcript_to_hgvs_protein['error'] = 'Unable to map %s to %s' % (
                hgvs_transcript.ac, associated_protein_accession)
            return hgvs_transcript_to_hgvs_protein

