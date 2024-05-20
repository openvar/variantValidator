import os
from configparser import ConfigParser
import vvhgvs
import vvhgvs.parser
import vvhgvs.dataproviders.uta
import vvhgvs.dataproviders.seqfetcher
import vvhgvs.assemblymapper
import vvhgvs.variantmapper
import vvhgvs.sequencevariant
import vvhgvs.validator
import vvhgvs.exceptions
import vvhgvs.location
import vvhgvs.posedit
import vvhgvs.edit
import vvhgvs.normalizer
from Bio.Seq import Seq

import re
import copy
from .vvDatabase import Database
from . import utils
from VariantValidator.settings import CONFIG_DIR
from VariantValidator.version import __version__


class InitialisationError(Exception):
    pass


class Mixin:
    """
    This mixin is the first for the validator object, which is instantiated in order to perform validator functions.
    The validator contains configuration information and permanent copies of database links and the like.
    Much of the validator's inner workings are stored in special one-off function container objects:
    validator.db : The validator's MySQL database access functions

    The validator configuration is loaded when the validator object is initialized.

    Running variant validator should hopefully be as simple as writing a script like this:
    import VariantValidator

    val=Validator()
    val.validate("some kind of gene situation","The genome version","the transcripts to use")

    """
    def __init__(self):
        """
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
        """

        # Load the configuration file.
        config = ConfigParser()
        config.read(CONFIG_DIR)

        # Handle databases
        self.entrez_email = config["Entrez"]["email"]
        self.entrez_api_key = None
        if config['Entrez']['api_key'] != 'YOUR_API_KEY':
            self.entrez_api_key = config['Entrez']['api_key']

        self.seqrepoVersion = config["seqrepo"]["version"]
        self.check_same_thread = config["seqrepo"]["require_threading"]
        if self.check_same_thread == "True":  # This is because the question asked is oposed to the required action
            self.check_same_thread = False
        elif self.check_same_thread == "False":
            self.check_same_thread = True
        self.seqrepoPath = os.path.join(config["seqrepo"]["location"], self.seqrepoVersion)
        self.vvdbVersion = config["mysql"]["version"]

        os.environ['HGVS_SEQREPO_DIR'] = self.seqrepoPath

        os.environ['UTA_DB_URL'] = "postgresql://%s:%s@%s:%s/%s/%s" % (
            config["postgres"]["user"],
            config["postgres"]["password"],
            config['postgres']['host'],
            config['postgres']['port'],
            config['postgres']['database'],
            config['postgres']['version']
        )
        self.utaPath = os.environ.get('UTA_DB_URL')

        self.dbConfig = {
            'user':     config["mysql"]["user"],
            'password': config["mysql"]["password"],
            'host':     config["mysql"]["host"],
            'port':     int(config["mysql"]["port"]),
            'database': config["mysql"]["database"],
            'raise_on_warnings': True
        }
        # Create database access objects
        self.db = Database(self.dbConfig)
        db_version = self.db.get_db_version()
        if db_version[0] != config["mysql"]["version"]:
            raise InitialisationError("Config error: VVDb version in config file is incorrect. VDb version is "
                                      + db_version[0])

        # Set up versions
        self.version = __version__
        if re.match(r'^\d+\.\d+\.\d+$', __version__) is not None:
            self.releasedVersion = True
            _is_released_version = True
        else:
            self.releasedVersion = False
        self.hgvsVersion = vvhgvs.__version__

        # Set up for test mode
        self.testing = False

        # Set up HGVS
        # Configure hgvs package global settings
        vvhgvs.global_config.uta.pool_max = 25
        vvhgvs.global_config.formatting.max_ref_length = 1000000

        # Create HGVS objects
        self.hdp = vvhgvs.dataproviders.uta.connect(pooling=True)
        self.hp = vvhgvs.parser.Parser()  # Parser
        self.vr = vvhgvs.validator.Validator(self.hdp)  # Validator
        self.vm = vvhgvs.variantmapper.VariantMapper(self.hdp)  # Variant mapper
        self.primary_assembly = 'GRCh38'  # Primary assembly defaults to GRCh38

        # Create a lose vm instance
        self.lose_vm = vvhgvs.variantmapper.VariantMapper(self.hdp,
                                                          replace_reference=True,
                                                          prevalidation_level=None
                                                          )

        self.nr_vm = vvhgvs.variantmapper.VariantMapper(self.hdp, replace_reference=False)  # No reverse variant mapper
        self.sf = vvhgvs.dataproviders.seqfetcher.SeqFetcher(self.check_same_thread)  # Seqfetcher

        # Set standard genome builds
        self.genome_builds = ['GRCh37', 'hg19', 'GRCh38']
        self.utaSchema = str(self.hdp.data_version())

        # When we are able to access Ensembl data we will need to use these normalizer instances
        # These are currently implemented in VF
        self.splign_normalizer = vvhgvs.normalizer.Normalizer(
            self.hdp,
            cross_boundaries=False,
            shuffle_direction=vvhgvs.global_config.normalizer.shuffle_direction,
            alt_aln_method='splign'  # RefSeq
            )

        self.genebuild_normalizer = vvhgvs.normalizer.Normalizer(
            self.hdp,
            cross_boundaries=False,
            shuffle_direction=vvhgvs.global_config.normalizer.shuffle_direction,
            alt_aln_method='genebuild'  # Ensembl
            )

        self.genebuild_normalizer_cross = vvhgvs.normalizer.Normalizer(
            self.hdp,
            cross_boundaries=True,
            shuffle_direction=vvhgvs.global_config.normalizer.shuffle_direction,
            alt_aln_method='genebuild'  # Ensembl
            )

        self.reverse_splign_normalizer = vvhgvs.normalizer.Normalizer(self.hdp,
                                                                      cross_boundaries=False,
                                                                      shuffle_direction=5,
                                                                      alt_aln_method='splign' # RefSeq
                                                                      )

        self.reverse_genebuild_normalizer = vvhgvs.normalizer.Normalizer(self.hdp,
                                                                         cross_boundaries=False,
                                                                         shuffle_direction=5,
                                                                         alt_aln_method='genebuild' # Ensembl
                                                                         )

        # Created during validate method
        self.selected_assembly = None
        self.select_transcripts = None
        self.alt_aln_method = None
        self.batch_list = []

    # Create additional normalizers
    def create_additional_normalizers_and_mappers(self):
        self.reverse_hn = vvhgvs.normalizer.Normalizer(self.hdp,
                                                       cross_boundaries=False,
                                                       shuffle_direction=5,
                                                       alt_aln_method=self.alt_aln_method
                                                       )

        self.merge_normalizer = vvhgvs.normalizer.Normalizer(
           self.hdp,
           cross_boundaries=False,
           shuffle_direction=vvhgvs.global_config.normalizer.shuffle_direction,
           alt_aln_method=self.alt_aln_method,
           validate=False
        )

        self.reverse_merge_normalizer = vvhgvs.normalizer.Normalizer(
           self.hdp,
           cross_boundaries=False,
           shuffle_direction=5,
           alt_aln_method=self.alt_aln_method,
           validate=False
        )

        self.no_norm_evm = vvhgvs.assemblymapper.AssemblyMapper(self.hdp,
                                                                assembly_name=self.primary_assembly,
                                                                alt_aln_method=self.alt_aln_method,
                                                                normalize=False,
                                                                replace_reference=True
                                                                )

    def __del__(self):
        try:
            del self.db
        except AttributeError:
            pass

    def my_config(self):
        """
        Returns configuration:
        version, hgvs version, uta schema, seqrepo db.
        """
        return {
            'variantvalidator_version': self.version,
            'variantvalidator_hgvs_version': self.hgvsVersion,
            'vvta_version': self.utaSchema,
            'vvseqrepo_db': self.seqrepoPath,
            'vvdb_version': self.vvdbVersion
        }

    def myc_to_p(self, hgvs_transcript, evm, re_to_p, hn):

        # Create dictionary to store the information
        hgvs_transcript_to_hgvs_protein = {'error': '', 'hgvs_protein': '', 'ref_residues': ''}

        associated_protein_accession = ''
        # Collect the associated protein
        if hgvs_transcript.type == 'c':
            associated_protein_accession = self.hdp.get_pro_ac_for_tx_ac(hgvs_transcript.ac)
            # This method sometimes fails
            if associated_protein_accession is None:
                cod = str(hgvs_transcript)
                cod = cod.replace('inv', 'del')
                cod = self.hp.parse(cod)  # Changed from parse
                p = evm.c_to_p(cod)
                associated_protein_accession = p.ac

        if hgvs_transcript.type == 'c':

            # Handle non inversions with simple c_to_p mapping
            if (hgvs_transcript.posedit.edit.type != 'inv') and (hgvs_transcript.posedit.edit.type != 'dup') and \
                    (hgvs_transcript.posedit.edit.type != 'delins') and (hgvs_transcript.posedit.edit.type != 'sub') and (hgvs_transcript.posedit.edit.type != 'identity') \
                    and (re_to_p is False):
                hgvs_protein = None
                # Does the edit affect the start codon?
                if ((1 <= hgvs_transcript.posedit.pos.start.base <= 3 and hgvs_transcript.posedit.pos.start.offset == 0)
                    or (1 <= hgvs_transcript.posedit.pos.end.base <= 3 and hgvs_transcript.posedit.pos.end.offset
                        == 0)) and '*' not in str(hgvs_transcript.posedit.pos):
                    residue_one = self.sf.fetch_seq(associated_protein_accession, start_i=1 - 1, end_i=1)
                    threed_residue_one = utils.one_to_three(residue_one)
                    r_one_report = '(%s1?)' % threed_residue_one  # was (MET1?)
                    hgvs_protein = vvhgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                          type='p',
                                                                          posedit=r_one_report
                                                                          )

                else:
                    try:
                        hgvs_protein = evm.c_to_p(hgvs_transcript)
                    except IndexError as e:
                        error = str(e)
                        if 'string index out of range' in error and 'dup' in str(hgvs_transcript):
                            hgvs_ins = self.hp.parse(str(hgvs_transcript))
                            hgvs_ins = hn.normalize(hgvs_ins)
                            inst = hgvs_ins.ac + ':c.' + str(hgvs_ins.posedit.pos.start.base - 1) + '_' + \
                                str(hgvs_ins.posedit.pos.start.base) + 'ins' + hgvs_ins.posedit.edit.ref
                            hgvs_transcript = self.hp.parse(inst)
                            hgvs_protein = evm.c_to_p(hgvs_transcript)

                if hgvs_protein:
                    hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                    # Replace Ter<pos>Ter with Ter=
                    if re.search('Ter\d+Ter', str(hgvs_transcript_to_hgvs_protein['hgvs_protein'].posedit)):
                        posedit = str(hgvs_transcript_to_hgvs_protein['hgvs_protein'].posedit)
                        posedit = posedit[:-4] + '=)'
                        hgvs_transcript_to_hgvs_protein['hgvs_protein'].posedit = posedit
                    try:
                        # Sometimes ins create an inline Ter in the alt. Needs to be terminated after the ter
                        if re.search("\*[A-Z]+", hgvs_transcript_to_hgvs_protein['hgvs_protein'].posedit.edit.alt):
                            pr_alt_ter_stp = hgvs_transcript_to_hgvs_protein['hgvs_protein'].posedit.edit.alt
                            pr_alt_ter_stp = pr_alt_ter_stp.split('*')[0] + '*'
                            hgvs_transcript_to_hgvs_protein['hgvs_protein'].posedit.edit.alt = pr_alt_ter_stp
                    except Exception:
                        pass
                    return hgvs_transcript_to_hgvs_protein
                else:
                    hgvs_transcript_to_hgvs_protein = self.myc_to_p(hgvs_transcript, evm, re_to_p=True, hn=hn)
                    return hgvs_transcript_to_hgvs_protein

            else:
                # Additional code required to process inversions
                # Note, this code was developed for VariantValidator and is not native to the biocommons hgvs
                # Python package
                # Convert positions to n. position
                hgvs_naughty = self.vm.c_to_n(hgvs_transcript)

                # Collect the deleted sequence using fetch_seq
                del_seq = self.sf.fetch_seq(str(hgvs_naughty.ac), start_i=hgvs_naughty.posedit.pos.start.base - 1,
                                            end_i=hgvs_naughty.posedit.pos.end.base)

                # Make the inverted sequence
                my_seq = Seq(del_seq)

                if hgvs_transcript.posedit.edit.type == 'inv':
                    inv_seq = my_seq.reverse_complement()
                elif 'del' in hgvs_transcript.posedit.edit.type:
                    inv_seq = hgvs_transcript.posedit.edit.alt
                    if inv_seq is None:
                        inv_seq = ''
                elif 'dup' in hgvs_transcript.posedit.edit.type:
                    inv_seq = del_seq + del_seq
                elif 'sub' in hgvs_transcript.posedit.edit.type:
                    inv_seq = hgvs_transcript.posedit.edit.alt
                elif 'identity' in hgvs_transcript.posedit.edit.type:
                    inv_seq = hgvs_transcript.posedit.edit.ref

                shifts = ''
                # Look for p. delins or del
                not_delins = True
                if hgvs_transcript.posedit.edit.type != 'inv':
                    try:
                        shifts = evm.c_to_p(hgvs_transcript)
                        if "identity" in shifts.posedit.edit.type:
                            not_delins = False
                        if 'del' in shifts.posedit.edit.type or 'dup' in shifts.posedit.edit.type:
                            not_delins = False
                        if "fs" in shifts.posedit.edit.type:
                            not_delins = True
                    except Exception:
                        not_delins = False
                else:
                    not_delins = False

                # Use inv delins code?
                if not not_delins:

                    # Collect the associated protein
                    associated_protein_accession = self.hdp.get_pro_ac_for_tx_ac(hgvs_transcript.ac)

                    # Intronic inversions are marked as uncertain i.e. p.?
                    if re.search(r'\d+-', str(hgvs_transcript.posedit.pos)) \
                            or re.search(r'\d+\+', str(hgvs_transcript.posedit.pos)) \
                            or re.search(r'\*', str(hgvs_transcript.posedit.pos)) \
                            or (re.search(r'[cn].-', str(hgvs_transcript)
                                          ) and "dup" not in hgvs_transcript.posedit.edit.type) or (
                            ("dup" in hgvs_transcript.posedit.edit.type and
                             "-" in str(hgvs_transcript.posedit.pos.end))
                            or
                            ("dup" in hgvs_transcript.posedit.edit.type and
                             "*" in str(hgvs_transcript.posedit.pos.start))
                            ):

                        if ((1 <= hgvs_transcript.posedit.pos.start.base <= 3 and
                            hgvs_transcript.posedit.pos.start.offset == 0) or (1 <=
                            hgvs_transcript.posedit.pos.end.base <= 3 and hgvs_transcript.posedit.pos.end.offset == 0))\
                                and '*' not in str(hgvs_transcript.posedit.pos):

                            residue_one = self.sf.fetch_seq(associated_protein_accession, start_i=1 - 1, end_i=1)
                            threed_residue_one = utils.one_to_three(residue_one)
                            r_one_report = '(%s1?)' % threed_residue_one  # was (MET1?)
                            hgvs_protein = vvhgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                                  type='p', posedit=r_one_report)
                        else:
                            # Make the variant
                            hgvs_protein = vvhgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                                  type='p', posedit='?')
                        hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                        return hgvs_transcript_to_hgvs_protein
                    else:
                        # Need to obtain the cds_start
                        inf = self.hdp.get_tx_identity_info(hgvs_transcript.ac)
                        cds_start = inf[3]
                        cds_end = inf[4]

                        # Extract the reference coding sequence from SeqRepo
                        try:
                            ref_seq = self.sf.fetch_seq(str(hgvs_naughty.ac))
                        except Exception as e:
                            error = str(e)
                            hgvs_transcript_to_hgvs_protein['error'] = error
                            return hgvs_transcript_to_hgvs_protein

                        # Create the variant coding sequence
                        var_seq = utils.n_inversion(ref_seq, del_seq, inv_seq,
                                                    hgvs_naughty.posedit.pos.start.base,
                                                    hgvs_naughty.posedit.pos.end.base)

                        # Check for modified amino acids
                        prot_seq = self.sf.fetch_seq(associated_protein_accession)
                        if "U" in prot_seq:
                            modified_aa = "Sec"
                        else:
                            modified_aa = None

                        # Translate the reference and variant proteins
                        prot_ref_seq = utils.translate(ref_seq, cds_start, modified_aa)
                        try:
                            prot_var_seq = utils.translate(var_seq, cds_start, modified_aa)
                        except IndexError:
                            hgvs_transcript_to_hgvs_protein['error'] = \
                                'Cannot identify an in-frame Termination codon in the variant mRNA sequence'
                            hgvs_protein = vvhgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                                  type='p', posedit='?')
                            hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                            return hgvs_transcript_to_hgvs_protein

                        if prot_ref_seq == 'error':
                            error = 'Unable to generate protein variant description'
                            hgvs_transcript_to_hgvs_protein['error'] = error
                            return hgvs_transcript_to_hgvs_protein
                        elif prot_var_seq == 'error':
                            # Does the edit affect the start codon?
                            if ((1 <= hgvs_transcript.posedit.pos.start.base <= 3 and
                                 hgvs_transcript.posedit.pos.start.offset == 0) or (
                                    1 <= hgvs_transcript.posedit.pos.end.base <= 3 and
                                    hgvs_transcript.posedit.pos.end.offset == 0)) \
                                    and '*' not in str(hgvs_transcript.posedit.pos):
                                residue_one = self.sf.fetch_seq(associated_protein_accession, start_i=1 - 1, end_i=1)
                                threed_residue_one = utils.one_to_three(residue_one)
                                r_one_report = '(%s1?)' % threed_residue_one  # was (MET1?)
                                hgvs_protein = vvhgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                                      type='p', posedit=r_one_report)

                                hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                                return hgvs_transcript_to_hgvs_protein
                            else:
                                error = 'Unable to generate protein variant description'
                                hgvs_transcript_to_hgvs_protein['error'] = error
                                return hgvs_transcript_to_hgvs_protein
                        elif ((1 <= hgvs_transcript.posedit.pos.start.base <= 3 and
                            hgvs_transcript.posedit.pos.start.offset == 0) or (1 <=
                            hgvs_transcript.posedit.pos.end.base <= 3 and hgvs_transcript.posedit.pos.end.offset == 0))\
                                and '*' not in str(hgvs_transcript.posedit.pos):

                            residue_one = self.sf.fetch_seq(associated_protein_accession, start_i=1 - 1, end_i=1)
                            threed_residue_one = utils.one_to_three(residue_one)
                            r_one_report = '(%s1?)' % threed_residue_one  # was (MET1?)
                            hgvs_protein = vvhgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                                  type='p', posedit=r_one_report)
                            hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                            return hgvs_transcript_to_hgvs_protein
                        else:

                            # Gather the required information regarding variant interval and sequences
                            if hgvs_transcript.posedit.edit.type != 'delins' and \
                                    hgvs_transcript.posedit.edit.type != 'dup':
                                pro_inv_info = utils.pro_inv_info(prot_ref_seq, prot_var_seq)
                            else:
                                # Test whether the length of the deletion, plus the insertion can be divided by 3
                                # This is trying to spot the difference between amino acid deletions
                                # and early terminations

                                # Get the cds length
                                cds_len = cds_end - cds_start

                                # Calculate the variant cds length
                                minus = False
                                plus = False

                                try:
                                    if len(hgvs_naughty.posedit.edit.ref) > len(hgvs_naughty.posedit.edit.alt):
                                        var_cds_len = cds_len - (len(hgvs_naughty.posedit.edit.ref)
                                                                 - len(hgvs_naughty.posedit.edit.alt))
                                        minus = True
                                    elif len(hgvs_naughty.posedit.edit.ref) < len(hgvs_naughty.posedit.edit.alt):
                                        var_cds_len = cds_len + (len(hgvs_naughty.posedit.edit.alt)
                                                                 - len(hgvs_naughty.posedit.edit.ref))
                                        plus = True
                                except AttributeError as e:
                                    if "'Dup' object has no attribute 'alt'" in str(e):
                                        var_cds_len = cds_len + (len(var_seq)
                                                                 - len(ref_seq))
                                        plus = True

                                # Do we have an in-frame variant i.e. divisible by 3?
                                in_frame = False
                                if minus is True:
                                    loss_gain = (cds_len) - (var_cds_len)
                                    if loss_gain % 3 == 0:
                                        loss_gain = loss_gain / 3
                                        loss_gain = 0 - loss_gain
                                        in_frame = loss_gain
                                elif plus is True:
                                    loss_gain = var_cds_len - cds_len
                                    if loss_gain % 3 == 0:
                                        loss_gain = loss_gain / 3
                                        in_frame = loss_gain

                                # Get the sequence info
                                pro_inv_info = utils.pro_delins_info(prot_ref_seq,
                                                                     prot_var_seq,
                                                                     in_frame)

                            # Error has occurred
                            if pro_inv_info['error'] == 'true':
                                error = 'Translation error occurred, please contact admin'
                                hgvs_transcript_to_hgvs_protein['error'] = error
                                return hgvs_transcript_to_hgvs_protein

                            # The Nucleotide variant has not affected the protein sequence i.e. synonymous
                            elif pro_inv_info['variant'] != 'true':

                                # Make the variant
                                hgvs_protein = vvhgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                                      type='p', posedit='=')
                                # Where possible, identify the exact positions of the amino acids
                                if isinstance(hgvs_transcript.posedit.pos.start.base, int) and isinstance(
                                        hgvs_transcript.posedit.pos.end.base, int):

                                    aa_start_pos = float(hgvs_transcript.posedit.pos.start.base / 3)
                                    aa_end_pos = float(hgvs_transcript.posedit.pos.end.base / 3)

                                    # end pos may be in the next amino acid i.e. float>0
                                    if not aa_end_pos.is_integer():
                                        aa_end_pos = int(aa_end_pos + 1)
                                    else:
                                        aa_end_pos = int(aa_end_pos)
                                    if not aa_start_pos.is_integer():
                                        aa_start_pos = int(aa_start_pos + 1)
                                    else:
                                        aa_start_pos = int(aa_start_pos)

                                    aa_seq = self.sf.fetch_seq(associated_protein_accession, start_i=aa_start_pos - 1,
                                                               end_i=aa_end_pos)

                                    # Handle Termination unaffected (note, * does not appear in the reference sequence)
                                    if aa_seq == "":
                                        ck_aa_seq = self.sf.fetch_seq(associated_protein_accession)
                                        length = len(ck_aa_seq)
                                        if aa_start_pos == length + 1 and aa_end_pos == length + 1:
                                            aa_seq = "*"

                                    start_aa = utils.one_to_three(aa_seq[0])
                                    end_aa = utils.one_to_three(aa_seq[-1])

                                    # create edit
                                    if aa_start_pos != aa_end_pos:
                                        posedit = '(%s%s_%s%s=)' % (start_aa,
                                                                    str(aa_start_pos),
                                                                    end_aa,
                                                                    str(aa_end_pos)
                                                                    )

                                        hgvs_protein = vvhgvs.sequencevariant.SequenceVariant(
                                            ac=associated_protein_accession, type='p', posedit=posedit)
                                    else:
                                        posedit = '(%s%s=)' % (start_aa, str(aa_start_pos))
                                        hgvs_protein = vvhgvs.sequencevariant.SequenceVariant(
                                            ac=associated_protein_accession, type='p', posedit=posedit)

                                hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein
                                return hgvs_transcript_to_hgvs_protein

                            else:

                                # Adjust extended aas if necessary
                                if modified_aa == "Sec":
                                    if "U" in pro_inv_info['prot_ins_seq'] and "U" not in pro_inv_info['prot_del_seq']:
                                        pro_inv_info['prot_ins_seq'] = pro_inv_info['prot_ins_seq'].replace("U", "*")
                                        pro_inv_info['ter_pos'] = pro_inv_info['edit_start'] + len(
                                            pro_inv_info['prot_ins_seq'].split("*")[0])

                                # Early termination i.e. stop gained
                                if pro_inv_info['terminate'] == 'true' and \
                                        (hgvs_transcript.posedit.edit.type == 'delins' or
                                         hgvs_transcript.posedit.edit.type == 'dup' or
                                         hgvs_transcript.posedit.edit.type == 'inv'):

                                    # This deals with early terminating delins in-frame prventing the format
                                    # NP_733765.1:p.(Gln259_Ser1042delinsProAla*) in issue #214 also #282
                                    if len(pro_inv_info['prot_del_seq']) + \
                                            int(pro_inv_info['edit_start'] - 1) == int(pro_inv_info['ter_pos']):
                                        end = 'Ter' + str(pro_inv_info['ter_pos'])
                                        pro_inv_info['prot_ins_seq'].replace('*', end)
                                        pro_inv_info['prot_ins_seq'] = pro_inv_info['prot_ins_seq']
                                        pro_inv_info['prot_del_seq'] = pro_inv_info['prot_del_seq'][0]
                                        pro_inv_info['edit_end'] = pro_inv_info['edit_start']
                                    elif hgvs_transcript.posedit.edit.type == 'dup' and pro_inv_info["prot_del_seq"] \
                                            == "" and (int(pro_inv_info["edit_end"]) < int(pro_inv_info["edit_start"])):

                                        # Handles in-frame dups only
                                        dup_len = (int(hgvs_transcript.posedit.pos.end.base) - int(
                                            hgvs_transcript.posedit.pos.start.base) + 1) / 3
                                        pro_inv_info['prot_del_seq'] = pro_inv_info['prot_ins_seq']
                                        pro_inv_info['edit_start'] = pro_inv_info['edit_end'] - \
                                                                     len(pro_inv_info['prot_del_seq']) + 1
                                        start_aa = self.sf.fetch_seq(associated_protein_accession,
                                                                     int(pro_inv_info['edit_start']-1),
                                                                     int(pro_inv_info['edit_start']) + (dup_len -1))
                                        pro_inv_info['prot_del_seq'] = start_aa
                                        pro_inv_info['prot_ins_seq'] = start_aa + \
                                                                       pro_inv_info['prot_ins_seq']

                                # Complete variant description
                                # Recode the single letter del and ins sequences into three letter amino acid codes
                                del_thr = utils.one_to_three(pro_inv_info['prot_del_seq'])
                                ins_thr = utils.one_to_three(pro_inv_info['prot_ins_seq'])

                                # Write the HGVS position and edit
                                del_len = len(del_thr)
                                from_aa = del_thr[0:3]
                                to_aa = del_thr[del_len - 3:]

                                # Handle a range of amino acids
                                if pro_inv_info['edit_start'] != pro_inv_info['edit_end']:

                                    # Handle duplications
                                    if pro_inv_info["prot_ins_seq"] == (pro_inv_info["prot_del_seq"]
                                                                          + pro_inv_info["prot_del_seq"]):

                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + \
                                                  str(pro_inv_info['edit_end']) + 'dup)'
                                    elif len(ins_thr) > 0:
                                        if 'Ter' in del_thr and ins_thr[-3:] != 'Ter':
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + \
                                                      str(pro_inv_info['edit_end']) + 'delins' + ins_thr + '?)'

                                        elif len(pro_inv_info["prot_ins_seq"]) > len(pro_inv_info["prot_del_seq"]) \
                                                and pro_inv_info["prot_ins_seq"] != (pro_inv_info["prot_del_seq"]
                                                                                     + pro_inv_info["prot_del_seq"]) and \
                                                pro_inv_info["prot_del_seq"] == "" and (pro_inv_info["edit_start"]
                                                > pro_inv_info["edit_end"]):

                                            from_aa = self.sf.fetch_seq(associated_protein_accession,
                                                                      int(pro_inv_info['edit_end']-len(pro_inv_info['prot_ins_seq'])),
                                                                      int(pro_inv_info['edit_end']-len(pro_inv_info['prot_ins_seq']))+1)

                                            to_aa = self.sf.fetch_seq(associated_protein_accession,
                                                                      int(pro_inv_info['edit_start']-2),
                                                                      int(pro_inv_info['edit_start']-1))

                                            posedit = '(' + from_aa + str(pro_inv_info['edit_end']-len(pro_inv_info['prot_ins_seq'])+1) + "_" + \
                                                      to_aa + str(pro_inv_info['edit_start']-1) + "dup)"

                                        else:
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + \
                                                      str(pro_inv_info['edit_end']) + 'delins' + ins_thr + ')'

                                    else:
                                        if 'Ter' in del_thr and ins_thr[-3:] != 'Ter':
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + \
                                                      str(pro_inv_info['edit_end']) + 'del?)'
                                        else:
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + '_' + to_aa + \
                                                      str(pro_inv_info['edit_end']) + 'del)'
                                else:
                                    # Handle duplications
                                    if pro_inv_info["prot_ins_seq"] == (pro_inv_info["prot_del_seq"]
                                                                          + pro_inv_info["prot_del_seq"]):
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + 'dup)'

                                    # Handle insertions
                                    elif len(pro_inv_info["prot_ins_seq"]) > len(pro_inv_info["prot_del_seq"]) \
                                        and pro_inv_info["prot_ins_seq"] != (pro_inv_info["prot_del_seq"]
                                                                          + pro_inv_info["prot_del_seq"]) and \
                                            (pro_inv_info["prot_ins_seq"][0] == pro_inv_info["prot_del_seq"][0]):

                                        to_aa = self.sf.fetch_seq(associated_protein_accession,
                                                                     int(pro_inv_info['edit_start']),
                                                                     int(pro_inv_info['edit_start'] + 1))
                                        posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + "_" + \
                                            to_aa + str(pro_inv_info['edit_start'] + 1) + "ins" + \
                                                  pro_inv_info["prot_ins_seq"][1:] + ")"

                                    # Handle extended proteins i.e. stop_lost
                                    elif del_thr == 'Ter' and (len(ins_thr) > len(del_thr)):
                                        # Nucleotide variant range aligns to the Termination codon
                                        if ins_thr[-3:] == 'Ter':
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + str(
                                                ins_thr[:3]) + 'ext' + str(ins_thr[-3:]) + str(int((len(ins_thr) / 3))
                                                                                               - 1) + ')'
                                        # Nucleotide variant range spans the Termination codon
                                        else:
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + str(
                                                ins_thr[:3]) + 'ext?)'

                                    # Nucleotide variation has not affected the length of the protein thus
                                    # substitution or del
                                    else:
                                        if len(ins_thr) == 3:
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + ins_thr + ')'
                                        elif len(ins_thr) == 0:
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + 'del)'
                                        else:
                                            posedit = '(' + from_aa + str(pro_inv_info['edit_start']) + 'delins' + \
                                                      ins_thr + ')'

                                # Complete the variant
                                hgvs_protein = vvhgvs.sequencevariant.SequenceVariant(ac=associated_protein_accession,
                                                                                      type='p',
                                                                                      posedit=posedit
                                                                                      )

                                hgvs_transcript_to_hgvs_protein['hgvs_protein'] = hgvs_protein

                else:
                    hgvs_transcript_to_hgvs_protein['hgvs_protein'] = shifts

                # Replace Ter<pos>Ter with Ter=
                if re.search('Ter\d+Ter', str(hgvs_transcript_to_hgvs_protein['hgvs_protein'].posedit)):
                    posedit = str(hgvs_transcript_to_hgvs_protein['hgvs_protein'].posedit)
                    posedit = posedit[:-4] + '=)'
                    hgvs_transcript_to_hgvs_protein['hgvs_protein'].posedit = posedit

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
