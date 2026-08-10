# -*- coding: utf-8 -*-

import logging
import os
from functools import lru_cache
from configparser import ConfigParser

import vvhgvs
import vvhgvs.assemblymapper
import vvhgvs.dataproviders.seqfetcher
import vvhgvs.dataproviders.uta
import vvhgvs.edit
import vvhgvs.exceptions
import vvhgvs.location
import vvhgvs.normalizer
import vvhgvs.parser
import vvhgvs.posedit
import vvhgvs.sequencevariant
import vvhgvs.validator
import vvhgvs.variantmapper
from Bio.Seq import Seq
from vvhgvs.edit import AAExt, AAFs, AARefAlt, Dup
from vvhgvs.location import AAPosition, Interval
from vvhgvs.posedit import PosEdit

from VariantValidator import settings, version
from VariantValidator.modules.hgvs_utils import (
    VVPosEdit,
    hgvs_delins_parts_to_hgvs_obj,
)

from . import hgvs_position_utils
from . import utils
from .vvDatabase import Database


logger = logging.getLogger(__name__)


class InitialisationError(Exception):
    """Raised when the Validator cannot be initialised."""
    pass


class CachedSeqFetcher:
    """
    LRU cache wrapper for SeqFetcher.

    This class implements the Decorator pattern around SeqFetcher.
    Only the most frequently repeated sequence lookup is cached using
    functools.lru_cache(); all other methods and attributes are
    transparently delegated to the wrapped SeqFetcher.

    The cache exists solely within the current Python process and is
    discarded when the process exits.
    """

    def __init__(self, sf):
        """
        Wrap an existing SeqFetcher.

        Parameters
        ----------
        sf
            An instantiated SeqFetcher object.
        """
        self._sf = sf

    def __getattr__(self, name):
        """
        Delegate all uncached methods and attributes to the wrapped
        SeqFetcher.

        Python only calls __getattr__ when an attribute is not found on
        CachedSeqFetcher itself. Consequently, only the methods
        explicitly implemented below are intercepted and cached;
        everything else behaves exactly as if the original SeqFetcher
        were being used directly.
        """
        return getattr(self._sf, name)

    @lru_cache(maxsize=settings.SEQFETCHER_CACHE_SIZE)
    def fetch_seq(self, ac, start_i=None, end_i=None):
        """
        Retrieve a sequence or sequence slice.

        The cache key is formed from the accession, start coordinate
        and end coordinate. Repeated requests for the same sequence
        slice are therefore served directly from memory rather than
        performing another SeqRepo lookup.
        """
        return self._sf.fetch_seq(ac, start_i, end_i)


class Mixin:
    """
    Initialise the persistent VariantValidator infrastructure.

    The Validator owns configuration information, database access,
    sequence access and HGVS objects that persist across individual
    variant validations.
    """

    def __init__(self):
        """
        Initialise Validator configuration and persistent infrastructure.

        Historical variable names
        -------------------------
        seqrepo_directory:
            HGVS_SEQREPO_DIR
        uta_url:
            UTA_DB_URL
        py_liftover_directory:
            PYLIFTOVER_DIR
        variantvalidator_data_url:
            VALIDATOR_DB_URL
        entrez_id:
            ENTREZ_ID
        variantvalidator_version:
            VERSION
        variantvalidator_hgvs_version:
            hgvs_version
        uta_schema:
            hdp.data_version()
        seqrepo_db:
            HGVS_SEQREPO_DIR.split('/')[-1]
        """

        # --------------------------------------------------------------
        # HGVS global configuration
        # --------------------------------------------------------------

        vvhgvs.global_config.uta.pool_max = 25
        vvhgvs.global_config.formatting.max_ref_length = 1000000

        if settings.vvHGVS_HDP_CACHE:
            vvhgvs.global_config.lru_cache.maxsize = settings.vvHGVS_HDP_CACHE_SIZE
        else:
            vvhgvs.global_config.lru_cache.maxsize = 200  # Default vvHGVS cache size.

        # --------------------------------------------------------------
        # Configuration
        # --------------------------------------------------------------

        config_path = settings.get_config_dir()

        if not os.path.exists(config_path):
            logger.error(
                "Configuration file not found at %s",
                config_path,
            )
            raise InitialisationError(
                "Configuration file not found, please create a new one at "
                f"{config_path}"
            )

        config = ConfigParser()
        config.read(config_path)

        logger.info(
            "Configuration file loaded from %s",
            config_path,
        )

        # --------------------------------------------------------------
        # Entrez configuration
        # --------------------------------------------------------------

        self.entrez_email = config["Entrez"]["email"]
        self.entrez_api_key = None

        api_key = config["Entrez"]["api_key"]

        if api_key != "YOUR_API_KEY":
            self.entrez_api_key = api_key

        # --------------------------------------------------------------
        # SeqRepo configuration
        # --------------------------------------------------------------

        self.seqrepoVersion = config["seqrepo"]["version"]

        require_threading = config["seqrepo"]["require_threading"]

        # The configuration question is opposite to the action required
        # by SeqFetcher.
        if require_threading == "True":
            self.check_same_thread = False
        elif require_threading == "False":
            self.check_same_thread = True
        else:
            self.check_same_thread = require_threading

        self.seqrepoPath = os.path.join(
            config["seqrepo"]["location"],
            self.seqrepoVersion,
        )

        os.environ["HGVS_SEQREPO_DIR"] = self.seqrepoPath

        # --------------------------------------------------------------
        # UTA configuration
        # --------------------------------------------------------------

        psql_host_or_socketfile = (
            config["postgres"]["host"].replace("/", "%2F")
        )

        os.environ["UTA_DB_URL"] = (
            "postgresql://%s:%s@%s:%s/%s/%s"
            % (
                config["postgres"]["user"],
                config["postgres"]["password"],
                psql_host_or_socketfile,
                config["postgres"]["port"],
                config["postgres"]["database"],
                config["postgres"]["version"],
            )
        )

        self.utaPath = os.environ["UTA_DB_URL"]

        # --------------------------------------------------------------
        # VariantValidator database
        # --------------------------------------------------------------

        self.vvdbVersion = config["mysql"]["version"]

        self.dbConfig = {
            "user": config["mysql"]["user"],
            "password": config["mysql"]["password"],
            "host": config["mysql"]["host"],
            "port": int(config["mysql"]["port"]),
            "database": config["mysql"]["database"],
            "raise_on_warnings": True,
        }

        mysql_unix_socket = config.get(
            "mysql",
            "unix_socket",
            fallback=False,
        )

        if mysql_unix_socket:
            self.dbConfig["unix_socket"] = mysql_unix_socket

        self.db = Database(self.dbConfig)

        db_version = self.db.get_db_version()

        if db_version[0] != self.vvdbVersion:
            raise InitialisationError(
                "Config error: VVDb version in config file is incorrect. "
                f"VDb version is {db_version[0]}"
            )

        # --------------------------------------------------------------
        # Version information
        # --------------------------------------------------------------

        self.version = version.__version__
        self.releasedVersion = version._is_released_version
        self.hgvsVersion = vvhgvs.__version__

        # --------------------------------------------------------------
        # Persistent Validator state
        # --------------------------------------------------------------

        self.testing = False
        self.primary_assembly = "GRCh38"

        self.genome_builds = [
            "GRCh37",
            "hg19",
            "GRCh38",
        ]

        # Populated during validation.
        self.selected_assembly = None
        self.select_transcripts = None
        self.alt_aln_method = None

        # Populated by create_additional_normalizers_and_mappers().
        self.reverse_hn = None
        self.hn = None
        self.merge_normalizer = None
        self.reverse_merge_normalizer = None
        self.no_norm_evm = None

        # --------------------------------------------------------------
        # HGVS data providers
        # --------------------------------------------------------------

        # Create the HGVS data provider.
        self.hdp = vvhgvs.dataproviders.uta.connect(pooling=True)

        self.utaSchema = str(
            self.hdp.data_version()
        )

        # --------------------------------------------------------------
        # HGVS parser and validator
        # --------------------------------------------------------------

        self.hp = vvhgvs.parser.Parser(
            expose_all_rules=True,
        )

        self.vr = vvhgvs.validator.Validator(
            self.hdp,
        )

        # --------------------------------------------------------------
        # Variant mappers
        # --------------------------------------------------------------

        self.vm = vvhgvs.variantmapper.VariantMapper(
            self.hdp,
        )

        self.lose_vm = vvhgvs.variantmapper.VariantMapper(
            self.hdp,
            replace_reference=True,
            prevalidation_level=None,
        )

        self.nr_vm = vvhgvs.variantmapper.VariantMapper(
            self.hdp,
            replace_reference=False,
        )

        # --------------------------------------------------------------
        # Sequence provider
        # --------------------------------------------------------------

        self.sf = vvhgvs.dataproviders.seqfetcher.SeqFetcher(
            self.check_same_thread,
        )

        # Wrap the SeqFetcher with the LRU cache layer.
        #
        # CachedSeqFetcher intercepts sequence fetches and caches the returned
        # sequence slices. All other methods and attributes are transparently
        # delegated to the original SeqFetcher via __getattr__().
        #
        # To disable all SeqFetcher caching, simply comment out the line below.
        if settings.SEQFETCHER_CACHE:
            self.sf = CachedSeqFetcher(self.sf)

        # --------------------------------------------------------------
        # Persistent alignment-specific normalizers
        # --------------------------------------------------------------

        shuffle_direction = (
            vvhgvs.global_config.normalizer.shuffle_direction
        )

        self.splign_normalizer = vvhgvs.normalizer.Normalizer(
            self.hdp,
            cross_boundaries=False,
            shuffle_direction=shuffle_direction,
            alt_aln_method="splign",
        )

        self.genebuild_normalizer = vvhgvs.normalizer.Normalizer(
            self.hdp,
            cross_boundaries=False,
            shuffle_direction=shuffle_direction,
            alt_aln_method="genebuild",
        )

        self.genebuild_normalizer_cross = (
            vvhgvs.normalizer.Normalizer(
                self.hdp,
                cross_boundaries=True,
                shuffle_direction=shuffle_direction,
                alt_aln_method="genebuild",
            )
        )

        self.reverse_splign_normalizer = (
            vvhgvs.normalizer.Normalizer(
                self.hdp,
                cross_boundaries=False,
                shuffle_direction=5,
                alt_aln_method="splign",
            )
        )

        self.reverse_genebuild_normalizer = (
            vvhgvs.normalizer.Normalizer(
                self.hdp,
                cross_boundaries=False,
                shuffle_direction=5,
                alt_aln_method="genebuild",
            )
        )

    def create_additional_normalizers_and_mappers(self):
        """
        Create HGVS objects that depend on the selected alignment method
        and primary assembly.
        """
        self.reverse_hn = vvhgvs.normalizer.Normalizer(
            self.hdp,
            cross_boundaries=False,
            shuffle_direction=5,
            alt_aln_method=self.alt_aln_method,
        )

        self.hn = vvhgvs.normalizer.Normalizer(
            self.hdp,
            cross_boundaries=False,
            shuffle_direction=3,
            alt_aln_method=self.alt_aln_method,
        )

        self.merge_normalizer = vvhgvs.normalizer.Normalizer(
            self.hdp,
            cross_boundaries=False,
            shuffle_direction=(
                vvhgvs.global_config.normalizer.shuffle_direction
            ),
            alt_aln_method=self.alt_aln_method,
            validate=False,
        )

        self.reverse_merge_normalizer = (
            vvhgvs.normalizer.Normalizer(
                self.hdp,
                cross_boundaries=False,
                shuffle_direction=5,
                alt_aln_method=self.alt_aln_method,
                validate=False,
            )
        )

        self.no_norm_evm = vvhgvs.assemblymapper.AssemblyMapper(
            self.hdp,
            assembly_name=self.primary_assembly,
            alt_aln_method=self.alt_aln_method,
            normalize=False,
            replace_reference=True,
        )

    def __del__(self):
        if getattr(self, "pool", None):
            self.pool = None

    def my_config(self):
        """
        Return VariantValidator configuration/version information.
        """
        return {
            "variantvalidator_version": self.version,
            "variantvalidator_hgvs_version": self.hgvsVersion,
            "vvta_version": self.utaSchema,
            "vvseqrepo_db": self.seqrepoPath,
            "vvdb_version": self.vvdbVersion,
        }

    def myc_to_p(self, hgvs_transcript, evm, re_to_p, hn):
        logger.info(
            "Translating %s to with myc_to_p",
            hgvs_transcript,
        )

        hgvs_transcript_to_hgvs_protein = {
            "error": "",
            "hgvs_protein": "",
            "ref_residues": "",
        }

        # Handle non-coding transcript and non-transcript descriptions.
        if hgvs_transcript.type == "n":
            return hgvs_transcript_to_hgvs_protein

        if hgvs_transcript.type != "c":
            hgvs_transcript_to_hgvs_protein["error"] = (
                f"Unable to map {hgvs_transcript.ac} "
                "to an associated protein"
            )
            return hgvs_transcript_to_hgvs_protein

        edit = hgvs_transcript.posedit.edit
        pos = hgvs_transcript.posedit.pos
        edit_type = edit.type

        associated_protein_accession = (
            self.hdp.get_pro_ac_for_tx_ac(
                hgvs_transcript.ac
            )
        )

        # This method sometimes fails.
        if associated_protein_accession is None:
            cod = hgvs_delins_parts_to_hgvs_obj(
                hgvs_transcript.ac,
                hgvs_transcript.type,
                pos,
                "",
                "",
            )
            p = evm.c_to_p(cod)
            associated_protein_accession = p.ac

        nucleotide_not_equal = edit_type != "identity"

        def _fb_unc(prot, base):
            return vvhgvs.sequencevariant.SequenceVariant(
                ac=prot,
                type="p",
                posedit=VVPosEdit(
                    pos=Interval(
                        start=AAPosition(
                            base=1,
                            aa=base,
                        )
                    ),
                    edit="",
                    uncertain=True,
                ),
            )

        def _tot_unc(prot):
            return vvhgvs.sequencevariant.SequenceVariant(
                ac=prot,
                type="p",
                posedit=VVPosEdit(
                    pos=Interval(),
                    edit="",
                    uncertain=True,
                ),
            )

        def _remake_unc(
                prot,
                nucleotide_not_equal=False,
        ):
            if prot.posedit is None:
                return prot

            return vvhgvs.sequencevariant.SequenceVariant(
                ac=prot.ac,
                type="p",
                posedit=VVPosEdit(
                    pos=prot.posedit.pos,
                    edit=prot.posedit.edit,
                    uncertain=True,
                    nucleotide_not_equal=nucleotide_not_equal,
                ),
            )

        if (
                edit_type not in (
                    "inv",
                    "dup",
                    "delins",
                    "sub",
                    "identity",
                    "del",
                    "ins",
                )
                and not re_to_p
        ):
            logger.info(
                "Passing %s into simple c_to_p mapping",
                hgvs_transcript,
            )

            hgvs_protein = None

            if (
                    (
                        1 <= pos.start.base <= 3
                        and not (
                            hgvs_position_utils
                            .start_position_is_intronic(
                                hgvs_transcript
                            )
                        )
                    )
                    or (
                        1 <= pos.end.base <= 3
                        and not (
                            hgvs_position_utils
                            .end_position_is_intronic(
                                hgvs_transcript
                            )
                        )
                    )
            ) and not (
                    hgvs_position_utils.start_is_3_prime_utr(
                        hgvs_transcript
                    )
                    or hgvs_position_utils.end_is_3_prime_utr(
                        hgvs_transcript
                    )
            ):
                residue_one = self.sf.fetch_seq(
                    associated_protein_accession,
                    start_i=0,
                    end_i=1,
                )

                hgvs_protein = _fb_unc(
                    associated_protein_accession,
                    residue_one,
                )

            else:
                try:
                    hgvs_protein = evm.c_to_p(
                        hgvs_transcript
                    )

                    hgvs_protein = _remake_unc(
                        hgvs_protein,
                        nucleotide_not_equal=(
                            nucleotide_not_equal
                        ),
                    )

                except IndexError as e:
                    if (
                            "string index out of range" in str(e)
                            and edit_type == "dup"
                    ):
                        hgvs_ins = hn.normalize(
                            hgvs_transcript
                        )

                        hgvs_transcript = (
                            hgvs_delins_parts_to_hgvs_obj(
                                hgvs_transcript.ac,
                                hgvs_transcript.type,
                                pos.start.base - 1,
                                "",
                                hgvs_ins.posedit.edit.ref,
                            )
                        )

                        hgvs_protein = evm.c_to_p(
                            hgvs_transcript
                        )

                        hgvs_protein = _remake_unc(
                            hgvs_protein,
                            nucleotide_not_equal=(
                                nucleotide_not_equal
                            ),
                        )

            if (
                    hgvs_protein
                    and hgvs_protein.posedit is None
            ):
                hgvs_protein = _tot_unc(
                    hgvs_protein.ac
                )

            if hgvs_protein:
                hgvs_transcript_to_hgvs_protein[
                    "hgvs_protein"
                ] = hgvs_protein

                try:
                    protein_alt = (
                        hgvs_protein.posedit.edit.alt
                    )

                    if "*" in protein_alt:
                        head, _, tail = protein_alt.partition("*")
                        if tail[:1].isupper():
                            protein_alt = head + "*"
                            hgvs_protein.posedit.edit.alt = protein_alt

                except Exception:
                    pass

            else:
                hgvs_transcript_to_hgvs_protein = (
                    self.myc_to_p(
                        hgvs_transcript,
                        evm,
                        re_to_p=True,
                        hn=hn,
                    )
                )

            return hgvs_transcript_to_hgvs_protein

        logger.info(
            "Passing %s into VV handled c_to_p mapping",
            hgvs_transcript,
        )

        hgvs_naughty = self.vm.c_to_n(
            hgvs_transcript
        )
        naughty_edit = hgvs_naughty.posedit.edit
        naughty_pos = hgvs_naughty.posedit.pos

        del_seq = self.sf.fetch_seq(
            hgvs_naughty.ac,
            start_i=naughty_pos.start.base - 1,
            end_i=naughty_pos.end.base,
        )

        if edit_type == "inv":
            inv_seq = str(
                Seq(del_seq).reverse_complement()
            )

        elif edit_type in ("del", "delins"):
            inv_seq = edit.alt or ""

        elif edit_type == "dup":
            inv_seq = del_seq + del_seq

        elif edit_type == "sub":
            inv_seq = edit.alt

        elif edit_type == "identity":
            inv_seq = edit.ref

        elif edit_type == "ins":
            inv_seq = (
                f"{del_seq[0]}"
                f"{edit.alt}"
                f"{del_seq[-1]}"
            )

        logger.info(
            "delSeq: %s and insSeq: %s extracted from %s",
            del_seq,
            inv_seq,
            hgvs_transcript,
        )

        shifts = ""
        not_delins = False

        try:
            shifts = evm.c_to_p(
                hgvs_transcript
            )

            shifts = _remake_unc(
                shifts,
                nucleotide_not_equal=nucleotide_not_equal,
            )

            if edit_type in ("inv", "delins"):
                if shifts.posedit.edit.type in (
                        "ins",
                        "fs",
                        "ext",
                ):
                    not_delins = True

            elif shifts.posedit.edit.type in (
                    "ins",
                    "sub",
                    "fs",
                    "ext",
            ):
                not_delins = True

        except Exception:
            not_delins = False

        if not_delins:
            hgvs_transcript_to_hgvs_protein[
                "hgvs_protein"
            ] = shifts

            return hgvs_transcript_to_hgvs_protein

        associated_protein_accession = (
            self.hdp.get_pro_ac_for_tx_ac(
                hgvs_transcript.ac
            )
        )

        logger.info(
            "Test for intronic and UTR trapping"
        )

        structural_edit = edit_type in (
            "dup",
            "del",
            "inv",
            "ins",
            "delins",
        )

        if (
                hgvs_position_utils.either_position_is_intronic(
                    hgvs_transcript
                )
                or (
                    structural_edit
                    and hgvs_position_utils.end_is_5_prime_utr(
                        hgvs_transcript
                    )
                )
                or (
                    structural_edit
                    and hgvs_position_utils.start_is_3_prime_utr(
                        hgvs_transcript
                    )
                )
                or (
                    hgvs_position_utils.start_is_3_prime_utr(
                        hgvs_transcript
                    )
                    and hgvs_position_utils.end_is_3_prime_utr(
                        hgvs_transcript
                    )
                )
                or (
                    hgvs_position_utils.start_is_5_prime_utr(
                        hgvs_transcript
                    )
                    and hgvs_position_utils.end_is_5_prime_utr(
                        hgvs_transcript
                    )
                )
        ):
            logger.info(
                "Translation passed into intronic handling code"
            )

            if (
                    (
                        1 <= pos.start.base <= 3
                        and not (
                            hgvs_position_utils
                            .start_position_is_intronic(
                                hgvs_transcript
                            )
                        )
                    )
                    or (
                        1 <= pos.end.base <= 3
                        and not (
                            hgvs_position_utils
                            .end_position_is_intronic(
                                hgvs_transcript
                            )
                        )
                        and (
                            hgvs_position_utils
                            .start_is_3_prime_utr(
                                hgvs_transcript
                            )
                        )
                        and (
                            hgvs_position_utils
                            .end_is_3_prime_utr(
                                hgvs_transcript
                            )
                        )
                    )
            ):
                residue_one = self.sf.fetch_seq(
                    associated_protein_accession,
                    start_i=0,
                    end_i=1,
                )

                hgvs_protein = _fb_unc(
                    associated_protein_accession,
                    residue_one,
                )

            else:
                hgvs_protein = _tot_unc(
                    associated_protein_accession
                )

            hgvs_transcript_to_hgvs_protein[
                "hgvs_protein"
            ] = hgvs_protein

            return hgvs_transcript_to_hgvs_protein

        logger.info(
            "Variant is not intronic and is not fully UTR, "
            "translate %s converted to %s",
            hgvs_transcript,
            hgvs_naughty,
        )

        inf = self.hdp.get_tx_identity_info(
            hgvs_transcript.ac
        )

        cds_start = inf[3]
        cds_end = inf[4]

        try:
            ref_seq = self.sf.fetch_seq(
                hgvs_naughty.ac
            )

        except Exception as e:
            hgvs_transcript_to_hgvs_protein[
                "error"
            ] = str(e)

            return hgvs_transcript_to_hgvs_protein

        var_seq = utils.n_inversion(
            ref_seq,
            del_seq,
            inv_seq,
            naughty_pos.start.base,
            naughty_pos.end.base,
        )

        logger.info(
            "Reference sequence:\n%s\n"
            "Deletion sequence:\n%s\n"
            "Inserted sequence:\n%s\n"
            "Var sequence:\n%s",
            ref_seq,
            del_seq,
            inv_seq,
            var_seq,
        )

        prot_seq = self.sf.fetch_seq(
            associated_protein_accession
        )

        if "U" in prot_seq:
            modified_aa = "Sec"

            hgvs_transcript_to_hgvs_protein["error"] = (
                "ProteinTranslationInfo: Selenocysteine detected "
                "in the original protein sequnce it may be "
                "incorporated instead of terminating at TGA/UGA "
                "termination codons"
            )

            logger.info(
                "Modified amino acid %s identified, "
                "update translation dict",
                modified_aa,
            )

        else:
            modified_aa = None

            logger.info(
                "No modified amino acid identified, "
                "use standard translation dict"
            )

        logger.info(
            "Translating reference and variant CDS outcomes"
        )

        try:
            prot_ref_seq = utils.translate(
                ref_seq,
                cds_start,
                modified_aa,
            )

        except IndexError:
            hgvs_transcript_to_hgvs_protein["error"] = (
                "ProteinTranslationError: Cannot generate a "
                "protein without an identifiable in-frame "
                "Termination codon in the reference mRNA sequence, "
                "this transcript may be subject to non-stop "
                "mediated decay"
            )

            hgvs_transcript_to_hgvs_protein[
                "hgvs_protein"
            ] = _tot_unc(
                associated_protein_accession
            )

            return hgvs_transcript_to_hgvs_protein

        except KeyError:
            hgvs_transcript_to_hgvs_protein["error"] = (
                "ProteinTranslationError: Unable to build protein "
                "sequence due to a non-CATG base included in the "
                "reference mRNA sequence, only standard "
                "unambiguous bases are accepted input for protein "
                "generation."
            )

            hgvs_transcript_to_hgvs_protein[
                "hgvs_protein"
            ] = _tot_unc(
                associated_protein_accession
            )

            return hgvs_transcript_to_hgvs_protein

        try:
            prot_var_seq = utils.translate(
                var_seq,
                cds_start,
                modified_aa,
            )

        except IndexError:
            hgvs_transcript_to_hgvs_protein["error"] = (
                "ProteinTranslationError: Cannot generate a "
                "protein without an identifiable in-frame "
                "Termination codon in the variant mRNA sequence, "
                "this transcript may be subject to non-stop decay"
            )

            hgvs_transcript_to_hgvs_protein[
                "hgvs_protein"
            ] = _tot_unc(
                associated_protein_accession
            )

            return hgvs_transcript_to_hgvs_protein

        except KeyError:
            hgvs_transcript_to_hgvs_protein["error"] = (
                "ProteinTranslationError: Unable to build protein "
                "sequence due to a non-CATG base included in the "
                "variant mRNA sequence, only standard unambiguous "
                "bases are accepted input for protein generation."
            )

            hgvs_transcript_to_hgvs_protein[
                "hgvs_protein"
            ] = _tot_unc(
                associated_protein_accession
            )

            return hgvs_transcript_to_hgvs_protein

        no_start_err = (
            "ProteinTranslationError: Unable to generate protein "
            "variant description due to the sequence missing an "
            "accepted start codon."
        )

        posedit = PosEdit(
            pos=Interval(),
            edit="?",
            uncertain=False,
        )

        hgvs_protein = (
            vvhgvs.sequencevariant.SequenceVariant(
                ac=associated_protein_accession,
                type="p",
                posedit=posedit,
            )
        )

        hgvs_transcript_to_hgvs_protein[
            "hgvs_protein"
        ] = hgvs_protein

        if prot_ref_seq == "error":
            hgvs_transcript_to_hgvs_protein[
                "error"
            ] = no_start_err.replace(
                "the sequence",
                "the reference sequence",
            )

            hgvs_transcript_to_hgvs_protein[
                "hgvs_protein"
            ] = _tot_unc(
                associated_protein_accession
            )

            return hgvs_transcript_to_hgvs_protein

        if prot_var_seq == "error":
            if (
                    (
                        1 <= pos.start.base <= 3
                        and not (
                            hgvs_position_utils
                            .start_position_is_intronic(
                                hgvs_transcript
                            )
                        )
                    )
                    or (
                        1 <= pos.end.base <= 3
                        and not (
                            hgvs_position_utils
                            .end_position_is_intronic(
                                hgvs_transcript
                            )
                        )
                    )
            ) and not (
                    hgvs_position_utils.start_is_3_prime_utr(
                        hgvs_transcript
                    )
                    or hgvs_position_utils.end_is_3_prime_utr(
                        hgvs_transcript
                    )
            ):
                residue_one = self.sf.fetch_seq(
                    associated_protein_accession,
                    start_i=0,
                    end_i=1,
                )

                hgvs_transcript_to_hgvs_protein[
                    "hgvs_protein"
                ] = _fb_unc(
                    associated_protein_accession,
                    residue_one,
                )

            else:
                hgvs_transcript_to_hgvs_protein[
                    "error"
                ] = no_start_err

            return hgvs_transcript_to_hgvs_protein

        if (
                (
                    1 <= pos.start.base <= 3
                    and not (
                        hgvs_position_utils
                        .start_position_is_intronic(
                            hgvs_transcript
                        )
                    )
                )
                or (
                    1 <= pos.end.base <= 3
                    and not (
                        hgvs_position_utils
                        .end_position_is_intronic(
                            hgvs_transcript
                        )
                    )
                )
        ) and not (
                hgvs_position_utils.start_is_3_prime_utr(
                    hgvs_transcript
                )
                or hgvs_position_utils.end_is_3_prime_utr(
                    hgvs_transcript
                )
        ):
            residue_one = self.sf.fetch_seq(
                associated_protein_accession,
                start_i=0,
                end_i=1,
            )

            hgvs_transcript_to_hgvs_protein[
                "hgvs_protein"
            ] = _fb_unc(
                associated_protein_accession,
                residue_one,
            )

            return hgvs_transcript_to_hgvs_protein

        if edit_type not in (
                "delins",
                "dup",
                "del",
                "ins",
        ):
            logger.info(
                "passing %s translations to pro_inv_info "
                "function",
                hgvs_transcript,
            )

            pro_inv_info = utils.pro_inv_info(
                prot_ref_seq,
                prot_var_seq,
            )

        else:
            logger.info(
                "passing %s translations to pro_delins_info "
                "function",
                hgvs_transcript,
            )

            cds_len = cds_end - cds_start
            minus = False
            plus = False

            if naughty_edit.type == "del":
                naughty_edit.alt = ""

            if naughty_edit.type == "ins":
                naughty_edit.ref = del_seq
                naughty_edit.alt = (
                    f"{del_seq[0]}"
                    f"{naughty_edit.alt}"
                    f"{del_seq[-1]}"
                )

            try:
                ref_len = len(naughty_edit.ref)
                alt_len = len(naughty_edit.alt)

                if ref_len > alt_len:
                    var_cds_len = (
                        cds_len - (ref_len - alt_len)
                    )
                    minus = True

                elif ref_len < alt_len:
                    var_cds_len = (
                        cds_len + (alt_len - ref_len)
                    )
                    plus = True

            except AttributeError as e:
                if (
                        "'Dup' object has no attribute 'alt'"
                        in str(e)
                ):
                    var_cds_len = (
                        cds_len
                        + (len(var_seq) - len(ref_seq))
                    )
                    plus = True

            in_frame = False

            if minus:
                loss_gain = cds_len - var_cds_len

                if loss_gain % 3 == 0:
                    in_frame = -(loss_gain / 3)

            elif plus:
                loss_gain = var_cds_len - cds_len

                if loss_gain % 3 == 0:
                    in_frame = loss_gain / 3

            pro_inv_info = utils.pro_delins_info(
                prot_ref_seq,
                prot_var_seq,
                in_frame,
            )

        logger.info(
            "RefSeq: %s",
            prot_ref_seq,
        )
        logger.info(
            "VarSeq: %s",
            prot_var_seq,
        )
        logger.info(
            "pro_inv_info: %s",
            pro_inv_info,
        )

        if (
                edit_type in ("del", "delins", "inv")
                and (
                    naughty_pos.start.base
                    < cds_end
                    <= naughty_pos.end.base
                )
                and pro_inv_info["prot_del_seq"][0] != "*"
        ):
            logger.info(
                "Variant %s starts upstream of the stop codon, "
                "and ends in or after the stop codon, could be "
                "a frame-shift",
                hgvs_transcript,
            )

            if (
                    (
                        "*" not in pro_inv_info["prot_ins_seq"]
                        and "*" not in pro_inv_info["prot_del_seq"]
                    )
                    or (
                        pro_inv_info["prot_ins_seq"][-1] == "*"
                        and pro_inv_info["prot_del_seq"][-1] == "*"
                    )
                    or edit_type == "del"
            ):
                ref = pro_inv_info["prot_del_seq"][0]
                alt = pro_inv_info["prot_ins_seq"][0]

                length = (
                    pro_inv_info["prot_ins_seq"].find("*")
                )
                length = length if length >= 0 else None

                if (
                        length is None
                        and pro_inv_info["terminate"] == "true"
                ):
                    length = (
                        pro_inv_info["ter_pos"]
                        - pro_inv_info["edit_start"]
                    )

                length += 1

                posedit = PosEdit(
                    pos=Interval(
                        start=AAPosition(
                            base=pro_inv_info["edit_start"],
                            aa=ref,
                        ),
                        end=AAPosition(
                            base=pro_inv_info["edit_start"],
                            aa=ref,
                        ),
                    ),
                    edit=AAFs(
                        ref=ref,
                        alt=alt,
                        length=length,
                    ),
                    uncertain=True,
                )

                hgvs_protein = (
                    vvhgvs.sequencevariant.SequenceVariant(
                        ac=associated_protein_accession,
                        type="p",
                        posedit=posedit,
                    )
                )

                hgvs_transcript_to_hgvs_protein[
                    "hgvs_protein"
                ] = hgvs_protein

                return hgvs_transcript_to_hgvs_protein

        if pro_inv_info["error"] == "true":
            hgvs_transcript_to_hgvs_protein["error"] = (
                "Translation error occurred, please contact admin"
            )

            return hgvs_transcript_to_hgvs_protein

        if pro_inv_info["variant"] != "true":
            posedit = VVPosEdit(
                pos=Interval(),
                edit=AARefAlt(),
                uncertain=True,
                nucleotide_not_equal=nucleotide_not_equal,
            )

            hgvs_protein = (
                vvhgvs.sequencevariant.SequenceVariant(
                    ac=associated_protein_accession,
                    type="p",
                    posedit=posedit,
                )
            )

            if (
                    isinstance(pos.start.base, int)
                    and isinstance(pos.end.base, int)
            ):
                aa_start_pos = (
                    pos.start.base + 2
                ) // 3

                aa_end_pos = (
                    pos.end.base + 2
                ) // 3

                aa_seq = self.sf.fetch_seq(
                    associated_protein_accession,
                    start_i=aa_start_pos - 1,
                    end_i=aa_end_pos,
                )

                if not aa_seq:
                    protein_sequence = self.sf.fetch_seq(
                        associated_protein_accession
                    )

                    if (
                            aa_start_pos
                            == len(protein_sequence) + 1
                            and aa_end_pos
                            == len(protein_sequence) + 1
                    ):
                        aa_seq = "*"

                start_aa = aa_seq[0]
                end_aa = aa_seq[-1]

                posedit = VVPosEdit(
                    pos=Interval(
                        start=AAPosition(
                            base=aa_start_pos,
                            aa=start_aa,
                        ),
                        end=AAPosition(
                            base=aa_end_pos,
                            aa=end_aa,
                        ),
                    ),
                    edit=AARefAlt(),
                    uncertain=True,
                    nucleotide_not_equal=nucleotide_not_equal,
                )

                hgvs_protein = (
                    vvhgvs.sequencevariant.SequenceVariant(
                        ac=associated_protein_accession,
                        type="p",
                        posedit=posedit,
                    )
                )

            hgvs_transcript_to_hgvs_protein[
                "hgvs_protein"
            ] = hgvs_protein

            return hgvs_transcript_to_hgvs_protein

        if (
                modified_aa == "Sec"
                and "U" in pro_inv_info["prot_ins_seq"]
                and "U" not in pro_inv_info["prot_del_seq"]
        ):
            pro_inv_info["ter_pos"] = (
                pro_inv_info["edit_start"]
                + len(pro_inv_info["prot_ins_seq"])
            )

        posedit = False

        if (
                pro_inv_info["terminate"] == "true"
                and edit_type in (
                    "delins",
                    "dup",
                    "inv",
                    "ins",
                )
        ):
            frameshift = False

            if naughty_edit.type == "dup":
                length = (
                    naughty_pos.end.base
                    - naughty_pos.start.base
                    + 1
                )
                frameshift = length % 3 != 0

            elif naughty_edit.type == "del":
                frameshift = (
                    len(naughty_edit.ref or "") % 3 != 0
                )

            elif naughty_edit.type == "ins":
                frameshift = (
                    len(naughty_edit.alt or "") % 3 != 0
                )

            elif naughty_edit.type == "delins":
                frameshift = (
                    len(naughty_edit.alt or "")
                    - len(naughty_edit.ref or "")
                ) % 3 != 0

            if frameshift:
                ref = pro_inv_info["prot_del_seq"][0]
                alt = pro_inv_info["prot_ins_seq"][0]

                length = (
                    pro_inv_info["prot_ins_seq"].find("*")
                )
                length = length if length >= 0 else None

                posedit = PosEdit(
                    pos=Interval(
                        start=AAPosition(
                            base=pro_inv_info["edit_start"],
                            aa=ref,
                        ),
                        end=AAPosition(
                            base=pro_inv_info["edit_start"],
                            aa=ref,
                        ),
                    ),
                    edit=AAFs(
                        ref=ref,
                        alt=alt,
                        length=length,
                    ),
                    uncertain=True,
                )

                hgvs_protein = (
                    vvhgvs.sequencevariant.SequenceVariant(
                        ac=associated_protein_accession,
                        type="p",
                        posedit=posedit,
                    )
                )

                hgvs_transcript_to_hgvs_protein[
                    "hgvs_protein"
                ] = hgvs_protein

                return hgvs_transcript_to_hgvs_protein

            if (
                    len(pro_inv_info["prot_del_seq"])
                    + pro_inv_info["edit_start"] - 1
                    == pro_inv_info["ter_pos"]
            ):
                pro_inv_info["prot_del_seq"] = (
                    pro_inv_info["prot_del_seq"][0]
                )
                pro_inv_info["edit_end"] = (
                    pro_inv_info["edit_start"]
                )

            elif (
                    edit_type == "dup"
                    and not pro_inv_info["prot_del_seq"]
                    and (
                        pro_inv_info["edit_end"]
                        < pro_inv_info["edit_start"]
                    )
            ):
                dup_len = (
                    pos.end.base
                    - pos.start.base
                    + 1
                ) / 3

                pro_inv_info["prot_del_seq"] = (
                    pro_inv_info["prot_ins_seq"]
                )

                pro_inv_info["edit_start"] = (
                    pro_inv_info["edit_end"]
                    - len(pro_inv_info["prot_del_seq"])
                    + 1
                )

                start_aa = self.sf.fetch_seq(
                    associated_protein_accession,
                    int(pro_inv_info["edit_start"] - 1),
                    int(
                        pro_inv_info["edit_start"]
                        + dup_len - 1
                    ),
                )

                pro_inv_info["prot_del_seq"] = start_aa
                pro_inv_info["prot_ins_seq"] = (
                    start_aa
                    + pro_inv_info["prot_ins_seq"]
                )

        prot_del_seq = pro_inv_info["prot_del_seq"]
        prot_ins_seq = pro_inv_info["prot_ins_seq"]
        edit_start = pro_inv_info["edit_start"]
        edit_end = pro_inv_info["edit_end"]

        if not prot_del_seq:
            assert edit_start != edit_end
            from_aa = prot_ref_seq[edit_start]
            to_aa = prot_ref_seq[edit_end]

        else:
            from_aa = prot_del_seq[0]
            to_aa = prot_del_seq[-1]

        if edit_start != edit_end:
            if prot_ins_seq == prot_del_seq + prot_del_seq:
                posedit = VVPosEdit(
                    pos=Interval(
                        start=AAPosition(
                            base=edit_start,
                            aa=from_aa,
                        ),
                        end=AAPosition(
                            base=edit_end,
                            aa=to_aa,
                        ),
                    ),
                    edit=Dup(
                        ref=prot_del_seq,
                    ),
                    uncertain=True,
                    nucleotide_not_equal=nucleotide_not_equal,
                )

            elif (
                    len(prot_ins_seq) > len(prot_del_seq)
                    and prot_ins_seq
                    != prot_del_seq + prot_del_seq
                    and not prot_del_seq
                    and edit_start > edit_end
            ):
                from_aa = self.sf.fetch_seq(
                    associated_protein_accession,
                    edit_end - len(prot_ins_seq),
                    edit_end - len(prot_ins_seq) + 1,
                )

                to_aa = self.sf.fetch_seq(
                    associated_protein_accession,
                    edit_start - 2,
                    edit_start - 1,
                )

                posedit = VVPosEdit(
                    pos=Interval(
                        start=AAPosition(
                            base=(
                                edit_end
                                - len(prot_ins_seq)
                                + 1
                            ),
                            aa=from_aa,
                        ),
                        end=AAPosition(
                            base=edit_start - 1,
                            aa=to_aa,
                        ),
                    ),
                    edit=Dup(
                        ref=prot_del_seq,
                    ),
                    uncertain=True,
                    nucleotide_not_equal=nucleotide_not_equal,
                )

            elif prot_ins_seq:
                if (
                        "*" in prot_del_seq
                        and not prot_ins_seq.endswith("*")
                ):
                    posedit = VVPosEdit(
                        pos=Interval(
                            start=AAPosition(
                                base=edit_start,
                                aa=from_aa,
                            ),
                            end=AAPosition(
                                base=edit_end,
                                aa=to_aa,
                            ),
                        ),
                        edit=AARefAlt(
                            ref="",
                            alt=prot_ins_seq + "?",
                        ),
                        uncertain=True,
                        nucleotide_not_equal=nucleotide_not_equal,
                    )

                else:
                    posedit = VVPosEdit(
                        pos=Interval(
                            start=AAPosition(
                                base=edit_start,
                                aa=from_aa,
                            ),
                            end=AAPosition(
                                base=edit_end,
                                aa=to_aa,
                            ),
                        ),
                        edit=AARefAlt(
                            ref="",
                            alt=prot_ins_seq,
                        ),
                        uncertain=True,
                        nucleotide_not_equal=nucleotide_not_equal,
                    )

            else:
                posedit = VVPosEdit(
                    pos=Interval(
                        start=AAPosition(
                            base=edit_start,
                            aa=from_aa,
                        ),
                        end=AAPosition(
                            base=edit_end,
                            aa=to_aa,
                        ),
                    ),
                    edit=AARefAlt(
                        ref=prot_del_seq,
                        alt=None,
                    ),
                    uncertain=True,
                    nucleotide_not_equal=nucleotide_not_equal,
                )

        else:
            if prot_ins_seq == prot_del_seq + prot_del_seq:
                posedit = VVPosEdit(
                    pos=Interval(
                        start=AAPosition(
                            base=edit_start,
                            aa=from_aa,
                        )
                    ),
                    edit=Dup(
                        ref=from_aa,
                    ),
                    uncertain=True,
                    nucleotide_not_equal=nucleotide_not_equal,
                )

            elif (
                    len(prot_ins_seq) > len(prot_del_seq)
                    and prot_ins_seq
                    != prot_del_seq + prot_del_seq
                    and prot_ins_seq.startswith(prot_del_seq[0])
            ):
                to_aa = self.sf.fetch_seq(
                    associated_protein_accession,
                    edit_start,
                    edit_start + 1,
                )

                posedit = VVPosEdit(
                    pos=Interval(
                        start=AAPosition(
                            base=edit_start,
                            aa=from_aa,
                        ),
                        end=AAPosition(
                            base=edit_end + 1,
                            aa=to_aa,
                        ),
                    ),
                    edit=AARefAlt(
                        alt=prot_ins_seq[1:],
                    ),
                    uncertain=True,
                    nucleotide_not_equal=nucleotide_not_equal,
                )

            elif (
                    prot_del_seq == "*"
                    and len(prot_ins_seq) > len(prot_del_seq)
            ):
                if prot_ins_seq.endswith("*"):
                    posedit = VVPosEdit(
                        pos=Interval(
                            start=AAPosition(
                                base=edit_start,
                                aa=from_aa,
                            )
                        ),
                        edit=AAExt(
                            alt=prot_ins_seq[0],
                            length=len(prot_ins_seq) - 1,
                            aaterm="*",
                        ),
                        uncertain=True,
                        nucleotide_not_equal=nucleotide_not_equal,
                    )

                else:
                    posedit = VVPosEdit(
                        pos=Interval(
                            start=AAPosition(
                                base=edit_start,
                                aa=from_aa,
                            )
                        ),
                        edit=AAExt(
                            alt=prot_ins_seq[-1],
                            length="?",
                        ),
                        uncertain=True,
                    )

            elif len(prot_ins_seq) == 1:
                posedit = VVPosEdit(
                    pos=Interval(
                        start=AAPosition(
                            base=edit_start,
                            aa=from_aa,
                        )
                    ),
                    edit=AARefAlt(
                        alt=prot_ins_seq,
                        ref=from_aa,
                    ),
                    uncertain=True,
                    nucleotide_not_equal=nucleotide_not_equal,
                )

            elif not prot_ins_seq:
                posedit = VVPosEdit(
                    pos=Interval(
                        start=AAPosition(
                            base=edit_start,
                            aa=from_aa,
                        )
                    ),
                    edit=AARefAlt(
                        ref=prot_del_seq,
                    ),
                    uncertain=True,
                    nucleotide_not_equal=nucleotide_not_equal,
                )

            else:
                posedit = VVPosEdit(
                    pos=Interval(
                        start=AAPosition(
                            base=edit_start,
                            aa=from_aa,
                        )
                    ),
                    edit=AARefAlt(
                        alt=prot_ins_seq,
                        ref=from_aa,
                    ),
                    uncertain=True,
                    nucleotide_not_equal=nucleotide_not_equal,
                )

        hgvs_protein = (
            vvhgvs.sequencevariant.SequenceVariant(
                ac=associated_protein_accession,
                type="p",
                posedit=posedit,
            )
        )

        hgvs_transcript_to_hgvs_protein[
            "hgvs_protein"
        ] = hgvs_protein

        return hgvs_transcript_to_hgvs_protein

    def revcomp(self, bases):
        """
        Return the reverse complement of a nucleotide sequence.
        """
        return utils.simple_dna_revcomp(bases)


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later