import copy
import logging

logger = logging.getLogger(__name__)

"""
A module for holding the TranscriptMapData specific to an individual transcript
in order to avoid having to repeatedly re-pull from the database.
"""
class TranscriptMapData():
    """
    This object is used to hold the transcript mapping data used in all tx
    variants, but in a way that does not force VariantFormatter to pull in the
    whole VV variant object. We may add more such functions later to reduce the
    amount sundry transcript mapping objects passed around in bits.
    """
    def __init__(self,hdp=None):
        """
        Initialisation of TranscriptMapData object, can keep the last validator
        hgvs data provider, to avoid constantly having to re-specify or pass
        around the validator separately.
        """
        self.hdp = hdp # current validator hgvs data provider in use
        self.mapping_opts = {} # used to only fetch once from vv/uta per tx
        self.gap_status = {} # made using above map data, for dict fetch
        self.mapped_strands = {} # made using above map data, for dict fetch
        self.mapping_types = {} # the mapping type for each tx->alt map
        self.exon_data = {} # the set of exon data per mapping
        self._known_bad_alignments = [
            ["NM_001009944.3", "NT_187607.1"]
        ]

    def mapping_options(self,tx_ac,hdp=None):
        """
        Get mapping options for the main variant transcript.

        These were originally fetched repeatedly from the database, this
        function centralises the fetch and caches the first result, it can't
        change without a server restart anyway.
        """
        if tx_ac in self.mapping_opts:
            return self.mapping_opts[tx_ac]
        if hdp:
            cur_hdp = hdp
            self.hdp = hdp
        elif self.hdp:
            cur_hdp = self.hdp
        else:
            raise ValueError(
                    "Either the creation, or the first data fetch, of each"
                    "TranscriptMapData object must specify a hgvs data "
                    "provider (hdp) for use as a data source")
        self.mapping_opts[tx_ac] = cur_hdp.get_tx_mapping_options(
                tx_ac,gap_warn=True)

        self.mapping_opts[tx_ac] = [
            option
            for option in self.mapping_opts[tx_ac]
            if [option[0], option[1]] not in self._known_bad_alignments
        ]

        logger.info(f"Mapping options for {tx_ac}: {self.mapping_opts[tx_ac]}")
        return self.mapping_opts[tx_ac]

    def map_strand(self,tx_ac,alt_ac,hdp=None):
        """
        Return the strand for a transcript, as mapped to a specific alt
        accession to allow strand fetch without exon fetch, or repeated
        database queries.
        """
        if tx_ac in self.mapped_strands:
            if alt_ac in self.mapped_strands[tx_ac]:
                return self.mapped_strands[tx_ac][alt_ac]
            return False
        map_opts = self.mapping_options(tx_ac,hdp)
        # opts layout is tx_ac,alt_ac,alt_aln_method, gapped, strand
        for opt in map_opts:
            if not opt[0] in self.mapped_strands:
                self.mapped_strands[opt[0]] = {}
            self.mapped_strands[opt[0]][opt[1]] = opt[4]

        if alt_ac in self.mapped_strands[tx_ac]:
            return self.mapped_strands[tx_ac][alt_ac]
        return False

    def is_gapped_map(self,tx_ac,alt_ac,hdp=None):
        """
        Returns true if the mapping in question is gapped, otherwise false.
        Allows us to do without the fixed bad gene list (which also impacts
        un-gapped mappings of problem transcripts.
        """
        if tx_ac in self.gap_status:
            if alt_ac in self.gap_status[tx_ac] \
                    and self.gap_status[tx_ac][alt_ac]:
                return True
            return False
        self.gap_status[tx_ac] = {}
        map_opts = self.mapping_options(tx_ac,hdp)
        for opt in map_opts:
            self.gap_status[tx_ac][opt[1]] = opt[3]
        if alt_ac in self.gap_status[tx_ac]:
            return self.gap_status[tx_ac][alt_ac]
        return False

    def map_type(self,tx_ac,alt_ac,hdp=None):
        """
        Look up the map type of the requested tx->alt ac pair, prefer not blat
        """
        if tx_ac in self.mapping_types:
            if alt_ac not in self.mapping_types[tx_ac]:
                return []
            return self.mapping_types[tx_ac][alt_ac]
        self.mapping_types[tx_ac] = {}
        map_opts = self.mapping_options(tx_ac,hdp)
        for opt in map_opts:
            if opt[2] == 'blat' and opt[1] in self.mapping_types[tx_ac]:
                continue
            self.mapping_types[tx_ac][opt[1]] = opt[2]
        if alt_ac not in self.mapping_types[tx_ac]:
            return False
        return self.mapping_types[tx_ac][alt_ac]

    def mapped_exons(self,tx_ac,alt_ac,alt_aln_method=None,hdp=None):
        """
        Return exon set data for mapping in question, without repeated db fetch.
        This version just returns the cached results in the same order as
        fetched, and presumes that the data is not edited on use in problematic
        ways.
        """
        if hdp:
            cur_hdp = hdp
            self.hdp = hdp
        elif self.hdp:
            cur_hdp = self.hdp
        else:
            raise ValueError(
                    "Either the creation, or the first data fetch, of each"
                    "TranscriptMapData object must specify a hgvs data "
                    "provider (hdp) for use as a data source")
        if tx_ac not in self.exon_data:
            self.exon_data[tx_ac] = {}
        if alt_ac not in self.exon_data[tx_ac]:
            if alt_aln_method:
                aln_method = alt_aln_method
            else:
                aln_method = self.map_type(tx_ac,alt_ac)
            self.exon_data[tx_ac][alt_ac] = cur_hdp.get_tx_exons(
                    tx_ac,alt_ac,aln_method)
        if alt_ac not in self.exon_data[tx_ac]:
            return []
        return self.exon_data[tx_ac][alt_ac]

    def tx_exons(self, tx_ac, alt_ac, alt_aln_method, hdp=None):
        """
        Alt strand flipped exon data function, cached version of orig validator
        function, some users of this method seem to edit this data.
        Returns exon information for a given transcript
        e.g. how the exons align to the genomic reference
        see vvhgvs.dataproviders.uta.py for details
        """
        # Interface with the UTA database via get_tx_exons in uta.py
        tx_exons = self.mapped_exons(
                tx_ac, alt_ac,
                hdp=hdp,
                alt_aln_method=alt_aln_method)
        tx_exons = copy.copy(tx_exons)
        try:
            tx_exons[0]['alt_strand']
        except TypeError:
            tx_exons = 'error'
            return tx_exons
        # If on the reverse strand, reverse the order of elements
        if tx_exons[0]['alt_strand'] == -1:
            tx_exons = tx_exons[::-1]
            logger.debug(f"Exon data for {tx_ac}: {tx_exons}")
            return tx_exons
        else:
            logger.debug(f"Exon data for {tx_ac}: {tx_exons}")
            return tx_exons

