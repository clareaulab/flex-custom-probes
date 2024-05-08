import os
from typing import Literal

import pandas as pd
from pyensembl import EnsemblRelease

from utils import transcriptome, rank_and_filter_transcripts
from probe_designer import *


class HumanProbeDesigner(ProbeDesigner):
    """
    A wrapper for ProbeDesigner that uses the human genome as the target.
    """

    def __init__(self, *args, **kwargs):
        # Latest Ensembl version for hg38
        ENSEMBL_RELEASE = 111
        os.environ['PYENSEMBL_CACHE_DIR'] = './ensembl_cache'
        # Overrides for canonical isoform info
        overrides = pd.read_csv(
            "https://raw.githubusercontent.com/genome-nexus/genome-nexus-importer/master/data/common_input/isoform_overrides_at_mskcc_grch38.txt",
            sep="\t")
        ensembl_genome = EnsemblRelease(ENSEMBL_RELEASE)
        ensembl_genome.download(overwrite=False)
        ensembl_genome.index(overwrite=False)
        kwargs['reference_database'] = transcriptome(ensembl_genome)  # ensembl_genome.transcript_sequences.fasta_dictionary
        super().__init__(*args, **kwargs)
        self.ensembl_genome = ensembl_genome
        self.overrides = overrides

    def create_human_probes(self,
                            gene_name: str,
                            n_probes: int = 1) -> list[Probe]:
        """
        Create probes for the human genome.
        :param gene_name: The name of the gene.
        :param n_probes: The ideal number of probes to create.
        :return: A list of probes. If none are found, an empty list is returned.
        """

        try:
            strand, original_sequence = self.get_gene_sequence(gene_name)
        except ValueError as e:
            print(e)
            return []

        return self.create_probes(
            name=gene_name,
            sequence=original_sequence,
            n_probes=n_probes
        )

    def get_gene_sequence(self,
                          gene_name: str) -> tuple[str, str]:
        """
        Get the sequence of a gene from the hg38 genome.
        :param chromosome: The chromosome the gene is on.
        :param gene_name: The name of the gene.
        :return: The strand, gene sequence.
        """
        gene_objs = self.ensembl_genome.genes_by_name(gene_name)
        if len(gene_objs) != 1:
            print(f"WARNING: Found {len(gene_objs)} genes for {gene_name}")
        gene_obj = gene_objs[0]
        transcript_candidates = rank_and_filter_transcripts(
            gene_obj.transcripts
        )
        # transcript_candidates = self.ensembl_genome.transcripts_at_locus(contig=chromosome, position=hg38_start, end=hg38_end)

        if len(transcript_candidates) == 0:
            raise ValueError("No transcripts found for the specified gene.")

        transcript = None
        # Prioritize msk canonical transcript, else use support level and length
        # ,
        for candidate in sorted(transcript_candidates, key=lambda x: (1e8 if self.overrides.enst_id.isin([x.transcript_id]).any() else 0, -(x.support_level or 0), x.length), reverse=True):  # Sort by support level and length
            if candidate.biotype == 'protein_coding' and candidate.complete:
                transcript = candidate
                original_sequence = transcript.coding_sequence

                original_sequence = [char for char in original_sequence if char.isalpha()]
                original_sequence = "".join(original_sequence)
                return transcript.strand, original_sequence

        raise ValueError("No correct protein coding transcripts found for the specified gene.")
