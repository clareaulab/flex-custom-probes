import itertools
import os
import subprocess
from collections import namedtuple
from io import StringIO

import pandas as pd
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from pyensembl import Genome, Transcript


Probe = namedtuple("Probe", ["name", "lhs", "rhs", "lhs_gene", "rhs_gene", "score", "lhs_GC_content", "rhs_GC_content"])


# From table 2: https://kb.10xgenomics.com/hc/en-us/articles/17623693026445-How-do-I-design-custom-probes-for-a-Single-Cell-Gene-Expression-Flex-i-e-Fixed-RNA-Profiling-for-multiplexed-samples-experiment
barcode_seqs = [
    "ACTTTAGG",
    "AACGGGAA",
    "AGTAGGCT",
    "ATGTTGAC",
    "ACAGACCT",
    "ATCCCAAC",
    "AAGTAGAG",
    "AGCTGTGA",
    "ACAGTCTG",
    "AGTGAGTG",
    "AGAGGCAA",
    "ACTACTCA",
    "ATACGTCA",
    "ATCATGTG",
    "AACGCCGA",
    "ATTCGGTT"
]


def parse_chromosome_location(hg38_location_string: str) -> tuple[str, int, int]:
    """
    Parse the chromosome location from the input string.
    :param hg38_location_string: The hg38 location string (ex. 9:133256215-133256215).
    :return: The chromosome, hg38 start, and hg38 end.
    """
    hg38_location_string = hg38_location_string.replace(" ", "").replace("..", "-")
    split_string = hg38_location_string.split(":")
    chromosome = split_string[0]
    location_range = split_string[1].split("-")
    hg38_start = int(location_range[0])
    hg38_end = int(location_range[1]) if len(location_range) == 2 else hg38_start

    return chromosome, hg38_start, hg38_end


def max_homopolymer_length(sequence: str) -> int:
    groups = itertools.groupby(sequence)
    return max(len(list(group)) for key, group in groups)


def setup_blast(blast_db: str, sequences: dict[str, str]):
    """
    Set up the BLAST database for a genome. Note that the name of each sequence should begin with the gene name.
    :param sequences: The sequences to add to the BLAST database. Keys are the names of the sequences. Values are the sequences.
    """
    with open("temp.fasta", "w") as fasta_file:
        for name, sequence in sequences.items():
            fasta_file.write(f">{name}\n{sequence}\n")

    subprocess.run(["makeblastdb", "-in", "temp.fasta", "-dbtype", "nucl", "-out", str(blast_db)])
    os.remove("temp.fasta")


def blast_search(blast_db: str,
                 sequence: str,
                 evalue: float = 1.0) -> NCBIXML:
    """
    Perform a BLAST search to identify the number of genes this matches to.
    :param blast_db: The BLAST database path to search.
    :param sequence: The sequence to search.
    :param evalue: The maximum E-value to consider.
    :return: The BLAST results.
    """
    blastn_command = ["blastn", "-db", str(blast_db), "-evalue", str(evalue), "-outfmt", "5", "-num_threads", "4", "-task", "blastn-short"]
    blast_out = subprocess.run(blastn_command, input=sequence, text=True, stdout=subprocess.PIPE).stdout
    if len(blast_out) == 0:
        return None

    # Parse the BLAST results
    blast_hits = NCBIXML.read(StringIO(blast_out))

    return blast_hits


def probe_blast_search(blast_db: str,
                       name: str,
                       lhs_sequence: str,
                       rhs_sequence: str,
                       evalue: float = 1.0
                       ) -> tuple[int, int]:
    """
    Perform a BLAST search to identify the number of genes this matches to.
    :param blast_db: The BLAST database path to search.
    :param name: The name of the gene the probe is targeting.
    :param lhs_sequence: The left-hand side probe sequence.
    :param rhs_sequence: The right-hand side probe sequence.
    :param evalue: The maximum E-value to consider.
    :return: The number of hits for the LHS and RHS probes.
    """
    lhs_hits = blast_search(blast_db, lhs_sequence, evalue)
    rhs_hits = blast_search(blast_db, rhs_sequence, evalue)

    if lhs_hits is None or rhs_hits is None:
        return 0, 0

    # 10X recommends at least 5 mismatches for a probe to be considered unique, so select alignments with < 5 mismatches to be "hits"
    # Also compute hits wrt gene name only
    n_lhs_hits = len(set([aln.hit_def.split(" ")[1] for aln in lhs_hits.alignments if sum([hsp.identities for hsp in aln.hsps]) >= len(lhs_sequence) - 5]))
    n_rhs_hits = len(set([aln.hit_def.split(" ")[1] for aln in rhs_hits.alignments if sum([hsp.identities for hsp in aln.hsps]) >= len(rhs_sequence) - 5]))

    return n_lhs_hits, n_rhs_hits


def rank_and_filter_transcripts(transcripts: list[Transcript], use_support=True) -> list[Transcript]:
    """
    Rank the transcripts by support level and length. Additionally only retain complete protein coding transcripts.
    :param transcripts: The transcripts to rank.
    :param use_support: Whether to use the support level in ranking.
    :return: The ranked transcripts.
    """
    filtered = list(filter(lambda x: x.complete, sorted(transcripts, key=lambda x: (x.length, -(x.support_level or 0) if use_support else 1), reverse=True)))
    return filtered


def transcriptome(genome: Genome) -> dict[str, str]:
    """
    Select the largest transcript for each gene to create a non-redundant transcriptome.
    :param genome: The genome to work with.
    :return: The  transcriptome. Key is transcript ID, value is the coding sequence.
    """
    transcriptome = dict()
    for gene in genome.genes():
        transcripts = rank_and_filter_transcripts(gene.transcripts)
        if len(transcripts) == 0:
            continue
        for transcript in transcripts:
            gene_name = gene.gene_name
            transcriptome[gene_name + " " + transcript.transcript_id] = transcript.coding_sequence
    return transcriptome


def export_probes_to_xlsx(probes: list[Probe], filename: str, probe_order_format: bool = False, n_barcodes: int = 1) -> pd.DataFrame:
    results = dict(
        name=[],
        lhs_probe=[],
        rhs_probe=[],
        lhs_seq=[],
        rhs_seq=[],
        lhs_gene=[],
        rhs_gene=[],
        lhs_GC_content=[],
        rhs_GC_content=[],
        GC_difference=[],
        score=[],
        barcode=[]
    )
    # Note that the LHS
    rhs_prefix = "/5Phos/"
    rhs_probe_barcode_bridge = "ACGCGGTTAGCACGTANN"
    rhs_suffix = "CGGTCCTAGCAA"

    lhs_prefix = "CCTTGGCACCCGAGAATTCCA"

    for probe in probes:
        if probe is None:
            continue
        for barcode in barcode_seqs[:n_barcodes]:
            results["name"].append(probe.name)
            results["lhs_probe"].append(lhs_prefix + probe.lhs)
            results["rhs_probe"].append(rhs_prefix + probe.rhs + rhs_probe_barcode_bridge + barcode + rhs_suffix)
            results["lhs_seq"].append(probe.lhs)
            results["rhs_seq"].append(probe.rhs)
            results["lhs_gene"].append(probe.lhs_gene)
            results["rhs_gene"].append(probe.rhs_gene)
            results["lhs_GC_content"].append(probe.lhs_GC_content)
            results["rhs_GC_content"].append(probe.rhs_GC_content)
            results["GC_difference"].append(abs(probe.lhs_GC_content - probe.rhs_GC_content))
            results["score"].append(probe.score)
            results["barcode"].append(barcode)

    df = pd.DataFrame(results)
    if not probe_order_format:
        df.to_excel(filename, index=False)
    else:
        exported_df = dict(
            name=[],
            value=[]
        )
        for i, row in df.iterrows():
            for j, barcode in enumerate(barcode_seqs[:n_barcodes]):
                exported_df['name'].append(row['name'] + " " + f"#{barcode}")
                exported_df['value'].append("5'-3'")

                exported_df["name"].append("RHSprobe")
                exported_df["value"].append(rhs_prefix + row['rhs_seq'] + rhs_probe_barcode_bridge + barcode + rhs_suffix)

                exported_df["name"].append("LHSprobe")
                exported_df["value"].append(lhs_prefix + row['lhs_seq'])
        pd.DataFrame(exported_df).to_excel(filename, index=False, header=False)
    return df
