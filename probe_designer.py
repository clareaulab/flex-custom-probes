import functools
from pathlib import Path
import subprocess

import numpy as np
from scipy.optimize import dual_annealing
from Bio.Seq import reverse_complement

from utils import max_homopolymer_length, probe_blast_search, Probe, setup_blast

# BLAST constants
COLLECTED_FASTA_PATH = "collected.fasta"


class ProbeDesigner:

    def __init__(self,
                 working_dir: str,
                 reference_database: dict[str, str],
                 lhs_probe_length: int = 25,
                 rhs_probe_length: int = 25,
                 invalid_score: float = 1e7):
        self.LHS_PROBE_LENGTH = lhs_probe_length
        self.RHS_PROBE_LENGTH = rhs_probe_length
        self.INVALID_SCORE = invalid_score

        self.working_dir = Path(working_dir)
        self.working_dir.mkdir(exist_ok=True)
        self.BLAST_DB_PATH = self.working_dir / "blast_db"
        setup_blast(self.BLAST_DB_PATH, reference_database)

    def test_probe(self,
                   args: np.ndarray,
                   name: str,
                   sequence: str,
                   visited: dict[tuple[int, int], float] = None,
                   strict_gc: bool = False) -> float:
        """
        Calculate a score for a given probe. Note the following requirements:
            - Avoid overlap with annotated repeat or low complexity sequences  # TODO: if run locally, run bwa to annotate low complexity
            - Avoid homopolymer repeats
            - Each probe should have GC content between 44-72%
            - Matches to off-target genes should have at least five mismatches

        :param args: The arguments to be optimized for the probe (LHS relative start position, RHS relative start position).
        :param name: The name of the gene.
        :param sequence: The original gene sequence.
        :param visited: The visited positions to avoid.
        :param strict_gc: If True, the GC content of the probes must be between 44-72%.
        :return: The score to minimize for the probe.
        """
        if visited is None:
            visited = dict()

        lhs_start = args[0]  # The start of the left-hand side probe
        lhs_start = int(np.round(lhs_start))
        lhs_end = lhs_start + self.LHS_PROBE_LENGTH  # The end of the left-hand side probe

        lhs_sequence = reverse_complement(sequence[lhs_start:lhs_end])  # The sequence of the left-hand side probe

        rhs_start = lhs_end

        rhs_end = rhs_start + self.RHS_PROBE_LENGTH  # The end of the right-hand side probe

        if rhs_end > len(sequence):
            visited[(lhs_start, rhs_start)] = self.INVALID_SCORE
            return self.INVALID_SCORE

        rhs_sequence = reverse_complement(sequence[rhs_start:rhs_end])  # The sequence of the right-hand side probe

        if (lhs_start, rhs_start) in visited:
            return visited[(lhs_start, rhs_start)]

        score = 0.0  # The score for the probe, start at perfect score and penalize as necessary

        # Test if low complexity sequence is already annotated by counting lower case letters
        n_lower_lhs = sum([1 for c in lhs_sequence if c.islower()])
        n_lower_rhs = sum([1 for c in rhs_sequence if c.islower()])
        score += n_lower_lhs * 10  # Penalize for each lower case letter
        score += n_lower_rhs * 10

        # the RHS probe must start with an A
        if rhs_sequence[-1] != "T":
            score += self.INVALID_SCORE

        # Calculate the GC content of each probe
        if len(lhs_sequence) == 0:
            lhs_gc = 0
        else:
            lhs_gc = (lhs_sequence.count("G") + lhs_sequence.count("C")) / len(lhs_sequence)
        if len(rhs_sequence) == 0:
            rhs_gc = 0
        else:
            rhs_gc = (rhs_sequence.count("G") + rhs_sequence.count("C")) / len(rhs_sequence)
        if lhs_gc < 0.44:
            if strict_gc:
                score += self.INVALID_SCORE
            else:
                score += (0.44 - lhs_gc) * 100  # Penalize for each percentage point below 44%
        elif lhs_gc > 0.72:
            if strict_gc:
                score += self.INVALID_SCORE
            else:
                score += (lhs_gc - 0.72) * 100
        if rhs_gc < 0.44:
            if strict_gc:
                score += self.INVALID_SCORE
            else:
                score += (0.44 - rhs_gc) * 100
        elif rhs_gc > 0.72:
            if strict_gc:
                score += self.INVALID_SCORE
            else:
                score += (rhs_gc - 0.72) * 100

        # Hard cutoffs regardless:
        if lhs_gc < 0.2:
            score += self.INVALID_SCORE
        elif lhs_gc > 0.8:
            score += self.INVALID_SCORE
        if rhs_gc < 0.2:
            score += self.INVALID_SCORE
        elif rhs_gc > 0.8:
            score += self.INVALID_SCORE

        # We want GC to be as close as possible between the two probes
        score += abs((lhs_gc - rhs_gc)*100)

        # Finally, identify off-target genes and penalize for each mismatch by using a BLAST search
        # But BLAST searches are slow, so we will only do this if the probe passes the other tests
        if score >= self.INVALID_SCORE:  # If the probe is already invalid, don't bother with the BLAST search
            visited[(lhs_start, rhs_start)] = score
            return score

        # Calculate the maximum homopolymer length of each probe
        lhs_homopolymer = max_homopolymer_length(lhs_sequence)
        rhs_homopolymer = max_homopolymer_length(rhs_sequence)

        score += lhs_homopolymer * 10  # Penalize for each base in the homopolymer
        score += rhs_homopolymer * 10

        # TODO: Calculate the repeats and low complexity sequences

        # When performing a BLAST search, we need to reverse complement the probes to match the transcriptome
        print(f"Found valid candidate ({lhs_start}, {rhs_start}), performing BLAST search...")
        n_hits_lhs, n_hits_rhs = probe_blast_search(self.BLAST_DB_PATH,
                                                    name,
                                                    reverse_complement(lhs_sequence),
                                                    reverse_complement(rhs_sequence))
        # We expect at least one hit (the gene we are targeting)
        print("LHS hits", n_hits_lhs, "RHS hits", n_hits_rhs)
        if n_hits_lhs < 1 or n_hits_rhs < 1 or n_hits_lhs > 10 or n_hits_rhs > 10:
            score += self.INVALID_SCORE
        else:
            # Calculate a penalty for each off-target hit
            score += (100 ** (n_hits_lhs - 1))
            score += (100 ** (n_hits_rhs - 1))

        visited[(lhs_start, rhs_start)] = score
        return score

    def create_probes(self,
                      name: str,
                      sequence: str,
                      n_probes: int = 1) -> list[Probe]:
        """
        Create potential probes for a sequence.
        :param name: The name of the gene.
        :param sequence: The transcript sequence being targeted.
        :param n_probes: The ideal number of probes to create.
        :return: A list of potential probes. If none are found, the list will be empty.
        """
        sequence = sequence
        probes = []
        visited = dict()
        for i in range(n_probes):
            max_start = len(sequence) - self.LHS_PROBE_LENGTH - self.RHS_PROBE_LENGTH - 1
            if max_start <= 0:
                return []
            result = dual_annealing(func=functools.partial(self.test_probe, strict_gc=True),
                                    bounds=[(0, max_start)],
                                    x0=np.array([len(sequence) // 2]),
                                    args=(name, sequence, visited),
                                    maxiter=1_000//n_probes,
                                    initial_temp=500
                                    )
            result_args = result.x
            result_score = result.fun

            # Note that we invert lhs and rhs since we were looking at lhs from the perspective of the transcript
            # But the probes are reverse complemented
            if result_args is None:
                return []
            rhs_start, = np.round(result_args).astype(int)
            visited.update({(int(rhs_start), int(rhs_start+self.LHS_PROBE_LENGTH)): self.INVALID_SCORE*10 for p in probes})
            lhs_start = rhs_start + self.LHS_PROBE_LENGTH
            lhs_probe = reverse_complement(sequence[lhs_start:lhs_start + self.LHS_PROBE_LENGTH])
            rhs_probe = reverse_complement(sequence[rhs_start:rhs_start + self.RHS_PROBE_LENGTH])
            lhs_gc = (lhs_probe.count("G") + lhs_probe.count("C")) / len(lhs_probe)
            rhs_gc = (rhs_probe.count("G") + rhs_probe.count("C")) / len(rhs_probe)
            probes.append(Probe(name, lhs_probe, rhs_probe,
                                reverse_complement(lhs_probe), reverse_complement(rhs_probe), result_score, lhs_gc, rhs_gc))

        # Remove invalid probes
        probes = [p for p in probes if p.score < self.INVALID_SCORE]
        return probes
