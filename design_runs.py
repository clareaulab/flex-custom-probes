import pandas as pd
from pyensembl import EnsemblRelease
from Bio import SeqIO

from utils import transcriptome, export_probes_to_xlsx
from probe_designer import ProbeDesigner


def design_car_probes(human_background=False):
    transcriptome_dict = dict()
    transcriptome_fasta = "CAR/WL_CARconstruct.fasta"
    for record in SeqIO.parse(transcriptome_fasta, "fasta"):
        transcriptome_dict[record.description] = str(record.seq)
    design_me = transcriptome_dict.copy()

    if human_background:
        ensembl_genome = EnsemblRelease(111)
        ensembl_genome.download(overwrite=False)
        ensembl_genome.index(overwrite=False)
        transcriptome_dict = transcriptome_dict | transcriptome(ensembl_genome)
	
    probe_design = ProbeDesigner(working_dir="CAR", reference_database=transcriptome_dict)
    n_probes = 5
    results = []
    missing_probes = []
    for name, sequence in design_me.items():
        probes = probe_design.create_probes(name, sequence, n_probes=n_probes)
        if probes is None or len(probes) < 1:
            missing_probes.append(name)
            print(f"MISSING: {name}")
            continue
        print(f"Found {len(probes)} probes for {name}.")
        results += probes

    # Add Results to DataFrame
    df = export_probes_to_xlsx(results, f"HHV6_probes_with_human_background_info.xlsx", n_barcodes=4)
    # Add truseq barcodes
    df['lhs_TruSeq_probe'] = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT" + df['lhs_seq']
    df.to_excel(f"CAR/WL_CAR_probes_with_human_background_info.xlsx", index=False)
    df = export_probes_to_xlsx(results, f"CAR/WL_CARprobes_with_human_background_1barcode_info.xlsx", n_barcodes=1)
    # Add truseq barcodes
    df['lhs_TruSeq_probe'] = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT" + df['lhs_seq']
    df.to_excel(f"CAR/WL_CAR_probes_with_human_background_1barcode_info.xlsx", index=False)
    print(f"Missing {len([r for r in results if r is None])} probes.")

    return pd.DataFrame(results)


if __name__ == "__main__":
    design_car_probes(human_background=False)
