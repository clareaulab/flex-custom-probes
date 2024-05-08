# flex-custom-probes
Public facing repository to aid in Flex probe design

## Requirements

The following Python packages are required:
- scipy
- numpy
- pyensembl
- pandas

Additionally, NCBI BLAST+ must be installed on the system. 

Using a conda environment, the requirements (including BLAST) can be installed with: `conda install -c conda-forge -c bioconda scipy numpy pandas pyensembl blast`

## Usage
The `design_runs.py` script provides a usage example of running the pipeline to generate probes.
