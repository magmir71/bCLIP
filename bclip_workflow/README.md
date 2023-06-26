# bCLIP_workflow

This workflow is used for the analysis of the clip variation protocol, called bCLIP.
The following operations are performed on the specified samples:
- extraction of UMIs
- demultiplexing of fastq files
- removal of poly(A) tails
- creation of STAR index
- mapping of reads to genome with STAR
- sorting and indexing of `.bam` files
- collapsing reads according to UMIs
- sorting and indexing of deduplicated `.bam` files

Samples present in the multiplexed `.fastq` files that are not specified in `samples.tsv` will not be processed after demultiplexing. A list of those "abandonned" samples will be written to a text file. 

## Installation 
### 1. Clone the repository

Go to the desired directory/folder on your file system, then clone/get the 
repository and move into the respective directory with:

```bash
git clone ssh://git@git.scicore.unibas.ch:2222/katsanto/bclip_workflow.git
cd bclip_workflow
```

### 2. Conda installation

Workflow dependencies can be conveniently installed with the [Conda][conda]
package manager. We recommend that you install [Miniconda][miniconda-installation] 
for your system (Linux). Be sure to select Python 3 option. 

### 3. Dependencies installation

For improved reproducibility and reusability of the workflow,
each individual step of the workflow runs in its own [Singularity][singularity]container. 
As a consequence, running this workflow has very few individual dependencies. 
If you want to make use of **container execution**, please [install
Singularity][singularity-install] in privileged mode, depending
on your system. You may have to ask an authorized person (e.g., a systems
administrator) to do that. This will almost certainly be required if you want to run the workflow on a high-performance computing (HPC) cluster. 

After installing Singularity, or should you choose not to use containerization, install the remaining dependencies with:
```bash
conda env create -f install/environment.yml
```

### 4. Activate environment

Activate the Conda environment with:

```bash
conda activate bclip
```

## Execute Workflow
### 1. Configuration
Check the files `config.yml` and `samples.tsv` and adjust filepaths and parameters (Refer to comments in `config.yml` for some explanation).

### 2. Run
To run the pipeline on the [Slurm][slurm] compute environment, using Singularity containers, navigate to the root directory of this repository and start the provided bash script:
```bash
bash run_slurm_singularity.sh
```
Adapt this script accordingly if you want to run the workflow in a different environment.


[conda]: <https://docs.conda.io/projects/conda/en/latest/index.html>
[miniconda-installation]: <https://docs.conda.io/en/latest/miniconda.html>
[rule-graph]: images/dag.svg
[snakemake]: <https://snakemake.readthedocs.io/en/stable/>
[singularity]: <https://sylabs.io/singularity/>
[singularity-install]: <https://sylabs.io/guides/3.8/user-guide/quick_start.html>
[slurm]: <https://slurm.schedmd.com/documentation.html>