# APTARS - Analysis of Pacbio TARgeted Sequencing

<!-- badges: start -->
![Maintainer](https://img.shields.io/badge/maintainer-SidSethi-blue)
[![Generic badge](https://img.shields.io/badge/WMS-snakemake-blue.svg)]
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)]
[![Linux](https://svgshare.com/i/Zhy.svg)]
[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)]
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/sid-sethi/APTARS/blob/main/LICENSE)
<!-- badges: end -->

APTARS is a `snakemake` pipeline that takes PacBio subreads as input, generates consensus reads using CCS, demultiplexes using lima, refines and cluster using isoseq3, map the reads to the genome using minimap2, assembles gene transcripts using cDNA_cupcake and annotates transcripts using Sqant3. Below is the dag of the pipeline:  

<p align="center">
  <img src="dag/dag.png" width="300" height="800"/>  
</p>


# Getting Started

## Input

- PacBio subreads.bam
- Primers in fasta format
- Reference genome assembly in fasta format
- Gencode gtf
- Cage peaks
- Intropolis
- Poly(a) atlas

## Depedencies

- [miniconda](https://conda.io/miniconda.html)
- SQANTI3: v4.2 - https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v4.2.tar.gz
- The rest of the dependencies (including snakemake) are installed via conda through the `environment.yml` file

**SQANTI3 setup:**  
```bash
wget https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v4.2.tar.gz
tar -xvf v4.2.tar.gz
echo "  - cdna_cupcake=22.0.0" >> SQANTI3-4.2/SQANTI3.conda_env.yml
```


## Installation

Clone the directory:

```bash
git clone --recursive https://github.com/sid-sethi/APTARS.git
```

Create conda environment for the pipeline which will install all the dependencies:

```bash
cd APTARS
conda env create -f environment.yml
```

## Usage

Edit `config.yml` to set up the working directory and input files/directories. `snakemake` command should be issued from within the pipeline directory. Please note that before you run any of the `snakemake` commands, make sure to first activate the conda environment using the command `conda activate aptars`.

```bash
cd APTARS
conda activate aptars
snakemake --use-conda -j <num_cores> all
```
It is a good idea to do a dry run (using -n parameter) to view what would be done by the pipeline before executing the pipeline.

```bash
snakemake --use-conda -n all
```

You can visualise the processes to be executed in a DAG:

```bash
snakemake --dag | dot -Tpng > dag.png
```

To exit a running `snakemake` pipeline, hit `ctrl+c` on the terminal. If the pipeline is running in the background, you can send a `TERM` signal which will stop the scheduling of new jobs and wait for all running jobs to be finished.

```bash
killall -TERM snakemake
```

To deactivate the conda environment:
```bash
conda deactivate
```
