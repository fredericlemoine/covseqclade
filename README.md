# CovSeqClade


## Introduction

This [Nextflow](https://www.nextflow.io/) workflow takes as input paired-end SARS-CoV-2 fastq files (gz or [dsrc](https://pubmed.ncbi.nlm.nih.gov/24747219/) compressed). It computes the relative coverage on all the mutations that are specific to all Nextclade clades.

For each sample, it is thus easy to visually check:

1) The probable virus clade
2) Cases of co-infections with several clades
3) Cases of potential recombination between several clades


## Installation

### Dependencies

To run this workflow, you need to install:

1) Java
2) Nextflow
3) Singularity

### Configuration

You should also update the `nextflow.config` to match your HPC plateform configuration.

1) For slurm: You just have to change the variables  `slurmqueue`, `slurmqos`.
2) For running locally: You just have to set the variable `executor='local'`, in the `process`section.
3) For other HPC platforms, have a look at the [nextflow documentation](https://www.nextflow.io/docs/latest/executor.html). 

## Usage

To display workflow help:

```
nextflow run main.nf --help
```

To run the workflow, type:

```
nextflow run main.nf --results <result directory> \
                     --fastqs <input data directory> \
                     --clade1 <First clade to compare|-1> \
                     --clade2 <Second clade to compare|-1>
```

If clade1 or clade2 = -1, then the coverage on all clades are displayed. Otherwise, only clade1 and clade2 are compared.

## Examples

- Co-infection visualized with clade1=-1 and clade2=-1

![Example of Co-infection](imgs/sample1_coinfection.svg)

- Co-infection visualized with clade1=21J and clade2=21K

![Example of Co-infection](imgs/sample1_coinfection_pair.svg)

- Recombination visualized with clade1=-1 and clade2=-1

![Example of recombinant](imgs/sample2_recombination.svg)

- Recombination visualized with clade1=21J and clade2=21K

![Example of recombinant](imgs/sample2_recombination_pair.svg)
