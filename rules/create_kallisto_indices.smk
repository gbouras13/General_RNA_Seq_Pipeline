"""
Snakefile for building Kalliso index for encode v39

snakemake -s rules/create_kallisto_indices.smk --use-conda --config Salmon_dir='/hpcfs/users/a1667917/Salmon_Ref_Genomes' Kallisto_dir="/hpcfs/users/a1667917/Kallisto_Ref_Genomes"

snakemake -c 1 -s rules/create_kallisto_indices.smk --use-conda --conda-create-envs-only --conda-frontend conda

"""

# https://ycl6.gitbook.io/guide-to-rna-seq-analysis/preparations/softwares-and-databases

import os

# load default config
configfile: os.path.join(workflow.basedir, '..', 'config', 'config.yaml')

BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]

if config['Salmon_dir'] is None:
    Salmon_dir = "/hpcfs/users/a1667917/Salmon_Ref_Genomes"
else:
    Salmon_dir = config["Salmon_dir"]

if config['Kallisto_dir'] is None:
    Kallisto_dir = "/hpcfs/users/a1667917/Kallisto_Ref_Genomes"
else:
    Kallisto_dir = config["Kallisto_dir"]


# needs to be created before (should exist)
if not os.path.exists(os.path.join(Salmon_dir)):
  os.makedirs(os.path.join(Salmon_dir))

if not os.path.exists(os.path.join(Kallisto_dir)):
  os.makedirs(os.path.join(Kallisto_dir))

rule all:
    input:
        os.path.join(Kallisto_dir, 'kallisto_index.flag')

rule kallisto_index:
    """create kallisto index."""
    input:
        os.path.join(Salmon_dir, 'gencode.v39.transcripts.fa')
    output:
        os.path.join(Kallisto_dir, 'kallisto_index.flag')
    params:
        os.path.join(Kallisto_dir, 'gencode.v39.idx')
    threads:
        BigJobCpu
    conda:
        os.path.join('..', 'envs','kallisto.yaml')
    resources:
        mem_mb=BigJobMem
    shell:
        """
        kallisto index -i {params[0]} {input[0]}
        touch {output[0]}
        """
