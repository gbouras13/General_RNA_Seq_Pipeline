"""
Snakefile for building Salmon index for encode v39

snakemake -s rules/create_salmon_indices.smk --use-conda --config HG38_dir="/hpcfs/users/a1667917/STAR_Ref_Genomes" Salmon_dir='/hpcfs/users/a1667917/Salmon_Ref_Genomes'

snakemake -c 1 -s rules/create_salmon_indices.smk --use-conda --conda-create-envs-only --conda-frontend conda

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

if config['HG38_dir'] is None:
    HG38_dir = "/hpcfs/users/a1667917/STAR_Ref_Genomes"
else:
    HG38_dir = config["HG38_dir"]


# needs to be created before (should exist)
if not os.path.exists(os.path.join(Salmon_dir)):
  os.makedirs(os.path.join(Salmon_dir))

rule all:
    input:
        os.path.join(Salmon_dir, 'GRCh38.decoys.txt'),
        os.path.join(Salmon_dir, 'GRCh38.gentrome.fa.gz'),
        os.path.join(Salmon_dir, 'salmon_index.flag'),
        os.path.join(Salmon_dir, 'star_transcripts.fa')

rule decoys:
    """get decoys"""
    input:
        os.path.join(Salmon_dir, 'GRCh38.primary_assembly.genome.fa')
    output:
       os.path.join(Salmon_dir, 'GRCh38.decoys.txt')
    threads:
        1
    shell:
        """
        grep "^>" {input[0]} | cut -d " " -f 1 > {output[0]}
        sed -i 's/>//g' {output[0]}
        """

rule concat:
    """concat transcript and fasta"""
    input:
        os.path.join(Salmon_dir, 'gencode.v39.transcripts.fa'),
        os.path.join(Salmon_dir, 'GRCh38.primary_assembly.genome.fa')
    output:
       os.path.join(Salmon_dir, 'GRCh38.gentrome.fa.gz')
    threads:
        1
    shell:
        """
        cat {input[0]} {input[1]} | gzip > {output[0]}
        """

rule salmon_index:
    """create salmon index."""
    input:
        os.path.join(Salmon_dir, 'GRCh38.gentrome.fa.gz'),
        os.path.join(Salmon_dir, 'GRCh38.decoys.txt')
    output:
        os.path.join(Salmon_dir, 'salmon_index.flag')
    params:
        os.path.join(Salmon_dir, 'gencode.v39_decoys_salmon')
    threads:
        BigJobCpu
    conda:
        os.path.join('..', 'envs','salmon.yaml')
    resources:
        mem_mb=BigJobMem
    shell:
        """
        salmon index -p {threads} --gencode -t {input[0]}  -d {input[1]} -i {params[0]} 
        touch {output[0]}
        """

rule transcript_for_star_salmon:
    """create salmon index."""
    input:
        os.path.join(HG38_dir, 'GRCh38.primary_assembly.genome.fa'),
        os.path.join(HG38_dir, 'gencode.v39.primary_assembly.annotation.gtf')
    output:
        os.path.join(Salmon_dir, 'star_transcripts.fa')
    threads:
        BigJobCpu
    conda:
        os.path.join('..', 'envs','salmon.yaml')
    resources:
        mem_mb=BigJobMem
    shell:
        """
        gffread -w {output[0]} -g {input[0]} {input[1]}
        """



